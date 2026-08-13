"""Gate 0 (Task 1): structure_ranker serialisable bundle — round-trip equality.

The release package requires a deterministic save/load path so blind-time
scoring (after Gate 0 freeze) can re-score future rows from a frozen state
dict without ever refitting. These tests pin the contract:

- ``fit_structure_ranker`` returns a ``StructureRankerBundle``;
- ``bundle.save`` / ``StructureRankerBundle.load`` round-trip losslessly;
- scoring from a loaded bundle is bit-identical to scoring from a fresh fit
  with the same seed (and identical to the legacy ``structure_ranker_scores``
  wrapper), so the frozen model and its fit-time scores are one and the same.
"""

from __future__ import annotations

import pytest

from plantpersulf.models.structure_ranker import (
    AblationConfig,
    BranchFeatures,
    StructureRankerBundle,
    fit_structure_ranker,
    score_structure_ranker_bundle,
    structure_ranker_scores,
)


def _toy_branches(n: int, *, mask_all: bool = True) -> BranchFeatures:
    # Two separable clusters so a model can actually learn something.
    sequence = [[float(i % 2), float((i + 1) % 2)] for i in range(n)]
    esm = [[float(i % 2)] * 4 for i in range(n)]
    structure = [[float(i % 2) * 10.0, 50.0 + (i % 2) * 40.0] for i in range(n)]
    structure_mask = [mask_all or (i % 2 == 0) for i in range(n)]
    study_ids = ["PXD000001" if i % 2 == 0 else "PXD000002" for i in range(n)]
    return BranchFeatures(
        sequence=sequence,
        esm=esm,
        structure=structure,
        structure_mask=structure_mask,
        study_ids=study_ids,
    )


def _toy_labels(n: int) -> list[str]:
    return ["positive" if i % 2 == 0 else "unlabeled" for i in range(n)]


def test_bundle_roundtrip_scores_match_fresh_fit(tmp_path) -> None:
    train = _toy_branches(20)
    train_y = _toy_labels(20)
    predict = _toy_branches(7)

    fresh = score_structure_ranker_bundle(
        fit_structure_ranker(train, train_y, seed=0), predict
    )

    bundle = fit_structure_ranker(train, train_y, seed=0)
    path = tmp_path / "ranker.pt"
    bundle.save(path)
    loaded = StructureRankerBundle.load(path)
    restored = score_structure_ranker_bundle(loaded, predict)

    assert restored.scores == pytest.approx(fresh.scores, abs=0.0)
    assert restored.uncertainty == pytest.approx(fresh.uncertainty, abs=0.0)


def test_bundle_score_equals_legacy_wrapper(tmp_path) -> None:
    train = _toy_branches(20)
    train_y = _toy_labels(20)
    predict = _toy_branches(7)

    legacy = structure_ranker_scores(train, train_y, predict, seed=1)

    bundle = fit_structure_ranker(train, train_y, seed=1)
    path = tmp_path / "ranker.pt"
    bundle.save(path)
    loaded = StructureRankerBundle.load(path)
    bundled = score_structure_ranker_bundle(loaded, predict)

    assert bundled.scores == pytest.approx(legacy.scores, abs=0.0)
    assert bundled.uncertainty == pytest.approx(legacy.uncertainty, abs=0.0)


def test_bundle_preserves_frozen_metadata(tmp_path) -> None:
    train = _toy_branches(20)
    train_y = _toy_labels(20)

    bundle = fit_structure_ranker(
        train,
        train_y,
        seed=42,
        ablation=AblationConfig(use_esm=False, use_study_context=False),
        hidden=8,
        dropout=0.3,
        epochs=50,
        lr=0.1,
        holdout_fraction=0.25,
        n_mc_dropout=8,
    )
    path = tmp_path / "ranker.pt"
    bundle.save(path)
    loaded = StructureRankerBundle.load(path)

    assert loaded.seed == 42
    assert loaded.ablation == AblationConfig(
        use_esm=False, use_study_context=False
    )
    assert loaded.input_dims == {"d_seq": 2, "d_esm": 4, "d_str": 2, "d_study": 2}
    assert loaded.hyperparameters == {
        "hidden": 8,
        "dropout": 0.3,
        "epochs": 50,
        "lr": 0.1,
        "holdout_fraction": 0.25,
        "n_mc_dropout": 8,
    }
    assert set(loaded.net_state_dict) == {
        "enc_seq.weight",
        "enc_seq.bias",
        "enc_esm.weight",
        "enc_esm.bias",
        "enc_str.weight",
        "enc_str.bias",
        "enc_study.weight",
        "enc_study.bias",
        "gate.weight",
        "gate.bias",
        "head.weight",
        "head.bias",
    }
    # Scaler statistics must survive the pickle round-trip.
    assert loaded.scalers.seq.mean is not None


def test_two_fits_same_seed_are_bit_identical() -> None:
    train = _toy_branches(20)
    train_y = _toy_labels(20)
    predict = _toy_branches(7)

    first = structure_ranker_scores(train, train_y, predict, seed=7)
    second = structure_ranker_scores(train, train_y, predict, seed=7)

    assert first.scores == second.scores
    assert first.uncertainty == second.uncertainty


def test_bundle_rejects_mismatched_feature_dims(tmp_path) -> None:
    train = _toy_branches(20)
    train_y = _toy_labels(20)
    bundle = fit_structure_ranker(train, train_y, seed=0)
    path = tmp_path / "ranker.pt"
    bundle.save(path)
    loaded = StructureRankerBundle.load(path)

    bad_predict = _toy_branches(3)
    bad_predict = BranchFeatures(
        sequence=[[0.0, 1.0, 0.0] for _ in range(3)],
        esm=bad_predict.esm,
        structure=bad_predict.structure,
        structure_mask=bad_predict.structure_mask,
        study_ids=bad_predict.study_ids,
    )
    with pytest.raises(ValueError, match="dims do not match"):
        score_structure_ranker_bundle(loaded, bad_predict)


def test_bundle_score_empty_predict_returns_empty(tmp_path) -> None:
    train = _toy_branches(20)
    train_y = _toy_labels(20)
    bundle = fit_structure_ranker(train, train_y, seed=0)
    path = tmp_path / "ranker.pt"
    bundle.save(path)
    loaded = StructureRankerBundle.load(path)

    out = score_structure_ranker_bundle(loaded, _toy_branches(0))
    assert out.scores == []
    assert out.uncertainty == []


def test_load_rejects_corrupt_bundle_file(tmp_path) -> None:
    path = tmp_path / "junk.pt"
    path.write_text("not a torch bundle")

    with pytest.raises(RuntimeError, match="corrupt"):
        StructureRankerBundle.load(path)
