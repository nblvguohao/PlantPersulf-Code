"""RED (Task 9): structure-aware PU ranker — output shapes and value ranges.

Exercises the minimal gated-fusion ranker on small synthetic multi-branch
"features" (sequence / frozen-ESM / structure + a missingness mask). Pure
algorithm correctness, not a biological claim — the real frozen features are
tested end-to-end elsewhere.

Expected RED: ``plantpersulf.models.structure_ranker`` does not exist yet.
"""

from __future__ import annotations

import pytest

from plantpersulf.models.structure_ranker import (  # RED: module missing
    AblationConfig,
    BranchFeatures,
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


def test_returns_one_score_and_one_uncertainty_per_predict_row() -> None:
    train = _toy_branches(20)
    train_y = _toy_labels(20)
    predict = _toy_branches(6)

    out = structure_ranker_scores(train, train_y, predict, seed=0)

    assert len(out.scores) == 6
    assert len(out.uncertainty) == 6
    for score in out.scores:
        assert isinstance(score, float)
    for unc in out.uncertainty:
        assert isinstance(unc, float)
        assert unc >= 0.0


def test_rejects_empty_training_positives() -> None:
    train = _toy_branches(4)
    all_unlabeled = ["unlabeled"] * 4

    with pytest.raises(ValueError, match="positive"):
        structure_ranker_scores(train, all_unlabeled, _toy_branches(2), seed=0)


def test_rejects_mismatched_branch_lengths() -> None:
    bad = BranchFeatures(
        sequence=[[0.0, 1.0], [1.0, 0.0]],
        esm=[[0.0]],  # length mismatch
        structure=[[0.0, 0.0], [1.0, 1.0]],
        structure_mask=[True, True],
        study_ids=["A", "B"],
    )
    with pytest.raises(ValueError, match="length"):
        structure_ranker_scores(bad, ["positive", "unlabeled"], bad, seed=0)


def test_ablation_disabling_esm_still_produces_valid_scores() -> None:
    train = _toy_branches(20)
    train_y = _toy_labels(20)
    predict = _toy_branches(6)

    out = structure_ranker_scores(
        train, train_y, predict, seed=0,
        ablation=AblationConfig(use_esm=False),
    )
    assert len(out.scores) == 6
