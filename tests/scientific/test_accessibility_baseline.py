"""RED (Task 8): leakage-safe no-learning accessibility baseline.

Ranks cysteines by the structural accessibility proxy extracted in
``plantpersulf.features.structure`` (``contact_number_proxy`` — a labelled
contact-count proxy, not a rigorous SASA/RSA calculation). Standardisation
statistics (mean/std) are fit on train-split, structure-observed cysteines
only, exactly like the motif baseline's train-only discipline in
``tests/scientific/test_baselines.py``. Cysteines with no structure
(``has_structure=False``) are never imputed with an average contact value —
they must always be assigned the same fixed "missing" score and therefore
rank last, never in the middle of the distribution by chance.

Reasoning for score direction: the literature this project builds on treats
solvent-exposed cysteines as more chemically reactive to H2S-derived
persulfidation signalling. A *lower* contact-number proxy means fewer
neighbouring residues near a Cys's C-alpha, i.e. more exposure — so the
baseline's score must increase as contact_number_proxy decreases.

Expected RED: ``train_accessibility_baseline`` /
``evaluate_accessibility_ranking`` do not exist yet in
``plantpersulf.models.baselines``.
"""

from __future__ import annotations

from pathlib import Path

from plantpersulf.features.structure import (
    CysStructureFeature,
    write_structure_features_tsv,
)
from plantpersulf.models.baselines import (  # RED: symbols missing
    AccessibilityBaseline,
    evaluate_accessibility_ranking,
    train_accessibility_baseline,
)


def _write_fixtures(tmp_path: Path) -> tuple[Path, Path]:
    splits = tmp_path / "splits.csv"
    splits.write_text(
        "protein_accession,cys_position,label,split\n"
        "A,3,positive,train\n"
        "B,2,unlabeled,train\n"
        "C,4,unlabeled,train\n"
        "D,5,positive,validation\n"
        "E,6,positive,test\n"
        "F,7,unlabeled,test\n",
        encoding="utf-8",
    )
    structure_features = tmp_path / "structure_features.tsv"
    write_structure_features_tsv(
        [
            # train: exposed (low contact) positive
            CysStructureFeature("A", 3, True, 90.0, False, 1.0),
            # train: buried (high contact) unlabeled
            CysStructureFeature("B", 2, True, 88.0, False, 5.0),
            # train: no structure at all (isoform / 404 / mapping failure)
            CysStructureFeature("C", 4, False, None, None, None),
            # validation: exposed
            CysStructureFeature("D", 5, True, 91.0, False, 0.0),
            # test: buried positive
            CysStructureFeature("E", 6, True, 40.0, True, 6.0),
            # test: no structure
            CysStructureFeature("F", 7, False, None, None, None),
        ],
        structure_features,
    )
    return splits, structure_features


def test_accessibility_baseline_fits_only_on_train_structured_rows(
    tmp_path: Path,
) -> None:
    splits_path, structure_path = _write_fixtures(tmp_path)
    model = train_accessibility_baseline(splits_path, structure_path)

    assert isinstance(model, AccessibilityBaseline)
    # Only A (1.0) and B (5.0) are train rows with has_structure=True;
    # C is train but has no structure and must be excluded from fitting.
    assert model.train_mean == 3.0
    assert model.train_std > 0.0


def test_more_exposed_cysteine_scores_higher(tmp_path: Path) -> None:
    splits_path, structure_path = _write_fixtures(tmp_path)
    model = train_accessibility_baseline(splits_path, structure_path)

    exposed = CysStructureFeature("X", 1, True, 90.0, False, 0.0)
    buried = CysStructureFeature("X", 2, True, 90.0, False, 10.0)
    assert model.score(exposed) > model.score(buried)


def test_missing_structure_scores_at_the_fixed_floor_never_imputed(
    tmp_path: Path,
) -> None:
    splits_path, structure_path = _write_fixtures(tmp_path)
    model = train_accessibility_baseline(splits_path, structure_path)

    missing = CysStructureFeature("Y", 1, False, None, None, None)
    present = CysStructureFeature("Y", 2, True, 90.0, False, 100.0)  # very buried
    assert model.score(missing) == model.missing_score
    assert model.score(missing) < model.score(present)


def test_evaluate_accessibility_ranking_matches_split_schema(
    tmp_path: Path,
) -> None:
    splits_path, structure_path = _write_fixtures(tmp_path)
    model = train_accessibility_baseline(splits_path, structure_path)

    for split_name in ("train", "validation", "test"):
        result = evaluate_accessibility_ranking(
            model, splits_path, structure_path, split_name
        )
        assert result.split_name == split_name
        assert result.total_sites > 0


def test_evaluate_accessibility_ranking_recovers_exposed_positive_first(
    tmp_path: Path,
) -> None:
    splits_path, structure_path = _write_fixtures(tmp_path)
    model = train_accessibility_baseline(splits_path, structure_path)

    result = evaluate_accessibility_ranking(model, splits_path, structure_path, "train")
    assert result.recall_at_k is not None
    assert result.recall_at_k[1] == 1.0  # A (exposed, positive) ranks #1


def test_accessibility_baseline_is_deterministic(tmp_path: Path) -> None:
    splits_path, structure_path = _write_fixtures(tmp_path)
    first = train_accessibility_baseline(splits_path, structure_path)
    second = train_accessibility_baseline(splits_path, structure_path)

    assert first.train_mean == second.train_mean
    assert first.train_std == second.train_std
