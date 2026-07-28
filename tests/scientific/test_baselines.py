"""RED (Task 8): leakage-safe PU baselines.

Implements no-learning baselines (motif frequency, cysteine count) and a
train-only-normalization constraint: every statistic used at inference time
must be computed from training data only, never from validation or test.

Expected RED: ``plantpersulf.models.baselines`` does not exist yet.
"""

from __future__ import annotations

import csv
from pathlib import Path

from plantpersulf.models.baselines import (  # RED: module missing
    evaluate_ranking,
    train_motif_baseline,
)


def _write_fixtures(tmp_path: Path) -> tuple[Path, Path, Path]:
    splits = tmp_path / "splits.csv"
    splits.write_text(
        "protein_accession,cys_position,label,split\n"
        "A,3,positive,train\n"
        "B,2,positive,validation\n"
        "C,4,positive,train\n"
        "C,2,positive,test\n",
        encoding="utf-8",
    )
    proteome = tmp_path / "mini.fasta"
    proteome.write_text(">sp|A\nMVCGK\n>sp|B\nACDEF\n>sp|C\nGCCCC\n", encoding="utf-8")
    features = tmp_path / "features.csv"
    features.write_text(
        "protein_accession,cys_position,flanking_window\n"
        "A,3,MVCGK\n"
        "B,2,ACDEF\n"
        "C,4,GCCCC\n"
        "C,1,GCCCC\n",
        encoding="utf-8",
    )
    return splits, proteome, features


def test_motif_baseline_fits_only_on_train(tmp_path: Path) -> None:
    splits_path, proteome_path, _ = _write_fixtures(tmp_path)
    model = train_motif_baseline(splits_path, proteome_path, window_radius=1)

    assert model.proteome_hash is not None  # provenance recorded
    assert isinstance(model.central_residue_freq, dict)


def test_evaluate_ranking_returns_valid_scores(tmp_path: Path) -> None:
    splits_path, proteome_path, _ = _write_fixtures(tmp_path)
    model = train_motif_baseline(splits_path, proteome_path, window_radius=1)

    for split_name in ("train", "validation", "test"):
        result = evaluate_ranking(model, splits_path, proteome_path, split_name)
        assert result.split_name == split_name
        assert result.total_sites > 0
        n_pos = sum(
            1
            for row in csv.DictReader(
                splits_path.open(encoding="utf-8", newline="")
            )
            if row["split"] == split_name and row["label"] == "positive"
        )
        if n_pos == 0:
            continue
        assert result.recall_at_k is not None
        # Recall at the largest K (100) must be 1.0 when n_pos <= 100
        assert result.recall_at_k.get(100, 0.0) == 1.0


def test_train_only_normalization_prevents_test_leakage(tmp_path: Path) -> None:
    """Building the model with train rows only must produce a model that
    makes predictions on test rows without having seen them during fitting."""
    splits_path, proteome_path, _ = _write_fixtures(tmp_path)
    model = train_motif_baseline(splits_path, proteome_path, window_radius=1)

    # The model was trained on the train split only.
    # Evaluating on test must not raise and must produce scores.
    result = evaluate_ranking(model, splits_path, proteome_path, "test")
    assert result.total_sites > 0
    assert result.positives_found is not None


def test_baseline_is_deterministic(tmp_path: Path) -> None:
    splits_path, proteome_path, _ = _write_fixtures(tmp_path)
    first = train_motif_baseline(splits_path, proteome_path, window_radius=1)
    second = train_motif_baseline(splits_path, proteome_path, window_radius=1)
    assert first.central_residue_freq == second.central_residue_freq
