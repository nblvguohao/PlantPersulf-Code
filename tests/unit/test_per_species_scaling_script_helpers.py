"""Unit tests for the per-species structure-scaling diagnostic's helpers."""

from __future__ import annotations

from scripts.accuracy_campaign.eval_per_species_structure_scaling import (
    _species_of,
    per_species_apk_and_recall,
)


def test_species_of_takes_prefix_before_pipe() -> None:
    assert _species_of("arabidopsis|A0A068FPW3") == "arabidopsis"
    assert _species_of("tomato|P12345") == "tomato"
    assert _species_of("no-pipe") == "no-pipe"


def test_per_species_apk_and_recall_splits_by_species() -> None:
    rows = [
        ("arabidopsis|A", 1, "positive"),
        ("arabidopsis|B", 2, "unlabeled"),
        ("rice|C", 3, "positive"),
        ("rice|D", 4, "unlabeled"),
        ("tomato|E", 5, "positive"),
    ]
    scores = [0.9, 0.5, 0.8, 0.4, 0.7]  # positive first within each species
    out = per_species_apk_and_recall(rows, scores, k=200)
    # Arabidopsis: 1 positive / 2 rows, ranked [pos, unlab] -> recall 1.0.
    assert out["arabidopsis"]["recall_k"] == 1.0
    assert out["rice"]["recall_k"] == 1.0
    assert out["tomato"]["recall_k"] == 1.0
    assert out["arabidopsis"]["n_positive"] == 1


def test_per_species_recall_k_respects_top_k_cutoff() -> None:
    rows = [
        ("arabidopsis|A", 1, "positive"),
        ("arabidopsis|B", 2, "positive"),
        ("arabidopsis|C", 3, "unlabeled"),
    ]
    scores = [0.1, 0.2, 0.9]  # only the unlabeled ranks above both positives
    out = per_species_apk_and_recall(rows, scores, k=2)
    # Top-2 are [unlabeled(0.9), positive B(0.2)] -> 1/2 positives reached.
    assert out["arabidopsis"]["recall_k"] == 0.5


def test_per_species_skips_species_without_positives() -> None:
    rows = [
        ("arabidopsis|A", 1, "unlabeled"),
        ("arabidopsis|B", 2, "unlabeled"),
    ]
    scores = [0.5, 0.6]
    assert per_species_apk_and_recall(rows, scores) == {}
