"""Unit tests for within-protein ranking metrics (diagnostic-only)."""

from __future__ import annotations

import pytest

from plantpersulf.evaluation.within_protein_ranking import (
    first_hit_rank,
    hit_at_k,
    mutagenesis_burden,
    rank_sites_desc,
    within_protein_metrics,
)


def test_mutagenesis_burden_matches_documented_random_baselines() -> None:
    """E = (n+1)/(k+1) — the closed form used in the methodology design."""
    assert mutagenesis_burden(n=7, k=1) == pytest.approx(4.0)  # SlWRKY6
    assert mutagenesis_burden(n=5, k=1) == pytest.approx(3.0)  # PyMYB10
    assert mutagenesis_burden(n=14, k=2) == pytest.approx(5.0)  # SlBRG3
    assert mutagenesis_burden(n=14, k=1) == pytest.approx(7.5)


def test_mutagenesis_burden_rejects_invalid_inputs() -> None:
    with pytest.raises(ValueError):
        mutagenesis_burden(n=0, k=1)
    with pytest.raises(ValueError):
        mutagenesis_burden(n=5, k=6)


def test_rank_sites_desc_sorts_by_score() -> None:
    assert rank_sites_desc({206: 0.3, 209: 0.9, 212: 0.1}) == [209, 206, 212]


def test_first_hit_rank_is_min_true_position_rank() -> None:
    # ranking [212, 206, 209, 173]; true = {206, 212} → first hit at rank 1
    assert first_hit_rank([212, 206, 209, 173], {206, 212}) == 1
    # true = {209} → rank 3
    assert first_hit_rank([212, 206, 209, 173], {209}) == 3
    # true = {173} → rank 4
    assert first_hit_rank([212, 206, 209, 173], {173}) == 4


def test_hit_at_k() -> None:
    ranking = [212, 206, 209, 173]
    assert hit_at_k(ranking, {206, 212}, k=2) is True
    assert hit_at_k(ranking, {209}, k=2) is False
    assert hit_at_k(ranking, {209}, k=3) is True


def test_within_protein_metrics_perfect_ranking() -> None:
    metrics = within_protein_metrics(
        positions=[173, 206, 209, 212],
        scores={173: 0.1, 206: 0.4, 209: 0.2, 212: 0.3},
        true_positions={206, 212},
    )
    assert metrics["n_cys"] == 4
    assert metrics["top1"] is True
    assert metrics["hit_at_2"] is True
    assert metrics["mrr"] == pytest.approx(0.75)  # (1/1 + 1/2) / 2 over two sites
    assert metrics["first_hit_burden"] == 1
    assert metrics["random_baseline_burden"] == pytest.approx((4 + 1) / (2 + 1))


def test_within_protein_metrics_worst_ranking() -> None:
    metrics = within_protein_metrics(
        positions=[173, 206, 209, 212],
        scores={173: 0.9, 206: 0.2, 209: 0.3, 212: 0.1},
        true_positions={206, 212},
    )
    assert metrics["top1"] is False
    assert metrics["hit_at_2"] is False
    # 206 ranks 3rd; building constructs in order hits it 3rd
    assert metrics["first_hit_burden"] == 3
    assert metrics["mrr"] == pytest.approx((1 / 3 + 1 / 4) / 2)


def test_within_protein_metrics_requires_all_positions_scored() -> None:
    with pytest.raises(ValueError):
        within_protein_metrics(
            positions=[173, 206, 209, 212],
            scores={173: 0.1, 206: 0.2, 209: 0.3},  # 212 missing
            true_positions={206, 212},
        )


def test_within_protein_metrics_requires_true_positions_present() -> None:
    with pytest.raises(ValueError):
        within_protein_metrics(
            positions=[173, 206, 209, 212],
            scores={173: 0.1, 206: 0.2, 209: 0.3, 212: 0.4},
            true_positions={206, 999},
        )
