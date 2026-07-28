"""RED (Task 10): protein-cluster bootstrap confidence intervals.

Bootstrap must resample whole *protein clusters*, not individual sites, so a
single high-homology cluster cannot manufacture a tight CI (Codex 10 conclusion
gate condition 5). Pure resampling algorithm — deterministic under a fixed seed.

Expected RED: ``plantpersulf.evaluation.bootstrap`` does not exist yet.
"""

from __future__ import annotations

from plantpersulf.evaluation.bootstrap import (  # RED: module missing
    cluster_bootstrap_ci,
)
from plantpersulf.evaluation.metrics import average_precision


def _scored() -> list[tuple[float, str, str]]:
    # (score, label, cluster_id): two clusters, clean separation.
    return [
        (0.9, "positive", "c1"),
        (0.8, "positive", "c1"),
        (0.2, "unlabeled", "c2"),
        (0.1, "unlabeled", "c2"),
        (0.85, "positive", "c3"),
        (0.15, "unlabeled", "c4"),
    ]


def test_ci_is_deterministic_for_fixed_seed() -> None:
    scored = _scored()
    a = cluster_bootstrap_ci(scored, average_precision, n_boot=200, seed=0)
    b = cluster_bootstrap_ci(scored, average_precision, n_boot=200, seed=0)
    assert a.point == b.point
    assert a.lower == b.lower
    assert a.upper == b.upper


def test_ci_brackets_point_estimate() -> None:
    scored = _scored()
    r = cluster_bootstrap_ci(scored, average_precision, n_boot=500, seed=1)
    assert r.lower <= r.point <= r.upper
    assert 0.0 <= r.lower <= 1.0
    assert 0.0 <= r.upper <= 1.0


def test_resamples_clusters_not_rows() -> None:
    """Cluster-level resampling leaves real between-cluster variance: with
    clusters that rank their positives differently, different cluster draws
    give different AP, so the CI has non-zero width."""
    # Imperfectly separated clusters: c_good ranks its positive top, c_bad
    # ranks its positive below an unlabeled decoy. Which clusters are drawn
    # changes the pooled AP.
    scored = [
        (0.95, "positive", "c_good"),
        (0.10, "unlabeled", "c_good"),
        (0.30, "positive", "c_bad"),
        (0.80, "unlabeled", "c_bad"),
        (0.05, "unlabeled", "c_extra"),
    ]
    r = cluster_bootstrap_ci(scored, average_precision, n_boot=500, seed=2)
    assert r.upper > r.lower


def test_rejects_empty_input() -> None:
    import pytest

    with pytest.raises(ValueError):
        cluster_bootstrap_ci([], average_precision, n_boot=10, seed=0)
