"""Unit tests for the paired effect-size statistics used by the Gate-2 wiring.

These statistics answer three distinct questions the conclusion gate asks:

* ``paired_delta_ci`` — is a run-level paired effect (e.g. structure-ablation
  AP deltas across fold x seed runs) distinguishable from zero?
* ``paired_cluster_bootstrap_delta_ci`` — is a per-site effect (release model
  vs baseline scored on the *same* rows) distinguishable from zero when
  resampling whole protein clusters, not rows?
* ``top_cluster_dominance`` — does the metric collapse when the single largest
  cluster is removed (i.e. is the result driven by one homology cluster)?
"""

from __future__ import annotations

import pytest

from plantpersulf.evaluation.effect_size import (
    paired_cluster_bootstrap_delta_ci,
    paired_delta_ci,
    top_cluster_dominance,
)


def _clustered(
    rows: list[tuple[float, str, str]],
) -> list[tuple[float, str, str]]:
    return rows


# ---------------------------------------------------------------------------
# paired_delta_ci (run-level)
# ---------------------------------------------------------------------------


def test_paired_delta_ci_all_positive_excludes_zero() -> None:
    deltas = [0.015, 0.028, 0.041, 0.022, 0.033, 0.019, 0.025, 0.031, 0.017, 0.024]
    result = paired_delta_ci(deltas, n_boot=2000, seed=0)
    assert result.point == pytest.approx(sum(deltas) / len(deltas))
    assert result.lower > 0.0
    assert result.upper > result.lower
    assert result.n_boot == 2000


def test_paired_delta_ci_mixed_sign_spans_zero() -> None:
    deltas = [-0.03, 0.02, -0.01, 0.04, -0.02, 0.01, -0.04, 0.03, -0.02, 0.02]
    result = paired_delta_ci(deltas, n_boot=2000, seed=0)
    assert result.lower < 0.0 < result.upper


def test_paired_delta_ci_is_deterministic() -> None:
    deltas = [0.01, -0.02, 0.03, 0.005, -0.01]
    a = paired_delta_ci(deltas, n_boot=500, seed=7)
    b = paired_delta_ci(deltas, n_boot=500, seed=7)
    assert (a.lower, a.upper) == (b.lower, b.upper)


def test_paired_delta_ci_empty_raises() -> None:
    with pytest.raises(ValueError, match="deltas"):
        paired_delta_ci([])


# ---------------------------------------------------------------------------
# paired_cluster_bootstrap_delta_ci (per-site, cluster-level)
# ---------------------------------------------------------------------------


def _separable_rows() -> tuple[
    list[tuple[float, str, str]], list[tuple[float, str, str]]
]:
    """8 clusters x 6 rows; 1 positive per cluster. The 'model' ranks the
    positive first in every cluster; the 'baseline' is near-flat."""
    model: list[tuple[float, str, str]] = []
    baseline: list[tuple[float, str, str]] = []
    for c in range(8):
        cluster = f"c{c}"
        model.append((0.9 + 0.01 * c, "positive", cluster))
        baseline.append((0.5, "positive", cluster))
        for u in range(5):
            model.append((0.05 + 0.01 * u, "unlabeled", cluster))
            baseline.append((0.5 + 0.001 * (u % 2), "unlabeled", cluster))
    return model, baseline


def test_paired_cluster_delta_better_model_excludes_zero() -> None:
    model, baseline = _separable_rows()
    result = paired_cluster_bootstrap_delta_ci(model, baseline, n_boot=500, seed=0)
    assert result.point > 0.5
    assert result.lower > 0.0


def test_paired_cluster_delta_identical_scores_is_zero() -> None:
    model, _ = _separable_rows()
    result = paired_cluster_bootstrap_delta_ci(model, model, n_boot=200, seed=0)
    assert result.point == pytest.approx(0.0)
    assert result.lower == pytest.approx(0.0)
    assert result.upper == pytest.approx(0.0)


def test_paired_cluster_delta_requires_matched_rows() -> None:
    model, baseline = _separable_rows()
    with pytest.raises(ValueError, match="same rows"):
        paired_cluster_bootstrap_delta_ci(model, baseline[:-1], n_boot=10)
    # labels must match row-by-row, not just in length
    mismatched = [(s, "unlabeled" if y == "positive" else y, c) for s, y, c in baseline]
    with pytest.raises(ValueError, match="same rows"):
        paired_cluster_bootstrap_delta_ci(model, mismatched, n_boot=10)


# ---------------------------------------------------------------------------
# top_cluster_dominance
# ---------------------------------------------------------------------------


def test_top_cluster_dominance_driven_when_ap_concentrated() -> None:
    """Every positive lives in one big cluster that the model ranks on top;
    all other clusters are unlabeled-only. Removing the big cluster leaves no
    positives at all -> the result is entirely cluster-driven."""
    rows: list[tuple[float, str, str]] = []
    for i in range(10):
        rows.append((0.9 - 0.01 * i, "positive", "big"))
    for i in range(40):
        rows.append((0.1 + 0.001 * i, "unlabeled", "big"))
    for c in range(5):
        for u in range(20):
            rows.append((0.2 + 0.001 * u, "unlabeled", f"small{c}"))
    result = top_cluster_dominance(rows)
    assert result.top_cluster_id == "big"
    assert result.driven is True
    assert result.retention_ratio < 0.5


def test_top_cluster_dominance_not_driven_when_spread() -> None:
    """Positives spread evenly across 10 equal clusters, each ranked first
    within its cluster: removing any single cluster barely dents AP."""
    rows: list[tuple[float, str, str]] = []
    for c in range(10):
        cluster = f"c{c}"
        rows.append((0.95, "positive", cluster))
        for u in range(9):
            rows.append((0.05 + 0.01 * u, "unlabeled", cluster))
    result = top_cluster_dominance(rows)
    assert result.driven is False
    assert result.retention_ratio >= 0.5


def test_top_cluster_dominance_no_positives_is_fail_closed() -> None:
    """No positives at all: there is no demonstrated effect to attribute, so
    the honest answer is 'robustness not demonstrated' (driven=True)."""
    rows = [(0.01 * i, "unlabeled", f"c{i % 4}") for i in range(40)]
    result = top_cluster_dominance(rows)
    assert result.driven is True


def test_top_cluster_dominance_empty_raises() -> None:
    with pytest.raises(ValueError, match="scored"):
        top_cluster_dominance([])
