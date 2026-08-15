"""Unit tests for the systematic in-regime structure separation helpers.

Covers ``evaluation.structure_regime`` — the pure-function layer behind the
known-control structure diagnostic (track ``structure_regime_separation_v1``).
Regime routing (proposition 3 of the methodology design) is tested at the
helper level here; the end-to-end 12-protein aggregation lives in
``scripts/evaluate_structure_regime_separation.py``.
"""

from __future__ import annotations

import numpy as np
import pytest

from plantpersulf.evaluation.structure_regime import (
    PLDDT_FOLDED,
    PLDDT_IDR,
    aggregate_site_z,
    composite_scores,
    loo_composite_burden,
    permutation_null,
    ranking_stats,
    regime_bucket,
    signs_from_aggregate,
    site_plddt_regime,
    subset_positions,
    within_protein_z,
)

# --- regime bucketing -------------------------------------------------------


def test_regime_bucket_thresholds() -> None:
    assert regime_bucket(90.0) == "folded"
    assert regime_bucket(PLDDT_FOLDED) == "folded"
    assert regime_bucket(PLDDT_FOLDED - 0.1) == "linker"
    assert regime_bucket(PLDDT_IDR) == "linker"
    assert regime_bucket(PLDDT_IDR - 0.1) == "disordered"
    assert regime_bucket(20.0) == "disordered"


def test_site_plddt_regime_reads_plddt_feature() -> None:
    features = {
        10: {"plddt": 92.0},
        20: {"plddt": 45.0},
        30: {"plddt": 61.0},
    }
    assert site_plddt_regime(features, 10) == "folded"
    assert site_plddt_regime(features, 20) == "disordered"
    assert site_plddt_regime(features, 30) == "linker"


# --- within-protein z -------------------------------------------------------


def test_within_protein_z_standardizes() -> None:
    # population std of {1,2,3} = sqrt(2/3)
    values = {1: 1.0, 2: 2.0, 3: 3.0}
    z1 = within_protein_z(values, 1)
    z3 = within_protein_z(values, 3)
    assert z1 == pytest.approx(-1.0 / (2.0 / 3.0) ** 0.5)
    assert z3 == pytest.approx(1.0 / (2.0 / 3.0) ** 0.5)


def test_within_protein_z_zero_std_guard() -> None:
    values = {1: 5.0, 2: 5.0, 3: 5.0}
    assert within_protein_z(values, 2) == 0.0


# --- subset_positions -------------------------------------------------------


def test_subset_positions_filters_by_predicate() -> None:
    feature_rows = {
        10: {"plddt": 92.0, "rsa_relative": 0.5},
        20: {"plddt": 45.0, "rsa_relative": 0.9},
        30: {"plddt": 88.0, "rsa_relative": 0.1},
    }
    folded = subset_positions(
        feature_rows, lambda p, row: regime_bucket(row["plddt"]) == "folded"
    )
    assert folded == [10, 30]


# --- aggregate_site_z -------------------------------------------------------


def _rows(*values: tuple[int, float]) -> dict[int, dict[str, float]]:
    """Nested feature_rows of one 'feature' for the positions given."""
    return {position: {"feature": value} for position, value in values}


def test_aggregate_site_z_mean_and_median_direction() -> None:
    # protein A: true site is the max of its distribution
    # protein B: true site is the min
    records = [
        (_rows((1, 1.0), (2, 2.0), (3, 3.0)), 3),  # z = +1.22
        (_rows((1, 1.0), (2, 2.0), (3, 3.0)), 1),  # z = -1.22
    ]
    agg = aggregate_site_z(records, "feature")
    assert agg["n"] == 2
    assert agg["mean_z"] == pytest.approx(0.0)
    assert agg["above_median"] == 1
    assert agg["per_protein_z"] == pytest.approx(
        [1.0 / (2.0 / 3.0) ** 0.5, -1.0 / (2.0 / 3.0) ** 0.5],
        rel=1e-3,  # per_protein_z is rounded to 4 decimals for JSON output
    )


# --- ranking_stats ----------------------------------------------------------


def test_ranking_stats_aggregates_burden_and_hits() -> None:
    records = [
        (_rows((1, 3.0), (2, 2.0), (3, 1.0)), 1),  # rank 1/3, top1 hit
        (_rows((1, 3.0), (2, 2.0), (3, 1.0)), 3),  # rank 3/3
    ]
    stats = ranking_stats(records, "feature")
    assert stats["n"] == 2
    assert stats["top1"] == 1
    assert stats["hit_at_2"] == 1
    assert stats["total_first_hit_burden"] == 4
    assert stats["total_random_burden"] == pytest.approx(2 * 2.0)  # (3+1)/(1+1)=2 each


def test_ranking_stats_requires_true_position_present() -> None:
    with pytest.raises(ValueError):
        ranking_stats([(_rows((1, 1.0), (2, 2.0)), 5)], "feature")


# --- composite_scores -------------------------------------------------------


def test_composite_scores_sign_flip() -> None:
    feature_rows = {
        1: {"a": 1.0, "b": 2.0},
        2: {"a": 2.0, "b": 4.0},
        3: {"a": 3.0, "b": 6.0},
    }
    # positively correlated features, both with the "high = true" direction
    std3 = (2.0 / 3.0) ** 0.5
    scores = composite_scores(feature_rows, {"a": 1.0, "b": 1.0})
    assert scores[1] == pytest.approx(-1.0 / std3)
    assert scores[2] == pytest.approx(0.0)
    assert scores[3] == pytest.approx(1.0 / std3)
    # flipping both signs inverts the ranking exactly (the sign semantics)
    flipped = composite_scores(feature_rows, {"a": -1.0, "b": -1.0})
    assert flipped[1] == pytest.approx(1.0 / std3)
    assert flipped[3] == pytest.approx(-1.0 / std3)


def test_composite_scores_restricted_scope() -> None:
    feature_rows = {
        1: {"a": 1.0},
        2: {"a": 2.0},
        3: {"a": 3.0},
        4: {"a": 100.0},
    }
    # z computed only within the given scope; position 4 is both out-of-scope
    # for normalization and absent from the returned scores (regime-local rank)
    scores = composite_scores(feature_rows, {"a": 1.0}, positions=[1, 2, 3])
    assert list(scores) == [1, 2, 3]
    std3 = (2.0 / 3.0) ** 0.5  # population std of {1,2,3}
    assert scores[1] == pytest.approx(-1.0 / std3)
    assert scores[3] == pytest.approx(1.0 / std3)


# --- LOO composite + permutation null ---------------------------------------


def test_signs_from_aggregate_follows_direction() -> None:
    records = [
        (_rows((1, 1.0), (2, 2.0), (3, 3.0)), 3),  # true = max
        (_rows((1, 1.0), (2, 2.0), (3, 3.0)), 3),  # true = max
    ]
    signs = signs_from_aggregate(records, ["feature"])
    assert signs["feature"] == 1.0
    inverted = [
        (_rows((1, 1.0), (2, 2.0), (3, 3.0)), 1),  # true = min
    ]
    assert signs_from_aggregate(inverted, ["feature"])["feature"] == -1.0


def test_loo_composite_burden_global_signs() -> None:
    # protein A true=max, protein B true=min: LOO sign for A comes from B
    # (-1), so A's true (max) is ranked LAST; for B the sign from A (+1)
    # makes B's true (min) ranked LAST too — the honest cross-validation cost
    # of sign learning on n=2.
    records = [
        (_rows((1, 1.0), (2, 2.0), (3, 3.0)), 3),
        (_rows((1, 1.0), (2, 2.0), (3, 3.0)), 1),
    ]
    res = loo_composite_burden(records, ["feature"], scope_mode="protein")
    assert res["n"] == 2
    assert res["per_protein"][0]["rank"] == 3
    assert res["per_protein"][1]["rank"] == 3
    assert res["total_first_hit_burden"] == 6
    assert res["total_random_burden"] == pytest.approx(2 * 2.0)


def test_loo_composite_burden_regime_grouped_signs() -> None:
    # two groups, 2 proteins each: held-out protein trains only on its group
    records = [
        (_rows((1, 1.0), (2, 2.0), (3, 3.0)), 3),  # g0: true=max
        (_rows((1, 1.0), (2, 2.0), (3, 3.0)), 3),  # g0
        (_rows((1, 1.0), (2, 2.0), (3, 3.0)), 1),  # g1: true=min
        (_rows((1, 1.0), (2, 2.0), (3, 3.0)), 1),  # g1
    ]
    res = loo_composite_burden(
        records, ["feature"], scope_mode="protein", group_key=lambda i: i // 2
    )
    assert res["per_protein"][0]["rank"] == 1  # trains on g0 -> sign +1 -> max=1
    assert res["per_protein"][2]["rank"] == 1  # trains on g1 -> sign -1 -> min=1
    assert res["total_first_hit_burden"] == 4


def test_loo_composite_burden_scope_mode_validation() -> None:
    with pytest.raises(ValueError):
        loo_composite_burden(
            [(_rows((1, 1.0)), 1)], ["feature"], scope_mode="bogus"
        )


def test_permutation_null_deterministic_and_sorted() -> None:
    records = [
        (_rows((1, 1.0), (2, 2.0), (3, 3.0)), 3),
        (_rows((1, 1.0), (2, 2.0), (3, 3.0)), 3),
    ]
    null1 = permutation_null(
        records, ["feature"], "protein", 20, np.random.RandomState(7)
    )
    null2 = permutation_null(
        records, ["feature"], "protein", 20, np.random.RandomState(7)
    )
    assert null1 == null2
    assert null1 == sorted(null1)
    assert len(null1) == 20
    assert all(2 <= value <= 6 for value in null1)
