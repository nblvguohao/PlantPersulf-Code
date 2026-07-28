"""RED (Task 10): label-permutation significance test.

Permuting the positive/unlabeled labels and recomputing the metric gives a null
distribution; the p-value is the fraction of permutations reaching the observed
metric. Deterministic under a fixed seed.

Expected RED: ``plantpersulf.evaluation.permutation`` does not exist yet.
"""

from __future__ import annotations

from plantpersulf.evaluation.metrics import average_precision
from plantpersulf.evaluation.permutation import (  # RED: module missing
    permutation_test,
)


def _separated() -> list[tuple[float, str]]:
    return [
        (0.95, "positive"),
        (0.90, "positive"),
        (0.85, "positive"),
        (0.20, "unlabeled"),
        (0.10, "unlabeled"),
        (0.05, "unlabeled"),
    ]


def _random_like() -> list[tuple[float, str]]:
    return [
        (0.5, "positive"),
        (0.5, "unlabeled"),
        (0.5, "positive"),
        (0.5, "unlabeled"),
    ]


def test_deterministic_for_fixed_seed() -> None:
    scored = _separated()
    a = permutation_test(scored, average_precision, n_perm=200, seed=0)
    b = permutation_test(scored, average_precision, n_perm=200, seed=0)
    assert a.observed == b.observed
    assert a.p_value == b.p_value


def test_strong_separation_is_significant() -> None:
    scored = _separated()
    r = permutation_test(scored, average_precision, n_perm=500, seed=0)
    # Perfect ranking should be rare under the null.
    assert r.p_value < 0.05
    assert 0.0 <= r.p_value <= 1.0


def test_p_value_in_unit_range_and_bounded_below() -> None:
    scored = _separated()
    r = permutation_test(scored, average_precision, n_perm=100, seed=3)
    # p is (1 + #>=obs) / (1 + n_perm): never zero, never above one.
    assert r.p_value >= 1.0 / (1.0 + 100)
    assert r.p_value <= 1.0
