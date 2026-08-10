"""Deterministic release-score helpers for the frozen tomato workflow."""

from plantpersulf.workflows.tomato_ranker_v2 import _percentiles


def test_percentiles_rank_descending_with_stable_tie_breaking() -> None:
    assert _percentiles((0.5, 0.8, 0.5)) == (0.5, 1.0, 0.0)
