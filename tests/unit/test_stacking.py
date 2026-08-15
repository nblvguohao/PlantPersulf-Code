"""Unit tests for deterministic model stacking / rank aggregation helpers."""

from __future__ import annotations

import pytest

from plantpersulf.evaluation.stacking import (
    align_site_scores,
    fit_score_combiner,
    format_site_key,
    rank_aggregate_scores,
)


def test_format_site_key_joins_protein_and_position() -> None:
    assert (
        format_site_key("arabidopsis|A0A068FJK1", 519)
        == "arabidopsis|A0A068FJK1|519"
    )


def test_align_site_scores_sorts_and_aligns() -> None:
    primary = {"b": 0.2, "a": 0.8}
    secondary = {"a": 0.1, "b": 0.9}
    keys, scores_a, scores_b = align_site_scores(primary, secondary)
    assert keys == ("a", "b")
    assert scores_a == (0.8, 0.2)
    assert scores_b == (0.1, 0.9)


def test_align_site_scores_rejects_mismatched_key_sets() -> None:
    with pytest.raises(ValueError):
        align_site_scores({"a": 0.5, "b": 0.4}, {"a": 0.5})


def test_rank_aggregate_scores_orders_by_mean_rank() -> None:
    # Model A ranks best the first site, model B the second; both agree the
    # third is worst.
    keys = ("a", "b", "c")
    scores_a = (0.9, 0.5, 0.1)
    scores_b = (0.2, 0.8, 0.05)
    combined = rank_aggregate_scores(keys, (scores_a, scores_b))
    # a: ranks (1, 2) -> mean 1.5 ; b: (2, 1) -> 1.5 ; c: (3, 3) -> 3
    assert combined == (-1.5, -1.5, -3.0)
    order = sorted(keys, key=lambda k: combined[keys.index(k)], reverse=True)
    assert order[:2] == ["a", "b"]  # ties keep stable input order
    assert order[2] == "c"


def test_rank_aggregate_scores_rejects_misaligned_runs() -> None:
    with pytest.raises(ValueError):
        rank_aggregate_scores(("a", "b"), ((0.5, 0.1), (0.3,)))


def test_fit_score_combiner_prefers_informative_first_model() -> None:
    # Model A is informative, model B is shuffled noise: the optimal convex
    # combination must lean toward A (alpha > 0.5).
    scores_a = (0.9, 0.8, 0.2, 0.1, 0.7, 0.3)
    scores_b = (0.3, 0.7, 0.5, 0.4, 0.2, 0.8)
    labels = ("positive", "positive", "negative", "negative", "positive", "negative")
    alpha = fit_score_combiner(scores_a, scores_b, labels, n_grid=11)
    assert alpha > 0.5


def test_fit_score_combiner_prefers_informative_second_model() -> None:
    scores_a = (0.3, 0.7, 0.5, 0.4, 0.2, 0.8)
    scores_b = (0.9, 0.8, 0.2, 0.1, 0.7, 0.3)
    labels = ("positive", "positive", "negative", "negative", "positive", "negative")
    alpha = fit_score_combiner(scores_a, scores_b, labels, n_grid=11)
    assert alpha < 0.5


def test_fit_score_combiner_validates_alignment() -> None:
    with pytest.raises(ValueError):
        fit_score_combiner((0.5, 0.1), (0.3, 0.4, 0.2), ("positive", "negative"))
    with pytest.raises(ValueError):
        fit_score_combiner((), (), ())
    with pytest.raises(ValueError):
        fit_score_combiner((0.5,), (0.3,), ("positive",), n_grid=1)
