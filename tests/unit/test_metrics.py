"""RED (Task 8): PU-ranking evaluation metrics.

Recall@K, Mean Reciprocal Rank, and Average Precision for a
positive-vs-unlabeled ranking. All metrics operate on (score, label)
pairs sorted by descending score; ties are broken deterministically by a
stable sort on the input order.

Expected RED: ``plantpersulf.evaluation.metrics`` does not exist yet.
"""

from __future__ import annotations

import pytest

from plantpersulf.evaluation.metrics import (  # RED: module missing
    average_precision,
    mean_reciprocal_rank,
    recall_at_k,
)


def test_recall_at_k_perfect_ranking() -> None:
    # 2 positives, both ranked first and second.
    scored = [(0.9, "positive"), (0.8, "positive"), (0.1, "unlabeled")]
    assert recall_at_k(scored, k=1) == pytest.approx(0.5)
    assert recall_at_k(scored, k=2) == pytest.approx(1.0)
    assert recall_at_k(scored, k=3) == pytest.approx(1.0)


def test_recall_at_k_worst_ranking() -> None:
    scored = [(0.9, "unlabeled"), (0.8, "unlabeled"), (0.1, "positive")]
    assert recall_at_k(scored, k=1) == pytest.approx(0.0)
    assert recall_at_k(scored, k=3) == pytest.approx(1.0)


def test_recall_at_k_no_positives_returns_none() -> None:
    scored = [(0.9, "unlabeled"), (0.1, "unlabeled")]
    assert recall_at_k(scored, k=1) is None


def test_mean_reciprocal_rank_first_positive_rank() -> None:
    # First positive at rank 2 -> MRR contribution 1/2 for this query.
    scored = [(0.9, "unlabeled"), (0.8, "positive"), (0.1, "positive")]
    assert mean_reciprocal_rank(scored) == pytest.approx(0.5)


def test_mean_reciprocal_rank_no_positives_returns_none() -> None:
    scored = [(0.9, "unlabeled")]
    assert mean_reciprocal_rank(scored) is None


def test_average_precision_perfect_ranking_is_one() -> None:
    scored = [(0.9, "positive"), (0.8, "positive"), (0.1, "unlabeled")]
    assert average_precision(scored) == pytest.approx(1.0)


def test_average_precision_interleaved_ranking() -> None:
    # positive at rank1 (p=1/1), unlabeled at rank2, positive at rank3 (p=2/3)
    scored = [(0.9, "positive"), (0.8, "unlabeled"), (0.7, "positive")]
    expected = (1 / 1 + 2 / 3) / 2
    assert average_precision(scored) == pytest.approx(expected)


def test_average_precision_no_positives_returns_none() -> None:
    scored = [(0.9, "unlabeled")]
    assert average_precision(scored) is None


def test_metrics_reject_forbidden_negative_label() -> None:
    scored = [(0.9, "negative"), (0.1, "positive")]
    with pytest.raises(ValueError, match="negative"):
        recall_at_k(scored, k=1)
