"""RED (Task 8): ESM embedding + linear head baseline.

Fits a leakage-safe logistic-regression head on top of frozen ESM-2
per-cysteine embeddings, reusing TrainOnlyScaler and logistic_regression_scores
from plantpersulf.models.traditional. Tested here with small synthetic
numeric embedding vectors (pure algorithm correctness, not a biological
claim) — the real frozen embeddings themselves are tested end-to-end in
tests/scientific/test_esm2_features.py and the integration test in
tests/scientific/test_esm_baseline_integration.py.

Expected RED: ``plantpersulf.models.esm_baseline`` does not exist yet.
"""

from __future__ import annotations

import numpy as np
import pytest

from plantpersulf.models.esm_baseline import (  # RED: module missing
    esm_linear_head_scores,
)


def _toy_embeddings() -> tuple[list[tuple[float, ...]], list[str]]:
    # Small synthetic 3-dim "embeddings" standing in for 1280-dim ESM-2
    # vectors; only used to check the fitting/scaling logic is correct.
    train_embeddings = [
        (1.0, 2.0, 3.0),
        (1.1, 2.1, 3.1),
        (10.0, 20.0, 30.0),
        (10.1, 20.1, 30.1),
    ]
    train_y = ["positive", "positive", "unlabeled", "unlabeled"]
    return train_embeddings, train_y


def test_returns_one_score_per_predict_row_in_unit_range() -> None:
    train_embeddings, train_y = _toy_embeddings()
    predict_embeddings = [(1.05, 2.05, 3.05), (10.05, 20.05, 30.05)]

    scores = esm_linear_head_scores(
        train_embeddings, train_y, predict_embeddings, seed=0
    )

    assert len(scores) == len(predict_embeddings)
    for score in scores:
        assert 0.0 <= score <= 1.0


def test_deterministic_given_fixed_seed() -> None:
    train_embeddings, train_y = _toy_embeddings()
    predict_embeddings = [(1.05, 2.05, 3.05), (10.05, 20.05, 30.05)]

    first = esm_linear_head_scores(
        train_embeddings, train_y, predict_embeddings, seed=42
    )
    second = esm_linear_head_scores(
        train_embeddings, train_y, predict_embeddings, seed=42
    )

    assert first == second


def test_predict_only_inputs_do_not_change_train_statistics() -> None:
    """Regression-style check on the TrainOnlyScaler reuse: scoring wildly
    different predict-only embeddings must not change the scores produced
    for a fixed predict row (i.e. statistics come from train only)."""
    train_embeddings, train_y = _toy_embeddings()
    baseline_predict = [(1.05, 2.05, 3.05)]

    baseline_scores = esm_linear_head_scores(
        train_embeddings, train_y, baseline_predict, seed=0
    )

    extreme_predict = [(1.05, 2.05, 3.05), (-9999.0, 1e6, -1e6)]
    extreme_scores = esm_linear_head_scores(
        train_embeddings, train_y, extreme_predict, seed=0
    )

    assert extreme_scores[0] == pytest.approx(baseline_scores[0])


def test_rejects_empty_train_data() -> None:
    with pytest.raises(ValueError, match="empty"):
        esm_linear_head_scores([], [], [(1.0, 2.0, 3.0)], seed=0)


def test_rejects_mismatched_embedding_dimensions() -> None:
    train_embeddings, train_y = _toy_embeddings()
    mismatched_predict = [(1.0, 2.0)]  # 2-dim, but train is 3-dim

    with pytest.raises(ValueError, match="dimension"):
        esm_linear_head_scores(train_embeddings, train_y, mismatched_predict, seed=0)


def test_accepts_memory_mapped_float32_embeddings(tmp_path) -> None:  # type: ignore[no-untyped-def]
    path = tmp_path / "embeddings.npy"
    np.save(
        path,
        np.asarray(
            [
                [1.0, 2.0, 3.0],
                [1.1, 2.1, 3.1],
                [10.0, 20.0, 30.0],
                [10.1, 20.1, 30.1],
                [1.05, 2.05, 3.05],
            ],
            dtype=np.float32,
        ),
    )
    matrix = np.load(path, mmap_mode="r")

    scores = esm_linear_head_scores(
        matrix[:4],
        ["positive", "positive", "unlabeled", "unlabeled"],
        matrix[4:],
        seed=0,
    )

    assert len(scores) == 1
    assert 0.0 <= scores[0] <= 1.0
