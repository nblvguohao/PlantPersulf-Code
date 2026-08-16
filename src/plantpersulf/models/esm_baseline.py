"""Task 8 — ESM embedding + linear head baseline.

A thin wrapper that standardizes frozen ESM-2 per-cysteine embeddings with
the same leakage-safe ``TrainOnlyScaler`` used by the other traditional
baselines, then fits a plain logistic-regression head on top via
``logistic_regression_scores``. No new model-fitting logic is invented here
— this module only wires together already-tested scaffolding.
"""

from __future__ import annotations

from typing import Any


def esm_linear_head_scores(
    train_embeddings: Any,
    train_y: list[str],
    predict_embeddings: Any,
    seed: int,
) -> list[float]:
    """Fit a logistic-regression head on standardized ESM-2 embeddings.

    ``TrainOnlyScaler`` is fit on ``train_embeddings`` only; the frozen
    statistics are then applied to both train and predict embeddings, so
    nothing about the predict distribution can leak into normalization.
    """
    import numpy as np
    from sklearn.linear_model import LogisticRegression  # type: ignore[import-untyped]

    train = np.asarray(train_embeddings, dtype=np.float32)
    predict = np.asarray(predict_embeddings, dtype=np.float32)
    if train.ndim != 2 or train.shape[0] == 0:
        raise ValueError("train_embeddings must not be empty")
    if predict.ndim != 2 or predict.shape[1] != train.shape[1]:
        predict_dim = predict.shape[1] if predict.ndim == 2 else "invalid"
        raise ValueError(
            "predict_embeddings dimension "
            f"{predict_dim} does not match train_embeddings dimension "
            f"{train.shape[1]}"
        )
    if train.shape[0] != len(train_y):
        raise ValueError("train embeddings and labels must have equal length")
    mean = train.mean(axis=0, dtype=np.float32)
    std = train.std(axis=0, dtype=np.float32)
    std = np.where(std > 0, std, np.float32(1.0))
    scaled_train = (train - mean) / std
    scaled_predict = (predict - mean) / std
    y_binary = np.asarray(
        [1 if label == "positive" else 0 for label in train_y], dtype=np.int8
    )
    model = LogisticRegression(random_state=seed, max_iter=1000)
    model.fit(scaled_train, y_binary)
    return [float(value) for value in model.predict_proba(scaled_predict)[:, 1]]
