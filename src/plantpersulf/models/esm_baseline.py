"""Task 8 — ESM embedding + linear head baseline.

A thin wrapper that standardizes frozen ESM-2 per-cysteine embeddings with
the same leakage-safe ``TrainOnlyScaler`` used by the other traditional
baselines, then fits a plain logistic-regression head on top via
``logistic_regression_scores``. No new model-fitting logic is invented here
— this module only wires together already-tested scaffolding.
"""

from __future__ import annotations

from plantpersulf.models.traditional import TrainOnlyScaler, logistic_regression_scores


def esm_linear_head_scores(
    train_embeddings: list[tuple[float, ...]],
    train_y: list[str],
    predict_embeddings: list[tuple[float, ...]],
    seed: int,
) -> list[float]:
    """Fit a logistic-regression head on standardized ESM-2 embeddings.

    ``TrainOnlyScaler`` is fit on ``train_embeddings`` only; the frozen
    statistics are then applied to both train and predict embeddings, so
    nothing about the predict distribution can leak into normalization.
    """
    if not train_embeddings:
        raise ValueError("train_embeddings must not be empty")

    train_dim = len(train_embeddings[0])
    for row in train_embeddings:
        if len(row) != train_dim:
            raise ValueError("all train_embeddings must have the same dimension")
    for row in predict_embeddings:
        if len(row) != train_dim:
            raise ValueError(
                "predict_embeddings dimension "
                f"{len(row)} does not match train_embeddings dimension {train_dim}"
            )

    train_X = [list(row) for row in train_embeddings]
    predict_X = [list(row) for row in predict_embeddings]

    scaler = TrainOnlyScaler.fit(train_X)
    scaled_train_X = scaler.transform(train_X)
    scaled_predict_X = scaler.transform(predict_X)

    return logistic_regression_scores(
        scaled_train_X, train_y, scaled_predict_X, seed=seed
    )
