"""Task 9 — score calibration + uncertainty aggregation for the PU ranker.

Two small, independently testable pieces:

1. ``PlattCalibrator`` — a one-parameter-pair logistic recalibration of raw
   ranker scores, fit **on validation only** (never on test). This maps raw
   scores to calibrated probabilities without touching the ranking order's
   fit; statistics come from validation so nothing about test leaks in.

2. ``aggregate_uncertainty`` — combine MC-dropout per-sample scores into a
   point estimate + standard deviation, so callers share one definition of the
   uncertainty head rather than re-deriving it.

No new model-fitting logic beyond a plain logistic on a single feature (the raw
score) is invented here.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class PlattCalibrator:
    """Logistic recalibration: ``p = sigmoid(a * raw_score + b)``."""

    a: float
    b: float

    @classmethod
    def fit(
        cls, scores: list[float], labels: list[str], seed: int = 0
    ) -> PlattCalibrator:
        """Fit ``a, b`` by logistic regression of the binary label on the raw
        score. Intended to be called on a **validation** split only."""
        if not scores:
            raise ValueError("scores must not be empty")
        if len(scores) != len(labels):
            raise ValueError("scores and labels must have equal length")
        y = [1 if label == "positive" else 0 for label in labels]
        if len(set(y)) < 2:
            # Degenerate validation split (all one class): identity calibration.
            return cls(a=1.0, b=0.0)

        from sklearn.linear_model import (  # type: ignore[import-untyped]
            LogisticRegression,
        )

        model = LogisticRegression(random_state=seed, max_iter=1000)
        model.fit([[s] for s in scores], y)
        return cls(a=float(model.coef_[0][0]), b=float(model.intercept_[0]))

    def transform(self, scores: list[float]) -> list[float]:
        import math

        return [1.0 / (1.0 + math.exp(-(self.a * s + self.b))) for s in scores]


def aggregate_uncertainty(
    samples: list[list[float]],
) -> tuple[list[float], list[float]]:
    """Reduce MC-dropout samples ``[pass][row]`` to per-row (mean, std)."""
    if not samples:
        return [], []
    n_pass = len(samples)
    n_rows = len(samples[0])
    for row in samples:
        if len(row) != n_rows:
            raise ValueError("all MC-dropout passes must have equal length")
    means: list[float] = []
    stds: list[float] = []
    for j in range(n_rows):
        col = [samples[i][j] for i in range(n_pass)]
        mean = sum(col) / n_pass
        var = sum((x - mean) ** 2 for x in col) / n_pass
        means.append(mean)
        stds.append(var**0.5)
    return means, stds
