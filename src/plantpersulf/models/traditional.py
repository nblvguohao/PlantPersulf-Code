"""Task 8 — traditional PU baselines with leakage-safe scaffolding.

Three structural guarantees, each independently testable:

1. ``TrainOnlyScaler`` freezes mean/std at ``.fit()`` time from the given
   rows only; ``.transform()`` never recomputes statistics, so applying it
   to validation/test data cannot leak their distribution into training.
2. ``select_hyperparameters`` has no test parameter in its signature at
   all — it is structurally impossible to select hyperparameters using
   test labels, not just discouraged by convention.
3. ``TestEvaluationGuard`` wraps the test set and raises on a second
   evaluation call, enforcing "test runs exactly once".

Model wrappers (Logistic, PU-Logistic via Elkan-Noto weights, Random
Forest, XGBoost) share this scaffolding and the same frozen split.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass
from typing import Any, TypeVar

from plantpersulf.models.pu_risk import pu_example_weights

T = TypeVar("T")


@dataclass(frozen=True)
class TrainOnlyScaler:
    """Feature standardizer whose statistics are frozen at fit time."""

    mean: tuple[float, ...]
    std: tuple[float, ...]

    @classmethod
    def fit(cls, train_features: list[list[float]]) -> TrainOnlyScaler:
        if not train_features:
            raise ValueError("train_features must not be empty")
        n_features = len(train_features[0])
        means = []
        stds = []
        for j in range(n_features):
            column = [row[j] for row in train_features]
            mean = sum(column) / len(column)
            variance = sum((x - mean) ** 2 for x in column) / len(column)
            std = math.sqrt(variance) if variance > 0 else 1.0
            means.append(mean)
            stds.append(std)
        return cls(mean=tuple(means), std=tuple(stds))

    def transform(self, features: list[list[float]]) -> list[list[float]]:
        return [
            [
                (value - m) / s
                for value, m, s in zip(row, self.mean, self.std, strict=True)
            ]
            for row in features
        ]


def select_hyperparameters(
    candidates: list[dict[str, Any]],
    train_X: list[list[float]],
    train_y: list[str],
    val_X: list[list[float]],
    val_y: list[str],
    score_fn: Callable[
        [dict[str, Any], list[list[float]], list[str], list[list[float]], list[str]],
        float,
    ],
) -> dict[str, Any]:
    """Pick the candidate with the highest validation score.

    No test parameter exists in this signature — hyperparameter selection
    is structurally unable to read test data.
    """
    if not candidates:
        raise ValueError("candidates must not be empty")
    best_candidate = candidates[0]
    best_score = float("-inf")
    for candidate in candidates:
        score = score_fn(candidate, train_X, train_y, val_X, val_y)
        if score > best_score:
            best_score = score
            best_candidate = candidate
    return best_candidate


class TestEvaluationGuard:
    """Wraps a test set and permits exactly one evaluation call."""

    __test__ = False  # not a pytest test class despite the name

    def __init__(self, test_X: list[list[float]], test_y: list[str]) -> None:
        self._test_X = test_X
        self._test_y = test_y
        self._evaluated = False

    def evaluate(
        self, score_fn: Callable[[list[list[float]], list[str]], T]
    ) -> T:
        if self._evaluated:
            raise RuntimeError("test set may be evaluated only once")
        self._evaluated = True
        return score_fn(self._test_X, self._test_y)


def logistic_regression_scores(
    train_X: list[list[float]],
    train_y: list[str],
    predict_X: list[list[float]],
    seed: int,
    sample_weight: list[float] | None = None,
) -> list[float]:
    """Fit sklearn LogisticRegression on train and score predict_X."""
    from sklearn.linear_model import LogisticRegression  # type: ignore[import-untyped]

    y_binary = [1 if label == "positive" else 0 for label in train_y]
    model = LogisticRegression(random_state=seed, max_iter=1000)
    model.fit(train_X, y_binary, sample_weight=sample_weight)
    return [float(p[1]) for p in model.predict_proba(predict_X)]


def pu_logistic_regression_scores(
    train_X: list[list[float]],
    train_y: list[str],
    predict_X: list[list[float]],
    seed: int,
    holdout_fraction: float = 0.2,
) -> list[float]:
    """PU-Logistic via Elkan-Noto: fit a non-traditional classifier on
    (positive vs unlabeled-as-negative), estimate label frequency c from a
    held-out slice of the labeled positives, then reweight and refit."""
    import random

    positive_indices = [i for i, label in enumerate(train_y) if label == "positive"]
    if len(positive_indices) < 2:
        raise ValueError("PU-Logistic requires at least 2 labeled positives")
    rng = random.Random(seed)
    shuffled = positive_indices[:]
    rng.shuffle(shuffled)
    n_holdout = max(1, int(len(shuffled) * holdout_fraction))
    holdout_idx = set(shuffled[:n_holdout])
    fit_idx = [i for i in range(len(train_y)) if i not in holdout_idx]

    fit_X = [train_X[i] for i in fit_idx]
    fit_y = [train_y[i] for i in fit_idx]
    nontraditional_scores = logistic_regression_scores(
        fit_X, fit_y, fit_X, seed=seed
    )

    holdout_X = [train_X[i] for i in holdout_idx]
    holdout_scores = logistic_regression_scores(fit_X, fit_y, holdout_X, seed=seed)

    from plantpersulf.models.pu_risk import estimate_label_frequency

    c = estimate_label_frequency(holdout_scores)
    is_labeled_positive = [fit_y[i] == "positive" for i in range(len(fit_y))]
    weights = pu_example_weights(nontraditional_scores, is_labeled_positive, c)

    return logistic_regression_scores(
        fit_X, fit_y, predict_X, seed=seed, sample_weight=weights
    )


def random_forest_scores(
    train_X: list[list[float]],
    train_y: list[str],
    predict_X: list[list[float]],
    seed: int,
    n_estimators: int = 100,
) -> list[float]:
    from sklearn.ensemble import RandomForestClassifier  # type: ignore[import-untyped]

    y_binary = [1 if label == "positive" else 0 for label in train_y]
    model = RandomForestClassifier(n_estimators=n_estimators, random_state=seed)
    model.fit(train_X, y_binary)
    return [float(p[1]) for p in model.predict_proba(predict_X)]


def xgboost_scores(
    train_X: list[list[float]],
    train_y: list[str],
    predict_X: list[list[float]],
    seed: int,
    n_estimators: int = 100,
) -> list[float]:
    from xgboost import XGBClassifier

    y_binary = [1 if label == "positive" else 0 for label in train_y]
    model = XGBClassifier(
        n_estimators=n_estimators,
        random_state=seed,
        eval_metric="logloss",
    )
    model.fit(train_X, y_binary)
    return [float(p[1]) for p in model.predict_proba(predict_X)]
