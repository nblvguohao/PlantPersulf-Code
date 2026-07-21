"""RED (Task 8): normalization statistics must come from the train fold only.

TrainOnlyScaler.fit() must be called with train rows only; its mean/std are
frozen at fit time and are never recomputed by .transform(), so applying it
to validation/test features cannot leak their distribution into the
normalization statistics.

Expected RED: ``plantpersulf.models.traditional`` does not exist yet.
"""

from __future__ import annotations

import pytest

from plantpersulf.models.traditional import (  # RED: module missing
    TrainOnlyScaler,
)


def test_scaler_mean_and_std_come_only_from_fit_data() -> None:
    train = [[1.0, 10.0], [2.0, 20.0], [3.0, 30.0]]
    scaler = TrainOnlyScaler.fit(train)

    assert scaler.mean == pytest.approx([2.0, 20.0])


def test_transform_does_not_change_frozen_statistics() -> None:
    train = [[1.0, 10.0], [2.0, 20.0], [3.0, 30.0]]
    scaler = TrainOnlyScaler.fit(train)
    mean_before = scaler.mean
    std_before = scaler.std

    # Transforming wildly different (e.g. test-distribution) data must not
    # alter the frozen train statistics.
    extreme_other_split = [[1000.0, -500.0], [-9999.0, 42.0]]
    scaler.transform(extreme_other_split)

    assert scaler.mean == mean_before
    assert scaler.std == std_before


def test_transform_uses_frozen_train_statistics_not_input_data() -> None:
    train = [[0.0], [10.0]]
    scaler = TrainOnlyScaler.fit(train)  # mean=5, std=5

    transformed = scaler.transform([[5.0]])  # (5 - 5) / 5 = 0.0
    assert transformed[0][0] == pytest.approx(0.0)


def test_scaler_is_immutable_after_fit() -> None:
    train = [[1.0], [2.0], [3.0]]
    scaler = TrainOnlyScaler.fit(train)
    with pytest.raises((AttributeError, TypeError)):
        scaler.mean = (999.0,)  # type: ignore[misc]


def test_fit_rejects_empty_train_data() -> None:
    with pytest.raises(ValueError, match="empty"):
        TrainOnlyScaler.fit([])
