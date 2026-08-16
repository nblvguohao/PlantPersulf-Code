"""RED (Task 8): PU learning risk weighting (Elkan-Noto).

Implements the standard positive-unlabeled weighting scheme (Elkan & Noto,
2008): given a "non-traditional" classifier trained to distinguish labeled
positives from unlabeled examples (as if unlabeled were negative), the
probability it outputs, f(x) = P(s=1|x), is a scaled version of the true
posterior P(y=1|x) = f(x) / c, where c = P(s=1|y=1) is estimated from the
classifier's scores on held-out labeled positives.

Expected RED: ``plantpersulf.models.pu_risk`` does not exist yet.
"""

from __future__ import annotations

import pytest

from plantpersulf.models.pu_risk import (  # RED: module missing
    estimate_label_frequency,
    pu_example_weights,
)


def test_estimate_label_frequency_is_mean_of_holdout_positive_scores() -> None:
    # c = P(s=1 | y=1) estimated as the mean non-traditional score on
    # held-out labeled positives.
    holdout_positive_scores = [0.8, 0.6, 0.7, 0.9]
    c = estimate_label_frequency(holdout_positive_scores)
    assert c == pytest.approx(0.75)


def test_estimate_label_frequency_rejects_empty_input() -> None:
    with pytest.raises(ValueError, match="empty"):
        estimate_label_frequency([])


def test_estimate_label_frequency_rejects_out_of_range_scores() -> None:
    with pytest.raises(ValueError, match="range"):
        estimate_label_frequency([0.5, 1.4])


def test_pu_weights_known_positive_gets_weight_one() -> None:
    # A labeled positive (s=1) always has true-positive weight 1.0.
    weights = pu_example_weights(
        scores=[0.9], is_labeled_positive=[True], label_frequency=0.75
    )
    assert weights[0] == pytest.approx(1.0)


def test_pu_weights_unlabeled_example_uses_elkan_noto_formula() -> None:
    # weight = (1 - c) / c * f(x) / (1 - f(x)), clipped to [0, 1].
    c = 0.75
    f = 0.6
    expected = min(1.0, max(0.0, (1 - c) / c * f / (1 - f)))
    weights = pu_example_weights(
        scores=[f], is_labeled_positive=[False], label_frequency=c
    )
    assert weights[0] == pytest.approx(expected)


def test_pu_weights_are_clipped_to_valid_probability_range() -> None:
    # A very high f(x) with a low c would exceed 1.0 without clipping.
    weights = pu_example_weights(
        scores=[0.99], is_labeled_positive=[False], label_frequency=0.1
    )
    assert 0.0 <= weights[0] <= 1.0


def test_pu_weights_rejects_mismatched_lengths() -> None:
    with pytest.raises(ValueError, match="length"):
        pu_example_weights(
            scores=[0.5, 0.6],
            is_labeled_positive=[True],
            label_frequency=0.5,
        )


def test_pu_weights_rejects_invalid_label_frequency() -> None:
    with pytest.raises(ValueError, match="label_frequency"):
        pu_example_weights(
            scores=[0.5], is_labeled_positive=[False], label_frequency=0.0
        )
