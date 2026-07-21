"""Task 8 — PU learning label-frequency estimation and example weighting.

Implements the Elkan & Noto (2008) "Learning classifiers from only positive
and unlabeled data" weighting scheme. A non-traditional classifier is first
trained to distinguish labeled positives (s=1) from unlabeled examples (as
if unlabeled were negative); its output f(x) approximates P(s=1|x), which is
a *scaled* version of the true posterior P(y=1|x): f(x) = c * P(y=1|x),
where c = P(s=1|y=1) is the label frequency, estimated as the mean f(x) on
held-out labeled positives (never on the same positives used to fit f).

This module only computes weights; it does not fit any classifier. Callers
are responsible for the train/held-out split used to estimate ``c`` and for
never estimating it from validation or test labels.
"""

from __future__ import annotations


def estimate_label_frequency(holdout_positive_scores: list[float]) -> float:
    """Estimate c = P(s=1|y=1) from a non-traditional classifier's scores
    on held-out labeled positives (Elkan-Noto)."""
    if not holdout_positive_scores:
        raise ValueError("holdout_positive_scores must not be empty")
    if any(not (0.0 <= score <= 1.0) for score in holdout_positive_scores):
        raise ValueError("holdout_positive_scores must be in range [0, 1]")
    return sum(holdout_positive_scores) / len(holdout_positive_scores)


def pu_example_weights(
    scores: list[float],
    is_labeled_positive: list[bool],
    label_frequency: float,
) -> list[float]:
    """Compute Elkan-Noto positive-class weights for every example.

    A labeled positive always gets weight 1.0 (it is certainly positive).
    An unlabeled example gets weight (1-c)/c * f(x)/(1-f(x)), clipped to
    [0, 1], the estimated probability it is a hidden true positive.
    """
    if len(scores) != len(is_labeled_positive):
        raise ValueError("scores and is_labeled_positive must have equal length")
    if not (0.0 < label_frequency <= 1.0):
        raise ValueError("label_frequency must be in (0, 1]")

    weights: list[float] = []
    for score, labeled in zip(scores, is_labeled_positive, strict=True):
        if labeled:
            weights.append(1.0)
            continue
        denom = 1.0 - score
        if denom <= 0.0:
            weights.append(1.0)
            continue
        raw = (1.0 - label_frequency) / label_frequency * score / denom
        weights.append(min(1.0, max(0.0, raw)))
    return weights
