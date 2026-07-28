"""Task 8 — positive-vs-unlabeled ranking metrics.

Every metric here operates on a list of ``(score, label)`` pairs where
``label`` is exactly ``"positive"`` or ``"unlabeled"``. The forbidden label
``"negative"`` is rejected explicitly (Codex red line: nondetection is never
a negative). Ties are broken by a stable sort on input order, so results are
deterministic given a fixed score list.

These are PU-appropriate metrics: "positive-vs-unlabeled" AP/MRR/Recall@K,
not standard positive-vs-negative classification metrics.
"""

from __future__ import annotations

ALLOWED_LABELS = frozenset({"positive", "unlabeled"})
Scored = list[tuple[float, str]]


def _validate(scored: Scored) -> None:
    labels = {label for _, label in scored}
    if not labels <= ALLOWED_LABELS:
        forbidden = labels - ALLOWED_LABELS
        raise ValueError(
            f"scored pairs contain forbidden label(s): {sorted(forbidden)}"
        )


def _ranked(scored: Scored) -> Scored:
    """Stable sort by descending score; ties keep input order."""
    return sorted(scored, key=lambda pair: -pair[0])


def recall_at_k(scored: Scored, k: int) -> float | None:
    """Fraction of all positives found in the top-k ranked items."""
    _validate(scored)
    total_positive = sum(1 for _, label in scored if label == "positive")
    if total_positive == 0:
        return None
    ranked = _ranked(scored)
    hits = sum(1 for _, label in ranked[:k] if label == "positive")
    return hits / total_positive


def mean_reciprocal_rank(scored: Scored) -> float | None:
    """Reciprocal rank of the first positive in the ranking."""
    _validate(scored)
    ranked = _ranked(scored)
    for rank, (_, label) in enumerate(ranked, start=1):
        if label == "positive":
            return 1.0 / rank
    return None


def average_precision(scored: Scored) -> float | None:
    """Average of precision-at-k evaluated at each positive's rank."""
    _validate(scored)
    ranked = _ranked(scored)
    total_positive = sum(1 for _, label in scored if label == "positive")
    if total_positive == 0:
        return None
    hits = 0
    precision_sum = 0.0
    for rank, (_, label) in enumerate(ranked, start=1):
        if label == "positive":
            hits += 1
            precision_sum += hits / rank
    return precision_sum / total_positive
