"""Task 10 — label-permutation significance test.

Shuffling the positive/unlabeled labels while holding the scores fixed gives a
null distribution for any ranking metric. The p-value is the (add-one
smoothed) fraction of permutations that reach the observed metric, so it is
never exactly zero and is bounded in ``[1/(n+1), 1]``.

Pure and deterministic under a fixed ``seed``.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

Scored = list[tuple[float, str]]
MetricFn = Callable[[Scored], float | None]


@dataclass(frozen=True)
class PermutationResult:
    observed: float
    p_value: float
    n_perm: int


def permutation_test(
    scored: Scored,
    metric_fn: MetricFn,
    n_perm: int = 1000,
    seed: int = 0,
) -> PermutationResult:
    """Permutation test for a ranking metric.

    Scores stay fixed; labels are shuffled ``n_perm`` times. ``p_value`` is
    ``(1 + #{permuted >= observed}) / (1 + n_perm)``.
    """
    import random

    if not scored:
        raise ValueError("scored must not be empty")

    scores = [s for s, _ in scored]
    labels = [y for _, y in scored]

    observed = metric_fn(scored)
    observed_val = 0.0 if observed is None else observed

    rng = random.Random(seed)
    at_least = 0
    for _ in range(n_perm):
        permuted_labels = labels[:]
        rng.shuffle(permuted_labels)
        value = metric_fn(list(zip(scores, permuted_labels, strict=True)))
        if value is not None and value >= observed_val:
            at_least += 1

    p_value = (1.0 + at_least) / (1.0 + n_perm)
    return PermutationResult(observed_val, p_value, n_perm)
