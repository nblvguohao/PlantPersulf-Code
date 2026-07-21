"""Task 10 — protein-cluster bootstrap confidence intervals.

Resampling is done over whole **protein clusters**, not individual sites. A
persulfidation benchmark is dominated by a handful of large homology clusters;
resampling rows would let one cluster's members appear many times and
manufacture an artificially tight interval. Cluster-level resampling keeps the
between-cluster variance that the Codex conclusion gate (condition 5: "results
are not driven by one high-homology cluster") requires the CI to reflect.

Pure and deterministic: a fixed ``seed`` fully determines the interval.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

Scored = list[tuple[float, str]]
ClusteredScored = list[tuple[float, str, str]]
MetricFn = Callable[[Scored], float | None]


@dataclass(frozen=True)
class BootstrapResult:
    point: float
    lower: float
    upper: float
    n_boot: int


def _percentile(sorted_values: list[float], q: float) -> float:
    """Linear-interpolation percentile (``q`` in [0, 1]) of a sorted list."""
    if not sorted_values:
        raise ValueError("cannot take a percentile of an empty sequence")
    if len(sorted_values) == 1:
        return sorted_values[0]
    pos = q * (len(sorted_values) - 1)
    lo = int(pos)
    hi = min(lo + 1, len(sorted_values) - 1)
    frac = pos - lo
    return sorted_values[lo] * (1.0 - frac) + sorted_values[hi] * frac


def cluster_bootstrap_ci(
    scored: ClusteredScored,
    metric_fn: MetricFn,
    n_boot: int = 1000,
    alpha: float = 0.05,
    seed: int = 0,
) -> BootstrapResult:
    """Bootstrap a metric's CI by resampling clusters with replacement.

    ``scored`` rows are ``(score, label, cluster_id)``. The point estimate is
    the metric on the full data; each bootstrap replicate draws ``n_clusters``
    clusters with replacement and pools their rows. Replicates whose metric is
    undefined (e.g. no positive in the draw) are skipped.
    """
    import random

    if not scored:
        raise ValueError("scored must not be empty")

    by_cluster: dict[str, Scored] = {}
    for score, label, cluster in scored:
        by_cluster.setdefault(cluster, []).append((score, label))
    clusters = sorted(by_cluster)

    point = metric_fn([(s, y) for s, y, _ in scored])
    point_val = 0.0 if point is None else point

    rng = random.Random(seed)
    replicates: list[float] = []
    for _ in range(n_boot):
        drawn = [clusters[rng.randrange(len(clusters))] for _ in clusters]
        pooled: Scored = []
        for c in drawn:
            pooled.extend(by_cluster[c])
        value = metric_fn(pooled)
        if value is not None:
            replicates.append(value)

    if not replicates:
        return BootstrapResult(point_val, point_val, point_val, n_boot)

    replicates.sort()
    lower = _percentile(replicates, alpha / 2.0)
    upper = _percentile(replicates, 1.0 - alpha / 2.0)
    return BootstrapResult(point_val, lower, upper, n_boot)
