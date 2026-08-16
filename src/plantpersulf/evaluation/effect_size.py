"""Task 10 (P2) — paired effect-size statistics for the Gate-2 wiring.

Three complementary measurements, all deterministic under a fixed ``seed``:

* ``paired_delta_ci`` — a plain paired bootstrap over run-level effect sizes
  (e.g. structure-ablation AP deltas across fold x seed runs). Used for the
  Gate-2 structure-gain condition, where the unit of evidence is one trained
  run.
* ``paired_cluster_bootstrap_delta_ci`` — the per-site effect (release model
  vs baseline scored on the *same* held-out rows) with whole protein clusters
  resampled, never individual rows. Resampling rows would manufacture a tight
  interval from one large homology cluster; cluster-level resampling is what
  makes the interval honest about cluster composition.
* ``top_cluster_dominance`` — does the metric survive removal of the single
  largest cluster? If more than half the metric is attributable to one
  cluster (retention ratio below ``min_retention``), the result is declared
  cluster-driven. Fail-closed: with no measurable effect at all, robustness
  is treated as not demonstrated.
"""

from __future__ import annotations

import random
from dataclasses import dataclass

from plantpersulf.evaluation.bootstrap import BootstrapResult, ClusteredScored
from plantpersulf.evaluation.metrics import average_precision

DEFAULT_MIN_RETENTION = 0.5


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


def _ci(
    point: float, replicates: list[float], n_boot: int, alpha: float
) -> BootstrapResult:
    if not replicates:
        return BootstrapResult(point, point, point, n_boot)
    replicates.sort()
    return BootstrapResult(
        point,
        _percentile(replicates, alpha / 2.0),
        _percentile(replicates, 1.0 - alpha / 2.0),
        n_boot,
    )


def paired_delta_ci(
    deltas: list[float],
    n_boot: int = 2000,
    alpha: float = 0.05,
    seed: int = 0,
) -> BootstrapResult:
    """Bootstrap CI over run-level paired deltas (resample deltas, take means).

    Each delta must already be a paired effect (same fold/seed, two arms);
    the bootstrap asks how stable the *mean effect* is across runs.
    """
    if not deltas:
        raise ValueError("deltas must not be empty")
    point = sum(deltas) / len(deltas)
    rng = random.Random(seed)
    replicates: list[float] = []
    n = len(deltas)
    for _ in range(n_boot):
        draw = [deltas[rng.randrange(n)] for _ in range(n)]
        replicates.append(sum(draw) / n)
    return _ci(point, replicates, n_boot, alpha)


def paired_cluster_bootstrap_delta_ci(
    model_scored: ClusteredScored,
    baseline_scored: ClusteredScored,
    n_boot: int = 1000,
    alpha: float = 0.05,
    seed: int = 0,
) -> BootstrapResult:
    """Cluster-bootstrap CI of the AP delta between two models scored on the
    same rows.

    Both inputs must describe the same sites in the same order — same label
    and same cluster at every index (only the scores may differ). Each
    replicate resamples whole clusters with replacement and evaluates both
    models on the *same* draw, so the pairing cancels site-composition
    variance and isolates the model difference.
    """
    if len(model_scored) != len(baseline_scored):
        raise ValueError(
            "model and baseline must be scored on the same rows "
            f"(got {len(model_scored)} vs {len(baseline_scored)})"
        )
    for (_, y_m, c_m), (_, y_b, c_b) in zip(model_scored, baseline_scored, strict=True):
        if (y_m, c_m) != (y_b, c_b):
            raise ValueError(
                "model and baseline rows diverge — both arms must be scored "
                "on the same rows in the same order"
            )
    if not model_scored:
        raise ValueError("scored must not be empty")

    by_cluster: dict[str, list[int]] = {}
    for i, (_, _, cluster) in enumerate(model_scored):
        by_cluster.setdefault(cluster, []).append(i)
    clusters = sorted(by_cluster)

    def _ap(rows: ClusteredScored, idx: list[int]) -> float | None:
        return average_precision([(rows[i][0], rows[i][1]) for i in idx])

    all_idx = list(range(len(model_scored)))
    ap_m = _ap(model_scored, all_idx)
    ap_b = _ap(baseline_scored, all_idx)
    point = (ap_m or 0.0) - (ap_b or 0.0)

    rng = random.Random(seed)
    replicates: list[float] = []
    for _ in range(n_boot):
        drawn = [clusters[rng.randrange(len(clusters))] for _ in clusters]
        idx = [i for c in drawn for i in by_cluster[c]]
        rep_m = _ap(model_scored, idx)
        rep_b = _ap(baseline_scored, idx)
        if rep_m is None or rep_b is None:
            continue
        replicates.append(rep_m - rep_b)
    return _ci(point, replicates, n_boot, alpha)


@dataclass(frozen=True)
class DominanceResult:
    top_cluster_id: str
    top_cluster_rows: int
    ap_full: float
    ap_without_top: float
    retention_ratio: float
    driven: bool


def top_cluster_dominance(
    scored: ClusteredScored,
    min_retention: float = DEFAULT_MIN_RETENTION,
) -> DominanceResult:
    """Measure how much of the average precision survives removal of the
    single largest cluster.

    ``retention_ratio = AP(without top cluster) / AP(full)``; the result is
    declared cluster-driven when the ratio falls below ``min_retention``
    (default 0.5 — more than half the metric attributable to one cluster).
    Fail-closed: when the full-data AP is zero/undefined (no positives), the
    ratio is 0 and the result is driven=True.
    """
    if not scored:
        raise ValueError("scored must not be empty")

    counts: dict[str, int] = {}
    for _, _, cluster in scored:
        counts[cluster] = counts.get(cluster, 0) + 1
    top_cluster = sorted(counts, key=lambda c: (-counts[c], c))[0]

    ap_full_val = average_precision([(s, y) for s, y, _ in scored])
    ap_full = ap_full_val if ap_full_val is not None else 0.0
    rest = [(s, y) for s, y, c in scored if c != top_cluster]
    ap_rest_val = average_precision(rest)
    ap_rest = ap_rest_val if ap_rest_val is not None else 0.0

    ratio = ap_rest / ap_full if ap_full > 0.0 else 0.0
    return DominanceResult(
        top_cluster_id=top_cluster,
        top_cluster_rows=counts[top_cluster],
        ap_full=ap_full,
        ap_without_top=ap_rest,
        retention_ratio=ratio,
        driven=ratio < min_retention,
    )
