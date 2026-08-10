"""Conservative paired out-of-fold model-admission decisions."""

from __future__ import annotations

import random
from dataclasses import dataclass


@dataclass(frozen=True)
class AdmissionDecision:
    admitted: bool
    reason: str
    paired_intervals: dict[str, tuple[float, float]]


def _percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    return ordered[int((len(ordered) - 1) * fraction)]


def _paired_block_interval(
    deltas_by_block: list[float], n_boot: int, seed: int
) -> tuple[float, float]:
    if len(deltas_by_block) < 5:
        raise ValueError("at least five repeated-CV blocks are required")
    if n_boot < 500:
        raise ValueError("n_boot must be at least 500")
    rng = random.Random(seed)
    means = [
        sum(rng.choice(deltas_by_block) for _ in deltas_by_block)
        / len(deltas_by_block)
        for _ in range(n_boot)
    ]
    return _percentile(means, 0.025), _percentile(means, 0.975)


def admit_candidate_model(
    candidate_runs: dict[str, tuple[float, float, float]],
    baseline_runs: dict[str, tuple[float, float, float]],
    run_blocks: dict[str, str],
    n_boot: int = 2_000,
    seed: int = 20_260_811,
) -> AdmissionDecision:
    """Admit only if all paired OOF metrics improve conservatively."""
    run_keys = set(candidate_runs)
    if run_keys != set(baseline_runs) or run_keys != set(run_blocks):
        raise ValueError("candidate, baseline, and block keys must match")
    metric_names = ("recall_at_k", "enrichment", "mrr")
    intervals: dict[str, tuple[float, float]] = {}
    for metric_index, metric_name in enumerate(metric_names):
        block_deltas = [
            sum(
                candidate_runs[key][metric_index] - baseline_runs[key][metric_index]
                for key in sorted(run_keys)
                if run_blocks[key] == block
            )
            / sum(run_blocks[key] == block for key in run_keys)
            for block in sorted(set(run_blocks.values()))
        ]
        interval = _paired_block_interval(block_deltas, n_boot, seed + metric_index)
        intervals[metric_name] = interval
        leave_one_block_out = [
            sum(value for index, value in enumerate(block_deltas) if index != held_out)
            / (len(block_deltas) - 1)
            for held_out in range(len(block_deltas))
        ]
        if interval[0] <= 0.0 or min(leave_one_block_out) <= 0.0:
            return AdmissionDecision(
                False,
                f"{metric_name}_conservative_interval_or_stability_failed",
                intervals,
            )
    return AdmissionDecision(
        True,
        "all_required_paired_intervals_exclude_zero",
        intervals,
    )
