"""Deterministic summaries of repeated out-of-fold rank runs."""

from __future__ import annotations

from dataclasses import dataclass
from statistics import median


@dataclass(frozen=True)
class RankSummary:
    site_key: str
    median_score: float
    median_rank: float
    best_rank: int
    worst_rank: int


def summarize_rank_runs(
    site_keys: tuple[str, ...], run_scores: tuple[tuple[float, ...], ...]
) -> tuple[RankSummary, ...]:
    """Summarize deterministic ranks across matched OOF score runs."""
    if not site_keys or not run_scores:
        raise ValueError("site_keys and run_scores must not be empty")
    if any(len(scores) != len(site_keys) for scores in run_scores):
        raise ValueError("every run must score every site")
    ranks_by_run: list[dict[int, int]] = []
    for scores in run_scores:
        order = sorted(
            range(len(scores)),
            key=lambda index: (-scores[index], site_keys[index]),
        )
        ranks_by_run.append(
            {index: rank for rank, index in enumerate(order, start=1)}
        )
    return tuple(
        RankSummary(
            site_key=key,
            median_score=median(scores[index] for scores in run_scores),
            median_rank=median(ranks[index] for ranks in ranks_by_run),
            best_rank=min(ranks[index] for ranks in ranks_by_run),
            worst_rank=max(ranks[index] for ranks in ranks_by_run),
        )
        for index, key in enumerate(site_keys)
    )


@dataclass(frozen=True)
class ApplicabilityEnvelope:
    """Feature-wise range fitted on train-fold rows only."""

    lower: tuple[float, ...]
    upper: tuple[float, ...]

    @classmethod
    def fit(cls, train_features: list[list[float]]) -> ApplicabilityEnvelope:
        if not train_features:
            raise ValueError("train_features must not be empty")
        columns = list(zip(*train_features, strict=True))
        return cls(
            lower=tuple(min(column) for column in columns),
            upper=tuple(max(column) for column in columns),
        )

    def contains(self, row: list[float]) -> bool:
        return all(
            lower <= value <= upper
            for value, lower, upper in zip(row, self.lower, self.upper, strict=True)
        )
