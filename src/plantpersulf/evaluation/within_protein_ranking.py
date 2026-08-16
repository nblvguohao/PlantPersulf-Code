"""Within-protein ranking metrics for registered known controls.

Diagnostic-only helpers: for each control protein, the frozen bundle's
scores of ALL its Cys are ranked within the protein, and the true published
site(s) are located in that ranking. This answers the task the four
published workflows actually perform — "which Cys of this already-selected
protein" — in contrast to proteome-wide Top-K percentiles. Metrics follow
the methodology design (Mutagenesis Burden = expected number of C->A
constructs until the first true site is hit, random baseline
``E = (n+1)/(k+1)``).

Nothing here fits, tunes, or mutates frozen artifacts.
"""

from __future__ import annotations

from collections.abc import Collection, Sequence
from typing import Any


def mutagenesis_burden(n: int, k: int) -> float:
    """Random-sorting baseline: expected constructs to first hit, (n+1)/(k+1)."""
    if n < 1 or not 1 <= k <= n:
        raise ValueError(f"invalid (n={n}, k={k}): need 1 <= k <= n")
    return (n + 1) / (k + 1)


def rank_sites_desc(scores: dict[int, float]) -> list[int]:
    """Protein positions sorted by score, descending."""
    return sorted(scores, key=lambda position: scores[position], reverse=True)


def first_hit_rank(ranking: Sequence[int], true_positions: Collection[int]) -> int:
    """1-based rank of the first true site when constructs are built in
    ranking order (the Mutagenesis Burden of a given ranking)."""
    for rank, position in enumerate(ranking, start=1):
        if position in true_positions:
            return rank
    raise ValueError("no true position present in ranking")


def hit_at_k(ranking: Sequence[int], true_positions: Collection[int], k: int) -> bool:
    return any(position in true_positions for position in ranking[:k])


def within_protein_metrics(
    *,
    positions: Sequence[int],
    scores: dict[int, float],
    true_positions: Collection[int],
) -> dict[str, Any]:
    """Ranking metrics for one protein.

    ``positions`` enumerates every Cys of the protein; ``scores`` maps each
    position to its frozen-bundle score; ``true_positions`` are the published
    site(s) of the control lineage.
    """
    missing = [position for position in positions if position not in scores]
    if missing:
        raise ValueError(f"positions without scores: {missing}")
    missing_true = [
        position for position in true_positions if position not in positions
    ]
    if missing_true:
        raise ValueError(f"true positions not in protein Cys set: {missing_true}")

    ranking = rank_sites_desc(scores)
    ranks = {position: rank for rank, position in enumerate(ranking, start=1)}
    n = len(positions)
    k = len(true_positions)

    return {
        "n_cys": n,
        "ranking": ranking,
        "true_site_ranks": {position: ranks[position] for position in true_positions},
        "top1": (
            ranks[min(true_positions, key=lambda position: ranks[position])] == 1
            if k
            else False
        ),
        "hit_at_2": hit_at_k(ranking, true_positions, k=2),
        "hit_at_3": hit_at_k(ranking, true_positions, k=3),
        "mrr": sum(1.0 / ranks[position] for position in true_positions) / k,
        "first_hit_burden": first_hit_rank(ranking, true_positions),
        "random_baseline_burden": mutagenesis_burden(n, k),
    }
