"""Systematic in-regime structure separation for known persulfidation controls.

The P2 early gate (``evaluate_p2_structure_separation.py``) showed on two
proteins that structure features separate co-peptide Cys *within a regime*:
BRG3's RING domain signed composite ranks the true site C206 first and the
gold-standard negative C209 last, while the same signs applied protein-wide
collapse because N-terminal disordered Cys have the same "exposed" shape.

This module is the pure-function layer behind the systematic extension of
that result to all 12 registered mapped controls (5 tomato + 7 Arabidopsis,
each with a registered AlphaFold DB model). The central claims tested:

- **Exposure signal**: do structure features rank the published true site
  above its own protein's other Cys, aggregated across all 12 controls?
  ``aggregate_site_z`` (mean within-protein z, above-median fraction) and
  ``ranking_stats`` (Top-1 / Hit@2 / first-hit burden vs random) answer this
  per feature, without fitting anything.
- **Regime routing (proposition 3)**: the signal is expected to differ by
  structural regime. ``regime_bucket`` / ``site_plddt_regime`` partition Cys
  by the structure's own pLDDT (an orthogonal prior — AF disorder confidence,
  never the site labels), and ``subset_positions`` / ``composite_scores``
  support ranking within a regime scope rather than the whole protein.

``composite_scores`` applies *signed* z-scores. Signs come from the
aggregate direction (a ceiling: what perfect direction learning at n=12
could do), so composite rankings are reported as an upper bound, not as
unbiased evidence — the per-feature aggregation is the evidence.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from typing import Any

import numpy as np

from plantpersulf.evaluation.within_protein_ranking import rank_sites_desc

PLDDT_FOLDED = 70.0
PLDDT_IDR = 50.0


def regime_bucket(plddt: float) -> str:
    """Structural regime of a site from its own AF pLDDT.

    ``folded`` >= 70 (confident structure), ``linker`` 50-70 (domain
    boundary / weakly structured), ``disordered`` < 50 (AF's disorder
    confidence band). Structural prior only — never touches site labels.
    """
    if plddt >= PLDDT_FOLDED:
        return "folded"
    if plddt >= PLDDT_IDR:
        return "linker"
    return "disordered"


def site_plddt_regime(
    feature_rows: Mapping[int, Mapping[str, float]], position: int
) -> str:
    """Regime of one Cys position from its ``plddt`` structure feature."""
    return regime_bucket(float(feature_rows[position]["plddt"]))


def _in_regime_bucket(
    bucket: str,
) -> Callable[[int, Mapping[str, float]], bool]:
    """Predicate factory for Cys whose pLDDT falls in ``bucket``.

    A module-level factory (rather than a closure in a loop) so mypy can
    infer the predicate type for ``subset_positions``.
    """

    def _predicate(position: int, row: Mapping[str, float]) -> bool:
        return regime_bucket(float(row["plddt"])) == bucket

    return _predicate


def within_protein_z(values: Mapping[int, float], position: int) -> float:
    """z-score of one site against the mean/std of the supplied distribution."""
    array = np.fromiter(values.values(), dtype=float, count=len(values))
    mean = float(array.mean())
    std = float(array.std()) or 1.0
    return (values[position] - mean) / std


def subset_positions(
    feature_rows: Mapping[int, Mapping[str, float]],
    predicate: Callable[[int, Mapping[str, float]], bool],
) -> list[int]:
    """Positions (1-based, sorted) whose feature row satisfies ``predicate``."""
    return sorted(
        position for position, row in feature_rows.items() if predicate(position, row)
    )


def _feature_values(
    feature_rows: Mapping[int, Mapping[str, float]], feature: str
) -> dict[int, float]:
    return {position: float(row[feature]) for position, row in feature_rows.items()}


def aggregate_site_z(
    records: Sequence[tuple[Mapping[int, Mapping[str, float]], int]],
    feature: str,
) -> dict[str, Any]:
    """Mean within-protein z of the true site across proteins, one feature.

    Each record is ``(feature_rows, true_position)`` using the registered
    k=1 representative site. Unbiased evidence level: computed per protein,
    then averaged — no sign or fit is learned from the aggregate.
    """
    zs = [
        within_protein_z(_feature_values(rows, feature), true_position)
        for rows, true_position in records
    ]
    return {
        "n": len(zs),
        "mean_z": float(np.mean(zs)) if zs else 0.0,
        "above_median": sum(
            1
            for rows, true_position in records
            if _feature_values(rows, feature)[true_position]
            > float(
                np.median(
                    np.fromiter(
                        _feature_values(rows, feature).values(),
                        dtype=float,
                        count=len(rows),
                    )
                )
            )
        ),
        "per_protein_z": [round(z, 4) for z in zs],
    }


def ranking_stats(
    records: Sequence[tuple[Mapping[int, Mapping[str, float]], int]],
    feature: str,
) -> dict[str, Any]:
    """Aggregate within-protein rank of true sites for one feature.

    For each protein, rank all Cys by the feature (descending) and locate the
    true site. Aggregates: Top-1 count, Hit@2 count, total first-hit burden
    and the total random baseline ``(n+1)/(k+1)`` per protein.
    """
    top1 = 0
    hit2 = 0
    total_burden = 0
    total_random = 0.0
    n = 0
    for rows, true_position in records:
        values = _feature_values(rows, feature)
        ranking = sorted(values, key=lambda position: values[position], reverse=True)
        rank = ranking.index(true_position) + 1
        top1 += rank == 1
        hit2 += rank <= 2
        total_burden += rank
        total_random += (len(rows) + 1) / 2.0  # k=1: (n+1)/(1+1)
        n += 1
    return {
        "n": n,
        "top1": top1,
        "hit_at_2": hit2,
        "total_first_hit_burden": total_burden,
        "total_random_burden": round(total_random, 3),
    }


def composite_scores(
    feature_rows: Mapping[int, Mapping[str, float]],
    signs: Mapping[str, float],
    positions: Sequence[int] | None = None,
) -> dict[int, float]:
    """Signed z-scored composite for each position (upper-bound ranking).

    ``signs`` maps each feature name to +1/-1. z-scores are computed within
    ``positions`` when given (regime-local scope), else within all rows.
    Only ``positions`` (or all positions) are returned, so a regime-local
    ranking ranks strictly within its regime scope. Because the signs are an
    input learned upstream, consumers must label the result a ceiling, not
    unbiased evidence.
    """
    scope = list(positions) if positions is not None else list(feature_rows)
    z_by_feature: dict[str, dict[int, float]] = {}
    for feature, sign in signs.items():
        values = {
            position: float(feature_rows[position][feature]) for position in scope
        }
        z_by_feature[feature] = {
            position: sign * within_protein_z(values, position) for position in scope
        }
    return {
        position: float(np.mean([z_by_feature[feature][position] for feature in signs]))
        for position in scope
    }


# --- leave-one-out composite + permutation null -----------------------------


def signs_from_aggregate(
    records: Sequence[tuple[Mapping[int, Mapping[str, float]], int]],
    features: Sequence[str],
) -> dict[str, float]:
    """Per-feature signs from the aggregate direction over training records.

    Sign = +1 when the mean within-protein z of the true site is >= 0 (the
    feature favours high values at true sites), else -1. Deterministic; used
    for LOO sign learning so the signs never see the held-out protein.
    """
    return {
        feature: (1.0 if aggregate_site_z(records, feature)["mean_z"] >= 0 else -1.0)
        for feature in features
    }


def loo_composite_burden(
    records: Sequence[tuple[Mapping[int, Mapping[str, float]], int]],
    features: Sequence[str],
    scope_mode: str,
    group_key: Callable[[int], object] | None = None,
) -> dict[str, Any]:
    """Leave-one-protein-out signed-composite ranking (unbiased estimate).

    For each held-out protein, signs are learned from the aggregate
    directions of the OTHER proteins only — optionally restricted to the
    held-out protein's ``group_key`` group (regime-grouped routing). The
    composite then ranks the held-out true site within its scope:

    - ``scope_mode="protein"``: rank among ALL Cys of the protein;
    - ``scope_mode="regime"``: rank among Cys in the true site's own pLDDT
      bucket (regime-local routing).

    This is the honest analogue of the in-sample signed-composite ceiling:
    signs never see the held-out protein, so the burden is a (nearly)
    unbiased estimate rather than an upper bound.
    """
    if scope_mode not in ("protein", "regime"):
        raise ValueError(f"unknown scope_mode: {scope_mode!r}")
    per_protein: dict[int, dict[str, object]] = {}
    total = 0
    top1 = 0
    for index, (rows, true_position) in enumerate(records):
        if group_key is None:
            training = [record for j, record in enumerate(records) if j != index]
        else:
            group = group_key(index)
            training = [
                record
                for j, record in enumerate(records)
                if j != index and group_key(j) == group
            ]
            if not training:  # single-member group: fall back to all others
                training = [record for j, record in enumerate(records) if j != index]
        signs = signs_from_aggregate(training, features)
        if scope_mode == "protein":
            scope: Sequence[int] | None = None
        else:
            bucket = site_plddt_regime(rows, true_position)
            scope = subset_positions(rows, _in_regime_bucket(bucket))
        scores = composite_scores(rows, signs, positions=scope)
        ranking = rank_sites_desc(scores)
        rank = ranking.index(true_position) + 1
        total += rank
        top1 += rank == 1
        per_protein[index] = {"rank": rank, "n_in_scope": len(scores)}
    return {
        "n": len(records),
        "top1": top1,
        "total_first_hit_burden": total,
        "total_random_burden": round(
            sum((len(record[0]) + 1) / 2.0 for record in records), 3
        ),
        "per_protein": per_protein,
    }


def permutation_null(
    records: Sequence[tuple[Mapping[int, Mapping[str, float]], int]],
    features: Sequence[str],
    scope_mode: str,
    n_perm: int,
    rng: np.random.RandomState,
) -> list[int]:
    """Permutation null for the LOO composite burden.

    Each permutation reassigns the true site to a uniform-random Cys of each
    protein (exchangeability null), re-learns LOO signs, and records the
    total first-hit burden. Returns the sorted null distribution; the
    observed burden's left-tail position is its permutation p-value.
    """
    null: list[int] = []
    for _ in range(n_perm):
        perm_records = [
            (rows, int(rng.choice(list(rows)))) for rows, _true in records
        ]
        burden = loo_composite_burden(perm_records, features, scope_mode)[
            "total_first_hit_burden"
        ]
        null.append(int(burden))
    return sorted(null)
