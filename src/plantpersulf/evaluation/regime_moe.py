"""Regime-routed MoE diagnostic helpers (diagnostic A of the v2 design).

The v2 structure branch design (``2026-08-15-regime-routing-moe-design.md``)
proposes replacing the frozen additive gated fusion with *regime routing*:
a gate driven by structural context (pLDDT) routes each site to a regime
expert (disordered / folded / RING-metal-cluster), and each expert scores
its own sites with its own feature direction.

This module is the pure-function layer behind the first diagnostic: **does
routing beat the single contact feature?** The contact_number_10a burial
signal is the only LOO-significant structure feature (§9.4 of the upgrade
diagnostics, p=0.018 protein-wide / p=0.032 regime-local). The routed
composite must at least match it to justify the architecture.

Key integrity property: **expert directions are written down as structural
priors, never learned from labeled data.** The frozen model's structure
branch is untouched; nothing here trains a model. Because the directions are
fixed a priori, there is no leave-one-out sign-learning bias — the
exchangeability permutation null is the honest reference (contrast
``structure_regime.loo_composite_burden``, which re-learns signs from the
other proteins).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

from plantpersulf.evaluation.structure_regime import (
    _in_regime_bucket,
    regime_bucket,
    subset_positions,
    within_protein_z,
)

# ---------------------------------------------------------------------------
# expert direction priors (structural hypotheses, NOT learned)
# ---------------------------------------------------------------------------
#
# Each regime expert assigns +1/-1/0 to each structure feature. The priors:
# - folded (pLDDT >= 70): burial/packing — contact_number is the ONLY
#   LOO-significant feature (§9.4), so it is the sole +1 here; the rest are
#   left at 0 (not equal-weight, exactly to avoid the composite dilution the
#   LOO showed for the 7-feature equal-weight composite).
# - disordered (pLDDT < 50): exposure/isolation — high RSA, far nearest-Sγ
#   (the P2/§9.2 hypothesis; marked a hypothesis, not a proven direction).
# - linker (50 <= pLDDT < 70): no prior — all 0. This is the honest default:
#   we have no directional evidence in the linker band, so it contributes
#   nothing to a routed score and is skipped by ``routed_burden``.

EXPERT_DIRECTIONS: dict[str, dict[str, float]] = {
    "folded": {
        "contact_number_10a": 1.0,
        "rsa_relative": 0.0,
        "cys_count_8a": 0.0,
        "nearest_sg_distance": 0.0,
        "positive_residue_count_6a": 0.0,
        "coulomb_potential_sg": 0.0,
        "plddt": 0.0,
    },
    "linker": {
        "contact_number_10a": 0.0,
        "rsa_relative": 0.0,
        "cys_count_8a": 0.0,
        "nearest_sg_distance": 0.0,
        "positive_residue_count_6a": 0.0,
        "coulomb_potential_sg": 0.0,
        "plddt": 0.0,
    },
    "disordered": {
        "contact_number_10a": 0.0,
        "rsa_relative": 1.0,
        "cys_count_8a": 0.0,
        "nearest_sg_distance": 1.0,
        "positive_residue_count_6a": 0.0,
        "coulomb_potential_sg": 0.0,
        "plddt": 0.0,
    },
}


def moe_route_scores(
    feature_rows: Mapping[int, Mapping[str, float]],
    directions: Mapping[str, Mapping[str, float]],
    positions: Sequence[int] | None = None,
) -> dict[int, float]:
    """Route each position to its pLDDT regime expert and score it.

    Each site's score is the mean of ``sign * within_regime_z`` over the
    features its regime expert activates (sign = expert direction). z-scores
    are computed within the site's OWN regime bucket only (regime-local scope,
    ``scope_mode="regime"``), so a disordered site is never compared to folded
    sites. A site whose regime expert activates no features (linker: all 0)
    gets score 0.0 — it contributes nothing to a routed ranking.
    """
    scope = list(positions) if positions is not None else list(feature_rows)
    by_bucket: dict[str, list[int]] = {}
    for position in scope:
        bucket = regime_bucket(float(feature_rows[position]["plddt"]))
        by_bucket.setdefault(bucket, []).append(position)

    scores: dict[int, float] = {}
    for position in scope:
        bucket = regime_bucket(float(feature_rows[position]["plddt"]))
        expert = directions[bucket]
        bucket_rows = {
            pos: feature_rows[pos] for pos in by_bucket[bucket]
        }
        terms: list[float] = []
        for feature, sign in expert.items():
            if sign == 0.0:
                continue
            values = {p: float(bucket_rows[p][feature]) for p in bucket_rows}
            terms.append(float(sign) * within_protein_z(values, position))
        scores[position] = float(np.mean(terms)) if terms else 0.0
    return scores


def routed_burden(
    records: Sequence[tuple[Mapping[int, Mapping[str, float]], int]],
    directions: Mapping[str, Mapping[str, float]],
    scope_mode: str,
) -> dict[str, Any]:
    """First-hit burden of true sites under the routed MoE score.

    For each protein, score ALL Cys with ``moe_route_scores`` and rank them
    descending; locate the true site. ``scope_mode`` selects which positions
    are scored/ranked:

    - ``"regime"``: score all positions but rank only the true site's own
      regime bucket (regime-local routing — the MoE claim);
    - ``"protein"``: rank among ALL Cys (routed scores across regimes compete
      on one scale — a stricter, arguably unfair bar since z-pools differ).

    Since expert directions are fixed priors (not learned), this needs no
    LOO sign learning — it is the routed analogue of
    ``structure_regime.loo_composite_burden`` with fixed signs.
    """
    if scope_mode not in ("protein", "regime"):
        raise ValueError(f"unknown scope_mode: {scope_mode!r}")
    per_protein: dict[int, dict[str, object]] = {}
    total = 0
    top1 = 0
    n = 0
    for index, (rows, true_position) in enumerate(records):
        scores = moe_route_scores(rows, directions)
        if scope_mode == "regime":
            bucket = regime_bucket(float(rows[true_position]["plddt"]))
            scope = subset_positions(rows, _in_regime_bucket(bucket))
            ranked_scope = sorted(scope, key=lambda p: scores[p], reverse=True)
            rank = ranked_scope.index(true_position) + 1
            n_in_scope = len(ranked_scope)
        else:
            ranked = sorted(scores, key=lambda p: scores[p], reverse=True)
            rank = ranked.index(true_position) + 1
            n_in_scope = len(scores)
        total += rank
        top1 += rank == 1
        n += 1
        per_protein[index] = {"rank": rank, "n_in_scope": n_in_scope}
    return {
        "n": n,
        "top1": top1,
        "total_first_hit_burden": total,
        "total_random_burden": round(
            sum((len(record[0]) + 1) / 2.0 for record in records), 3
        ),
        "per_protein": per_protein,
    }


def routed_permutation_null(
    records: Sequence[tuple[Mapping[int, Mapping[str, float]], int]],
    directions: Mapping[str, Mapping[str, float]],
    scope_mode: str,
    n_perm: int,
    rng: np.random.RandomState,
) -> list[int]:
    """Permutation null for the routed burden (exchangeability).

    Each permutation reassigns the true site to a uniform-random Cys of each
    protein and records the total routed burden. Because the directions are
    fixed priors, no sign re-learning happens inside the loop. Returns the
    sorted null distribution; the observed burden's left-tail position is its
    permutation p-value.
    """
    null: list[int] = []
    for _ in range(n_perm):
        perm_records = [
            (rows, int(rng.choice(list(rows)))) for rows, _true in records
        ]
        burden = routed_burden(perm_records, directions, scope_mode)[
            "total_first_hit_burden"
        ]
        null.append(int(burden))
    return sorted(null)
