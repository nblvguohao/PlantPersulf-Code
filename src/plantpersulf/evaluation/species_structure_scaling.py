"""Per-species structure feature scaling (diagnostic B, campaign-only).

The W1 structure-coverage evaluation surfaced a mechanism lead: the frozen
``_BranchScalers.struct`` is fit on *all* present-structure rows pooled, and
after release v3 tomato rows dominate (~78% of structure rows), so the global
scaler statistics are tomato-driven and arabidopsis/rice structure inputs get
re-scaled against tomato statistics — consistent with the large per-species
AP regression of those two crops (W1 report §3.2).

This module is the campaign-only transform that standardizes structure
features *within each species* before the ranker's global scaler runs:

- ``fit_species_struct_scalers`` fits one ``TrainOnlyScaler`` per species on
  that species' present-structure rows only (masked rows never contribute).
- ``transform_species_struct`` applies each species' scaler to its rows;
  rows whose species has no fitted scaler (no present rows in train) are
  returned unchanged (identity fallback); masked rows are never touched.

Why this survives the ranker's internal second global fit: after per-species
z-scoring each species has column mean 0 and variance 1, so a subsequent
global fit over merged rows computes mean ≈ 0 and variance ≈ 1 for every
column — i.e. the internal transform is approximately the identity and does
not re-introduce the cross-species pull. This is exactly the opposite of the
W1 failure mode (raw values scaled against tomato-dominated statistics).

Claim class ``diagnostic_only``; nothing here trains a model or modifies the
frozen bundle path (``_BranchScalers`` is untouched).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

from plantpersulf.models.traditional import TrainOnlyScaler


def fit_species_struct_scalers(
    rows: Sequence[Sequence[float]],
    masks: Sequence[bool],
    species: Sequence[str],
) -> dict[str, TrainOnlyScaler]:
    """Fit one structure scaler per species over its present rows.

    Species with no present-structure rows get no scaler (the caller falls
    back to identity). Masked rows are excluded from every species'
    statistics, matching the frozen ``_BranchScalers`` invariant.
    """
    by_species: dict[str, list[list[float]]] = {}
    for row, present, name in zip(rows, masks, species, strict=True):
        if present:
            by_species.setdefault(name, []).append(list(row))
    return {
        name: TrainOnlyScaler.fit(rows)
        for name, rows in by_species.items()
    }


def transform_species_struct(
    rows: Sequence[Sequence[float]],
    masks: Sequence[bool],
    species: Sequence[str],
    scalers: Mapping[str, TrainOnlyScaler],
) -> list[list[float]]:
    """Transform rows with their species' scaler (identity if none fitted).

    Masked rows are returned untouched (their structure values are zeroed by
    the branch mask downstream anyway, but we never mutate them here). Rows
    whose species is not in ``scalers`` are returned as-is — no fabricated
    statistics.
    """
    out: list[list[float]] = []
    for row, present, name in zip(rows, masks, species, strict=True):
        scaler = scalers.get(name)
        if not present or scaler is None:
            out.append(list(row))
        else:
            out.append(scaler.transform([list(row)])[0])
    return out
