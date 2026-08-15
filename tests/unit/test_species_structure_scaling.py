"""Unit tests for per-species structure feature scaling (diagnostic B).

Covers ``evaluation.species_structure_scaling`` — the campaign-only transform
that standardizes structure features *within each species* before the frozen
global ``_BranchScalers`` runs. Because each species is z-scored to (0, 1),
the subsequent global fit on merged rows is approximately the identity, so the
per-species alignment survives the ranker's internal second scaling.

This is the diagnosis for the W1 mechanism lead: the global structure scaler
is dominated by tomato rows (78% of structure rows after release v3), so
arabidopsis/rice structure inputs get re-scaled against tomato statistics.
Per-species scaling removes that cross-species pull.
"""

from __future__ import annotations

import pytest

from plantpersulf.evaluation.species_structure_scaling import (
    fit_species_struct_scalers,
    transform_species_struct,
)


def test_fit_returns_one_scaler_per_species_with_present_rows() -> None:
    rows = [
        [1.0, 50.0],  # tomato 1
        [2.0, 60.0],  # tomato 2
        [9.0, 80.0],  # rice 1
        [11.0, 90.0],  # rice 2
        [0.0, 0.0],  # magnaporthe (masked, no structure)
    ]
    masks = [True, True, True, True, False]
    species = ["tomato", "tomato", "rice", "rice", "magnaporthe"]
    scalers = fit_species_struct_scalers(rows, masks, species)
    assert set(scalers) == {"tomato", "rice"}
    # tomato mean over present rows: [(1+2)/2, (50+60)/2] = [1.5, 55]
    assert scalers["tomato"].mean == pytest.approx((1.5, 55.0))
    assert scalers["rice"].mean == pytest.approx((10.0, 85.0))


def test_transform_applies_per_species_scaler() -> None:
    rows = [
        [1.0, 50.0],
        [2.0, 60.0],
        [9.0, 80.0],
        [11.0, 90.0],
    ]
    masks = [True, True, True, True]
    species = ["tomato", "tomato", "rice", "rice"]
    scalers = fit_species_struct_scalers(rows, masks, species)
    out = transform_species_struct(rows, masks, species, scalers)
    # Each species' present rows are z-scored: within species, first value
    # below the species mean -> negative, second above -> positive.
    assert out[0][0] < 0.0 and out[1][0] > 0.0  # tomato
    assert out[2][0] < 0.0 and out[3][0] > 0.0  # rice


def test_transform_keeps_masked_rows_untouched() -> None:
    rows = [
        [1.0, 50.0],
        [2.0, 60.0],
        [0.0, 0.0],  # masked row
    ]
    masks = [True, True, False]
    species = ["tomato", "tomato", "magnaporthe"]
    scalers = fit_species_struct_scalers(rows, masks, species)
    out = transform_species_struct(rows, masks, species, scalers)
    assert out[2] == [0.0, 0.0]  # untouched raw zeros


def test_transform_falls_back_to_identity_for_species_without_scaler() -> None:
    rows = [
        [1.0, 50.0],
        [5.0, 70.0],  # species with no scaler fit on train
    ]
    masks = [True, True]
    species = ["tomato", "novel"]
    scalers = fit_species_struct_scalers(rows[:1], masks[:1], species[:1])
    out = transform_species_struct(rows, masks, species, scalers)
    assert out[1] == [5.0, 70.0]  # identity for unknown species


def test_all_masked_train_produces_no_scalers() -> None:
    rows = [[0.0, 0.0], [0.0, 0.0]]
    masks = [False, False]
    species = ["magnaporthe", "magnaporthe"]
    assert fit_species_struct_scalers(rows, masks, species) == {}
