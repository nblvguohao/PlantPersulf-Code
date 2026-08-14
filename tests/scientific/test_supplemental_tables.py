"""Tests for the supplemental-table generator (Supplemental Table S1).

S1 documents the release training panel's per-species composition; the
counts asserted here are pinned to the frozen release manifest
(multispecies-v2-candidate-release-v1, SHA256 30e432fb...).
"""

from pathlib import Path

import pytest

from scripts.generate_supplemental_tables import (
    SPECIES_ORDER,
    composition_from_seed_panels,
)

REPO = Path(__file__).resolve().parents[2]
SEED_GLOB = (
    REPO
    / "results/experiments/multispecies_v2_global_clusters_v11"
    / "literature_random_protein/sul_bertgru_seed*/shared_panel.tsv"
)


@pytest.fixture(scope="module")
def composition():
    return composition_from_seed_panels(str(SEED_GLOB))


def test_union_panel_totals_match_frozen_release(composition):
    assert composition["n_rows"] == 389609
    assert composition["n_positives"] == 2334


def test_per_species_counts_match_candidate_scan_table(composition):
    expected = {
        "arabidopsis": (111496, 327),
        "rice": (185379, 739),
        "tomato": (72051, 79),
        "magnaporthe": (20683, 1189),
    }
    for sp, (rows, pos) in expected.items():
        assert composition["species"][sp]["rows"] == rows
        assert composition["species"][sp]["positives"] == pos


def test_species_order_is_reporting_order(composition):
    assert list(composition["species"]) == SPECIES_ORDER
    assert sum(s["rows"] for s in composition["species"].values()) == 389609
    assert sum(s["positives"] for s in composition["species"].values()) == 2334
