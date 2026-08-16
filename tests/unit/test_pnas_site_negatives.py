"""Unit tests for pnas_site_negatives.py.

Covers the coordinate-shift-safe mapping of paper-confirmed modified Cys to
in-peptide positions, the S1 inventory cross-check, and the fail-closed
co-peptide group construction.
"""

from __future__ import annotations

import pytest

from plantpersulf.evidence.copeptide_negatives import classify_peptide_cys
from plantpersulf.evidence.pnas_site_negatives import (
    co_peptide_groups,
    cys_in_peptide_positions,
    map_modified_positions,
    parse_position_field,
)

# The Q5ZCB1 case: paper inventory off by +1 vs registered coordinates.
# SLPPICHCADEVASCAAACKECDMVNSSSEPPR, Cys at in-peptide 6, 8, 15, 19, 22.
PEPTIDE = "SLPPICHCADEVASCAAACKECDMVNSSSEPPR"
INVENTORY_PAPER = (144, 146, 153, 157, 160)  # S1 paper-space positions
CYS_IN_PEP = (6, 8, 15, 19, 22)


class TestParsePositionField:
    def test_comma_separated(self) -> None:
        assert parse_position_field("413,417") == (413, 417)

    def test_semicolon_separated(self) -> None:
        assert parse_position_field("144;146;153") == (144, 146, 153)

    def test_single_value(self) -> None:
        assert parse_position_field("21") == (21,)

    def test_numeric_cell(self) -> None:
        assert parse_position_field(21) == (21,)

    def test_empty_and_none(self) -> None:
        assert parse_position_field("") == ()
        assert parse_position_field(None) == ()


class TestCysInPeptidePositions:
    def test_counts_cys(self) -> None:
        assert cys_in_peptide_positions("ACPGVCR") == (2, 6)

    def test_no_cys(self) -> None:
        assert cys_in_peptide_positions("AVGLPR") == ()

    def test_real_peptide(self) -> None:
        assert cys_in_peptide_positions(PEPTIDE) == CYS_IN_PEP


class TestMapModifiedPositions:
    def test_rank_mapping_invariant_to_coordinate_shift(self) -> None:
        # modified paper positions 146,153,157 -> ranks 2,3,4 -> in-peptide 8,15,19
        mapped = map_modified_positions(
            (146, 153, 157), INVENTORY_PAPER, CYS_IN_PEP, PEPTIDE
        )
        assert mapped == (8, 15, 19)

    def test_single_modified_first_cys(self) -> None:
        # ACPTDVLEMIPWDGCK: inventory (21, 34), modified 21 -> rank 1 -> in-peptide 2
        assert map_modified_positions(
            (21,), (21, 34), (2, 15), "ACPTDVLEMIPWDGCK"
        ) == (2,)

    def test_single_modified_second_cys(self) -> None:
        # MAGPGVCINLLNGTTMHLSVGCVFR: inventory (7, 22), modified 22 -> in-peptide 22
        assert map_modified_positions(
            (22,), (7, 22), (7, 22), "MAGPGVCINLLNGTTMHLSVGCVFR"
        ) == (22,)

    def test_missing_position_fails_closed(self) -> None:
        with pytest.raises(ValueError, match="not in peptide inventory"):
            map_modified_positions((999,), (21, 34), (2, 15), "ACPTDVLEMIPWDGCK")

    def test_inventory_peptide_cys_mismatch_fails_closed(self) -> None:
        with pytest.raises(ValueError, match="inventory.*disagree"):
            map_modified_positions((21,), (21, 34), (2,), "ACPTDVLEMIPWDGCK")


class TestCoPeptideGroups:
    def _site(self, peptide, razor, positions, protein_ids=("P1",)):
        return {
            "peptide": peptide,
            "protein_ids": list(protein_ids),
            "leading_razor": razor,
            "modified_positions": tuple(positions),
            "domain": "",
        }

    def test_partial_peptide_produces_group_with_negative(self) -> None:
        sites = [self._site(PEPTIDE, "Q5ZCB1", (146, 153, 157))]
        inventory = {(PEPTIDE, "Q5ZCB1"): INVENTORY_PAPER}
        groups = co_peptide_groups(sites, inventory)
        assert len(groups) == 1
        group = groups[0]
        assert group["modified_in_peptide"] == (8, 15, 19)
        assert group["negative_in_peptide"] == (6, 22)

    def test_multi_row_sites_merge_modified_set(self) -> None:
        sites = [
            self._site(PEPTIDE, "Q5ZCB1", (146,)),
            self._site(PEPTIDE, "Q5ZCB1", (153,)),
            self._site(PEPTIDE, "Q5ZCB1", (157,)),
        ]
        inventory = {(PEPTIDE, "Q5ZCB1"): INVENTORY_PAPER}
        groups = co_peptide_groups(sites, inventory)
        assert len(groups) == 1
        assert groups[0]["modified_in_peptide"] == (8, 15, 19)

    def test_fully_modified_peptide_skipped(self) -> None:
        # AAVCEMPFATVASDDLGGVGGTCVLR: 2 Cys, both modified -> no negative
        peptide = "AAVCEMPFATVASDDLGGVGGTCVLR"
        sites = [
            self._site(peptide, "Q8S5T1", (107,)),
            self._site(peptide, "Q8S5T1", (126,)),
        ]
        inventory = {(peptide, "Q8S5T1"): (107, 126)}
        assert co_peptide_groups(sites, inventory) == []

    def test_single_cys_peptide_skipped(self) -> None:
        sites = [self._site("LVTLCDNAPATPIDVVR", "A0A0N7KE93", (106,))]
        inventory = {("LVTLCDNAPATPIDVVR", "A0A0N7KE93"): (106,)}
        assert co_peptide_groups(sites, inventory) == []

    def test_missing_inventory_fails_closed(self) -> None:
        sites = [self._site(PEPTIDE, "Q5ZCB1", (146,))]
        with pytest.raises(RuntimeError, match="no Dataset S1 inventory"):
            co_peptide_groups(sites, {})


class TestClassificationIntegration:
    def test_q5zcb1_three_state_labels(self) -> None:
        # Regression of the real Q5ZCB1 case with the shared classifier.
        states = classify_peptide_cys(
            peptide=PEPTIDE,
            modified_in_peptide_positions=(8, 15, 19),
            localization_confirmed=True,
            site_determining_ions_confirmed=True,
        )
        assert states == {6: "negative", 8: "positive", 15: "positive",
                          19: "positive", 22: "negative"}
