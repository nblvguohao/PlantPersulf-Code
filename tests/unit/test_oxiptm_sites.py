"""Unit tests for the oxiPTM site table and fail-closed coordinate audit."""

from __future__ import annotations

from plantpersulf.evidence.oxiptm_sites import (
    ACCESSION_UNRESOLVED,
    ARABIDOPSIS,
    PERSULFIDATION,
    POSITION_NOT_CYS,
    RICE,
    S_NITROSYLATION,
    SITE_TABLE,
    TOMATO,
    UNMAPPABLE,
    VERIFIED,
    OxiptmSite,
    audit_sites,
    verified_sites,
    verify_site,
)


def _seq(length: int, residue_at: dict[int, str] | None = None) -> str:
    """A length-``length`` string; residues at 1-based positions overridden."""
    residue_at = residue_at or {}
    return "".join(residue_at.get(i, "A") for i in range(1, length + 1))


# --- verify_site status branches ---------------------------------------------


def test_verify_verified_when_position_is_cys() -> None:
    proteomes = {ARABIDOPSIS: {"P1": _seq(20, {10: "C"})}}
    site = OxiptmSite("G1", ("P1",), ARABIDOPSIS, 10, PERSULFIDATION, "t")
    result = verify_site(site, proteomes=proteomes)
    assert result.status == VERIFIED
    assert result.accession == "P1"
    assert result.observed_residue == "C"


def test_verify_position_not_cys() -> None:
    proteomes = {ARABIDOPSIS: {"P1": _seq(20, {10: "M"})}}
    site = OxiptmSite("G1", ("P1",), ARABIDOPSIS, 10, S_NITROSYLATION, "t")
    result = verify_site(site, proteomes=proteomes)
    assert result.status == POSITION_NOT_CYS
    assert result.observed_residue == "M"


def test_verify_unmappable_when_known_gene() -> None:
    proteomes: dict[str, dict[str, str]] = {TOMATO: {}}
    site = OxiptmSite("APX1", (), TOMATO, 168, PERSULFIDATION, "t")
    result = verify_site(site, proteomes=proteomes)
    assert result.status == UNMAPPABLE


def test_verify_accession_unresolved_without_candidates() -> None:
    proteomes: dict[str, dict[str, str]] = {ARABIDOPSIS: {}}
    site = OxiptmSite("PRMT5", (), ARABIDOPSIS, 125, S_NITROSYLATION, "t")
    result = verify_site(site, proteomes=proteomes)
    assert result.status == ACCESSION_UNRESOLVED


def test_verify_accession_unresolved_when_candidates_absent() -> None:
    proteomes = {ARABIDOPSIS: {"OTHER": _seq(10, {3: "C"})}}
    site = OxiptmSite("G1", ("P1", "P2"), ARABIDOPSIS, 3, S_NITROSYLATION, "t")
    result = verify_site(site, proteomes=proteomes)
    assert result.status == ACCESSION_UNRESOLVED


def test_verify_position_out_of_range() -> None:
    proteomes = {ARABIDOPSIS: {"P1": _seq(5)}}
    site = OxiptmSite("G1", ("P1",), ARABIDOPSIS, 50, PERSULFIDATION, "t")
    result = verify_site(site, proteomes=proteomes)
    assert result.status == POSITION_NOT_CYS


def test_verify_tries_candidates_in_order_first_verified_wins() -> None:
    # First candidate fails (non-Cys at 5), second verifies.
    proteomes = {
        ARABIDOPSIS: {"P1": _seq(10, {5: "M"}), "P2": _seq(10, {5: "C"})}
    }
    site = OxiptmSite("G1", ("P1", "P2"), ARABIDOPSIS, 5, S_NITROSYLATION, "t")
    result = verify_site(site, proteomes=proteomes)
    assert result.status == VERIFIED
    assert result.accession == "P2"


# --- audit_sites / verified_sites over a small table --------------------------


def test_audit_sites_returns_json_safe_rows() -> None:
    sites = (
        OxiptmSite("G1", ("P1",), ARABIDOPSIS, 5, PERSULFIDATION, "t"),
        OxiptmSite("G2", (), ARABIDOPSIS, 9, S_NITROSYLATION, "t"),
    )
    proteomes = {ARABIDOPSIS: {"P1": _seq(20, {5: "C"})}}
    rows = audit_sites(sites, proteomes=proteomes)
    assert [r["status"] for r in rows] == [VERIFIED, ACCESSION_UNRESOLVED]
    assert rows[0]["cys_position"] == 5
    assert rows[1]["observed_residue"] is None


def test_verified_sites_filters_to_passing_coordinates() -> None:
    sites = (
        OxiptmSite("G1", ("P1",), ARABIDOPSIS, 5, PERSULFIDATION, "t"),
        OxiptmSite("G2", ("P2",), ARABIDOPSIS, 9, S_NITROSYLATION, "t"),
    )
    proteomes = {
        ARABIDOPSIS: {
            "P1": _seq(20, {5: "C"}),
            "P2": _seq(20, {9: "M"}),
        }
    }
    assert [s.gene for s in verified_sites(sites, proteomes=proteomes)] == ["G1"]


# --- pin the documented real-world audit outcomes (synthetic proteome) --------
# These mirror the 2026-08-16 primary-source resolution of the review's Table 1.
# Verified SNO sites with their correct species accession: MPK6 C201 (Q39026),
# RAB7/RABG3E C171 (Q9XI98), SlMEK1 C172 (O48616), SlP5CR C5 (A0A3Q7FME1),
# ACOh4 C172 (A0A3Q7FZA2 = NCBI 'ACO homolog 4' LOC101265426),
# HA2 C206 (Q9SPD5 = NCBI LHA2). Verified persulfidation: bZIP68 C171
# (rice A2YXP7, the paper's dual-Cys171/245 pattern). The remaining SNO
# sites are reference-proteome coverage gaps and must NOT verify.


def _pin_proteome() -> dict[str, dict[str, str]]:
    return {
        ARABIDOPSIS: {
            "Q39026": _seq(210, {201: "C"}),  # MPK6 C201
            "Q9XI98": _seq(200, {171: "C"}),  # RABG3E C171
        },
        TOMATO: {
            "O48616": _seq(180, {172: "C"}),  # SlMEK1 C172
            "A0A3Q7FME1": _seq(280, {5: "C"}),  # SlP5CR C5
            "A0A3Q7FZA2": _seq(370, {172: "C"}),  # SlACOh4 C172
            "Q9SPD5": _seq(956, {206: "C"}),  # LHA2 C206
        },
        RICE: {
            "A2YXP7": _seq(435, {171: "C", 245: "C"}),  # bZIP68 dual-Cys
        },
    }


def test_resolved_sno_sites_verify() -> None:
    proteomes = _pin_proteome()
    cases = (
        ("MPK6", "Q39026", ARABIDOPSIS, 201),
        ("RAB7", "Q9XI98", ARABIDOPSIS, 171),
        ("MEK1", "O48616", TOMATO, 172),
        ("P5CR", "A0A3Q7FME1", TOMATO, 5),
        ("ACOh4", "A0A3Q7FZA2", TOMATO, 172),
        ("HA2", "Q9SPD5", TOMATO, 206),
    )
    for gene, accession, species, pos in cases:
        site = OxiptmSite(gene, (accession,), species, pos, S_NITROSYLATION, "t")
        assert verify_site(site, proteomes=proteomes).status == VERIFIED


def test_bzip68_c171_verifies_in_rice() -> None:
    proteomes = _pin_proteome()
    site = OxiptmSite("bZIP68", ("A2YXP7",), RICE, 171, PERSULFIDATION, "t")
    assert verify_site(site, proteomes=proteomes).status == VERIFIED


def test_coverage_gap_sites_not_verified() -> None:
    """Reference-proteome coverage gaps must stay unverified."""
    proteomes = _pin_proteome()
    unresolved = {"GSNOR", "LCD", "GSNOR1", "PRMT5"}
    for site in SITE_TABLE:
        if site.gene in unresolved:
            assert verify_site(site, proteomes=proteomes).status != VERIFIED


def test_site_table_contains_both_columns() -> None:
    n_persulf = sum(1 for s in SITE_TABLE if s.oxiptm == PERSULFIDATION)
    n_sno = sum(1 for s in SITE_TABLE if s.oxiptm == S_NITROSYLATION)
    assert n_persulf == 13
    assert n_sno == 10
