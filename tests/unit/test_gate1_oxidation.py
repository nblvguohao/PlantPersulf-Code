"""Unit tests for the gate-1 (oxidizable-Cys) registry module."""

from __future__ import annotations

from plantpersulf.evidence.gate1_oxidation import (
    SNO,
    TIER1_SNO,
    Gate1Site,
    gate1_audit,
)
from plantpersulf.evidence.oxiptm_sites import (
    ARABIDOPSIS,
    S_NITROSYLATION,
    SITE_TABLE,
)


def _seq(length: int, residue_at: dict[int, str] | None = None) -> str:
    residue_at = residue_at or {}
    return "".join(residue_at.get(i, "A") for i in range(1, length + 1))


def test_tier1_contains_all_review_sno_sites() -> None:
    sno_from_table = [s for s in SITE_TABLE if s.oxiptm == S_NITROSYLATION]
    assert len(TIER1_SNO) == len(sno_from_table) == 10
    assert all(site.oxidation_type == SNO for site in TIER1_SNO)


def test_gate1_audit_verified_status() -> None:
    proteomes = {ARABIDOPSIS: {"Q39026": _seq(210, {201: "C"})}}
    site = Gate1Site(
        "MPK6", ("Q39026",), ARABIDOPSIS, 201, SNO, "review", "t"
    )
    rows = gate1_audit((site,), proteomes=proteomes)
    assert rows[0]["status"] == "verified"
    assert rows[0]["accession"] == "Q39026"


def test_gate1_audit_coverage_gap() -> None:
    proteomes: dict[str, dict[str, str]] = {ARABIDOPSIS: {}}
    site = Gate1Site("PRMT5", (), ARABIDOPSIS, 125, SNO, "review", "t")
    rows = gate1_audit((site,), proteomes=proteomes)
    assert rows[0]["status"] == "coverage_gap"
