"""Gate-1 (oxidizable-Cys) site registry — stage-1 labels for the two-step model.

Mechanistic basis (Corpas et al., in-review COPLBI-D-26-00068, Fig.1C): a
cysteine thiol must first be oxidized to sulfenic acid (–SOH), a disulfide
(–S–S–), or an S-nitrosated form (–SNO) before it can react with H2S. The
two-step decomposition P(persulfidation) = P(oxidizable) × P(persulfidated |
oxidized) requires a stage-1 label set of Cys known to be oxidizable.

This module holds the *registered* gate-1 site table with provenance and a
fail-closed coordinate check, mirroring the oxiPTM-site discipline: a site
is ``verified`` only if its accession is in a registered reference proteome
AND the reported position is a Cys; otherwise it is recorded with its status
(``coverage_gap`` when the protein is absent from the registered proteome,
``unresolved`` otherwise).

Tier-1 content (this file): the S-nitrosylation column of the COPLBI review
Table 1 / body text — the regulatory-oxidation class that best matches the
"oxiPTMs compete for the same reactive Cys" framing and is the cleanest
gate-1 proxy. Disulfide-bonded Cys (structural oxidation, mostly secreted
proteins) and cross-species oxidation databases (human/mouse sulfenylation,
SNO) are documented as the stage-1 volume extension in the data-strategy
note; dbSNO is defunct and the plant sulfenome literature is sparse.

Claim class ``diagnostic_only``; nothing here trains a model.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from plantpersulf.evidence.oxiptm_sites import (
    S_NITROSYLATION,
    SITE_TABLE,
    VERIFIED,
    OxiptmSite,
    verify_site,
)

# Oxidation-state labels used as stage-1 classes.
SNO = "s-nitrosylation"
SULFENYLATION = "sulfenylation"
DISULFIDE = "disulfide"


@dataclass(frozen=True)
class Gate1Site:
    """One oxidizable-Cys label with provenance."""

    gene: str
    candidate_accessions: tuple[str, ...]
    species: str
    cys_position: int
    oxidation_type: str
    source: str
    source_note: str


# Tier-1: the S-nitrosylation sites sourced from the COPLBI review (the same
# 10 sites audited in oxiptm_sites).  Accessions are the resolved ones; the
# audit's ``verify_site`` decides verified vs coverage_gap per registered
# proteome.
TIER1_SNO: tuple[Gate1Site, ...] = tuple(
    Gate1Site(
        gene=site.gene,
        candidate_accessions=site.candidate_accessions,
        species=site.species,
        cys_position=site.cys_position,
        oxidation_type=SNO,
        source="COPLBI-D-26-00068 Table 1 / body text (S-nitrosation column)",
        source_note=site.source_note,
    )
    for site in SITE_TABLE
    if site.oxiptm == S_NITROSYLATION
)


def gate1_audit(
    sites: tuple[Gate1Site, ...] = TIER1_SNO,
    *,
    proteomes: dict[str, dict[str, str]],
) -> list[dict[str, Any]]:
    """Coordinate-verify every gate-1 site against the registered proteomes.

    Reuses the oxiPTM fail-closed ``verify_site``; statuses are ``verified``
    (Cys at position in a registered proteome), ``coverage_gap`` (protein
    absent / fragment / non-Cys in the registered set), or ``unresolved``.
    """
    rows: list[dict[str, Any]] = []
    for site in sites:
        result = verify_site(
            OxiptmSite(
                gene=site.gene,
                candidate_accessions=site.candidate_accessions,
                species=site.species,
                cys_position=site.cys_position,
                oxiptm=site.oxidation_type,
                source_note=site.source_note,
            ),
            proteomes=proteomes,
        )
        status = (
            "verified"
            if result.status == VERIFIED
            else "coverage_gap"
            if result.status in ("position_not_cys", "accession_unresolved")
            else "unmappable"
        )
        rows.append(
            {
                "gene": site.gene,
                "species": site.species,
                "cys_position": site.cys_position,
                "oxidation_type": site.oxidation_type,
                "accession": result.accession,
                "observed_residue": result.observed_residue,
                "status": status,
                "source": site.source,
                "source_note": site.source_note,
            }
        )
    return rows
