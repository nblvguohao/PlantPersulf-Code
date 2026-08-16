"""OxiPTM site audit — the two oxiPTM columns of the in-review Corpas review.

The premise to test (methodology design: a multi-oxiPTM head is a registered
candidate; the review's own framing is that "oxiPTMs compete for the same
reactive Cys residues within proteins"): can a persulfidation-trained frozen
model separate persulfidation sites from S-nitrosylation sites, so that "the
other modification" becomes a more reliable negative than a putatively
unmodified site? Before any scoring, every proposed site must survive a
fail-closed coordinate check against the registered reference proteomes —
the same discipline that resolved BRG3, WRKY71, RNF144b and ERF.D3, where
cross-paper numbering or isoform offsets repeatedly broke naive transfer.

This module therefore holds the site table (both columns of the review's
Table 1) and the pure verification logic, and nothing more. Scoring lives in
``scripts/evaluate_oxiptm_sites.py``. Claim class ``diagnostic_only``; no
frozen artifact is touched and no site enters the formal registry until it
verifies and carries a primary-source DOI chain.

Source: Corpas et al., COPLBI-D-26-00068 (in-review Current Opinion in Plant
Biology), Table 1 — persulfidation column and S-nitrosation column, with the
persulfidation accessions cross-checked against ``known_controls.py`` where
already registered.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

# Species keys match the proteome dicts used by the release-scoring path.
ARABIDOPSIS = "arabidopsis"
TOMATO = "tomato"
RICE = "rice"

PERSULFIDATION = "persulfidation"
S_NITROSYLATION = "s-nitrosylation"

# Status values from verify_site().
VERIFIED = "verified"
POSITION_NOT_CYS = "position_not_cys"
ACCESSION_UNRESOLVED = "accession_unresolved"
UNMAPPABLE = "unmappable"

# Documented resolution failures from known_controls.py provenance: tomato
# APX1 (candidate accessions carry Ala, not Cys, at position 168) and POD5
# (zero Solanum lycopersicum hits) cannot be registered without the primary
# report's own accession-number section.
KNOWN_UNMAPPABLE = ("APX1", "POD5")


@dataclass(frozen=True)
class OxiptmSite:
    """One proposed oxiPTM site from the review's Table 1.

    ``candidate_accessions`` is ordered: verification accepts the first that
    survives; empty means no accession could be resolved in the registered
    reference proteome.
    """

    gene: str
    candidate_accessions: tuple[str, ...]
    species: str
    cys_position: int
    oxiptm: str
    source_note: str


@dataclass(frozen=True)
class Verification:
    """Outcome of a fail-closed coordinate check."""

    status: str
    accession: str | None
    observed_residue: str | None
    detail: str


SITE_TABLE: tuple[OxiptmSite, ...] = (
    # ---- persulfidation column (accessions cross-checked against
    #      known_controls.py where registered) ----------------------------
    OxiptmSite(
        "DES1", ("F4K5T2",), ARABIDOPSIS, 44, PERSULFIDATION,
        "review Table 1; Shen et al. 2020 Plant Cell (doi:10.1105/tpc.19.00826), "
        "registered control F4K5T2 C44",
    ),
    OxiptmSite(
        "DES1", ("F4K5T2",), ARABIDOPSIS, 205, PERSULFIDATION,
        "paired second site, same report (not scored separately in the control)",
    ),
    OxiptmSite(
        "RBOHD", ("Q9FIJ0",), ARABIDOPSIS, 825, PERSULFIDATION,
        "review Table 1; Shen et al. 2020 Plant Cell (doi:10.1105/tpc.19.00826), "
        "registered control Q9FIJ0 C825",
    ),
    OxiptmSite(
        "RBOHD", ("Q9FIJ0",), ARABIDOPSIS, 890, PERSULFIDATION,
        "paired second site, same report",
    ),
    OxiptmSite(
        "SnRK2.6", ("Q940H6",), ARABIDOPSIS, 131, PERSULFIDATION,
        "review Table 1; Chen et al. Mol Plant (doi:10.1016/j.molp.2021.07.002), "
        "registered control Q940H6 C131",
    ),
    OxiptmSite(
        "SnRK2.6", ("Q940H6",), ARABIDOPSIS, 137, PERSULFIDATION,
        "paired second site, same report",
    ),
    OxiptmSite(
        "ATG4a", ("Q8S929",), ARABIDOPSIS, 170, PERSULFIDATION,
        "review Table 1; Laureano-Marin et al. 2020 Plant Cell "
        "(doi:10.1105/tpc.20.00766), registered control Q8S929 C170",
    ),
    OxiptmSite(
        "ABI4", ("A0MES8",), ARABIDOPSIS, 250, PERSULFIDATION,
        "review Table 1; Zhou et al. 2021 Mol Plant (doi:10.1016/j.molp.2021.03.007), "
        "registered control A0MES8 C250",
    ),
    OxiptmSite(
        "PAD3", ("Q9LW27",), ARABIDOPSIS, 440, PERSULFIDATION,
        "review Table 1; Zhang et al. 2026 Plant Cell Environ "
        "(doi:10.1111/pce.70593), registered control Q9LW27 C440",
    ),
    OxiptmSite(
        "bZIP68", ("A2YXP7",), RICE, 171, PERSULFIDATION,
        "review Table 1 / body text; RICE japonica cv. Wuyungeng 7, Ma X et al. "
        "2026 (IJMS 27:3841, doi:10.3390/ijms27093841). RESOLVED (open-access "
        "full text): persulfidation at Cys171 + oxidation at Cys245; A2YXP7 "
        "(OsI_30119, 435aa, BZIP domain) is the registered proteome's only "
        "protein with BOTH cysteines at exactly 171 and 245 (= NCBI 125562410). "
        "Caveat: entry is indica-annotated; the japonica locus ID was not "
        "resolvable (Ensembl down, no NCBI gene record) but is the sole "
        "representative of this bZIP68 locus in the reference",
    ),
    OxiptmSite(
        "CAT1", ("P30264",), TOMATO, 234, PERSULFIDATION,
        "review Table 1; Li et al. 2020 Plant Physiol Biochem 156:257-266, "
        "registered control P30264 C234",
    ),
    OxiptmSite(
        "APX1", (), TOMATO, 168, PERSULFIDATION,
        "review Table 1; known-unmappable in tomato (known_controls.py "
        "provenance: candidate accessions carry Ala at 168)",
    ),
    OxiptmSite(
        "POD5", (), TOMATO, 61, PERSULFIDATION,
        "review Table 1; known-unmappable in tomato (zero Solanum "
        "lycopersicum hits)",
    ),
    # ---- S-nitrosation column.  Species attribution follows the review's
    #      body text (primary source per site); the earlier pass wrongly
    #      assumed Arabidopsis for every symbol (2026-08-16 resolution). ----
    OxiptmSite(
        "ACOh4", ("A0A3Q7FZA2",), TOMATO, 172, S_NITROSYLATION,
        "review Table 1 / body text; TOMATO cv. Ailsa Craig, Liu et al. 2023 "
        "(New Phytol 239:159-173, doi:10.1111/nph.18928). RESOLVED via NCBI "
        "Gene LOC101265426 '1-aminocyclopropane-1-carboxylate oxidase homolog "
        "4' (protein 460382410, 366aa): A0A3Q7FZA2 in the registered tomato "
        "proteome is an exact match and carries Cys at 172 (annotated as "
        "'Fe2OG dioxygenase', which is why the name scan missed it)",
    ),
    OxiptmSite(
        "MEK1", ("O48616", "A0A3Q7JS13"), TOMATO, 172, S_NITROSYLATION,
        "review Table 1 / body text; TOMATO cv. Micro-Tom, Fang et al. 2026 — "
        "tomato MAPKK, Cys verified at 172",
    ),
    OxiptmSite(
        "GSNOR1", ("Q0WM36",), ARABIDOPSIS, 10, S_NITROSYLATION,
        "review Table 1 / body text; Arabidopsis (gsnor1-3), Zhan et al. 2018 "
        "(hypoxia). Coverage gap: canonical GSNOR1 absent from proteome v2; "
        "only 195-residue fragment Q0WM36 (Pro at 10)",
    ),
    OxiptmSite(
        "GSNOR", (), TOMATO, 47, S_NITROSYLATION,
        "review Table 1 / body text; TOMATO cv. Micro-Tom, Huang et al. 2026 "
        "(Cd stress). Coverage gap: tomato S-nitrosoglutathione reductase "
        "absent from registered tomato proteome",
    ),
    OxiptmSite(
        "LCD", (), TOMATO, 225, S_NITROSYLATION,
        "review Table 1 / body text; TOMATO cv. Micro-Tom, Huang et al. 2026. "
        "Coverage gap: tomato L-cysteine desulfhydrase absent (only "
        "D-cysteine desulfhydrase present)",
    ),
    OxiptmSite(
        "P5CR", ("A0A3Q7FME1",), TOMATO, 5, S_NITROSYLATION,
        "review Table 1 / body text; TOMATO SlP5CR (pyrroline-5-carboxylate "
        "reductase, SlP5CRC5S mutants), Liu et al. 2024 — Cys verified at 5",
    ),
    OxiptmSite(
        "HA2", ("Q9SPD5",), TOMATO, 206, S_NITROSYLATION,
        "review Table 1 / body text; TOMATO cv. Ailsa Craig, Wei et al. 2025 "
        "(Plant Cell 37:koaf035, doi:10.1093/plcell/koaf035). RESOLVED: "
        "Q9SPD5 (GN=LHA2, 956aa) is 100% identical to NCBI 'plasma membrane "
        "H+-ATPase isoform LHA2' (5901757) and carries Cys at 206; the "
        "704-residue P23980 fragment with Gly at 206 is the old truncated "
        "record of the same gene",
    ),
    OxiptmSite(
        "RAB7", ("Q9XI98",), ARABIDOPSIS, 171, S_NITROSYLATION,
        "review Table 1 / body text; RABG3E (aka RAB7, rab7-1/rab7-2), "
        "Lin et al. 2023 — Cys verified at 171",
    ),
    OxiptmSite(
        "PRMT5", (), ARABIDOPSIS, 125, S_NITROSYLATION,
        "review Table 1 / body text; Arabidopsis (prmt5-1), Hu et al. 2017. "
        "Coverage gap: PRMT5/SKB1 absent from reference proteome v2",
    ),
    OxiptmSite(
        "MPK6", ("Q39026",), ARABIDOPSIS, 201, S_NITROSYLATION,
        "review Table 1 / body text; Arabidopsis (mpk6-3), Wang et al. 2025 — "
        "Cys verified at 201",
    ),
)


def verify_site(
    site: OxiptmSite,
    *,
    proteomes: dict[str, dict[str, str]],
    known_unmappable: tuple[str, ...] = KNOWN_UNMAPPABLE,
) -> Verification:
    """Fail-closed coordinate check of one proposed site.

    Returns ``VERIFIED`` only when some candidate accession is present in the
    registered reference proteome AND the reported position is a Cys. All
    other outcomes carry the observed residue and a reason, so a failed audit
    row distinguishes "wrong residue at the position" (isoform/numbering
    offset, the BRG3-class failure) from "accession unresolvable" (absent
    entry / unindexed locus) from "documented unmappable" (gene-level
    resolution failure already recorded in known_controls.py provenance).
    """
    if not site.candidate_accessions:
        if site.gene in known_unmappable:
            return Verification(
                UNMAPPABLE, None, None,
                f"{site.gene}: known-unmappable (documented in known_controls.py)",
            )
        return Verification(
            ACCESSION_UNRESOLVED, None, None,
            f"{site.gene}: no accession resolved in reference proteome",
        )

    species_proteome = proteomes.get(site.species, {})
    # Every candidate isoform is scanned: verification succeeds if ANY
    # candidate carries a Cys at the reported position (the paper may have
    # meant that isoform); otherwise the FIRST non-Cys finding is reported,
    # so a multi-isoform gene whose candidates all fail is one finding, not
    # several.
    first_non_cys: tuple[str, str] | None = None
    for accession in site.candidate_accessions:
        sequence = species_proteome.get(accession)
        if sequence is None:
            continue
        if not 1 <= site.cys_position <= len(sequence):
            if first_non_cys is None:
                first_non_cys = (accession, "?")
            continue
        residue = sequence[site.cys_position - 1]
        if residue == "C":
            return Verification(
                VERIFIED,
                accession,
                "C",
                f"{accession} residue C at {site.cys_position}",
            )
        if first_non_cys is None:
            first_non_cys = (accession, residue)
    if first_non_cys is not None:
        accession, residue = first_non_cys
        return Verification(
            POSITION_NOT_CYS, accession, residue,
            f"{accession} has {residue} (not Cys) at {site.cys_position}",
        )
    return Verification(
        ACCESSION_UNRESOLVED, None, None,
        f"none of {list(site.candidate_accessions)} present in reference proteome",
    )


def audit_sites(
    sites: tuple[OxiptmSite, ...] = SITE_TABLE,
    *,
    proteomes: dict[str, dict[str, str]],
) -> list[dict[str, Any]]:
    """Audit every proposed site; returns one dict per site (JSON-safe)."""
    rows: list[dict[str, Any]] = []
    for site in sites:
        result = verify_site(site, proteomes=proteomes)
        rows.append(
            {
                "gene": site.gene,
                "oxiptm": site.oxiptm,
                "species": site.species,
                "cys_position": site.cys_position,
                "status": result.status,
                "accession": result.accession,
                "observed_residue": result.observed_residue,
                "detail": result.detail,
                "source_note": site.source_note,
            }
        )
    return rows


def verified_sites(
    sites: tuple[OxiptmSite, ...] = SITE_TABLE,
    *,
    proteomes: dict[str, dict[str, str]],
) -> list[OxiptmSite]:
    """The subset that passes the coordinate check (the scoreable set)."""
    return [
        site
        for site in sites
        if verify_site(site, proteomes=proteomes).status == VERIFIED
    ]
