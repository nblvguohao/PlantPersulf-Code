"""Compartment conditioning overlay — subcellular pH/redox of within-protein ranks.

Mechanistic motivation (Corpas et al., in-review COPLBI-D-26-00068, and the
two-step persulfidation constraint it states): NO/H2S signalling is
compartmentalized, and the thiolate fraction of a cysteine of fixed pKa
scales with the local pH — the illuminated chloroplast stroma (~8.0) vs the
cytosol (~7.2) spans roughly an order of magnitude in reactivity. A
compartment-conditioned score is a registered next-generation feature
candidate (model-improvement register L10).

This module provides ONLY the zero-cost, in-repo slice of that feature: the
compartment is taken from the registered reference-proteome headers and the
registered tomato GO annotation file (taxon 4081), never asserted from
memory. If neither source annotates a protein, the compartment is reported
as ``None`` (not_annotated) — a full overlay would require a localization
prediction (TargetP/DeepLoc), which is a new data input, not part of this
diagnostic.

Claim class ``diagnostic_only``; nothing fits, mutates, or touches frozen
artifacts.
"""

from __future__ import annotations

import re
from typing import Any

# --- curated compartment pH/redox (literature values, used only as an
#     interpretive overlay, not a model input) --------------------------------
# pH sources: cytosol ~7.2 (plant cytosolic pH), illuminated chloroplast
# stroma ~8.0 (dark ~7.4), mitochondrial matrix ~7.8, nucleus ~7.3,
# ER ~7.2, Golgi ~6.4, vacuole ~5.0, apoplast ~5.5-6.5, peroxisome ~7.4.
# "plasma_membrane" refers to the cytosolic face (pH ~7.2).
COMPARTMENT_PH: dict[str, float] = {
    "cytosol": 7.2,
    "cytoplasm": 7.2,
    "nucleus": 7.3,
    "chloroplast": 8.0,
    "plastid": 7.8,
    "mitochondrion": 7.8,
    "peroxisome": 7.4,
    "ER": 7.2,
    "Golgi": 6.4,
    "vacuole": 5.0,
    "apoplast": 5.8,
    "plasma_membrane": 7.2,
    "membrane": 7.2,
    "extracellular": 6.0,
}

COMPARTMENT_REDOX_NOTE: dict[str, str] = {
    "cytosol": "reducing (GSH/GRX), light-independent",
    "chloroplast": "strongly oxidising under illumination (Fd-TRX pool)",
    "plastid": "oxidising during photosynthesis",
    "mitochondrion": "matrix more reducing, IMS oxidising",
    "peroxisome": "high H2O2, oxidising",
    "apoplast": "oxidising (H2O2-generating)",
    "vacuole": "low pH, low redox activity",
    "nucleus": "mildly oxidising (Trx import)",
    "ER": "oxidising (ERO1, disulfide formation)",
    "Golgi": "mildly oxidising",
}

# --- GO cellular-component ids -> compartment --------------------------------
GO_CC_TO_COMPARTMENT: dict[str, str] = {
    "GO:0005634": "nucleus",
    "GO:0005737": "cytoplasm",
    "GO:0005829": "cytosol",
    "GO:0009507": "chloroplast",
    "GO:0009536": "plastid",
    "GO:0005739": "mitochondrion",
    "GO:0005777": "peroxisome",
    "GO:0005783": "ER",
    "GO:0005794": "Golgi",
    "GO:0005773": "vacuole",
    "GO:0005886": "plasma_membrane",
    "GO:0016020": "membrane",
    "GO:0048046": "apoplast",
    "GO:0005618": "apoplast",  # cell wall (apoplastic face)
    "GO:0005576": "extracellular",
}

# Primary-compartment tie-break when a protein carries several CC terms:
# the most pH-relevant (distinctive) compartment wins over generic ones.
_COMPARTMENT_PRIORITY: tuple[str, ...] = (
    "apoplast",
    "vacuole",
    "chloroplast",
    "plastid",
    "peroxisome",
    "mitochondrion",
    "Golgi",
    "ER",
    "nucleus",
    "cytosol",
    "cytoplasm",
    "plasma_membrane",
    "membrane",
    "extracellular",
)

# Ordered header keyword scan: specific compartments first, so e.g.
# "chloroplast membrane" reads as chloroplast, not membrane.
_HEADER_KEYWORDS: tuple[tuple[str, str], ...] = (
    ("chloroplast", "chloroplast"),
    ("mitochondri", "mitochondrion"),
    ("peroxisom", "peroxisome"),
    ("vacuol", "vacuole"),
    ("endoplasmic", "ER"),
    ("golgi", "Golgi"),
    ("nucle", "nucleus"),
    ("cytosolic", "cytosol"),
    ("cytoplasmic", "cytoplasm"),
    ("apoplast", "apoplast"),
    ("cell wall", "apoplast"),
    ("extracellular", "extracellular"),
    ("secreted", "extracellular"),
    ("membrane", "membrane"),
    ("transmembrane", "membrane"),
)

_GO_TERM_RE = re.compile(r"(?P<term>.+?)\s*\[(?P<go>GO:\d+)\]")


def parse_go_field(field: str | None) -> list[tuple[str, str]]:
    """Parse a ``"term [GO:id]; term [GO:id]; ..."`` GO field.

    Returns ``[(term, go_id), ...]`` in file order; empty for None/empty.
    """
    if not field:
        return []
    terms: list[tuple[str, str]] = []
    for chunk in field.split(";"):
        chunk = chunk.strip()
        match = _GO_TERM_RE.match(chunk)
        if match:
            terms.append((match.group("term").strip(), match.group("go")))
    return terms


def cellular_component_compartments(
    go_terms: list[tuple[str, str]],
) -> list[str]:
    """The compartments implied by the CC terms of a GO field, in file order."""
    return [
        GO_CC_TO_COMPARTMENT[go]
        for _, go in go_terms
        if go in GO_CC_TO_COMPARTMENT
    ]


def header_compartment(header: str) -> str | None:
    """Subcellular compartment from a reference-proteome FASTA header, or None."""
    lowered = header.lower()
    for keyword, compartment in _HEADER_KEYWORDS:
        if keyword in lowered:
            return compartment
    return None


def primary_compartment(compartments: list[str]) -> str | None:
    """The most pH-relevant compartment by ``_COMPARTMENT_PRIORITY``."""
    for candidate in _COMPARTMENT_PRIORITY:
        if candidate in compartments:
            return candidate
    return compartments[0] if compartments else None


def compartment_ph(compartment: str | None) -> float | None:
    """pH overlay for a compartment name, or None if unannotated/unknown."""
    if compartment is None:
        return None
    return COMPARTMENT_PH.get(compartment)


def compartment_overlay(
    *,
    go_field: str | None,
    header: str | None,
) -> dict[str, Any]:
    """Combine the two in-repo localization sources into one overlay row.

    GO CC terms (registered tomato annotation) take precedence over the
    header keyword scan. ``not_annotated`` is the honest outcome when neither
    source carries a compartment.
    """
    compartments: list[str] = []
    if go_field:
        compartments = cellular_component_compartments(parse_go_field(go_field))
    if not compartments and header:
        header_comp = header_compartment(header)
        if header_comp is not None:
            compartments = [header_comp]
    primary = primary_compartment(compartments)
    return {
        "compartments": compartments,
        "primary": primary,
        "ph": compartment_ph(primary),
        "redox_note": COMPARTMENT_REDOX_NOTE.get(primary) if primary else None,
        "annotated": primary is not None,
        "source": "go_cc" if compartments and go_field else (
            "header" if primary is not None else "not_annotated"
        ),
    }
