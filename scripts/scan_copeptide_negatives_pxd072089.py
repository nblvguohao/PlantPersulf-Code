#!/usr/bin/env python
"""Scan PXD072089 PNAS supplements for co-peptide negatives (rice).

PXD072089 (Huang et al., *PNAS* 2025, ``10.1073/pnas.2608150123``) is a rice
(*Oryza sativa*) persulfidome study — a 4th species for the co-peptide
negative axis and an independent laboratory from kiad070 / PXD024061, which
is exactly the Gate 2 condition-1 independence the negative dimension needs.

Two supplementary tables are combined (see ``pnas_site_negatives`` module):
Dataset S4 (paper-confirmed ``-SSH`` sites, one row per modified Cys) and
Dataset S1 (MaxQuant peptide inventory — every Cys of each detected
persulfidated peptide). Where a detected peptide carries more Cys than the
paper confirms modified, the confirmed positions are registered ``positive``
and the remaining Cys of the same detected peptide are registered
``negative`` — the same paper-reported basis as the registered kiad070
tomato entries (unambiguous published localisation on the same peptide; raw
spectra not locally available is recorded as a limitation).

Peptides are re-localised by unique string match in the registered rice
reference proteome v1 (UniProt Oryza sativa species-level); the paper's
search coordinates occasionally differ by a residue, so the paper-space
modified positions are mapped to in-peptide positions via the S1 Cys
inventory (coordinate-shift invariant) and the peptide's own Cys ranks.

New rows are appended to ``data/registry/copeptide_negatives_v1.tsv``
(idempotent by site_id) and a scan summary is written to
``results/diagnostics/``.

Usage:
    python scripts/scan_copeptide_negatives_pxd072089.py
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import pandas as pd

from plantpersulf.evidence.copeptide_negatives import (
    classify_peptide_cys,
    load_copeptide_registry,
    write_copeptide_registry,
)
from plantpersulf.evidence.maxquant_sites_negatives import locate_peptide
from plantpersulf.evidence.pnas_site_negatives import (
    co_peptide_groups,
    parse_position_field,
)
from plantpersulf.features.sequence import _load_proteome

_REPO_ROOT = Path(__file__).resolve().parents[1]
SUPPLEMENTS = _REPO_ROOT / "data" / "raw" / "supplements" / "PXD072089"
REGISTRY = _REPO_ROOT / "data" / "registry" / "copeptide_negatives_v1.tsv"
RICE_PROTEOME = (
    _REPO_ROOT
    / "data"
    / "raw"
    / "references"
    / "rice_proteome_v1"
    / "uniprot_rice_v1.fasta"
)
OUTPUT = _REPO_ROOT / "results" / "diagnostics" / "copeptide_negatives_pxd072089.json"

DOI = "10.1073/pnas.2608150123"
SITES_FILE = SUPPLEMENTS / "pnas.2608150123.sd04.xlsx"
INVENTORY_FILE = SUPPLEMENTS / "pnas.2608150123.sd01.xlsx"


def load_sites(path: Path) -> list[dict]:
    """Read Dataset S4 into site records (paper-confirmed modified Cys)."""
    frame = pd.read_excel(path, sheet_name=0, header=1)
    rows: list[dict] = []
    for _, raw in frame.iterrows():
        rows.append(
            {
                "peptide": str(raw["Peptide sequence"]),
                "protein_ids": [
                    part.strip()
                    for part in str(raw["Protein_ID"]).split(";")
                    if part.strip()
                ],
                "leading_razor": str(raw["leading razor protein"]),
                "modified_positions": parse_position_field(
                    raw["-SSH Cys position"]
                ),
                "domain": str(raw["Domain"]) if pd.notna(raw["Domain"]) else "",
            }
        )
    return rows


def load_inventory(path: Path) -> dict[tuple[str, str], tuple[int, ...]]:
    """Read Dataset S1 into ``{(peptide, leading_razor): cys_inventory}``."""
    frame = pd.read_excel(path, sheet_name=0, header=1)
    result: dict[tuple[str, str], tuple[int, ...]] = {}
    for _, raw in frame.iterrows():
        key = (str(raw["Sequence"]), str(raw["Leading razor protein"]))
        result[key] = parse_position_field(raw["Position"])
    return result


def _registry_rows(
    group: dict,
    proteome: dict[str, str],
    counters: Counter,
) -> list[dict[str, str]]:
    """Build registry rows for one co-peptide group, per registered protein."""
    peptide = group["peptide"]
    modified_in_peptide = group["modified_in_peptide"]
    rows: list[dict[str, str]] = []
    for accession in group["protein_ids"]:
        if accession.startswith("CON__"):
            counters["skipped_decoy_or_contaminant"] += 1
            continue
        if accession not in proteome:
            counters["skipped_protein_not_in_registered_proteome"] += 1
            continue
        start = locate_peptide(peptide, proteome[accession])
        if start is None:
            counters["skipped_peptide_not_uniquely_located"] += 1
            continue
        states = classify_peptide_cys(
            peptide=peptide,
            modified_in_peptide_positions=modified_in_peptide,
            localization_confirmed=True,
            site_determining_ions_confirmed=True,
        )
        negatives = {
            position for position, state in states.items() if state == "negative"
        }
        if not negatives:
            counters["skipped_peptide_without_determinate_negative"] += 1
            continue
        counters["peptides_registered"] += 1
        end = start + len(peptide) - 1
        by_state = {
            state: sorted(
                start + position - 1
                for position, label in states.items()
                if label == state
            )
            for state in ("positive", "negative", "undetermined")
        }
        domain = next(iter(group["domains"]), "")
        for position, state in sorted(states.items()):
            rows.append(
                {
                    "site_id": (
                        f"PXD072089_{accession}_"
                        f"{start + position - 1}_{state.upper()}"
                    ),
                    "species": "rice",
                    "protein_accession": accession,
                    "cys_position": start + position - 1,
                    "state": state,
                    "peptide_sequence": peptide,
                    "peptide_start": start,
                    "peptide_end": end,
                    "in_peptide_position": position,
                    "positive_positions": ";".join(map(str, by_state["positive"])),
                    "negative_positions": ";".join(map(str, by_state["negative"])),
                    "undetermined_positions": ";".join(
                        map(str, by_state["undetermined"])
                    ),
                    "localization_confirmed": "True",
                    "site_determining_ion_coverage": "True",
                    "source_doi": DOI,
                    "source_detail": "PNAS suppl. Dataset S4 (pnas.2608150123)",
                    "provenance": (
                        f"PXD072089 PNAS supplementary Dataset S4 "
                        f"(pnas.2608150123, Huang et al.), persulfidated sites "
                        f"with corresponding domain, row for peptide {peptide} "
                        f"(domain {domain or 'n/a'}); paper-confirmed -SSH Cys "
                        f"position(s) {sorted(modified_in_peptide)}; peptide "
                        f"re-localised by unique match at {start}-{end} of "
                        f"{accession} in registered rice reference proteome "
                        f"v1 (UniProt Oryza sativa species-level); other Cys "
                        f"of the same detected peptide observed unmodified -> "
                        f"same peptide, same spectrum, same enrichment -> "
                        f"detectability matched by construction (AGENTS.md "
                        f"explicit-negative registration); "
                        f"site_determining_ion_coverage=True as "
                        f"paper-reported unambiguous localisation of the "
                        f"co-modified sites on the same peptide; raw spectra "
                        f"not locally available (limitation)."
                    ),
                }
            )
    return rows


def main() -> None:
    proteome = _load_proteome(RICE_PROTEOME)
    existing = load_copeptide_registry(REGISTRY)
    existing_ids = {record["site_id"] for record in existing}

    sites = load_sites(SITES_FILE)
    inventory = load_inventory(INVENTORY_FILE)
    groups = co_peptide_groups(sites, inventory)

    counters: Counter = Counter()
    counters["sites_rows"] = len(sites)
    counters["inventory_peptides"] = len(inventory)
    counters["co_peptide_groups_found"] = len(groups)

    new_rows: list[dict[str, str]] = []
    for group in groups:
        rows = _registry_rows(group, proteome, counters)
        for row in rows:
            if row["site_id"] in existing_ids:
                counters["skipped_already_registered"] += 1
                continue
            new_rows.append(row)
            existing_ids.add(row["site_id"])

    if new_rows:
        write_copeptide_registry(REGISTRY, existing + new_rows)

    summary = {
        "registry": str(REGISTRY.relative_to(_REPO_ROOT)),
        "source": {
            "sites": str(SITES_FILE.relative_to(_REPO_ROOT)),
            "inventory": str(INVENTORY_FILE.relative_to(_REPO_ROOT)),
            "doi": DOI,
        },
        "rows_before": len(existing),
        "rows_added": len(new_rows),
        "rows_after": len(existing) + len(new_rows),
        "counters": dict(counters),
        "groups": [
            {
                "peptide": group["peptide"],
                "protein_ids": group["protein_ids"],
                "modified_in_peptide": list(group["modified_in_peptide"]),
                "negative_in_peptide": list(group["negative_in_peptide"]),
            }
            for group in groups
        ],
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
