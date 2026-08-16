"""Scan PXD024061 MaxQuant sites tables for co-peptide negatives.

PXD024061 (Aroca et al. 2021, *Antioxidants* 10:508) deposited two MaxQuant
sites tables: ``Sulfide(C)Sites.txt`` (82 rows) and ``CianoBiotin(C)Sites.txt``
(6 rows), registered in ``data/registry/supplementary_sources.tsv``. Each row
is one candidate localisation; rows sharing a mod-peptide ID are merged into a
PeptideEvidence whose modification state is determinate (or provably
indeterminate) — the unmodified Cys of the same detected peptide are then
co-peptide negatives (detectability matched by construction).

Conservative rules (see ``maxquant_sites_negatives`` module docstring):
strong candidate (probability >= 0.75) -> positive; weak candidate ->
undetermined (never negative); non-candidate Cys -> negative only when the
total modified count of the peptide form is determinate. Peptides with no
determinate negative are not registered (no co-peptide contrast).

New rows are appended to ``data/registry/copeptide_negatives_v1.tsv``
(idempotent by site_id) and a scan summary is written to
``results/diagnostics/``.

Usage:
    python scripts/scan_copeptide_negatives_pxd024061.py
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import pandas as pd

from plantpersulf.evidence.copeptide_negatives import (
    load_copeptide_registry,
    write_copeptide_registry,
)
from plantpersulf.evidence.maxquant_sites_negatives import (
    PeptideEvidence,
    classify_evidence,
    evidence_from_site_rows,
    locate_peptide,
)
from plantpersulf.features.sequence import _load_proteome

_REPO_ROOT = Path(__file__).resolve().parents[1]
SUPPLEMENTS = _REPO_ROOT / "data" / "raw" / "supplements" / "PXD024061"
REGISTRY = _REPO_ROOT / "data" / "registry" / "copeptide_negatives_v1.tsv"
PROTEOME = (
    _REPO_ROOT / "data" / "raw" / "references" / "arabidopsis_ref_proteome_v2.fasta"
)
OUTPUT = _REPO_ROOT / "results" / "diagnostics" / "copeptide_negatives_pxd024061.json"

DOI = "10.3390/antiox10040508"

PROBES = ("Sulfide(C)", "CianoBiotin(C)")


def _peptide_id(row: pd.Series) -> str:
    """Group key: mod-peptide ID with fallbacks to peptide ID and row index."""
    for column in ("Mod. peptide IDs", "Peptide IDs"):
        value = row.get(column)
        if value is not None and str(value) not in ("", "nan", "<NA>"):
            return f"{column}:{value}"
    return f"row:{row.name}"


def _evidence_groups(
    path: Path, probe: str
) -> list[tuple[str, PeptideEvidence, str]]:
    """Group sites-table rows by peptide; return (key, evidence, proteins)."""
    table = pd.read_csv(path, sep="\t", dtype=str)
    probability_column = f"{probe} Probabilities"
    number_column = f"Number of {probe}"
    groups: dict[str, list[dict[str, object]]] = {}
    proteins_by_key: dict[str, str] = {}
    order: list[str] = []
    for _, row in table.iterrows():
        key = _peptide_id(row)
        if key not in groups:
            groups[key] = []
            order.append(key)
            proteins_by_key[key] = str(row["Proteins"])
        groups[key].append(
            {
                probability_column: row[probability_column],
                number_column: row[number_column],
            }
        )
    results: list[tuple[str, PeptideEvidence, str]] = []
    for key in order:
        try:
            evidence = evidence_from_site_rows(
                groups[key],
                probability_column=probability_column,
                number_column=number_column,
            )
        except (RuntimeError, ValueError) as error:
            raise RuntimeError(f"{path.name} group {key}: {error}") from error
        results.append((key, evidence, proteins_by_key[key]))
    return results


def _registry_rows(
    probe: str,
    key: str,
    evidence: PeptideEvidence,
    proteins_field: str,
    proteome: dict[str, str],
    counters: Counter,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for protein in proteins_field.split(";"):
        accession = protein.strip()
        if accession.startswith("CON__"):
            counters["skipped_decoy_or_contaminant"] += 1
            continue
        if accession not in proteome:
            counters["skipped_protein_not_in_registered_proteome"] += 1
            continue
        start = locate_peptide(evidence.sequence, proteome[accession])
        if start is None:
            counters["skipped_peptide_not_uniquely_located"] += 1
            continue
        states = classify_evidence(evidence)
        negatives = {
            position for position, state in states.items() if state == "negative"
        }
        if not negatives:
            counters["skipped_peptide_without_determinate_negative"] += 1
            continue
        counters["peptides_registered"] += 1
        end = start + len(evidence.sequence) - 1
        by_state = {
            state: sorted(
                start + position - 1
                for position, label in states.items()
                if label == state
            )
            for state in ("positive", "negative", "undetermined")
        }
        pairs = ";".join(
            f"{position}:{probability:.3f}"
            for position, probability in evidence.candidates
        )
        for position, state in sorted(states.items()):
            rows.append(
                {
                    "site_id": (
                        f"PXD024061_{probe.split('(')[0]}_"
                        f"{accession}_{start + position - 1}_{state.upper()}"
                    ),
                    "species": "arabidopsis",
                    "protein_accession": accession,
                    "cys_position": start + position - 1,
                    "state": state,
                    "peptide_sequence": evidence.sequence,
                    "peptide_start": start,
                    "peptide_end": end,
                    "in_peptide_position": position,
                    "positive_positions": ";".join(map(str, by_state["positive"])),
                    "negative_positions": ";".join(map(str, by_state["negative"])),
                    "undetermined_positions": ";".join(
                        map(str, by_state["undetermined"])
                    ),
                    "localization_confirmed": (
                        "True" if by_state["positive"] else "False"
                    ),
                    "site_determining_ion_coverage": (
                        "True" if evidence.n_modified is not None else "False"
                    ),
                    "source_doi": DOI,
                    "source_detail": f"PXD024061 {probe.split('(')[0]}(C)Sites.txt",
                    "provenance": (
                        f"PXD024061 MaxQuant {probe.split('(')[0]}(C)Sites.txt "
                        f"(Aroca et al. 2021, Antioxidants 10:508), modified "
                        f"peptide group {key}, peptide "
                        f"{evidence.sequence} re-localised by unique match at "
                        f"{start}-{end} of {accession} in registered arabidopsis "
                        f"reference proteome v2 (MaxQuant search coordinates "
                        f"differ from the registered proteome space); total "
                        f"modified Cys {evidence.n_modified}; candidate "
                        f"localisation probabilities {pairs}; same peptide, "
                        f"same spectrum, same enrichment -> detectability "
                        f"matched by construction (AGENTS.md explicit-negative "
                        f"registration); raw spectra not locally available "
                        f"(limitation)."
                    ),
                }
            )
    return rows


def main() -> None:
    proteome = _load_proteome(PROTEOME)
    existing = load_copeptide_registry(REGISTRY)
    existing_ids = {record["site_id"] for record in existing}

    counters: Counter = Counter()
    new_rows: list[dict[str, str]] = []
    for probe in PROBES:
        path = SUPPLEMENTS / f"{probe.split('(')[0]}_C_Sites.txt"
        counters[f"probe:{probe.split('(')[0]}_files"] += 1
        for key, evidence, proteins_field in _evidence_groups(path, probe):
            counters["peptides_scanned"] += 1
            rows = _registry_rows(
                probe, key, evidence, proteins_field, proteome, counters
            )
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
        "rows_before": len(existing),
        "rows_added": len(new_rows),
        "rows_after": len(existing) + len(new_rows),
        "counters": dict(counters),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
