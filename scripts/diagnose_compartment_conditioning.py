#!/usr/bin/env python
"""Compartment conditioning overlay for the registered known controls.

Tabulates the subcellular pH/redox overlay (mechanistic two-step
persulfidation constraint) against the frozen-bundle within-protein ranks of
the registered controls and the one verified S-nitrosylation site (MPK6).
Compartment comes ONLY from in-repo sources — the registered tomato GO
annotation (taxon 4081) and the reference-proteome FASTA headers — and is
reported as not_annotated when neither source carries it. This is the
zero-cost slice of the compartment-conditioning feature candidate (register
L10): if coverage is low, the honest finding is that a full overlay needs a
localization prediction (TargetP/DeepLoc), which is a new data input.

diagnostic_only; no fit, no frozen-artifact change, no claim beyond the table.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from plantpersulf.evaluation.compartment_conditioning import (  # noqa: E402
    compartment_overlay,
)
from plantpersulf.evaluation.known_controls import REGISTERED_CONTROLS  # noqa: E402

TOMATO_GO = (
    _REPO_ROOT
    / "data"
    / "raw"
    / "references"
    / "go_annotations_v1"
    / "tomato_go_taxon4081.tsv"
)
ARABIDOPSIS_PROTEOME = (
    _REPO_ROOT / "data" / "raw" / "references" / "arabidopsis_ref_proteome_v2.fasta"
)
TOMATO_PROTEOME = (
    _REPO_ROOT / "data" / "raw" / "references" / "tomato_ref_proteome_v1.fasta"
)
RANKING_JSON = (
    _REPO_ROOT / "results" / "known_controls" / "within_protein_ranking_v1.json"
)
OUTPUT = _REPO_ROOT / "results" / "diagnostics" / "compartment_conditioning_v1.json"
OUTPUT_ROWS = (
    _REPO_ROOT / "results" / "diagnostics" / "compartment_conditioning_v1.rows.tsv"
)

# The single verified S-nitrosylation site (oxiptm_sites_v1 audit) appended as
# an out-of-group contrast row.
MPK6 = {
    "gene": "MPK6",
    "uniprot_accession": "Q39026",
    "species": "arabidopsis",
    "cys_position": 201,
    "control_site": 201,
    "note": "verified S-nitrosylation site (oxiptm_sites_v1)",
}


def _load_go(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    with path.open(encoding="utf-8") as handle:
        handle.readline()  # skip the "Entry\tGene Ontology (GO)" header
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                result[parts[0]] = "\t".join(parts[1:])
    return result


def _load_headers(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for line in path.open(encoding="utf-8"):
        if line.startswith(">"):
            header = line.strip()
            accession = header.split()[0].split("|")[1]
            result[accession] = header
    return result


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args(argv)

    go = _load_go(TOMATO_GO)
    headers: dict[str, dict[str, str]] = {
        "arabidopsis": _load_headers(ARABIDOPSIS_PROTEOME),
        "tomato": _load_headers(TOMATO_PROTEOME),
    }
    ranking = json.loads(RANKING_JSON.read_text(encoding="utf-8"))
    ranking_by_gene = {r["gene"]: r for r in ranking["controls"]}

    rows: list[dict[str, object]] = []
    controls = [c for c in REGISTERED_CONTROLS if c.status == "mapped"]
    entries: list[dict[str, Any]] = [
        *(
            {
                "gene": c.gene,
                "uniprot_accession": c.uniprot_accession,
                "species": ("tomato" if c.control_species == "Solanum lycopersicum"
                            else "arabidopsis"),
                "cys_position": c.cys_position,
                "control_site": c.cys_position,
                "note": "registered persulfidation control",
            }
            for c in controls
        ),
        MPK6,
    ]

    for entry in entries:
        accession = entry["uniprot_accession"]
        species = entry["species"]
        go_field = go.get(accession) if species == "tomato" else None
        header = headers[species].get(accession)
        overlay = compartment_overlay(go_field=go_field, header=header)

        rank: int | None = None
        burden: float | None = None
        baseline: float | None = None
        if entry["note"] == "registered persulfidation control":
            rec = ranking_by_gene.get(entry["gene"])
            if rec is not None:
                reg = rec["registered_metrics"]
                # JSON object keys are strings; true_site_ranks is keyed by
                # the control's registered position.
                rank = reg["true_site_ranks"][str(entry["cys_position"])]
                burden = reg["first_hit_burden"]
                baseline = reg["random_baseline_burden"]

        row = {
            "gene": entry["gene"],
            "species": species,
            "protein_accession": accession,
            "cys_position": entry["cys_position"],
            "note": entry["note"],
            "compartments": overlay["compartments"],
            "primary_compartment": overlay["primary"],
            "compartment_ph": overlay["ph"],
            "redox_note": overlay["redox_note"],
            "source": overlay["source"],
            "within_protein_rank": rank,
            "first_hit_burden": burden,
            "random_baseline_burden": baseline,
        }
        rows.append(row)

    n_annotated = sum(1 for r in rows if r["primary_compartment"] is not None)
    annotated = [r for r in rows if r["primary_compartment"] is not None]

    print("=== compartment overlay (zero-cost, in-repo sources only) ===")
    print(
        f"{'gene':<9}{'compartment':<15}{'pH':>5}  {'rank':>5}  "
        f"{'burden':>7}  source"
    )
    for r in rows:
        ph = f"{r['compartment_ph']:.1f}" if r["compartment_ph"] is not None else "-"
        rank_str = (
            f"{r['within_protein_rank']}"
            if r["within_protein_rank"] is not None else "-"
        )
        burden_str = (
            f"{r['first_hit_burden']}"
            if r["first_hit_burden"] is not None else "-"
        )
        print(
            f"{r['gene']:<9}{(r['primary_compartment'] or 'not_annotated'):<15}"
            f"{ph:>5}  {rank_str:>5}  {burden_str:>7}  {r['source']}"
        )

    document = {
        "track": "compartment_conditioning_v1",
        "claim_class": "diagnostic_only",
        "note": (
            "Zero-cost compartment overlay (registered tomato GO CC + "
            "reference-proteome headers) against frozen within-protein ranks "
            "of the registered controls and the sole verified S-nitrosylation "
            "site (MPK6). Compartment is never asserted from memory; a protein "
            "with no in-repo annotation is not_annotated. A full overlay "
            "requires a localization prediction (TargetP/DeepLoc) — a new data "
            "input for the compartment-conditioning feature candidate (L10)."
        ),
        "n_controls": len(controls),
        "n_annotated_in_repo": n_annotated,
        "coverage": f"{n_annotated}/{len(rows)}",
        "rows": rows,
        "reading": (
            "In-repo localization annotation covers only a minority of the "
            "controls (tomato GO file + a single header match), so no "
            "compartment-conditioned statement about the within-protein ranks "
            "is possible at this coverage — the reference headers do not carry "
            "subcellular annotation for most Arabidopsis controls. This "
            "confirms the compartment overlay is a new-data-input feature "
            "(TargetP/DeepLoc), not a free add-on, and the pH/redox table is "
            "ready to apply once localizations are predicted."
            if n_annotated < len(rows) // 2
            else "Coverage sufficient for a descriptive overlay."
        ),
        "annotated_rows": annotated,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(document, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    with OUTPUT_ROWS.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = list(rows[0])
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for r in rows:
            writer.writerow({k: (";".join(v) if isinstance(v, list) else v)
                            for k, v in r.items()})

    print(f"\ncoverage: {n_annotated}/{len(rows)} annotated in-repo")
    print(f"wrote {args.output}")
    print(f"wrote {OUTPUT_ROWS}")


if __name__ == "__main__":
    main()
