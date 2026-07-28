#!/usr/bin/env python
"""PXD072089 attrition compositional-bias check (self-review response,
2026-07-23).

Question: is the ~55% UniProt-accession-deletion attrition in the rice
persulfidome dataset (``docs/phase_z_evidence_audit.md`` §5.3) redox-biased
— i.e. does it disproportionately drop non-redox proteins (which would
artificially redox-enrich the surviving/verified set and confound the
§5.5.1 family-level and §5.5.2 structural-level convergence findings), or
disproportionately drop redox proteins (deflating those findings)?

Method: every unique protein accession referenced in SD01 or SD04 carries a
UniProt-style free-text "Protein description" field, captured by the
original MS search *before* any later accession deletion — this field
survives even for accessions no longer resolvable in the current reference
proteome. Each accession is classified:

* **kept**   — accession resolves in the current rice reference proteome
  (``RICE_PROTEOME``), i.e. survives to be coordinate-verifiable and usable
  downstream (as in ``pxd072089_sites.parse_pxd072089_sites``).
* **dropped** — accession does not resolve (deleted from UniProt since the
  original MS search).

A two-sided Fisher's exact test on the 2x2 (kept/dropped) x
(redox-related/not) table (``evaluation/attrition_bias_check.py``)
determines whether the two groups' redox-keyword composition differs.

Usage::

    python scripts/check_attrition_bias.py \
        --output results/cross_species_conservation/attrition_bias_v1.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

SD01 = Path("data/raw/supplements/PXD072089/pnas.2608150123.sd01.xlsx")
SD04 = Path("data/raw/supplements/PXD072089/pnas.2608150123.sd04.xlsx")
RICE_PROTEOME = Path("data/raw/references/rice_proteome_v1/uniprot_rice_v1.fasta")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _collect_accession_descriptions(
    sd01_path: Path, sd04_path: Path
) -> dict[str, str]:
    """One description per unique accession referenced in SD01 or SD04.

    SD01 uses "Leading razor protein"; SD04 uses "leading razor protein"
    (lowercase). Both carry a "Protein description" column. Accessions with
    an empty accession field are skipped (they carry no identity to
    classify); first non-empty description wins if an accession appears
    with more than one description string (should not happen for a stable
    UniProt entry, but the loop is defensive).
    """
    from plantpersulf.proteomics.pxd072089_sites import _read_xlsx

    descriptions: dict[str, str] = {}
    for path, accession_col in (
        (sd01_path, "Leading razor protein"),
        (sd04_path, "leading razor protein"),
    ):
        for row in _read_xlsx(path):
            accession = str(row.get(accession_col, "") or "").strip()
            if not accession:
                continue
            description = str(row.get("Protein description", "") or "").strip()
            if not description:
                continue
            descriptions.setdefault(accession, description)
    return descriptions


def run_attrition_bias_check(output_path: Path) -> dict:
    from plantpersulf.evaluation.attrition_bias_check import attrition_bias_test
    from plantpersulf.features.sequence import _load_proteome

    for path in (SD01, SD04, RICE_PROTEOME):
        if not path.is_file():
            raise RuntimeError(f"required input missing: {path}")

    proteome = _load_proteome(RICE_PROTEOME)
    descriptions = _collect_accession_descriptions(SD01, SD04)

    kept_descriptions = []
    dropped_descriptions = []
    kept_accessions = []
    dropped_accessions = []
    for accession, description in sorted(descriptions.items()):
        if accession in proteome:
            kept_descriptions.append(description)
            kept_accessions.append(accession)
        else:
            dropped_descriptions.append(description)
            dropped_accessions.append(accession)

    result = attrition_bias_test(kept_descriptions, dropped_descriptions)

    summary = {
        "question": (
            "Is the UniProt-accession-deletion attrition in PXD072089 "
            "redox-keyword-biased (kept vs. dropped proteins)?"
        ),
        "n_unique_accessions_referenced": len(descriptions),
        "n_kept": len(kept_accessions),
        "n_dropped": len(dropped_accessions),
        "attrition_fraction": (
            len(dropped_accessions) / len(descriptions) if descriptions else 0.0
        ),
        "kept_redox": result.kept_redox,
        "kept_non_redox": result.kept_non_redox,
        "dropped_redox": result.dropped_redox,
        "dropped_non_redox": result.dropped_non_redox,
        "kept_redox_fraction": result.kept_redox_fraction,
        "dropped_redox_fraction": result.dropped_redox_fraction,
        "odds_ratio_kept_vs_dropped": result.odds_ratio,
        "p_value_two_sided": result.p_value,
        "significant_0.05": result.p_value < 0.05,
        "interpretation": (
            "kept_redox_fraction > dropped_redox_fraction and significant "
            "would indicate the surviving/verified dataset is artificially "
            "redox-enriched by differential attrition (a real confound for "
            "the family/structural convergence findings); "
            "kept_redox_fraction < dropped_redox_fraction and significant "
            "would indicate the attrition removed redox proteins "
            "disproportionately (deflating those findings, not inflating "
            "them); p >= 0.05 indicates no detectable compositional bias."
        ),
        "inputs_sha256": {
            "sd01_xlsx": _sha256(SD01),
            "sd04_xlsx": _sha256(SD04),
            "rice_proteome": _sha256(RICE_PROTEOME),
        },
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print(f"\nwrote -> {output_path}")
    return summary


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="PXD072089 attrition bias check")
    p.add_argument(
        "--output",
        type=Path,
        default=Path("results/cross_species_conservation/attrition_bias_v1.json"),
    )
    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    run_attrition_bias_check(args.output)
