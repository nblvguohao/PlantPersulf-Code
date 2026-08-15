#!/usr/bin/env python
"""Diagnostic-only: within-protein ranking of every registered control.

The release's published headline was "0 of 99 published sites appear in the
candidate table's top 2,000" — a proteome-wide Top-K percentile. But the
four published workflows (kiad068/070/100/271) never blind-scan the
proteome: they first select a target protein (transcriptome / interaction /
phenotype), then ask "which Cys of THIS protein". Both statements can hold
simultaneously, because the model was trained on within-protein pairwise
ranking, not on cross-protein score comparability.

This diagnostic scores EVERY Cys of each mapped registered control protein
with the FROZEN candidate-release bundle (structure branch masked — the
release's tomato scoring semantics) and reports the within-protein rank of
the published site(s), with Top-1 / Hit@2 / Hit@3 / MRR / first-hit
Mutagenesis Burden vs the random baseline E=(n+1)/(k+1).

Pre-registered discipline: this is a mock-blind expectation-setting
diagnostic (same claim class as the kiae271 release-bundle recovery). It
changes no model, feature, Top-K, threshold or SAP; published positives are
excluded from the SAP primary endpoint.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from plantpersulf.evaluation.known_controls import REGISTERED_CONTROLS  # noqa: E402
from plantpersulf.evaluation.release_scoring import (  # noqa: E402
    cys_positions,
    score_sites,
)
from plantpersulf.evaluation.within_protein_ranking import (  # noqa: E402
    within_protein_metrics,
)
from plantpersulf.features.sequence import _load_proteome  # noqa: E402
from plantpersulf.models.structure_ranker import StructureRankerBundle  # noqa: E402

BUNDLE_PATH = (
    _REPO_ROOT
    / "results"
    / "candidates"
    / "multispecies_v2_candidate_release_v1"
    / "model_weights"
    / "structure_ranker_bundle.pt"
)
TOMATO_PROTEOME = (
    _REPO_ROOT / "data" / "raw" / "references" / "tomato_ref_proteome_v1.fasta"
)
ARABIDOPSIS_PROTEOME = (
    _REPO_ROOT / "data" / "raw" / "references" / "arabidopsis_ref_proteome_v2.fasta"
)
DEFAULT_OUTPUT = (
    _REPO_ROOT / "results" / "known_controls" / "within_protein_ranking_v1.json"
)

# Second sites of paired-mutant lineages documented in known_controls.py
# provenance ("second site in the same paired report, not scored separately"
# convention). Used for the k>=2 group-variant metrics only; the primary
# metrics use the registered representative position (k=1, release
# single-position-per-lineage convention).
GROUP_SITE_VARIANTS: dict[str, tuple[int, ...]] = {
    "BRG3_H2S_UBIQUITINATION": (206, 212),
    "WRKY71_H2S_UBIQUITINATION": (35, 40),
    "DES1_H2S_SELF_RBOHD_ABA": (44, 205),
    "RBOHD_H2S_ROS_ABA": (825, 890),
    "SNRK26_H2S_PHOSPHORYLATION": (131, 137),
}

SPECIES_PROTEOMES = {
    "Solanum lycopersicum": ("tomato", TOMATO_PROTEOME),
    "Arabidopsis thaliana": ("arabidopsis", ARABIDOPSIS_PROTEOME),
}


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)

    controls = [c for c in REGISTERED_CONTROLS if c.status == "mapped"]
    print(f"mapped controls: {len(controls)}")

    bundle = StructureRankerBundle.load(BUNDLE_PATH)
    proteomes: dict[str, dict[str, str]] = {}
    by_species: dict[str, list] = {}
    for control in controls:
        species, proteome_path = SPECIES_PROTEOMES[control.control_species]
        proteomes.setdefault(species, _load_proteome(proteome_path))
        by_species.setdefault(species, []).append(control)

    results: list[dict[str, object]] = []
    for species, species_controls in by_species.items():
        for control in species_controls:
            sequence = proteomes[species][control.uniprot_accession]
            positions = cys_positions(sequence)
            sites = [(species, control.uniprot_accession, p) for p in positions]
            scored = score_sites(sites, bundle=bundle, proteomes=proteomes)
            scores = {
                p: scored[(species, control.uniprot_accession, p)].score
                for p in positions
            }

            primary = within_protein_metrics(
                positions=positions,
                scores=scores,
                true_positions={control.cys_position},
            )
            variants = GROUP_SITE_VARIANTS.get(control.mechanism_lineage_id)
            group = None
            if variants is not None:
                group = within_protein_metrics(
                    positions=positions,
                    scores=scores,
                    true_positions=set(variants),
                )

            record = {
                "mechanism_lineage_id": control.mechanism_lineage_id,
                "gene": control.gene,
                "species": species,
                "protein_accession": control.uniprot_accession,
                "protein_length": len(sequence),
                "control_site": control.cys_position,
                "registered_metrics": primary,
                "group_variant_metrics": group,
            }
            results.append(record)
            reg = primary
            print(
                f"  {control.gene} Cys{control.cys_position} (n={reg['n_cys']}, "
                f"k=1): rank={reg['true_site_ranks'][control.cys_position]}, "
                f"top1={reg['top1']}, hit@2={reg['hit_at_2']}, "
                f"burden={reg['first_hit_burden']} "
                f"(random {reg['random_baseline_burden']})"
            )
            if group:
                print(
                    f"    group k={len(variants)}: ranks={group['true_site_ranks']}, "
                    f"burden={group['first_hit_burden']} "
                    f"(random {group['random_baseline_burden']})"
                )

    summary = {
        "track": "within_protein_ranking_v1",
        "claim_class": "pre_blind_expectation_setting_only",
        "note": (
            "MOCK-BLIND expectation-setting diagnostic: every registered "
            "mapped control's protein scored with the FROZEN candidate-release "
            "bundle, ranking the published site among ALL Cys of its own "
            "protein (the task the four published workflows actually perform). "
            "Structure branch masked (release tomato semantics). NOT blind "
            "evidence; published positives are excluded from the SAP primary "
            "endpoint. Does not modify the model, features, Top-K, thresholds "
            "or SAP."
        ),
        "bundle": {"path": str(BUNDLE_PATH), "seed": bundle.seed},
        "controls": results,
        "summary": {
            "n_controls": len(results),
            "n_top1_k1": sum(
                1 for r in results if r["registered_metrics"]["top1"]
            ),
            "n_hit2_k1": sum(
                1 for r in results if r["registered_metrics"]["hit_at_2"]
            ),
            "burden_k1": {
                "observed": [
                    r["registered_metrics"]["first_hit_burden"] for r in results
                ],
                "random_baseline": [
                    r["registered_metrics"]["random_baseline_burden"] for r in results
                ],
            },
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    rows_path = args.output.with_suffix(".rows.tsv")
    flat_rows = []
    for r in results:
        reg = r["registered_metrics"]
        flat_rows.append(
            {
                "gene": r["gene"],
                "species": r["species"],
                "protein_accession": r["protein_accession"],
                "protein_length": r["protein_length"],
                "control_site": r["control_site"],
                "n_cys": reg["n_cys"],
                "true_site_rank": reg["true_site_ranks"][r["control_site"]],
                "top1": reg["top1"],
                "hit_at_2": reg["hit_at_2"],
                "hit_at_3": reg["hit_at_3"],
                "mrr": round(reg["mrr"], 4),
                "first_hit_burden": reg["first_hit_burden"],
                "random_baseline_burden": round(reg["random_baseline_burden"], 4),
                "ranking": ";".join(str(p) for p in reg["ranking"]),
            }
        )
    with rows_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(flat_rows[0].keys()),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(flat_rows)
    print(f"\nwrote {args.output}")
    print(f"wrote {rows_path}")


if __name__ == "__main__":
    main()
