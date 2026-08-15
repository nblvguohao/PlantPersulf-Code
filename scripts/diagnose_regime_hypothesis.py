#!/usr/bin/env python
"""Diagnostic-only: the cys-density regime hypothesis, on frozen scores.

The 2026-08-14 increment proposal (PlantPersulf_增量建议 §4) states a
falsifiable mechanism for the model's flat tomato ranking: the SAME
``local_cys_density`` feature would need opposite weights in the two site
regimes —

* IDR-type true sites (SlWRKY6 C396, PyMYB10 C218) sit in Cys-poor
  disordered C-termini: the feature must act as a *negative* signal;
* metal-cluster true sites (BRG3 C206/C212) sit in a Cys-dense RING cluster:
  the feature must act as a *positive* signal.

A single global additive model must compromise between the two. This script
tests the data side of the hypothesis on the registered controls:

1. what the FROZEN bundle's actual sequence features are (hydrophobicity,
   protein-level cys_density, local positive-charge density — note that the
   model's ``cys_density`` is protein-level and therefore CONSTANT within a
   protein, so it cannot drive within-protein ranking at all);
2. whether the diagnostic ``local_cys_density`` (C count in the +/-10
   window) separates true sites from decoys in the predicted directions;
3. whether the frozen model's within-protein ranking correlates with
   ``local_cys_density``.

claim class: diagnostic_only. Changes no model, feature, Top-K, threshold
or SAP.
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
    make_site_row,
    score_sites,
)
from plantpersulf.features.sequence import (  # noqa: E402
    _flanking_window,
    _load_proteome,
)
from plantpersulf.models.structure_ranker import StructureRankerBundle  # noqa: E402
from plantpersulf.proteomics.multispecies_v2_sources import (  # noqa: E402
    sequence_feature_map,
)

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
DEFAULT_OUTPUT = _REPO_ROOT / "results" / "diagnostics" / "regime_hypothesis_v1.json"

SPECIES_PROTEOMES = {
    "Solanum lycopersicum": ("tomato", TOMATO_PROTEOME),
    "Arabidopsis thaliana": ("arabidopsis", ARABIDOPSIS_PROTEOME),
}

# Canonical regime labels for the headline cases (methodology design §3,
# verified in the increment proposal appendix):
# IDR-type: true site in a Cys-poor disordered region, decoys in a Cys-dense
#   metal/DBD cluster — hypothesis: local_cys_density(true) < local density(decoys)
# metal-cluster-type: true sites inside a Cys-dense cluster — hypothesis:
#   local_cys_density(true) ~= local density of its own cluster
REGIME_EXPECTATION: dict[str, dict[str, object]] = {
    "SLWRKY6_H2S_PHOSPHORYLATION": {
        "regime": "IDR",
        "expected": "true_below_decoys",
        "note": (
            "C396 in C-terminal IDR; five decoy Cys in the WRKY zinc-finger cluster"
        ),
    },
    "BRG3_H2S_UBIQUITINATION": {
        "regime": "metal_cluster",
        "expected": "true_within_cluster",
        "note": (
            "C206/C212 inside the RING cluster (9 Cys in 35 aa); "
            "C209 is the co-peptide negative"
        ),
    },
}


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)

    controls = [c for c in REGISTERED_CONTROLS if c.status == "mapped"]
    bundle = StructureRankerBundle.load(BUNDLE_PATH)
    proteomes: dict[str, dict[str, str]] = {}

    per_protein_rows: list[dict[str, object]] = []
    regime_results: list[dict[str, object]] = []

    for control in controls:
        species, proteome_path = SPECIES_PROTEOMES[control.control_species]
        proteomes.setdefault(species, _load_proteome(proteome_path))
        sequence = proteomes[species][control.uniprot_accession]
        positions = cys_positions(sequence)
        sites = [(species, control.uniprot_accession, p) for p in positions]
        scored = score_sites(sites, bundle=bundle, proteomes=proteomes)

        rows = [make_site_row(*site) for site in sites]
        raw = sequence_feature_map(rows, proteomes)

        cys_rows: list[dict[str, object]] = []
        for position in positions:
            flank = _flanking_window(sequence, position, 10)
            local_cys_density = flank.count("C") / len(flank) if flank else 0.0
            model_features = raw[(f"{species}|{control.uniprot_accession}", position)]
            cys_rows.append(
                {
                    "cys_position": position,
                    "is_true_site": position == control.cys_position,
                    "hydrophobicity": round(float(model_features[0]), 4),
                    "protein_cys_density": round(float(model_features[1]), 4),
                    "local_positive_charge_density": round(float(model_features[2]), 4),
                    "local_cys_density_diag": round(local_cys_density, 4),
                    "score": round(
                        float(
                            scored[
                                (species, control.uniprot_accession, position)
                            ].score
                        ),
                        6,
                    ),
                }
            )

        true_rows = [r for r in cys_rows if r["is_true_site"]]
        decoy_rows = [r for r in cys_rows if not r["is_true_site"]]
        mean_true = sum(r["local_cys_density_diag"] for r in true_rows) / len(true_rows)
        mean_decoy = (
            sum(r["local_cys_density_diag"] for r in decoy_rows) / len(decoy_rows)
        )

        expectation = REGIME_EXPECTATION.get(control.mechanism_lineage_id)
        observed = (
            "true_below_decoys"
            if mean_true < mean_decoy
            else "true_at_least_decoys"
        )
        matched = None
        if expectation is not None:
            expected = expectation["expected"]
            if expected == "true_below_decoys":
                matched = observed == "true_below_decoys"
            else:  # metal cluster: true sites inside the cluster, comparable density
                matched = observed == "true_at_least_decoys"

        regime_results.append(
            {
                "mechanism_lineage_id": control.mechanism_lineage_id,
                "gene": control.gene,
                "regime": expectation["regime"] if expectation else "unclassified",
                "expected": expectation["expected"] if expectation else None,
                "observed_direction": observed,
                "hypothesis_matched": matched,
                "mean_local_cys_density_true": round(mean_true, 4),
                "mean_local_cys_density_decoy": round(mean_decoy, 4),
            }
        )
        per_protein_rows.extend(
            {
                "gene": control.gene,
                "species": species,
                "protein_accession": control.uniprot_accession,
                **r,
            }
            for r in cys_rows
        )

        tag = expectation["regime"] if expectation else "---"
        print(
            f"  {control.gene} [{tag}] Cys{control.cys_position}: "
            f"local_density true={mean_true:.3f} decoys={mean_decoy:.3f} "
            f"-> {observed}"
        )

    summary = {
        "track": "regime_hypothesis_v1",
        "claim_class": "diagnostic_only",
        "note": (
            "Frozen-bundle feature/score export for the cys-density regime "
            "hypothesis (increment proposal 2026-08-14 §4). KEY MODEL FACT: the "
            "frozen bundle's only density feature is PROTEIN-level cys_density "
            "(count(C)/protein length), constant within a protein, so it cannot "
            "drive within-protein ranking; the local-window density that the "
            "hypothesis concerns is NOT in the model. This diagnostic computes "
            "it separately. Does not modify the model, features, Top-K, "
            "thresholds or SAP."
        ),
        "bundle": {"path": str(BUNDLE_PATH), "seed": bundle.seed},
        "model_sequence_features": [
            "hydrophobicity",
            "protein_cys_density (constant within protein)",
            "local_positive_charge_density",
        ],
        "regime_checks": regime_results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    rows_path = args.output.with_suffix(".rows.tsv")
    with rows_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(per_protein_rows[0].keys()),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(per_protein_rows)
    print(f"\nwrote {args.output}")
    print(f"wrote {rows_path}")


if __name__ == "__main__":
    main()
