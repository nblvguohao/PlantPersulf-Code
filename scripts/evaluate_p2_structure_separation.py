#!/usr/bin/env python
"""P2 early gate: can structure-derived features separate co-peptide Cys?

The methodology design's P2 gate ("three days, decides whether the project
is worth doing") asks whether L1-style structure features can rank
SlBRG3 C206/C212 ABOVE the gold-standard co-peptide negative C209 — a task
that is information-theoretically impossible for +/-10 sequence windows
(empirically confirmed on the frozen bundle: 8/9 registered peptides
unseparated).

This diagnostic computes the structure feature set
(``evaluation.structure_features``) for every Cys of the AFDB-covered
control proteins and reports three evidence levels:

1. per-feature separation direction on BRG3 (true {206,212} vs negative
   {209}) — is there ANY structure feature with a consistent direction?
2. per-feature within-protein ranking of the true sites (Top-1 / Hit@2 /
   first-hit burden vs random baseline) — which features carry signal;
3. a z-score composite ranking (all features, equal weight) — the
   "hand-built features + ranking" analogue of the P2 plan, plus a note
   on why supervised gradient boosting is not estimable at n=4 labelled
   sites (the honest sample-size finding of the gate itself).

PyMYB10 (A0A0U2QCQ1) has NO AlphaFold DB file (API entry exists, all model
versions 404) — folding (ESMFold/Boltz-2) is a separate decision; the
script records the gap instead of inventing a structure.

Diagnostic-only: no fitting of frozen artifacts, no SAP changes.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from plantpersulf.evaluation.structure_features import (
    STRUCTURE_FEATURE_NAMES,
    cys_structure_features,
    parse_pdb,
)
from plantpersulf.evaluation.within_protein_ranking import (
    first_hit_rank,
    hit_at_k,
    mutagenesis_burden,
    rank_sites_desc,
)
from plantpersulf.features.sequence import _load_proteome

_REPO_ROOT = Path(__file__).resolve().parents[1]
ALPHAFOLD_DIR = _REPO_ROOT / "data" / "raw" / "alphafold"
PROTEOMES = {
    "tomato": (
        _REPO_ROOT / "data" / "raw" / "references" / "tomato_ref_proteome_v1.fasta"
    ),
}
OUTPUT = _REPO_ROOT / "results" / "diagnostics" / "p2_structure_separation_v1.json"

# (species, accession, pdb file, true sites, negative sites)
PROTEINS = [
    (
        "tomato",
        "A0A3Q7EW23",  # SlBRG3 — RING domain, co-peptide gate case
        "AF-A0A3Q7EW23-F1-model.pdb",
        {206, 212},
        {209},
    ),
    (
        "tomato",
        "A0A3Q7F586",  # SlWRKY6 — C-terminal IDR case
        "AF-A0A3Q7F586-F1-model.pdb",
        {396},
        set(),
    ),
]
PYMYB10 = ("tomato", "A0A0U2QCQ1")  # no AFDB file (recorded, not fabricated)


def _ranking_summary(
    scores: dict[int, float], true_positions: set[int], n_cys: int
) -> dict[str, object]:
    ranking = rank_sites_desc(scores)
    ranks = {position: rank for rank, position in enumerate(ranking, start=1)}
    k = len(true_positions)
    return {
        "ranking": ranking,
        "true_ranks": {p: ranks[p] for p in sorted(true_positions)},
        "top1": ranks[min(true_positions)] == 1,
        "hit_at_2": hit_at_k(ranking, true_positions, k=2),
        "first_hit_burden": first_hit_rank(ranking, true_positions),
        "random_baseline_burden": mutagenesis_burden(n_cys, k),
    }


def main() -> None:
    proteomes = {"tomato": _load_proteome(PROTEOMES["tomato"])}

    per_protein: dict[str, dict[str, object]] = {}
    brg3_features: dict[int, dict[str, float]] | None = None
    brg3_true = set()
    brg3_negative = set()

    for species, accession, pdb_file, true_sites, negative_sites in PROTEINS:
        pdb_path = ALPHAFOLD_DIR / pdb_file
        if not pdb_path.exists():
            per_protein[accession] = {"skipped": f"missing {pdb_file}"}
            continue
        residues = parse_pdb(pdb_path.read_text(encoding="utf-8"))
        sequence = proteomes[species][accession]
        if len(residues) != len(sequence):
            raise RuntimeError(
                f"{accession}: PDB {len(residues)} residues != proteome "
                f"{len(sequence)}"
            )
        mismatches = [
            position
            for position in range(1, len(sequence) + 1)
            if residues[position - 1].type != sequence[position - 1]
        ]
        if mismatches:
            raise RuntimeError(
                f"{accession}: PDB/proteome mismatch at {mismatches[:10]}"
            )
        positions = [
            i + 1 for i, residue in enumerate(sequence) if residue == "C"
        ]
        features = cys_structure_features(residues, positions)
        per_protein[accession] = {
            "species": species,
            "pdb": pdb_file,
            "n_cys": len(positions),
            "true_sites": sorted(true_sites),
            "negative_sites": sorted(negative_sites),
            "features": features,
        }
        if accession == "A0A3Q7EW23":
            brg3_features = features
            brg3_true = true_sites
            brg3_negative = negative_sites

    assert brg3_features is not None
    assert brg3_true and brg3_negative

    # --- level 1: per-feature separation direction on BRG3 -----------------
    separation: dict[str, dict[str, object]] = {}
    for feature in STRUCTURE_FEATURE_NAMES:
        true_values = sorted(brg3_features[p][feature] for p in brg3_true)
        negative_values = sorted(brg3_features[p][feature] for p in brg3_negative)
        min_true, max_true = true_values[0], true_values[-1]
        separation[feature] = {
            "true_values": true_values,
            "negative_values": negative_values,
            "all_true_above_all_negative": min_true > negative_values[-1],
            "all_true_below_all_negative": max_true < negative_values[0],
            "direction": (
                "true_above"
                if min_true > negative_values[-1]
                else "true_below"
                if max_true < negative_values[0]
                else "no_clean_separation"
            ),
        }

    # --- level 2: single-feature within-protein ranking --------------------
    ranking_by_feature: dict[str, dict[str, object]] = {}
    for accession, record in per_protein.items():
        if "features" not in record:
            continue
        n_cys = record["n_cys"]
        true_sites = set(record["true_sites"])
        if not true_sites:
            continue
        for feature in STRUCTURE_FEATURE_NAMES:
            scores = {
                position: float(row[feature])
                for position, row in record["features"].items()
            }
            key = f"{accession}:{feature}"
            ranking_by_feature[key] = _ranking_summary(
                scores, true_sites, n_cys
            )

    # --- level 3: z-score composite ranking ---------------------------------
    composite: dict[str, dict[str, object]] = {}
    for accession, record in per_protein.items():
        if "features" not in record:
            continue
        true_sites = set(record["true_sites"])
        if not true_sites:
            continue
        n_cys = record["n_cys"]
        feature_rows = record["features"]
        z_by_feature: dict[str, dict[int, float]] = {}
        for feature in STRUCTURE_FEATURE_NAMES:
            values = [feature_rows[p][feature] for p in feature_rows]
            mean = float(np.mean(values))
            std = float(np.std(values)) or 1.0
            z_by_feature[feature] = {
                position: (feature_rows[position][feature] - mean) / std
                for position in feature_rows
            }
        scores = {
            position: float(
                np.mean(
                    [
                        z_by_feature[feature][position]
                        for feature in STRUCTURE_FEATURE_NAMES
                    ]
                )
            )
            for position in feature_rows
        }
        composite[accession] = _ranking_summary(scores, true_sites, n_cys)

    # --- level 4: signed upper bound (BRG3, overfit by construction) -------
    # Directions taken FROM the BRG3 separation table itself: this is the
    # ceiling of what the feature set could achieve with perfect direction
    # learning, NOT an unbiased estimate. Only features with a clean
    # direction participate; others contribute nothing.
    signed: dict[str, float] = {}
    for feature in STRUCTURE_FEATURE_NAMES:
        direction = separation[feature]["direction"]
        if direction == "true_above":
            signed[feature] = 1.0
        elif direction == "true_below":
            signed[feature] = -1.0
    signed_upper_bound: dict[str, object] = {}
    if signed:
        brg3_record = per_protein["A0A3Q7EW23"]
        feature_rows = brg3_record["features"]
        n_cys = brg3_record["n_cys"]
        signed_z: dict[str, dict[int, float]] = {}
        for feature in signed:
            values = [feature_rows[p][feature] for p in feature_rows]
            mean = float(np.mean(values))
            std = float(np.std(values)) or 1.0
            signed_z[feature] = {
                position: (feature_rows[position][feature] - mean) / std
                for position in feature_rows
            }
        scores = {
            position: float(
                np.mean(
                    [signed[f] * signed_z[f][position] for f in signed]
                )
            )
            for position in feature_rows
        }
        signed_upper_bound = _ranking_summary(
            scores, brg3_true, n_cys
        )
    signed_note = (
        "BRG3 signed upper bound: z-scores combined with signs taken from "
        "the BRG3 separation table itself (overfit ceiling, NOT an unbiased "
        "estimate). It quantifies how much ranking power the structure "
        "features could deliver if a direction learner were available — "
        "which requires MIL-scale labels (methodology L0), not n=4."
    )

    # --- level 5: regime-local signed ranking (BRG3 RING domain) -----------
    # Directional features are regime-local: the N-terminal disordered Cys
    # (C24/C30/C173/C184/C185, pLDDT 42-53) have the same "exposed" shape as
    # a true site but are not modified, so a protein-wide signed composite
    # loses the true sites to the disordered tail. Within the RING domain
    # (C197-C231, the methodology's metal-cluster regime) the same signs
    # should push the true sites up. This is the regime-routing hypothesis
    # (proposition 3) tested with hand features instead of a MoE gate.
    ring_domain: dict[str, object] = {}
    if signed:
        brg3_record = per_protein["A0A3Q7EW23"]
        feature_rows = brg3_record["features"]
        domain_positions = [
            p for p in feature_rows if 197 <= p <= 231
        ]
        ring_signed_z: dict[str, dict[int, float]] = {}
        for feature in signed:
            values = [feature_rows[p][feature] for p in domain_positions]
            mean = float(np.mean(values))
            std = float(np.std(values)) or 1.0
            ring_signed_z[feature] = {
                p: (feature_rows[p][feature] - mean) / std
                for p in domain_positions
            }
        domain_scores = {
            p: float(np.mean([signed[f] * ring_signed_z[f][p] for f in signed]))
            for p in domain_positions
        }
        ring_domain = _ranking_summary(domain_scores, brg3_true, len(domain_positions))
    ring_note = (
        "Regime-local signed ranking: identical signs to level 4, but "
        "z-normalised and ranked WITHIN the RING domain (C197-C231, "
        "metal-cluster regime). The regime prior comes from the methodology "
        "design's RING annotation, not from the data. Compares directly with "
        "level 4's protein-wide signed composite."
    )

    # --- tree note: supervised boosting is not estimable at n=4 ------------
    tree_note = (
        "Supervised gradient boosting is NOT estimable on this gate: only "
        "4 labelled sites exist (BRG3 true 206/212 + negative 209, WRKY6 "
        "true 396). A tree on 21 Cys with 4 labels cannot generalise; the "
        "honest gate finding is that structure evidence must be evaluated "
        "unsupervised (feature separation + ranking), which is what levels "
        "1-3 report. MIL-scale labelled data (methodology L0) is the "
        "prerequisite for the P2 tree plan as written."
    )

    summary = {
        "track": "p2_structure_separation_v1",
        "claim_class": "diagnostic_only",
        "note": (
            "P2 early gate: structure features (pLDDT, RSA, S-gamma cluster "
            "geometry, positive-residue microenvironment, Coulomb potential, "
            "contact number) computed from registered AlphaFold DB models; "
            "no fitting of frozen artifacts; metal coordination is a "
            "geometry proxy because AFDB models contain no metal ions; "
            "Coulomb potential is a simplified dielectric-scaled sum (APBS "
            "PB is the rigorous follow-up)."
        ),
        "structures": {
            "A0A3Q7EW23": "AFDB v6 (registered alphafold_structures.tsv)",
            "A0A3Q7F586": "AFDB v6 (registered alphafold_structures.tsv)",
            "A0A0U2QCQ1": (
                "NO AFDB file (API entry exists, model files 404); ESMFold/"
                "Boltz-2 folding is a separate decision — gap recorded, not "
                "fabricated"
            ),
        },
        "per_protein": per_protein,
        "brg3_per_feature_separation": separation,
        "ranking_by_feature": ranking_by_feature,
        "composite_zscore_ranking": composite,
        "signed_upper_bound_brg3": signed_upper_bound,
        "signed_upper_bound_note": signed_note,
        "ring_domain_signed_ranking": ring_domain,
        "ring_domain_note": ring_note,
        "supervised_tree_note": tree_note,
        "summary": {
            "n_proteins_with_structure": sum(
                1 for r in per_protein.values() if "features" in r
            ),
            "n_features": len(STRUCTURE_FEATURE_NAMES),
            "brg3_features_with_clean_direction": sum(
                1
                for s in separation.values()
                if s["direction"] != "no_clean_separation"
            ),
            "composite_top1": composite.get("A0A3Q7EW23", {}).get("top1"),
            "composite_first_hit_burden": composite.get("A0A3Q7EW23", {}).get(
                "first_hit_burden"
            ),
            "signed_upper_bound_burden": signed_upper_bound.get(
                "first_hit_burden"
            ),
            "ring_domain_signed_burden": ring_domain.get(
                "first_hit_burden"
            ),
            "ring_domain_signed_top1": ring_domain.get("top1"),
        },
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    # console digest
    print("BRG3 per-feature separation:")
    for feature, s in separation.items():
        print(
            f"  {feature:<28} true={s['true_values']} neg={s['negative_values']} "
            f"-> {s['direction']}"
        )
    print(f"\nBRG3 z-score composite: {composite['A0A3Q7EW23']}")
    print(f"WRKY6 z-score composite: {composite.get('A0A3Q7F586', {})}")
    if signed_upper_bound:
        print(
            f"BRG3 signed upper bound (whole protein): "
            f"{signed_upper_bound['first_hit_burden']} "
            f"(random {signed_upper_bound['random_baseline_burden']}) "
            f"true_ranks={signed_upper_bound['true_ranks']}"
        )
    if ring_domain:
        print(
            f"BRG3 signed ranking (RING domain only): "
            f"{ring_domain['first_hit_burden']} "
            f"(random {ring_domain['random_baseline_burden']}) "
            f"true_ranks={ring_domain['true_ranks']} "
            f"ranking={ring_domain['ranking']}"
        )
    print(f"\nwrote {OUTPUT}")


if __name__ == "__main__":
    main()
