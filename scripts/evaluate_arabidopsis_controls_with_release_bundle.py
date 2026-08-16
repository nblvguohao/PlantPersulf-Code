#!/usr/bin/env python
"""Mock-blind diagnostic for Arabidopsis: registered controls vs the release bundle.

Scores the 7 registered Arabidopsis known-mechanism controls (DES1, RBOHD,
SnRK2.6/OST1, ABI4, ATG4a, AtG6PD6, PAD3 — all from labs independent of the
Seville training network, all registered as 'mapped') with the FROZEN
candidate-release structure_ranker bundle and ranks them within the frozen-
bundle Arabidopsis candidate table (results/diagnostics/candidates_arabidopsis/).

Same discipline as the kiae271 tomato diagnostic: expectation-setting only,
no model/feature/Top-K/SAP change, results registered regardless of outcome.
The registered positions were verified against Arabidopsis proteome v1; this
script re-verifies them against the v2 proteome used by the v11 pipeline and
reports any coordinate mismatch instead of silently rescoring.
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

from scipy.stats import binomtest  # noqa: E402

from plantpersulf.evaluation.known_controls import REGISTERED_CONTROLS  # noqa: E402
from plantpersulf.features.sequence import _load_proteome  # noqa: E402
from plantpersulf.models.structure_ranker import (  # noqa: E402
    BranchFeatures,
    StructureRankerBundle,
    score_structure_ranker_bundle,
)
from plantpersulf.proteomics.multispecies_v2_dataset import MultispeciesV2SiteRow  # noqa: E402
from plantpersulf.proteomics.multispecies_v2_sources import sequence_feature_map  # noqa: E402

ARABIDOPSIS_PROTEOME_V2 = Path("data/raw/references/arabidopsis_ref_proteome_v2.fasta")
CANDIDATE_TABLE = (
    _REPO_ROOT / "results" / "diagnostics" / "candidates_arabidopsis" / "top_k_candidates.tsv"
)
BUNDLE_PATH = (
    _REPO_ROOT
    / "results"
    / "candidates"
    / "multispecies_v2_candidate_release_v1"
    / "model_weights"
    / "structure_ranker_bundle.pt"
)
DEFAULT_OUTPUT = (
    _REPO_ROOT / "results" / "known_controls" / "arabidopsis_release_bundle_recovery_v1.json"
)


def _percentile_rank(score: float, reference: list[float]) -> float:
    below = sum(1 for value in reference if value < score)
    return 100.0 * below / len(reference)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)

    if not CANDIDATE_TABLE.is_file():
        raise RuntimeError(
            f"arabidopsis candidate table not found: {CANDIDATE_TABLE}; "
            "run scripts/generate_tomato_topk_candidates.py --species arabidopsis first"
        )
    candidate_scores: dict[tuple[str, int], float] = {}
    with CANDIDATE_TABLE.open(encoding="utf-8", newline="") as handle:
        for record in csv.DictReader(handle, delimiter="\t"):
            candidate_scores[(record["site_key"], int(record["cys_position"]))] = float(
                record["score"]
            )
    reference = sorted(candidate_scores.values())
    print(f"arabidopsis candidate reference distribution: {len(reference)} sites")

    controls = [
        control
        for control in REGISTERED_CONTROLS
        if control.control_species == "Arabidopsis thaliana" and control.status == "mapped"
    ]
    print(f"registered arabidopsis mapped controls: {len(controls)}")

    proteome = _load_proteome(ARABIDOPSIS_PROTEOME_V2)

    # Re-verify coordinates against the v2 proteome used by the pipeline.
    verified: list[tuple] = []
    mismatched: list[dict[str, object]] = []
    for control in controls:
        sequence = proteome.get(control.uniprot_accession)
        position = control.cys_position
        if sequence is None:
            mismatched.append(
                {
                    "mechanism_lineage_id": control.mechanism_lineage_id,
                    "reason": "accession absent from proteome v2",
                }
            )
            continue
        if position < 1 or position > len(sequence) or sequence[position - 1] != "C":
            mismatched.append(
                {
                    "mechanism_lineage_id": control.mechanism_lineage_id,
                    "reason": f"position {position} is not Cys in proteome v2",
                }
            )
            continue
        verified.append((control, sequence))
    print(
        f"controls verified against proteome v2: {len(verified)}; "
        f"mismatched: {len(mismatched)}"
    )

    # Real structure features for the controls (same registry as the
    # candidate table, so scoring is apples-to-apples with the reference).
    from plantpersulf.evaluation.comparable_track import (
        build_registered_structure_features,
    )

    control_keys = {
        (f"arabidopsis|{c.uniprot_accession}", c.cys_position) for c, _ in verified
    }
    struct_values = build_registered_structure_features(
        control_keys,
        registry_path=_REPO_ROOT
        / "data"
        / "registry"
        / "releases"
        / "alphafold_structures_release_v2.tsv",
        registry_base=_REPO_ROOT / "data" / "registry",
    )

    rows: list[dict[str, object]] = []
    percentiles: list[float] = []
    for control, sequence in verified:
        global_id = f"arabidopsis|{control.uniprot_accession}"
        site_row = MultispeciesV2SiteRow(
            species="arabidopsis",
            protein_accession=control.uniprot_accession,
            cys_position=control.cys_position,
            label="unlabeled",
            study_accessions=(),
            global_protein_id=global_id,
            cluster_id="",
            split="development",
            development_fold=None,
        )
        features = sequence_feature_map([site_row], {"arabidopsis": proteome})
        key = (global_id, control.cys_position)
        values = features[key]
        if any(value is None for value in values):
            rows.append(
                {
                    "mechanism_lineage_id": control.mechanism_lineage_id,
                    "gene": control.gene,
                    "site_key": global_id,
                    "cys_position": control.cys_position,
                    "score": None,
                    "percentile_rank": None,
                    "recovery_status": "feature_failed",
                    "in_candidate_table": key in candidate_scores,
                }
            )
            continue
        branch = BranchFeatures(
            sequence=[list(float(v) for v in values)],
            esm=[[0.0]],
            structure=[list(struct_values[key][0])],
            structure_mask=[struct_values[key][1]],
            study_ids=None,
        )
        bundle = StructureRankerBundle.load(BUNDLE_PATH)
        output = score_structure_ranker_bundle(bundle, branch, device_name="cpu")
        score = output.scores[0]
        pct = _percentile_rank(score, reference)
        percentiles.append(pct)
        recovered = pct > 50
        rows.append(
            {
                "mechanism_lineage_id": control.mechanism_lineage_id,
                "gene": control.gene,
                "site_key": global_id,
                "cys_position": control.cys_position,
                "score": round(score, 8),
                "percentile_rank": round(pct, 2),
                "recovery_status": "recovered" if recovered else "not_recovered",
                "in_candidate_table": key in candidate_scores,
            }
        )

    n_scored = len(percentiles)
    n_recovered = sum(1 for r in rows if r["recovery_status"] == "recovered")
    mean_pct = sum(percentiles) / n_scored if percentiles else None
    median_pct = sorted(percentiles)[n_scored // 2] if percentiles else None
    two_sided_p = (
        binomtest(n_recovered, n_scored, 0.5, alternative="two-sided").pvalue
        if n_scored
        else None
    )
    print(
        f"\nrecovery: {n_recovered}/{n_scored} (pct>50), "
        f"mean percentile={mean_pct:.1f}, median={median_pct:.1f}, "
        f"two-sided binomial p={two_sided_p:.4f}"
    )
    for row in rows:
        print(
            f"  {row['gene']}: pct={row['percentile_rank']} "
            f"status={row['recovery_status']} in_candidates={row['in_candidate_table']}"
        )

    summary = {
        "track": "arabidopsis_controls_release_bundle_mock_blind",
        "claim_class": "pre_blind_expectation_setting_only",
        "note": (
            "Mock-blind diagnostic: 7 registered Arabidopsis controls "
            "(independent labs, never training positives in the multispecies "
            "v2 panel) scored with the FROZEN release bundle, percentiles "
            "within the frozen-bundle arabidopsis candidate table. "
            "Expectation-setting only; no model/feature/Top-K/SAP change."
        ),
        "controls": [
            {"mechanism_lineage_id": c.mechanism_lineage_id, "gene": c.gene}
            for c in controls
        ],
        "coordinate_mismatches": mismatched,
        "scoring": {
            "n_scored": n_scored,
            "n_recovered_pct_gt_50": n_recovered,
            "recovery_rate": (n_recovered / n_scored) if n_scored else None,
            "mean_percentile": mean_pct,
            "median_percentile": median_pct,
            "two_sided_binomial_p_vs_50": two_sided_p,
            "candidate_reference_n": len(reference),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    rows_path = args.output.with_suffix(".rows.tsv")
    with rows_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {args.output}")
    print(f"wrote {rows_path}")


if __name__ == "__main__":
    main()
