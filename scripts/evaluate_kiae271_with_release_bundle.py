#!/usr/bin/env python
"""Pre-blind expectation-setting diagnostic: kiae271 sites vs the release bundle.

The 99 coordinate-verified tomato persulfidation sites published in kiae271
(Zhang et al. 2024, doi:10.1093/plphys/kiae271) are scored with the FROZEN
candidate-release structure_ranker bundle and their percentile ranks are
computed within the frozen candidate table (top_k_candidates.tsv, 179,736
tomato sites).

This is a MOCK-BLIND DIAGNOSTIC for expectation-setting ONLY. It:

- uses the release bundle byte-for-byte (no refit, no feature change);
- does NOT change the model, features, Top-K, thresholds or SAP (constraint 1);
- records its result in the negative-results register regardless of outcome;
- must never be reported as blind-validation evidence (kiae271 sites are
  published positives; the SAP excludes them from the primary endpoint).

Also reports the overlap between the kiae271 sites and the release training
panel: a site present in the candidate table is by construction NOT in the
training panel (candidates exclude panel sites), so its absence from the
candidate table flags training contamination of that site.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from scipy.stats import binomtest  # noqa: E402

from plantpersulf.features.sequence import _load_proteome  # noqa: E402
from plantpersulf.models.structure_ranker import (  # noqa: E402
    BranchFeatures,
    StructureRankerBundle,
    score_structure_ranker_bundle,
)
from plantpersulf.proteomics.kiae271_sites import parse_kiae271_sites  # noqa: E402
from plantpersulf.proteomics.multispecies_v2_dataset import MultispeciesV2SiteRow  # noqa: E402
from plantpersulf.proteomics.multispecies_v2_sources import sequence_feature_map  # noqa: E402

KIAE271_XLSX = Path("data/raw/supplements/KIAE271_SUPPL/kiae271_DSs.xlsx")
TOMATO_PROTEOME = Path("data/raw/references/tomato_ref_proteome_v1.fasta")
RELEASE_DIR = (
    _REPO_ROOT / "results" / "candidates" / "multispecies_v2_candidate_release_v1"
)
BUNDLE_PATH = RELEASE_DIR / "model_weights" / "structure_ranker_bundle.pt"
CANDIDATE_TABLE = RELEASE_DIR / "top_k_candidates.tsv"
DEFAULT_OUTPUT = (
    _REPO_ROOT / "results" / "known_controls" / "kiae271_release_bundle_recovery_v1.json"
)


def _percentile_rank(score: float, reference: list[float]) -> float:
    """Fraction of reference scores strictly below ``score``, in percent."""
    below = sum(1 for value in reference if value < score)
    return 100.0 * below / len(reference)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--localization-min", type=float, default=0.75)
    args = parser.parse_args(argv)

    # --- candidate-table reference distribution ------------------------------
    candidate_scores: dict[tuple[str, int], float] = {}
    with CANDIDATE_TABLE.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for record in reader:
            candidate_scores[(record["site_key"], int(record["cys_position"]))] = float(
                record["score"]
            )
    reference = sorted(candidate_scores.values())
    print(f"candidate reference distribution: {len(reference)} sites")

    # --- kiae271 verified sites ----------------------------------------------
    tomato_proteome = _load_proteome(TOMATO_PROTEOME)
    table = parse_kiae271_sites(
        KIAE271_XLSX, tomato_proteome, localization_min=args.localization_min
    )
    print(
        f"kiae271: {table.rows_total} rows -> {table.total_verified_sites} "
        f"verified sites"
    )

    # --- overlap vs training panel (via candidate-table membership) ----------
    in_candidates = [
        site
        for site in table.sites
        if (f"tomato|{site.protein_accession}", site.cys_position) in candidate_scores
    ]
    in_panel = [
        site
        for site in table.sites
        if (f"tomato|{site.protein_accession}", site.cys_position) not in candidate_scores
    ]
    print(
        f"kiae271 sites in candidate table (NOT in training panel): {len(in_candidates)}; "
        f"absent from candidate table (IN training panel): {len(in_panel)}"
    )

    # --- features for all kiae271 sites (sequence only; structure masked) ----
    scan_rows = [
        MultispeciesV2SiteRow(
            species="tomato",
            protein_accession=site.protein_accession,
            cys_position=site.cys_position,
            label="unlabeled",
            study_accessions=(),
            global_protein_id=f"tomato|{site.protein_accession}",
            cluster_id="",
            split="development",
            development_fold=None,
        )
        for site in table.sites
    ]
    raw_features = sequence_feature_map(scan_rows, {"tomato": tomato_proteome})
    sequence_features = {
        key: tuple(float(value) for value in values if value is not None)
        for key, values in raw_features.items()
    }
    ordered_keys = sorted(sequence_features)
    train = BranchFeatures(
        sequence=[list(sequence_features[key]) for key in ordered_keys],
        esm=[[0.0] for _ in ordered_keys],
        structure=[[0.0, 0.0] for _ in ordered_keys],
        structure_mask=[False for _ in ordered_keys],
        study_ids=None,
    )
    bundle = StructureRankerBundle.load(BUNDLE_PATH)
    output = score_structure_ranker_bundle(bundle, train, device_name="cpu")

    # --- assemble rows ---------------------------------------------------------
    rows: list[dict[str, object]] = []
    percentiles: list[float] = []
    n_recovered = 0
    for site, key, score in zip(table.sites, ordered_keys, output.scores, strict=True):
        pct = _percentile_rank(score, reference)
        recovered = pct > 50
        percentiles.append(pct)
        if recovered:
            n_recovered += 1
        rows.append(
            {
                "protein_accession": site.protein_accession,
                "cys_position": site.cys_position,
                "regulation": site.regulation,
                "localization_prob": site.localization_prob,
                "score": round(float(score), 8),
                "percentile_rank": round(pct, 2),
                "recovery_status": "recovered" if recovered else "not_recovered",
                "in_candidate_table": (
                    f"tomato|{site.protein_accession}", site.cys_position
                ) in candidate_scores,
            }
        )

    n_scored = len(percentiles)
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

    summary = {
        "track": "kiae271_release_bundle_mock_blind_diagnostic",
        "claim_class": "pre_blind_expectation_setting_only",
        "note": (
            "MOCK-BLIND diagnostic for expectation-setting: the 99 published "
            "kiae271 tomato sites scored with the FROZEN candidate-release "
            "bundle, percentiles within the frozen candidate table. NOT blind "
            "evidence; kiae271 sites are published positives and are excluded "
            "from the SAP primary endpoint. Does not modify the model, "
            "features, Top-K, thresholds or SAP."
        ),
        "bundle": {
            "path": "model_weights/structure_ranker_bundle.pt",
            "seed": bundle.seed,
        },
        "source": {
            "doi": "10.1093/plphys/kiae271",
            "supplementary_file": str(KIAE271_XLSX),
            "localization_min": args.localization_min,
        },
        "parsed": {
            "rows_total": table.rows_total,
            "total_verified_sites": table.total_verified_sites,
            "regulation_breakdown": dict(Counter(s.regulation for s in table.sites)),
        },
        "overlap": {
            "n_in_candidate_table": len(in_candidates),
            "n_absent_from_candidate_table": len(in_panel),
            "interpretation": (
                "sites in the candidate table are by construction NOT in the "
                "release training panel; absent sites are training-panel sites"
            ),
        },
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
