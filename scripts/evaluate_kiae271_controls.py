#!/usr/bin/env python
"""Gate 2 condition 3 (recovery_no_leakage) expansion: kiae271 differential
persulfidome as a large cross-species control panel.

Extends the retrospective known-control recovery evaluation
(``scripts/evaluate_known_controls.py``) from the ~9 hand-curated controls in
``known_controls.py`` to the ~99 coordinate-verified sites parsed from
kiae271's own published Supplementary Dataset S1 (see
``plantpersulf.proteomics.kiae271_sites``) — a genuinely independent lab,
species, and chemistry relative to the Arabidopsis Seville tag-switch
training benchmark. This does NOT flip Gate 2 condition 1 (independent
studies within the same benchmark; cross-species data structurally cannot,
same conclusion as the rice/fungal transfer tracks) — it only increases the
statistical power of condition 3 (recovery rate) by roughly an order of
magnitude, using the exact same shared-PU-scorer methodology as
``evaluate_known_controls.py`` so the two are directly comparable.

Usage::

    python scripts/evaluate_kiae271_controls.py \
        --output results/known_controls/kiae271_differential_recovery_v1.json
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any

from evaluate_known_controls import (  # type: ignore[import-not-found]
    _fit_pu_scorer,
    _percentile_rank,
    _score_control,
)

from plantpersulf.features.sequence import _load_proteome
from plantpersulf.proteomics.kiae271_sites import parse_kiae271_sites

KIAE271_XLSX = Path("data/raw/supplements/KIAE271_SUPPL/kiae271_DSs.xlsx")
ARABIDOPSIS_BENCHMARK = Path("data/processed/benchmark_v1/sites.tsv")
ARABIDOPSIS_PROTEOME = Path("data/raw/references/arabidopsis_ref_proteome_v1.fasta")
TOMATO_PROTEOME = Path("data/raw/references/tomato_ref_proteome_v1.fasta")


def run_kiae271_control_evaluation(
    output_path: Path | None,
    xlsx_path: Path = KIAE271_XLSX,
    benchmark_path: Path = ARABIDOPSIS_BENCHMARK,
    proteome_path: Path = ARABIDOPSIS_PROTEOME,
    tomato_proteome_path: Path = TOMATO_PROTEOME,
    localization_min: float = 0.75,
) -> dict[str, Any]:
    for path in (xlsx_path, benchmark_path, proteome_path, tomato_proteome_path):
        if not path.is_file():
            raise RuntimeError(f"required input missing: {path}")

    tomato_proteome = _load_proteome(tomato_proteome_path)
    table = parse_kiae271_sites(
        xlsx_path, tomato_proteome, localization_min=localization_min
    )
    print(
        f"kiae271 Dataset S1: {table.rows_total} rows -> "
        f"{table.total_verified_sites} verified sites / "
        f"{table.total_verified_proteins} proteins"
    )

    # Same shared PU scorer as evaluate_known_controls.py, trained ONCE on
    # the Arabidopsis benchmark — these are cross-species controls, so
    # nothing here is excluded from the training pool (unlike same-species
    # Arabidopsis controls, which must be held out of it).
    scorer, unlabeled_X = _fit_pu_scorer(benchmark_path, proteome_path)
    unlabeled_scores = scorer(unlabeled_X)
    print(f"Unlabeled reference distribution: {len(unlabeled_scores)} samples")

    rows: list[dict[str, Any]] = []
    for site in table.sites:
        score, detail = _score_control(
            scorer, site.protein_accession, site.cys_position, tomato_proteome_path
        )
        pct = _percentile_rank(score, unlabeled_scores) if detail == "ok" else None
        rows.append(
            {
                "protein_accession": site.protein_accession,
                "cys_position": site.cys_position,
                "regulation": site.regulation,
                "localization_prob": site.localization_prob,
                "intensity_lcd": site.intensity_lcd,
                "intensity_wt": site.intensity_wt,
                "score": score if detail == "ok" else None,
                "percentile_rank": round(pct, 1) if pct is not None else None,
                "recovery_status": (
                    "recovered" if pct is not None and pct > 50 else "not_recovered"
                ),
                "detail": detail,
            }
        )

    percentiles = [
        r["percentile_rank"] for r in rows if r["percentile_rank"] is not None
    ]
    n_recovered = sum(1 for r in rows if r["recovery_status"] == "recovered")
    n_scored = len(percentiles)

    def _stratum_stats(regulation: str) -> dict[str, Any]:
        vals = [
            r["percentile_rank"]
            for r in rows
            if r["regulation"] == regulation and r["percentile_rank"] is not None
        ]
        n_rec = sum(
            1
            for r in rows
            if r["regulation"] == regulation and r["recovery_status"] == "recovered"
        )
        return {
            "n": sum(1 for r in rows if r["regulation"] == regulation),
            "n_scored": len(vals),
            "n_recovered": n_rec,
            "recovery_rate": (n_rec / len(vals)) if vals else None,
            "mean_percentile": (sum(vals) / len(vals)) if vals else None,
            "median_percentile": (sorted(vals)[len(vals) // 2]) if vals else None,
        }

    summary: dict[str, Any] = {
        "track": "kiae271_differential_control_panel",
        "claim_class": "gate2_condition3_recovery_expansion_only",
        "note": (
            "Cross-species (tomato) data; does NOT satisfy Gate 2 condition 1 "
            "(independent studies within the same benchmark) — only expands "
            "condition 3 (recovery_no_leakage) statistical power, same "
            "conclusion already established for the rice/fungal transfer "
            "tracks. Comparable in methodology to "
            "results/known_controls/recovery_v1.json (same shared PU scorer, "
            "same percentile-rank definition)."
        ),
        "source": {
            "doi": "10.1093/plphys/kiae271",
            "supplementary_file": str(xlsx_path),
            "registry_row": "data/registry/supplementary_sources.tsv (KIAE271_SUPPL)",
            "localization_min": localization_min,
        },
        "parsed": {
            "rows_total": table.rows_total,
            "dropped_not_cysteine": table.dropped_not_cysteine,
            "dropped_no_position_alignment": table.dropped_no_position_alignment,
            "dropped_missing_accession": table.dropped_missing_accession,
            "dropped_coordinate_mismatch": table.dropped_coordinate_mismatch,
            "dropped_low_localization": table.dropped_low_localization,
            "dropped_no_quantified_intensity": table.dropped_no_quantified_intensity,
            "dropped_duplicate": table.dropped_duplicate,
            "total_verified_sites": table.total_verified_sites,
            "total_verified_proteins": table.total_verified_proteins,
            "regulation_breakdown": dict(Counter(s.regulation for s in table.sites)),
        },
        "scoring": {
            "n_scored": n_scored,
            "n_feature_extraction_failed": sum(1 for r in rows if r["detail"] != "ok"),
            "n_recovered_pct_gt_50": n_recovered,
            "recovery_rate": (n_recovered / n_scored) if n_scored else None,
            "mean_percentile": (sum(percentiles) / len(percentiles))
            if percentiles
            else None,
            "median_percentile": sorted(percentiles)[len(percentiles) // 2]
            if percentiles
            else None,
        },
        "scoring_by_regulation": {
            REGULATION: _stratum_stats(REGULATION)
            for REGULATION in ("lcd_gain", "wt_only", "both")
        },
    }

    if output_path:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(
            json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        rows_path = output_path.with_suffix(".rows.tsv")
        with rows_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=list(rows[0].keys()) if rows else [],
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)
        print(f"\nSummary written to {output_path}")
        print(f"Per-row detail written to {rows_path}")

    print(
        f"\nOverall: {n_recovered}/{n_scored} recovered (pct>50), "
        f"mean percentile={summary['scoring']['mean_percentile']:.1f}"
        if percentiles
        else "no rows scored"
    )
    for reg, stats in summary["scoring_by_regulation"].items():
        if stats["n_scored"]:
            print(
                f"  {reg}: {stats['n_recovered']}/{stats['n_scored']} recovered, "
                f"mean percentile={stats['mean_percentile']:.1f}"
            )

    return summary


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="kiae271 differential persulfidome control-panel recovery evaluation"
    )
    p.add_argument(
        "--output",
        type=Path,
        default=Path("results/known_controls/kiae271_differential_recovery_v1.json"),
    )
    p.add_argument("--localization-min", type=float, default=0.75)
    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    run_kiae271_control_evaluation(
        output_path=args.output, localization_min=args.localization_min
    )
