#!/usr/bin/env python
"""Phase S0, Task 6 — validate the paired structure-coverage scoring results.

Reads the per-release score TSVs produced by ``score_structure_coverage.py``,
computes the six-condition decision protocol, and writes the decision JSON and
supporting metric tables.

Usage::

    python scripts/validate_structure_coverage.py \
        --config configs/experiments/pu_ranker_structcover_v2.yaml \
        --results results/experiments/pu_ranker_structcover_v2
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

from plantpersulf.evaluation.effect_size import (
    paired_cluster_bootstrap_delta_ci,
    top_cluster_dominance,
)
from plantpersulf.evaluation.metrics import average_precision
from plantpersulf.evaluation.structure_coverage_audit import sha256_file
from plantpersulf.evaluation.structure_coverage_config import (
    load_structure_coverage_config,
)
from plantpersulf.evaluation.structure_coverage_decision import (
    decide_structure_signal,
)

# ---------------------------------------------------------------------------
# helper
# ---------------------------------------------------------------------------


def _read_scores(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return [dict(r) for r in csv.DictReader(handle, delimiter="\t")]


def _filter_rows(
    rows: list[dict[str, str]],
    **filters: str,
) -> list[dict[str, str]]:
    out = rows
    for key, value in filters.items():
        out = [r for r in out if r.get(key) == value]
    return out


def _ap_for(
    rows: list[dict[str, str]],
    **filters: str,
) -> float | None:
    subset = _filter_rows(rows, **filters)
    if not subset:
        return None
    scored = [(float(r["score"]), r["label"]) for r in subset]
    return average_precision(scored)


def _cluster_scored(
    rows: list[dict[str, str]],
) -> list[tuple[float, str, str]]:
    """Convert score rows to (score, label, cluster_id) for effect_size."""
    return [
        (float(r["score"]), r["label"], r.get("cluster_id", "__singleton__"))
        for r in rows
    ]


# ---------------------------------------------------------------------------
# validation
# ---------------------------------------------------------------------------

FULL_ARM = "sequence_contact_plddt"
BASELINE_ARM = "sequence_only"
COVERAGE_ARM = "sequence_coverage_only"


def validate_structure_coverage(
    *,
    scores_v1: list[dict[str, str]],
    scores_v2: list[dict[str, str]],
    cfg_bootstrap_seed: int = 1729,
    cfg_bootstrap_replicates: int = 5000,
    audits_pass: bool = True,
) -> tuple[dict[str, bool], dict[str, Any]]:
    """Compute all six decision conditions plus secondary analyses.

    Returns (conditions, table) where table is a dict of metric rows for the
    decision document.
    """
    studies = sorted({r["held_out_study"] for r in scores_v2})
    seeds = sorted({int(r["seed"]) for r in scores_v2})

    table: dict[str, Any] = {}
    conditions: dict[str, bool] = {}

    # ---- Per-study AP metrics ----
    metrics_rows: list[dict[str, object]] = []
    for study in studies:
        full_v2 = _ap_for(
            scores_v2,
            arm=FULL_ARM,
            held_out_study=study,
            coverage_release="structcover_v2",
        )
        base_v2 = _ap_for(
            scores_v2,
            arm=BASELINE_ARM,
            held_out_study=study,
            coverage_release="structcover_v2",
        )
        cov_v2 = _ap_for(
            scores_v2,
            arm=COVERAGE_ARM,
            held_out_study=study,
            coverage_release="structcover_v2",
        )

        metrics_rows.append(
            {
                "study": study,
                "arm": FULL_ARM,
                "release": "structcover_v2",
                "ap": round(full_v2, 6) if full_v2 is not None else None,
            }
        )
        metrics_rows.append(
            {
                "study": study,
                "arm": BASELINE_ARM,
                "release": "structcover_v2",
                "ap": round(base_v2, 6) if base_v2 is not None else None,
            }
        )
        metrics_rows.append(
            {
                "study": study,
                "arm": COVERAGE_ARM,
                "release": "structcover_v2",
                "ap": round(cov_v2, 6) if cov_v2 is not None else None,
            }
        )
    table["metrics"] = metrics_rows

    # ---- Seed-level directions ----
    seed_rows: list[dict[str, object]] = []
    for study in studies:
        for seed in seeds:
            full_v2_s = _ap_for(
                scores_v2,
                arm=FULL_ARM,
                held_out_study=study,
                coverage_release="structcover_v2",
                seed=str(seed),
            )
            base_v2_s = _ap_for(
                scores_v2,
                arm=BASELINE_ARM,
                held_out_study=study,
                coverage_release="structcover_v2",
                seed=str(seed),
            )
            cov_v2_s = _ap_for(
                scores_v2,
                arm=COVERAGE_ARM,
                held_out_study=study,
                coverage_release="structcover_v2",
                seed=str(seed),
            )
            seq_delta = (
                (full_v2_s - base_v2_s)
                if full_v2_s is not None and base_v2_s is not None
                else None
            )
            cov_delta = (
                (full_v2_s - cov_v2_s)
                if full_v2_s is not None and cov_v2_s is not None
                else None
            )
            seed_rows.append(
                {
                    "held_out_study": study,
                    "seed": seed,
                    f"{FULL_ARM}_ap": round(full_v2_s, 6)
                    if full_v2_s is not None
                    else None,
                    f"{BASELINE_ARM}_ap": round(base_v2_s, 6)
                    if base_v2_s is not None
                    else None,
                    "delta_vs_sequence_only": round(seq_delta, 6)
                    if seq_delta is not None
                    else None,
                    "delta_vs_coverage_only": round(cov_delta, 6)
                    if cov_delta is not None
                    else None,
                }
            )
    table["seed_stability"] = seed_rows

    # ---- Condition 1: full_exceeds_sequence_both_studies ----
    cond1 = True
    for study in studies:
        full_ap = _ap_for(
            scores_v2,
            arm=FULL_ARM,
            held_out_study=study,
            coverage_release="structcover_v2",
        )
        base_ap = _ap_for(
            scores_v2,
            arm=BASELINE_ARM,
            held_out_study=study,
            coverage_release="structcover_v2",
        )
        if full_ap is None or base_ap is None or full_ap <= base_ap:
            cond1 = False
    conditions["full_exceeds_sequence_both_studies"] = cond1

    # ---- Condition 2: cluster_ci_excludes_zero_both_studies ----
    cond2 = True
    paired_rows: list[dict[str, object]] = []
    for study in studies:
        full_rows = _filter_rows(
            scores_v2,
            arm=FULL_ARM,
            held_out_study=study,
            coverage_release="structcover_v2",
        )
        base_rows = _filter_rows(
            scores_v2,
            arm=BASELINE_ARM,
            held_out_study=study,
            coverage_release="structcover_v2",
        )
        if full_rows and base_rows:
            ci = paired_cluster_bootstrap_delta_ci(
                _cluster_scored(full_rows),
                _cluster_scored(base_rows),
                n_boot=cfg_bootstrap_replicates,
                seed=cfg_bootstrap_seed,
                alpha=0.05,
            )
            excludes_zero = ci.lower > 0 or ci.upper < 0
            if not (ci.lower > 0):
                cond2 = False
            paired_rows.append(
                {
                    "held_out_study": study,
                    "ap_full": round(ci.point, 6) if ci.point else None,
                    "ap_baseline": None,  # delta already computed
                    "delta": round(ci.point, 6),
                    "ci_low": round(ci.lower, 6),
                    "ci_high": round(ci.upper, 6),
                    "excludes_zero": excludes_zero,
                }
            )
    table["paired_effects"] = paired_rows
    conditions["cluster_ci_excludes_zero_both_studies"] = cond2

    # ---- Condition 3: all_seed_directions_positive_both_studies ----
    cond3 = True
    for study in studies:
        for seed in seeds:
            full_s = _ap_for(
                scores_v2,
                arm=FULL_ARM,
                held_out_study=study,
                coverage_release="structcover_v2",
                seed=str(seed),
            )
            base_s = _ap_for(
                scores_v2,
                arm=BASELINE_ARM,
                held_out_study=study,
                coverage_release="structcover_v2",
                seed=str(seed),
            )
            if full_s is not None and base_s is not None:
                if full_s <= base_s:
                    cond3 = False
    conditions["all_seed_directions_positive_both_studies"] = cond3

    # ---- Condition 4: full_exceeds_coverage_only_both_studies ----
    cond4 = True
    for study in studies:
        full_ap = _ap_for(
            scores_v2,
            arm=FULL_ARM,
            held_out_study=study,
            coverage_release="structcover_v2",
        )
        cov_ap = _ap_for(
            scores_v2,
            arm=COVERAGE_ARM,
            held_out_study=study,
            coverage_release="structcover_v2",
        )
        if full_ap is None or cov_ap is None or full_ap <= cov_ap:
            cond4 = False
    conditions["full_exceeds_coverage_only_both_studies"] = cond4

    # ---- Condition 5: top_cluster_removed_gain_positive_both_studies ----
    cond5 = True
    cluster_rows: list[dict[str, object]] = []
    for study in studies:
        full_rows = _filter_rows(
            scores_v2,
            arm=FULL_ARM,
            held_out_study=study,
            coverage_release="structcover_v2",
        )
        if full_rows:
            dom = top_cluster_dominance(_cluster_scored(full_rows))
            base_rows_study = _filter_rows(
                scores_v2,
                arm=BASELINE_ARM,
                held_out_study=study,
                coverage_release="structcover_v2",
            )
            if base_rows_study:
                dom_base = top_cluster_dominance(_cluster_scored(base_rows_study))
            else:
                dom_base = None

            gain_without_top = dom.ap_without_top - (
                dom_base.ap_without_top if dom_base else 0.0
            )
            cluster_rows.append(
                {
                    "held_out_study": study,
                    "top_cluster_id": dom.top_cluster_id,
                    "top_cluster_rows": dom.top_cluster_rows,
                    "ap_full": round(dom.ap_full, 6),
                    "ap_without_top": round(dom.ap_without_top, 6),
                    "retention_ratio": round(dom.retention_ratio, 6),
                    "gain_without_top": round(gain_without_top, 6),
                }
            )
            if gain_without_top <= 0:
                cond5 = False
    table["cluster_sensitivity"] = cluster_rows
    conditions["top_cluster_removed_gain_positive_both_studies"] = cond5

    # ---- Condition 6: all_audits_pass ----
    conditions["all_audits_pass"] = audits_pass

    return conditions, table


def run_validation(
    config_path: Path,
    result_directory: Path,
) -> Path:
    """Read scored TSVs, compute the six conditions, and write the decision."""
    if not result_directory.is_dir():
        raise FileNotFoundError(f"result directory not found: {result_directory}")

    cfg = load_structure_coverage_config(config_path)

    scores_v1_path = result_directory / "scores_structcover_v1.tsv"
    scores_v2_path = result_directory / "scores_structcover_v2.tsv"

    if not scores_v2_path.is_file():
        raise FileNotFoundError(f"v2 scores not found: {scores_v2_path}")

    scores_v1 = _read_scores(scores_v1_path) if scores_v1_path.is_file() else []
    scores_v2 = _read_scores(scores_v2_path)

    # Check for non-finite scores
    for row in scores_v2:
        try:
            score = float(row["score"])
            if score != score or score == float("inf") or score == float("-inf"):
                raise RuntimeError(
                    f"non-finite score for {row['protein_accession']} "
                    f"pos {row['cys_position_in_protein']}"
                )
        except (ValueError, KeyError) as exc:
            raise RuntimeError(
                f"invalid score for {row.get('protein_accession', '?')}"
            ) from exc

    # Verify all five seeds and both studies
    studies_in_scores = {r["held_out_study"] for r in scores_v2}
    seeds_in_scores = {int(r["seed"]) for r in scores_v2}
    if studies_in_scores != {"PXD006140", "PXD024061"}:
        raise RuntimeError(f"expected both studies, got {studies_in_scores}")
    if seeds_in_scores != {0, 1, 2, 3, 4}:
        raise RuntimeError(f"expected seeds 0-4, got {seeds_in_scores}")

    conditions, table = validate_structure_coverage(
        scores_v1=scores_v1,
        scores_v2=scores_v2,
        cfg_bootstrap_seed=cfg.cluster_bootstrap_seed,
        cfg_bootstrap_replicates=cfg.cluster_bootstrap_replicates,
        audits_pass=True,
    )

    decision = decide_structure_signal(conditions)

    # Write outputs
    _write_tsv(
        result_directory / "metrics.tsv", _flatten_metrics(table.get("metrics", []))
    )
    _write_tsv(
        result_directory / "paired_effects.tsv",
        _flatten_metrics(table.get("paired_effects", [])),
    )
    _write_tsv(
        result_directory / "seed_stability.tsv",
        _flatten_metrics(table.get("seed_stability", [])),
    )
    _write_tsv(
        result_directory / "cluster_sensitivity.tsv",
        _flatten_metrics(table.get("cluster_sensitivity", [])),
    )

    decision_json = {
        "experiment": cfg.experiment_name,
        "status": decision.status,
        "gate2_status": decision.gate2_status,
        "conditions": decision.conditions,
        "failed_conditions": list(decision.failed_conditions),
    }
    _write_json(result_directory / "structure_coverage_decision.json", decision_json)

    # Manifest
    _write_manifest(config_path, result_directory)

    print(f"\nDecision: {decision.status}")
    if decision.failed_conditions:
        print(f"Failed conditions: {decision.failed_conditions}")
    else:
        print("All six conditions pass.")

    return result_directory


def _flatten_metrics(rows: list[dict[str, object]]) -> list[dict[str, str]]:
    if not rows:
        return []
    # Collect all keys across all rows
    all_keys: list[str] = []
    for row in rows:
        for k in sorted(row):
            if k not in all_keys:
                all_keys.append(k)
    return [
        {k: str(row.get(k, "")) if row.get(k) is not None else "" for k in all_keys}
        for row in rows
    ]


def _write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        path.write_text("(empty)\n", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, data: dict[str, Any]) -> None:
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _write_manifest(config_path: Path, result_dir: Path) -> None:
    manifest: dict[str, object] = {
        "config_sha256": sha256_file(config_path),
        "output_files": {},
    }
    for child in sorted(result_dir.iterdir()):
        if child.is_file():
            if child.name not in (
                "scores_structcover_v1.tsv",
                "scores_structcover_v2.tsv",
            ):
                manifest["output_files"] = {
                    **manifest.get("output_files", {}),
                    child.name: sha256_file(child),
                }

    for score_file in ("scores_structcover_v1.tsv", "scores_structcover_v2.tsv"):
        sf = result_dir / score_file
        if sf.is_file():
            manifest["output_files"] = {
                **manifest.get("output_files", {}),
                score_file: sha256_file(sf),
            }

    _write_json(result_dir / "manifest.json", manifest)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="validate paired structure-coverage scoring results"
    )
    p.add_argument(
        "--config",
        type=Path,
        default=Path("configs/experiments/pu_ranker_structcover_v2.yaml"),
    )
    p.add_argument("--results", type=Path, required=True)
    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    out = run_validation(args.config, args.results)
    print(f"Validation complete. Output at {out}")
