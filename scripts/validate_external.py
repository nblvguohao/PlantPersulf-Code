#!/usr/bin/env python
"""Task 10 — strict external validation orchestrator.

Stitches together the pieces a cross-study predictive claim must survive and
writes one serialisable report under ``results/external_validation/<release>/``:

* leave-study-out fold metrics for the model release and a traditional
  baseline, keeping **every** seed (never only the best);
* a protein-cluster bootstrap CI and a permutation p-value on any per-site
  scored predictions supplied (``--scored``); skipped-with-a-note otherwise;
* known-mechanism control recovery: leakage check, independent-unit count, and
  a report that includes failed/unmappable controls — reusing the registered
  Zhang-lab controls and the tested integrity rules;
* a manifest recording the SHA256 of every input and the mandatory same-lab /
  same-chemistry / same-species limitation.

The report is deliberately fail-closed and honest: missing inputs are recorded
as gaps, never silently filled. It computes evidence for the Gate-2 conclusion
gate (``plantpersulf.evaluation.conclusion_gate``) but does not itself claim
predictive value.

Usage::

    python scripts/validate_external.py --model-release pu_ranker_v1
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

from plantpersulf.evaluation.bootstrap import cluster_bootstrap_ci
from plantpersulf.evaluation.external_validation import (
    ControlRecord,
    ExternalValidation,
    StudyFoldMetrics,
    build_recovery_report,
    count_independent_validation_units,
    find_control_training_leakage,
)
from plantpersulf.evaluation.metrics import average_precision
from plantpersulf.evaluation.permutation import permutation_test

LIMITATION = (
    "Both training studies are from the same laboratory (Romero/Gotor, "
    "Universidad de Sevilla), the same tag-switch persulfidation chemistry, "
    "and the same species (Arabidopsis thaliana). Leave-study-out is "
    "method-level validation only; it does NOT demonstrate chemically, "
    "laboratory-, or species-independent predictive value."
)


def _sha256(path: Path) -> str:
    from plantpersulf.provenance.hashing import hash_file

    return hash_file(path, "sha256")


def _parse_fold_study(model_field: str) -> str | None:
    """The study-split runner encodes the fold as ``leave_<study>_out|...``."""
    head = model_field.split("|", 1)[0]
    if head.startswith("leave_") and head.endswith("_out"):
        return head[len("leave_") : -len("_out")]
    return None


def _read_metrics(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return [dict(r) for r in csv.DictReader(handle, delimiter="\t")]


def _collect_fold_metrics(
    release_rows: list[dict[str, str]],
    baseline_rows: list[dict[str, str]],
    model_tag: str,
    baseline_tag: str,
) -> list[StudyFoldMetrics]:
    """Group test AP by held-out study for the release model and a baseline,
    keeping all seeds."""

    def _by_study(
        rows: list[dict[str, str]], tag: str
    ) -> dict[str, list[tuple[int, float]]]:
        out: dict[str, list[tuple[int, float]]] = {}
        for r in rows:
            model = r["model"]
            if tag not in model:
                continue
            study = _parse_fold_study(model)
            if study is None or not r.get("test_ap"):
                continue
            out.setdefault(study, []).append((int(r["seed"]), float(r["test_ap"])))
        return out

    release_by = _by_study(release_rows, model_tag)
    baseline_by = _by_study(baseline_rows, baseline_tag)

    folds: list[StudyFoldMetrics] = []
    for study in sorted(release_by):
        rel = sorted(release_by[study])
        base = dict(baseline_by.get(study, []))
        seeds = [s for s, _ in rel]
        folds.append(
            StudyFoldMetrics(
                holdout_study=study,
                model=model_tag,
                seeds=seeds,
                test_ap=[ap for _, ap in rel],
                baseline_test_ap=[base.get(s, float("nan")) for s in seeds],
            )
        )
    return folds


def _control_records(recovery_path: Path | None) -> list[ControlRecord]:
    """Build ControlRecords from the registered Zhang-lab controls, attaching
    percentile ranks from a prior recovery run if one is available."""
    from evaluate_known_controls import (
        REGISTERED_CONTROLS,  # type: ignore[import-not-found]
    )

    ranks: dict[str, float] = {}
    if recovery_path is not None and recovery_path.is_file():
        loaded = json.loads(recovery_path.read_text(encoding="utf-8"))
        for row in loaded:
            pct = row.get("percentile_rank")
            if pct is not None:
                ranks[str(row["mechanism_lineage_id"])] = float(pct)

    records: list[ControlRecord] = []
    for ctrl in REGISTERED_CONTROLS:
        lineage = str(ctrl["mechanism_lineage_id"])
        records.append(
            ControlRecord(
                mechanism_lineage_id=lineage,
                gene=str(ctrl["gene"]),
                uniprot_accession=str(ctrl["uniprot_accession"]),
                cys_position=int(ctrl["cys_position"]),
                status=str(ctrl["status"]),
                percentile_rank=ranks.get(lineage),
            )
        )
    return records


def _training_positives(benchmark_path: Path) -> set[tuple[str, int]]:
    positives: set[tuple[str, int]] = set()
    with benchmark_path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["label"] == "positive":
                positives.add(
                    (row["protein_accession"], int(row["cys_position_in_protein"]))
                )
    return positives


def _read_scored(path: Path) -> list[tuple[float, str, str]]:
    """Read per-site scored predictions: columns score,label,cluster_id."""
    rows: list[tuple[float, str, str]] = []
    with path.open(encoding="utf-8", newline="") as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            rows.append((float(r["score"]), r["label"], r["cluster_id"]))
    return rows


def run_external_validation(
    release: str,
    benchmark_path: Path,
    cluster_path: Path,
    baseline_tag: str,
    model_tag: str,
    recovery_path: Path | None,
    scored_path: Path | None,
    output_dir: Path,
) -> ExternalValidation:
    release_metrics = Path("results/experiments") / release / "metrics.tsv"
    baseline_metrics = (
        Path("results/experiments") / "baseline_leave_study_out_v1" / "metrics.tsv"
    )

    fold_metrics: list[StudyFoldMetrics] = []
    inputs: dict[str, str] = {}
    if release_metrics.is_file() and baseline_metrics.is_file():
        fold_metrics = _collect_fold_metrics(
            _read_metrics(release_metrics),
            _read_metrics(baseline_metrics),
            model_tag,
            baseline_tag,
        )
        inputs["release_metrics_sha256"] = _sha256(release_metrics)
        inputs["baseline_metrics_sha256"] = _sha256(baseline_metrics)
    else:
        print(
            f"NOTE: release metrics {release_metrics} not found — "
            "fold metrics skipped."
        )

    # --- bootstrap + permutation on per-site scores if provided ---
    bootstrap: dict[str, Any] = {"status": "skipped_no_per_site_scores"}
    permutation: dict[str, Any] = {"status": "skipped_no_per_site_scores"}
    if scored_path is not None and scored_path.is_file():
        scored = _read_scored(scored_path)
        ci = cluster_bootstrap_ci(scored, average_precision, n_boot=1000, seed=0)
        bootstrap = {
            "metric": "average_precision", "point": ci.point,
            "lower": ci.lower, "upper": ci.upper, "n_boot": ci.n_boot,
        }
        perm = permutation_test(
            [(s, y) for s, y, _ in scored], average_precision, n_perm=1000, seed=0
        )
        permutation = {
            "metric": "average_precision", "observed": perm.observed,
            "p_value": perm.p_value, "n_perm": perm.n_perm,
        }
        inputs["scored_sha256"] = _sha256(scored_path)

    # --- known-mechanism control recovery ---
    controls = _control_records(recovery_path)
    training_positives = _training_positives(benchmark_path)
    leaks = find_control_training_leakage(controls, training_positives)
    if leaks:
        raise RuntimeError(
            "CONTROL LEAKAGE: "
            + ", ".join(c.mechanism_lineage_id for c in leaks)
            + " is both a training positive and a reported control"
        )
    report = build_recovery_report(controls)
    units = count_independent_validation_units(controls)

    inputs["benchmark_sha256"] = _sha256(benchmark_path)
    if cluster_path.is_file():
        inputs["cluster_file_sha256"] = _sha256(cluster_path)
    else:
        print(f"NOTE: cluster file {cluster_path} not found — SHA256 not recorded.")

    validation = ExternalValidation(
        release=release,
        fold_metrics=fold_metrics,
        bootstrap=bootstrap,
        permutation=permutation,
        control_report=report,
        independent_units=units,
        control_leakage=[c.mechanism_lineage_id for c in leaks],
        limitation=LIMITATION,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "external_validation.json").write_text(
        json.dumps(validation.to_dict(), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    with (output_dir / "control_recovery.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(report[0]), delimiter="\t", lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(report)
    manifest = {
        "release": release,
        "inputs": inputs,
        "independent_validation_units": units,
        "limitation": LIMITATION,
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(
        f"External validation written to {output_dir} "
        f"(folds={len(fold_metrics)}, independent_units={units})"
    )
    return validation


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="strict external validation report")
    p.add_argument("--model-release", default="pu_ranker_v1")
    p.add_argument(
        "--benchmark", type=Path,
        default=Path("data/processed/benchmark_v1/sites.tsv"),
    )
    p.add_argument(
        "--clusters", type=Path,
        default=Path("data/processed/clusters/protein_clusters_v1.tsv"),
    )
    p.add_argument("--baseline-tag", default="pu_logistic")
    p.add_argument("--model-tag", default="structure_ranker:full")
    p.add_argument("--recovery", type=Path, default=None)
    p.add_argument("--scored", type=Path, default=None)
    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    out = Path("results/external_validation") / str(args.model_release)
    run_external_validation(
        release=str(args.model_release),
        benchmark_path=args.benchmark,
        cluster_path=args.clusters,
        baseline_tag=str(args.baseline_tag),
        model_tag=str(args.model_tag),
        recovery_path=args.recovery,
        scored_path=args.scored,
        output_dir=out,
    )
