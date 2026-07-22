#!/usr/bin/env python
"""Task 10 — strict external validation orchestrator.

Stitches together the pieces a cross-study predictive claim must survive and
writes one serialisable report under ``results/external_validation/<release>/``:

* leave-study-out fold metrics for the model release and a traditional
  baseline, keeping **every** seed (never only the best);
* a protein-cluster bootstrap CI and a permutation p-value on any per-site
  scored predictions supplied (``--scored``); skipped-with-a-note otherwise;
* a **paired effect CI** (release vs baseline scored on the same held-out
  rows, clusters resampled) when ``--scored-baseline`` is also supplied;
* a structure-gain measurement (release vs structure-ablated arm, paired per
  fold x seed run) when ``--scored-ablated`` is supplied, plus a
  top-cluster-dominance measurement for the single-cluster condition;
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

    python scripts/validate_external.py --model-release pu_ranker_v1 \
        --scored results/external_validation/pu_ranker_v1/scored/model.tsv \
        --scored-baseline results/external_validation/pu_ranker_v1/scored/baseline.tsv \
        --scored-ablated results/external_validation/pu_ranker_v1/scored/ablated.tsv \
        --recovery results/known_controls/recovery_v1.json
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

from plantpersulf.evaluation.bootstrap import cluster_bootstrap_ci
from plantpersulf.evaluation.effect_size import (
    paired_cluster_bootstrap_delta_ci,
    paired_delta_ci,
    top_cluster_dominance,
)
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
    """Build ControlRecords from the registered known-mechanism controls,
    attaching percentile ranks from a prior recovery run if one is available."""
    from plantpersulf.evaluation.known_controls import REGISTERED_CONTROLS

    ranks: dict[str, float] = {}
    if recovery_path is not None and recovery_path.is_file():
        loaded = json.loads(recovery_path.read_text(encoding="utf-8"))
        for row in loaded:
            pct = row.get("percentile_rank")
            if pct is not None:
                ranks[str(row["mechanism_lineage_id"])] = float(pct)

    records: list[ControlRecord] = []
    for ctrl in REGISTERED_CONTROLS:
        records.append(
            ControlRecord(
                mechanism_lineage_id=ctrl.mechanism_lineage_id,
                gene=ctrl.gene,
                uniprot_accession=ctrl.uniprot_accession,
                cys_position=ctrl.cys_position,
                status=ctrl.status,
                percentile_rank=ranks.get(ctrl.mechanism_lineage_id),
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


# ---------------------------------------------------------------------------
# rich per-site scored files (fold/seed/has_structure aware)
# ---------------------------------------------------------------------------

_SCORED_RICH_COLUMNS = (
    "fold",
    "seed",
    "protein_accession",
    "cys_position",
    "label",
    "cluster_id",
    "has_structure",
    "score",
)


def _read_scored_rich(path: Path) -> list[dict[str, Any]]:
    """Read a rich scored file written by ``scripts/score_release.py``."""
    rows: list[dict[str, Any]] = []
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != _SCORED_RICH_COLUMNS:
            raise RuntimeError(f"scored file has invalid columns: {path}")
        for r in reader:
            rows.append(
                {
                    "fold": r["fold"],
                    "seed": int(r["seed"]),
                    "key": (r["protein_accession"], int(r["cys_position"])),
                    "label": r["label"],
                    "cluster_id": r["cluster_id"],
                    "has_structure": r["has_structure"] == "1",
                    "score": float(r["score"]),
                }
            )
    return rows


def _ensemble_by_fold(
    rows: list[dict[str, Any]],
) -> dict[str, list[dict[str, Any]]]:
    """Mean-ensemble scores across seeds, per fold.

    The test subsample is deterministic, so every seed of one fold scores the
    identical row set; a mismatch means the input files are inconsistent and
    ensembling must fail loudly rather than pair the wrong rows."""
    by_fold_seed: dict[str, dict[int, list[dict[str, Any]]]] = {}
    for r in rows:
        by_fold_seed.setdefault(r["fold"], {}).setdefault(r["seed"], []).append(r)

    out: dict[str, list[dict[str, Any]]] = {}
    for fold, seed_groups in sorted(by_fold_seed.items()):
        seed_keys = {
            seed: [r["key"] for r in group] for seed, group in seed_groups.items()
        }
        ref_seed = sorted(seed_keys)[0]
        ref = seed_keys[ref_seed]
        for seed, keys in seed_keys.items():
            if keys != ref:
                raise RuntimeError(
                    f"scored rows differ across seeds in fold {fold} "
                    f"(seed {seed} vs seed {ref_seed}) — cannot ensemble"
                )
        ensembles: list[dict[str, Any]] = []
        for i, key in enumerate(ref):
            scores = [g[i]["score"] for g in (seed_groups[s] for s in seed_keys)]
            first = seed_groups[ref_seed][i]
            ensembles.append(
                {
                    "key": key,
                    "label": first["label"],
                    "cluster_id": first["cluster_id"],
                    "has_structure": first["has_structure"],
                    "score": sum(scores) / len(scores),
                }
            )
        out[fold] = ensembles
    return out


def _as_clustered(
    rows: list[dict[str, Any]],
) -> list[tuple[float, str, str]]:
    return [(r["score"], r["label"], r["cluster_id"]) for r in rows]


def _align_arms(
    model_rows: list[dict[str, Any]],
    other_rows: list[dict[str, Any]],
    arm_name: str,
) -> tuple[list[tuple[float, str, str]], list[tuple[float, str, str]]]:
    """Pair two arms row-by-row; both must cover the same sites."""
    keys_model = [r["key"] for r in model_rows]
    keys_other = [r["key"] for r in other_rows]
    if keys_model != keys_other:
        raise RuntimeError(
            f"model and {arm_name} scored rows differ — both arms must be "
            "scored on the same held-out rows"
        )
    return _as_clustered(model_rows), _as_clustered(other_rows)


def _per_seed_fold_aps(
    rows: list[dict[str, Any]],
) -> dict[tuple[str, int], float]:
    """AP per (fold, seed) run — the unit for run-level paired deltas."""
    from plantpersulf.evaluation.metrics import average_precision

    groups: dict[tuple[str, int], list[tuple[float, str]]] = {}
    for r in rows:
        groups.setdefault((r["fold"], r["seed"]), []).append((r["score"], r["label"]))
    out: dict[tuple[str, int], float] = {}
    for key, scored in groups.items():
        ap = average_precision(scored)
        if ap is not None:
            out[key] = ap
    return out


def _load_gate2_config(path: Path) -> tuple[bool, dict[str, Any]]:
    import yaml

    if not path.is_file():
        return False, {}
    loaded = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    return (
        bool(loaded.get("studies_are_independent", False)),
        dict(loaded.get("thresholds", {})),
    )


def _effect_evidence(
    model_rich: list[dict[str, Any]],
    baseline_rich: list[dict[str, Any]],
) -> dict[str, Any]:
    """Paired release-vs-baseline AP delta with a cluster bootstrap CI, per
    fold (seed-ensembled); the Gate-2 effect CI is the *minimum* per-fold
    lower bound — every held-out study must show the effect."""
    model_folds = _ensemble_by_fold(model_rich)
    baseline_folds = _ensemble_by_fold(baseline_rich)
    per_fold: list[dict[str, Any]] = []
    for fold in sorted(model_folds):
        if fold not in baseline_folds:
            continue
        model_cl, base_cl = _align_arms(
            model_folds[fold], baseline_folds[fold], "baseline"
        )
        ci = paired_cluster_bootstrap_delta_ci(model_cl, base_cl, n_boot=1000, seed=0)
        per_fold.append(
            {
                "fold": fold,
                "delta_point": ci.point,
                "delta_lower": ci.lower,
                "delta_upper": ci.upper,
                "n_boot": ci.n_boot,
                "n_rows": len(model_cl),
            }
        )
    lowers = [f["delta_lower"] for f in per_fold]
    return {
        "arms": "model vs baseline (paired, cluster bootstrap, seed-ensembled)",
        "per_fold": per_fold,
        "delta_ci_lower": min(lowers) if lowers else None,
    }


def _structure_gain_evidence(
    model_rich: list[dict[str, Any]],
    ablated_rich: list[dict[str, Any]],
) -> dict[str, Any]:
    """Structure-branch gain: paired per-(fold, seed) AP deltas of the release
    vs structure-ablated arm, bootstrapped over runs; plus a structured-subset
    diagnostic (the subset is small — 14/390 benchmark positives carry an
    AlphaFold structure — so it is reported, not gated on)."""
    model_aps = _per_seed_fold_aps(model_rich)
    ablated_aps = _per_seed_fold_aps(ablated_rich)
    deltas = [
        model_aps[k] - ablated_aps[k] for k in sorted(model_aps) if k in ablated_aps
    ]
    ci = paired_delta_ci(deltas, n_boot=2000, seed=0) if deltas else None

    subset: dict[str, Any]
    model_folds = _ensemble_by_fold(model_rich)
    ablated_folds = _ensemble_by_fold(ablated_rich)
    subset_rows: list[tuple[float, str, str]] = []
    subset_abl: list[tuple[float, str, str]] = []
    for fold in sorted(model_folds):
        if fold not in ablated_folds:
            continue
        m_struct = [r for r in model_folds[fold] if r["has_structure"]]
        a_struct = [r for r in ablated_folds[fold] if r["has_structure"]]
        if m_struct and a_struct:
            ms_cl, as_cl = _align_arms(m_struct, a_struct, "ablated")
            subset_rows.extend(ms_cl)
            subset_abl.extend(as_cl)
    n_pos = sum(1 for _, y, _ in subset_rows if y == "positive")
    if len(subset_rows) >= 10 and n_pos >= 2:
        sub_ci = paired_cluster_bootstrap_delta_ci(
            subset_rows, subset_abl, n_boot=1000, seed=0
        )
        subset = {
            "status": "measured",
            "n_rows": len(subset_rows),
            "n_positives": n_pos,
            "delta_point": sub_ci.point,
            "delta_lower": sub_ci.lower,
            "delta_upper": sub_ci.upper,
        }
    else:
        subset = {
            "status": "subset_too_small",
            "n_rows": len(subset_rows),
            "n_positives": n_pos,
        }
    return {
        "arms": "model vs structure-ablated (paired per fold x seed run)",
        "deltas": deltas,
        "gain_point": ci.point if ci else None,
        "gain_lower": ci.lower if ci else None,
        "gain_upper": ci.upper if ci else None,
        "structured_subset_diagnostic": subset,
    }


def _dominance_evidence(model_rich: list[dict[str, Any]]) -> dict[str, Any]:
    """Top-cluster dominance on the pooled, seed-ensembled model scores."""
    folds = _ensemble_by_fold(model_rich)
    pooled = [r for fold in sorted(folds) for r in folds[fold]]
    dom = top_cluster_dominance(_as_clustered(pooled))
    return {
        "top_cluster_id": dom.top_cluster_id,
        "top_cluster_rows": dom.top_cluster_rows,
        "ap_full": dom.ap_full,
        "ap_without_top": dom.ap_without_top,
        "retention_ratio": dom.retention_ratio,
        "single_cluster_driven": dom.driven,
    }


def _evaluate_gate2(
    validation: ExternalValidation,
    effect: dict[str, Any],
    structure: dict[str, Any],
    dominance: dict[str, Any],
    permutation: dict[str, Any],
) -> Any:
    """Build the Gate-2 evidence dict from the report and judge it mechanically.

    Evidence not (yet) demonstrated is passed as ``None``/conservative
    defaults, so unproven conditions fail rather than being assumed."""
    from plantpersulf.evaluation.conclusion_gate import evaluate_gate2

    independent, thresholds = _load_gate2_config(Path("configs/gate2_v1.yaml"))

    per_study: list[dict[str, Any]] = []
    for fold in validation.fold_metrics:
        full = [x for x in fold.test_ap if x == x]
        base = [x for x in fold.baseline_test_ap if x == x]
        if not full or not base:
            continue
        per_study.append(
            {
                "study": fold.holdout_study,
                "full_ap": sum(full) / len(full),
                "baseline_ap": sum(base) / len(base),
            }
        )

    delta_ci_lower = effect.get("delta_ci_lower")
    perm_p = permutation.get("p_value") if "p_value" in permutation else None
    single_cluster_driven = bool(dominance.get("single_cluster_driven", True))
    structure_gain = structure.get("gain_lower")

    evidence: dict[str, Any] = {
        "studies_are_independent": independent,
        "per_study": per_study,
        "delta_ci_lower": delta_ci_lower,
        "permutation_p": perm_p,
        "control_leakage": validation.control_leakage,
        "independent_units": validation.independent_units,
        "structure_gain": structure_gain,
        "single_cluster_driven": single_cluster_driven,
    }
    return evaluate_gate2(evidence, thresholds)


def run_external_validation(
    release: str,
    benchmark_path: Path,
    cluster_path: Path,
    baseline_tag: str,
    model_tag: str,
    recovery_path: Path | None,
    scored_path: Path | None,
    output_dir: Path,
    scored_baseline_path: Path | None = None,
    scored_ablated_path: Path | None = None,
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
            f"NOTE: release metrics {release_metrics} not found — fold metrics skipped."
        )

    # --- bootstrap + permutation on per-site scores if provided ---
    bootstrap: dict[str, Any] = {"status": "skipped_no_per_site_scores"}
    permutation: dict[str, Any] = {"status": "skipped_no_per_site_scores"}
    effect: dict[str, Any] = {"status": "skipped_no_paired_baseline_scores"}
    structure: dict[str, Any] = {"status": "skipped_no_ablated_scores"}
    dominance: dict[str, Any] = {"status": "skipped_no_per_site_scores"}
    if scored_path is not None and scored_path.is_file():
        model_rich = _read_scored_rich(scored_path)
        # Raw AP CI and permutation run on seed-ENSEMBLED pooled rows — the
        # five seeds score the identical sites, so pooling raw per-seed rows
        # would count every site five times.
        ens_folds = _ensemble_by_fold(model_rich)
        pooled = [r for fold in sorted(ens_folds) for r in ens_folds[fold]]
        scored = _as_clustered(pooled)
        ci = cluster_bootstrap_ci(scored, average_precision, n_boot=1000, seed=0)
        bootstrap = {
            "metric": "average_precision",
            "point": ci.point,
            "lower": ci.lower,
            "upper": ci.upper,
            "n_boot": ci.n_boot,
        }
        perm = permutation_test(
            [(s, y) for s, y, _ in scored], average_precision, n_perm=1000, seed=0
        )
        permutation = {
            "metric": "average_precision",
            "observed": perm.observed,
            "p_value": perm.p_value,
            "n_perm": perm.n_perm,
        }
        inputs["scored_sha256"] = _sha256(scored_path)

        dominance = _dominance_evidence(model_rich)
        if scored_baseline_path is not None and scored_baseline_path.is_file():
            effect = _effect_evidence(
                model_rich, _read_scored_rich(scored_baseline_path)
            )
            inputs["scored_baseline_sha256"] = _sha256(scored_baseline_path)
        else:
            print(
                "NOTE: no --scored-baseline — paired effect CI recorded as a "
                "gap (Gate-2 condition 2 stays unproven)."
            )
        if scored_ablated_path is not None and scored_ablated_path.is_file():
            structure = _structure_gain_evidence(
                model_rich, _read_scored_rich(scored_ablated_path)
            )
            inputs["scored_ablated_sha256"] = _sha256(scored_ablated_path)
        else:
            print(
                "NOTE: no --scored-ablated — structure gain recorded as a gap "
                "(Gate-2 condition 4 stays unproven)."
            )

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
        effect=effect,
        structure_gain=structure,
        cluster_dominance=dominance,
    )

    # --- Gate 2 conclusion gate (mechanical, honest) ---
    gate2 = _evaluate_gate2(validation, effect, structure, dominance, permutation)

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "external_validation.json").write_text(
        json.dumps(validation.to_dict(), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    (output_dir / "gate2_decision.json").write_text(
        json.dumps(gate2.to_dict(), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    with (output_dir / "control_recovery.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(report[0]),
            delimiter="\t",
            lineterminator="\n",
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
    print(f"Gate 2 decision: {gate2.decision}")
    for cond in gate2.conditions:
        mark = "PASS" if cond.passed else "FAIL"
        print(f"  [{mark}] {cond.name}: {cond.detail}")
    return validation


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="strict external validation report")
    p.add_argument("--model-release", default="pu_ranker_v1")
    p.add_argument(
        "--benchmark",
        type=Path,
        default=Path("data/processed/benchmark_v1/sites.tsv"),
    )
    p.add_argument(
        "--clusters",
        type=Path,
        default=Path("data/processed/clusters/protein_clusters_v1.tsv"),
    )
    p.add_argument("--baseline-tag", default="pu_logistic")
    p.add_argument("--model-tag", default="structure_ranker:full")
    p.add_argument("--recovery", type=Path, default=None)
    p.add_argument("--scored", type=Path, default=None)
    p.add_argument("--scored-baseline", type=Path, default=None)
    p.add_argument("--scored-ablated", type=Path, default=None)
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
        scored_baseline_path=args.scored_baseline,
        scored_ablated_path=args.scored_ablated,
    )
