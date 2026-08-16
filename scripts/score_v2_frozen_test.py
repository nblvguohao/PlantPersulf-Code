#!/usr/bin/env python
"""Gate 1 — one-shot strict-track frozen test scoring (unlock-gated).

This is the ONLY code path that materialises and scores the frozen 20%
homology-cluster test partition of the multispecies v2 strict track. It is
gated by ``assert_test_unlocked``: the unlock record must exactly match the
config SHA256, the code revision, and the frozen-split SHA256, so the test
partition can only be scored against the exact frozen configuration.

One-shot discipline (NC entry-gate roadmap, Task 2 / Gate 1):

- the model arm (``structure_ranker``) is trained ONCE on the strict
  development panel (1:20 PU sampling, config seed) with the frozen v11
  hyperparameters verbatim and the freeze-day seed 20260813 (never
  performance-selected);
- the baseline arm (``pu_logistic``) is trained on the SAME development
  rows; its PU holdout fraction is 0.1, the modal selection of the five
  strict development folds (0.10/0.10/0.20/0.20/0.10) — a development-side
  choice, recorded here, not tuned on test;
- test rows are ALL cysteine sites on test-partition proteins (no PU
  subsampling — subsampling is a training-panel concept);
- the statistical report is composed by the frozen Task 9.6 orchestrator
  ``build_multispecies_statistical_report`` (per-species metrics, pooled
  primary AP, cluster bootstrap delta with the frozen policy
  n_boot=10,000 / seed=20260811, literature-track lead check, claim gate);
- outputs: per-site score TSVs for both arms, the statistical report JSON,
  a copy of the unlock record, and a run manifest that is audited with
  ``audit_v2_run_manifest`` before the script exits.

Usage::

    python scripts/score_v2_frozen_test.py \
        --config configs/experiments/multispecies_v2_global_clusters_v11.yaml \
        --test-unlock configs/gates/v2_frozen_test_unlock_20260813.json
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import yaml  # noqa: E402

from plantpersulf.benchmark.multispecies_splits import (  # noqa: E402
    TEST_SPLIT,
    load_global_cluster_table,
)
from plantpersulf.evaluation.comparable_track import (  # noqa: E402
    build_registered_structure_features,
)
from plantpersulf.evaluation.multispecies_reporting import (  # noqa: E402
    build_multispecies_statistical_report,
    compute_species_metric,
    strict_track_bootstrap_delta,
)
from plantpersulf.models.structure_ranker import (  # noqa: E402
    AblationConfig,
    BranchFeatures,
    fit_structure_ranker,
    score_structure_ranker_bundle,
)
from plantpersulf.models.traditional import pu_logistic_regression_scores  # noqa: E402
from plantpersulf.proteomics.multispecies_v2_dataset import (  # noqa: E402
    MultispeciesV2SiteRow,
)
from plantpersulf.proteomics.multispecies_v2_sources import (  # noqa: E402
    development_positive_rows,
    reference_cysteines,
    sequence_feature_map,
)
from plantpersulf.workflows.multispecies_v2 import (  # noqa: E402
    FROZEN_TEST_MODELS,
    FROZEN_TEST_SEED,
    FROZEN_TEST_TRACK,
    TrainingEvent,
    append_training_event,
    audit_v2_run_manifest,
    build_task_fingerprint,
    collect_runtime_environment,
    fit_v2_development_pipeline,
    run_multispecies_experiment,
    write_run_manifest,
    write_task_checkpoint,
)
from scripts.train_multispecies_v2 import (  # noqa: E402
    _code_revision,
    _development_rows,
    _dirty_code_paths,
    _dirty_paths,
    _strict_runtime_input_paths,
)

PRIMARY_SPECIES = ("arabidopsis", "rice", "tomato")
PRESSURE_SPECIES = ("magnaporthe",)
BASELINE_HOLDOUT = 0.1  # modal selection of the 5 strict development folds

DEFAULT_OUTPUT_DIR = (
    _REPO_ROOT
    / "results"
    / "experiments"
    / "multispecies_v2_global_clusters_v11"
    / "frozen_test"
)


def _sha256_bytes(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("configs/experiments/multispecies_v2_global_clusters_v11.yaml"),
    )
    parser.add_argument("--test-unlock", type=Path, required=True)
    parser.add_argument(
        "--evidence-config",
        type=Path,
        default=Path("configs/experiments/multispecies_v2_global_clusters_v4.yaml"),
        help=(
            "config carrying the raw positive_evidence paths; v11 shares the "
            "v4 split/cluster/proteome lineage, so v4 is the canonical source"
        ),
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--device", type=str, default=None)
    parser.add_argument("--batch-size", type=int, default=16384)
    return parser


def _build_test_rows(
    evidence: dict[str, object],
    frozen,
    clusters,
    proteomes: dict[str, dict[str, str]],
) -> list[MultispeciesV2SiteRow]:
    """Materialise frozen-test site rows through the unlocked boundary.

    ``build_split_bound_v2_rows`` deliberately omits test rows, so this
    boundary builds them directly: same parsers, same frozen split, same
    cluster table — the test analogue of the development row builder.
    """
    split_by_protein = {row.global_protein_id: row for row in frozen.rows}
    cluster_by_protein = {row.global_protein_id: row for row in clusters}
    test_allowed: dict[str, set[str]] = {}
    for row in frozen.rows:
        if row.split == TEST_SPLIT:
            test_allowed.setdefault(row.species, set()).add(
                row.global_protein_id.split("|", 1)[1]
            )
    if not test_allowed:
        raise RuntimeError("frozen split contains no test proteins")
    n_test_proteins = sum(len(value) for value in test_allowed.values())
    print(f"test proteins: {n_test_proteins}")

    evidence_paths = {
        key: Path(value) for key, value in evidence.items() if key != "registry"
    }
    positives = development_positive_rows(
        benchmark=evidence_paths["arabidopsis_benchmark"],
        tomato_xlsx=evidence_paths["tomato_kiae271"],
        rice_sd01=evidence_paths["rice_sd01"],
        rice_sd04=evidence_paths["rice_sd04"],
        rice_ss_all=evidence_paths["rice_ss_all"],
        magnaporthe_tsv=evidence_paths["magnaporthe_sites"],
        proteomes=proteomes,
        allowed_proteins=test_allowed,
    )
    positive_studies: dict[tuple[str, str, int], set[str]] = {}
    for row in positives:
        protein_id = f"{row.species}|{row.protein_accession}"
        split_row = split_by_protein.get(protein_id)
        if split_row is None or split_row.split != TEST_SPLIT:
            raise RuntimeError(
                f"test positive outside the test partition: {protein_id}"
            )
        positive_studies.setdefault(
            (row.species, row.protein_accession, row.cys_position), set()
        ).add(row.study_accession)

    rows: list[MultispeciesV2SiteRow] = []
    for species, accession, position in sorted(
        reference_cysteines(proteomes, test_allowed)
    ):
        protein_id = f"{species}|{accession}"
        split_row = split_by_protein[protein_id]
        cluster_row = cluster_by_protein.get(protein_id)
        if cluster_row is None or cluster_row.cluster_id != split_row.cluster_id:
            raise RuntimeError(
                f"cluster/split disagreement on test protein: {protein_id}"
            )
        studies = tuple(
            sorted(positive_studies.get((species, accession, position), set()))
        )
        rows.append(
            MultispeciesV2SiteRow(
                species=species,
                protein_accession=accession,
                cys_position=position,
                label="positive" if studies else "unlabeled",
                study_accessions=studies,
                global_protein_id=protein_id,
                cluster_id=split_row.cluster_id,
                split=TEST_SPLIT,
                development_fold=None,
            )
        )
    n_pos = sum(1 for row in rows if row.label == "positive")
    print(f"test rows: {len(rows)} sites, {n_pos} positives")
    if n_pos == 0:
        raise RuntimeError("frozen test partition contains no positives")
    return rows


def _branch_features(
    rows: list[MultispeciesV2SiteRow],
    sequence_features: dict[tuple[str, int], tuple[float, ...]],
    structure_features: dict[tuple[str, int], tuple[tuple[float, ...], bool]],
    keys: list[tuple[str, int]],
) -> BranchFeatures:
    return BranchFeatures(
        sequence=[list(sequence_features[key]) for key in keys],
        esm=[[0.0] for _ in keys],
        structure=[list(structure_features[key][0]) for key in keys],
        structure_mask=[structure_features[key][1] for key in keys],
        study_ids=None,
    )


def main(argv: list[str] | None = None) -> None:
    args = build_argument_parser().parse_args(argv)
    code_revision = _code_revision()

    # --- the unlock gate: everything below runs only with an exact record ----
    prepared = run_multispecies_experiment(
        args.config,
        score_frozen_test=True,
        test_unlock_path=args.test_unlock,
        code_revision=code_revision,
    )
    frozen = prepared.split
    cfg = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    if not isinstance(cfg, dict):
        raise RuntimeError("multispecies v2 config is required")
    global_cfg = cfg.get("global_mmseqs2")
    clusters = load_global_cluster_table(Path(str(global_cfg["cluster_table"])))
    print(f"unlock gate passed; frozen split sha256={frozen.sha256}")

    # --- rows ----------------------------------------------------------------
    dev_rows, _frozen_check, proteomes = _development_rows(args.config)
    if _frozen_check.sha256 != frozen.sha256:
        raise RuntimeError("development row builder disagrees with the frozen split")
    dev_rows = sorted(
        dev_rows, key=lambda row: (row.global_protein_id, row.cys_position)
    )
    n_dev_pos = sum(1 for row in dev_rows if row.label == "positive")
    print(f"development panel: {len(dev_rows)} rows, {n_dev_pos} positives (1:20 PU)")

    evidence_cfg = yaml.safe_load(args.evidence_config.read_text(encoding="utf-8"))
    if not isinstance(evidence_cfg, dict) or "positive_evidence" not in evidence_cfg:
        raise RuntimeError("evidence config must carry positive_evidence paths")
    test_rows = _build_test_rows(
        evidence_cfg["positive_evidence"], frozen, clusters, proteomes
    )
    for species in (*PRIMARY_SPECIES, *PRESSURE_SPECIES):
        n_pos = sum(
            1
            for row in test_rows
            if row.species == species and row.label == "positive"
        )
        if n_pos == 0:
            raise RuntimeError(f"frozen test has no positives for {species}")
        print(f"test positives: {species}={n_pos}")
    dev_proteins = {row.global_protein_id for row in dev_rows}
    if any(row.global_protein_id in dev_proteins for row in test_rows):
        raise RuntimeError("development/test protein overlap detected")

    # --- features -------------------------------------------------------------
    all_rows = list(dev_rows) + list(test_rows)
    raw_features = sequence_feature_map(all_rows, proteomes)
    if any(any(value is None for value in values) for values in raw_features.values()):
        raise RuntimeError("sequence features contain missing values")
    sequence_features = {
        key: tuple(float(value) for value in values if value is not None)
        for key, values in raw_features.items()
    }
    inputs_cfg = cfg.get("comparison_inputs", {})
    build_cfg = cfg.get("comparison_feature_build")
    structure_features = build_registered_structure_features(
        set(sequence_features),
        registry_path=Path(str(inputs_cfg["structure_features"])),
        registry_base=Path(str(build_cfg["structure_registry_base"])),
    )

    dev_keys = [(row.global_protein_id, row.cys_position) for row in dev_rows]
    test_keys = [(row.global_protein_id, row.cys_position) for row in test_rows]

    # --- model arm: structure_ranker, frozen hyperparameters, one fit --------
    params = cfg["literature_baseline_parameters"]["structure_ranker"]
    device = args.device or "cpu"
    print(
        f"model arm: fit structure_ranker on {len(dev_keys)} dev rows "
        f"(seed={FROZEN_TEST_SEED}, device={device})"
    )
    model_start = time.monotonic()
    bundle = fit_structure_ranker(
        _branch_features(dev_rows, sequence_features, structure_features, dev_keys),
        [row.label for row in dev_rows],
        seed=FROZEN_TEST_SEED,
        ablation=AblationConfig(use_esm=False, use_study_context=False),
        hidden=int(params["hidden"]),
        dropout=float(params["dropout"]),
        epochs=int(params["epochs"]),
        lr=float(params["lr"]),
        holdout_fraction=float(params["holdout_fraction"]),
        n_mc_dropout=int(params["n_mc_dropout"]),
        device_name=args.device,
        batch_size=None,
    )
    model_output = score_structure_ranker_bundle(
        bundle,
        _branch_features(test_rows, sequence_features, structure_features, test_keys),
        device_name=args.device,
        batch_size=args.batch_size,
    )
    model_scores = dict(zip(test_keys, model_output.scores, strict=True))
    model_wall = time.monotonic() - model_start

    # --- baseline arm: pu_logistic on the same dev rows ----------------------
    print(
        f"baseline arm: pu_logistic (holdout={BASELINE_HOLDOUT}, "
        "modal dev-fold choice)"
    )
    baseline_start = time.monotonic()
    pipeline = fit_v2_development_pipeline(dev_rows, sequence_features)
    baseline_scores_list = pu_logistic_regression_scores(
        pipeline.transform(dev_rows, sequence_features),
        [row.label for row in dev_rows],
        pipeline.transform(test_rows, sequence_features),
        seed=FROZEN_TEST_SEED,
        holdout_fraction=BASELINE_HOLDOUT,
    )
    baseline_scores = dict(zip(test_keys, baseline_scores_list, strict=True))
    baseline_wall = time.monotonic() - baseline_start

    # --- statistics (frozen Task 9.6 policy) ----------------------------------
    def scored_by_species(scores: dict[tuple[str, int], float]) -> dict[str, list]:
        out: dict[str, list] = {}
        for row in test_rows:
            key = (row.global_protein_id, row.cys_position)
            out.setdefault(row.species, []).append((scores[key], row.label))
        return out

    model_by_species = scored_by_species(model_scores)
    baseline_by_species = scored_by_species(baseline_scores)
    species_deltas = {
        species: (
            compute_species_metric(model_by_species[species], species).average_precision
            - compute_species_metric(
                baseline_by_species[species], species
            ).average_precision
        )
        for species in PRIMARY_SPECIES
    }
    primary_rows = [row for row in test_rows if row.species in PRIMARY_SPECIES]
    model_clustered = [
        (
            model_scores[(row.global_protein_id, row.cys_position)],
            row.label,
            row.cluster_id,
        )
        for row in primary_rows
    ]
    baseline_clustered = [
        (
            baseline_scores[(row.global_protein_id, row.cys_position)],
            row.label,
            row.cluster_id,
        )
        for row in primary_rows
    ]
    bootstrap = strict_track_bootstrap_delta(model_clustered, baseline_clustered)
    print(
        f"strict bootstrap delta (pooled primary species): "
        f"point={bootstrap.point:.4f} "
        f"CI=[{bootstrap.lower:.4f}, {bootstrap.upper:.4f}] "
        f"(n_boot={bootstrap.n_boot})"
    )

    # literature-track per-seed macro APs for the lead check
    summary_path = (
        _REPO_ROOT
        / "results"
        / "experiments"
        / "multispecies_v2_global_clusters_v11"
        / "literature_random_protein"
        / "summary.json"
    )
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    by_seed: dict[int, dict[str, float]] = {}
    for run in summary["runs"]:
        by_seed.setdefault(int(run["seed"]), {})[run["model"]] = float(
            run["three_crop_macro_average_precision"]
        )
    seeds = sorted(by_seed)
    candidate_series = [by_seed[seed]["structure_ranker"] for seed in seeds]
    baseline_series = [
        max(
            value
            for model, value in by_seed[seed].items()
            if model != "structure_ranker"
        )
        for seed in seeds
    ]

    report = build_multispecies_statistical_report(
        scored_by_species=model_by_species,
        primary_species=PRIMARY_SPECIES,
        pressure_species=PRESSURE_SPECIES,
        strict_ci_lower=bootstrap.lower,
        species_deltas=species_deltas,
        same_frozen_inputs=True,
        candidate_literature_macro_aps=candidate_series,
        baseline_literature_macro_aps=baseline_series,
    )
    report["baseline_arm"] = {
        "model": "pu_logistic",
        "holdout_fraction": BASELINE_HOLDOUT,
        "per_species": {
            species: {
                "average_precision": compute_species_metric(
                    baseline_by_species[species], species
                ).average_precision
            }
            for species in (*PRIMARY_SPECIES, *PRESSURE_SPECIES)
        },
    }
    report["strict_bootstrap"] = {
        "point": bootstrap.point,
        "lower": bootstrap.lower,
        "upper": bootstrap.upper,
        "n_boot": bootstrap.n_boot,
        "pooled_species": list(PRIMARY_SPECIES),
    }
    report["frozen_test"] = {
        "n_rows": len(test_rows),
        "n_positives": sum(1 for row in test_rows if row.label == "positive"),
        "seed": FROZEN_TEST_SEED,
        "device": device,
        "model_arm": (
            "structure_ranker (frozen v11 hyperparameters, "
            "fit on strict development panel only)"
        ),
    }
    print(f"claim_class: {report['claim_class']}")

    # --- outputs ---------------------------------------------------------------
    args.output_dir.mkdir(parents=True, exist_ok=True)
    score_files = (
        ("scores_model.tsv", model_scores),
        ("scores_baseline.tsv", baseline_scores),
    )
    for name, scores in score_files:
        path = args.output_dir / name
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(
                ["site_key", "cys_position", "species", "label", "cluster_id", "score"]
            )
            for row in test_rows:
                key = (row.global_protein_id, row.cys_position)
                writer.writerow(
                    [
                        row.global_protein_id,
                        row.cys_position,
                        row.species,
                        row.label,
                        row.cluster_id,
                        f"{scores[key]:.10g}",
                    ]
                )

    report_path = args.output_dir / "statistical_report.json"
    report_path.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    unlock_copy = args.output_dir / "unlock_record.json"
    unlock_copy.write_text(
        args.test_unlock.read_text(encoding="utf-8"), encoding="utf-8"
    )

    input_paths = _strict_runtime_input_paths(cfg)
    input_paths["test_unlock"] = args.test_unlock
    input_paths["evidence_config"] = args.evidence_config
    input_paths["literature_summary"] = summary_path
    dirty_paths = _dirty_paths()

    # one fingerprint + checkpoint + training event per arm (the audit requires
    # exact roster/log/checkpoint coverage)
    arm_walls = {"structure_ranker": model_wall, "pu_logistic": baseline_wall}
    fingerprints = {}
    checkpoints = {}
    checkpoint_dir = args.output_dir / "checkpoints"
    log_path = args.output_dir / "training.jsonl"
    if log_path.is_file():
        log_path.write_text("", encoding="utf-8")
    for model in FROZEN_TEST_MODELS:
        fingerprint = build_task_fingerprint(
            track=FROZEN_TEST_TRACK,
            fold=0,
            seed=FROZEN_TEST_SEED,
            model=model,
            input_paths=input_paths,
            config_path=args.config,
            code_revision=code_revision,
            dirty_paths=_dirty_code_paths(dirty_paths),
        )
        fingerprints[model] = fingerprint
        checkpoint_path = checkpoint_dir / f"{model}.checkpoint.json"
        write_task_checkpoint(
            checkpoint_path,
            fingerprint,
            {
                "arm": model,
                "n_dev_rows": len(dev_rows),
                "n_test_rows": len(test_rows),
            },
        )
        checkpoints[model] = checkpoint_path
        append_training_event(
            log_path,
            TrainingEvent(
                timestamp=datetime.now(timezone.utc).isoformat(),
                track=FROZEN_TEST_TRACK,
                fold=0,
                seed=FROZEN_TEST_SEED,
                model=model,
                epoch=0,
                loss=None,
                val_ap=None,
                lr=None,
                device=device,
                gpu_memory=None,
                wall_seconds=arm_walls[model],
            ),
        )
    first_fingerprint = fingerprints[FROZEN_TEST_MODELS[0]]
    manifest_path = args.output_dir / "manifest.json"
    write_run_manifest(
        manifest_path,
        config_sha256=prepared.config_sha256,
        config_path=args.config,
        split_sha256=frozen.sha256,
        code_revision=code_revision,
        code_sha256=first_fingerprint.code_sha256,
        input_sha256=first_fingerprint.input_sha256,
        command=[sys.executable, *sys.argv],
        device=device,
        environment=collect_runtime_environment(),
        artifacts={
            "training_log": log_path,
            "statistical_report": report_path,
            "scores_model": args.output_dir / "scores_model.tsv",
            "scores_baseline": args.output_dir / "scores_baseline.tsv",
            "unlock_record": unlock_copy,
        },
        checkpoint_paths={
            "model_arm": checkpoints["structure_ranker"],
            "baseline_arm": checkpoints["pu_logistic"],
        },
        task_roster=tuple(fingerprints[model] for model in FROZEN_TEST_MODELS),
        dirty=bool(dirty_paths),
        dirty_paths=dirty_paths,
    )
    audit_v2_run_manifest(manifest_path)
    print(f"manifest audited: {manifest_path}")
    print(f"report: {report_path}")


if __name__ == "__main__":
    main()
