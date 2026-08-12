#!/usr/bin/env python
"""Preflight and execute the development-only multispecies v2 path.

The default command validates the registered global MMseqs2 cluster table and
frozen split. ``--prepare-development`` builds the split-bound 1:20 PU panel,
prepares five homology-grouped folds, and fits every preprocessing and model
selection step on each training fold only. Passing ``--score-test`` requires a
matching unlock record and does not by itself materialize test labels.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import tempfile
import time
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import cast

import yaml

from plantpersulf.benchmark.literature_random_track import (
    bind_model_to_comparison_panel,
    build_literature_random_track,
    run_direct_comparison_roster,
)
from plantpersulf.benchmark.multispecies_splits import (
    DEVELOPMENT_SPLIT,
    FrozenMultispeciesSplit,
    load_frozen_multispecies_split,
    load_global_cluster_table,
)
from plantpersulf.evaluation.comparable_track import (
    build_registered_structure_features,
    comparator_statuses,
    run_sul_bertgru_adapter,
    summarize_comparable_scores,
    validate_sul_environment_manifest,
)
from plantpersulf.features.esm2_windows import (
    build_cysteine_windows,
    extract_esm2_window_artifact,
    load_window_embedding_artifact,
)
from plantpersulf.proteomics.multispecies_v2_dataset import (
    MultispeciesV2SiteRow,
    build_split_bound_v2_rows,
    prepare_v2_development_fold,
    sample_v2_unlabeled_panels,
)
from plantpersulf.proteomics.multispecies_v2_sources import (
    load_frozen_development_positives,
    reference_cysteines,
    sequence_feature_map,
)
from plantpersulf.provenance.audit import assert_registered_input
from plantpersulf.workflows.multispecies_v2 import (
    TaskFingerprint,
    TrainingEvent,
    append_training_event,
    audit_v2_run_manifest,
    build_task_fingerprint,
    collect_runtime_environment,
    load_resumable_checkpoint,
    run_multispecies_experiment,
    run_task_group,
    run_with_oom_batch_retry,
    select_v2_development_hyperparameters,
    train_v2_development_fold,
    write_run_manifest,
    write_task_checkpoint,
)


def _code_revision() -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"], check=True, capture_output=True, text=True
    )
    return completed.stdout.strip()


def _dirty_paths() -> tuple[str, ...]:
    completed = subprocess.run(
        ["git", "status", "--porcelain"], check=True, capture_output=True, text=True
    )
    return tuple(
        sorted(line[3:] for line in completed.stdout.splitlines() if len(line) >= 4)
    )


def _dirty_code_paths(dirty_paths: tuple[str, ...]) -> tuple[Path, ...]:
    """Return changed executable files whose bytes affect the active code state."""
    code_paths = tuple(
        Path(path)
        for path in dirty_paths
        if Path(path).suffix == ".py" and Path(path).is_file()
    )
    missing = [
        path
        for path in dirty_paths
        if Path(path).suffix == ".py" and not Path(path).is_file()
    ]
    if missing:
        raise RuntimeError(f"cannot fingerprint deleted dirty code paths: {missing}")
    return code_paths


def _strict_runtime_input_paths(cfg: dict[str, object]) -> dict[str, Path]:
    """List every physical input read by strict development task preparation."""
    references = cfg.get("reference_proteomes")
    strict = cfg.get("strict_cluster_holdout")
    global_clusters = cfg.get("global_mmseqs2")
    positives = cfg.get("development_positive_manifest")
    if not all(
        isinstance(value, dict) for value in (strict, global_clusters, positives)
    ) or not isinstance(references, list):
        raise RuntimeError("strict runtime inputs are not configured")
    paths: dict[str, Path] = {
        "development_positive_manifest": Path(str(positives["path"])),
        "frozen_split": Path(str(strict["split_path"])),
        "global_cluster_table": Path(str(global_clusters["cluster_table"])),
    }
    for reference in references:
        if not isinstance(reference, dict):
            raise RuntimeError("reference proteome entry is invalid")
        species = reference.get("species")
        source_path = reference.get("path")
        if not isinstance(species, str) or not isinstance(source_path, str):
            raise RuntimeError("reference proteome entry is invalid")
        paths[f"reference_proteome:{species}"] = Path(source_path)
    return paths


def _development_rows_sha256(rows: tuple[MultispeciesV2SiteRow, ...]) -> str:
    lines = [
        "\t".join(
            (
                row.global_protein_id,
                str(row.cys_position),
                row.label,
                row.cluster_id,
                str(row.development_fold),
            )
        )
        for row in sorted(
            rows, key=lambda row: (row.global_protein_id, row.cys_position)
        )
    ]
    return hashlib.sha256("\n".join(lines).encode()).hexdigest()


def _registered(path: Path) -> None:
    registries = (
        Path("data/registry/cross_crop_target_label_free_inputs_v1.tsv"),
        Path("data/registry/supplementary_sources.tsv"),
        Path("data/registry/model_inputs.tsv"),
    )
    for registry in registries:
        try:
            assert_registered_input(path, registry)
            return
        except RuntimeError:
            pass
    raise RuntimeError(f"input is absent from every approved registry: {path}")


def _development_rows(
    config: Path, *, sample_unlabeled: bool = True
) -> tuple[
    tuple[MultispeciesV2SiteRow, ...],
    FrozenMultispeciesSplit,
    dict[str, dict[str, str]],
]:
    """Build the only default v2 training dataset: development rows only."""
    cfg = yaml.safe_load(config.read_text(encoding="utf-8"))
    assert isinstance(cfg, dict)
    ref = {item["species"]: Path(item["path"]) for item in cfg["reference_proteomes"]}
    manifest_cfg = cfg.get("development_positive_manifest")
    if not isinstance(manifest_cfg, dict):
        raise RuntimeError("development_positive_manifest configuration is required")
    manifest_path = Path(str(manifest_cfg["path"]))
    for path in [*ref.values(), manifest_path]:
        _registered(path)
    manifest_sha256 = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    if manifest_sha256 != manifest_cfg.get("sha256"):
        raise RuntimeError("development-positive manifest SHA256 mismatch")
    from plantpersulf.features.sequence import _load_proteome
    from plantpersulf.proteomics.pxd063170_sites import load_ensembl_fungi_proteome

    proteomes = {
        "arabidopsis": _load_proteome(ref["arabidopsis"]),
        "tomato": _load_proteome(ref["tomato"]),
        "rice": _load_proteome(ref["rice"]),
        "magnaporthe": load_ensembl_fungi_proteome(ref["magnaporthe"]),
    }
    _registered(Path(cfg["global_mmseqs2"]["cluster_table"]))
    clusters = load_global_cluster_table(Path(cfg["global_mmseqs2"]["cluster_table"]))
    _registered(Path(cfg["strict_cluster_holdout"]["split_path"]))
    frozen = load_frozen_multispecies_split(
        Path(cfg["strict_cluster_holdout"]["split_path"])
    )
    development_ids = {
        row.global_protein_id for row in frozen.rows if row.split == DEVELOPMENT_SPLIT
    }
    allowed = {
        species: {
            row.global_protein_id.split("|", 1)[1]
            for row in frozen.rows
            if row.species == species and row.global_protein_id in development_ids
        }
        for species in proteomes
    }
    site_rows = load_frozen_development_positives(manifest_path)
    all_cys = list(reference_cysteines(proteomes, allowed))
    rows = build_split_bound_v2_rows(
        positives=site_rows,
        all_cysteines=all_cys,
        clusters=clusters,
        frozen_split=frozen,
    )
    panel = (
        sample_v2_unlabeled_panels(
            rows,
            per_positive=cfg["unlabeled_panel"]["per_positive"],
            seed=cfg["unlabeled_panel"]["seed"],
        )
        if sample_unlabeled
        else rows
    )
    return panel, frozen, proteomes


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path, default=Path("configs/experiments/multispecies_v2.yaml")
    )
    parser.add_argument("--score-test", action="store_true")
    parser.add_argument("--test-unlock", type=Path)
    parser.add_argument("--prepare-development", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--prepare-literature-track", action="store_true")
    parser.add_argument("--build-comparison-features", action="store_true")
    parser.add_argument("--run-literature-baselines", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_argument_parser().parse_args(argv)
    prepared = run_multispecies_experiment(
        args.config,
        score_frozen_test=args.score_test,
        test_unlock_path=args.test_unlock,
        code_revision=_code_revision(),
    )
    print(
        f"prepared {prepared.cluster_count} clusters; "
        f"test_scoring_enabled={prepared.test_scoring_enabled}"
    )
    if args.prepare_development:
        rows, frozen, proteomes = _development_rows(args.config)
        features = sequence_feature_map(rows, proteomes)
        cfg = yaml.safe_load(args.config.read_text(encoding="utf-8"))
        if not isinstance(cfg, dict):
            raise RuntimeError("multispecies v2 config is required")
        output_cfg = cfg.get("output")
        compute_cfg = cfg.get("compute")
        if not isinstance(output_cfg, dict) or not isinstance(compute_cfg, dict):
            raise RuntimeError("output and compute configuration are required")
        output_directory = Path(str(output_cfg["directory"]))
        log_path = output_directory / str(output_cfg["jsonl_log"])
        checkpoint_directory = output_directory / "checkpoints"
        default_device = str(compute_cfg["device"])
        raw_gpu_map = compute_cfg.get("gpu_map", [])
        if not isinstance(raw_gpu_map, list) or not all(
            isinstance(value, str) for value in raw_gpu_map
        ):
            raise RuntimeError("compute.gpu_map must be a list of device names")
        gpu_map = tuple(raw_gpu_map)
        max_parallel_tasks = compute_cfg.get("max_parallel_tasks")
        if not isinstance(max_parallel_tasks, int) or max_parallel_tasks < 1:
            raise RuntimeError("compute.max_parallel_tasks must be a positive integer")
        code_revision = _code_revision()
        dirty_paths = _dirty_paths()
        input_paths = _strict_runtime_input_paths(cfg)
        dirty_code_paths = _dirty_code_paths(dirty_paths)

        def run_fold(fold: int, device: str) -> dict[str, object]:
            fingerprint = build_task_fingerprint(
                track="strict_cluster_holdout",
                fold=fold,
                seed=fold,
                model="pu_logistic",
                input_paths=input_paths,
                config_path=args.config,
                code_revision=code_revision,
                dirty_paths=dirty_code_paths,
            )
            checkpoint_path = checkpoint_directory / f"fold{fold}.checkpoint.json"
            prepared_fold = prepare_v2_development_fold(
                rows,
                validation_fold=fold,
                n_development_folds=frozen.n_development_folds,
            )
            if not prepared_fold.validation_rows:
                raise RuntimeError(
                    f"development fold has no v2 rows: {fold}; "
                    "frozen split requires reviewer remediation"
                )
            has_validation_positive = any(
                row.label == "positive" for row in prepared_fold.validation_rows
            )
            if not has_validation_positive:
                raise RuntimeError(
                    f"development fold has no positive v2 rows: {fold}; "
                    "frozen split requires reviewer remediation"
                )
            start = time.monotonic()
            if args.resume:
                state = load_resumable_checkpoint(checkpoint_path, fingerprint)
                selected_holdout = float(state["pu_holdout"])
                validation_ap = float(state["validation_ap"])
            else:
                selected_holdout = select_v2_development_hyperparameters(
                    prepared_fold.fit_rows,
                    prepared_fold.validation_rows,
                    features,
                    candidates=(0.1, 0.2),
                    seed=fold,
                )
                result = train_v2_development_fold(
                    prepared_fold.fit_rows,
                    prepared_fold.validation_rows,
                    features,
                    seed=fold,
                    holdout_fraction=selected_holdout,
                )
                validation_ap = result.validation_ap
                write_task_checkpoint(
                    checkpoint_path,
                    fingerprint,
                    {
                        "pu_holdout": selected_holdout,
                        "validation_ap": validation_ap,
                    },
                )
            if not args.resume:
                event = TrainingEvent(
                    timestamp=datetime.now(timezone.utc).isoformat(),
                    track=fingerprint.track,
                    fold=fold,
                    seed=fold,
                    model=fingerprint.model,
                    epoch=0,
                    loss=None,
                    val_ap=validation_ap,
                    lr=None,
                    device=device,
                    gpu_memory=None,
                    wall_seconds=time.monotonic() - start,
                )
            else:
                event = None
            return {
                "fingerprint": fingerprint,
                "checkpoint_path": checkpoint_path,
                "event": event,
                "fit_count": len(prepared_fold.fit_rows),
                "validation_count": len(prepared_fold.validation_rows),
                "validation_ap": validation_ap,
                "selected_holdout": selected_holdout,
            }

        task_records = run_task_group(
            tuple(range(frozen.n_development_folds)),
            max_parallel_tasks=max_parallel_tasks,
            gpu_map=gpu_map,
            default_device=default_device,
            runner=run_fold,
        )
        checkpoint_paths: dict[str, Path] = {}
        task_roster = []
        for record in task_records:
            fingerprint = cast(TaskFingerprint, record["fingerprint"])
            checkpoint_path = cast(Path, record["checkpoint_path"])
            event = cast(TrainingEvent | None, record["event"])
            checkpoint_paths[f"fold{fingerprint.fold}"] = checkpoint_path
            task_roster.append(fingerprint)
            if event is not None:
                append_training_event(log_path, event)
            print(
                f"development fold {fingerprint.fold}: fit={record['fit_count']} "
                f"validation={record['validation_count']} "
                f"val_ap={record['validation_ap']:.6f} "
                f"pu_holdout={record['selected_holdout']:.2f}"
            )
        first_fingerprint = task_roster[0]
        manifest_path = output_directory / str(output_cfg["manifest"])
        write_run_manifest(
            manifest_path,
            config_sha256=prepared.config_sha256,
            config_path=args.config,
            split_sha256=frozen.sha256,
            code_revision=code_revision,
            code_sha256=first_fingerprint.code_sha256,
            input_sha256=first_fingerprint.input_sha256,
            command=[
                sys.executable,
                *(argument for argument in sys.argv if argument != "--resume"),
            ],
            device=default_device,
            environment=collect_runtime_environment(),
            gpu_map=gpu_map,
            artifacts={"training_log": log_path},
            checkpoint_paths=checkpoint_paths,
            task_roster=tuple(task_roster),
            dirty=bool(dirty_paths),
            dirty_paths=dirty_paths,
        )
        audit_v2_run_manifest(manifest_path)
    if (
        args.prepare_literature_track
        or args.build_comparison_features
        or args.run_literature_baselines
    ):
        rows, frozen, proteomes = _development_rows(args.config, sample_unlabeled=False)
        cfg = yaml.safe_load(args.config.read_text(encoding="utf-8"))
        random_cfg = cfg["literature_random_protein"]
        seeds = tuple(random_cfg["seeds"])
        if seeds != tuple(range(10)):
            raise RuntimeError("literature random track requires seeds 0-9")
        if (
            random_cfg["test_fraction"] != 0.2
            or random_cfg["validation_fraction_of_remaining"] != 0.2
        ):
            raise RuntimeError("literature random track requires 64/16/20 geometry")
        runs = build_literature_random_track(
            rows,
            seeds=seeds,
            test_fraction=random_cfg["test_fraction"],
            validation_fraction_of_remaining=random_cfg[
                "validation_fraction_of_remaining"
            ],
            unlabeled_per_positive=cfg["unlabeled_panel"]["per_positive"],
            sampling_seed=cfg["unlabeled_panel"]["seed"],
        )
        for run in runs:
            assigned = dict(run.protein_partitions)
            assignment_counts = {
                name: sum(value == name for value in assigned.values())
                for name in ("train", "validation", "test")
            }
            panel_counts = {
                name: len({row.global_protein_id for row in run.partitions[name]})
                for name in ("train", "validation", "test")
            }
            print(
                f"literature seed {run.seed}: assigned_proteins={assignment_counts} "
                f"panel_proteins={panel_counts} "
                f"panel_sha256={run.panel_sha256}"
            )
        shared_site_keys = {
            (row.global_protein_id, row.cys_position)
            for run in runs
            for partition in run.partitions.values()
            for row in partition
        }
        if args.build_comparison_features:
            build_cfg = cfg.get("comparison_feature_build")
            if not isinstance(build_cfg, dict):
                raise RuntimeError("comparison_feature_build configuration is required")
            checkpoint = Path(str(build_cfg["esm_checkpoint"]))
            _registered(checkpoint)
            global_sequences = {
                f"{species}|{accession}": sequence
                for species, species_sequences in proteomes.items()
                for accession, sequence in species_sequences.items()
            }
            windows = build_cysteine_windows(shared_site_keys, global_sequences)
            compute_cfg = cfg.get("compute")
            if not isinstance(compute_cfg, dict):
                raise RuntimeError("compute configuration is required")
            score_batch_size = compute_cfg.get("score_batch_size")
            if not isinstance(score_batch_size, int) or score_batch_size < 1:
                raise RuntimeError("compute.score_batch_size must be positive")
            configured_batch = int(build_cfg["esm_batch_size"])

            def extract(batch_size: int) -> None:
                extract_esm2_window_artifact(
                    windows,
                    Path(str(build_cfg["esm_output_directory"])),
                    model_checkpoint=checkpoint,
                    input_sha256={
                        str(item["species"]): str(item["sha256"])
                        for item in cfg["reference_proteomes"]
                    },
                    batch_size=batch_size,
                )

            _, effective_batch, deviations = run_with_oom_batch_retry(
                extract,
                initial_batch_size=min(score_batch_size, configured_batch),
            )
            if deviations:
                print(
                    "comparison ESM OOM recovery: "
                    f"effective_batch_size={effective_batch} deviations={deviations}"
                )
            print(
                "comparison ESM artifact: "
                f"{build_cfg['esm_output_directory']} sites={len(windows)}"
            )
        if args.run_literature_baselines:
            selected_rows = {
                (row.global_protein_id, row.cys_position): row
                for run in runs
                for partition in run.partitions.values()
                for row in partition
            }
            raw_features = sequence_feature_map(selected_rows.values(), proteomes)
            if any(
                any(value is None for value in values)
                for values in raw_features.values()
            ):
                raise RuntimeError(
                    "sequence comparison features contain missing values"
                )
            features: dict[tuple[str, int], tuple[float, ...]] = {}
            for key, values in raw_features.items():
                features[key] = tuple(
                    float(value) for value in values if value is not None
                )
            inputs_cfg = cfg.get("comparison_inputs", {})
            if not isinstance(inputs_cfg, dict):
                raise RuntimeError("comparison_inputs must be a mapping")

            def optional_path(name: str) -> Path | None:
                value = inputs_cfg.get(name)
                return Path(str(value)) if value else None

            prerequisite_paths = {
                "esm_features": optional_path("esm_features"),
                "structure_features": optional_path("structure_features"),
                "sul_environment_manifest": optional_path("sul_environment_manifest"),
                "pcysmod_scores": optional_path("pcysmod_scores"),
            }
            registered_paths: set[Path] = set()
            for path in prerequisite_paths.values():
                if path is not None and path.is_file():
                    _registered(path)
                    registered_paths.add(path.resolve())
            statuses = comparator_statuses(
                **prerequisite_paths,
                registered_input_paths=frozenset(registered_paths),
            )
            status_by_model = {status.model: status for status in statuses}
            if any(
                status_by_model[model].status != "ready"
                for model in ("esm_linear_head", "structure_ranker")
            ):
                raise RuntimeError(
                    "all five direct baselines require registered ESM and "
                    "structure inputs"
                )
            esm_manifest = prerequisite_paths["esm_features"]
            structure_registry = prerequisite_paths["structure_features"]
            assert esm_manifest is not None
            assert structure_registry is not None
            esm_features = load_window_embedding_artifact(
                esm_manifest.parent,
                required_site_keys=shared_site_keys,
            )
            build_cfg = cfg.get("comparison_feature_build")
            if not isinstance(build_cfg, dict):
                raise RuntimeError("comparison_feature_build configuration is required")
            structure_features = build_registered_structure_features(
                shared_site_keys,
                registry_path=structure_registry,
                registry_base=Path(str(build_cfg["structure_registry_base"])),
            )
            reports: list[dict[str, object]] = []
            sul_manifest = prerequisite_paths["sul_environment_manifest"]
            sul_environment = None
            if status_by_model["sul_bertgru"].status == "ready":
                assert sul_manifest is not None
                sul_environment = validate_sul_environment_manifest(sul_manifest)
            output_dir = Path(cfg["output"]["directory"]) / "literature_random_protein"
            log_path = output_dir / "training.jsonl"
            checkpoint_directory = output_dir / "checkpoints"
            compute_cfg = cfg.get("compute")
            if not isinstance(compute_cfg, dict):
                raise RuntimeError("compute configuration is required")
            default_device = str(compute_cfg.get("device", "cpu"))
            raw_gpu_map = compute_cfg.get("gpu_map", [])
            if not isinstance(raw_gpu_map, list) or not all(
                isinstance(device, str) for device in raw_gpu_map
            ):
                raise RuntimeError("compute.gpu_map must be a list of device names")
            gpu_map = tuple(raw_gpu_map)
            code_revision = _code_revision()
            dirty_paths = _dirty_paths()
            dirty_code_paths = _dirty_code_paths(dirty_paths)
            runtime_inputs = {
                **_strict_runtime_input_paths(cfg),
                "esm_features": esm_manifest,
                "structure_features": structure_registry,
            }
            if prerequisite_paths["sul_environment_manifest"] is not None:
                runtime_inputs["sul_environment_manifest"] = prerequisite_paths[
                    "sul_environment_manifest"
                ]
            if sul_environment is not None:
                runtime_inputs.update(
                    {
                        "sul_python_executable": sul_environment.python_executable,
                        "sul_environment_lock": sul_environment.environment_lock,
                        "sul_adapter": sul_environment.adapter_path,
                    }
                )
                runtime_inputs.update(
                    {
                        f"sul_input_artifact_{index}": artifact
                        for index, (artifact, _sha256) in enumerate(
                            sul_environment.input_artifacts
                        )
                    }
                )
            task_roster: list[TaskFingerprint] = []
            checkpoint_paths: dict[str, Path] = {}

            def record_comparator_task(
                *,
                model: str,
                seed: int,
                report: dict[str, object],
                device: str,
                wall_seconds: float,
                deviations: list[dict[str, int]] | None = None,
                extra_input_paths: dict[str, Path] | None = None,
            ) -> None:
                fingerprint = build_task_fingerprint(
                    track="literature_random_protein",
                    fold=0,
                    seed=seed,
                    model=model,
                    input_paths={**runtime_inputs, **(extra_input_paths or {})},
                    config_path=args.config,
                    code_revision=code_revision,
                    dirty_paths=dirty_code_paths,
                )
                checkpoint_path = checkpoint_directory / f"{model}_seed{seed}.json"
                write_task_checkpoint(
                    checkpoint_path,
                    fingerprint,
                    {"report": report, "oom_batch_deviations": deviations or []},
                )
                task_roster.append(fingerprint)
                checkpoint_paths[f"{model}_seed{seed}"] = checkpoint_path
                append_training_event(
                    log_path,
                    TrainingEvent(
                        timestamp=datetime.now(timezone.utc).isoformat(),
                        track=fingerprint.track,
                        fold=fingerprint.fold,
                        seed=seed,
                        model=model,
                        epoch=0,
                        loss=None,
                        val_ap=float(report["three_crop_macro_average_precision"]),
                        lr=None,
                        device=device,
                        gpu_memory=None,
                        wall_seconds=wall_seconds,
                    ),
                )

            max_parallel_tasks = compute_cfg.get("max_parallel_tasks")
            if not isinstance(max_parallel_tasks, int) or max_parallel_tasks < 1:
                raise RuntimeError("compute.max_parallel_tasks must be positive")
            raw_sul_timeout = compute_cfg.get("sul_timeout_seconds", 21600.0)
            if not isinstance(raw_sul_timeout, int | float) or raw_sul_timeout <= 0.0:
                raise RuntimeError("compute.sul_timeout_seconds must be positive")
            sul_timeout_seconds = float(raw_sul_timeout)

            def run_direct_task(
                run: object, device: str
            ) -> tuple[object, str, object, float]:
                start = time.monotonic()
                return (
                    run,
                    device,
                    run_direct_comparison_roster(
                        run,
                        sequence_features=features,
                        esm_features=esm_features,
                        structure_features=structure_features,
                        parameters=cfg["literature_baseline_parameters"],
                        device=device,
                    ),
                    time.monotonic() - start,
                )

            def complete_direct_task(
                _task_index: int, result: tuple[object, str, object, float]
            ) -> None:
                run, assigned_device, direct_scores, wall_seconds = result
                for scores in direct_scores:
                    model_input = bind_model_to_comparison_panel(run, scores.model)
                    report = summarize_comparable_scores(
                        model_input,
                        scores,
                        primary_species=tuple(cfg["species"]["primary"]),
                        pressure_species=tuple(cfg["species"]["pressure"]),
                    )
                    reports.append(report)
                    record_comparator_task(
                        model=scores.model,
                        seed=run.seed,
                        report=report,
                        device=(
                            assigned_device
                            if scores.model == "structure_ranker"
                            else "cpu"
                        ),
                        wall_seconds=wall_seconds,
                    )
                    print(
                        f"literature seed {run.seed} model={scores.model} "
                        "three_crop_macro_ap="
                        f"{report['three_crop_macro_average_precision']:.6f}"
                    )

            run_task_group(
                tuple(runs),
                max_parallel_tasks=max_parallel_tasks,
                gpu_map=gpu_map,
                default_device=default_device,
                runner=run_direct_task,
                on_complete=complete_direct_task,
            )
            if status_by_model["sul_bertgru"].status == "ready":
                assert sul_environment is not None
                global_sequences = {
                    f"{species}|{accession}": sequence
                    for species, species_sequences in proteomes.items()
                    for accession, sequence in species_sequences.items()
                }
                output_dir = (
                    Path(cfg["output"]["directory"]) / "literature_random_protein"
                )

                def run_sul_task(
                    run: object, device: str
                ) -> tuple[object, str, object, float, list[dict[str, int]]]:
                    model_input = bind_model_to_comparison_panel(run, "sul_bertgru")
                    start = time.monotonic()
                    deviations: list[dict[str, int]] = []
                    scores = run_sul_bertgru_adapter(
                        model_input,
                        global_sequences,
                        sul_environment,
                        work_directory=output_dir / f"sul_bertgru_seed{run.seed}",
                        timeout_seconds=sul_timeout_seconds,
                        device=device,
                        oom_deviations=deviations,
                    )
                    return run, device, scores, time.monotonic() - start, deviations

                sul_results = run_task_group(
                    tuple(runs),
                    max_parallel_tasks=max_parallel_tasks,
                    gpu_map=gpu_map,
                    default_device=default_device,
                    runner=run_sul_task,
                )
                for run, device, scores, wall_seconds, deviations in sul_results:
                    model_input = bind_model_to_comparison_panel(run, "sul_bertgru")
                    report = summarize_comparable_scores(
                        model_input,
                        scores,
                        primary_species=tuple(cfg["species"]["primary"]),
                        pressure_species=tuple(cfg["species"]["pressure"]),
                    )
                    reports.append(report)
                    record_comparator_task(
                        model="sul_bertgru",
                        seed=run.seed,
                        report=report,
                        device=device,
                        wall_seconds=wall_seconds,
                        deviations=deviations,
                        extra_input_paths={
                            "shared_panel": output_dir
                            / f"sul_bertgru_seed{run.seed}"
                            / "shared_panel.tsv"
                        },
                    )
                    print(
                        f"literature seed {run.seed} model=sul_bertgru "
                        "three_crop_macro_ap="
                        f"{report['three_crop_macro_average_precision']:.6f}"
                    )
            reports.sort(key=lambda item: (int(item["seed"]), str(item["model"])))
            task_roster.sort(key=lambda item: (item.seed, item.model))
            payload = {
                "track": "literature_random_protein_development_zone",
                "config_sha256": hashlib.sha256(args.config.read_bytes()).hexdigest(),
                "code_revision": _code_revision(),
                "baseline_parameters": cfg["literature_baseline_parameters"],
                "statuses": [asdict(status) for status in statuses],
                "runs": reports,
                "limitation": "within_dataset_literature_comparable_not_gate2",
            }
            output_dir.mkdir(parents=True, exist_ok=True)
            output_path = output_dir / "summary.json"
            with tempfile.NamedTemporaryFile(
                mode="w",
                encoding="utf-8",
                dir=output_dir,
                delete=False,
            ) as handle:
                temporary = Path(handle.name)
                handle.write(json.dumps(payload, indent=2, sort_keys=True) + "\n")
            temporary.replace(output_path)
            first_fingerprint = task_roster[0]
            manifest_path = output_dir / "manifest.json"
            write_run_manifest(
                manifest_path,
                config_sha256=first_fingerprint.config_sha256,
                config_path=args.config,
                split_sha256=frozen.sha256,
                code_revision=code_revision,
                code_sha256=first_fingerprint.code_sha256,
                input_sha256=first_fingerprint.input_sha256,
                command=[sys.executable, *sys.argv],
                device=default_device,
                environment=collect_runtime_environment(),
                external_environments=(
                    {
                        "sul_bertgru": collect_runtime_environment(
                            sul_environment.python_executable
                        )
                    }
                    if sul_environment is not None
                    else {}
                ),
                gpu_map=gpu_map,
                artifacts={"training_log": log_path, "summary": output_path},
                checkpoint_paths=checkpoint_paths,
                task_roster=tuple(task_roster),
                dirty=bool(dirty_paths),
                dirty_paths=dirty_paths,
            )
            audit_v2_run_manifest(manifest_path)
            print(f"literature comparison summary: {output_path}")


if __name__ == "__main__":
    main()
