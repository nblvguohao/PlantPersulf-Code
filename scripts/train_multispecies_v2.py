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
import tempfile
from dataclasses import asdict
from pathlib import Path

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
    run_multispecies_experiment,
    select_v2_development_hyperparameters,
    train_v2_development_fold,
)


def _code_revision() -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"], check=True, capture_output=True, text=True
    )
    return completed.stdout.strip()


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
        for fold in range(frozen.n_development_folds):
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
            print(
                f"development fold {fold}: fit={len(prepared_fold.fit_rows)} "
                f"validation={len(prepared_fold.validation_rows)} "
                f"val_ap={result.validation_ap:.6f} "
                f"pu_holdout={selected_holdout:.2f}"
            )
    if (
        args.prepare_literature_track
        or args.build_comparison_features
        or args.run_literature_baselines
    ):
        rows, _, proteomes = _development_rows(
            args.config, sample_unlabeled=False
        )
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
            extract_esm2_window_artifact(
                windows,
                Path(str(build_cfg["esm_output_directory"])),
                model_checkpoint=checkpoint,
                input_sha256={
                    str(item["species"]): str(item["sha256"])
                    for item in cfg["reference_proteomes"]
                },
                batch_size=int(build_cfg["esm_batch_size"]),
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
                "sul_environment_manifest": optional_path(
                    "sul_environment_manifest"
                ),
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
            for run in runs:
                direct_scores = run_direct_comparison_roster(
                    run,
                    sequence_features=features,
                    esm_features=esm_features,
                    structure_features=structure_features,
                    parameters=cfg["literature_baseline_parameters"],
                )
                for scores in direct_scores:
                    model_input = bind_model_to_comparison_panel(run, scores.model)
                    report = summarize_comparable_scores(
                        model_input,
                        scores,
                        primary_species=tuple(cfg["species"]["primary"]),
                        pressure_species=tuple(cfg["species"]["pressure"]),
                    )
                    reports.append(report)
                    print(
                        f"literature seed {run.seed} model={scores.model} "
                        "three_crop_macro_ap="
                        f"{report['three_crop_macro_average_precision']:.6f}"
                    )
            sul_manifest = prerequisite_paths["sul_environment_manifest"]
            if status_by_model["sul_bertgru"].status == "ready":
                assert sul_manifest is not None
                sul_environment = validate_sul_environment_manifest(sul_manifest)
                global_sequences = {
                    f"{species}|{accession}": sequence
                    for species, species_sequences in proteomes.items()
                    for accession, sequence in species_sequences.items()
                }
                output_dir = (
                    Path(cfg["output"]["directory"])
                    / "literature_random_protein"
                )
                for run in runs:
                    model_input = bind_model_to_comparison_panel(run, "sul_bertgru")
                    scores = run_sul_bertgru_adapter(
                        model_input,
                        global_sequences,
                        sul_environment,
                        work_directory=output_dir / f"sul_bertgru_seed{run.seed}",
                    )
                    report = summarize_comparable_scores(
                        model_input,
                        scores,
                        primary_species=tuple(cfg["species"]["primary"]),
                        pressure_species=tuple(cfg["species"]["pressure"]),
                    )
                    reports.append(report)
                    print(
                        f"literature seed {run.seed} model=sul_bertgru "
                        "three_crop_macro_ap="
                        f"{report['three_crop_macro_average_precision']:.6f}"
                    )
            payload = {
                "track": "literature_random_protein_development_zone",
                "config_sha256": hashlib.sha256(args.config.read_bytes()).hexdigest(),
                "code_revision": _code_revision(),
                "baseline_parameters": cfg["literature_baseline_parameters"],
                "statuses": [asdict(status) for status in statuses],
                "runs": reports,
                "limitation": "within_dataset_literature_comparable_not_gate2",
            }
            output_dir = Path(cfg["output"]["directory"]) / "literature_random_protein"
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
            print(f"literature comparison summary: {output_path}")


if __name__ == "__main__":
    main()
