"""End-to-end development-only wiring for the v2 CLI entrypoint."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

from plantpersulf.benchmark.literature_random_track import ComparisonModelScores
from plantpersulf.benchmark.multispecies_splits import (
    DEVELOPMENT_SPLIT,
    FrozenMultispeciesSplit,
    FrozenSplitRow,
    GlobalClusterRow,
    MultispeciesSiteRow,
)
from plantpersulf.evaluation.comparable_track import ComparatorStatus
from plantpersulf.proteomics.multispecies_v2_dataset import MultispeciesV2SiteRow


def _load_cli_module():
    path = Path("scripts/train_multispecies_v2.py")
    spec = importlib.util.spec_from_file_location("train_multispecies_v2_test", path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_development_cli_path_passes_only_development_proteins_to_source_parser(
    tmp_path: Path, monkeypatch
) -> None:
    """A test protein must be excluded before the positive-source adapter runs."""
    cli = _load_cli_module()
    config = tmp_path / "config.yaml"
    development_manifest = tmp_path / "development_positives.tsv"
    development_manifest.write_text("policy marker\n", encoding="utf-8")
    development_manifest_sha256 = hashlib.sha256(
        development_manifest.read_bytes()
    ).hexdigest()
    config.write_text(
        "reference_proteomes:\n"
        "  - species: arabidopsis\n    path: marker.fa\n    sha256: 'a'\n"
        "  - species: tomato\n    path: marker.fa\n    sha256: 'a'\n"
        "  - species: rice\n    path: marker.fa\n    sha256: 'a'\n"
        "  - species: magnaporthe\n    path: marker.fa\n    sha256: 'a'\n"
        "positive_evidence:\n"
        "  registry: positive_registry.tsv\n"
        "  arabidopsis_benchmark: raw_positive.tsv\n"
        "  tomato_kiae271: raw_positive.tsv\n  rice_sd01: raw_positive.tsv\n"
        "  rice_sd04: raw_positive.tsv\n  rice_ss_all: raw_positive.tsv\n"
        "  magnaporthe_sites: raw_positive.tsv\n"
        "development_positive_manifest:\n"
        f"  path: {development_manifest.as_posix()}\n"
        f"  sha256: {development_manifest_sha256}\n"
        "global_mmseqs2: {cluster_table: clusters.tsv}\n"
        "strict_cluster_holdout: {split_path: split.tsv}\n"
        "unlabeled_panel: {per_positive: 1, seed: 7}\n",
        encoding="utf-8",
    )
    frozen = FrozenMultispeciesSplit(
        rows=(
            FrozenSplitRow(
                "arabidopsis|DEV", "C1", DEVELOPMENT_SPLIT, 0, "arabidopsis", "a" * 64
            ),
            FrozenSplitRow(
                "arabidopsis|TEST", "C2", "test", None, "arabidopsis", "a" * 64
            ),
        ),
        seed=1,
        test_fraction=0.2,
        n_development_folds=2,
        status="OK",
        test_cluster_fraction=0.2,
        positive_test_fraction_by_stratum={},
    )
    clusters = (
        GlobalClusterRow("arabidopsis", "DEV", "arabidopsis|DEV", "C1", 0.3, "a" * 64),
        GlobalClusterRow(
            "arabidopsis", "TEST", "arabidopsis|TEST", "C2", 0.3, "a" * 64
        ),
    )
    audited: list[Path] = []

    def frozen_source(path: Path):
        assert path == development_manifest
        return (MultispeciesSiteRow("arabidopsis", "DEV", 3, "positive", "S"),)

    monkeypatch.setattr(cli, "_registered", lambda path: audited.append(path))
    monkeypatch.setattr(cli, "load_global_cluster_table", lambda path: clusters)
    monkeypatch.setattr(cli, "load_frozen_multispecies_split", lambda path: frozen)
    monkeypatch.setattr(cli, "load_frozen_development_positives", frozen_source)
    monkeypatch.setattr(
        cli,
        "reference_cysteines",
        lambda proteomes, allowed: [
            ("arabidopsis", "DEV", 3),
            ("arabidopsis", "DEV", 5),
        ],
    )
    monkeypatch.setattr(
        "plantpersulf.features.sequence._load_proteome",
        lambda path: {"DEV": "MCMCC", "TEST": "MCMCC"},
    )
    monkeypatch.setattr(
        "plantpersulf.proteomics.pxd063170_sites.load_ensembl_fungi_proteome",
        lambda path: {},
    )

    rows, _, _ = cli._development_rows(config, sample_unlabeled=False)

    assert {row.global_protein_id for row in rows} == {"arabidopsis|DEV"}
    assert Path("clusters.tsv") in audited
    assert Path("split.tsv") in audited
    assert development_manifest in audited
    assert Path("raw_positive.tsv") not in audited
    assert not hasattr(cli, "development_positive_rows")


def test_v2_cli_exposes_development_only_literature_track() -> None:
    cli = _load_cli_module()

    arguments = cli.build_argument_parser().parse_args(
        [
            "--config",
            "configs/experiments/multispecies_v2_global_clusters_v4.yaml",
            "--prepare-literature-track",
        ]
    )

    assert arguments.prepare_literature_track is True
    assert arguments.score_test is False


def test_v2_cli_exposes_shared_panel_baseline_execution() -> None:
    cli = _load_cli_module()

    arguments = cli.build_argument_parser().parse_args(
        [
            "--config",
            "configs/experiments/multispecies_v2_global_clusters_v4.yaml",
            "--run-literature-baselines",
        ]
    )

    assert arguments.run_literature_baselines is True
    assert arguments.score_test is False


def test_v2_cli_exposes_registered_comparison_feature_build() -> None:
    cli = _load_cli_module()

    arguments = cli.build_argument_parser().parse_args(
        [
            "--config",
            "configs/experiments/multispecies_v2_global_clusters_v4.yaml",
            "--build-comparison-features",
        ]
    )

    assert arguments.build_comparison_features is True
    assert arguments.score_test is False


def test_v2_cli_exposes_hash_bound_development_resume() -> None:
    """Development resume is explicit and must not imply frozen-test scoring."""
    cli = _load_cli_module()

    arguments = cli.build_argument_parser().parse_args(
        [
            "--config",
            "configs/experiments/multispecies_v2_global_clusters_v5.yaml",
            "--prepare-development",
            "--resume",
        ]
    )

    assert arguments.prepare_development is True
    assert arguments.resume is True
    assert arguments.score_test is False


def test_cli_executes_complete_comparator_roster_on_shared_random_panel(
    tmp_path: Path, monkeypatch
) -> None:
    """The production flag must invoke all direct models and Sul, not preflight."""
    cli = _load_cli_module()
    output_dir = tmp_path / "output"
    esm_manifest = tmp_path / "esm_manifest.json"
    structure_registry = tmp_path / "structures.tsv"
    sul_manifest = tmp_path / "sul.json"
    for path in (esm_manifest, structure_registry, sul_manifest):
        path.write_text("registered policy marker\n", encoding="utf-8")
    config = tmp_path / "config.yaml"
    config.write_text(
        "species: {primary: [arabidopsis, rice, tomato], pressure: [magnaporthe]}\n"
        "literature_random_protein:\n"
        "  test_fraction: 0.2\n"
        "  validation_fraction_of_remaining: 0.2\n"
        "  seeds: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]\n"
        "unlabeled_panel: {per_positive: 1, seed: 3}\n"
        "comparison_inputs:\n"
        f"  esm_features: {esm_manifest.as_posix()}\n"
        f"  structure_features: {structure_registry.as_posix()}\n"
        f"  sul_environment_manifest: {sul_manifest.as_posix()}\n"
        "  pcysmod_scores: null\n"
        "comparison_feature_build: {structure_registry_base: data/registry}\n"
        "literature_baseline_parameters:\n"
        "  pu_logistic: {}\n"
        "  random_forest: {}\n"
        "  xgboost: {}\n"
        "  esm_linear_head: {}\n"
        "  structure_ranker: {}\n"
        f"output: {{directory: {output_dir.as_posix()}}}\n",
        encoding="utf-8",
    )
    rows = tuple(
        MultispeciesV2SiteRow(
            species="arabidopsis",
            protein_accession=f"P{index}",
            cys_position=position,
            label=label,
            study_accessions=("REAL_SOURCE",) if label == "positive" else (),
            global_protein_id=f"arabidopsis|P{index}",
            cluster_id=f"C{index}",
            split="development",
            development_fold=0,
        )
        for index in range(3)
        for position, label in ((2, "positive"), (4, "unlabeled"))
    )
    monkeypatch.setattr(
        cli,
        "run_multispecies_experiment",
        lambda *args, **kwargs: SimpleNamespace(
            cluster_count=3, test_scoring_enabled=False
        ),
    )
    monkeypatch.setattr(
        cli,
        "_development_rows",
        lambda *args, **kwargs: (
            rows,
            None,
            {"arabidopsis": {"P0": "MCAMC", "P1": "MCAMC", "P2": "MCAMC"}},
        ),
    )
    monkeypatch.setattr(
        cli,
        "sequence_feature_map",
        lambda selected, proteomes: {
            (row.global_protein_id, row.cys_position): (0.1, 0.2) for row in selected
        },
    )
    monkeypatch.setattr(cli, "_registered", lambda path: None)
    ready = tuple(
        ComparatorStatus(model, "direct_baseline", "ready", "registered")
        for model in (
            "pu_logistic",
            "random_forest",
            "xgboost",
            "esm_linear_head",
            "structure_ranker",
        )
    ) + (ComparatorStatus("sul_bertgru", "external_comparator", "ready", "registered"),)
    monkeypatch.setattr(cli, "comparator_statuses", lambda **kwargs: ready)
    monkeypatch.setattr(
        cli, "load_window_embedding_artifact", lambda *args, **kwargs: {}
    )
    monkeypatch.setattr(
        cli, "build_registered_structure_features", lambda *args, **kwargs: {}
    )
    monkeypatch.setattr(cli, "validate_sul_environment_manifest", lambda path: object())
    calls: list[tuple[int, tuple[str, ...]]] = []

    def score_models(run, **kwargs):
        calls.append((run.seed, tuple(sorted(kwargs))))
        partitions = tuple(
            (name, {(row.global_protein_id, row.cys_position): 0.5 for row in values})
            for name, values in run.partitions.items()
        )
        return tuple(
            ComparisonModelScores(model, run.seed, run.panel_sha256, partitions)
            for model in (
                "pu_logistic",
                "random_forest",
                "xgboost",
                "esm_linear_head",
                "structure_ranker",
            )
        )

    monkeypatch.setattr(cli, "run_direct_comparison_roster", score_models)
    monkeypatch.setattr(
        cli,
        "run_sul_bertgru_adapter",
        lambda model_input, *args, **kwargs: ComparisonModelScores(
            "sul_bertgru",
            model_input.seed,
            model_input.panel_sha256,
            tuple(
                (name, {(row.site_key): 0.5 for row in values})
                for name, values in model_input.partition_rows
            ),
        ),
    )
    monkeypatch.setattr(
        cli,
        "summarize_comparable_scores",
        lambda model_input, scores, **kwargs: {
            "model": scores.model,
            "seed": scores.seed,
            "panel_sha256": scores.panel_sha256,
            "three_crop_macro_average_precision": 0.5,
        },
    )

    cli.main(["--config", str(config), "--run-literature-baselines"])

    payload = json.loads(
        (output_dir / "literature_random_protein" / "summary.json").read_text()
    )
    assert len(calls) == 10
    assert {report["model"] for report in payload["runs"]} == {
        "pu_logistic",
        "random_forest",
        "xgboost",
        "esm_linear_head",
        "structure_ranker",
        "sul_bertgru",
    }
    assert len(payload["runs"]) == 60
