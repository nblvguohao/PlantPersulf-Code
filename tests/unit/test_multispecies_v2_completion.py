"""Remaining Task 9.3 fail-closed policy tests."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from plantpersulf.benchmark.multispecies_splits import (
    FrozenMultispeciesSplit,
    FrozenSplitRow,
    audit_frozen_multispecies_split,
    write_frozen_multispecies_split,
)
from plantpersulf.cli import main
from plantpersulf.evaluation.multispecies_ranking_audit import (
    LabelFreeRankedCysteine,
    audit_full_proteome_ranking,
)
from plantpersulf.evaluation.multispecies_reporting import (
    SpeciesMetric,
    summarize_species_metrics,
)


def _split() -> FrozenMultispeciesSplit:
    return FrozenMultispeciesSplit(
        rows=(
            FrozenSplitRow(
                "arabidopsis|P1",
                "C1",
                "development",
                0,
                "arabidopsis",
                "a" * 64,
            ),
            FrozenSplitRow("rice|P2", "C2", "test", None, "rice", "b" * 64),
        ),
        seed=20260811,
        test_fraction=0.2,
        n_development_folds=5,
        status="OK",
        test_cluster_fraction=0.5,
        positive_test_fraction_by_stratum={},
    )


def test_audit_leakage_cli_loads_requested_immutable_split(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    split_path = tmp_path / "multispecies_strict_v4.tsv"
    write_frozen_multispecies_split(split_path, _split())

    assert main([
        "audit-leakage",
        "--split-version", "multispecies_strict_v4",
        "--split-root", str(tmp_path),
    ]) == 0

    payload = json.loads(capsys.readouterr().out)
    assert payload["split_version"] == "multispecies_strict_v4"
    assert payload["protein_count"] == 2
    assert payload["test_protein_count"] == 1
    assert payload["development_fold_count"] == 5


def test_audit_leakage_cli_rejects_path_like_version(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="invalid split version"):
        main([
            "audit-leakage",
            "--split-version", "../multispecies_strict_v4",
            "--split-root", str(tmp_path),
        ])


def test_leakage_audit_rejects_unknown_partition() -> None:
    invalid = _split()
    invalid = FrozenMultispeciesSplit(
        rows=(
            FrozenSplitRow(
                "arabidopsis|P1", "C1", "discarded", None, "arabidopsis", "a" * 64
            ),
        ),
        seed=invalid.seed,
        test_fraction=invalid.test_fraction,
        n_development_folds=invalid.n_development_folds,
        status=invalid.status,
        test_cluster_fraction=invalid.test_cluster_fraction,
        positive_test_fraction_by_stratum={},
    )

    with pytest.raises(RuntimeError, match="invalid frozen split assignment"):
        audit_frozen_multispecies_split(invalid)


def _write_freeze_manifest(
    tmp_path: Path, *, model: Path, config: Path, model_sha256: str | None = None
) -> Path:
    manifest = tmp_path / "model_selection_freeze.json"
    manifest.write_text(json.dumps({
        "status": "frozen",
        "model_sha256": model_sha256 or hashlib.sha256(model.read_bytes()).hexdigest(),
        "config_sha256": hashlib.sha256(config.read_bytes()).hexdigest(),
        "code_revision": "abc123",
    }), encoding="utf-8")
    return manifest


def test_full_proteome_ranking_audit_requires_matching_frozen_model(
    tmp_path: Path,
) -> None:
    rows = (
        LabelFreeRankedCysteine("arabidopsis", "P1", 3, 0.8),
        LabelFreeRankedCysteine("magnaporthe", "M1", 7, 0.7),
    )
    model = tmp_path / "model.bin"
    config = tmp_path / "config.yaml"
    model.write_bytes(b"software-policy-model-marker")
    config.write_text("version: 2\n", encoding="utf-8")
    manifest = _write_freeze_manifest(
        tmp_path, model=model, config=config, model_sha256="0" * 64
    )

    with pytest.raises(RuntimeError, match="model artifact hash mismatch"):
        audit_full_proteome_ranking(
            rows,
            expected_cysteines=(("arabidopsis", "P1", 3), ("magnaporthe", "M1", 7)),
            model_artifact=model,
            config_path=config,
            model_freeze_manifest=manifest,
            code_revision="abc123",
            expected_primary_species=("arabidopsis", "rice", "tomato"),
            pressure_species=("magnaporthe",),
        )


def test_full_proteome_ranking_audit_is_label_free_and_separates_fungus(
    tmp_path: Path,
) -> None:
    assert "label" not in LabelFreeRankedCysteine.__dataclass_fields__
    rows = (
        LabelFreeRankedCysteine("arabidopsis", "P1", 3, 0.8),
        LabelFreeRankedCysteine("rice", "R1", 5, 0.6),
        LabelFreeRankedCysteine("tomato", "T1", 9, 0.9),
        LabelFreeRankedCysteine("magnaporthe", "M1", 7, 0.7),
    )
    model = tmp_path / "model.bin"
    config = tmp_path / "config.yaml"
    model.write_bytes(b"software-policy-model-marker")
    config.write_text("version: 2\n", encoding="utf-8")
    manifest = _write_freeze_manifest(tmp_path, model=model, config=config)

    report = audit_full_proteome_ranking(
        rows,
        expected_cysteines=tuple(
            (row.species, row.protein_accession, row.cys_position) for row in rows
        ),
        model_artifact=model,
        config_path=config,
        model_freeze_manifest=manifest,
        code_revision="abc123",
        expected_primary_species=("arabidopsis", "rice", "tomato"),
        pressure_species=("magnaporthe",),
    )

    assert set(report["three_crop_primary_ranking"]) == {
        "arabidopsis", "rice", "tomato"
    }
    assert set(report["magnaporthe_pressure_test_ranking"]) == {"magnaporthe"}
    assert "magnaporthe" not in report["three_crop_primary_ranking"]
    assert report["output_limitation"] == "label_free_ranking_audit_not_test_metrics"


def test_full_proteome_ranking_rejects_incomplete_reference_cys(
    tmp_path: Path,
) -> None:
    model = tmp_path / "model.bin"
    config = tmp_path / "config.yaml"
    model.write_bytes(b"software-policy-model-marker")
    config.write_text("version: 2\n", encoding="utf-8")
    manifest = _write_freeze_manifest(tmp_path, model=model, config=config)

    with pytest.raises(
        RuntimeError, match="ranking does not cover reference cysteines"
    ):
        audit_full_proteome_ranking(
            (LabelFreeRankedCysteine("arabidopsis", "P1", 3, 0.8),),
            expected_cysteines=(
                ("arabidopsis", "P1", 3),
                ("arabidopsis", "P1", 9),
            ),
            model_artifact=model,
            config_path=config,
            model_freeze_manifest=manifest,
            code_revision="abc123",
            expected_primary_species=("arabidopsis",),
            pressure_species=(),
        )


def test_species_metric_report_uses_non_interchangeable_panel_titles() -> None:
    def metric(species: str) -> SpeciesMetric:
        return SpeciesMetric(species, 0.5, 0.1, 0.2, 0.3, 0.4)

    report = summarize_species_metrics(
        [
            metric("arabidopsis"),
            metric("rice"),
            metric("tomato"),
            metric("magnaporthe"),
        ],
        primary_species=("arabidopsis", "rice", "tomato"),
        pressure_species=("magnaporthe",),
    )

    assert report["primary_panel_title"] == "three_crop_primary_analysis"
    assert report["pressure_panel_title"] == "magnaporthe_pressure_test"
