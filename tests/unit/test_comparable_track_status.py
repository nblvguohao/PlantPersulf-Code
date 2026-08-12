"""Task 9.4 comparator eligibility and fail-closed status tests."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest

from plantpersulf.benchmark.literature_random_track import (
    ComparisonModelInput,
    ComparisonModelScores,
    ComparisonSite,
)
from plantpersulf.evaluation.comparable_track import (
    SulEnvironment,
    comparator_statuses,
    load_complete_external_scores,
    run_sul_bertgru_adapter,
    summarize_comparable_scores,
    validate_sul_environment_manifest,
)


def test_missing_registered_features_block_only_dependent_comparators() -> None:
    statuses = {
        status.model: status
        for status in comparator_statuses(
            esm_features=None,
            structure_features=None,
            sul_environment_manifest=None,
            pcysmod_scores=None,
        )
    }

    assert statuses["pu_logistic"].status == "ready"
    assert statuses["random_forest"].status == "ready"
    assert statuses["xgboost"].status == "ready"
    assert statuses["esm_linear_head"].status == "blocked"
    assert statuses["structure_ranker"].status == "blocked"
    assert statuses["sul_bertgru"].status == "blocked"
    assert statuses["pcysmod"].status == "qualitative_only"
    assert statuses["tree"].status == "architecture_reference_only"
    assert statuses["graft"].status == "architecture_reference_only"


def test_missing_inputs_are_never_replaced_by_sequence_features() -> None:
    statuses = {
        status.model: status
        for status in comparator_statuses(
            esm_features=None,
            structure_features=None,
            sul_environment_manifest=None,
            pcysmod_scores=None,
        )
    }

    assert statuses["esm_linear_head"].reason == "missing_registered_esm_features"
    assert statuses["structure_ranker"].reason == (
        "missing_registered_structure_features"
    )


def test_comparable_report_separates_three_crops_from_fungal_pressure() -> None:
    species = ("arabidopsis", "rice", "tomato", "magnaporthe")
    test_rows = tuple(
        ComparisonSite(name, f"{name}_P", position, label, ())
        for name in species
        for position, label in ((3, "positive"), (7, "unlabeled"))
    )
    model_input = ComparisonModelInput(
        model="pu_logistic",
        seed=0,
        panel_sha256="a" * 64,
        partition_rows=(
            ("train", test_rows),
            ("validation", test_rows),
            ("test", test_rows),
        ),
    )
    scores = {row.site_key: float(row.label == "positive") for row in test_rows}
    model_scores = ComparisonModelScores(
        model="pu_logistic",
        seed=0,
        panel_sha256="a" * 64,
        partition_scores=(("validation", scores), ("test", scores)),
    )

    report = summarize_comparable_scores(
        model_input,
        model_scores,
        primary_species=("arabidopsis", "rice", "tomato"),
        pressure_species=("magnaporthe",),
    )

    assert set(report["three_crop_primary_test"]) == {
        "arabidopsis",
        "rice",
        "tomato",
    }
    assert set(report["magnaporthe_pressure_test"]) == {"magnaporthe"}
    assert "magnaporthe" not in report["three_crop_primary_test"]
    assert set(report["three_crop_primary_validation"]) == {
        "arabidopsis",
        "rice",
        "tomato",
    }


def test_pcysmod_import_rejects_partial_frozen_test_scores(tmp_path: Path) -> None:
    rows = (
        ComparisonSite("arabidopsis", "P1", 3, "positive", ()),
        ComparisonSite("arabidopsis", "P1", 7, "unlabeled", ()),
    )
    model_input = ComparisonModelInput(
        model="pcysmod",
        seed=0,
        panel_sha256="b" * 64,
        partition_rows=(("train", rows), ("validation", rows), ("test", rows)),
    )
    path = tmp_path / "pcysmod.tsv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("global_protein_id", "cys_position", "score"))
        writer.writerow(("arabidopsis|P1", 3, 0.9))

    with pytest.raises(RuntimeError, match="complete frozen input"):
        load_complete_external_scores(path, model_input, partition="test")


def test_sul_environment_must_use_a_separate_python(tmp_path: Path) -> None:
    lock = tmp_path / "environment.lock"
    adapter = tmp_path / "adapter.py"
    lock.write_text("policy lock\n", encoding="utf-8")
    adapter.write_text("# policy adapter\n", encoding="utf-8")
    manifest = tmp_path / "sul_environment.json"
    manifest.write_text(
        "{\n"
        '  "isolated": true,\n'
        f'  "python_executable": "{Path(sys.executable).as_posix()}",\n'
        f'  "environment_lock": "{lock.as_posix()}",\n'
        f'  "adapter_path": "{adapter.as_posix()}",\n'
        '  "repo_commit": "52c030c3b23e20ff8170d73a417ce852bd46c627"\n'
        "}\n",
        encoding="utf-8",
    )

    with pytest.raises(RuntimeError, match="separate Python environment"):
        validate_sul_environment_manifest(manifest)


def test_sul_adapter_returns_complete_shared_panel_scores(tmp_path: Path) -> None:
    rows = (
        ComparisonSite("arabidopsis", "P1", 3, "positive", ()),
        ComparisonSite("arabidopsis", "P1", 7, "unlabeled", ()),
    )
    model_input = ComparisonModelInput(
        model="sul_bertgru",
        seed=0,
        panel_sha256="c" * 64,
        partition_rows=(("train", rows), ("validation", rows), ("test", rows)),
    )
    adapter = tmp_path / "adapter.py"
    adapter.write_text(
        "import argparse, csv\n"
        "p=argparse.ArgumentParser(); p.add_argument('--input'); "
        "p.add_argument('--output'); p.add_argument('--seed'); "
        "p.add_argument('--epochs', required=True); a=p.parse_args()\n"
        "rows=list(csv.DictReader(open(a.input), delimiter='\\t'))\n"
        "f=open(a.output,'w',newline=''); "
        "w=csv.writer(f,delimiter='\\t',lineterminator='\\n'); "
        "w.writerow(['partition','global_protein_id','cys_position','score'])\n"
        "[w.writerow([r['partition'],r['global_protein_id'],"
        "r['cys_position'],str(float(a.epochs)/2)]) "
        "for r in rows if r['partition']!='train']; "
        "f.close()\n",
        encoding="utf-8",
    )
    audited = SulEnvironment(
        python_executable=Path(sys.executable),
        environment_lock=tmp_path / "lock",
        adapter_path=adapter,
        repo_commit="52c030c3b23e20ff8170d73a417ce852bd46c627",
        parameters=(("epochs", "1"),),
    )
    sequences = {"arabidopsis|P1": "MMC" + "AAA" + "C" + "M" * 20}

    scores = run_sul_bertgru_adapter(
        model_input,
        sequences,
        audited,
        work_directory=tmp_path / "run",
    )

    assert set(dict(scores.partition_scores)) == {"validation", "test"}
