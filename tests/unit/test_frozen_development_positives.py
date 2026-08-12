"""Fail-closed frozen development-positive manifest tests."""

from __future__ import annotations

import importlib.util
from pathlib import Path

from plantpersulf.benchmark.multispecies_splits import MultispeciesSiteRow
from plantpersulf.proteomics.multispecies_v2_sources import (
    load_frozen_development_positives,
    write_frozen_development_positives,
)


def test_frozen_development_positive_loader_uses_label_free_fixed_schema(
    tmp_path: Path,
) -> None:
    manifest = tmp_path / "development_positives.tsv"
    manifest.write_text(
        "species\tprotein_accession\tcys_position\tstudy_accessions\t"
        "source_sha256\n"
        f"arabidopsis\tPOLICY_P1\t3\tPOLICY_STUDY\t{'a' * 64}\n",
        encoding="utf-8",
    )

    rows = load_frozen_development_positives(manifest)

    assert len(rows) == 1
    assert rows[0].label == "positive"
    assert rows[0].protein_accession == "POLICY_P1"
    assert rows[0].study_accession == "POLICY_STUDY"


def test_frozen_development_positive_writer_is_immutable_and_merges_studies(
    tmp_path: Path,
) -> None:
    path = tmp_path / "development_positives.tsv"
    rows = (
        MultispeciesSiteRow("arabidopsis", "POLICY_P1", 3, "positive", "S2"),
        MultispeciesSiteRow("arabidopsis", "POLICY_P1", 3, "positive", "S1"),
    )

    write_frozen_development_positives(
        path,
        rows,
        source_sha256_by_study={"S1": "a" * 64, "S2": "a" * 64},
    )

    text = path.read_text(encoding="utf-8")
    assert "S1;S2" in text
    try:
        write_frozen_development_positives(
            path,
            rows,
            source_sha256_by_study={"S1": "a" * 64, "S2": "a" * 64},
        )
    except FileExistsError:
        pass
    else:
        raise AssertionError("frozen development-positive manifest was overwritten")


def test_freeze_command_is_separate_from_development_entry() -> None:
    script = Path("scripts/freeze_multispecies_v2_development_positives.py")
    spec = importlib.util.spec_from_file_location("freeze_dev_positives_test", script)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    arguments = module.build_argument_parser().parse_args(
        ["--config", "policy.yaml", "--output", "frozen.tsv"]
    )
    assert arguments.config == Path("policy.yaml")
    assert arguments.output == Path("frozen.tsv")
