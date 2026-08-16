"""Unit tests for the structure-registry release snapshot helper."""

from __future__ import annotations

import csv

from scripts.accuracy_campaign.snapshot_structure_release import (
    MODEL_INPUTS_FIELDS,
    build_snapshot_path,
    model_inputs_has_input_id,
    sorted_registry_rows,
)


def test_build_snapshot_path_lowercases_input_id() -> None:
    path = build_snapshot_path("ALPHAFOLD_STRUCTURES_RELEASE_V3")
    assert path.name == "alphafold_structures_release_v3.tsv"


def test_sorted_registry_rows_sorts_by_accession() -> None:
    rows = [
        {"accession": "B1", "sha256": "a"},
        {"accession": "A2", "sha256": "b"},
    ]
    assert [row["accession"] for row in sorted_registry_rows(rows)] == ["A2", "B1"]


def test_model_inputs_has_input_id_true_and_false(tmp_path) -> None:
    import scripts.accuracy_campaign.snapshot_structure_release as module

    path = tmp_path / "model_inputs.tsv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=MODEL_INPUTS_FIELDS, delimiter="\t"
        )
        writer.writeheader()
        writer.writerow({"input_id": "EXISTING_ID"})
    monkeypatch_target = module
    original = monkeypatch_target.MODEL_INPUTS
    monkeypatch_target.MODEL_INPUTS = path
    try:
        assert model_inputs_has_input_id("EXISTING_ID") is True
        assert model_inputs_has_input_id("NEW_ID") is False
    finally:
        monkeypatch_target.MODEL_INPUTS = original
