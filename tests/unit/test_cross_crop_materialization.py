"""RED tests for Task 9A's deterministic, PU-safe materialization policy."""

import hashlib
from pathlib import Path

import pytest

from plantpersulf.workflows.cross_crop_materialization import (
    MaterializationSite,
    _load_registered_inputs,
    _source_rows_from_proteome,
    build_source_payload,
    build_target_payload,
    deterministic_pu_source_sample,
    validate_label_free_target_payload,
)


def test_source_sampling_is_order_invariant_and_preserves_pu_labels() -> None:
    """A frozen hash seed selects real unlabeled rows without relabelling them."""
    rows = (
        MaterializationSite("marker-a", 1, "positive"),
        MaterializationSite("marker-b", 2, "positive"),
        MaterializationSite("marker-c", 3, "unlabeled"),
        MaterializationSite("marker-d", 4, "unlabeled"),
        MaterializationSite("marker-e", 5, "unlabeled"),
        MaterializationSite("marker-f", 6, "unlabeled"),
        MaterializationSite("marker-g", 7, "unlabeled"),
    )

    first = deterministic_pu_source_sample(rows, unlabeled_per_positive=2, seed=17)
    second = deterministic_pu_source_sample(
        tuple(reversed(rows)), unlabeled_per_positive=2, seed=17
    )

    assert first == second
    assert [row.label for row in first].count("positive") == 2
    assert [row.label for row in first].count("unlabeled") == 4
    assert {row.label for row in first} == {"positive", "unlabeled"}


def test_target_payload_rejects_any_target_label_field() -> None:
    try:
        validate_label_free_target_payload(
            {"site_keys": ["marker-a:C1"], "features": [[0.0]], "labels": []}
        )
    except RuntimeError as exc:
        assert "target label" in str(exc)
    else:
        raise AssertionError("target labels must fail closed")


def test_source_payload_has_only_contract_features_and_never_negative_labels() -> None:
    rows = (
        MaterializationSite("marker-a", 2, "positive"),
        MaterializationSite("marker-b", 2, "unlabeled"),
        MaterializationSite("marker-c", 2, "unlabeled"),
    )
    payload = build_source_payload(
        species="Marker species",
        study_accession="MARKER001",
        rows=rows,
        proteome={"marker-a": "ACA", "marker-b": "ACA", "marker-c": "ACA"},
        source_input_sha256s=("marker-input-sha",),
        unlabeled_per_positive=2,
        seed=17,
        panther_subfamilies={},
    )

    assert payload["feature_names"] == [
        "hydrophobicity",
        "protein_cys_density",
        "local_positive_charge_density",
        "local_negative_charge_density",
        "local_cys_density",
        "local_sequence_entropy",
    ]
    assert payload["labels"] == ["positive", "unlabeled", "unlabeled"]
    assert "negative" not in payload["labels"]
    assert "cluster_ids" not in payload  # missing orthology calls fall back to protein.


def test_target_payload_enumerates_all_reference_cysteines_without_labels() -> None:
    payload = build_target_payload(
        {"tomato-marker-a": "ACA", "tomato-marker-b": "CCC"},
        source_sha256="marker-fasta-sha",
    )

    assert payload["site_keys"] == [
        "tomato-marker-a:C2",
        "tomato-marker-b:C1",
        "tomato-marker-b:C2",
        "tomato-marker-b:C3",
    ]
    assert "labels" not in payload


def test_supplemental_isoform_sequence_supplies_only_its_observed_positive() -> None:
    rows = _source_rows_from_proteome(
        {"canonical-marker": "ACA"},
        {("isoform-marker", 2)},
        supplemental_proteome={"isoform-marker": "ACACA"},
    )

    assert rows == (
        MaterializationSite("canonical-marker", 2, "unlabeled"),
        MaterializationSite("isoform-marker", 2, "positive"),
    )


def test_input_registry_rejects_duplicate_ids_before_a_source_is_read(
    tmp_path: Path,
) -> None:
    marker = tmp_path / "marker.txt"
    marker.write_text("non-biological marker", encoding="utf-8")
    digest = hashlib.sha256(marker.read_bytes()).hexdigest()
    registry = tmp_path / "inputs.tsv"
    row = (
        "MARKER\tmarker\tpolicy_test\tmarker.txt\t"
        f"{digest}\tnon-biological marker\n"
    )
    registry.write_text(
        "input_id\tsource_accession\tsource_kind\tpath\tsha256\tscientific_use\n"
        + row
        + row,
        encoding="utf-8",
    )

    with pytest.raises(RuntimeError, match="duplicate"):
        _load_registered_inputs(registry)
