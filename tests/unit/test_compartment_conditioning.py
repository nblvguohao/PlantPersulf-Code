"""Unit tests for the compartment conditioning overlay helpers."""

from __future__ import annotations

from plantpersulf.evaluation.compartment_conditioning import (
    cellular_component_compartments,
    compartment_overlay,
    compartment_ph,
    header_compartment,
    parse_go_field,
    primary_compartment,
)


def test_parse_go_field_returns_term_and_go_id() -> None:
    field = (
        "nucleus [GO:0005634]; "
        "DNA-binding transcription factor activity [GO:0003700]"
    )
    assert parse_go_field(field) == [
        ("nucleus", "GO:0005634"),
        ("DNA-binding transcription factor activity", "GO:0003700"),
    ]


def test_parse_go_field_empty() -> None:
    assert parse_go_field("") == []
    assert parse_go_field(None) == []


def test_cellular_component_compartments_filters_to_cc() -> None:
    terms = [
        ("nucleus", "GO:0005634"),
        ("DNA-binding transcription factor activity", "GO:0003700"),  # MF, ignored
    ]
    assert cellular_component_compartments(terms) == ["nucleus"]


def test_header_compartment_matches_specific_first() -> None:
    assert header_compartment(
        "Glucose-6-phosphate 1-dehydrogenase 6, cytoplasmic OS=..."
    ) == "cytoplasm"
    assert header_compartment(
        "GTP diphosphokinase RSH1, chloroplastic"
    ) == "chloroplast"
    assert header_compartment(
        "Respiratory burst oxidase, plasma membrane"
    ) == "membrane"


def test_header_compartment_returns_none_without_keyword() -> None:
    assert header_compartment("Cysteine protease ATG4a OS=...") is None


def test_primary_compartment_follows_priority() -> None:
    # Peroxisome outranks the generic cytoplasm when both are annotated.
    assert primary_compartment(["cytoplasm", "peroxisome"]) == "peroxisome"
    assert primary_compartment(["cytoplasm"]) == "cytoplasm"
    assert primary_compartment([]) is None


def test_compartment_ph() -> None:
    assert compartment_ph("chloroplast") == 8.0
    assert compartment_ph("vacuole") == 5.0
    assert compartment_ph("not_a_compartment") is None
    assert compartment_ph(None) is None


def test_overlay_go_takes_precedence_over_header() -> None:
    overlay = compartment_overlay(
        go_field="peroxisome [GO:0005777]",
        header="Catalase isozyme 1, cytoplasmic",
    )
    assert overlay["primary"] == "peroxisome"
    assert overlay["source"] == "go_cc"


def test_overlay_header_fallback() -> None:
    overlay = compartment_overlay(
        go_field=None,
        header="Glucose-6-phosphate 1-dehydrogenase 6, cytoplasmic",
    )
    assert overlay["primary"] == "cytoplasm"
    assert overlay["source"] == "header"


def test_overlay_not_annotated_when_neither_source() -> None:
    overlay = compartment_overlay(go_field=None, header="Cysteine protease ATG4a")
    assert overlay["primary"] is None
    assert overlay["annotated"] is False
    assert overlay["source"] == "not_annotated"
