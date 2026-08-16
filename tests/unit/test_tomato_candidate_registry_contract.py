"""Contract tests for the label-free tomato candidate registry."""

from dataclasses import fields

from plantpersulf.proteomics.tomato_candidate_registry import TomatoCandidate


def test_candidate_contract_has_no_label_field() -> None:
    assert tuple(field.name for field in fields(TomatoCandidate)) == (
        "protein_accession",
        "cys_position",
        "site_key",
    )
