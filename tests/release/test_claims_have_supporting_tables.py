"""RED (Task 10): every outward claim must cite an existing supporting table.

A claim in a collaboration/outward document may only be made if the result
table it cites actually exists on disk. This prevents "predictive value"
statements that are not backed by a shipped, reproducible table.

Expected RED: ``plantpersulf.evaluation.conclusion_gate`` does not exist yet.
"""

from __future__ import annotations

from plantpersulf.evaluation.conclusion_gate import (  # RED: module missing
    find_unsupported_claims,
)


def test_claim_citing_existing_table_is_supported() -> None:
    claims = [
        {"text": "Leave-study-out AP is reported", "table": "metrics.tsv"},
    ]
    available = {"metrics.tsv", "external_validation.json"}
    assert find_unsupported_claims(claims, available) == []


def test_claim_citing_missing_table_is_flagged() -> None:
    claims = [
        {"text": "Model predicts novel sites", "table": "predictions.tsv"},
        {"text": "Recovery table", "table": "control_recovery.tsv"},
    ]
    available = {"control_recovery.tsv"}
    flagged = find_unsupported_claims(claims, available)
    assert [c["table"] for c in flagged] == ["predictions.tsv"]


def test_claim_without_any_table_is_flagged() -> None:
    claims = [{"text": "Bold unbacked claim", "table": ""}]
    assert find_unsupported_claims(claims, {"metrics.tsv"})
