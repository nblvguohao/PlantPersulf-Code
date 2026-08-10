"""RED (Task 10): integrity rules for known-mechanism control recovery.

Covers the Codex "future RED" rules that do not require a fictional site — they
operate on already-registered control records and pure software policy:

* a control cannot be both a training positive and a reported independent
  recovery (no self-recovery);
* duplicate ``mechanism_lineage_id`` controls collapse to one independent unit;
* model/ablation selection cannot see control ranks (structural);
* failed and unmappable controls must appear in the recovery report.

Expected RED: ``plantpersulf.evaluation.external_validation`` does not exist yet.
"""

from __future__ import annotations

import inspect

import pytest

from plantpersulf.evaluation.external_validation import (  # RED: module missing
    ControlRecord,
    build_recovery_report,
    count_independent_validation_units,
    find_control_training_leakage,
    select_ablation_by_validation,
)


def _mapped(lineage: str, acc: str, pos: int) -> ControlRecord:
    return ControlRecord(
        mechanism_lineage_id=lineage,
        gene=lineage,
        uniprot_accession=acc,
        cys_position=pos,
        status="mapped",
        percentile_rank=90.0,
    )


def test_control_that_is_a_training_positive_is_flagged() -> None:
    controls = [_mapped("LIN_A", "Q1", 5), _mapped("LIN_B", "Q2", 9)]
    training_positives = {("Q1", 5)}  # LIN_A leaked into training

    leaks = find_control_training_leakage(controls, training_positives)

    assert [c.mechanism_lineage_id for c in leaks] == ["LIN_A"]


def test_no_leakage_when_controls_are_cross_species_holdouts() -> None:
    controls = [_mapped("LIN_A", "Q1", 5)]
    training_positives = {("ARABIDOPSIS1", 12)}
    assert find_control_training_leakage(controls, training_positives) == []


def test_duplicate_mechanism_lineage_is_not_counted_twice() -> None:
    controls = [
        _mapped("LIN_A", "Q1", 5),
        _mapped("LIN_A", "Q1b", 6),  # same lineage, different mapping
        _mapped("LIN_B", "Q2", 9),
    ]
    assert count_independent_validation_units(controls) == 2


def test_unmappable_control_is_not_an_independent_unit() -> None:
    controls = [
        _mapped("LIN_A", "Q1", 5),
        ControlRecord(
            mechanism_lineage_id="LIN_C",
            gene="BRG3",
            uniprot_accession="",
            cys_position=0,
            status="unmappable",
            percentile_rank=None,
        ),
    ]
    assert count_independent_validation_units(controls) == 1


def test_selection_cannot_read_control_ranks() -> None:
    # Structural guarantee: the selection function's signature has no control
    # parameter — it is impossible to select an ablation using control ranks.
    params = set(inspect.signature(select_ablation_by_validation).parameters)
    assert not any("control" in p or "rank" in p for p in params)

    best = select_ablation_by_validation({"full": 0.4, "seq_only": 0.1})
    assert best == "full"


def test_failed_and_unmappable_controls_are_reported() -> None:
    controls = [
        _mapped("LIN_A", "Q1", 5),
        ControlRecord(
            mechanism_lineage_id="LIN_C",
            gene="BRG3",
            uniprot_accession="",
            cys_position=0,
            status="unmappable",
            percentile_rank=None,
        ),
        ControlRecord(
            mechanism_lineage_id="LIN_D",
            gene="X",
            uniprot_accession="Q9",
            cys_position=3,
            status="feature_extraction_failed",
            percentile_rank=None,
        ),
    ]
    report = build_recovery_report(controls)
    reported_statuses = {row["status"] for row in report}
    assert "unmappable" in reported_statuses
    assert "feature_extraction_failed" in reported_statuses
    # every control appears exactly once
    assert len(report) == len(controls)


def test_report_rejects_empty_controls() -> None:
    with pytest.raises(ValueError):
        build_recovery_report([])
