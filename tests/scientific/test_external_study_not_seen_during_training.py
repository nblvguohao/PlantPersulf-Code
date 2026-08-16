"""RED (Task 10): leave-study-out folds never leak the held-out study.

The canonical leave-study-out partition must place *all* positives of the
held-out study in the test fold and *none* in train. Unlabeled cysteines are
the shared comparison distribution and may appear in both folds (they are not
training labels), but no held-out positive may be trained on.

Expected RED: ``plantpersulf.evaluation.external_validation`` does not exist yet.
"""

from __future__ import annotations

from plantpersulf.evaluation.external_validation import (  # RED: module missing
    partition_leave_study_out,
)


def _rows() -> list[dict[str, str]]:
    return [
        {
            "protein_accession": "P1",
            "cys_position_in_protein": "10",
            "label": "positive",
            "study_accession": "PXD_A",
        },
        {
            "protein_accession": "P2",
            "cys_position_in_protein": "20",
            "label": "positive",
            "study_accession": "PXD_B",
        },
        {
            "protein_accession": "P3",
            "cys_position_in_protein": "30",
            "label": "unlabeled",
            "study_accession": "",
        },
    ]


def test_heldout_study_positives_are_not_in_train() -> None:
    train, test = partition_leave_study_out(_rows(), holdout_study="PXD_A")

    train_positives = [r for r in train if r["label"] == "positive"]
    assert all(r["study_accession"] != "PXD_A" for r in train_positives)

    test_positives = [r for r in test if r["label"] == "positive"]
    assert {r["study_accession"] for r in test_positives} == {"PXD_A"}


def test_unlabeled_rows_are_shared_comparison_distribution() -> None:
    train, test = partition_leave_study_out(_rows(), holdout_study="PXD_A")
    train_unlabeled = [r for r in train if r["label"] != "positive"]
    test_unlabeled = [r for r in test if r["label"] != "positive"]
    # The unlabeled comparison rows are present in both folds.
    assert train_unlabeled and test_unlabeled
