"""Task 10 — strict external-validation aggregation and integrity rules.

This module holds the *logic* the external-validation report is built from,
factored out of the runner and the control-evaluation script so it is unit
testable without heavy data:

* ``partition_leave_study_out`` — the canonical no-leakage fold split (all
  held-out-study positives to test, none to train; unlabeled shared).
* ``ControlRecord`` + ``find_control_training_leakage`` — a mapped known-
  mechanism control may never be both a training positive and a reported
  independent recovery.
* ``count_independent_validation_units`` — duplicate ``mechanism_lineage_id``
  controls collapse to one unit; unmappable/failed controls count for none.
* ``select_ablation_by_validation`` — model/ablation selection whose signature
  structurally cannot see control ranks (validation metrics only).
* ``build_recovery_report`` — every control is reported, including failed and
  unmappable ones (never only the successes).

The higher-level ``build_external_validation`` stitches per-study
leave-study-out metrics, a cluster bootstrap CI, a permutation p-value and the
control recovery report into one serialisable structure, recording all seeds.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

# ---------------------------------------------------------------------------
# leave-study-out partition (canonical, reused by the runner)
# ---------------------------------------------------------------------------

Row = dict[str, str]


def partition_leave_study_out(
    rows: list[Row],
    holdout_study: str,
) -> tuple[list[Row], list[Row]]:
    """Split benchmark rows for one leave-study-out fold.

    Positives of ``holdout_study`` go to test only; positives of other studies
    go to train only. Unlabeled rows are the shared comparison distribution and
    are appended to both folds (they are never training labels, so this leaks
    nothing).
    """
    train: list[Row] = []
    test: list[Row] = []
    unlabeled: list[Row] = []
    for row in rows:
        if row.get("label") != "positive":
            unlabeled.append(row)
        elif row.get("study_accession") == holdout_study:
            test.append(row)
        else:
            train.append(row)
    train.extend(unlabeled)
    test.extend(unlabeled)
    return train, test


def partition_random_protein(
    rows: list[Row],
    seed: int,
    test_ratio: float = 0.2,
) -> tuple[list[Row], list[Row]]:
    """Random protein-level split — literature-comparable regime, NOT
    cross-study validation.

    Mirrors the evaluation regime of published cysteine-PTM predictors
    (Sul-BertGRU, Bioinformatics 2025, btaf078: 20% of proteins held out as
    an independent test set, 10 repetitions). ``test_ratio`` of *proteins*
    are held out entirely — all their Cys rows (positive and unlabeled) go
    to test; the remaining proteins go to train. Deterministic per ``seed``.

    Deliberately NO homology control: this matches the published regime so
    the resulting numbers are comparable with the literature. It therefore
    does NOT cross a study, laboratory, chemistry, or species boundary, and
    its output must never enter Gate 2 — scripts/validate_external.py only
    parses ``leave_<study>_out``-tagged model names, and the runner tags
    this track's folds ``protein_split_seed<k>|...`` so they are structurally
    invisible to Gate 2 (locked in by
    tests/release/test_gate2_ignores_within_dataset_split_metrics.py).
    """
    import random

    proteins = sorted({r["protein_accession"] for r in rows})
    rng = random.Random(seed)
    rng.shuffle(proteins)
    n_test = max(1, round(len(proteins) * test_ratio))
    test_proteins = set(proteins[:n_test])
    train: list[Row] = []
    test: list[Row] = []
    for row in rows:
        if row["protein_accession"] in test_proteins:
            test.append(row)
        else:
            train.append(row)
    return train, test


# ---------------------------------------------------------------------------
# known-mechanism control integrity
# ---------------------------------------------------------------------------

RECOVERED_STATUS = "mapped"


@dataclass(frozen=True)
class ControlRecord:
    mechanism_lineage_id: str
    gene: str
    uniprot_accession: str
    cys_position: int
    status: str
    percentile_rank: float | None = None


def find_control_training_leakage(
    controls: list[ControlRecord],
    training_positives: set[tuple[str, int]],
) -> list[ControlRecord]:
    """Return mapped controls whose (accession, position) is a training
    positive — such a control cannot be an independent recovery."""
    return [
        c
        for c in controls
        if c.status == RECOVERED_STATUS
        and (c.uniprot_accession, c.cys_position) in training_positives
    ]


def count_independent_validation_units(controls: list[ControlRecord]) -> int:
    """Number of distinct mechanism lineages among *recovered* controls.

    Duplicate lineages collapse to one; unmappable/failed controls do not
    count as independent validation units.
    """
    lineages = {
        c.mechanism_lineage_id for c in controls if c.status == RECOVERED_STATUS
    }
    return len(lineages)


def select_ablation_by_validation(validation_ap: dict[str, float]) -> str:
    """Pick the ablation with the best *validation* AP.

    The signature deliberately admits only validation metrics — there is no
    control-rank parameter, so selection cannot be tuned on control recovery.
    """
    if not validation_ap:
        raise ValueError("validation_ap must not be empty")
    return max(sorted(validation_ap), key=lambda name: validation_ap[name])


def build_recovery_report(controls: list[ControlRecord]) -> list[dict[str, Any]]:
    """One row per control — successes, failures, and unmappables alike."""
    if not controls:
        raise ValueError("controls must not be empty")
    return [
        {
            "mechanism_lineage_id": c.mechanism_lineage_id,
            "gene": c.gene,
            "uniprot_accession": c.uniprot_accession,
            "cys_position": c.cys_position,
            "status": c.status,
            "percentile_rank": c.percentile_rank,
            "is_independent_unit": c.status == RECOVERED_STATUS,
        }
        for c in controls
    ]


# ---------------------------------------------------------------------------
# top-level aggregation
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class StudyFoldMetrics:
    holdout_study: str
    model: str
    seeds: list[int]
    test_ap: list[float]
    baseline_test_ap: list[float]


@dataclass(frozen=True)
class ExternalValidation:
    release: str
    fold_metrics: list[StudyFoldMetrics] = field(default_factory=list)
    bootstrap: dict[str, Any] = field(default_factory=dict)
    permutation: dict[str, Any] = field(default_factory=dict)
    control_report: list[dict[str, Any]] = field(default_factory=list)
    independent_units: int = 0
    control_leakage: list[str] = field(default_factory=list)
    limitation: str = ""
    effect: dict[str, Any] = field(default_factory=dict)
    structure_gain: dict[str, Any] = field(default_factory=dict)
    cluster_dominance: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "release": self.release,
            "fold_metrics": [
                {
                    "holdout_study": f.holdout_study,
                    "model": f.model,
                    "seeds": f.seeds,
                    "test_ap": f.test_ap,
                    "baseline_test_ap": f.baseline_test_ap,
                }
                for f in self.fold_metrics
            ],
            "bootstrap": self.bootstrap,
            "permutation": self.permutation,
            "control_report": self.control_report,
            "independent_units": self.independent_units,
            "control_leakage": self.control_leakage,
            "limitation": self.limitation,
            "effect": self.effect,
            "structure_gain": self.structure_gain,
            "cluster_dominance": self.cluster_dominance,
        }
