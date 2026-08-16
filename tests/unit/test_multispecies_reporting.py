"""Reporting rules that prevent fungal data or weak tracks inflating claims."""

from __future__ import annotations

import pytest

from plantpersulf.evaluation.multispecies_reporting import (
    HOMOLOGY_ROBUST_CLAIM,
    LITERATURE_COMPARABLE_CLAIM,
    NO_SUPPORTED_CLAIM,
    STRICT_BOOTSTRAP_N_BOOT,
    STRICT_BOOTSTRAP_SEED,
    SpeciesMetric,
    build_multispecies_statistical_report,
    claim_class_for_multispecies_result,
    claim_statement_for_class,
    compute_species_metric,
    find_forbidden_external_claims,
    literature_track_leads,
    pooled_average_precision,
    strict_track_bootstrap_delta,
    summarize_species_metrics,
)


def test_primary_macro_excludes_fungal_pressure_species() -> None:
    """Including the fungal value would change the advertised three-crop result."""
    report = summarize_species_metrics(
        [
            SpeciesMetric("arabidopsis", 0.10, 0.05, 0.1, 0.2, 0.3),
            SpeciesMetric("rice", 0.20, 0.05, 0.1, 0.2, 0.3),
            SpeciesMetric("tomato", 0.30, 0.05, 0.1, 0.2, 0.3),
            SpeciesMetric("magnaporthe", 0.90, 0.05, 0.1, 0.2, 0.3),
        ],
        primary_species=("arabidopsis", "rice", "tomato"),
        pressure_species=("magnaporthe",),
    )

    assert report["primary_macro_average_precision"] == 0.2
    assert report["pressure_species"]["magnaporthe"]["average_precision"] == 0.9


def test_claim_requires_strict_track_ci_and_all_three_positive_deltas() -> None:
    """A random-split win must not be presented as homology-robust evidence."""
    assert claim_class_for_multispecies_result(
        strict_ci_lower=0.01,
        species_deltas={"arabidopsis": 0.01, "rice": 0.02, "tomato": 0.03},
        same_frozen_inputs=True,
        literature_leads=False,
    ) == HOMOLOGY_ROBUST_CLAIM
    assert claim_class_for_multispecies_result(
        strict_ci_lower=-0.01,
        species_deltas={"arabidopsis": 0.01, "rice": 0.02, "tomato": 0.03},
        same_frozen_inputs=True,
        literature_leads=True,
    ) == LITERATURE_COMPARABLE_CLAIM


def test_claim_gate_returns_no_supported_claim_without_either_track_win() -> None:
    """Previously this silently fell back to the literature-improvement claim
    even though the literature track never won anything — a fail-open bug.
    With neither track's evidence clearing its bar, no positive claim is
    allowed at all."""
    assert claim_class_for_multispecies_result(
        strict_ci_lower=-0.02,
        species_deltas={"arabidopsis": -0.01, "rice": 0.0, "tomato": -0.02},
        same_frozen_inputs=True,
        literature_leads=False,
    ) == NO_SUPPORTED_CLAIM


def test_claim_statement_matches_mandated_wording() -> None:
    assert "文献同口径" in claim_statement_for_class(LITERATURE_COMPARABLE_CLAIM)
    assert "对未见同源家族稳健" in claim_statement_for_class(HOMOLOGY_ROBUST_CLAIM)
    assert "两条轨道均未证明" in claim_statement_for_class(NO_SUPPORTED_CLAIM)
    with pytest.raises(ValueError, match="unknown"):
        claim_statement_for_class("made_up_class")


def test_forbidden_external_claim_phrases_are_flagged() -> None:
    """No blind tomato cohort exists yet, so these phrases are unconditionally
    forbidden in any multispecies v2 outward text."""
    assert find_forbidden_external_claims("这是外部泛化结果") == ["外部泛化"]
    assert find_forbidden_external_claims("within-dataset improvement only") == []


def test_statistical_report_composes_species_pooled_and_claim() -> None:
    """The Task 9.6 orchestrator must wire raw scored rows all the way
    through to one claim class and its mandated wording."""
    scored_by_species = {
        "arabidopsis": [(0.9, "positive"), (0.1, "unlabeled")],
        "rice": [(0.8, "positive"), (0.2, "unlabeled")],
        "tomato": [(0.7, "positive"), (0.3, "unlabeled")],
        "magnaporthe": [(0.6, "positive"), (0.4, "unlabeled")],
    }

    report = build_multispecies_statistical_report(
        scored_by_species=scored_by_species,
        primary_species=("arabidopsis", "rice", "tomato"),
        pressure_species=("magnaporthe",),
        strict_ci_lower=0.01,
        species_deltas={"arabidopsis": 0.01, "rice": 0.02, "tomato": 0.03},
        same_frozen_inputs=True,
        candidate_literature_macro_aps=[0.5, 0.6],
        baseline_literature_macro_aps=[0.4, 0.4],
    )

    assert report["claim_class"] == HOMOLOGY_ROBUST_CLAIM
    assert report["claim_statement"] == claim_statement_for_class(
        HOMOLOGY_ROBUST_CLAIM
    )
    assert report["gate2_eligible"] is False
    assert report["literature_track_leads"] is True
    assert "primary_pooled_average_precision" in report
    assert report["primary_macro_average_precision"] == pytest.approx(1.0)


def test_compute_species_metric_matches_underlying_pu_metrics() -> None:
    """The frozen-test report must derive Recall@K/MRR from the same
    PU-appropriate metric functions used everywhere else, not reimplement
    them ad hoc."""
    scored = [
        (0.9, "positive"),
        (0.8, "unlabeled"),
        (0.7, "positive"),
        (0.1, "unlabeled"),
    ]

    metric = compute_species_metric(scored, "arabidopsis")

    assert metric.species == "arabidopsis"
    assert metric.base_rate == 0.5
    assert metric.recall_at_50 == 1.0
    assert metric.mean_reciprocal_rank == 1.0


def test_pooled_average_precision_differs_from_macro_average() -> None:
    """Pooling ranks across species is not the same statistic as macro-
    averaging their independent APs; the report must expose both."""
    scored_by_species = {
        "arabidopsis": [(0.9, "positive"), (0.1, "unlabeled")],
        "rice": [(0.2, "positive")] + [(0.8, "unlabeled")] * 9,
    }

    pooled = pooled_average_precision(scored_by_species, ("arabidopsis", "rice"))

    assert 0.0 < pooled < 1.0


def test_literature_track_leads_compares_paired_seed_means() -> None:
    """The literature-track lead check is a plain paired mean over the same
    seeds, not a CI — the spec only mandates a CI for the strict track."""
    assert literature_track_leads([0.5, 0.6], [0.4, 0.4]) is True
    assert literature_track_leads([0.3, 0.3], [0.5, 0.5]) is False
    with pytest.raises(ValueError, match="equal length"):
        literature_track_leads([0.5], [0.4, 0.4])


def test_strict_track_bootstrap_delta_uses_frozen_policy() -> None:
    """The strict-track claim depends on exactly this bootstrap policy; a
    silently different n_boot/seed would make the CI non-reproducible."""
    model = [(0.9, "positive", "C1"), (0.2, "unlabeled", "C2")]
    baseline = [(0.6, "positive", "C1"), (0.5, "unlabeled", "C2")]

    result = strict_track_bootstrap_delta(model, baseline)

    assert STRICT_BOOTSTRAP_N_BOOT == 10_000
    assert STRICT_BOOTSTRAP_SEED == 20260811
    assert result.n_boot == STRICT_BOOTSTRAP_N_BOOT
