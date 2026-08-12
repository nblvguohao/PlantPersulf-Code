"""Reporting rules that prevent fungal data or weak tracks inflating claims."""

from __future__ import annotations

from plantpersulf.evaluation.multispecies_reporting import (
    SpeciesMetric,
    claim_class_for_multispecies_result,
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
    ) == "homology_robust_multiplant_within_dataset_ranking"
    assert claim_class_for_multispecies_result(
        strict_ci_lower=-0.01,
        species_deltas={"arabidopsis": 0.01, "rice": 0.02, "tomato": 0.03},
        same_frozen_inputs=True,
    ) == "literature_comparable_within_dataset_improvement_only"
