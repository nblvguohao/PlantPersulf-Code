"""Species-stratified reporting and claim gates for multispecies v2."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class SpeciesMetric:
    species: str
    average_precision: float
    base_rate: float
    recall_at_50: float
    recall_at_200: float
    mean_reciprocal_rank: float

    @property
    def enrichment_over_base_rate(self) -> float:
        if self.base_rate <= 0.0:
            raise ValueError("base_rate must be positive")
        return self.average_precision / self.base_rate


def _metric_dict(metric: SpeciesMetric) -> dict[str, float | str]:
    result: dict[str, float | str] = asdict(metric)
    result["enrichment_over_base_rate"] = metric.enrichment_over_base_rate
    return result


def summarize_species_metrics(
    metrics: Iterable[SpeciesMetric],
    *,
    primary_species: tuple[str, ...],
    pressure_species: tuple[str, ...],
) -> dict[str, object]:
    """Summarize plants and pressure species in non-interchangeable panels."""
    if set(primary_species) & set(pressure_species):
        raise ValueError("primary and pressure species must be disjoint")
    by_species = {metric.species: metric for metric in metrics}
    required = set(primary_species) | set(pressure_species)
    missing = required - set(by_species)
    if missing:
        raise RuntimeError(f"missing species metrics: {sorted(missing)}")
    primary_ap = [by_species[species].average_precision for species in primary_species]
    return {
        "primary_panel_title": "three_crop_primary_analysis",
        "pressure_panel_title": "magnaporthe_pressure_test",
        "primary_species": {
            species: _metric_dict(by_species[species]) for species in primary_species
        },
        "pressure_species": {
            species: _metric_dict(by_species[species]) for species in pressure_species
        },
        "primary_macro_average_precision": round(sum(primary_ap) / len(primary_ap), 12),
    }


def claim_class_for_multispecies_result(
    *,
    strict_ci_lower: float,
    species_deltas: dict[str, float],
    same_frozen_inputs: bool,
) -> str:
    """Return only the claim class allowed by strict-track evidence."""
    required = {"arabidopsis", "rice", "tomato"}
    if (
        same_frozen_inputs
        and required <= set(species_deltas)
        and strict_ci_lower > 0.0
        and all(species_deltas[species] > 0.0 for species in required)
    ):
        return "homology_robust_multiplant_within_dataset_ranking"
    return "literature_comparable_within_dataset_improvement_only"
