"""Species-stratified reporting and claim gates for multispecies v2."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import asdict, dataclass

from plantpersulf.evaluation.bootstrap import BootstrapResult, ClusteredScored
from plantpersulf.evaluation.metrics import (
    average_precision,
    mean_reciprocal_rank,
    recall_at_k,
)

RECALL_AT_50 = 50
RECALL_AT_200 = 200

# Task 9.6 frozen strict-track statistical policy: the same seed already used
# to freeze the strict 20% split, so the whole strict-track story reproduces
# from one documented constant. Never overridden per call site.
STRICT_BOOTSTRAP_N_BOOT = 10_000
STRICT_BOOTSTRAP_SEED = 20260811


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


def compute_species_metric(
    scored: list[tuple[float, str]], species: str
) -> SpeciesMetric:
    """Derive one species' frozen-test metric bundle from raw scored pairs.

    Reuses the PU-appropriate metric functions everywhere else in the
    project rather than reimplementing AP/Recall@K/MRR ad hoc.
    """
    if not scored:
        raise ValueError("scored must not be empty")
    ap = average_precision(scored)
    if ap is None:
        raise RuntimeError(f"no positive sites scored for species: {species}")
    positive_count = sum(1 for _, label in scored if label == "positive")
    base_rate = positive_count / len(scored)
    recall_50 = recall_at_k(scored, RECALL_AT_50)
    recall_200 = recall_at_k(scored, RECALL_AT_200)
    mrr = mean_reciprocal_rank(scored)
    assert recall_50 is not None and recall_200 is not None and mrr is not None
    return SpeciesMetric(species, ap, base_rate, recall_50, recall_200, mrr)


def pooled_average_precision(
    scored_by_species: dict[str, list[tuple[float, str]]],
    species: tuple[str, ...],
) -> float:
    """Average precision over rows *pooled* across ``species`` into one
    ranking — a different statistic from macro-averaging independent APs."""
    if not species:
        raise ValueError("species must not be empty")
    pooled: list[tuple[float, str]] = []
    for name in species:
        pooled.extend(scored_by_species[name])
    value = average_precision(pooled)
    if value is None:
        raise RuntimeError("pooled scored rows contain no positives")
    return value


def literature_track_leads(
    candidate_macro_aps: list[float], baseline_macro_aps: list[float]
) -> bool:
    """Did the candidate beat the best baseline's mean macro AP, seed for
    seed, over the literature-comparable random-protein track?"""
    if not candidate_macro_aps or len(candidate_macro_aps) != len(
        baseline_macro_aps
    ):
        raise ValueError(
            "candidate and baseline macro AP lists must be equal length and "
            "non-empty"
        )
    candidate_mean = sum(candidate_macro_aps) / len(candidate_macro_aps)
    baseline_mean = sum(baseline_macro_aps) / len(baseline_macro_aps)
    return candidate_mean > baseline_mean


def strict_track_bootstrap_delta(
    model_scored: ClusteredScored, baseline_scored: ClusteredScored
) -> BootstrapResult:
    """Frozen Task 9.6 statistical policy: 10,000-replicate protein-cluster
    bootstrap of the model-vs-baseline AP delta, seed 20260811."""
    from plantpersulf.evaluation.effect_size import (
        paired_cluster_bootstrap_delta_ci,
    )

    return paired_cluster_bootstrap_delta_ci(
        model_scored,
        baseline_scored,
        n_boot=STRICT_BOOTSTRAP_N_BOOT,
        alpha=0.05,
        seed=STRICT_BOOTSTRAP_SEED,
    )


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


HOMOLOGY_ROBUST_CLAIM = "homology_robust_multiplant_within_dataset_ranking"
LITERATURE_COMPARABLE_CLAIM = "literature_comparable_within_dataset_improvement_only"
NO_SUPPORTED_CLAIM = "no_supported_multispecies_claim"

_CLAIM_STATEMENTS = {
    LITERATURE_COMPARABLE_CLAIM: (
        "文献同口径的数据集内性能提升。 "
        "Literature-comparable within-dataset performance improvement."
    ),
    HOMOLOGY_ROBUST_CLAIM: (
        "对未见同源家族稳健的多植物数据集内排序。 "
        "Multiplant within-dataset ranking robust to unseen homology "
        "families."
    ),
    NO_SUPPORTED_CLAIM: (
        "两条轨道均未证明相对基线的优势，不得声称任何跨物种性能提升。 "
        "Neither track demonstrates an advantage over baselines; no "
        "multispecies performance-improvement claim may be made."
    ),
}

# Task 9.6: multispecies v2 is internal/within-dataset only until a tomato
# blind cohort returns real prospective data. These phrases are therefore
# unconditionally forbidden at this project stage — there is no flag to
# unlock them; that decision belongs to a future task once a blind cohort
# actually exists.
FORBIDDEN_EXTERNAL_CLAIM_PHRASES = (
    "external generalization",
    "external generalisation",
    "cross-crop generalization",
    "cross-crop generalisation",
    "prospective validation",
    "外部泛化",
    "跨作物泛化",
    "前瞻验证",
)


def claim_class_for_multispecies_result(
    *,
    strict_ci_lower: float,
    species_deltas: dict[str, float],
    same_frozen_inputs: bool,
    literature_leads: bool,
) -> str:
    """Return only the claim class allowed by strict- and literature-track
    evidence.

    Fails closed to ``NO_SUPPORTED_CLAIM`` when neither track's evidence
    clears its bar — a result with zero demonstrated advantage must never be
    presented as a literature-comparable improvement.
    """
    required = {"arabidopsis", "rice", "tomato"}
    if (
        same_frozen_inputs
        and required <= set(species_deltas)
        and strict_ci_lower > 0.0
        and all(species_deltas[species] > 0.0 for species in required)
    ):
        return HOMOLOGY_ROBUST_CLAIM
    if literature_leads:
        return LITERATURE_COMPARABLE_CLAIM
    return NO_SUPPORTED_CLAIM


def claim_statement_for_class(claim_class: str) -> str:
    """Map a claim class to its exact mandated bilingual wording."""
    try:
        return _CLAIM_STATEMENTS[claim_class]
    except KeyError as exc:
        raise ValueError(f"unknown multispecies claim class: {claim_class}") from exc


def find_forbidden_external_claims(text: str) -> list[str]:
    """Return the forbidden external/prospective phrases present in ``text``."""
    lowered = text.lower()
    return [
        phrase
        for phrase in FORBIDDEN_EXTERNAL_CLAIM_PHRASES
        if phrase.lower() in lowered
    ]


def build_multispecies_statistical_report(
    *,
    scored_by_species: dict[str, list[tuple[float, str]]],
    primary_species: tuple[str, ...],
    pressure_species: tuple[str, ...],
    strict_ci_lower: float,
    species_deltas: dict[str, float],
    same_frozen_inputs: bool,
    candidate_literature_macro_aps: list[float],
    baseline_literature_macro_aps: list[float],
) -> dict[str, object]:
    """Compose the full Task 9.6 frozen-test statistical report and claim.

    This is the single orchestrator a future frozen-test run calls: per-
    species metrics, the pooled primary-species AP, the literature-track lead
    check, the resulting claim class, and its mandated wording. It never
    unlocks or scores a frozen test itself — every scored row is supplied by
    the caller, which is responsible for having obtained it through the
    already-existing `assert_test_unlocked` gate.
    """
    metrics = [
        compute_species_metric(scored_by_species[species], species)
        for species in (*primary_species, *pressure_species)
    ]
    report = summarize_species_metrics(
        metrics, primary_species=primary_species, pressure_species=pressure_species
    )
    report["primary_pooled_average_precision"] = pooled_average_precision(
        scored_by_species, primary_species
    )
    leads = literature_track_leads(
        candidate_literature_macro_aps, baseline_literature_macro_aps
    )
    claim_class = claim_class_for_multispecies_result(
        strict_ci_lower=strict_ci_lower,
        species_deltas=species_deltas,
        same_frozen_inputs=same_frozen_inputs,
        literature_leads=leads,
    )
    report["strict_ci_lower"] = strict_ci_lower
    report["species_deltas"] = dict(species_deltas)
    report["literature_track_leads"] = leads
    report["claim_class"] = claim_class
    report["claim_statement"] = claim_statement_for_class(claim_class)
    report["gate2_eligible"] = False
    return report
