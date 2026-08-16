"""Target-label-free cross-crop scoring with a fail-closed source firewall."""

from __future__ import annotations

import json
import math
import random
import statistics
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from plantpersulf.evaluation.metrics import average_precision
from plantpersulf.evaluation.rank_ensemble import ApplicabilityEnvelope
from plantpersulf.features.site_biology import BIOLOGY_FEATURE_NAMES
from plantpersulf.models.additive_pu_ranker import (
    AdditivePuConfig,
    fit_additive_pu_ranker,
)
from plantpersulf.models.traditional import (
    TrainOnlyScaler,
    pu_logistic_regression_scores,
)
from plantpersulf.provenance.hashing import hash_file
from plantpersulf.workflows.tomato_ranker_v2 import (
    _lower_quantile,
    _percentiles,
    configure_run_logger,
    sensitivity_priors,
)


@dataclass(frozen=True)
class SourceBatch:
    species: str
    study_accession: str
    feature_names: tuple[str, ...]
    features: tuple[tuple[float, ...], ...]
    labels: tuple[str, ...]
    protein_ids: tuple[str, ...]
    source_sha256: str
    cluster_ids: tuple[str, ...] = ()
    batch_condition_ids: tuple[str, ...] = ()


@dataclass(frozen=True)
class TargetCandidateBatch:
    site_keys: tuple[str, ...]
    feature_names: tuple[str, ...]
    features: tuple[tuple[float, ...], ...]
    source_sha256: str


@dataclass(frozen=True)
class SourceDomainEvidence:
    source_domain: str
    candidate_ap: float
    baseline_ap: float
    preprocessing_fit_scope: str = "training_fold_only"

    def to_dict(self) -> dict[str, object]:
        return {
            "source_domain": self.source_domain,
            "candidate_ap": self.candidate_ap,
            "baseline_ap": self.baseline_ap,
            "preprocessing_fit_scope": self.preprocessing_fit_scope,
        }


@dataclass(frozen=True)
class SourceTransferDecision:
    admitted: bool
    reason: str
    domain_evidence: tuple[SourceDomainEvidence, ...] = ()


@dataclass(frozen=True)
class SourceEvidenceProfile:
    minimum_coverage_met: bool
    within_species_replication_supported: bool
    study_counts: dict[str, int]
    within_species_replication_by_species: dict[str, bool]
    pure_species_effect_supported: bool = False
    universal_cross_crop_generalization_supported: bool = False


PLANT_MAIN_SPECIES = frozenset({"Arabidopsis thaliana", "Oryza sativa"})
MAGNAPORTHE_PRESSURE_SPECIES = "Magnaporthe oryzae"


def partition_plant_main_and_pressure_sources(
    sources: tuple[SourceBatch, ...],
) -> tuple[tuple[SourceBatch, ...], tuple[SourceBatch, ...]]:
    """Keep fungal stress data out of the plant model and tomato scoring fit."""
    plant_main = tuple(
        source for source in sources if source.species in PLANT_MAIN_SPECIES
    )
    pressure = tuple(
        source for source in sources if source.species == MAGNAPORTHE_PRESSURE_SPECIES
    )
    if {source.species for source in plant_main} != PLANT_MAIN_SPECIES:
        raise RuntimeError("plant main analysis lacks a required source species")
    if len(pressure) != 1:
        raise RuntimeError(
            "Magnaporthe pressure test requires exactly one source batch"
        )
    if len(plant_main) + len(pressure) != len(sources):
        raise RuntimeError("source partition contains an unsupported non-plant domain")
    return plant_main, pressure


def cross_crop_arm_admissible(
    plant_transfer_passed: bool,
    pressure_test_passed: bool,
    target_applicability_passed: bool,
) -> bool:
    """Admit no X-arm candidates unless every frozen cross-crop check passes."""
    return (
        plant_transfer_passed
        and pressure_test_passed
        and target_applicability_passed
    )


@dataclass(frozen=True)
class TrainFoldPreprocessor:
    medians: tuple[float, ...]
    scaler: TrainOnlyScaler

    @classmethod
    def fit(cls, features: tuple[tuple[float, ...], ...]) -> TrainFoldPreprocessor:
        if not features or not features[0]:
            raise ValueError("training features must not be empty")
        width = len(features[0])
        if any(len(row) != width for row in features):
            raise ValueError("training feature width mismatch")
        medians = []
        for column in range(width):
            observed = [row[column] for row in features if math.isfinite(row[column])]
            if not observed:
                raise ValueError("training feature column has no observed value")
            medians.append(float(statistics.median(observed)))
        imputed = [
            [
                value if math.isfinite(value) else medians[index]
                for index, value in enumerate(row)
            ]
            for row in features
        ]
        return cls(tuple(medians), TrainOnlyScaler.fit(imputed))

    def transform(
        self, features: tuple[tuple[float, ...], ...]
    ) -> tuple[tuple[float, ...], ...]:
        imputed = [
            [
                value if math.isfinite(value) else self.medians[index]
                for index, value in enumerate(row)
            ]
            for row in features
        ]
        return tuple(tuple(row) for row in self.scaler.transform(imputed))


@dataclass(frozen=True)
class FoldEvaluation:
    candidate_ap: float
    baseline_ap: float
    preprocessing_fit_scope: str = "training_fold_only"


def evaluate_pu_fold(
    train_features: tuple[tuple[float, ...], ...],
    train_labels: tuple[str, ...],
    train_protein_ids: tuple[str, ...],
    test_features: tuple[tuple[float, ...], ...],
    test_labels: tuple[str, ...],
    feature_names: tuple[str, ...],
    model_seeds: tuple[int, ...] = (0, 1, 2, 3, 4),
    additive_epochs: int = 200,
) -> FoldEvaluation:
    processor = TrainFoldPreprocessor.fit(train_features)
    train_x = [list(row) for row in processor.transform(train_features)]
    test_x = [list(row) for row in processor.transform(test_features)]
    priors = sensitivity_priors(
        train_labels.count("positive"), len(train_labels), (1.0, 1.5, 2.0), 0.5
    )
    baseline_columns = [
        pu_logistic_regression_scores(train_x, list(train_labels), test_x, model_seed)
        for model_seed in model_seeds
    ]
    representative_seed = model_seeds[0]
    candidate_columns = [
        fit_additive_pu_ranker(
            train_x,
            list(train_labels),
            list(train_protein_ids),
            feature_names,
            AdditivePuConfig(
                class_prior=prior,
                seed=representative_seed,
                epochs=additive_epochs,
            ),
        ).score(test_x)
        for prior in priors
    ]
    candidate_scores = [
        statistics.median(column[index] for column in candidate_columns)
        for index in range(len(test_x))
    ]
    baseline_scores = [
        statistics.median(column[index] for column in baseline_columns)
        for index in range(len(test_x))
    ]
    candidate_ap = average_precision(
        list(zip(candidate_scores, test_labels, strict=True))
    )
    baseline_ap = average_precision(
        list(zip(baseline_scores, test_labels, strict=True))
    )
    if candidate_ap is None or baseline_ap is None:
        raise RuntimeError("held-out fold has no positive")
    return FoldEvaluation(candidate_ap, baseline_ap)


@dataclass(frozen=True)
class WithinSourceSplitPlan:
    grouping_unit: str
    repeated_fold_assignments: tuple[tuple[int, ...], ...]
    batch_condition_holdouts: tuple[tuple[str, tuple[int, ...], tuple[int, ...]], ...]


@dataclass(frozen=True)
class WithinSourceFoldEvidence:
    split_type: str
    split_id: str
    candidate_ap: float
    baseline_ap: float

    def to_dict(self) -> dict[str, object]:
        return {
            "split_type": self.split_type,
            "split_id": self.split_id,
            "candidate_ap": self.candidate_ap,
            "baseline_ap": self.baseline_ap,
        }


@dataclass(frozen=True)
class WithinSourceEvidence:
    species: str
    study_accession: str
    grouping_unit: str
    grouped_fold_count: int
    batch_condition_holdout_count: int
    median_candidate_ap: float
    median_baseline_ap: float
    preprocessing_fit_scope: str = "training_fold_only"
    fold_evaluations: tuple[WithinSourceFoldEvidence, ...] = ()

    def to_dict(self) -> dict[str, object]:
        return {
            "species": self.species,
            "study_accession": self.study_accession,
            "grouping_unit": self.grouping_unit,
            "grouped_fold_count": self.grouped_fold_count,
            "batch_condition_holdout_count": self.batch_condition_holdout_count,
            "median_candidate_ap": self.median_candidate_ap,
            "median_baseline_ap": self.median_baseline_ap,
            "preprocessing_fit_scope": self.preprocessing_fit_scope,
            "fold_evaluations": [item.to_dict() for item in self.fold_evaluations],
        }


def source_evidence_profile(
    sources: tuple[SourceBatch, ...], ideal_studies_per_species: int
) -> SourceEvidenceProfile:
    counts = {
        species: len(
            {source.study_accession for source in sources if source.species == species}
        )
        for species in sorted({source.species for source in sources})
    }
    return SourceEvidenceProfile(
        minimum_coverage_met=bool(counts)
        and all(count >= 1 for count in counts.values()),
        within_species_replication_supported=bool(counts)
        and all(count >= ideal_studies_per_species for count in counts.values()),
        study_counts=counts,
        within_species_replication_by_species={
            species: count >= ideal_studies_per_species
            for species, count in counts.items()
        },
    )


def repeated_grouped_folds(
    group_ids: tuple[str, ...],
    labels: tuple[str, ...],
    n_folds: int,
    repetitions: int,
    seed: int,
) -> tuple[tuple[int, ...], ...]:
    if len(group_ids) != len(labels):
        raise ValueError("group and label row count mismatch")
    groups = sorted(set(group_ids))
    if len(groups) < n_folds:
        raise ValueError(
            f"repeated {n_folds}-fold CV requires at least {n_folds} groups"
        )
    positive_groups = sorted(
        {
            group
            for group, label in zip(group_ids, labels, strict=True)
            if label == "positive"
        }
    )
    if len(positive_groups) < n_folds:
        raise ValueError(
            f"repeated {n_folds}-fold CV requires at least {n_folds} "
            "positive-bearing groups"
        )
    other_groups = sorted(set(groups) - set(positive_groups))
    assignments = []
    for repetition in range(repetitions):
        rng = random.Random(seed + repetition)
        positive, other = positive_groups[:], other_groups[:]
        rng.shuffle(positive)
        rng.shuffle(other)
        fold_of = {group: index % n_folds for index, group in enumerate(positive)}
        for index, group in enumerate(other):
            fold_of[group] = index % n_folds
        assignments.append(tuple(fold_of[group] for group in group_ids))
    return tuple(assignments)


def build_within_source_split_plan(
    source: SourceBatch,
    n_folds: int,
    repetitions: int,
    seed: int,
) -> WithinSourceSplitPlan:
    use_clusters = bool(source.cluster_ids) and all(source.cluster_ids)
    if use_clusters:
        if len(source.cluster_ids) != len(source.labels):
            raise ValueError("cluster ID row count mismatch")
        group_ids = source.cluster_ids
        grouping_unit = "homology_cluster"
    else:
        group_ids = source.protein_ids
        grouping_unit = "protein"
    assignments = repeated_grouped_folds(
        group_ids, source.labels, n_folds, repetitions, seed
    )
    condition_holdouts = []
    if source.batch_condition_ids:
        if len(source.batch_condition_ids) != len(source.labels):
            raise ValueError("batch/condition ID row count mismatch")
        conditions = sorted(set(source.batch_condition_ids))
        if len(conditions) > 1:
            for condition in conditions:
                test = tuple(
                    index
                    for index, value in enumerate(source.batch_condition_ids)
                    if value == condition
                )
                train = tuple(
                    index
                    for index, value in enumerate(source.batch_condition_ids)
                    if value != condition
                )
                condition_holdouts.append((condition, train, test))
    return WithinSourceSplitPlan(
        grouping_unit=grouping_unit,
        repeated_fold_assignments=assignments,
        batch_condition_holdouts=tuple(condition_holdouts),
    )


def evaluate_within_source(
    source: SourceBatch,
    n_folds: int,
    repetitions: int,
    seed: int,
    model_seeds: tuple[int, ...] = (0, 1, 2, 3, 4),
    additive_epochs: int = 200,
    on_split_complete: Callable[[WithinSourceFoldEvidence], None] | None = None,
) -> WithinSourceEvidence:
    plan = build_within_source_split_plan(source, n_folds, repetitions, seed)
    splits: list[
        tuple[str, str, tuple[int, ...], tuple[int, ...]]
    ] = []
    for repetition, assignment in enumerate(plan.repeated_fold_assignments):
        for fold in range(n_folds):
            test = tuple(
                index for index, value in enumerate(assignment) if value == fold
            )
            train = tuple(
                index for index, value in enumerate(assignment) if value != fold
            )
            splits.append(
                ("repeated_grouped_cv", f"r{repetition}_f{fold}", train, test)
            )
    splits.extend(
        ("leave_batch_condition_out", condition, train, test)
        for condition, train, test in plan.batch_condition_holdouts
    )

    evaluations = []
    fold_evaluations = []
    for split_type, split_id, train, test in splits:
        evaluation = evaluate_pu_fold(
            train_features=tuple(source.features[index] for index in train),
            train_labels=tuple(source.labels[index] for index in train),
            train_protein_ids=tuple(source.protein_ids[index] for index in train),
            test_features=tuple(source.features[index] for index in test),
            test_labels=tuple(source.labels[index] for index in test),
            feature_names=source.feature_names,
            model_seeds=model_seeds,
            additive_epochs=additive_epochs,
        )
        evaluations.append(evaluation)
        fold_evidence = WithinSourceFoldEvidence(
            split_type=split_type,
            split_id=split_id,
            candidate_ap=evaluation.candidate_ap,
            baseline_ap=evaluation.baseline_ap,
        )
        fold_evaluations.append(fold_evidence)
        if on_split_complete is not None:
            on_split_complete(fold_evidence)
    return WithinSourceEvidence(
        species=source.species,
        study_accession=source.study_accession,
        grouping_unit=plan.grouping_unit,
        grouped_fold_count=n_folds * repetitions,
        batch_condition_holdout_count=len(plan.batch_condition_holdouts),
        median_candidate_ap=statistics.median(
            evaluation.candidate_ap for evaluation in evaluations
        ),
        median_baseline_ap=statistics.median(
            evaluation.baseline_ap for evaluation in evaluations
        ),
        fold_evaluations=tuple(fold_evaluations),
    )


def leave_one_source_domain_out(
    sources: tuple[SourceBatch, ...],
) -> tuple[tuple[str, tuple[SourceBatch, ...], tuple[SourceBatch, ...]], ...]:
    holdouts = []
    for held_out in sorted(
        sources, key=lambda item: (item.species, item.study_accession)
    ):
        name = f"{held_out.species}|{held_out.study_accession}"
        training = tuple(source for source in sources if source is not held_out)
        holdouts.append((name, training, (held_out,)))
    return tuple(holdouts)


def validate_target_label_free_sources(
    sources: tuple[SourceBatch, ...], target_species: str
) -> None:
    if not sources:
        raise RuntimeError("at least one registered source batch is required")
    keys = [(source.species, source.study_accession) for source in sources]
    if len(keys) != len(set(keys)):
        raise RuntimeError("duplicate source species-study batch")
    for source in sources:
        if source.species == target_species:
            raise RuntimeError("target label leakage: target species in source batches")
        if not source.source_sha256 or not source.study_accession:
            raise RuntimeError("source batch lacks registered accession or SHA256")
        if source.feature_names != BIOLOGY_FEATURE_NAMES:
            raise RuntimeError("source feature order differs from biological contract")
        if set(source.labels) - {"positive", "unlabeled"}:
            raise RuntimeError("source batch contains a forbidden label")
        if not source.features or not (
            len(source.features) == len(source.labels) == len(source.protein_ids)
        ):
            raise RuntimeError("source batch row count mismatch")
        for metadata in (source.cluster_ids, source.batch_condition_ids):
            if metadata and len(metadata) != len(source.labels):
                raise RuntimeError("optional source metadata row count mismatch")
            if metadata and not all(metadata):
                raise RuntimeError("optional source metadata contains an empty ID")
        if any(len(row) != len(source.feature_names) for row in source.features):
            raise RuntimeError("source feature width mismatch")


def validate_source_coverage(
    sources: tuple[SourceBatch, ...],
    minimum_source_species: int,
    minimum_total_source_studies: int,
    minimum_source_studies_per_species: int,
) -> None:
    """Enforce multi-species coverage while permitting singleton species."""
    species = {source.species for source in sources}
    if len(species) < minimum_source_species:
        raise RuntimeError("source species coverage is below the frozen minimum")
    studies = {(source.species, source.study_accession) for source in sources}
    if len(studies) < minimum_total_source_studies:
        raise RuntimeError("total source study coverage is below the frozen minimum")
    for name in species:
        species_studies = {
            source.study_accession for source in sources if source.species == name
        }
        if len(species_studies) < minimum_source_studies_per_species:
            raise RuntimeError(f"source species lacks study support: {name}")


def source_transfer_gate(sources: tuple[SourceBatch, ...]) -> SourceTransferDecision:
    if len({source.species for source in sources}) < 2:
        return SourceTransferDecision(False, "fewer_than_two_source_species")
    domain_evidence = []
    failed_domains = []
    for name, training, held_out in leave_one_source_domain_out(sources):
        train_x = tuple(row for batch in training for row in batch.features)
        train_y = tuple(label for batch in training for label in batch.labels)
        train_proteins = tuple(
            f"{batch.species}|{batch.study_accession}|{protein}"
            for batch in training
            for protein in batch.protein_ids
        )
        test_x = tuple(row for batch in held_out for row in batch.features)
        test_y = tuple(label for batch in held_out for label in batch.labels)
        evaluation = evaluate_pu_fold(
            train_features=train_x,
            train_labels=train_y,
            train_protein_ids=train_proteins,
            test_features=test_x,
            test_labels=test_y,
            feature_names=training[0].feature_names,
        )
        domain_evidence.append(
            SourceDomainEvidence(
                source_domain=name,
                candidate_ap=evaluation.candidate_ap,
                baseline_ap=evaluation.baseline_ap,
            )
        )
        if evaluation.candidate_ap <= evaluation.baseline_ap:
            failed_domains.append(name)
    evidence = tuple(domain_evidence)
    if failed_domains:
        return SourceTransferDecision(
            False,
            "source_domain_transfer_not_better:" + ",".join(failed_domains),
            evidence,
        )
    return SourceTransferDecision(True, "all_source_domains_better", evidence)


def magnaporthe_pressure_test(
    plant_main_sources: tuple[SourceBatch, ...],
    pressure_sources: tuple[SourceBatch, ...],
) -> SourceTransferDecision:
    """Stress-test the plant-only fit on Magnaporthe without scoring tomato from it."""
    if not plant_main_sources or len(pressure_sources) != 1:
        raise RuntimeError("pressure test requires plant sources and one fungal source")
    if any(source.species not in PLANT_MAIN_SPECIES for source in plant_main_sources):
        raise RuntimeError("pressure-test training sources must be plant-only")
    if pressure_sources[0].species != MAGNAPORTHE_PRESSURE_SPECIES:
        raise RuntimeError("pressure-test holdout must be Magnaporthe")
    train_x = tuple(row for batch in plant_main_sources for row in batch.features)
    train_y = tuple(label for batch in plant_main_sources for label in batch.labels)
    train_proteins = tuple(
        f"{batch.study_accession}|{protein}"
        for batch in plant_main_sources
        for protein in batch.protein_ids
    )
    held_out = pressure_sources[0]
    evaluation = evaluate_pu_fold(
        train_features=train_x,
        train_labels=train_y,
        train_protein_ids=train_proteins,
        test_features=held_out.features,
        test_labels=held_out.labels,
        feature_names=held_out.feature_names,
    )
    evidence = (
        SourceDomainEvidence(
            source_domain=f"{held_out.species}|{held_out.study_accession}",
            candidate_ap=evaluation.candidate_ap,
            baseline_ap=evaluation.baseline_ap,
        ),
    )
    if evaluation.candidate_ap <= evaluation.baseline_ap:
        return SourceTransferDecision(
            False, "pressure_domain_transfer_not_better", evidence
        )
    return SourceTransferDecision(True, "all_pressure_domains_better", evidence)


def _is_frozen_v2_partition(cfg: dict[str, Any]) -> bool:
    partition = cfg.get("analysis_partition")
    if partition is None:
        return False
    if partition != {
        "plant_main_source_species": ["Arabidopsis thaliana", "Oryza sativa"],
        "plant_main_minimum_source_species": 2,
        "plant_main_minimum_total_source_studies": 3,
        "pressure_test_species": "Magnaporthe oryzae",
        "pressure_test_role": "held_out_stress_test_only",
        "target_scoring_sources": "plant_main_only",
    }:
        raise RuntimeError("plant-main and pressure-test partition differs from v2")
    if cfg.get("x_arm_policy") != {
        "required_checks": [
            "plant_source_domain_transfer",
            "magnaporthe_pressure_test",
            "target_applicability",
        ],
        "on_any_failure": "cross_crop_k_zero",
    }:
        raise RuntimeError("cross-crop candidate-arm policy differs from v2")
    return True


def run_cross_crop_target_label_free(
    config_path: Path,
    sources: tuple[SourceBatch, ...],
    target_candidates: TargetCandidateBatch,
    output_dir: Path,
) -> dict[str, object]:
    cfg: dict[str, Any] = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if cfg.get("claim_class") != "target_label_free_transfer_not_gate2":
        raise RuntimeError("cross-crop workflow cannot emit a Gate 2 claim")
    target_species = str(cfg["target_species"])
    v2_partitioned_analysis = _is_frozen_v2_partition(cfg)
    evaluation = dict(cfg["evaluation"])
    within_source_validation = dict(cfg["within_source_validation"])
    claims = dict(cfg["claims"])
    wetlab = dict(cfg["wetlab"])
    if (
        evaluation["source_selection"] != "leave_one_source_domain_out"
        or evaluation["source_domain_unit"] != "species_study"
        or evaluation["target_labels_visible_during_fit"] is not False
        or evaluation["target_labels_visible_during_tuning"] is not False
        or within_source_validation
        != {
            "split": "repeated_grouped_protein_or_homology_cluster_cv",
            "folds": 10,
            "repetitions": 5,
            "grouping_preference": "homology_cluster_then_protein",
            "preprocessing_fit_scope": "training_fold_only",
            "missing_value_policy": "training_fold_median",
            "leave_batch_condition_out_when_available": True,
        }
        or claims
        != {
            "pure_species_effect_supported": False,
            "universal_cross_crop_generalization_supported": False,
            "final_generalization_requires_tomato_blind_validation": True,
        }
        or wetlab["enabled"] is not False
        or float(wetlab["max_noncontrol_fraction"]) > 0.25
    ):
        raise RuntimeError("target-label or wet-lab policy differs from frozen config")
    validate_target_label_free_sources(sources, target_species)
    if {source.species for source in sources} != set(cfg["source_species"]):
        raise RuntimeError("source species differ from frozen config")
    validate_source_coverage(
        sources,
        int(cfg["minimum_source_species"]),
        int(cfg["minimum_total_source_studies"]),
        int(cfg["minimum_source_studies_per_species"]),
    )
    if (
        int(cfg["minimum_source_studies_per_species"]) != 1
        or int(cfg["ideal_source_studies_per_species"]) < 2
    ):
        raise RuntimeError("source evidence grading differs from frozen config")
    if v2_partitioned_analysis:
        scoring_sources, pressure_sources = partition_plant_main_and_pressure_sources(
            sources
        )
        validate_source_coverage(
            scoring_sources,
            minimum_source_species=2,
            minimum_total_source_studies=3,
            minimum_source_studies_per_species=1,
        )
    else:
        scoring_sources = sources
        pressure_sources = ()
    if output_dir.exists():
        raise FileExistsError(output_dir)
    if (
        not target_candidates.source_sha256
        or target_candidates.feature_names != BIOLOGY_FEATURE_NAMES
    ):
        raise RuntimeError(
            "target candidate batch lacks registered biological contract"
        )
    if (
        not target_candidates.features
        or len(target_candidates.site_keys) != len(target_candidates.features)
        or len(set(target_candidates.site_keys)) != len(target_candidates.site_keys)
    ):
        raise RuntimeError("target candidate row count or keys are invalid")
    if any(
        len(row) != len(BIOLOGY_FEATURE_NAMES) for row in target_candidates.features
    ):
        raise RuntimeError("target feature width mismatch")
    output_dir.mkdir(parents=True, exist_ok=False)
    logger = configure_run_logger(output_dir / "run.log")
    logger.info("phase=run_start config=%s", config_path)
    profile = source_evidence_profile(
        scoring_sources, int(cfg["ideal_source_studies_per_species"])
    )
    within_source_results = []
    for source in sorted(
        scoring_sources, key=lambda item: (item.species, item.study_accession)
    ):
        logger.info(
            "phase=within_source_start species=%s study=%s",
            source.species,
            source.study_accession,
        )

        def _log_split_complete(
            fold: WithinSourceFoldEvidence,
            source: SourceBatch = source,
        ) -> None:
            logger.info(
                "phase=within_source_split_complete species=%s study=%s "
                "split_type=%s split_id=%s candidate_ap=%.6f baseline_ap=%.6f",
                source.species,
                source.study_accession,
                fold.split_type,
                fold.split_id,
                fold.candidate_ap,
                fold.baseline_ap,
            )

        evidence = evaluate_within_source(
            source,
            n_folds=int(within_source_validation["folds"]),
            repetitions=int(within_source_validation["repetitions"]),
            seed=0,
            on_split_complete=_log_split_complete,
        )
        within_source_results.append(evidence)
        logger.info(
            "phase=within_source_complete species=%s study=%s "
            "grouped_folds=%s condition_holdouts=%s candidate_ap=%.6f "
            "baseline_ap=%.6f",
            evidence.species,
            evidence.study_accession,
            evidence.grouped_fold_count,
            evidence.batch_condition_holdout_count,
            evidence.median_candidate_ap,
            evidence.median_baseline_ap,
        )
    within_source_evidence = tuple(within_source_results)
    decision = source_transfer_gate(scoring_sources)
    logger.info(
        "phase=source_domain_gate_complete admitted=%s domains=%s reason=%s",
        decision.admitted,
        len(decision.domain_evidence),
        decision.reason,
    )
    pressure_decision = (
        magnaporthe_pressure_test(scoring_sources, pressure_sources)
        if v2_partitioned_analysis
        else None
    )
    if pressure_decision is not None:
        logger.info(
            "phase=magnaporthe_pressure_test_complete admitted=%s reason=%s",
            pressure_decision.admitted,
            pressure_decision.reason,
        )
    raw_train_x = tuple(row for batch in scoring_sources for row in batch.features)
    train_y = [label for batch in scoring_sources for label in batch.labels]
    proteins = [
        f"{batch.species}|{batch.study_accession}|{protein}"
        for batch in scoring_sources
        for protein in batch.protein_ids
    ]
    processor = TrainFoldPreprocessor.fit(raw_train_x)
    train_x = [list(row) for row in processor.transform(raw_train_x)]
    target_x = [
        list(row) for row in processor.transform(target_candidates.features)
    ]
    priors = sensitivity_priors(
        train_y.count("positive"), len(train_y), (1.0, 1.5, 2.0), 0.5
    )
    columns = [
        _percentiles(
            fit_additive_pu_ranker(
                train_x,
                train_y,
                proteins,
                BIOLOGY_FEATURE_NAMES,
                AdditivePuConfig(class_prior=prior, seed=seed),
            ).score(target_x)
        )
        for seed in range(5)
        for prior in priors
    ]
    envelope = ApplicabilityEnvelope.fit(train_x)
    in_domain = [envelope.contains(row) for row in target_x]
    fraction = sum(in_domain) / len(in_domain)
    scores = {
        key: {
            "median_percentile": statistics.median(column[i] for column in columns),
            "conservative_percentile": _lower_quantile(
                [column[i] for column in columns], 0.2
            ),
            "in_domain": in_domain[i],
        }
        for i, key in enumerate(target_candidates.site_keys)
    }
    (output_dir / "scores.json").write_text(
        json.dumps({"scores": scores}, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    applicability_passed = fraction >= float(
        dict(cfg["applicability"])["min_target_in_domain_fraction"]
    )
    cross_crop_eligible = (
        cross_crop_arm_admissible(
            decision.admitted,
            pressure_decision.admitted,
            applicability_passed,
        )
        if pressure_decision is not None
        else decision.admitted and applicability_passed
    )
    result: dict[str, object] = {
        "complete": True,
        "claim_class": cfg["claim_class"],
        "target_species": target_species,
        "sources": [
            {
                "species": s.species,
                "study_accession": s.study_accession,
                "sha256": s.source_sha256,
            }
            for s in sources
        ],
        "target_candidate_sha256": target_candidates.source_sha256,
        "target_label_hash": None,
        "source_study_counts": profile.study_counts,
        "minimum_source_coverage_met": profile.minimum_coverage_met,
        "within_species_replication_supported": (
            profile.within_species_replication_supported
        ),
        "within_species_replication_by_species": (
            profile.within_species_replication_by_species
        ),
        "within_source_evidence": [
            evidence.to_dict() for evidence in within_source_evidence
        ],
        "source_design_limitation": (
            "A species represented by one study confounds species with study. "
            "Leave-one-source-domain-out therefore supports cross-source domain "
            "transfer, not a pure species effect."
        ),
        "source_gate_admitted": decision.admitted,
        "source_gate_reason": decision.reason,
        "source_domain_evidence": [
            evidence.to_dict() for evidence in decision.domain_evidence
        ],
        "cross_source_domain_transfer_supported": decision.admitted,
        "pure_species_effect_supported": profile.pure_species_effect_supported,
        "universal_cross_crop_generalization_supported": (
            profile.universal_cross_crop_generalization_supported
        ),
        "final_generalization_requires_tomato_blind_validation": claims[
            "final_generalization_requires_tomato_blind_validation"
        ],
        "target_in_domain_fraction": fraction,
        "applicability_passed": fraction
        >= float(dict(cfg["applicability"])["min_target_in_domain_fraction"]),
        "wetlab_eligible": False,
        "scores_sha256": hash_file(output_dir / "scores.json", "sha256"),
    }
    if pressure_decision is not None:
        result["analysis_partition"] = {
            "plant_main_source_species": sorted(PLANT_MAIN_SPECIES),
            "plant_main_studies": [
                source.study_accession for source in scoring_sources
            ],
            "pressure_test_species": MAGNAPORTHE_PRESSURE_SPECIES,
            "pressure_test_studies": [
                source.study_accession for source in pressure_sources
            ],
            "target_scoring_sources": "plant_main_only",
        }
        result["magnaporthe_pressure_test"] = {
            "admitted": pressure_decision.admitted,
            "reason": pressure_decision.reason,
            "domain_evidence": [
                evidence.to_dict() for evidence in pressure_decision.domain_evidence
            ],
        }
        result["cross_crop_arm_admissible"] = cross_crop_eligible
        result["cross_crop_k_policy"] = "result_dependent_else_zero"
    result["wetlab_eligible"] = cross_crop_eligible
    (output_dir / "manifest.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    logger.info("phase=run_complete output_dir=%s", output_dir)
    return result
