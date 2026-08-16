"""Graded cross-source evidence and leakage-safe split behavior."""

from math import nan
from pathlib import Path

import pytest
import yaml

import plantpersulf.workflows.cross_crop_target_label_free as workflow
from plantpersulf.features.site_biology import BIOLOGY_FEATURE_NAMES
from plantpersulf.workflows.cross_crop_target_label_free import (
    FoldEvaluation,
    SourceBatch,
    SourceDomainEvidence,
    SourceTransferDecision,
    TargetCandidateBatch,
    TrainFoldPreprocessor,
    WithinSourceEvidence,
    build_within_source_split_plan,
    cross_crop_arm_admissible,
    evaluate_pu_fold,
    evaluate_within_source,
    leave_one_source_domain_out,
    magnaporthe_pressure_test,
    partition_plant_main_and_pressure_sources,
    repeated_grouped_folds,
    run_cross_crop_target_label_free,
    source_evidence_profile,
    source_transfer_gate,
    validate_target_label_free_sources,
)
from scripts.run_cross_crop_target_label_free_v1 import source_batch_from_json


def _source(species: str, study: str) -> SourceBatch:
    return SourceBatch(
        species=species,
        study_accession=study,
        feature_names=BIOLOGY_FEATURE_NAMES,
        features=((0.0,) * 6, (1.0,) * 6),
        labels=("positive", "unlabeled"),
        protein_ids=(f"{study}-protein", f"{study}-protein"),
        source_sha256=f"marker-{species}-{study}",
    )


def test_single_study_species_is_admitted_but_not_called_replicated() -> None:
    sources = (
        _source("Arabidopsis thaliana", "at-1"),
        _source("Arabidopsis thaliana", "at-2"),
        _source("Oryza sativa", "os-1"),
        _source("Magnaporthe oryzae", "mo-1"),
    )
    profile = source_evidence_profile(sources, ideal_studies_per_species=2)
    assert profile.minimum_coverage_met is True
    assert profile.within_species_replication_supported is False
    assert profile.study_counts == {
        "Arabidopsis thaliana": 2,
        "Magnaporthe oryzae": 1,
        "Oryza sativa": 1,
    }
    assert profile.within_species_replication_by_species == {
        "Arabidopsis thaliana": True,
        "Magnaporthe oryzae": False,
        "Oryza sativa": False,
    }
    assert profile.pure_species_effect_supported is False
    assert profile.universal_cross_crop_generalization_supported is False


def test_plant_main_and_magnaporthe_pressure_sources_are_disjoint() -> None:
    sources = (
        _source("Arabidopsis thaliana", "at-1"),
        _source("Arabidopsis thaliana", "at-2"),
        _source("Oryza sativa", "os-1"),
        _source("Magnaporthe oryzae", "mo-1"),
    )

    plant_main, pressure = partition_plant_main_and_pressure_sources(sources)

    assert tuple(source.study_accession for source in plant_main) == (
        "at-1",
        "at-2",
        "os-1",
    )
    assert tuple(source.study_accession for source in pressure) == ("mo-1",)
    assert not {id(source) for source in plant_main} & {
        id(source) for source in pressure
    }


def test_cross_crop_arm_requires_plant_transfer_pressure_and_applicability() -> None:
    assert cross_crop_arm_admissible(True, True, True) is True
    assert cross_crop_arm_admissible(False, True, True) is False
    assert cross_crop_arm_admissible(True, False, True) is False
    assert cross_crop_arm_admissible(True, True, False) is False


def test_magnaporthe_pressure_test_trains_only_on_plant_sources(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    plant_main = (
        _source("Arabidopsis thaliana", "at-1"),
        _source("Oryza sativa", "os-1"),
    )
    pressure = (_source("Magnaporthe oryzae", "mo-1"),)
    calls: list[object] = []

    def _evaluate(**kwargs: object) -> FoldEvaluation:
        calls.append(kwargs)
        return FoldEvaluation(candidate_ap=0.8, baseline_ap=0.5)

    monkeypatch.setattr(workflow, "evaluate_pu_fold", _evaluate)
    decision = magnaporthe_pressure_test(plant_main, pressure)

    assert decision.admitted is True
    assert decision.reason == "all_pressure_domains_better"
    assert tuple(item.source_domain for item in decision.domain_evidence) == (
        "Magnaporthe oryzae|mo-1",
    )
    assert len(calls) == 1
    assert calls[0]["train_protein_ids"] == (
        "at-1|at-1-protein",
        "at-1|at-1-protein",
        "os-1|os-1-protein",
        "os-1|os-1-protein",
    )


def test_repeated_ten_fold_split_never_separates_a_protein_group() -> None:
    groups = tuple(f"protein-{index // 2}" for index in range(40))
    labels = tuple("positive" if index % 4 == 0 else "unlabeled" for index in range(40))
    assignments = repeated_grouped_folds(groups, labels, 10, 3, seed=17)
    assert len(assignments) == 3
    for assignment in assignments:
        assert set(assignment) == set(range(10))
        for group in set(groups):
            assert (
                len({assignment[i] for i, value in enumerate(groups) if value == group})
                == 1
            )


def test_preprocessing_statistics_are_fitted_on_train_rows_only() -> None:
    processor = TrainFoldPreprocessor.fit(((0.0, nan), (2.0, 4.0)))
    assert processor.transform(((100.0, nan),)) == ((99.0, 0.0),)


@pytest.mark.parametrize("field", ["cluster_ids", "batch_condition_ids"])
def test_optional_grouping_metadata_must_match_source_rows(field: str) -> None:
    kwargs = {field: ("only-one-row",)}
    source = SourceBatch(
        species="Oryza sativa",
        study_accession="os-1",
        feature_names=BIOLOGY_FEATURE_NAMES,
        features=((0.0,) * 6, (1.0,) * 6),
        labels=("positive", "unlabeled"),
        protein_ids=("protein-1", "protein-2"),
        source_sha256="marker-os-1",
        **kwargs,
    )
    with pytest.raises(RuntimeError, match="row count mismatch"):
        validate_target_label_free_sources((source,), "Solanum lycopersicum")


def test_request_parser_preserves_optional_grouping_metadata_as_tuples() -> None:
    parsed = source_batch_from_json(
        {
            "species": "Oryza sativa",
            "study_accession": "os-1",
            "feature_names": list(BIOLOGY_FEATURE_NAMES),
            "features": [[0.0] * 6, [1.0] * 6],
            "labels": ["positive", "unlabeled"],
            "protein_ids": ["protein-1", "protein-2"],
            "cluster_ids": ["cluster-1", "cluster-2"],
            "batch_condition_ids": ["batch-a", "batch-b"],
            "source_sha256": "marker-os-1",
        }
    )
    assert parsed.cluster_ids == ("cluster-1", "cluster-2")
    assert parsed.batch_condition_ids == ("batch-a", "batch-b")


def test_leave_one_source_domain_out_holds_entire_species_study_batch() -> None:
    sources = (
        _source("Arabidopsis thaliana", "at-1"),
        _source("Arabidopsis thaliana", "at-2"),
        _source("Oryza sativa", "os-1"),
    )
    holdouts = leave_one_source_domain_out(sources)
    assert tuple(name for name, _, _ in holdouts) == (
        "Arabidopsis thaliana|at-1",
        "Arabidopsis thaliana|at-2",
        "Oryza sativa|os-1",
    )
    for name, training, held_out in holdouts:
        assert len(held_out) == 1
        assert f"{held_out[0].species}|{held_out[0].study_accession}" == name
        assert held_out[0] not in training


def test_source_transfer_gate_evaluates_only_species_study_domains(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    sources = (
        _source("Arabidopsis thaliana", "at-1"),
        _source("Arabidopsis thaliana", "at-2"),
        _source("Oryza sativa", "os-1"),
        _source("Magnaporthe oryzae", "mo-1"),
    )
    calls: list[object] = []

    def _evaluate(**kwargs: object) -> FoldEvaluation:
        calls.append(kwargs)
        return FoldEvaluation(candidate_ap=0.8, baseline_ap=0.5)

    monkeypatch.setattr(workflow, "evaluate_pu_fold", _evaluate)
    decision = source_transfer_gate(sources)
    assert decision.admitted is True
    assert decision.reason == "all_source_domains_better"
    assert len(calls) == len(sources)
    assert tuple(item.source_domain for item in decision.domain_evidence) == (
        "Arabidopsis thaliana|at-1",
        "Arabidopsis thaliana|at-2",
        "Magnaporthe oryzae|mo-1",
        "Oryza sativa|os-1",
    )
    assert all(
        item.preprocessing_fit_scope == "training_fold_only"
        for item in decision.domain_evidence
    )


def test_frozen_task9_config_uses_graded_evidence_and_source_domains() -> None:
    cfg = yaml.safe_load(
        Path("configs/experiments/cross_crop_target_label_free_v1.yaml").read_text(
            encoding="utf-8"
        )
    )
    assert cfg["minimum_source_studies_per_species"] == 1
    assert cfg["ideal_source_studies_per_species"] == 2
    assert cfg["within_source_validation"] == {
        "split": "repeated_grouped_protein_or_homology_cluster_cv",
        "folds": 10,
        "repetitions": 5,
        "grouping_preference": "homology_cluster_then_protein",
        "preprocessing_fit_scope": "training_fold_only",
        "missing_value_policy": "training_fold_median",
        "leave_batch_condition_out_when_available": True,
    }
    assert cfg["evaluation"]["source_selection"] == "leave_one_source_domain_out"
    assert cfg["evaluation"]["source_domain_unit"] == "species_study"
    assert cfg["claims"] == {
        "pure_species_effect_supported": False,
        "universal_cross_crop_generalization_supported": False,
        "final_generalization_requires_tomato_blind_validation": True,
    }


def test_task9_manifest_reports_graded_evidence_without_species_claims(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = (
        _source("Arabidopsis thaliana", "at-1"),
        _source("Arabidopsis thaliana", "at-2"),
        _source("Oryza sativa", "os-1"),
        _source("Magnaporthe oryzae", "mo-1"),
    )

    def _within(source: SourceBatch, **_: object) -> WithinSourceEvidence:
        return WithinSourceEvidence(
            species=source.species,
            study_accession=source.study_accession,
            grouping_unit="protein",
            grouped_fold_count=50,
            batch_condition_holdout_count=0,
            median_candidate_ap=0.8,
            median_baseline_ap=0.5,
        )

    class _Model:
        def score(self, features: list[list[float]]) -> tuple[float, ...]:
            return tuple(row[0] for row in features)

    monkeypatch.setattr(workflow, "evaluate_within_source", _within)
    monkeypatch.setattr(
        workflow,
        "source_transfer_gate",
        lambda _: SourceTransferDecision(
            True,
            "all_source_domains_better",
            (
                SourceDomainEvidence(
                    source_domain="Arabidopsis thaliana|at-1",
                    candidate_ap=0.8,
                    baseline_ap=0.5,
                ),
            ),
        ),
    )
    monkeypatch.setattr(workflow, "fit_additive_pu_ranker", lambda *_, **__: _Model())
    target = TargetCandidateBatch(
        site_keys=("tomato-protein|C1", "tomato-protein|C2"),
        feature_names=BIOLOGY_FEATURE_NAMES,
        features=((0.25,) * 6, (0.75,) * 6),
        source_sha256="marker-tomato-candidates",
    )
    result = run_cross_crop_target_label_free(
        Path("configs/experiments/cross_crop_target_label_free_v1.yaml"),
        sources,
        target,
        tmp_path / "task9",
    )
    assert result["source_study_counts"] == {
        "Arabidopsis thaliana": 2,
        "Magnaporthe oryzae": 1,
        "Oryza sativa": 1,
    }
    assert result["within_species_replication_supported"] is False
    assert result["within_species_replication_by_species"] == {
        "Arabidopsis thaliana": True,
        "Magnaporthe oryzae": False,
        "Oryza sativa": False,
    }
    assert len(result["within_source_evidence"]) == 4
    assert result["source_domain_evidence"] == [
        {
            "source_domain": "Arabidopsis thaliana|at-1",
            "candidate_ap": 0.8,
            "baseline_ap": 0.5,
            "preprocessing_fit_scope": "training_fold_only",
        }
    ]
    assert result["cross_source_domain_transfer_supported"] is True
    assert result["pure_species_effect_supported"] is False
    assert result["universal_cross_crop_generalization_supported"] is False
    assert result["final_generalization_requires_tomato_blind_validation"] is True
    log = (tmp_path / "task9" / "run.log").read_text(encoding="utf-8")
    assert "phase=within_source_complete" in log
    assert "phase=source_domain_gate_complete" in log
    assert "phase=run_complete" in log


def test_v2_scores_tomato_from_plant_main_and_reports_fungal_pressure_only(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = (
        _source("Arabidopsis thaliana", "at-1"),
        _source("Arabidopsis thaliana", "at-2"),
        _source("Oryza sativa", "os-1"),
        _source("Magnaporthe oryzae", "mo-1"),
    )
    within_calls: list[str] = []
    transfer_calls: list[tuple[str, ...]] = []
    pressure_calls: list[tuple[str, ...]] = []

    def _within(source: SourceBatch, **_: object) -> WithinSourceEvidence:
        within_calls.append(source.study_accession)
        return WithinSourceEvidence(
            species=source.species,
            study_accession=source.study_accession,
            grouping_unit="protein",
            grouped_fold_count=50,
            batch_condition_holdout_count=0,
            median_candidate_ap=0.8,
            median_baseline_ap=0.5,
        )

    def _transfer(items: tuple[SourceBatch, ...]) -> SourceTransferDecision:
        transfer_calls.append(tuple(item.study_accession for item in items))
        return SourceTransferDecision(True, "all_source_domains_better")

    def _pressure(
        plant_main: tuple[SourceBatch, ...], pressure: tuple[SourceBatch, ...]
    ) -> SourceTransferDecision:
        pressure_calls.append(
            tuple(item.study_accession for item in (*plant_main, *pressure))
        )
        return SourceTransferDecision(True, "all_pressure_domains_better")

    class _Model:
        def score(self, features: list[list[float]]) -> tuple[float, ...]:
            return tuple(row[0] for row in features)

    monkeypatch.setattr(workflow, "evaluate_within_source", _within)
    monkeypatch.setattr(workflow, "source_transfer_gate", _transfer)
    monkeypatch.setattr(workflow, "magnaporthe_pressure_test", _pressure)
    monkeypatch.setattr(workflow, "fit_additive_pu_ranker", lambda *_, **__: _Model())
    target = TargetCandidateBatch(
        site_keys=("tomato-protein|C1", "tomato-protein|C2"),
        feature_names=BIOLOGY_FEATURE_NAMES,
        features=((0.25,) * 6, (0.75,) * 6),
        source_sha256="marker-tomato-candidates",
    )

    result = run_cross_crop_target_label_free(
        Path("configs/experiments/cross_crop_target_label_free_v2.yaml"),
        sources,
        target,
        tmp_path / "task9-v2",
    )

    assert within_calls == ["at-1", "at-2", "os-1"]
    assert transfer_calls == [("at-1", "at-2", "os-1")]
    assert pressure_calls == [("at-1", "at-2", "os-1", "mo-1")]
    assert result["cross_crop_arm_admissible"] is True
    assert result["cross_crop_k_policy"] == "result_dependent_else_zero"
    assert result["magnaporthe_pressure_test"]["admitted"] is True


def test_ten_fold_split_rejects_fewer_than_ten_groups() -> None:
    with pytest.raises(ValueError, match="10 groups"):
        repeated_grouped_folds(
            tuple(f"protein-{i}" for i in range(9)),
            ("positive",) * 9,
            10,
            1,
            seed=0,
        )


def test_ten_fold_split_rejects_fewer_than_ten_positive_groups() -> None:
    with pytest.raises(ValueError, match="positive-bearing groups"):
        repeated_grouped_folds(
            tuple(f"protein-{i}" for i in range(20)),
            tuple("positive" if i < 9 else "unlabeled" for i in range(20)),
            10,
            1,
            seed=0,
        )


def test_within_source_plan_prefers_clusters_and_holds_out_conditions() -> None:
    rows = 40
    source = SourceBatch(
        species="Oryza sativa",
        study_accession="study-os-1",
        feature_names=BIOLOGY_FEATURE_NAMES,
        features=tuple((float(index),) * 6 for index in range(rows)),
        labels=tuple(
            "positive" if index % 4 == 0 else "unlabeled" for index in range(rows)
        ),
        protein_ids=tuple(f"protein-{index}" for index in range(rows)),
        source_sha256="marker-rice",
        cluster_ids=tuple(f"cluster-{index // 2}" for index in range(rows)),
        batch_condition_ids=tuple(
            "condition-a" if index < 20 else "condition-b" for index in range(rows)
        ),
    )
    plan = build_within_source_split_plan(source, 10, 3, seed=11)
    assert plan.grouping_unit == "homology_cluster"
    assert len(plan.repeated_fold_assignments) == 3
    for assignment in plan.repeated_fold_assignments:
        for cluster in set(source.cluster_ids):
            assert (
                len(
                    {
                        assignment[index]
                        for index, value in enumerate(source.cluster_ids)
                        if value == cluster
                    }
                )
                == 1
            )
    assert tuple(name for name, _, _ in plan.batch_condition_holdouts) == (
        "condition-a",
        "condition-b",
    )


def test_within_source_plan_falls_back_to_protein_groups() -> None:
    source = SourceBatch(
        species="Oryza sativa",
        study_accession="study-os-1",
        feature_names=BIOLOGY_FEATURE_NAMES,
        features=tuple((float(index),) * 6 for index in range(20)),
        labels=tuple(
            "positive" if index % 2 == 0 else "unlabeled" for index in range(20)
        ),
        protein_ids=tuple(f"protein-{index // 2}" for index in range(20)),
        source_sha256="marker-rice",
    )
    plan = build_within_source_split_plan(source, 10, 1, seed=0)
    assert plan.grouping_unit == "protein"


def test_fold_evaluation_reports_train_only_preprocessing() -> None:
    train_features = tuple(
        (float(index), float(index % 3), nan if index == 0 else 1.0)
        for index in range(20)
    )
    evaluation = evaluate_pu_fold(
        train_features=train_features,
        train_labels=("positive",) * 10 + ("unlabeled",) * 10,
        train_protein_ids=tuple(f"protein-{index}" for index in range(20)),
        test_features=(
            (100.0, 1.0, nan),
            (101.0, 2.0, 1.0),
            (102.0, 0.0, 1.0),
            (103.0, 1.0, 1.0),
        ),
        test_labels=("positive", "unlabeled", "positive", "unlabeled"),
        feature_names=("f1", "f2", "f3"),
        model_seeds=(0,),
        additive_epochs=2,
    )
    assert evaluation.preprocessing_fit_scope == "training_fold_only"
    assert 0.0 <= evaluation.candidate_ap <= 1.0
    assert 0.0 <= evaluation.baseline_ap <= 1.0


def test_within_source_evaluation_executes_grouped_and_condition_holdouts() -> None:
    rows = 16
    source = SourceBatch(
        species="Oryza sativa",
        study_accession="study-os-1",
        feature_names=BIOLOGY_FEATURE_NAMES,
        features=tuple(
            (float(index), float(index % 3), 1.0, 0.0, 0.5, 2.0)
            for index in range(rows)
        ),
        labels=tuple(
            "positive" if (index // 2) % 2 == 0 else "unlabeled"
            for index in range(rows)
        ),
        protein_ids=tuple(f"protein-{index // 2}" for index in range(rows)),
        source_sha256="marker-rice",
        batch_condition_ids=tuple(
            "condition-a" if index < 8 else "condition-b" for index in range(rows)
        ),
    )
    completed: list[object] = []
    evidence = evaluate_within_source(
        source,
        n_folds=2,
        repetitions=1,
        seed=0,
        model_seeds=(0,),
        additive_epochs=2,
        on_split_complete=completed.append,
    )
    assert evidence.grouping_unit == "protein"
    assert evidence.grouped_fold_count == 2
    assert evidence.batch_condition_holdout_count == 2
    assert evidence.preprocessing_fit_scope == "training_fold_only"
    assert len(evidence.fold_evaluations) == 4
    assert len(completed) == 4
    assert {item.split_type for item in evidence.fold_evaluations} == {
        "repeated_grouped_cv",
        "leave_batch_condition_out",
    }
