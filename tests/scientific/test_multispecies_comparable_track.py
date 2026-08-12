"""Software-policy tests for the Task 9.4 literature-comparable track."""

from __future__ import annotations

from pathlib import Path

from plantpersulf.benchmark.literature_random_track import (
    assert_complete_model_scores,
    bind_model_to_comparison_panel,
    build_literature_random_track,
    competitor_role,
    prepare_sul_bertgru_adapter_rows,
    run_direct_comparison_roster,
    run_structure_direct_baseline,
    run_tabular_direct_baseline,
)
from plantpersulf.evaluation.comparable_track import (
    build_registered_structure_features,
)
from plantpersulf.proteomics.multispecies_v2_dataset import MultispeciesV2SiteRow


def _rows(n_proteins: int = 100) -> tuple[MultispeciesV2SiteRow, ...]:
    rows: list[MultispeciesV2SiteRow] = []
    for index in range(n_proteins):
        accession = f"POLICY_{index:03d}"
        protein_id = f"arabidopsis|{accession}"
        for position, label in ((3, "positive"), (7, "unlabeled")):
            rows.append(
                MultispeciesV2SiteRow(
                    species="arabidopsis",
                    protein_accession=accession,
                    cys_position=position,
                    label=label,
                    study_accessions=("POLICY_STUDY",) if label == "positive" else (),
                    global_protein_id=protein_id,
                    cluster_id=f"POLICY_CLUSTER_{index:03d}",
                    split="development",
                    development_fold=index % 5,
                )
            )
    return tuple(rows)


def test_random_track_has_ten_exact_protein_grouped_repetitions() -> None:
    runs = build_literature_random_track(
        _rows(),
        seeds=tuple(range(10)),
        test_fraction=0.2,
        validation_fraction_of_remaining=0.2,
        unlabeled_per_positive=1,
        sampling_seed=20260811,
    )

    assert tuple(run.seed for run in runs) == tuple(range(10))
    for run in runs:
        assigned = dict(run.protein_partitions)
        assert {
            name: sum(value == name for value in assigned.values())
            for name in ("train", "validation", "test")
        } == {"train": 64, "validation": 16, "test": 20}
        proteins = {
            name: {row.global_protein_id for row in partition}
            for name, partition in run.partitions.items()
        }
        assert {name: len(ids) for name, ids in proteins.items()} == {
            "train": 64,
            "validation": 16,
            "test": 20,
        }
        assert proteins["train"].isdisjoint(proteins["validation"])
        assert proteins["train"].isdisjoint(proteins["test"])
        assert proteins["validation"].isdisjoint(proteins["test"])


def test_random_track_samples_unlabeled_only_after_each_partition() -> None:
    runs = build_literature_random_track(
        _rows(),
        seeds=(0,),
        test_fraction=0.2,
        validation_fraction_of_remaining=0.2,
        unlabeled_per_positive=1,
        sampling_seed=20260811,
    )
    run = runs[0]

    site_owners: dict[tuple[str, int], str] = {}
    for partition_name, partition in run.partitions.items():
        assert sum(row.label == "positive" for row in partition) == sum(
            row.label == "unlabeled" for row in partition
        )
        for row in partition:
            key = (row.global_protein_id, row.cys_position)
            assert key not in site_owners
            site_owners[key] = partition_name


def test_random_track_rejects_frozen_test_rows() -> None:
    rows = list(_rows(10))
    row = rows[0]
    rows[0] = MultispeciesV2SiteRow(
        species=row.species,
        protein_accession=row.protein_accession,
        cys_position=row.cys_position,
        label=row.label,
        study_accessions=row.study_accessions,
        global_protein_id=row.global_protein_id,
        cluster_id=row.cluster_id,
        split="test",
        development_fold=None,
    )

    try:
        build_literature_random_track(
            rows,
            seeds=(0,),
            test_fraction=0.2,
            validation_fraction_of_remaining=0.2,
            unlabeled_per_positive=1,
            sampling_seed=20260811,
        )
    except RuntimeError as exc:
        assert "strict frozen test row" in str(exc)
    else:
        raise AssertionError("strict frozen test row entered comparison track")


def test_all_models_bind_to_identical_sites_labels_and_panel_hash() -> None:
    run = build_literature_random_track(
        _rows(20),
        seeds=(0,),
        test_fraction=0.2,
        validation_fraction_of_remaining=0.2,
        unlabeled_per_positive=1,
        sampling_seed=20260811,
    )[0]
    models = (
        "pu_logistic",
        "random_forest",
        "xgboost",
        "esm_linear_head",
        "structure_ranker",
        "sul_bertgru",
    )

    inputs = [bind_model_to_comparison_panel(run, model) for model in models]

    assert {item.panel_sha256 for item in inputs} == {run.panel_sha256}
    assert len({item.partition_rows for item in inputs}) == 1
    assert {item.model for item in inputs} == set(models)


def test_numeric_comparator_requires_one_score_for_every_frozen_test_site() -> None:
    run = build_literature_random_track(
        _rows(20),
        seeds=(0,),
        test_fraction=0.2,
        validation_fraction_of_remaining=0.2,
        unlabeled_per_positive=1,
        sampling_seed=20260811,
    )[0]
    model_input = bind_model_to_comparison_panel(run, "pcysmod")
    expected = dict(model_input.partition_rows)["test"]
    incomplete = {row.site_key: 0.5 for row in expected[:-1]}

    try:
        assert_complete_model_scores(model_input, "test", incomplete)
    except RuntimeError as exc:
        assert "complete frozen input" in str(exc)
    else:
        raise AssertionError("incomplete pCysMod scores entered numeric table")


def test_tree_and_graft_are_architecture_references_not_competitors() -> None:
    assert competitor_role("tree") == "architecture_reference_only"
    assert competitor_role("graft") == "architecture_reference_only"
    assert competitor_role("pcysmod") == "quantitative_only_with_complete_scores"
    assert competitor_role("sul_bertgru") == "isolated_external_adapter"


def test_sul_adapter_uses_putative_negative_without_rewriting_pu_labels() -> None:
    source_rows = _rows(20)
    run = build_literature_random_track(
        source_rows,
        seeds=(0,),
        test_fraction=0.2,
        validation_fraction_of_remaining=0.2,
        unlabeled_per_positive=1,
        sampling_seed=20260811,
    )[0]
    model_input = bind_model_to_comparison_panel(run, "sul_bertgru")
    sequences = {
        row.global_protein_id: "M" * 2 + "C" + "A" * 3 + "C" + "M" * 20
        for row in source_rows
    }

    adapter_rows = prepare_sul_bertgru_adapter_rows(model_input, sequences)

    assert {row.adapter_label for row in adapter_rows} == {
        "positive",
        "putative_negative",
    }
    assert all(len(row.sequence_window) == 31 for row in adapter_rows)
    assert {row.label for row in source_rows} == {"positive", "unlabeled"}
    assert all(row.sequence_window[15] == "C" for row in adapter_rows)


def test_tabular_direct_baseline_scores_complete_shared_validation_and_test() -> None:
    run = build_literature_random_track(
        _rows(20),
        seeds=(0,),
        test_fraction=0.2,
        validation_fraction_of_remaining=0.2,
        unlabeled_per_positive=1,
        sampling_seed=20260811,
    )[0]
    model_input = bind_model_to_comparison_panel(run, "pu_logistic")
    features = {
        (row.global_protein_id, row.cys_position): (
            float(row.cys_position),
            float(row.cys_position == 3),
        )
        for row in _rows(20)
    }

    scored = run_tabular_direct_baseline(model_input, features)

    for partition in ("validation", "test"):
        assert_complete_model_scores(
            model_input, partition, dict(scored.partition_scores)[partition]
        )


def test_structure_ranker_scores_complete_shared_validation_and_test() -> None:
    rows = _rows(20)
    run = build_literature_random_track(
        rows,
        seeds=(0,),
        test_fraction=0.2,
        validation_fraction_of_remaining=0.2,
        unlabeled_per_positive=1,
        sampling_seed=20260811,
    )[0]
    model_input = bind_model_to_comparison_panel(run, "structure_ranker")
    sequence = {
        (row.global_protein_id, row.cys_position): (
            float(row.cys_position),
            float(row.cys_position == 3),
        )
        for row in rows
    }
    structure = {
        key: ((0.25, 80.0), True) for key in sequence
    }

    scored = run_structure_direct_baseline(
        model_input,
        sequence,
        None,
        structure,
        parameters={
            "epochs": 1,
            "hidden": 4,
            "n_mc_dropout": 1,
            "use_esm": 0,
        },
    )

    for partition in ("validation", "test"):
        assert_complete_model_scores(
            model_input, partition, dict(scored.partition_scores)[partition]
        )


def test_comparable_structure_features_use_frozen_registered_models() -> None:
    features = build_registered_structure_features(
        {
            ("arabidopsis|A0A067Y2H7", 31),
            ("rice|NO_REGISTERED_MODEL", 7),
        },
        registry_path=Path(
            "data/registry/releases/alphafold_structures_release_v2.tsv"
        ),
        registry_base=Path("data/registry"),
    )

    values, available = features[("arabidopsis|A0A067Y2H7", 31)]
    assert available is True
    assert values[1] == 95.62
    assert features[("rice|NO_REGISTERED_MODEL", 7)] == (
        (0.0, 0.0),
        False,
    )


def test_direct_roster_executes_all_five_models_on_one_shared_panel() -> None:
    rows = _rows(20)
    run = build_literature_random_track(
        rows,
        seeds=(0,),
        test_fraction=0.2,
        validation_fraction_of_remaining=0.2,
        unlabeled_per_positive=1,
        sampling_seed=20260811,
    )[0]
    sequence = {
        (row.global_protein_id, row.cys_position): (
            float(row.cys_position),
            float(row.cys_position == 3),
        )
        for row in rows
    }
    esm = {
        key: (values[0], values[1], values[0] + values[1])
        for key, values in sequence.items()
    }
    structure = {key: ((0.25, 80.0), True) for key in sequence}

    scores = run_direct_comparison_roster(
        run,
        sequence_features=sequence,
        esm_features=esm,
        structure_features=structure,
        parameters={
            "pu_logistic": {"holdout_fraction": 0.2},
            "random_forest": {"n_estimators": 2},
            "xgboost": {"n_estimators": 2},
            "structure_ranker": {
                "epochs": 1,
                "hidden": 4,
                "n_mc_dropout": 1,
                "use_esm": 0,
            },
        },
    )

    assert {result.model for result in scores} == {
        "pu_logistic",
        "random_forest",
        "xgboost",
        "esm_linear_head",
        "structure_ranker",
    }
    assert {result.panel_sha256 for result in scores} == {run.panel_sha256}
