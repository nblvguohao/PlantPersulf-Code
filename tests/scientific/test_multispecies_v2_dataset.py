"""Task 9.3 split-bound, positive-unlabelled v2 dataset policy."""

from __future__ import annotations

import pytest

from plantpersulf.benchmark.multispecies_splits import (
    DEVELOPMENT_SPLIT,
    TEST_SPLIT,
    FrozenMultispeciesSplit,
    FrozenSplitRow,
    GlobalClusterRow,
    MultispeciesSiteRow,
    global_protein_id,
)
from plantpersulf.proteomics.multispecies_v2_dataset import (
    MultispeciesV2SiteRow,
    build_split_bound_v2_rows,
    prepare_v2_development_fold,
    sample_v2_unlabeled_panels,
)


def _cluster(species: str, accession: str, cluster: str) -> GlobalClusterRow:
    return GlobalClusterRow(
        species=species,
        protein_accession=accession,
        global_protein_id=global_protein_id(species, accession),
        cluster_id=cluster,
        identity_threshold=0.3,
        source_proteome_sha256="a" * 64,
    )


def _frozen() -> FrozenMultispeciesSplit:
    return FrozenMultispeciesSplit(
        rows=(
            FrozenSplitRow(
                global_protein_id="arabidopsis|P1",
                cluster_id="C1",
                split=DEVELOPMENT_SPLIT,
                development_fold=0,
                species="arabidopsis",
                source_sha256="b" * 64,
            ),
            FrozenSplitRow(
                global_protein_id="tomato|P1",
                cluster_id="C1",
                split=DEVELOPMENT_SPLIT,
                development_fold=0,
                species="tomato",
                source_sha256="b" * 64,
            ),
            FrozenSplitRow(
                global_protein_id="rice|P2",
                cluster_id="C2",
                split=TEST_SPLIT,
                development_fold=None,
                species="rice",
                source_sha256="b" * 64,
            ),
        ),
        seed=20260811,
        test_fraction=0.2,
        n_development_folds=5,
        status="OK",
        test_cluster_fraction=0.2,
        positive_test_fraction_by_stratum={},
    )


def test_split_bound_rows_merge_duplicate_study_evidence_and_keep_all_cys_together(
) -> None:
    """A site and every unlabeled Cys of its protein inherit its one split."""
    rows = build_split_bound_v2_rows(
        positives=(
            MultispeciesSiteRow("arabidopsis", "P1", 3, "positive", "PXD_A"),
            MultispeciesSiteRow("arabidopsis", "P1", 3, "positive", "PXD_B"),
        ),
        all_cysteines=(
            ("arabidopsis", "P1", 3),
            ("arabidopsis", "P1", 9),
            ("tomato", "P1", 4),
            ("rice", "P2", 7),
            ("rice", "P2", 11),
        ),
        clusters=(
            _cluster("arabidopsis", "P1", "C1"),
            _cluster("tomato", "P1", "C1"),
            _cluster("rice", "P2", "C2"),
        ),
        frozen_split=_frozen(),
    )

    by_key = {
        (row.species, row.protein_accession, row.cys_position): row for row in rows
    }
    assert by_key[("arabidopsis", "P1", 3)].study_accessions == ("PXD_A", "PXD_B")
    assert by_key[("arabidopsis", "P1", 9)].label == "unlabeled"
    assert {row.split for row in rows if row.global_protein_id == "arabidopsis|P1"} == {
        DEVELOPMENT_SPLIT
    }
    assert "rice|P2" not in {row.global_protein_id for row in rows}
    assert all(
        row.cluster_id == "C1" for row in rows if row.global_protein_id == "tomato|P1"
    )


def test_split_bound_rows_reject_duplicate_global_cluster_records() -> None:
    """A duplicate key must not be silently overwritten before assignment."""
    with pytest.raises(RuntimeError, match="duplicate global cluster"):
        build_split_bound_v2_rows(
            positives=(),
            all_cysteines=(("arabidopsis", "P1", 3),),
            clusters=(
                _cluster("arabidopsis", "P1", "C1"),
                _cluster("arabidopsis", "P1", "C2"),
            ),
            frozen_split=_frozen(),
        )


def test_development_builder_rejects_any_frozen_test_positive_label() -> None:
    """Supplying a test positive proves that a locked label was read."""
    with pytest.raises(RuntimeError, match="frozen test positive"):
        build_split_bound_v2_rows(
            positives=(
                MultispeciesSiteRow("rice", "P2", 7, "positive", "PXD_R"),
            ),
            all_cysteines=(("rice", "P2", 7),),
            clusters=(_cluster("rice", "P2", "C2"),),
            frozen_split=_frozen(),
        )


def test_split_bound_rows_fail_when_any_cys_protein_lacks_cluster_membership() -> None:
    """No singleton fallback is allowed for an unlabeled protein either."""
    with pytest.raises(RuntimeError, match="missing global cluster"):
        build_split_bound_v2_rows(
            positives=(),
            all_cysteines=(("arabidopsis", "NO_CLUSTER", 3),),
            clusters=(),
            frozen_split=_frozen(),
        )


def test_v2_panel_sampling_is_split_local_and_deterministic() -> None:
    """An unlabeled Cys cannot be drawn by both development and frozen test."""
    rows = tuple(
        MultispeciesV2SiteRow(
            species="arabidopsis",
            protein_accession=f"P{i}",
            cys_position=3,
            label="positive",
            study_accessions=("PXD_A",),
            global_protein_id=f"arabidopsis|P{i}",
            cluster_id=f"C{i}",
            split=DEVELOPMENT_SPLIT if i < 2 else TEST_SPLIT,
            development_fold=0 if i < 2 else None,
        )
        for i in range(3)
    ) + tuple(
        MultispeciesV2SiteRow(
            species="arabidopsis",
            protein_accession=f"U{i}",
            cys_position=position,
            label="unlabeled",
            study_accessions=(),
            global_protein_id=f"arabidopsis|U{i}",
            cluster_id=f"UC{i}",
            split=DEVELOPMENT_SPLIT if i < 40 else TEST_SPLIT,
            development_fold=0 if i < 40 else None,
        )
        for i in range(60)
        for position in (3,)
    )

    first = sample_v2_unlabeled_panels(rows, per_positive=20, seed=12345)
    second = sample_v2_unlabeled_panels(rows, per_positive=20, seed=12345)

    assert first == second
    assert sum(row.label == "positive" for row in first) == 3
    development_unlabeled = [
        row
        for row in first
        if row.split == DEVELOPMENT_SPLIT and row.label == "unlabeled"
    ]
    test_unlabeled = [
        row for row in first if row.split == TEST_SPLIT and row.label == "unlabeled"
    ]
    assert len(development_unlabeled) == 40
    assert len(test_unlabeled) == 20
    assert {row.global_protein_id for row in development_unlabeled}.isdisjoint(
        {row.global_protein_id for row in test_unlabeled}
    )


def test_development_fold_preparation_excludes_test_and_validation_from_fit_rows() -> (
    None
):
    """Fold-local fitting receives only development proteins outside validation."""
    rows = (
        MultispeciesV2SiteRow(
            "arabidopsis",
            "TRAIN",
            3,
            "positive",
            ("A",),
            "arabidopsis|TRAIN",
            "C1",
            DEVELOPMENT_SPLIT,
            1,
        ),
        MultispeciesV2SiteRow(
            "arabidopsis",
            "VALID",
            5,
            "positive",
            ("A",),
            "arabidopsis|VALID",
            "C2",
            DEVELOPMENT_SPLIT,
            0,
        ),
        MultispeciesV2SiteRow(
            "arabidopsis",
            "TEST",
            7,
            "positive",
            ("A",),
            "arabidopsis|TEST",
            "C3",
            TEST_SPLIT,
            None,
        ),
    )

    prepared = prepare_v2_development_fold(
        rows, validation_fold=0, n_development_folds=5
    )

    assert [row.global_protein_id for row in prepared.fit_rows] == ["arabidopsis|TRAIN"]
    assert [row.global_protein_id for row in prepared.validation_rows] == [
        "arabidopsis|VALID"
    ]
