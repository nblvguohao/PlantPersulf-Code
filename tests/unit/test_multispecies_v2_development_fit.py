"""Fold-local fitting contract for the v2 development runner."""

from __future__ import annotations

from plantpersulf.proteomics.multispecies_v2_dataset import MultispeciesV2SiteRow
from plantpersulf.workflows.multispecies_v2 import (
    fit_v2_development_pipeline,
    select_v2_development_hyperparameters,
    train_v2_development_fold,
)


def _row(accession: str, label: str) -> MultispeciesV2SiteRow:
    return MultispeciesV2SiteRow(
        species="arabidopsis",
        protein_accession=accession,
        cys_position=3,
        label=label,
        study_accessions=("STUDY",) if label == "positive" else (),
        global_protein_id=f"arabidopsis|{accession}",
        cluster_id=accession,
        split="development",
        development_fold=0,
    )


def test_v2_pipeline_fits_imputation_selection_scaling_and_pu_prior_from_fit_rows() -> (
    None
):
    """Changing validation features cannot alter any fitted development state."""
    fit_rows = (_row("FIT_POS", "positive"), _row("FIT_UNL", "unlabeled"))
    validation_rows = (_row("VALID", "unlabeled"),)
    features = {
        ("arabidopsis|FIT_POS", 3): (1.0, None, 5.0),
        ("arabidopsis|FIT_UNL", 3): (3.0, 7.0, 5.0),
        ("arabidopsis|VALID", 3): (999.0, 999.0, 5.0),
    }

    pipeline = fit_v2_development_pipeline(fit_rows, features)
    transformed = pipeline.transform(validation_rows, features)

    assert pipeline.imputation_values == (2.0, 7.0, 5.0)
    assert pipeline.selected_feature_indices == (0,)
    assert pipeline.pu_prior == 0.5
    assert transformed != [[999.0, 999.0]]


def test_v2_pipeline_preserves_nonzero_standard_deviation_below_one() -> None:
    fit_rows = (_row("LOW", "positive"), _row("HIGH", "unlabeled"))
    features = {
        ("arabidopsis|LOW", 3): (0.0,),
        ("arabidopsis|HIGH", 3): (0.5,),
    }

    pipeline = fit_v2_development_pipeline(fit_rows, features)

    assert pipeline.feature_std == (0.25,)
    assert pipeline.transform(fit_rows, features) == [[-1.0], [1.0]]


def test_v2_fold_trainer_accepts_only_fit_and_validation_rows() -> None:
    """The development model-selection API has no frozen-test argument."""
    fit_rows = (
        _row("FIT_POS_A", "positive"),
        _row("FIT_POS_B", "positive"),
        _row("FIT_UNL_A", "unlabeled"),
        _row("FIT_UNL_B", "unlabeled"),
    )
    validation_rows = (
        _row("VALID_POS", "positive"),
        _row("VALID_UNL", "unlabeled"),
    )
    features = {
        (row.global_protein_id, row.cys_position): (
            float(index),
            float(index % 2),
        )
        for index, row in enumerate((*fit_rows, *validation_rows), start=1)
    }

    result = train_v2_development_fold(fit_rows, validation_rows, features, seed=0)

    assert 0.0 <= result.validation_ap <= 1.0
    assert result.pipeline.pu_prior == 0.5


def test_v2_hyperparameter_selection_has_no_frozen_test_argument() -> None:
    """Selection must be structurally limited to fit and validation rows."""
    fit_rows = (_row("P1", "positive"), _row("P2", "positive"), _row("U1", "unlabeled"))
    validation_rows = (_row("VP", "positive"), _row("VU", "unlabeled"))
    features = {
        (row.global_protein_id, row.cys_position): (float(index), float(index + 1))
        for index, row in enumerate((*fit_rows, *validation_rows))
    }

    selected = select_v2_development_hyperparameters(
        fit_rows, validation_rows, features, candidates=(0.1, 0.2), seed=0
    )

    assert selected in (0.1, 0.2)
