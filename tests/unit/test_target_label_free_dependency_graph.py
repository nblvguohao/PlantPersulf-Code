import pytest

from plantpersulf.features.site_biology import BIOLOGY_FEATURE_NAMES
from plantpersulf.workflows.cross_crop_target_label_free import (
    SourceBatch,
    validate_source_coverage,
    validate_target_label_free_sources,
)


def _source(species: str, study: str) -> SourceBatch:
    return SourceBatch(
        species=species,
        study_accession=study,
        feature_names=BIOLOGY_FEATURE_NAMES,
        features=((0.0,) * 6, (1.0,) * 6),
        labels=("positive", "unlabeled"),
        protein_ids=("group-a", "group-a"),
        source_sha256=f"marker-{species}-{study}",
    )


def test_tomato_source_batch_is_rejected() -> None:
    source = SourceBatch(
        species="Solanum lycopersicum",
        study_accession="policy-marker-study",
        feature_names=BIOLOGY_FEATURE_NAMES,
        features=((0.0,) * 6, (1.0,) * 6),
        labels=("positive", "unlabeled"),
        protein_ids=("group-a", "group-a"),
        source_sha256="marker-sha",
    )
    with pytest.raises(RuntimeError, match="target label leakage"):
        validate_target_label_free_sources((source,), "Solanum lycopersicum")


def test_three_source_species_allow_single_study_species_with_four_batches() -> None:
    sources = (
        _source("Arabidopsis thaliana", "study-at-1"),
        _source("Arabidopsis thaliana", "study-at-2"),
        _source("Oryza sativa", "study-os-1"),
        _source("Magnaporthe oryzae", "study-mo-1"),
    )
    validate_source_coverage(
        sources,
        minimum_source_species=3,
        minimum_total_source_studies=4,
        minimum_source_studies_per_species=1,
    )


def test_two_source_species_fail_the_multispecies_coverage_gate() -> None:
    sources = (
        _source("Arabidopsis thaliana", "study-at-1"),
        _source("Arabidopsis thaliana", "study-at-2"),
        _source("Oryza sativa", "study-os-1"),
        _source("Oryza sativa", "study-os-2"),
    )
    with pytest.raises(RuntimeError, match="source species"):
        validate_source_coverage(sources, 3, 4, 1)
