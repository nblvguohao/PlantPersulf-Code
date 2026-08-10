import pytest

from plantpersulf.features.site_biology import BIOLOGY_FEATURE_NAMES
from plantpersulf.workflows.cross_crop_target_label_free import (
    SourceBatch,
    validate_target_label_free_sources,
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
