from pathlib import Path

import pytest
import yaml


def test_dataset_record_requires_traceable_source() -> None:
    from plantpersulf.provenance.schema import DatasetRecord

    with pytest.raises(ValueError, match="accession"):
        DatasetRecord(
            accession="",
            repository="PRIDE",
            source_url="",
            scientific_role="training",
            metadata_retrieved_at="",
            metadata_sha256="",
        )


@pytest.mark.parametrize(
    "field",
    [
        "repository",
        "source_url",
        "scientific_role",
        "metadata_retrieved_at",
        "metadata_sha256",
    ],
)
def test_dataset_record_rejects_each_missing_traceability_field(field: str) -> None:
    from plantpersulf.provenance.schema import DatasetRecord

    values = {
        "accession": "PXD006140",
        "repository": "PRIDE",
        "source_url": "https://www.ebi.ac.uk/pride/archive/projects/PXD006140",
        "scientific_role": "training",
        "metadata_retrieved_at": "2026-07-20T00:00:00Z",
        "metadata_sha256": "a" * 64,
    }
    values[field] = ""

    with pytest.raises(ValueError, match=field):
        DatasetRecord(**values)


def test_dataset_record_requires_valid_metadata_sha256() -> None:
    from plantpersulf.provenance.schema import DatasetRecord

    with pytest.raises(ValueError, match="metadata_sha256"):
        DatasetRecord(
            accession="PXD006140",
            repository="PRIDE",
            source_url="https://www.ebi.ac.uk/pride/archive/projects/PXD006140",
            scientific_role="training",
            metadata_retrieved_at="2026-07-20T00:00:00Z",
            metadata_sha256="not-a-sha256",
        )


def test_scientific_integrity_config_is_fail_closed() -> None:
    config_path = Path("configs/scientific_integrity.yaml")
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))

    assert config["biological_data"]["allow_synthetic"] is False
    assert config["biological_data"]["allow_random_values"] is False
    assert config["biological_data"]["allow_simulated_labels"] is False
    assert config["biological_data"]["require_registered_input"] is True
    assert config["labels"]["unobserved_cysteine"] == "unlabeled"
    assert config["splits"]["allow_test_tuning"] is False
    assert config["splits"]["allow_frozen_split_overwrite"] is False
    assert config["known_mechanisms"]["allow_train_and_independent_validation"] is False
    assert config["prospective_validation"]["allow_post_wetlab_reranking"] is False
    assert config["release"]["allow_placeholder_results"] is False
    assert set(config["provenance"]["required_dataset_fields"]) == {
        "accession",
        "repository",
        "source_url",
        "scientific_role",
        "metadata_retrieved_at",
        "metadata_sha256",
    }
    assert set(config["release"]["forbidden_filename_terms"]) >= {
        "synthetic",
        "fake",
        "dummy",
        "mock_biology",
    }
