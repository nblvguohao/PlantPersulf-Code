from pathlib import Path

import yaml

SOURCE_CONFIG = Path("configs/evidence_method_sources_v1.yaml")
MAPPING_CONFIG = Path("configs/evidence_method_mappings_v1.yaml")
METHOD_REGISTRY = Path("data/registry/evidence_methods.tsv")


def test_pxd035795_method_source_must_be_registered_and_hashed() -> None:
    from plantpersulf.evidence.methods import audit_method_sources

    sources = audit_method_sources(SOURCE_CONFIG, METHOD_REGISTRY)
    source = next(row for row in sources if row.study_accession == "PXD035795")

    assert source.identifier_type == "DOI"
    assert source.identifier == "10.1111/nph.18838"
    assert source.repository == "University of Seville IDUS"
    assert source.repository_record_id == ("ce5a8e5e-6bf0-4564-a0d5-c6dcb4eeaa6f")
    assert source.size_bytes == 11_657_462
    assert len(source.sha256) == 64
    assert source.local_path.is_file()
    assert source.method_scope == "dimedone_switch_lc_ms_ms"
    assert source.data_use_status == "method_audit_only"


def test_version_one_has_no_reviewed_site_upgrade_mapping() -> None:
    mapping = yaml.safe_load(MAPPING_CONFIG.read_text(encoding="utf-8"))

    assert mapping == {"version": 1, "mappings": []}


def test_cli_exposes_method_acquisition_and_audit() -> None:
    from plantpersulf.cli import build_parser

    acquire = build_parser().parse_args(
        ["acquire-evidence-method", "--accession", "PXD035795"]
    )
    audit = build_parser().parse_args(["audit-evidence-methods"])

    assert acquire.command == "acquire-evidence-method"
    assert audit.command == "audit-evidence-methods"
