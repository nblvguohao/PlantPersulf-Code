"""RED (Phase A2 pipeline): deterministic publish/audit of persulfidation sites.

Turns the two site parsers into registered, reproducible on-disk outputs
(sites/low_confidence/conflicts/excluded TSV + manifest) that the readiness v2
gate can consume. Both studies are driven from
``configs/persulfidation_sites_v1.yaml`` against the registered supplement
sources and a SHA256-pinned UniProt reference bundle.

Expected RED: ``plantpersulf.proteomics.persulfidation_publish`` does not exist.
"""

from __future__ import annotations

import json
from pathlib import Path

from plantpersulf.proteomics.persulfidation_publish import (  # RED: module missing
    OUTPUT_NAMES,
    PERSULF_SITE_FIELDS,
    audit_persulfidation_sites,
    publish_persulfidation_sites,
)

CONFIG = Path("configs/persulfidation_sites_v1.yaml")
REGISTRY = Path("data/registry")


def _publish(study: str, root: Path) -> object:
    return publish_persulfidation_sites(
        study_accession=study,
        config_path=CONFIG,
        registry_dir=REGISTRY,
        output_root=root,
    )


def test_pxd024061_publish_counts_match_full_scale(tmp_path: Path) -> None:
    summary = _publish("PXD024061", tmp_path)

    assert summary.site_count == 73
    assert summary.distinct_site_count == 73
    assert summary.low_confidence_count == 8
    assert summary.conflict_count == 5
    assert summary.excluded_count == 2


def test_pxd006140_publish_counts_match_full_scale(tmp_path: Path) -> None:
    summary = _publish("PXD006140", tmp_path)

    assert summary.site_count == 405  # spectra (not deduped)
    assert summary.distinct_site_count == 320
    assert summary.conflict_count == 54


def test_published_outputs_are_site_ms_and_label_free(tmp_path: Path) -> None:
    _publish("PXD024061", tmp_path)
    output = tmp_path / "PXD024061" / "persulfidation_sites_v1"

    assert {p.name for p in output.iterdir()} == set(OUTPUT_NAMES)
    header = (output / "sites.tsv").read_text(encoding="utf-8").splitlines()[0]
    assert header.split("\t") == list(PERSULF_SITE_FIELDS)
    assert "label" not in header.lower()
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["labels_created"] is False
    assert manifest["nondetection_labeled_negative"] is False
    assert manifest["biological_values_modified"] is False


def test_publish_is_deterministic_and_audits(tmp_path: Path) -> None:
    first = _publish("PXD024061", tmp_path)
    # Re-publishing into the existing directory must reproduce byte-identical
    # outputs (the publisher verifies the tree hashes), then audit re-verifies.
    second = _publish("PXD024061", tmp_path)
    audited = audit_persulfidation_sites(
        "PXD024061", tmp_path / "PXD024061" / "persulfidation_sites_v1"
    )

    assert first == second == audited
