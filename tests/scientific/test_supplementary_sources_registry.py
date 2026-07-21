"""RED (Phase A2 pipeline): register persulfidation-site supplement sources.

The two Gate-1 site-level studies draw evidence from supplement/search-output
files that are not the primary PRIDE result files: PXD006140 Dataset S3 (paper
supplement) and PXD024061's MaxQuant Sulfide(C)/CianoBiotin(C) site tables
(deposited search output). This registry records each with full provenance so the
site parsers only ever read hash-verified, registered inputs.

Expected RED: ``plantpersulf.provenance.supplementary`` does not exist yet.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from plantpersulf.provenance.supplementary import (  # RED: module missing
    SUPPLEMENTARY_FIELDS,
    audit_supplementary_sources,
)

CONFIG = Path("configs/supplementary_sources_v1.yaml")
REGISTRY = Path("data/registry/supplementary_sources.tsv")

EXPECTED = {
    ("PXD006140", "erx294_suppl_supplementary_data_set_s3.xlsx"): (
        "b4831646643204fca48dc6b13e76e68f5a8fe88fa999fc00884e5c0d4e1ded87"
    ),
    ("PXD024061", "txt_persulfproject/Sulfide(C)Sites.txt"): (
        "e9a1abfb5b38a75fe64626f707a0fc04917d628d1f3437d807982065891003d8"
    ),
    ("PXD024061", "txt_persulfproject/CianoBiotin(C)Sites.txt"): (
        "a98a49e91ad41dbfa53f91a5de32629a7591b51f59d209f2f696e6b9365d00fb"
    ),
}


def test_registry_schema_has_provenance_fields() -> None:
    for field in (
        "study_accession",
        "member",
        "publication_doi",
        "sha256",
        "data_level",
        "scientific_use",
    ):
        assert field in SUPPLEMENTARY_FIELDS


def test_all_site_supplement_sources_audit_against_local_sha256() -> None:
    sources = audit_supplementary_sources(CONFIG, REGISTRY)

    keyed = {(s.study_accession, s.member): s for s in sources}
    assert set(keyed) == set(EXPECTED)
    for key, sha256 in EXPECTED.items():
        assert keyed[key].sha256 == sha256
        assert keyed[key].data_level == "C"
        assert keyed[key].local_path.is_file()


def test_both_pxd024061_members_are_registered_distinctly() -> None:
    sources = audit_supplementary_sources(CONFIG, REGISTRY)
    pxd024061 = [s for s in sources if s.study_accession == "PXD024061"]

    # Two distinct members of the same study must both survive (not deduped).
    assert len(pxd024061) == 2
    assert {s.member for s in pxd024061} == {
        "txt_persulfproject/Sulfide(C)Sites.txt",
        "txt_persulfproject/CianoBiotin(C)Sites.txt",
    }


def test_tampered_registry_sha256_is_rejected(tmp_path: Path) -> None:
    tampered = tmp_path / "supplementary_sources.tsv"
    original = REGISTRY.read_text(encoding="utf-8")
    tampered.write_text(original.replace("b48316466", "0" * 9), encoding="utf-8")

    with pytest.raises(RuntimeError):
        audit_supplementary_sources(CONFIG, tampered)
