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
    ("PXD063170", "41467_2025_61582_MOESM3_ESM.xlsx"): (
        "5da579af1560819863b45c1af6966a53746725b20d5546b115212381b86bec49"
    ),
    ("PXD063170", "Magnaporthe_oryzae.MG8.pep.all.fa"): (
        "0466bf3d2af6fa6a6d44cf74d72a6aa717dc948963ef8b52f7224e70f0470293"
    ),
    ("PXD072089", "pnas.2608150123.sd01.xlsx"): (
        "6318477160ed6b0ca335aff97505d0deab899fcba261dcd76f8c06af614a7cae"
    ),
    ("PXD072089", "pnas.2608150123.sd02.xlsx"): (
        "9c6520e7dca452e8e473b092c4307e0a06263a3a15566201fde90956430e7446"
    ),
    ("PXD072089", "pnas.2608150123.sd03.xlsx"): (
        "f8e8e02925207092590e09893533b3b11dc57926ae6f69d1df73606232addc52"
    ),
    ("PXD072089", "pnas.2608150123.sd04.xlsx"): (
        "834c21a60f950216b39c151b95e576b42c02edb4e9fdfd4af4a22f548833c8fe"
    ),
    ("PXD072089", "pnas.2608150123.sd05.xlsx"): (
        "5d339c68442d40acf938bb04d4ff017ab0b833b94f202bbb1460f6f764160c87"
    ),
    ("PXD072089", "pnas.2608150123.sd06.xlsx"): (
        "4961449856f0ff0b6684648d2a477e5b5545d4cf8b0e68efe0f09fd3d4a5944f"
    ),
    ("PXD072089", "uniprot_rice_v1.fasta"): (
        "fc20c76a58c6e5695e57cf2f389b91cff2af15bd7d23182b1dc1f99123b9322b"
    ),
    ("PXD072089", "SS-all-peptides.tsv"): (
        "57c4ce64b5747eb591e879565916c04e0901ea228c20e511b2e1804fb384890c"
    ),
    ("PXD072300", "Persulfidation_MDHAR5.csv"): (
        "33d1c750f3b17e4140822275b5acb96f328ae5f4fd15ce1f0b02752da665891d"
    ),
    ("PXD072300", "Persulfidation_ALDP.csv"): (
        "5d0d51424cd7860a92d52b8ace52abd96ea51d505be4d9c3162959fb487aff16"
    ),
    ("PXD072300", "Persulfidation_PFP.csv"): (
        "776b1f322787d1e7643e9c089b506363547d72b06c155acdcf98bfc3750ceb95"
    ),
    ("PXD072300", "Persulfidation_FBA1.csv"): (
        "83f4df4750b354ecaeb07f4b434f0b9e238b30079bc67fff0bab92fc9c10501d"
    ),
    ("PXD072300", "Persulfidation_MDHAR4.csv"): (
        "c8eb247b2a94152c94d952ef8f93879eea6834e701c11f40494344ff12906d94"
    ),
    ("PXD072300", "Persulfidation_TKT.csv"): (
        "9ab4628bb514eadfbdcde9fd47d513d89cf32cc845a5e9a8d6e191377bf3adc2"
    ),
    ("PXD072300", "Persulfidation_RPI.csv"): (
        "3c333a3044106876b0df307b68250737376918e47d53d31b5aa1f8de5601ab12"
    ),
    ("PXD072300", "Persulfidation_FBA3.csv"): (
        "508ee10ac5919fdd3c01ded704eae3c7b79cebff360549fc9df979f4134bbdf7"
    ),
    ("PXD072300", "Persulfidation_MDHAR3.csv"): (
        "af4efe6c87af381aad71de87e028b94ea64cfce0e67dfd3298aa1e66156666cd"
    ),
    ("PXD072300", "Persulfidation_TAL.csv"): (
        "67fea3ed43243ecd31f1871f5ffa5aa196a3c1db0b7d8fb8cd46c06d5054b761"
    ),
    ("KIAE271_SUPPL", "kiae271_DSs.xlsx"): (
        "a437f1a941a348ccf15d538bcea193d9feefe29019c5162f0427551416b1b98a"
    ),
    ("GOANNOT_TOMATO", "tomato_go_taxon4081.tsv"): (
        "aef1eb8cd827bb6abd5591fe6b7a68a437963c5756b4cc3b7fc6f66caeea294c"
    ),
    ("PANTHERDB_ORTHOLOGY", "panther_tomato_taxon4081.tsv"): (
        "58cd0b4b3f610e8f54c6a8063f0c914c317e21e226be095b5b2f19de3d03c432"
    ),
    ("PANTHERDB_ORTHOLOGY", "panther_arabidopsis_taxon3702.tsv"): (
        "859156896af8d70d2c3262df31af6da7cd7411e56564b23b2958428348203dde"
    ),
    ("PANTHERDB_ORTHOLOGY", "panther_rice_taxon4530.tsv"): (
        "d362ccf6a109c6e2eca337b3a1799cfba489bb324aeee1e0b45180e659be7f50"
    ),
    ("PANTHERDB_ORTHOLOGY", "panther_magnaporthe_taxon242507.tsv"): (
        "79ff37494f3e5fc42589a863ae062c5d5a1b7d43ef95c997c47ce5559ac55c9d"
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
