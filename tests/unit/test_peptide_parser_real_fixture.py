import csv
import hashlib
import json
from pathlib import Path

import pytest

from plantpersulf.provenance.hashing import hash_file

FIXTURE_MANIFEST = Path("tests/fixtures/real/PXD006140/source_manifest.json")
DOWNLOAD_REGISTRY = Path("data/registry/downloads.tsv")
REFERENCE_REGISTRY = Path("data/registry/reference_sequences.tsv")
REFERENCE_PATHS = (
    Path("data/registry/cache/uniprot/Q93VK9.fasta"),
    Path("data/registry/cache/uniprot/Q9ZW96.fasta"),
)
OUTPUT_FILES = {
    "psms.tsv",
    "sites.tsv",
    "conflicts.tsv",
    "excluded.tsv",
    "missing_metadata.tsv",
    "manifest.json",
}


def _write_registry(path: Path, source_path: Path, sha256: str) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("path", "sha256"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow(
            {
                "path": source_path.resolve().as_posix(),
                "sha256": sha256,
            }
        )


def _fixture_source(tmp_path: Path):  # type: ignore[no-untyped-def]
    from plantpersulf.proteomics.metadata import ProteomicsSource

    manifest = json.loads(FIXTURE_MANIFEST.read_text(encoding="utf-8"))
    fixture_path = FIXTURE_MANIFEST.parent / manifest["fixture_file"]
    registry_path = tmp_path / "fixture_registry.tsv"
    _write_registry(registry_path, fixture_path, manifest["fixture_sha256"])
    return ProteomicsSource(
        study_accession="PXD006140",
        source_file=fixture_path,
        registry_file=registry_path,
        source_sha256=manifest["fixture_sha256"],
    )


def _tree_hashes(root: Path) -> dict[str, str]:
    return {
        path.relative_to(root).as_posix(): hashlib.sha256(
            path.read_bytes()
        ).hexdigest()
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


def test_real_omssa_fixture_preserves_peptide_and_modification_bytes(
    tmp_path: Path,
) -> None:
    from plantpersulf.proteomics.peptide_parser import parse_omssa

    parsed = parse_omssa(_fixture_source(tmp_path))
    q93vk9 = next(record for record in parsed.psms if record.spectrum_id == "29")

    assert q93vk9.modified_sequence == "ETmASLGLICEK"
    assert q93vk9.peptide_sequence == "ETMASLGLICEK"
    assert q93vk9.modification_name_raw == "oxidation of M:3"
    assert q93vk9.protein_accession_raw == "Q93VK9"
    assert q93vk9.protein_accession_canonical == "Q93VK9"
    assert q93vk9.start_one_based == 15
    assert q93vk9.stop_one_based == 26
    assert q93vk9.source_sha256 == (
        "10fce0b8ac26f564c7c16879b0d6e32cfa14a8ee1116f0f7911cbc758650b43f"
    )


def test_real_decoy_and_missing_sample_are_audited(tmp_path: Path) -> None:
    from plantpersulf.proteomics.peptide_parser import parse_omssa

    parsed = parse_omssa(_fixture_source(tmp_path))

    assert any(
        issue.spectrum_id == "32" and issue.reason == "decoy"
        for issue in parsed.excluded
    )
    assert all(record.spectrum_id != "32" for record in parsed.psms)
    assert all(record.sample_id == "" for record in parsed.psms)
    assert [issue.reason for issue in parsed.missing] == [
        "missing_sample_metadata"
    ]


def test_real_multi_protein_peptide_mappings_are_preserved() -> None:
    from plantpersulf.proteomics.metadata import ProteomicsSource
    from plantpersulf.proteomics.peptide_parser import parse_omssa

    source_path = Path(
        "data/raw/PXD006140/"
        "omssa.ne.20150821_01_AAroca_TMT6plex.cmpd.mgf.txt"
    )
    parsed = parse_omssa(
        ProteomicsSource(
            study_accession="PXD006140",
            source_file=source_path,
            registry_file=DOWNLOAD_REGISTRY,
            source_sha256=(
                "9bbf2f88f51ab9a80b62315dd2b61f86f0d4290d8457591d28912a252af971ea"
            ),
        )
    )

    accessions = {
        record.protein_accession_raw
        for record in parsed.psms
        if record.peptide_sequence == "GNESYEDAIEALKK"
    }
    assert accessions == {"P42737-2", "F4K875"}


def test_invalid_omssa_header_fails_before_biological_output(
    tmp_path: Path,
) -> None:
    from plantpersulf.proteomics.metadata import ProteomicsSource
    from plantpersulf.proteomics.peptide_parser import parse_omssa

    source_path = tmp_path / "invalid.csv"
    source_path.write_text("not_an_omssa_header\n", encoding="utf-8")
    source_sha256 = hash_file(source_path, "sha256")
    registry_path = tmp_path / "registry.tsv"
    _write_registry(registry_path, source_path, source_sha256)
    source = ProteomicsSource(
        study_accession="POLICYTEST",
        source_file=source_path,
        registry_file=registry_path,
        source_sha256=source_sha256,
    )

    with pytest.raises(RuntimeError, match="OMSSA source has invalid columns"):
        parse_omssa(source)


def test_cli_exposes_task4_commands() -> None:
    from plantpersulf.cli import build_parser

    parse_args = build_parser().parse_args(
        ["parse-proteomics", "--accession", "PXD006140"]
    )
    audit_args = build_parser().parse_args(
        ["audit-sites", "--accession", "PXD006140"]
    )

    assert parse_args.command == "parse-proteomics"
    assert audit_args.command == "audit-sites"


def test_fixture_pipeline_publication_is_deterministic(tmp_path: Path) -> None:
    from plantpersulf.proteomics.peptide_parser import (
        publish_parsed_proteomics,
    )

    source = _fixture_source(tmp_path)
    first_output = tmp_path / "first/PXD006140/proteomics_parser_v1"
    second_output = tmp_path / "second/PXD006140/proteomics_parser_v1"

    first = publish_parsed_proteomics(
        (source,),
        REFERENCE_PATHS,
        REFERENCE_REGISTRY,
        first_output,
    )
    second = publish_parsed_proteomics(
        (source,),
        REFERENCE_PATHS,
        REFERENCE_REGISTRY,
        second_output,
    )

    assert first == second
    assert first.site_count == 3
    assert {path.name for path in first_output.iterdir()} == OUTPUT_FILES
    assert _tree_hashes(first_output) == _tree_hashes(second_output)


def test_fixture_pipeline_rejects_changed_source_sha_without_output(
    tmp_path: Path,
) -> None:
    from plantpersulf.proteomics.metadata import ProteomicsSource
    from plantpersulf.proteomics.peptide_parser import (
        publish_parsed_proteomics,
    )

    manifest = json.loads(FIXTURE_MANIFEST.read_text(encoding="utf-8"))
    fixture_path = FIXTURE_MANIFEST.parent / manifest["fixture_file"]
    registry_path = tmp_path / "wrong_registry.tsv"
    _write_registry(registry_path, fixture_path, "0" * 64)
    source = ProteomicsSource(
        study_accession="PXD006140",
        source_file=fixture_path,
        registry_file=registry_path,
        source_sha256=manifest["fixture_sha256"],
    )
    output = tmp_path / "output"

    with pytest.raises(RuntimeError, match="SHA256 mismatch"):
        publish_parsed_proteomics(
            (source,),
            REFERENCE_PATHS,
            REFERENCE_REGISTRY,
            output,
        )

    assert not output.exists()
