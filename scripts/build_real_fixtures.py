"""Build deterministic tests-only fixtures from audited real public files."""

from __future__ import annotations

import argparse
import json
import tempfile
from collections.abc import Mapping
from dataclasses import asdict, dataclass
from pathlib import Path
from types import MappingProxyType

from plantpersulf.provenance.audit import assert_registered_input
from plantpersulf.provenance.hashing import hash_file

EXTRACTION_COMMAND = (
    "python scripts/build_real_fixtures.py --accessions PXD006140,PXD051570,GSE163745"
)


@dataclass(frozen=True)
class FixtureSource:
    accession: str
    source_repository: str
    source_registry: Path
    source_file: Path
    expected_source_sha256: str
    fixture_name: str
    extraction_method: str
    physical_line_count: int | None = None


@dataclass(frozen=True)
class FixtureBuildSummary:
    accession_count: int
    fixture_count: int


SOURCES: Mapping[str, FixtureSource] = MappingProxyType(
    {
        "PXD006140": FixtureSource(
            accession="PXD006140",
            source_repository="PRIDE",
            source_registry=Path("data/registry/downloads.tsv"),
            source_file=Path(
                "data/raw/PXD006140/omssa.20150821_01_AAroca_TMT6plex.cmpd.mgf.txt"
            ),
            expected_source_sha256=(
                "f8d052626f9f792495c785b7b42d657e2effb3a982c2f8ba561c7c320f37bdcb"
            ),
            fixture_name="omssa_head_32.csv",
            extraction_method="binary_header_plus_32_lines",
            physical_line_count=33,
        ),
        "PXD051570": FixtureSource(
            accession="PXD051570",
            source_repository="iProX",
            source_registry=Path("data/registry/files.tsv"),
            source_file=Path("data/registry/cache/iprox/PXD051570.json"),
            expected_source_sha256=(
                "84a5d83e099910a780721a3984c753a35f1e25f781675f5a2e5097174a336168"
            ),
            fixture_name="PXD051570.metadata.json",
            extraction_method="byte_for_byte_copy",
        ),
        "GSE163745": FixtureSource(
            accession="GSE163745",
            source_repository="GEO",
            source_registry=Path("data/registry/files.tsv"),
            source_file=Path("data/registry/cache/geo/GSE163745.json"),
            expected_source_sha256=(
                "20f52bd62689003b8d7185167c6e18ec70b9e1178da3342134a447c138c5f228"
            ),
            fixture_name="GSE163745.metadata.json",
            extraction_method="byte_for_byte_copy",
        ),
    }
)


def _validate_existing_manifest(path: Path) -> None:
    if not path.exists():
        return
    try:
        parsed: object = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"invalid existing fixture manifest: {path}") from exc
    if not isinstance(parsed, dict) or parsed.get("schema_version") != 1:
        raise RuntimeError(f"unsupported existing fixture manifest: {path}")


def _copy_source(source: FixtureSource, temporary_path: Path) -> None:
    with source.source_file.open("rb") as input_handle:
        with temporary_path.open("wb") as output_handle:
            if source.extraction_method == "byte_for_byte_copy":
                for chunk in iter(lambda: input_handle.read(1024 * 1024), b""):
                    output_handle.write(chunk)
                return
            if source.extraction_method == "binary_header_plus_32_lines":
                if source.physical_line_count is None:
                    raise RuntimeError("physical line count is required")
                for _ in range(source.physical_line_count):
                    line = input_handle.readline()
                    if not line:
                        raise RuntimeError(
                            f"fixture source has too few lines: {source.accession}"
                        )
                    output_handle.write(line)
                return
    raise RuntimeError(f"unsupported extraction method: {source.extraction_method}")


def _temporary_path(directory: Path, name: str) -> Path:
    with tempfile.NamedTemporaryFile(
        dir=directory,
        prefix=f".{name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        return Path(handle.name)


def build_fixture(
    source: FixtureSource,
    output_root: Path = Path("tests/fixtures/real"),
    extraction_command: str = EXTRACTION_COMMAND,
) -> Path:
    """Build one fixture only after its source passes registered SHA256 audit."""
    assert_registered_input(source.source_file, source.source_registry)
    source_sha256 = hash_file(source.source_file, "sha256")
    if source_sha256 != source.expected_source_sha256:
        raise RuntimeError(
            "fixture source SHA256 differs from approved Task 3 source: "
            f"{source.accession}"
        )

    output_directory = output_root / source.accession
    fixture_path = output_directory / source.fixture_name
    manifest_path = output_directory / "source_manifest.json"
    _validate_existing_manifest(manifest_path)
    output_directory.mkdir(parents=True, exist_ok=True)
    temporary_fixture_path = _temporary_path(
        output_directory,
        source.fixture_name,
    )
    temporary_fixture: Path | None = temporary_fixture_path
    temporary_manifest: Path | None = None
    try:
        _copy_source(source, temporary_fixture_path)
        fixture_sha256 = hash_file(temporary_fixture_path, "sha256")
        manifest: dict[str, object] = {
            "schema_version": 1,
            "accession": source.accession,
            "source_repository": source.source_repository,
            "source_registry": source.source_registry.as_posix(),
            "source_file": source.source_file.as_posix(),
            "source_sha256": source_sha256,
            "fixture_file": source.fixture_name,
            "fixture_sha256": fixture_sha256,
            "extraction_method": source.extraction_method,
            "extraction_command": extraction_command,
            "biological_values_modified": False,
            "use": "tests_only",
        }
        temporary_manifest_path = _temporary_path(
            output_directory,
            "source_manifest.json",
        )
        temporary_manifest = temporary_manifest_path
        serialized = json.dumps(
            manifest,
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        )
        temporary_manifest_path.write_bytes(f"{serialized}\n".encode())
        temporary_fixture_path.replace(fixture_path)
        temporary_fixture = None
        temporary_manifest_path.replace(manifest_path)
        temporary_manifest = None
        if hash_file(fixture_path, "sha256") != fixture_sha256:
            raise RuntimeError(f"published fixture SHA256 mismatch: {fixture_path}")
        return manifest_path
    finally:
        if temporary_fixture is not None:
            temporary_fixture.unlink(missing_ok=True)
        if temporary_manifest is not None:
            temporary_manifest.unlink(missing_ok=True)


def build_real_fixtures(
    accessions: tuple[str, ...],
    output_root: Path = Path("tests/fixtures/real"),
) -> FixtureBuildSummary:
    """Build the requested approved accessions in the supplied order."""
    normalized = tuple(accession.strip().upper() for accession in accessions)
    if any(not accession for accession in normalized):
        raise RuntimeError("fixture accession cannot be empty")
    if len(set(normalized)) != len(normalized):
        raise RuntimeError("duplicate fixture accession")
    unsupported = [accession for accession in normalized if accession not in SOURCES]
    if unsupported:
        raise RuntimeError(f"unsupported fixture accession: {unsupported[0]}")
    for accession in normalized:
        build_fixture(SOURCES[accession], output_root=output_root)
    return FixtureBuildSummary(
        accession_count=len(normalized),
        fixture_count=len(normalized),
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accessions", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = build_parser().parse_args(argv)
    accessions = tuple(arguments.accessions.split(","))
    summary = build_real_fixtures(accessions)
    print(json.dumps(asdict(summary), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
