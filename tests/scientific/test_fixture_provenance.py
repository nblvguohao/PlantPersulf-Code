import csv
import hashlib
import json
import subprocess
from pathlib import Path

import pytest

APPROVED_ACCESSIONS = ("PXD006140", "PXD051570", "GSE163745")
FIXTURE_ROOT = Path("tests/fixtures/real")
EXTRACTION_COMMAND = (
    "python scripts/build_real_fixtures.py --accessions PXD006140,PXD051570,GSE163745"
)
MANIFEST_FIELDS = {
    "schema_version",
    "accession",
    "source_repository",
    "source_registry",
    "source_file",
    "source_sha256",
    "fixture_file",
    "fixture_sha256",
    "extraction_method",
    "extraction_command",
    "biological_values_modified",
    "use",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _standard_manifests() -> dict[str, tuple[Path, dict[str, object]]]:
    """Manifests produced by the canonical fixture builder only.

    Additional scientifically scoped fixtures have their own schemas and
    dedicated parser tests, so they are not subject to this builder contract.
    """
    result: dict[str, tuple[Path, dict[str, object]]] = {}
    for accession in APPROVED_ACCESSIONS:
        path = FIXTURE_ROOT / accession / "source_manifest.json"
        parsed: object = json.loads(path.read_text(encoding="utf-8"))
        assert isinstance(parsed, dict)
        accession = parsed.get("accession")
        assert isinstance(accession, str)
        result[accession] = (path, parsed)
    return result


def _tree_hashes(root: Path) -> dict[str, str]:
    return {
        path.relative_to(root).as_posix(): _sha256(path)
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


def test_real_fixture_bytes_are_exempt_from_git_text_normalization() -> None:
    result = subprocess.run(
        [
            "git",
            "check-attr",
            "text",
            "diff",
            "--",
            "tests/fixtures/real/PXD006140/omssa_head_32.csv",
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    assert "text: unset" in result.stdout
    assert "diff: unset" in result.stdout


def test_real_fixture_manifests_use_canonical_lf_bytes() -> None:
    for manifest_path, _ in _standard_manifests().values():
        manifest_bytes = manifest_path.read_bytes()
        assert manifest_bytes.endswith(b"\n")
        assert b"\r\n" not in manifest_bytes


def test_every_real_fixture_has_complete_manifest() -> None:
    manifests = _standard_manifests()

    assert set(manifests) == set(APPROVED_ACCESSIONS)
    for accession, (manifest_path, manifest) in manifests.items():
        assert set(manifest) == MANIFEST_FIELDS
        assert manifest["schema_version"] == 1
        assert manifest["accession"] == accession
        assert len(str(manifest["source_sha256"])) == 64
        assert len(str(manifest["fixture_sha256"])) == 64
        assert manifest["extraction_command"] == EXTRACTION_COMMAND
        assert manifest["biological_values_modified"] is False
        assert manifest["use"] == "tests_only"
        fixture_path = manifest_path.parent / str(manifest["fixture_file"])
        assert fixture_path.is_file()
        assert sorted(path.name for path in manifest_path.parent.iterdir()) == [
            fixture_path.name,
            "source_manifest.json",
        ]


def test_real_fixtures_preserve_registered_source_bytes() -> None:
    manifests = _standard_manifests()

    for manifest_path, manifest in manifests.values():
        source_path = Path(str(manifest["source_file"]))
        fixture_path = manifest_path.parent / str(manifest["fixture_file"])
        assert _sha256(source_path) == manifest["source_sha256"]
        assert _sha256(fixture_path) == manifest["fixture_sha256"]

    pxd_manifest_path, pxd_manifest = manifests["PXD006140"]
    pxd_source = Path(str(pxd_manifest["source_file"]))
    with pxd_source.open("rb") as handle:
        expected_prefix = b"".join(handle.readline() for _ in range(33))
    pxd_fixture = pxd_manifest_path.parent / str(pxd_manifest["fixture_file"])
    assert pxd_fixture.read_bytes() == expected_prefix
    assert len(pxd_fixture.read_bytes().splitlines()) == 33

    for accession in ("PXD051570", "GSE163745"):
        manifest_path, manifest = manifests[accession]
        fixture_path = manifest_path.parent / str(manifest["fixture_file"])
        source_path = Path(str(manifest["source_file"]))
        assert fixture_path.read_bytes() == source_path.read_bytes()


def test_real_fixture_build_is_deterministic(tmp_path: Path) -> None:
    from scripts.build_real_fixtures import build_real_fixtures

    output_root = tmp_path / "real"
    build_real_fixtures(APPROVED_ACCESSIONS, output_root=output_root)
    first = _tree_hashes(output_root)
    build_real_fixtures(APPROVED_ACCESSIONS, output_root=output_root)
    second = _tree_hashes(output_root)

    assert first == second
    assert len(first) == 6


def test_unregistered_fixture_source_is_rejected(tmp_path: Path) -> None:
    from scripts.build_real_fixtures import FixtureSource, build_fixture

    source_path = tmp_path / "policy.bin"
    source_path.write_bytes(b"software-policy-marker")
    registry_path = tmp_path / "files.tsv"
    with registry_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("path", "sha256"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow({"path": source_path.name, "sha256": "0" * 64})
    source = FixtureSource(
        accession="POLICYTEST",
        source_repository="POLICY",
        source_registry=registry_path,
        source_file=source_path,
        expected_source_sha256=_sha256(source_path),
        fixture_name="policy.bin",
        extraction_method="byte_for_byte_copy",
    )
    output_root = tmp_path / "fixtures"

    with pytest.raises(RuntimeError, match="SHA256 mismatch"):
        build_fixture(source, output_root=output_root)

    assert not (output_root / "POLICYTEST/policy.bin").exists()
    assert not (output_root / "POLICYTEST/source_manifest.json").exists()


@pytest.mark.parametrize(
    "accessions",
    [
        ("PXD006140", "PXD006140"),
        ("UNREGISTERED",),
    ],
)
def test_duplicate_or_unsupported_accessions_are_rejected(
    tmp_path: Path,
    accessions: tuple[str, ...],
) -> None:
    from scripts.build_real_fixtures import build_real_fixtures

    with pytest.raises(RuntimeError, match="duplicate|unsupported"):
        build_real_fixtures(accessions, output_root=tmp_path / "fixtures")


def test_unsupported_existing_manifest_schema_is_rejected(
    tmp_path: Path,
) -> None:
    from scripts.build_real_fixtures import SOURCES, build_fixture

    output_root = tmp_path / "fixtures"
    output_dir = output_root / "GSE163745"
    output_dir.mkdir(parents=True)
    manifest_path = output_dir / "source_manifest.json"
    manifest_path.write_text('{"schema_version": 2}\n', encoding="utf-8")

    with pytest.raises(RuntimeError, match="unsupported existing fixture manifest"):
        build_fixture(SOURCES["GSE163745"], output_root=output_root)

    assert not (output_dir / "GSE163745.metadata.json").exists()
    assert manifest_path.read_text(encoding="utf-8") == '{"schema_version": 2}\n'
