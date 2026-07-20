import hashlib
from pathlib import Path

import pytest


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_unregistered_file_is_rejected(tmp_path: Path) -> None:
    from plantpersulf.provenance.audit import assert_registered_input

    unknown = tmp_path / "unknown.tsv"
    unknown.write_text("column\tvalue\nentry\t1\n", encoding="utf-8")

    with pytest.raises(RuntimeError, match="not registered"):
        assert_registered_input(unknown, registry_path=tmp_path / "files.tsv")


def test_file_absent_from_registry_is_rejected(tmp_path: Path) -> None:
    from plantpersulf.provenance.audit import assert_registered_input

    unknown = tmp_path / "unknown.tsv"
    unknown.write_text("column\tvalue\nentry\t1\n", encoding="utf-8")
    registry = tmp_path / "files.tsv"
    registry.write_text(
        "path\tsha256\nother.tsv\t" + "a" * 64 + "\n",
        encoding="utf-8",
    )

    with pytest.raises(RuntimeError, match="not registered"):
        assert_registered_input(unknown, registry_path=registry)


def test_registered_file_with_wrong_sha256_is_rejected(tmp_path: Path) -> None:
    from plantpersulf.provenance.audit import assert_registered_input

    source = tmp_path / "source.tsv"
    source.write_text("column\tvalue\nentry\t1\n", encoding="utf-8")
    registry = tmp_path / "files.tsv"
    registry.write_text(
        "path\tsha256\nsource.tsv\t" + "0" * 64 + "\n",
        encoding="utf-8",
    )

    with pytest.raises(RuntimeError, match="SHA256 mismatch"):
        assert_registered_input(source, registry_path=registry)


def test_registered_file_with_matching_sha256_is_accepted(tmp_path: Path) -> None:
    from plantpersulf.provenance.audit import assert_registered_input

    source = tmp_path / "source.tsv"
    source.write_text("column\tvalue\nentry\t1\n", encoding="utf-8")
    registry = tmp_path / "files.tsv"
    registry.write_text(
        f"path\tsha256\nsource.tsv\t{_sha256(source)}\n",
        encoding="utf-8",
    )

    assert_registered_input(source, registry_path=registry)
