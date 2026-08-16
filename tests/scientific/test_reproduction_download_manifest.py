import csv
from pathlib import Path, PurePosixPath

DOWNLOAD_REGISTRY = Path("data/registry/downloads.tsv")
METHOD_REGISTRY = Path("data/registry/evidence_methods.tsv")
REPRODUCTION_MANIFEST = Path("data/registry/reproduction_downloads.tsv")


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_reproduction_manifest_exactly_covers_registered_inputs() -> None:
    downloads = _read_tsv(DOWNLOAD_REGISTRY)
    methods = _read_tsv(METHOD_REGISTRY)
    manifest = _read_tsv(REPRODUCTION_MANIFEST)

    assert len(manifest) == len(downloads) + len(methods) == 12

    by_key = {(row["study_accession"], row["file_name"]): row for row in manifest}
    for source in downloads:
        row = by_key[(source["dataset_accession"], source["file_name"])]
        assert row["source_type"] == "repository_file"
        assert row["official_url"] == source["download_url"]
        assert row["registry_size_bytes"] == source["registry_size_bytes"]
        assert row["download_size_bytes"] == source["size_bytes"]
        assert row["remote_checksum_algorithm"] == source["remote_checksum_algorithm"]
        assert row["remote_checksum"] == source["remote_checksum"]
        assert row["sha256"] == source["sha256"]

    for source in methods:
        file_name = PurePosixPath(source["local_path"].replace("\\", "/")).name
        row = by_key[(source["study_accession"], file_name)]
        assert row["source_type"] == "method_publication"
        assert row["official_url"] == source["official_url"]
        assert row["registry_size_bytes"] == source["size_bytes"]
        assert row["download_size_bytes"] == source["size_bytes"]
        assert row["remote_checksum_algorithm"] == ""
        assert row["remote_checksum"] == ""
        assert row["sha256"] == source["sha256"]


def test_reproduction_manifest_paths_are_safe_and_relative() -> None:
    for row in _read_tsv(REPRODUCTION_MANIFEST):
        relative_path = PurePosixPath(row["relative_path"])

        assert not relative_path.is_absolute()
        assert ".." not in relative_path.parts
        assert row["official_url"].startswith("https://")
        assert len(row["sha256"]) == 64
