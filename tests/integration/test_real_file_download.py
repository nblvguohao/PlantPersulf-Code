import csv
from pathlib import Path

import pytest


@pytest.mark.network
def test_registered_pride_sdrf_download_matches_size_and_checksum(
    tmp_path: Path,
) -> None:
    from plantpersulf.download.base import DownloadRequest, download_verified_file
    from plantpersulf.provenance.hashing import hash_file

    with Path("data/registry/files.tsv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    row = next(
        row
        for row in rows
        if row["dataset_accession"] == "PXD035795" and row["file_name"] == "SDRF.txt"
    )
    destination = tmp_path / row["file_name"]

    result = download_verified_file(
        DownloadRequest(
            source_url=row["source_url"],
            destination=destination,
            expected_size=int(row["size_bytes"]),
            expected_checksum=row["remote_checksum"],
            expected_checksum_algorithm=row["remote_checksum_algorithm"],
        )
    )

    assert result.size_bytes == 4671
    assert result.sha256 == hash_file(destination, "sha256")
    assert hash_file(destination, "sha1") == row["remote_checksum"]
