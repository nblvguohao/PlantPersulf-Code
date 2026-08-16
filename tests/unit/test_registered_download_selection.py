import csv
import hashlib
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from threading import Thread

import pytest


class _RegisteredPayloadHandler(BaseHTTPRequestHandler):
    payload = b"registered-download-policy-marker"

    def do_GET(self) -> None:  # noqa: N802
        self.send_response(200)
        self.send_header("Content-Length", str(len(self.payload)))
        self.end_headers()
        self.wfile.write(self.payload)

    def log_message(self, format: str, *args: object) -> None:
        return


class _PartialFailureHandler(BaseHTTPRequestHandler):
    payload = b"first-registered-policy-marker"

    def do_GET(self) -> None:  # noqa: N802
        if self.path.endswith("/second.txt"):
            self.send_response(503)
            self.end_headers()
            return
        self.send_response(200)
        self.send_header("Content-Length", str(len(self.payload)))
        self.end_headers()
        self.wfile.write(self.payload)

    def log_message(self, format: str, *args: object) -> None:
        return


def test_pxd006140_selection_is_exact_and_explicit() -> None:
    from plantpersulf.download.registered import resolve_registered_selection

    selected = resolve_registered_selection(
        accession="PXD006140",
        file_classes=("metadata", "results"),
        selection_path=Path("configs/download_selection.yaml"),
        files_registry_path=Path("data/registry/files.tsv"),
    )

    assert [(item.file_class, item.file_name) for item in selected] == [
        ("metadata", "PXD006140.json"),
        ("results", "omssa.20150821_01_AAroca_TMT6plex.cmpd.mgf.txt"),
        ("results", "omssa.ne.20150821_01_AAroca_TMT6plex.cmpd.mgf.txt"),
    ]


def test_unapproved_file_class_is_fatal() -> None:
    from plantpersulf.download.registered import resolve_registered_selection

    with pytest.raises(RuntimeError, match="not approved"):
        resolve_registered_selection(
            accession="PXD006140",
            file_classes=("raw",),
            selection_path=Path("configs/download_selection.yaml"),
            files_registry_path=Path("data/registry/files.tsv"),
        )


def test_registered_download_writes_complete_provenance_manifest(
    tmp_path: Path,
) -> None:
    from plantpersulf.download.registered import (
        audit_downloaded_files,
        download_registered_files,
    )

    server = ThreadingHTTPServer(("127.0.0.1", 0), _RegisteredPayloadHandler)
    thread = Thread(target=server.serve_forever, daemon=True)
    thread.start()
    source_url = f"http://127.0.0.1:{server.server_port}/policy.txt"
    registry_dir = tmp_path / "data/registry"
    registry_dir.mkdir(parents=True)
    selection_path = tmp_path / "configs/selection.yaml"
    selection_path.parent.mkdir()
    selection_path.write_text(
        "version: 1\napproved:\n  POLICYTEST:\n    metadata:\n"
        "      - policy.json\n    results:\n"
        "      - policy.txt\n",
        encoding="utf-8",
    )
    metadata_path = registry_dir / "cache/policy.json"
    metadata_path.parent.mkdir()
    metadata_path.write_bytes(b"official-metadata-policy-marker")
    metadata_sha256 = hashlib.sha256(metadata_path.read_bytes()).hexdigest()
    files_path = registry_dir / "files.tsv"
    files_path.write_text(
        "dataset_accession\trepository\trecord_type\tfile_name\tfile_category"
        "\tsource_url\tsize_bytes\tremote_checksum"
        "\tremote_checksum_algorithm\tpath\tsha256\tstatus\tretrieved_at\n"
        "POLICYTEST\tPOLICY\tmetadata_cache\tpolicy.json\tofficial_metadata\t"
        f"https://example.invalid/policy\t{metadata_path.stat().st_size}\t\t\t"
        f"cache/policy.json\t{metadata_sha256}\tcached\t"
        "2026-07-20T00:00:00+00:00\n"
        f"POLICYTEST\tPOLICY\tsource_file\tpolicy.txt\tPOLICY\t{source_url}\t"
        f"{len(_RegisteredPayloadHandler.payload) + 1}\t"
        f"{hashlib.sha1(_RegisteredPayloadHandler.payload).hexdigest()}\tSHA1\t"
        "\t\tremote_only\t2026-07-20T00:00:00+00:00\n",
        encoding="utf-8",
    )
    datasets_path = registry_dir / "datasets.tsv"
    datasets_path.write_text(
        "accession\trepository\tlicense_or_usage\n"
        "POLICYTEST\tPOLICY\tSoftware policy test only\n",
        encoding="utf-8",
    )
    manifest_path = registry_dir / "downloads.tsv"
    summary = download_registered_files(
        accession="POLICYTEST",
        file_classes=("metadata", "results"),
        selection_path=selection_path,
        files_registry_path=files_path,
        datasets_registry_path=datasets_path,
        downloads_registry_path=manifest_path,
        raw_dir=tmp_path / "data/raw",
    )

    assert summary.selected_count == 2
    assert summary.downloaded_count == 1
    assert summary.cached_count == 1
    with manifest_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    row = rows[0]
    assert row["dataset_accession"] == "POLICYTEST"
    assert row["source_url"] == source_url
    assert row["download_url"] == source_url
    assert row["remote_checksum_algorithm"] == "SHA1"
    assert row["registry_size_bytes"] == str(len(_RegisteredPayloadHandler.payload) + 1)
    assert row["size_bytes"] == str(len(_RegisteredPayloadHandler.payload))
    assert len(row["sha256"]) == 64
    assert row["license_or_usage"] == "Software policy test only"
    assert row["downloader_version"] == "plantpersulf/0.0.0"
    assert (manifest_path.parent / row["path"]).read_bytes() == (
        _RegisteredPayloadHandler.payload
    )
    audit_summary = audit_downloaded_files(
        accession="POLICYTEST",
        selection_path=selection_path,
        files_registry_path=files_path,
        datasets_registry_path=datasets_path,
        downloads_registry_path=manifest_path,
    )
    assert audit_summary.selected_count == 2
    assert audit_summary.downloaded_count == 1
    assert audit_summary.cached_count == 1
    reused_summary = download_registered_files(
        accession="POLICYTEST",
        file_classes=("metadata", "results"),
        selection_path=selection_path,
        files_registry_path=files_path,
        datasets_registry_path=datasets_path,
        downloads_registry_path=manifest_path,
        raw_dir=tmp_path / "data/raw",
    )
    assert reused_summary.selected_count == 2
    assert reused_summary.downloaded_count == 0
    assert reused_summary.cached_count == 2
    local_path = manifest_path.parent / row["path"]
    local_path.write_bytes(b"x" * len(_RegisteredPayloadHandler.payload))
    with pytest.raises(RuntimeError, match="SHA256 mismatch"):
        audit_downloaded_files(
            accession="POLICYTEST",
            selection_path=selection_path,
            files_registry_path=files_path,
            datasets_registry_path=datasets_path,
            downloads_registry_path=manifest_path,
        )
    recovered_summary = download_registered_files(
        accession="POLICYTEST",
        file_classes=("metadata", "results"),
        selection_path=selection_path,
        files_registry_path=files_path,
        datasets_registry_path=datasets_path,
        downloads_registry_path=manifest_path,
        raw_dir=tmp_path / "data/raw",
    )
    server.shutdown()
    server.server_close()
    thread.join()

    assert recovered_summary.selected_count == 2
    assert recovered_summary.downloaded_count == 1
    assert recovered_summary.cached_count == 1
    repaired_summary = audit_downloaded_files(
        accession="POLICYTEST",
        selection_path=selection_path,
        files_registry_path=files_path,
        datasets_registry_path=datasets_path,
        downloads_registry_path=manifest_path,
    )
    assert repaired_summary.downloaded_count == 1


def test_each_successful_file_is_registered_before_a_later_failure(
    tmp_path: Path,
) -> None:
    from plantpersulf.download.registered import download_registered_files

    server = ThreadingHTTPServer(("127.0.0.1", 0), _PartialFailureHandler)
    thread = Thread(target=server.serve_forever, daemon=True)
    thread.start()
    registry_dir = tmp_path / "data/registry"
    registry_dir.mkdir(parents=True)
    selection_path = tmp_path / "selection.yaml"
    selection_path.write_text(
        "version: 1\napproved:\n  POLICYTEST:\n    results:\n"
        "      - first.txt\n      - second.txt\n",
        encoding="utf-8",
    )
    checksum = hashlib.sha1(_PartialFailureHandler.payload).hexdigest()
    source_prefix = f"http://127.0.0.1:{server.server_port}"
    files_path = registry_dir / "files.tsv"
    header = (
        "dataset_accession\trepository\trecord_type\tfile_name\tfile_category"
        "\tsource_url\tsize_bytes\tremote_checksum"
        "\tremote_checksum_algorithm\tpath\tsha256\tstatus\tretrieved_at\n"
    )
    rows = "".join(
        "POLICYTEST\tPOLICY\tsource_file\t"
        f"{name}\tPOLICY\t{source_prefix}/{name}\t"
        f"{len(_PartialFailureHandler.payload)}\t{checksum}\tSHA1\t\t\t"
        "remote_only\t2026-07-20T00:00:00+00:00\n"
        for name in ("first.txt", "second.txt")
    )
    files_path.write_text(header + rows, encoding="utf-8")
    datasets_path = registry_dir / "datasets.tsv"
    datasets_path.write_text(
        "accession\trepository\tlicense_or_usage\n"
        "POLICYTEST\tPOLICY\tSoftware policy test only\n",
        encoding="utf-8",
    )
    manifest_path = registry_dir / "downloads.tsv"
    try:
        with pytest.raises(RuntimeError, match="download failed"):
            download_registered_files(
                accession="POLICYTEST",
                file_classes=("results",),
                selection_path=selection_path,
                files_registry_path=files_path,
                datasets_registry_path=datasets_path,
                downloads_registry_path=manifest_path,
                raw_dir=tmp_path / "data/raw",
            )
    finally:
        server.shutdown()
        server.server_close()
        thread.join()

    with manifest_path.open(encoding="utf-8", newline="") as handle:
        manifest_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["file_name"] for row in manifest_rows] == ["first.txt"]
    assert (tmp_path / "data/raw/POLICYTEST/first.txt").is_file()
    assert not (tmp_path / "data/raw/POLICYTEST/second.txt").exists()
