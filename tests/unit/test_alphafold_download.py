"""RED (Task 7): fail-closed AlphaFold DB structure download + registry.

AlphaFold DB predicted structures live at a stable per-accession URL. Not
every UniProt accession has a model: isoform accessions (containing "-")
are never hosted, and some canonical accessions 404. Both are real
"no predicted structure" facts and must come back as an explicit, clean
result — never an exception, and never a substituted/fabricated file. Any
other transport failure (5xx, timeout) must still fail closed by raising.

Expected RED: ``plantpersulf.download.alphafold`` does not exist yet.
"""

from __future__ import annotations

from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from threading import Thread

import pytest

from plantpersulf.download.alphafold import (  # RED: module missing
    STRUCTURE_FIELDS,
    AlphaFoldFetchResult,
    alphafold_pdb_url,
    audit_alphafold_structures,
    fetch_alphafold_structure,
    is_isoform_accession,
    load_approved_accessions,
    register_alphafold_structure,
)

CONFIG = Path("configs/alphafold_sources_v1.yaml")


class _NotFoundHandler(BaseHTTPRequestHandler):
    def do_GET(self) -> None:  # noqa: N802
        self.send_response(404)
        self.end_headers()

    def log_message(self, format: str, *args: object) -> None:
        return


class _ServerErrorHandler(BaseHTTPRequestHandler):
    def do_GET(self) -> None:  # noqa: N802
        self.send_response(503)
        self.end_headers()

    def log_message(self, format: str, *args: object) -> None:
        return


class _PdbPayloadHandler(BaseHTTPRequestHandler):
    # Not a real structure: a short marker payload for exercising code paths
    # only (mirrors tests/unit/test_download_failure_is_fatal.py style).
    payload = b"alphafold-download-policy-marker"

    def do_GET(self) -> None:  # noqa: N802
        self.send_response(200)
        self.send_header("Content-Length", str(len(self.payload)))
        self.end_headers()
        self.wfile.write(self.payload)

    def log_message(self, format: str, *args: object) -> None:
        return


def _run_server(handler: type) -> ThreadingHTTPServer:
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = Thread(target=server.serve_forever, daemon=True)
    thread.start()
    return server


def _stop_server(server: ThreadingHTTPServer) -> None:
    server.shutdown()
    server.server_close()


def test_alphafold_pdb_url_uses_official_v4_pattern() -> None:
    assert alphafold_pdb_url("O03042") == (
        "https://alphafold.ebi.ac.uk/files/AF-O03042-F1-model_v4.pdb"
    )


def test_isoform_accession_is_detected() -> None:
    assert is_isoform_accession("P27140-2") is True
    assert is_isoform_accession("O03042") is False


def test_isoform_accession_returns_clean_result_without_network(
    tmp_path: Path,
) -> None:
    result = fetch_alphafold_structure(
        "P27140-2",
        tmp_path / "P27140-2.pdb",
    )

    assert result.status == "isoform_not_applicable"
    assert result.destination is None
    assert not (tmp_path / "P27140-2.pdb").exists()


def test_404_response_is_explicit_missing_result(tmp_path: Path) -> None:
    server = _run_server(_NotFoundHandler)
    try:
        result = fetch_alphafold_structure(
            "O03042",
            tmp_path / "O03042.pdb",
            source_url=f"http://127.0.0.1:{server.server_port}/missing.pdb",
        )
    finally:
        _stop_server(server)

    assert result.status == "not_found"
    assert result.destination is None
    assert not (tmp_path / "O03042.pdb").exists()


def test_server_error_fails_closed_and_leaves_no_file(tmp_path: Path) -> None:
    server = _run_server(_ServerErrorHandler)
    try:
        with pytest.raises(RuntimeError, match="AlphaFold download failed"):
            fetch_alphafold_structure(
                "O03042",
                tmp_path / "O03042.pdb",
                source_url=f"http://127.0.0.1:{server.server_port}/error.pdb",
            )
    finally:
        _stop_server(server)

    assert not (tmp_path / "O03042.pdb").exists()
    assert list(tmp_path.iterdir()) == []


def test_successful_download_registers_full_provenance(tmp_path: Path) -> None:
    server = _run_server(_PdbPayloadHandler)
    try:
        destination = tmp_path / "raw" / "O03042.pdb"
        result = fetch_alphafold_structure(
            "O03042",
            destination,
            source_url=f"http://127.0.0.1:{server.server_port}/AF-O03042.pdb",
        )
    finally:
        _stop_server(server)

    assert result.status == "downloaded"
    assert isinstance(result, AlphaFoldFetchResult)
    assert destination.is_file()
    assert result.size_bytes == len(_PdbPayloadHandler.payload)

    registry_path = tmp_path / "registry" / "alphafold_structures.tsv"
    register_alphafold_structure(result, registry_path)

    audited = audit_alphafold_structures(registry_path)
    assert len(audited) == 1
    assert audited[0].accession == "O03042"
    assert audited[0].sha256 == result.sha256

    # header round-trips through the declared schema
    header = registry_path.read_text(encoding="utf-8").splitlines()[0].split("\t")
    assert tuple(header) == STRUCTURE_FIELDS


def test_only_a_downloaded_result_can_be_registered(tmp_path: Path) -> None:
    missing = fetch_alphafold_structure("P27140-2", tmp_path / "x.pdb")

    with pytest.raises(RuntimeError, match="only a downloaded"):
        register_alphafold_structure(missing, tmp_path / "registry.tsv")


def test_tampered_registry_sha256_is_rejected(tmp_path: Path) -> None:
    server = _run_server(_PdbPayloadHandler)
    try:
        destination = tmp_path / "O03042.pdb"
        result = fetch_alphafold_structure(
            "O03042",
            destination,
            source_url=f"http://127.0.0.1:{server.server_port}/AF-O03042.pdb",
        )
    finally:
        _stop_server(server)

    registry_path = tmp_path / "alphafold_structures.tsv"
    register_alphafold_structure(result, registry_path)
    tampered = registry_path.read_text(encoding="utf-8").replace(
        result.sha256 or "", "0" * 64
    )
    registry_path.write_text(tampered, encoding="utf-8")

    with pytest.raises(RuntimeError):
        audit_alphafold_structures(registry_path)


def test_approved_accessions_config_lists_the_pilot_batch() -> None:
    accessions = load_approved_accessions(CONFIG)

    assert accessions == (
        "O03042",
        "Q9SLA0",
        "A0A1P8B1I9",
        "A0A1P8APX1",
        "Q8VZF1",
        "Q944G9",
        "P27140-2",
    )
    assert is_isoform_accession(accessions[-1]) is True
