from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from threading import Thread

import pytest


class _UnavailableHandler(BaseHTTPRequestHandler):
    def do_GET(self) -> None:  # noqa: N802
        self.send_response(503)
        self.end_headers()

    def log_message(self, format: str, *args: object) -> None:
        return


class _PayloadHandler(BaseHTTPRequestHandler):
    payload = b"software-policy-marker"

    def do_GET(self) -> None:  # noqa: N802
        self.send_response(200)
        self.send_header("Content-Length", str(len(self.payload)))
        self.end_headers()
        self.wfile.write(self.payload)

    def log_message(self, format: str, *args: object) -> None:
        return


def test_pride_ftp_url_uses_same_host_and_path_over_https() -> None:
    from plantpersulf.download.base import canonical_download_url

    source_url = (
        "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2023/02/PXD035795/SDRF.txt"
    )

    assert canonical_download_url(source_url) == source_url.replace(
        "ftp://", "https://", 1
    )


def test_http_failure_is_fatal_and_leaves_no_file(tmp_path: Path) -> None:
    from plantpersulf.download.base import DownloadRequest, download_verified_file

    server = ThreadingHTTPServer(("127.0.0.1", 0), _UnavailableHandler)
    thread = Thread(target=server.serve_forever, daemon=True)
    thread.start()
    destination = tmp_path / "result.txt"
    try:
        with pytest.raises(RuntimeError, match="download failed"):
            download_verified_file(
                DownloadRequest(
                    source_url=f"http://127.0.0.1:{server.server_port}/result.txt",
                    destination=destination,
                )
            )
    finally:
        server.shutdown()
        server.server_close()
        thread.join()

    assert not destination.exists()
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize(
    ("expected_size", "expected_checksum", "expected_algorithm", "message"),
    [
        (1, "", "", "size mismatch"),
        (
            len(_PayloadHandler.payload),
            "0" * 40,
            "SHA1",
            "checksum mismatch",
        ),
    ],
)
def test_integrity_failure_is_fatal_and_leaves_no_file(
    tmp_path: Path,
    expected_size: int,
    expected_checksum: str,
    expected_algorithm: str,
    message: str,
) -> None:
    from plantpersulf.download.base import DownloadRequest, download_verified_file

    server = ThreadingHTTPServer(("127.0.0.1", 0), _PayloadHandler)
    thread = Thread(target=server.serve_forever, daemon=True)
    thread.start()
    destination = tmp_path / "result.txt"
    try:
        with pytest.raises(RuntimeError, match=message):
            download_verified_file(
                DownloadRequest(
                    source_url=f"http://127.0.0.1:{server.server_port}/result.txt",
                    destination=destination,
                    expected_size=expected_size,
                    expected_checksum=expected_checksum,
                    expected_checksum_algorithm=expected_algorithm,
                )
            )
    finally:
        server.shutdown()
        server.server_close()
        thread.join()

    assert not destination.exists()
    assert list(tmp_path.iterdir()) == []
