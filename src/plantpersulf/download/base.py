"""Atomic transport for registered scientific files."""

from __future__ import annotations

import tempfile
from dataclasses import dataclass
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from plantpersulf.provenance.hashing import hash_file


@dataclass(frozen=True)
class DownloadRequest:
    source_url: str
    destination: Path
    expected_size: int | None = None
    expected_checksum: str = ""
    expected_checksum_algorithm: str = ""
    timeout_seconds: float = 60.0


@dataclass(frozen=True)
class DownloadResult:
    destination: Path
    size_bytes: int
    sha256: str


def canonical_download_url(source_url: str) -> str:
    """Use PRIDE's official HTTPS transport for its registered FTP paths."""
    pride_ftp_prefix = "ftp://ftp.pride.ebi.ac.uk/pride/"
    if source_url.startswith(pride_ftp_prefix):
        return source_url.replace("ftp://", "https://", 1)
    return source_url


def download_verified_file(request: DownloadRequest) -> DownloadResult:
    """Stream a URL to a temporary sibling and publish it atomically."""
    request.destination.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None
    try:
        source_request = Request(
            canonical_download_url(request.source_url),
            headers={"User-Agent": "PlantPersulf-Code/0.0.0"},
        )
        with urlopen(source_request, timeout=request.timeout_seconds) as response:
            with tempfile.NamedTemporaryFile(
                mode="wb",
                dir=request.destination.parent,
                prefix=f".{request.destination.name}.",
                suffix=".part",
                delete=False,
            ) as handle:
                temporary_path = Path(handle.name)
                size_bytes = 0
                while chunk := response.read(1024 * 1024):
                    handle.write(chunk)
                    size_bytes += len(chunk)
        if request.expected_size is not None and size_bytes != request.expected_size:
            raise RuntimeError(
                "download size mismatch: "
                f"expected {request.expected_size}, observed {size_bytes}"
            )
        if request.expected_checksum:
            if not request.expected_checksum_algorithm:
                raise RuntimeError("download checksum algorithm is required")
            observed_checksum = hash_file(
                temporary_path,
                request.expected_checksum_algorithm,
            )
            if observed_checksum.lower() != request.expected_checksum.lower():
                raise RuntimeError(
                    "download checksum mismatch: "
                    f"expected {request.expected_checksum}, "
                    f"observed {observed_checksum}"
                )
        sha256 = hash_file(temporary_path, "sha256")
        temporary_path.replace(request.destination)
        temporary_path = None
        return DownloadResult(request.destination, size_bytes, sha256)
    except (HTTPError, URLError, TimeoutError, OSError) as exc:
        raise RuntimeError(f"download failed: {request.source_url}") from exc
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
