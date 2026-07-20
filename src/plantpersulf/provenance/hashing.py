"""Streaming cryptographic hashes for provenance records."""

from __future__ import annotations

import hashlib
from pathlib import Path


def hash_file(path: Path, algorithm: str) -> str:
    """Return a lowercase SHA1 or SHA256 digest for a local file."""
    normalized = algorithm.strip().lower()
    if normalized not in {"sha1", "sha256"}:
        raise ValueError(f"unsupported checksum algorithm: {algorithm}")
    digest = hashlib.new(normalized)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
