"""Fail-closed AlphaFold DB structure download and provenance registry.

AlphaFold DB predicted structures are served at a stable per-accession URL:
``https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v4.pdb``. Not
every UniProt accession has a model: UniProt isoform accessions (those
containing ``-``, e.g. ``P27140-2``) are never hosted by AlphaFold DB, and
some canonical accessions legitimately return HTTP 404. Both cases are real
"no predicted structure" facts, not download failures — batch callers must
record them as an explicit missing-structure result and never raise, and
never substitute a placeholder or fabricated file. Any other transport or
integrity failure (timeout, 5xx, malformed response) still fails closed by
raising ``RuntimeError``, exactly like ``plantpersulf.download.base``.
"""

from __future__ import annotations

import csv
import json as _json
import os
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request, urlopen

import yaml

from plantpersulf.download.base import DownloadRequest, download_verified_file
from plantpersulf.provenance.hashing import hash_file

ALPHAFOLD_API_URL = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"
ALPHAFOLD_MODEL_VERSION = "v4"  # fallback; overridden by API response when available
DOWNLOADER_VERSION = "plantpersulf/0.0.0"

STRUCTURE_FIELDS = (
    "accession",
    "model_version",
    "source_url",
    "local_path",
    "retrieved_at",
    "size_bytes",
    "sha256",
    "downloader_version",
)

_STATUS_DOWNLOADED = "downloaded"
_STATUS_ISOFORM = "isoform_not_applicable"
_STATUS_NOT_FOUND = "not_found"


@dataclass(frozen=True)
class AlphaFoldFetchResult:
    accession: str
    status: str  # "downloaded" | "isoform_not_applicable" | "not_found"
    source_url: str
    destination: Path | None
    size_bytes: int | None
    sha256: str | None
    retrieved_at: str | None
    model_version: str | None = None


@dataclass(frozen=True)
class AlphaFoldStructureSource:
    accession: str
    model_version: str
    source_url: str
    local_path: Path
    retrieved_at: str
    size_bytes: int
    sha256: str
    downloader_version: str


@dataclass(frozen=True)
class AlphaFoldMetadata:
    accession: str
    model_entity_id: str
    latest_version: int
    pdb_url: str
    cif_url: str
    global_metric_value: float | None
    fraction_plddt_very_high: float | None
    fraction_plddt_low: float | None


def query_alphafold_api(
    accession: str,
    timeout_seconds: float = 30.0,
) -> AlphaFoldMetadata | None:
    """Query the AlphaFold DB prediction API for per-accession metadata.

    Returns ``None`` when the accession has no prediction (genuine absence);
    raises only on network or parse errors (fail-closed: the caller must
    decide whether a transient failure becomes a retry, a clean
    not_found, or a raised RuntimeError).
    """
    url = ALPHAFOLD_API_URL.format(accession=accession)
    req = Request(url, headers={"Accept": "application/json"})
    with urlopen(req, timeout=timeout_seconds) as response:
        data = _json.load(response)
    if not data or not isinstance(data, list):
        return None
    entry = data[0]
    return AlphaFoldMetadata(
        accession=accession,
        model_entity_id=entry.get("modelEntityId", f"AF-{accession}-F1"),
        latest_version=int(entry.get("latestVersion", 4)),
        pdb_url=str(entry.get("pdbUrl", "")),
        cif_url=str(entry.get("cifUrl", "")),
        global_metric_value=(
            float(entry["globalMetricValue"]) if "globalMetricValue" in entry else None
        ),
        fraction_plddt_very_high=(
            float(entry["fractionPlddtVeryHigh"])
            if "fractionPlddtVeryHigh" in entry
            else None
        ),
        fraction_plddt_low=(
            float(entry["fractionPlddtLow"]) if "fractionPlddtLow" in entry else None
        ),
    )


def alphafold_pdb_url(accession: str) -> str:
    """Convenience: return the API URL (not the final PDB URL — callers that
    need the actual download URL should use ``query_alphafold_api`` or pass
    the source_url through ``fetch_alphafold_structure``)."""
    return ALPHAFOLD_API_URL.format(accession=accession)


def is_isoform_accession(accession: str) -> bool:
    """UniProt isoform accessions (e.g. ``P27140-2``) are not in AlphaFold DB."""
    return "-" in accession


def fetch_alphafold_structure(
    accession: str,
    destination: Path,
    source_url: str | None = None,
    timeout_seconds: float = 60.0,
) -> AlphaFoldFetchResult:
    """Download one AlphaFold model, or return an explicit missing result.

    Isoform accessions and HTTP 404 responses are real "no structure" facts
    and are returned as a clean, non-raising result. Any other transport or
    integrity failure raises ``RuntimeError`` (fail closed) and never leaves
    a partial or placeholder file behind.
    """
    # Use the provided URL, or discover the latest version via API.
    # Do NOT pre-filter on "-" (isoform syntax): AlphaFold DB does host some
    # UniProt isoform accessions (e.g. P27140-2 → AF-P27140-2-F1).
    if source_url is not None:
        resolved_url = source_url
        model_version = "custom"
    else:
        try:
            meta = query_alphafold_api(accession, timeout_seconds=timeout_seconds)
        except Exception as exc:
            raise RuntimeError(f"AlphaFold API query failed: {accession}") from exc
        if meta is None:
            return AlphaFoldFetchResult(
                accession=accession,
                status=_STATUS_NOT_FOUND,
                source_url=ALPHAFOLD_API_URL.format(accession=accession),
                destination=None,
                size_bytes=None,
                sha256=None,
                retrieved_at=None,
            )
        resolved_url = meta.pdb_url
        model_version = f"v{meta.latest_version}"

    try:
        result = download_verified_file(
            DownloadRequest(
                source_url=resolved_url,
                destination=destination,
                timeout_seconds=timeout_seconds,
            )
        )
    except RuntimeError as exc:
        cause = exc.__cause__
        if isinstance(cause, HTTPError) and cause.code == 404:
            return AlphaFoldFetchResult(
                accession=accession,
                status=_STATUS_NOT_FOUND,
                source_url=resolved_url,
                destination=None,
                size_bytes=None,
                sha256=None,
                retrieved_at=None,
            )
        raise RuntimeError(f"AlphaFold download failed: {resolved_url}") from exc
    return AlphaFoldFetchResult(
        accession=accession,
        status=_STATUS_DOWNLOADED,
        source_url=resolved_url,
        destination=result.destination,
        size_bytes=result.size_bytes,
        sha256=result.sha256,
        retrieved_at=datetime.now(timezone.utc).isoformat(),
        model_version=model_version,
    )


def _read_registry(registry_path: Path) -> list[dict[str, str]]:
    if not registry_path.exists():
        return []
    with registry_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != STRUCTURE_FIELDS:
            raise RuntimeError("AlphaFold structure registry has invalid columns")
        return [dict(row) for row in reader]


def _write_registry(registry_path: Path, rows: list[dict[str, str]]) -> None:
    registry_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = registry_path.with_suffix(".tsv.tmp")
    with temporary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=STRUCTURE_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    temporary_path.replace(registry_path)


def register_alphafold_structure(
    result: AlphaFoldFetchResult,
    registry_path: Path,
) -> None:
    """Append a successfully downloaded structure with full provenance.

    Fails closed: only a ``status == "downloaded"`` result with a complete
    size/sha256/retrieved_at triple may be registered.
    """
    if (
        result.status != _STATUS_DOWNLOADED
        or result.destination is None
        or result.size_bytes is None
        or result.sha256 is None
        or result.retrieved_at is None
    ):
        raise RuntimeError(
            f"only a downloaded AlphaFold structure can be registered: "
            f"{result.accession}"
        )
    rows = _read_registry(registry_path)
    rows = [row for row in rows if row["accession"] != result.accession]
    rows.append(
        {
            "accession": result.accession,
            "model_version": result.model_version or ALPHAFOLD_MODEL_VERSION,
            "source_url": result.source_url,
            "local_path": Path(
                os.path.relpath(result.destination, registry_path.parent)
            ).as_posix(),
            "retrieved_at": result.retrieved_at,
            "size_bytes": str(result.size_bytes),
            "sha256": result.sha256,
            "downloader_version": DOWNLOADER_VERSION,
        }
    )
    rows.sort(key=lambda row: row["accession"])
    _write_registry(registry_path, rows)


def audit_alphafold_structures(
    registry_path: Path = Path("data/registry/alphafold_structures.tsv"),
    *,
    base_directory: Path | None = None,
) -> tuple[AlphaFoldStructureSource, ...]:
    """Re-verify every registered structure's SHA256 against local bytes.

    ``local_path`` is recorded relative to the directory holding the registry
    (see ``register_alphafold_structure``), so that is the default resolution
    root. A frozen *snapshot* of a registry may legitimately be stored
    somewhere else (e.g. ``data/registry/releases/``) while still describing
    the same downloaded files; such a caller must state the original registry
    directory via ``base_directory``. It is passed explicitly and never
    guessed: silently searching parent directories for a structure file would
    be exactly the kind of quiet repair this auditor exists to prevent.
    """
    root = registry_path.parent if base_directory is None else base_directory
    rows = _read_registry(registry_path)
    sources: list[AlphaFoldStructureSource] = []
    for row in rows:
        local_path = (root / row["local_path"]).resolve()
        try:
            size = int(row["size_bytes"])
        except ValueError as exc:
            raise RuntimeError(
                f"AlphaFold structure has invalid size: {row['accession']}"
            ) from exc
        if (
            not row["retrieved_at"]
            or re.fullmatch(r"[0-9a-f]{64}", row["sha256"]) is None
            or not local_path.is_file()
            or local_path.stat().st_size != size
            or hash_file(local_path, "sha256") != row["sha256"]
        ):
            raise RuntimeError(
                f"AlphaFold structure failed local audit: {row['accession']}"
            )
        sources.append(
            AlphaFoldStructureSource(
                accession=row["accession"],
                model_version=row["model_version"],
                source_url=row["source_url"],
                local_path=local_path,
                retrieved_at=row["retrieved_at"],
                size_bytes=size,
                sha256=row["sha256"],
                downloader_version=row["downloader_version"],
            )
        )
    return tuple(sources)


def load_approved_accessions(
    config_path: Path = Path("configs/alphafold_sources_v1.yaml"),
) -> tuple[str, ...]:
    """Load the allow-listed candidate UniProt accessions for this cycle."""
    try:
        loaded: object = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise RuntimeError(f"invalid AlphaFold source config: {config_path}") from exc
    if not isinstance(loaded, dict) or loaded.get("version") != 1:
        raise RuntimeError("AlphaFold source config requires version 1")
    accessions = loaded.get("accessions")
    if not isinstance(accessions, list) or not accessions:
        raise RuntimeError("AlphaFold source config requires accessions")
    if not all(isinstance(item, str) and item.strip() for item in accessions):
        raise RuntimeError("AlphaFold source config has an invalid accession")
    if len(set(accessions)) != len(accessions):
        raise RuntimeError("AlphaFold source config has a duplicate accession")
    return tuple(accessions)
