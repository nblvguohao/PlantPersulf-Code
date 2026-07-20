"""Official PRIDE Archive metadata client with provenance-preserving cache."""

from __future__ import annotations

import hashlib
import json
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, cast
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

PRIDE_API_BASE = "https://www.ebi.ac.uk/pride/ws/archive/v3"


@dataclass(frozen=True)
class PrideFile:
    file_name: str
    category: str
    size_bytes: int
    checksum: str
    public_urls: tuple[str, ...]


@dataclass(frozen=True)
class PrideReference:
    pubmed_id: str
    doi: str
    citation: str


@dataclass(frozen=True)
class PrideProject:
    accession: str
    title: str
    publication_date: str
    organisms: tuple[str, ...]
    files: tuple[PrideFile, ...]
    references: tuple[PrideReference, ...]
    license_name: str
    source_url: str
    retrieved_at: str
    cache_path: Path
    cache_sha256: str


JsonObject = dict[str, Any]
JsonValue = JsonObject | list[Any]


def _sha256_bytes(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def _mapping(value: object, context: str) -> JsonObject:
    if not isinstance(value, dict):
        raise RuntimeError(f"PRIDE response requires mapping at {context}")
    return cast(JsonObject, value)


def _items(value: object, context: str) -> list[Any]:
    if not isinstance(value, list):
        raise RuntimeError(f"PRIDE response requires list at {context}")
    return value


class PrideClient:
    """Fetch PRIDE project metadata and its official file listing."""

    def __init__(
        self,
        cache_dir: Path = Path("data/registry/cache/pride"),
        timeout_seconds: float = 30.0,
        retries: int = 3,
    ) -> None:
        self.cache_dir = cache_dir
        self.timeout_seconds = timeout_seconds
        self.retries = retries

    def get_project(self, accession: str) -> PrideProject:
        normalized_accession = accession.strip().upper()
        if not normalized_accession.startswith("PXD"):
            raise ValueError("PRIDE accession must start with PXD")
        cache_path = self.cache_dir / f"{normalized_accession}.json"
        if cache_path.exists():
            cache = self._load_cache(cache_path, normalized_accession)
        else:
            cache = self._fetch_cache(normalized_accession, cache_path)
        return self._project_from_cache(cache, cache_path)

    def _request_json(self, url: str) -> tuple[JsonValue, bytes]:
        request = Request(
            url,
            headers={
                "Accept": "application/json",
                "User-Agent": "PlantPersulf-Code/0.0.0",
            },
        )
        for attempt in range(self.retries):
            try:
                with urlopen(request, timeout=self.timeout_seconds) as response:
                    content = response.read()
                parsed: object = json.loads(content.decode("utf-8"))
                if not isinstance(parsed, (dict, list)):
                    raise RuntimeError("PRIDE returned non-container JSON")
                return cast(JsonValue, parsed), content
            except HTTPError as exc:
                retryable = exc.code == 429 or 500 <= exc.code < 600
                if not retryable or attempt == self.retries - 1:
                    raise RuntimeError(
                        f"PRIDE request failed: {url}: HTTP {exc.code}"
                    ) from exc
            except TimeoutError as exc:
                if attempt == self.retries - 1:
                    raise RuntimeError(
                        f"PRIDE request failed: {url}: timeout"
                    ) from exc
            except URLError as exc:
                if attempt == self.retries - 1:
                    raise RuntimeError(
                        f"PRIDE request failed: {url}: {exc.reason}"
                    ) from exc
            time.sleep(2**attempt)
        raise RuntimeError(f"PRIDE request exhausted retries: {url}")

    def _fetch_cache(self, accession: str, cache_path: Path) -> JsonObject:
        project_url = f"{PRIDE_API_BASE}/projects/{accession}"
        files_url = f"{project_url}/files"
        project_payload, project_content = self._request_json(project_url)
        files_payload, files_content = self._request_json(files_url)
        cache: JsonObject = {
            "schema_version": 1,
            "repository": "PRIDE",
            "accession": accession,
            "retrieved_at": datetime.now(timezone.utc).isoformat(),
            "project_url": project_url,
            "files_url": files_url,
            "project_response_sha256": _sha256_bytes(project_content),
            "files_response_sha256": _sha256_bytes(files_content),
            "project": project_payload,
            "files": files_payload,
        }
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        serialized = json.dumps(cache, ensure_ascii=False, indent=2, sort_keys=True)
        encoded = f"{serialized}\n".encode()
        temporary_path = cache_path.with_suffix(".json.tmp")
        temporary_path.write_bytes(encoded)
        temporary_path.replace(cache_path)
        return cache

    def _load_cache(self, cache_path: Path, accession: str) -> JsonObject:
        try:
            parsed: object = json.loads(cache_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise RuntimeError(f"invalid PRIDE metadata cache: {cache_path}") from exc
        cache = _mapping(parsed, "cache")
        if cache.get("repository") != "PRIDE" or cache.get("accession") != accession:
            raise RuntimeError(f"PRIDE cache identity mismatch: {cache_path}")
        return cache

    def _project_from_cache(self, cache: JsonObject, cache_path: Path) -> PrideProject:
        project = _mapping(cache.get("project"), "project")
        files = tuple(
            self._file_from_payload(_mapping(item, "files[]"))
            for item in _items(cache.get("files"), "files")
        )
        organisms = tuple(
            str(_mapping(item, "organisms[]").get("name", ""))
            for item in _items(project.get("organisms", []), "organisms")
        )
        references = tuple(
            self._reference_from_payload(_mapping(item, "references[]"))
            for item in _items(project.get("references", []), "references")
        )
        return PrideProject(
            accession=str(project.get("accession", "")),
            title=str(project.get("title", "")),
            publication_date=str(project.get("publicationDate", "")),
            organisms=organisms,
            files=files,
            references=references,
            license_name=str(project.get("license", "")),
            source_url=str(cache.get("project_url", "")),
            retrieved_at=str(cache.get("retrieved_at", "")),
            cache_path=cache_path,
            cache_sha256=_sha256_bytes(cache_path.read_bytes()),
        )

    @staticmethod
    def _file_from_payload(payload: JsonObject) -> PrideFile:
        category = _mapping(payload.get("fileCategory", {}), "fileCategory")
        locations = _items(
            payload.get("publicFileLocations", []),
            "publicFileLocations",
        )
        return PrideFile(
            file_name=str(payload.get("fileName", "")),
            category=str(category.get("value", "")),
            size_bytes=int(payload.get("fileSizeBytes", 0)),
            checksum=str(payload.get("checksum", "")),
            public_urls=tuple(
                str(_mapping(item, "publicFileLocations[]").get("value", ""))
                for item in locations
            ),
        )

    @staticmethod
    def _reference_from_payload(payload: JsonObject) -> PrideReference:
        return PrideReference(
            pubmed_id=str(payload.get("pubmedID", "")),
            doi=str(payload.get("doi", "")),
            citation=str(payload.get("referenceLine", "")),
        )
