"""Official NCBI SRA RunInfo client with a traceable JSON cache."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, cast
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

EUTILS_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


@dataclass(frozen=True)
class SraRun:
    run_accession: str
    experiment_accession: str
    biosample_accession: str
    sample_name: str
    bioproject_accession: str
    sra_study_accession: str


@dataclass(frozen=True)
class SraStudy:
    query_accession: str
    bioproject_accession: str
    sra_study_accession: str
    runs: tuple[SraRun, ...]
    source_url: str
    retrieved_at: str
    cache_path: Path
    cache_sha256: str


JsonObject = dict[str, Any]


def _sha256_bytes(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


class SraClient:
    """Resolve an official SRA study or BioProject to its RunInfo rows."""

    def __init__(
        self,
        cache_dir: Path = Path("data/registry/cache/sra"),
        timeout_seconds: float = 30.0,
        retries: int = 3,
    ) -> None:
        self.cache_dir = cache_dir
        self.timeout_seconds = timeout_seconds
        self.retries = retries

    def get_runs(self, accession: str) -> SraStudy:
        normalized = accession.strip().upper()
        if not normalized.startswith(("SRP", "PRJNA")):
            raise ValueError("SRA query accession must start with SRP or PRJNA")
        cache_path = self.cache_dir / f"{normalized}.json"
        if cache_path.exists():
            cache = self._load_cache(cache_path, normalized)
        else:
            cache = self._fetch_cache(normalized, cache_path)
        return self._study_from_cache(cache, cache_path)

    def _request_bytes(self, url: str, accept: str) -> bytes:
        request = Request(
            url,
            headers={"Accept": accept, "User-Agent": "PlantPersulf-Code/0.0.0"},
        )
        for attempt in range(self.retries):
            try:
                with urlopen(request, timeout=self.timeout_seconds) as response:
                    return cast(bytes, response.read())
            except HTTPError as exc:
                retryable = exc.code == 429 or 500 <= exc.code < 600
                if not retryable or attempt == self.retries - 1:
                    raise RuntimeError(
                        f"SRA request failed: {url}: HTTP {exc.code}"
                    ) from exc
            except TimeoutError as exc:
                if attempt == self.retries - 1:
                    raise RuntimeError(f"SRA request failed: {url}: timeout") from exc
            except URLError as exc:
                if attempt == self.retries - 1:
                    raise RuntimeError(
                        f"SRA request failed: {url}: {exc.reason}"
                    ) from exc
            time.sleep(2**attempt)
        raise RuntimeError(f"SRA request exhausted retries: {url}")

    def _fetch_cache(self, accession: str, cache_path: Path) -> JsonObject:
        field = "ACCN" if accession.startswith("SRP") else "BioProject"
        search_query = urlencode(
            {
                "db": "sra",
                "term": f"{accession}[{field}]",
                "retmax": "10000",
                "retmode": "json",
            }
        )
        search_url = f"{EUTILS_BASE_URL}/esearch.fcgi?{search_query}"
        search_content = self._request_bytes(search_url, "application/json")
        try:
            search_payload = json.loads(search_content.decode("utf-8"))
            identifiers = search_payload["esearchresult"]["idlist"]
        except (KeyError, TypeError, json.JSONDecodeError, UnicodeDecodeError) as exc:
            raise RuntimeError(f"invalid SRA search response for {accession}") from exc
        if not identifiers:
            raise RuntimeError(f"SRA search returned no records for {accession}")
        fetch_query = urlencode(
            {
                "db": "sra",
                "id": ",".join(identifiers),
                "rettype": "runinfo",
                "retmode": "text",
            }
        )
        fetch_url = f"{EUTILS_BASE_URL}/efetch.fcgi?{fetch_query}"
        runinfo_content = self._request_bytes(fetch_url, "text/csv")
        cache: JsonObject = {
            "schema_version": 1,
            "repository": "SRA",
            "query_accession": accession,
            "retrieved_at": datetime.now(timezone.utc).isoformat(),
            "search_url": search_url,
            "fetch_url": fetch_url,
            "search_response_sha256": _sha256_bytes(search_content),
            "runinfo_response_sha256": _sha256_bytes(runinfo_content),
            "search_response": search_payload,
            "runinfo_csv": runinfo_content.decode("utf-8-sig"),
        }
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        serialized = json.dumps(cache, ensure_ascii=False, indent=2, sort_keys=True)
        temporary_path = cache_path.with_suffix(".json.tmp")
        temporary_path.write_bytes(f"{serialized}\n".encode())
        temporary_path.replace(cache_path)
        return cache

    def _load_cache(self, cache_path: Path, accession: str) -> JsonObject:
        try:
            parsed: object = json.loads(cache_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise RuntimeError(f"invalid SRA metadata cache: {cache_path}") from exc
        if not isinstance(parsed, dict):
            raise RuntimeError(f"invalid SRA metadata cache: {cache_path}")
        cache = cast(JsonObject, parsed)
        if (
            cache.get("repository") != "SRA"
            or cache.get("query_accession") != accession
        ):
            raise RuntimeError(f"SRA cache identity mismatch: {cache_path}")
        return cache

    def _study_from_cache(self, cache: JsonObject, cache_path: Path) -> SraStudy:
        runinfo_csv = cache.get("runinfo_csv")
        if not isinstance(runinfo_csv, str):
            raise RuntimeError(f"SRA cache lacks RunInfo CSV: {cache_path}")
        rows = tuple(csv.DictReader(io.StringIO(runinfo_csv)))
        if not rows:
            raise RuntimeError(f"SRA cache has no runs: {cache_path}")
        runs = tuple(
            SraRun(
                run_accession=row["Run"],
                experiment_accession=row["Experiment"],
                biosample_accession=row["BioSample"],
                sample_name=row["SampleName"],
                bioproject_accession=row["BioProject"],
                sra_study_accession=row["SRAStudy"],
            )
            for row in rows
        )
        query_accession = str(cache["query_accession"])
        if query_accession.startswith("SRP") and {
            run.sra_study_accession for run in runs
        } != {query_accession}:
            raise RuntimeError(f"SRA study identity mismatch: {cache_path}")
        if query_accession.startswith("PRJNA") and {
            run.bioproject_accession for run in runs
        } != {query_accession}:
            raise RuntimeError(f"SRA BioProject identity mismatch: {cache_path}")
        return SraStudy(
            query_accession=query_accession,
            bioproject_accession=runs[0].bioproject_accession,
            sra_study_accession=runs[0].sra_study_accession,
            runs=tuple(sorted(runs, key=lambda run: run.run_accession)),
            source_url=str(cache["fetch_url"]),
            retrieved_at=str(cache["retrieved_at"]),
            cache_path=cache_path,
            cache_sha256=_sha256_bytes(cache_path.read_bytes()),
        )
