"""Official NCBI GEO SOFT metadata client with a traceable JSON cache."""

from __future__ import annotations

import hashlib
import json
import re
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, cast
from urllib.error import HTTPError, URLError
from urllib.parse import parse_qs, urlencode, urlparse
from urllib.request import Request, urlopen

GEO_ACCESSION_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"


@dataclass(frozen=True)
class GeoSample:
    accession: str
    title: str
    organism: str
    source_name: str
    characteristics: tuple[str, ...]
    biosample_accession: str
    sra_experiment_accession: str


@dataclass(frozen=True)
class GeoSeries:
    accession: str
    title: str
    organism: str
    publication_date: str
    samples: tuple[GeoSample, ...]
    pubmed_ids: tuple[str, ...]
    supplementary_files: tuple[str, ...]
    bioproject_accession: str
    sra_study_accession: str
    source_url: str
    samples_source_url: str
    retrieved_at: str
    cache_path: Path
    cache_sha256: str


@dataclass(frozen=True)
class _SoftEntity:
    entity_type: str
    accession: str
    attributes: dict[str, tuple[str, ...]]


JsonObject = dict[str, Any]


def _sha256_bytes(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def _soft_url(accession: str, target: str) -> str:
    query = urlencode(
        {
            "acc": accession,
            "targ": target,
            "form": "text",
            "view": "brief",
        }
    )
    return f"{GEO_ACCESSION_URL}?{query}"


def _parse_soft(content: str) -> tuple[_SoftEntity, ...]:
    entities: list[_SoftEntity] = []
    entity_type = ""
    accession = ""
    attributes: dict[str, list[str]] = {}

    def finish_entity() -> None:
        if not entity_type:
            return
        frozen_attributes = {key: tuple(values) for key, values in attributes.items()}
        entities.append(_SoftEntity(entity_type, accession, frozen_attributes))

    for raw_line in content.splitlines():
        line = raw_line.strip()
        if line.startswith("^") and " = " in line:
            finish_entity()
            label, accession_value = line[1:].split(" = ", 1)
            entity_type = label.upper()
            accession = accession_value.strip()
            attributes = {}
        elif line.startswith("!") and " = " in line and entity_type:
            label, value = line[1:].split(" = ", 1)
            attributes.setdefault(label, []).append(value.strip())
    finish_entity()
    return tuple(entities)


def _first(attributes: dict[str, tuple[str, ...]], key: str) -> str:
    values = attributes.get(key, ())
    return values[0] if values else ""


def _relation_accession(relations: tuple[str, ...], relation_type: str) -> str:
    prefix = f"{relation_type}:"
    for relation in relations:
        if not relation.startswith(prefix):
            continue
        url = relation.split(":", 1)[1].strip()
        parsed = urlparse(url)
        term = parse_qs(parsed.query).get("term", [])
        if term:
            return term[0]
        return parsed.path.rstrip("/").rsplit("/", 1)[-1]
    return ""


def _publication_date(status: str) -> str:
    match = re.fullmatch(r"Public on (.+)", status)
    if match is None:
        raise RuntimeError(f"GEO series has unsupported publication status: {status}")
    try:
        return datetime.strptime(match.group(1), "%b %d %Y").date().isoformat()
    except ValueError as exc:
        raise RuntimeError(
            f"GEO series has invalid publication date: {status}"
        ) from exc


class GeoClient:
    """Fetch a GEO Series and its submitter-supplied Sample metadata."""

    def __init__(
        self,
        cache_dir: Path = Path("data/registry/cache/geo"),
        timeout_seconds: float = 30.0,
        retries: int = 3,
    ) -> None:
        self.cache_dir = cache_dir
        self.timeout_seconds = timeout_seconds
        self.retries = retries

    def get_series(self, accession: str) -> GeoSeries:
        normalized_accession = accession.strip().upper()
        if not normalized_accession.startswith("GSE"):
            raise ValueError("GEO Series accession must start with GSE")
        cache_path = self.cache_dir / f"{normalized_accession}.json"
        if cache_path.exists():
            cache = self._load_cache(cache_path, normalized_accession)
        else:
            cache = self._fetch_cache(normalized_accession, cache_path)
        return self._series_from_cache(cache, cache_path)

    def _request_text(self, url: str) -> tuple[str, bytes]:
        request = Request(
            url,
            headers={
                "Accept": "text/plain",
                "User-Agent": "PlantPersulf-Code/0.0.0",
            },
        )
        for attempt in range(self.retries):
            try:
                with urlopen(request, timeout=self.timeout_seconds) as response:
                    content = response.read()
                return content.decode("utf-8"), content
            except HTTPError as exc:
                retryable = exc.code == 429 or 500 <= exc.code < 600
                if not retryable or attempt == self.retries - 1:
                    raise RuntimeError(
                        f"GEO request failed: {url}: HTTP {exc.code}"
                    ) from exc
            except TimeoutError as exc:
                if attempt == self.retries - 1:
                    raise RuntimeError(f"GEO request failed: {url}: timeout") from exc
            except URLError as exc:
                if attempt == self.retries - 1:
                    raise RuntimeError(
                        f"GEO request failed: {url}: {exc.reason}"
                    ) from exc
            time.sleep(2**attempt)
        raise RuntimeError(f"GEO request exhausted retries: {url}")

    def _fetch_cache(self, accession: str, cache_path: Path) -> JsonObject:
        series_url = _soft_url(accession, "self")
        samples_url = _soft_url(accession, "gsm")
        series_soft, series_content = self._request_text(series_url)
        samples_soft, samples_content = self._request_text(samples_url)
        cache: JsonObject = {
            "schema_version": 1,
            "repository": "GEO",
            "accession": accession,
            "retrieved_at": datetime.now(timezone.utc).isoformat(),
            "series_url": series_url,
            "samples_url": samples_url,
            "series_response_sha256": _sha256_bytes(series_content),
            "samples_response_sha256": _sha256_bytes(samples_content),
            "series_soft": series_soft,
            "samples_soft": samples_soft,
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
            raise RuntimeError(f"invalid GEO metadata cache: {cache_path}") from exc
        if not isinstance(parsed, dict):
            raise RuntimeError(f"GEO metadata cache must be a mapping: {cache_path}")
        cache = cast(JsonObject, parsed)
        if cache.get("repository") != "GEO" or cache.get("accession") != accession:
            raise RuntimeError(f"GEO cache identity mismatch: {cache_path}")
        return cache

    def _series_from_cache(self, cache: JsonObject, cache_path: Path) -> GeoSeries:
        series_entities = _parse_soft(str(cache.get("series_soft", "")))
        sample_entities = _parse_soft(str(cache.get("samples_soft", "")))
        series_matches = [
            entity for entity in series_entities if entity.entity_type == "SERIES"
        ]
        if len(series_matches) != 1:
            raise RuntimeError("GEO response must contain exactly one Series")
        series = series_matches[0]
        if series.accession != cache.get("accession"):
            raise RuntimeError("GEO response accession does not match cache")

        expected_samples = set(series.attributes.get("Series_sample_id", ()))
        samples = tuple(
            self._sample_from_entity(entity)
            for entity in sample_entities
            if entity.entity_type == "SAMPLE"
        )
        if {sample.accession for sample in samples} != expected_samples:
            raise RuntimeError("GEO Sample response does not match Series sample list")

        relations = series.attributes.get("Series_relation", ())
        return GeoSeries(
            accession=series.accession,
            title=_first(series.attributes, "Series_title"),
            organism=_first(series.attributes, "Series_sample_organism"),
            publication_date=_publication_date(
                _first(series.attributes, "Series_status")
            ),
            samples=samples,
            pubmed_ids=series.attributes.get("Series_pubmed_id", ()),
            supplementary_files=series.attributes.get("Series_supplementary_file", ()),
            bioproject_accession=_relation_accession(relations, "BioProject"),
            sra_study_accession=_relation_accession(relations, "SRA"),
            source_url=str(cache.get("series_url", "")),
            samples_source_url=str(cache.get("samples_url", "")),
            retrieved_at=str(cache.get("retrieved_at", "")),
            cache_path=cache_path,
            cache_sha256=_sha256_bytes(cache_path.read_bytes()),
        )

    @staticmethod
    def _sample_from_entity(entity: _SoftEntity) -> GeoSample:
        relations = entity.attributes.get("Sample_relation", ())
        return GeoSample(
            accession=entity.accession,
            title=_first(entity.attributes, "Sample_title"),
            organism=_first(entity.attributes, "Sample_organism_ch1"),
            source_name=_first(entity.attributes, "Sample_source_name_ch1"),
            characteristics=entity.attributes.get("Sample_characteristics_ch1", ()),
            biosample_accession=_relation_accession(relations, "BioSample"),
            sra_experiment_accession=_relation_accession(relations, "SRA"),
        )
