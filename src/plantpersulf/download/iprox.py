"""ProteomeXchange metadata client for datasets hosted by iProX."""

from __future__ import annotations

import hashlib
import json
import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, cast
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

PROTEOMEXCHANGE_DATASET_URL = (
    "https://proteomecentral.proteomexchange.org/cgi/GetDataset"
)


@dataclass(frozen=True)
class IproxFile:
    file_name: str
    category: str
    source_url: str


@dataclass(frozen=True)
class IproxReference:
    pubmed_id: str
    citation: str


@dataclass(frozen=True)
class IproxDataset:
    accession: str
    title: str
    publication_date: str
    organisms: tuple[str, ...]
    files: tuple[IproxFile, ...]
    references: tuple[IproxReference, ...]
    source_url: str
    retrieved_at: str
    cache_path: Path
    cache_sha256: str


JsonObject = dict[str, Any]


def _sha256_bytes(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def _metadata_url(accession: str) -> str:
    query = urlencode(
        {
            "ID": f"{accession}.0-1",
            "outputMode": "XML",
            "test": "no",
        }
    )
    return f"{PROTEOMEXCHANGE_DATASET_URL}?{query}"


def _parse_xml(source_xml: str, accession: str) -> ET.Element:
    try:
        root = ET.fromstring(source_xml)
    except ET.ParseError as exc:
        raise RuntimeError("ProteomeXchange returned invalid XML") from exc
    if root.tag != "ProteomeXchangeDataset" or root.get("id") != accession:
        raise RuntimeError("ProteomeXchange XML accession mismatch")
    summary = root.find("DatasetSummary")
    if summary is None or summary.get("hostingRepository") != "iProX":
        raise RuntimeError(
            f"ProteomeXchange dataset is not hosted by iProX: {accession}"
        )
    return root


class IproxClient:
    """Fetch ProteomeXchange announcement metadata for an iProX dataset."""

    def __init__(
        self,
        cache_dir: Path = Path("data/registry/cache/iprox"),
        timeout_seconds: float = 30.0,
        retries: int = 3,
    ) -> None:
        self.cache_dir = cache_dir
        self.timeout_seconds = timeout_seconds
        self.retries = retries

    def get_dataset(self, accession: str) -> IproxDataset:
        normalized_accession = accession.strip().upper()
        if not normalized_accession.startswith("PXD"):
            raise ValueError("ProteomeXchange accession must start with PXD")
        cache_path = self.cache_dir / f"{normalized_accession}.json"
        if cache_path.exists():
            cache = self._load_cache(cache_path, normalized_accession)
        else:
            cache = self._fetch_cache(normalized_accession, cache_path)
        return self._dataset_from_cache(cache, cache_path)

    def _request_xml(self, url: str) -> tuple[str, bytes]:
        request = Request(
            url,
            headers={
                "Accept": "application/xml",
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
                        f"ProteomeXchange request failed: {url}: HTTP {exc.code}"
                    ) from exc
            except TimeoutError as exc:
                if attempt == self.retries - 1:
                    raise RuntimeError(
                        f"ProteomeXchange request failed: {url}: timeout"
                    ) from exc
            except URLError as exc:
                if attempt == self.retries - 1:
                    raise RuntimeError(
                        f"ProteomeXchange request failed: {url}: {exc.reason}"
                    ) from exc
            time.sleep(2**attempt)
        raise RuntimeError(f"ProteomeXchange request exhausted retries: {url}")

    def _fetch_cache(self, accession: str, cache_path: Path) -> JsonObject:
        source_url = _metadata_url(accession)
        source_xml, source_content = self._request_xml(source_url)
        _parse_xml(source_xml, accession)
        cache: JsonObject = {
            "schema_version": 1,
            "repository": "iProX",
            "accession": accession,
            "retrieved_at": datetime.now(timezone.utc).isoformat(),
            "source_url": source_url,
            "source_response_sha256": _sha256_bytes(source_content),
            "source_xml": source_xml,
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
            raise RuntimeError(f"invalid iProX metadata cache: {cache_path}") from exc
        if not isinstance(parsed, dict):
            raise RuntimeError(f"iProX metadata cache must be a mapping: {cache_path}")
        cache = cast(JsonObject, parsed)
        if cache.get("repository") != "iProX" or cache.get("accession") != accession:
            raise RuntimeError(f"iProX cache identity mismatch: {cache_path}")
        return cache

    def _dataset_from_cache(
        self,
        cache: JsonObject,
        cache_path: Path,
    ) -> IproxDataset:
        accession = str(cache.get("accession", ""))
        root = _parse_xml(str(cache.get("source_xml", "")), accession)
        summary = root.find("DatasetSummary")
        if summary is None:
            raise RuntimeError("ProteomeXchange XML is missing DatasetSummary")
        organisms = tuple(
            parameter.get("value", "")
            for parameter in root.findall("./SpeciesList/Species/cvParam")
            if parameter.get("name") == "taxonomy: scientific name"
        )
        files = tuple(
            self._file_from_element(element)
            for element in root.findall("./DatasetFileList/DatasetFile")
        )
        references = tuple(
            self._reference_from_element(element)
            for element in root.findall("./PublicationList/Publication")
        )
        return IproxDataset(
            accession=accession,
            title=summary.get("title", ""),
            publication_date=summary.get("announceDate", ""),
            organisms=organisms,
            files=files,
            references=references,
            source_url=str(cache.get("source_url", "")),
            retrieved_at=str(cache.get("retrieved_at", "")),
            cache_path=cache_path,
            cache_sha256=_sha256_bytes(cache_path.read_bytes()),
        )

    @staticmethod
    def _file_from_element(element: ET.Element) -> IproxFile:
        parameter = element.find("cvParam")
        if parameter is None:
            raise RuntimeError("ProteomeXchange DatasetFile lacks a URI")
        return IproxFile(
            file_name=element.get("name", ""),
            category=parameter.get("name", ""),
            source_url=parameter.get("value", ""),
        )

    @staticmethod
    def _reference_from_element(element: ET.Element) -> IproxReference:
        pubmed_id = ""
        citation = ""
        for parameter in element.findall("cvParam"):
            if parameter.get("name") == "PubMed identifier":
                pubmed_id = parameter.get("value", "")
            elif parameter.get("name") == "Reference":
                citation = parameter.get("value", "")
        return IproxReference(pubmed_id=pubmed_id, citation=citation)
