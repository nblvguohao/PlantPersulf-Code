"""Versioned scientific policy for evidence preflight."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

import yaml

EXPECTED_EVIDENCE_CLASSES = (
    "site_ms",
    "site_mutagenesis",
    "site_biochemical",
    "protein_level_only",
    "identification_only",
    "unresolved",
)
EXPECTED_PROHIBITED_TERMS = ("positive", "unlabeled", "negative")
LARGE_FILE_THRESHOLD_BYTES = 104_857_600


@dataclass(frozen=True)
class PreflightPolicy:
    version: int
    studies: tuple[str, ...]
    evidence_classes: tuple[str, ...]
    large_file_threshold_bytes: int
    prohibited_label_terms: tuple[str, ...]


def _string_tuple(value: object, context: str) -> tuple[str, ...]:
    if not isinstance(value, list) or any(
        not isinstance(item, str) or not item for item in value
    ):
        raise RuntimeError(f"evidence preflight policy has invalid {context}")
    return tuple(value)


def load_preflight_policy(path: Path) -> PreflightPolicy:
    """Load the exact reviewed evidence-preflight contract."""
    try:
        loaded: object = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise RuntimeError(f"invalid evidence preflight policy: {path}") from exc
    if not isinstance(loaded, dict):
        raise RuntimeError("evidence preflight policy must be a mapping")
    policy = cast(dict[str, Any], loaded)
    if policy.get("version") != 1:
        raise RuntimeError("evidence preflight policy requires version 1")
    studies = _string_tuple(policy.get("studies"), "studies")
    if len(set(studies)) != len(studies) or any(
        not accession.startswith("PXD") for accession in studies
    ):
        raise RuntimeError("evidence preflight policy has invalid studies")
    evidence_classes = _string_tuple(
        policy.get("evidence_classes"),
        "evidence classes",
    )
    if evidence_classes != EXPECTED_EVIDENCE_CLASSES:
        raise RuntimeError("evidence preflight policy has invalid evidence classes")
    threshold = policy.get("large_file_threshold_bytes")
    if threshold != LARGE_FILE_THRESHOLD_BYTES:
        raise RuntimeError("evidence preflight policy has invalid threshold")
    prohibited = _string_tuple(
        policy.get("prohibited_label_terms"),
        "prohibited label terms",
    )
    if prohibited != EXPECTED_PROHIBITED_TERMS:
        raise RuntimeError("evidence preflight policy has invalid prohibited terms")
    return PreflightPolicy(
        version=1,
        studies=studies,
        evidence_classes=evidence_classes,
        large_file_threshold_bytes=LARGE_FILE_THRESHOLD_BYTES,
        prohibited_label_terms=prohibited,
    )
