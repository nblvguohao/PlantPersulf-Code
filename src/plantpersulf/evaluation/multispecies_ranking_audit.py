"""Label-free full-proteome ranking audit for a frozen v2 model."""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from pathlib import Path


@dataclass(frozen=True)
class LabelFreeRankedCysteine:
    """A scored reference Cys; labels are intentionally not representable."""

    species: str
    protein_accession: str
    cys_position: int
    score: float


def audit_full_proteome_ranking(
    rows: Iterable[LabelFreeRankedCysteine],
    *,
    expected_cysteines: Iterable[tuple[str, str, int]],
    model_artifact: Path,
    config_path: Path,
    model_freeze_manifest: Path,
    code_revision: str,
    expected_primary_species: tuple[str, ...],
    pressure_species: tuple[str, ...],
) -> dict[str, object]:
    """Validate and separate a label-free ranking after model selection."""
    try:
        freeze = json.loads(model_freeze_manifest.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError("model selection freeze manifest is invalid") from exc
    if not isinstance(freeze, dict) or freeze.get("status") != "frozen":
        raise RuntimeError("model selection must be frozen before ranking audit")
    if not model_artifact.is_file():
        raise RuntimeError(f"model artifact is missing: {model_artifact}")
    if not config_path.is_file():
        raise RuntimeError(f"frozen config is missing: {config_path}")

    model_sha256 = hashlib.sha256(model_artifact.read_bytes()).hexdigest()
    config_sha256 = hashlib.sha256(config_path.read_bytes()).hexdigest()
    if freeze.get("model_sha256") != model_sha256:
        raise RuntimeError("model artifact hash mismatch")
    if freeze.get("config_sha256") != config_sha256:
        raise RuntimeError("frozen config hash mismatch")
    if freeze.get("code_revision") != code_revision:
        raise RuntimeError("frozen code revision mismatch")
    if set(expected_primary_species) & set(pressure_species):
        raise ValueError("primary and pressure species must be disjoint")

    materialized = tuple(rows)
    seen: set[tuple[str, str, int]] = set()
    by_species: dict[str, list[LabelFreeRankedCysteine]] = {}
    allowed_species = set(expected_primary_species) | set(pressure_species)
    for row in materialized:
        key = (row.species, row.protein_accession, row.cys_position)
        if key in seen:
            raise RuntimeError(f"duplicate ranked cysteine coordinate: {key}")
        seen.add(key)
        if row.species not in allowed_species:
            raise RuntimeError(f"unregistered ranking species: {row.species}")
        if row.cys_position < 1:
            raise ValueError("cysteine positions must be one-based")
        if not math.isfinite(row.score):
            raise ValueError("ranking scores must be finite")
        by_species.setdefault(row.species, []).append(row)

    expected_coordinates = set(expected_cysteines)
    if len(expected_coordinates) == 0:
        raise RuntimeError("reference cysteine inventory is empty")
    if seen != expected_coordinates:
        missing_coordinates = expected_coordinates - seen
        extra_coordinates = seen - expected_coordinates
        raise RuntimeError(
            "ranking does not cover reference cysteines exactly: "
            f"missing={len(missing_coordinates)} extra={len(extra_coordinates)}"
        )

    missing = allowed_species - set(by_species)
    if missing:
        raise RuntimeError(f"missing full-proteome ranking species: {sorted(missing)}")

    def ranked(species: str) -> list[dict[str, object]]:
        ordered = sorted(
            by_species[species],
            key=lambda row: (-row.score, row.protein_accession, row.cys_position),
        )
        return [asdict(row) for row in ordered]

    return {
        "three_crop_primary_ranking": {
            species: ranked(species) for species in expected_primary_species
        },
        "magnaporthe_pressure_test_ranking": {
            species: ranked(species) for species in pressure_species
        },
        "output_limitation": "label_free_ranking_audit_not_test_metrics",
        "model_sha256": model_sha256,
        "config_sha256": config_sha256,
        "code_revision": code_revision,
    }
