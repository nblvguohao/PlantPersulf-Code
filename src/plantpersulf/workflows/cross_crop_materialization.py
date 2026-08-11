"""Deterministic, target-label-free input materialization for Task 9A.

This module builds source positive--unlabeled (PU) rows only from registered
public inputs. It deliberately has no dependency on tomato experimental
labels: every cysteine in the frozen tomato reference is a scoring candidate.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from plantpersulf.evaluation.cross_species_conservation import (
    load_panther_annotations,
    load_panther_annotations_by_orf_gene,
    mg8_accession_to_gene,
)
from plantpersulf.features.sequence import _load_proteome
from plantpersulf.features.site_biology import (
    BIOLOGY_FEATURE_NAMES,
    build_site_biology_vector_from_sequence,
)
from plantpersulf.proteomics.pxd063170_sites import (
    load_ensembl_fungi_proteome,
    parse_pxd063170_sites,
)
from plantpersulf.proteomics.pxd072089_sites import parse_pxd072089_sites
from plantpersulf.provenance.audit import assert_registered_input
from plantpersulf.provenance.hashing import hash_file


@dataclass(frozen=True, order=True)
class MaterializationSite:
    """A coordinate-verified source site; ``unlabeled`` is never a negative."""

    protein_accession: str
    cys_position: int
    label: str

    @property
    def site_key(self) -> str:
        return f"{self.protein_accession}:C{self.cys_position}"


@dataclass(frozen=True)
class RegisteredMaterializationInput:
    input_id: str
    source_accession: str
    source_kind: str
    path: Path
    sha256: str
    scientific_use: str


def _sample_rank(site: MaterializationSite, seed: int) -> str:
    payload = f"{seed}|{site.site_key}".encode()
    return hashlib.sha256(payload).hexdigest()


def deterministic_pu_source_sample(
    rows: tuple[MaterializationSite, ...],
    unlabeled_per_positive: int,
    seed: int,
) -> tuple[MaterializationSite, ...]:
    """Keep all positives and a hash-selected real unlabeled background.

    Selection is independent of input order. It changes neither biological
    coordinates nor labels, and retains positive--unlabeled rather than
    artificial positive--negative semantics.
    """
    if unlabeled_per_positive <= 0:
        raise ValueError("unlabeled_per_positive must be positive")
    if set(row.label for row in rows) - {"positive", "unlabeled"}:
        raise ValueError("source rows must use positive or unlabeled labels")
    if len({row.site_key for row in rows}) != len(rows):
        raise ValueError("source rows contain duplicate site coordinates")
    positives = tuple(sorted(row for row in rows if row.label == "positive"))
    if not positives:
        raise ValueError("source rows must contain at least one positive")
    unlabeled = sorted(
        (row for row in rows if row.label == "unlabeled"),
        key=lambda row: (_sample_rank(row, seed), row.site_key),
    )
    selected = unlabeled[: unlabeled_per_positive * len(positives)]
    return tuple(sorted((*positives, *selected)))


def validate_label_free_target_payload(payload: dict[str, Any]) -> None:
    """Fail closed if an input request attempts to provide tomato labels."""
    if "labels" in payload:
        raise RuntimeError("target label leakage: target payload contains labels")


def _canonical_sha256(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def build_source_payload(
    *,
    species: str,
    study_accession: str,
    rows: tuple[MaterializationSite, ...],
    proteome: dict[str, str],
    source_input_sha256s: tuple[str, ...],
    unlabeled_per_positive: int,
    seed: int,
    panther_subfamilies: dict[str, frozenset[str]],
) -> dict[str, Any]:
    """Build a runner-compatible source batch from verified site coordinates.

    PANTHER subfamilies are used only when they cover every sampled protein;
    otherwise the runner's declared fallback is protein-grouped validation.
    """
    if not species or not study_accession or not source_input_sha256s:
        raise ValueError("source identity and registered input hashes are required")
    selected = deterministic_pu_source_sample(rows, unlabeled_per_positive, seed)
    protein_ids = [row.protein_accession for row in selected]
    features = []
    for row in selected:
        sequence = proteome.get(row.protein_accession)
        if sequence is None:
            raise RuntimeError(
                f"source site protein absent from proteome: {row.site_key}"
            )
        vector = build_site_biology_vector_from_sequence(
            row.protein_accession, row.cys_position, sequence
        )
        if vector.names != BIOLOGY_FEATURE_NAMES:
            raise RuntimeError("source feature contract mismatch")
        features.append(list(vector.values))
    payload: dict[str, Any] = {
        "species": species,
        "study_accession": study_accession,
        "feature_names": list(BIOLOGY_FEATURE_NAMES),
        "features": features,
        "labels": [row.label for row in selected],
        "protein_ids": protein_ids,
    }
    if all(panther_subfamilies.get(protein) for protein in protein_ids):
        payload["cluster_ids"] = [
            "panther:" + sorted(panther_subfamilies[protein])[0]
            for protein in protein_ids
        ]
    payload["source_sha256"] = _canonical_sha256(
        {
            "source_input_sha256s": sorted(source_input_sha256s),
            "sampling": {
                "unlabeled_per_positive": unlabeled_per_positive,
                "deterministic_hash_seed": seed,
            },
            "batch": payload,
        }
    )
    return payload


def build_target_payload(
    proteome: dict[str, str], *, source_sha256: str
) -> dict[str, Any]:
    """Enumerate every reference-proteome cysteine as a label-free target row."""
    if not source_sha256:
        raise ValueError("target reference SHA256 is required")
    site_keys: list[str] = []
    features: list[list[float]] = []
    for accession in sorted(proteome):
        sequence = proteome[accession]
        for position, residue in enumerate(sequence, start=1):
            if residue != "C":
                continue
            site_keys.append(f"{accession}:C{position}")
            vector = build_site_biology_vector_from_sequence(
                accession, position, sequence
            )
            if vector.names != BIOLOGY_FEATURE_NAMES:
                raise RuntimeError("target feature contract mismatch")
            features.append(list(vector.values))
    payload: dict[str, Any] = {
        "site_keys": site_keys,
        "feature_names": list(BIOLOGY_FEATURE_NAMES),
        "features": features,
        "source_sha256": source_sha256,
    }
    validate_label_free_target_payload(payload)
    if not site_keys:
        raise RuntimeError("target reference contains no cysteine candidates")
    return payload


def _load_registered_inputs(
    registry_path: Path,
) -> dict[str, RegisteredMaterializationInput]:
    with registry_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "input_id",
            "source_accession",
            "source_kind",
            "path",
            "sha256",
            "scientific_use",
        }
        if set(reader.fieldnames or ()) != required:
            raise RuntimeError("materialization input registry schema mismatch")
        records = [
            RegisteredMaterializationInput(
                input_id=row["input_id"],
                source_accession=row["source_accession"],
                source_kind=row["source_kind"],
                path=(registry_path.parent / row["path"]).resolve(),
                sha256=row["sha256"],
                scientific_use=row["scientific_use"],
            )
            for row in reader
        ]
    if not records or len({item.input_id for item in records}) != len(records):
        raise RuntimeError("materialization input registry has duplicate or no IDs")
    inputs = {item.input_id: item for item in records}
    for item in inputs.values():
        assert_registered_input(item.path, registry_path)
    return inputs


def _source_rows_from_proteome(
    proteome: dict[str, str],
    positive_keys: set[tuple[str, int]],
    supplemental_proteome: dict[str, str] | None = None,
) -> tuple[MaterializationSite, ...]:
    if not positive_keys:
        raise RuntimeError("source study has no coordinate-verified positives")
    rows = []
    seen_positive: set[tuple[str, int]] = set()
    for accession in sorted(proteome):
        sequence = proteome[accession]
        for position, residue in enumerate(sequence, start=1):
            if residue != "C":
                continue
            key = (accession, position)
            if key in positive_keys:
                seen_positive.add(key)
            rows.append(
                MaterializationSite(
                    accession,
                    position,
                    "positive" if key in positive_keys else "unlabeled",
                )
            )
    supplemental = supplemental_proteome or {}
    missing = positive_keys - seen_positive
    for accession, position in sorted(missing):
        supplemental_sequence = supplemental.get(accession)
        if (
            supplemental_sequence is None
            or position < 1
            or position > len(supplemental_sequence)
            or supplemental_sequence[position - 1] != "C"
        ):
            continue
        rows.append(MaterializationSite(accession, position, "positive"))
        seen_positive.add((accession, position))
    missing = positive_keys - seen_positive
    if missing:
        raise RuntimeError(
            "positive coordinates are absent from the reference proteome"
        )
    return tuple(rows)


def _benchmark_positive_keys(
    benchmark_path: Path, study_accession: str
) -> set[tuple[str, int]]:
    with benchmark_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected = (
            "protein_accession",
            "cys_position_in_protein",
            "label",
            "study_accession",
            "evidence_level",
            "source_sha256",
        )
        if tuple(reader.fieldnames or ()) != expected:
            raise RuntimeError("Arabidopsis benchmark schema mismatch")
        positives = {
            (row["protein_accession"], int(row["cys_position_in_protein"]))
            for row in reader
            if row["label"] == "positive"
            and row["study_accession"] == study_accession
        }
    if not positives:
        raise RuntimeError(
            f"benchmark has no positives for source study: {study_accession}"
        )
    return positives


def _checked_config(config_path: Path) -> dict[str, Any]:
    loaded = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if not isinstance(loaded, dict):
        raise RuntimeError("materialization config must be a mapping")
    sampling = loaded.get("source_unlabeled_sampling")
    target = loaded.get("target")
    sources = loaded.get("sources")
    if sampling != {
        "label_semantics": "positive_unlabeled_only",
        "unlabeled_per_positive": 20,
        "selection": "deterministic_sha256_site_key_order",
        "deterministic_hash_seed": 20260811,
    }:
        raise RuntimeError("materialization sampling policy differs from frozen config")
    if not isinstance(target, dict) or target.get("labels_visible") is not False:
        raise RuntimeError(
            "target labels must remain unavailable during materialization"
        )
    if (
        target.get("candidate_universe")
        != "every_cysteine_in_frozen_reference_proteome"
    ):
        raise RuntimeError("target candidate universe differs from frozen config")
    if not isinstance(sources, list) or len(sources) != 4:
        raise RuntimeError("materialization requires exactly four source domains")
    if not isinstance(loaded.get("input_registry"), str) or not isinstance(
        loaded.get("request_registry"), str
    ):
        raise RuntimeError("materialization registry paths are required")
    return loaded


def _source_payloads(
    cfg: dict[str, Any], inputs: dict[str, RegisteredMaterializationInput]
) -> list[dict[str, Any]]:
    sampling = dict(cfg["source_unlabeled_sampling"])
    sources: list[dict[str, Any]] = []
    for spec in cfg["sources"]:
        if not isinstance(spec, dict):
            raise RuntimeError("source specification must be a mapping")
        input_ids = tuple(str(item) for item in spec.get("input_ids", ()))
        if not input_ids or any(item not in inputs for item in input_ids):
            raise RuntimeError("source specification references an unregistered input")
        by_id = {item: inputs[item] for item in input_ids}
        parser = spec.get("parser")
        species = str(spec.get("species", ""))
        study = str(spec.get("study_accession", ""))
        if parser == "benchmark_v1_study_subset":
            proteome = _load_proteome(by_id["AT_PROTEOME_V1"].path)
            positives = _benchmark_positive_keys(by_id["AT_BENCHMARK_V1"].path, study)
            supplemental_proteome = (
                _load_proteome(by_id["AT_P42737_2_ISOFORM"].path)
                if "AT_P42737_2_ISOFORM" in by_id
                else {}
            )
            panther = load_panther_annotations(
                by_id["AT_PANTHER_V1"].path
            ).family_by_accession
        elif parser == "pxd072089_coordinate_verified_union":
            proteome = _load_proteome(by_id["RICE_PROTEOME_V1"].path)
            rice_table = parse_pxd072089_sites(
                by_id["RICE_SD01"].path,
                by_id["RICE_SD04"].path,
                proteome,
                ss_all_peptides_path=by_id["RICE_SS_ALL"].path,
            )
            positives = {
                (site.protein_accession, site.cys_position)
                for site in rice_table.sites
            }
            supplemental_proteome = {}
            panther = load_panther_annotations(
                by_id["RICE_PANTHER_V1"].path
            ).family_by_accession
        elif parser == "pxd063170_coordinate_verified_sites":
            proteome = load_ensembl_fungi_proteome(
                by_id["MAGNAPORTHE_PROTEOME_V1"].path
            )
            magnaporthe_table = parse_pxd063170_sites(
                by_id["MAGNAPORTHE_SITES_TSV"].path, proteome
            )
            positives = {
                (site.protein_accession, site.cys_position)
                for site in magnaporthe_table.sites
            }
            supplemental_proteome = {}
            panther = load_panther_annotations_by_orf_gene(
                by_id["MAGNAPORTHE_PANTHER_V1"].path,
                {accession: mg8_accession_to_gene(accession) for accession in proteome},
            ).family_by_accession
        else:
            raise RuntimeError(f"unsupported materialization source parser: {parser}")
        sources.append(
            build_source_payload(
                species=species,
                study_accession=study,
                rows=_source_rows_from_proteome(
                    proteome, positives, supplemental_proteome
                ),
                proteome={**proteome, **supplemental_proteome},
                source_input_sha256s=tuple(by_id[item].sha256 for item in input_ids),
                unlabeled_per_positive=int(sampling["unlabeled_per_positive"]),
                seed=int(sampling["deterministic_hash_seed"]),
                panther_subfamilies=panther,
            )
        )
    return sources


def _write_request_registry(
    request_path: Path,
    request_sha256: str,
    registry_path: Path,
) -> None:
    if registry_path.exists():
        raise FileExistsError(registry_path)
    registry_path.parent.mkdir(parents=True, exist_ok=True)
    relative = os.path.relpath(request_path.resolve(), registry_path.parent.resolve())
    registry_path.write_text(
        "input_id\tsource_accession\tsource_kind\tpath\tsha256\tscientific_use\n"
        "CROSS_CROP_TARGET_LABEL_FREE_REQUEST_V1\t"
        "task9a\tmaterialized_target_label_free_request\t"
        f"{Path(relative).as_posix()}\t{request_sha256}\t"
        "Task 9 target-label-free cross-crop runner input\n",
        encoding="utf-8",
    )


def materialize_cross_crop_request(
    config_path: Path,
    output_dir: Path,
    request_registry_path: Path | None = None,
) -> dict[str, object]:
    """Materialize the frozen Task 9A request and its hash registry once."""
    cfg = _checked_config(config_path)
    repo_root = config_path.resolve().parents[2]
    input_registry = (repo_root / str(cfg["input_registry"])).resolve()
    output_dir = output_dir.resolve()
    request_registry = (
        request_registry_path.resolve()
        if request_registry_path is not None
        else (repo_root / str(cfg["request_registry"])).resolve()
    )
    if output_dir.exists():
        raise FileExistsError(output_dir)
    if request_registry.exists():
        raise FileExistsError(request_registry)
    inputs = _load_registered_inputs(input_registry)
    target_spec = dict(cfg["target"])
    target_ids = tuple(str(item) for item in target_spec.get("input_ids", ()))
    if target_ids != ("TOMATO_PROTEOME_V1",):
        raise RuntimeError("target input policy differs from frozen config")
    target_input = inputs[target_ids[0]]
    forbidden = tuple(
        str(item).lower() for item in target_spec["forbidden_target_label_sources"]
    )
    if any(term in target_input.path.as_posix().lower() for term in forbidden):
        raise RuntimeError("target candidate input references a forbidden label source")
    target_payload = build_target_payload(
        _load_proteome(target_input.path), source_sha256=target_input.sha256
    )
    request: dict[str, Any] = {
        "sources": _source_payloads(cfg, inputs),
        "target_candidates": target_payload,
    }
    validate_label_free_target_payload(request["target_candidates"])
    output_dir.mkdir(parents=True, exist_ok=False)
    request_path = output_dir / "request.json"
    request_path.write_text(
        json.dumps(request, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    request_sha256 = hash_file(request_path, "sha256")
    _write_request_registry(request_path, request_sha256, request_registry)
    manifest: dict[str, object] = {
        "complete": True,
        "schema_version": 1,
        "task": "9A",
        "claim_class": "target_label_free_transfer_not_gate2",
        "materialization_config_sha256": hash_file(config_path, "sha256"),
        "input_registry_sha256": hash_file(input_registry, "sha256"),
        "request": {"path": request_path.name, "sha256": request_sha256},
        "source_domain_count": len(request["sources"]),
        "source_species": sorted({item["species"] for item in request["sources"]}),
        "source_row_counts": {
            f"{item['species']}|{item['study_accession']}": {
                "positive": item["labels"].count("positive"),
                "unlabeled": item["labels"].count("unlabeled"),
                "grouping": "homology_cluster"
                if "cluster_ids" in item
                else "protein",
            }
            for item in request["sources"]
        },
        "target": {
            "species": target_spec["species"],
            "candidate_count": len(target_payload["site_keys"]),
            "labels_visible": False,
            "source_sha256": target_input.sha256,
        },
        "labels_are_positive_or_unlabeled_only": True,
        "unobserved_sites_are_not_negatives": True,
        "tomato_label_input_used": False,
        "biological_values_modified": False,
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest
