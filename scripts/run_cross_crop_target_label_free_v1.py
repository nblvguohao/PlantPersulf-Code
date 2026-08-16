"""Run a hash-verified, target-label-free cross-crop request."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from plantpersulf.provenance.audit import assert_registered_input
from plantpersulf.provenance.hashing import hash_file
from plantpersulf.workflows.cross_crop_target_label_free import (
    SourceBatch,
    TargetCandidateBatch,
    run_cross_crop_target_label_free,
)


def source_batch_from_json(item: dict[str, Any]) -> SourceBatch:
    required = {
        "species",
        "study_accession",
        "feature_names",
        "features",
        "labels",
        "protein_ids",
        "source_sha256",
    }
    optional = {"cluster_ids", "batch_condition_ids"}
    if not required.issubset(item) or set(item) - required - optional:
        raise RuntimeError("source batch request schema mismatch")
    return SourceBatch(
        species=str(item["species"]),
        study_accession=str(item["study_accession"]),
        feature_names=tuple(str(value) for value in item["feature_names"]),
        features=tuple(
            tuple(float(value) for value in row) for row in item["features"]
        ),
        labels=tuple(str(value) for value in item["labels"]),
        protein_ids=tuple(str(value) for value in item["protein_ids"]),
        source_sha256=str(item["source_sha256"]),
        cluster_ids=tuple(str(value) for value in item.get("cluster_ids", ())),
        batch_condition_ids=tuple(
            str(value) for value in item.get("batch_condition_ids", ())
        ),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--request", type=Path, required=True)
    parser.add_argument("--request-registry", type=Path, required=True)
    parser.add_argument("--request-sha256", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    assert_registered_input(args.request, args.request_registry)
    if hash_file(args.request, "sha256") != args.request_sha256:
        raise RuntimeError("materialized request SHA256 mismatch")
    raw = json.loads(args.request.read_text(encoding="utf-8"))
    if set(raw) != {"sources", "target_candidates"}:
        raise RuntimeError("materialized request schema mismatch")
    sources = tuple(source_batch_from_json(item) for item in raw["sources"])
    target = raw["target_candidates"]
    if "labels" in target:
        raise RuntimeError("target label leakage or schema mismatch")
    batch = TargetCandidateBatch(
        tuple(target["site_keys"]),
        tuple(target["feature_names"]),
        tuple(tuple(float(x) for x in r) for r in target["features"]),
        str(target["source_sha256"]),
    )
    run_cross_crop_target_label_free(args.config, sources, batch, args.output_dir)


if __name__ == "__main__":
    main()
