#!/usr/bin/env python
"""Freeze the main AlphaFold structure registry into a release snapshot.

W1 (tomato structure coverage) needs a frozen, registered structure-registry
release (v3) that the campaign runners can point at. This script mirrors how
release v2 was produced:

1. audit every registered PDB (SHA256 + size) via the fail-closed auditor;
2. write a sorted-by-accession snapshot tsv to ``data/registry/releases/``;
3. compute the snapshot SHA256;
4. append one row to ``data/registry/model_inputs.tsv`` (refuses to duplicate
   an existing input_id, so re-runs are safe).

Usage::

    PYTHONPATH='src' python scripts/accuracy_campaign/snapshot_structure_release.py \\
        --input-id ALPHAFOLD_STRUCTURES_RELEASE_V3 \\
        --scientific-use multispecies_v2_accuracy_campaign_w1_tomato_structures

The snapshot is append-only with respect to the main registry: nothing in the
main registry or existing releases is modified.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
from datetime import datetime, timezone
from pathlib import Path

REGISTRY = Path("data/registry/alphafold_structures.tsv")
RELEASES_DIR = Path("data/registry/releases")
MODEL_INPUTS = Path("data/registry/model_inputs.tsv")

MODEL_INPUTS_FIELDS = (
    "input_id",
    "source_accession",
    "source_kind",
    "container_url",
    "path",
    "retrieved_at",
    "size_bytes",
    "sha256",
    "data_level",
    "scientific_use",
    "derivation",
)


def read_registry_rows(registry_path: Path) -> list[dict[str, str]]:
    with registry_path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sorted_registry_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Sort registry rows by accession (deterministic snapshot order)."""
    return sorted(rows, key=lambda row: row["accession"])


def build_snapshot_path(input_id: str) -> Path:
    return RELEASES_DIR / f"{input_id.lower()}.tsv"


def registry_fieldnames() -> tuple[str, ...]:
    from plantpersulf.download.alphafold import STRUCTURE_FIELDS

    return STRUCTURE_FIELDS


def model_inputs_has_input_id(input_id: str) -> bool:
    if not MODEL_INPUTS.is_file():
        return False
    with MODEL_INPUTS.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return any(row["input_id"] == input_id for row in reader)


def append_model_inputs_row(
    *,
    input_id: str,
    snapshot_path: Path,
    size_bytes: int,
    sha256: str,
    scientific_use: str,
    registry_row_count: int,
) -> None:
    row = {
        "input_id": input_id,
        "source_accession": "derived:ALPHAFOLD_STRUCTURE_REGISTRY",
        "source_kind": "frozen_structure_registry_snapshot",
        "container_url": "local:data/registry/alphafold_structures.tsv",
        "path": str(snapshot_path.relative_to("data/registry")),
        "retrieved_at": datetime.now(timezone.utc).isoformat(),
        "size_bytes": str(size_bytes),
        "sha256": sha256,
        "data_level": "C",
        "scientific_use": scientific_use,
        "derivation": (
            f"Frozen snapshot of {registry_row_count} SHA256-registered "
            "AlphaFold structure records; absent proteins remain explicitly "
            "masked."
        ),
    }
    with MODEL_INPUTS.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MODEL_INPUTS_FIELDS, delimiter="\t")
        writer.writerow(row)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-id", required=True)
    parser.add_argument("--scientific-use", required=True)
    args = parser.parse_args(argv)

    from plantpersulf.download.alphafold import audit_alphafold_structures

    # Verification pass only: fail closed if any registered PDB no longer
    # matches its recorded SHA256/size. The snapshot itself keeps the
    # registry-relative ``local_path`` values (release convention).
    audit = audit_alphafold_structures(REGISTRY)
    if not audit:
        raise RuntimeError("registry audit returned no structures")

    rows = sorted_registry_rows(read_registry_rows(REGISTRY))
    if not rows:
        raise RuntimeError("registry has no rows to snapshot")

    snapshot_path = build_snapshot_path(args.input_id)
    snapshot_path.parent.mkdir(parents=True, exist_ok=True)
    with snapshot_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=registry_fieldnames(), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)
    snapshot_sha256 = hashlib.sha256(snapshot_path.read_bytes()).hexdigest()
    print(
        f"snapshot {snapshot_path}: {len(rows)} rows, "
        f"size={snapshot_path.stat().st_size}, sha256={snapshot_sha256}"
    )

    if model_inputs_has_input_id(args.input_id):
        print(f"model_inputs.tsv already carries {args.input_id}; not appended")
        return
    append_model_inputs_row(
        input_id=args.input_id,
        snapshot_path=snapshot_path,
        size_bytes=snapshot_path.stat().st_size,
        sha256=snapshot_sha256,
        scientific_use=args.scientific_use,
        registry_row_count=len(rows),
    )
    print(f"appended {args.input_id} to {MODEL_INPUTS}")


if __name__ == "__main__":
    main()
