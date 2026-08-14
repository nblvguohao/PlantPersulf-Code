#!/usr/bin/env python
"""Bulk AlphaFold DB structure download for the cross-species structural
conservation analysis (2026-07-23).

Reuses ``plantpersulf.download.alphafold`` (fail-closed, per-accession
provenance, handles isoform/404 as clean "no structure" results) to fetch
AlphaFold models for every persulfidated protein across the three
independent species tracks (Arabidopsis benchmark_v1, PXD072089 rice,
PXD063170 Magnaporthe — the latter bridged from MG8 gene IDs to UniProt
accessions via the PANTHER gene-name join already used for the
conservation analysis).

Idempotent: accessions already present in the registry are skipped, so the
script can be safely re-run to resume an interrupted batch.

Usage::

    python scripts/download_alphafold_structures_bulk.py \
        --accessions data/raw/references/persulfidated_protein_accessions_v1.txt
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

REGISTRY = Path("data/registry/alphafold_structures.tsv")
DEST_DIR = Path("data/raw/alphafold")


def select_accessions_to_fetch(
    accessions: list[str],
    already: set[str],
    force: bool,
) -> list[str]:
    """De-duplicate the batch, preserving order, and drop already-registered
    accessions unless ``force`` is set.

    ``force`` exists for registry refreshes: ``register_alphafold_structure``
    replaces the row for an accession, so refetching an accession that is
    already registered upgrades it in place (used to lift structures pinned
    at a superseded AlphaFold model version onto the current one).
    """
    seen: set[str] = set()
    todo: list[str] = []
    for accession in accessions:
        if accession in seen:
            continue
        seen.add(accession)
        if not force and accession in already:
            continue
        todo.append(accession)
    return todo


def run_bulk_download(
    accessions_path: Path,
    sleep_seconds: float = 0.1,
    limit: int | None = None,
    force: bool = False,
) -> dict[str, int]:
    from plantpersulf.download.alphafold import (
        audit_alphafold_structures,
        fetch_alphafold_structure,
        register_alphafold_structure,
    )

    accessions = [
        line.strip()
        for line in accessions_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if limit is not None:
        accessions = accessions[:limit]

    already = (
        {s.accession for s in audit_alphafold_structures(REGISTRY)}
        if REGISTRY.is_file()
        else set()
    )
    todo = select_accessions_to_fetch(accessions, already, force)
    already_in_batch = len(accessions) - len(todo)
    print(
        f"{len(accessions)} total accessions in this batch, "
        f"{already_in_batch} already registered, {len(todo)} to fetch"
    )

    counts = {"downloaded": 0, "not_found": 0, "isoform": 0, "error": 0}
    for i, accession in enumerate(todo, start=1):
        dest = DEST_DIR / f"AF-{accession}-F1-model.pdb"
        try:
            result = fetch_alphafold_structure(accession, dest)
        except RuntimeError as exc:
            counts["error"] += 1
            print(f"  [{i}/{len(todo)}] {accession}: ERROR {exc}")
            continue

        if result.status == "downloaded":
            register_alphafold_structure(result, REGISTRY)
            counts["downloaded"] += 1
        elif result.status == "not_found":
            counts["not_found"] += 1
        else:
            counts["isoform"] += 1

        if i % 100 == 0 or i == len(todo):
            print(
                f"  [{i}/{len(todo)}] downloaded={counts['downloaded']} "
                f"not_found={counts['not_found']} "
                f"isoform={counts['isoform']} error={counts['error']}"
            )
        time.sleep(sleep_seconds)

    print(f"\nFinal: {counts}")
    return counts


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="bulk AlphaFold DB structure download")
    p.add_argument("--accessions", type=Path, required=True)
    p.add_argument("--sleep-seconds", type=float, default=0.1)
    p.add_argument("--limit", type=int, default=None)
    p.add_argument(
        "--force",
        action="store_true",
        help="refetch accessions already in the registry (in-place upgrade)",
    )
    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    run_bulk_download(
        accessions_path=args.accessions,
        sleep_seconds=args.sleep_seconds,
        limit=args.limit,
        force=args.force,
    )
