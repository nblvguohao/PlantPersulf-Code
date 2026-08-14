"""Bulk AlphaFold download: which accessions a batch actually fetches.

The batch is idempotent by default so an interrupted run can resume, but a
stale registry entry (e.g. a structure registered at a superseded AlphaFold
model version) has to be refetchable without hand-editing the registry —
``register_alphafold_structure`` already replaces the row for an accession,
so the only thing standing in the way is the skip rule.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

from download_alphafold_structures_bulk import (  # noqa: E402
    select_accessions_to_fetch,
)


def test_already_registered_accessions_are_skipped_by_default() -> None:
    todo = select_accessions_to_fetch(["A", "B", "C"], already={"B"}, force=False)
    assert todo == ["A", "C"]


def test_force_refetches_everything_in_the_batch() -> None:
    todo = select_accessions_to_fetch(["A", "B", "C"], already={"A", "B"}, force=True)
    assert todo == ["A", "B", "C"]


def test_selection_preserves_input_order_and_drops_duplicates() -> None:
    todo = select_accessions_to_fetch(["C", "A", "C", "B"], already=set(), force=False)
    assert todo == ["C", "A", "B"]


def test_empty_batch_is_empty() -> None:
    assert select_accessions_to_fetch([], already={"A"}, force=True) == []
