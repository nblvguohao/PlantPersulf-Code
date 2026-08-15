"""Unit tests for the parallel-batch bulk AlphaFold downloader path."""

from __future__ import annotations

from pathlib import Path

from scripts.download_alphafold_structures_bulk import (
    _fetch_batch,
    iter_batches,
)


class _FakeResult:
    def __init__(self, accession: str, status: str) -> None:
        self.accession = accession
        self.status = status


class _Boom(RuntimeError):
    pass


def test_iter_batches_chunks_exactly() -> None:
    assert list(iter_batches(["a", "b", "c", "d", "e"], 2)) == [
        ["a", "b"],
        ["c", "d"],
        ["e"],
    ]
    assert list(iter_batches([], 2)) == []


def test_fetch_batch_returns_input_order_with_workers() -> None:
    order: list[str] = []

    def fetch(accession: str, dest: Path):
        order.append(accession)
        if accession == "bad":
            raise _Boom("transport failure")
        return _FakeResult(accession, "downloaded")

    pairs = _fetch_batch(["a", "bad", "c", "d"], Path("."), fetch, workers=4)
    assert [accession for accession, _ in pairs] == ["a", "bad", "c", "d"]
    assert isinstance(pairs[1][1], _Boom)
    assert isinstance(pairs[0][1], _FakeResult)
    # completion order is not asserted; input order is the contract


def test_fetch_batch_sequential_path() -> None:
    calls: list[str] = []

    def fetch(accession: str, dest: Path):
        calls.append(accession)
        return _FakeResult(accession, "not_found")

    pairs = _fetch_batch(["x", "y"], Path("."), fetch, workers=1)
    assert calls == ["x", "y"]
    assert [result.status for _, result in pairs] == ["not_found", "not_found"]
