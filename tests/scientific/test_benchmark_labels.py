"""RED (Task 5): positive-unlabeled persulfidation benchmark v1.

Consumes the deterministic persulfidation site outputs and a SHA256-pinned
Arabidopsis reference proteome. The benchmark must label every cysteine as
positive (has published site-level persulfidation evidence) or unlabeled;
the label "negative" must never appear, and every positive record must be
traceable to its source study, file, and SHA256.

Expected RED: ``plantpersulf.benchmark.labels`` does not exist yet, so
collection fails with ImportError.
"""

from __future__ import annotations

import shutil
import tempfile
from pathlib import Path

import pytest

from plantpersulf.benchmark.labels import (  # RED: module missing
    BenchmarkLabels,
    build_benchmark_labels,
)

_SITE_CACHE: Path | None = None
_LABELS_CACHE: BenchmarkLabels | None = None


@pytest.fixture(scope="module")
def cached_data(
    tmp_path_factory: pytest.TempPathFactory,
) -> tuple[BenchmarkLabels, Path]:
    """Build once, return (labels, site_output_root)."""
    global _LABELS_CACHE, _SITE_CACHE
    if _LABELS_CACHE is not None:
        return _LABELS_CACHE, _SITE_CACHE  # type: ignore[return-value]
    root = tmp_path_factory.mktemp("bench")
    from plantpersulf.proteomics.persulfidation_publish import (
        publish_persulfidation_sites,
    )

    publish_persulfidation_sites("PXD006140", output_root=root)
    publish_persulfidation_sites("PXD024061", output_root=root)
    _SITE_CACHE = root
    _LABELS_CACHE = build_benchmark_labels(
        site_output_root=root,
        proteome_path=Path("data/raw/references/arabidopsis_ref_proteome_v1.fasta"),
        output_directory=root / "benchmark_v1",
    )
    return _LABELS_CACHE, _SITE_CACHE


@pytest.fixture(scope="module")
def real_benchmark(cached_data: tuple[BenchmarkLabels, Path]) -> BenchmarkLabels:
    return cached_data[0]


def test_labels_are_exclusively_positive_or_unlabeled(
    real_benchmark: BenchmarkLabels,
) -> None:
    labels = real_benchmark

    observed = {labels.sites[i].label for i in range(len(labels))}
    assert observed <= {"positive", "unlabeled"}
    assert "negative" not in observed
    assert "positive" in observed
    assert labels.positive_count >= 320  # PXD006140 alone has 320 distinct
    assert labels.unlabeled_count > 0


def test_every_positive_has_traceable_experimental_evidence(
    real_benchmark: BenchmarkLabels,
) -> None:
    labels = real_benchmark

    positives = [(i, s) for i, s in enumerate(labels) if s.label == "positive"]
    assert positives
    for _, site in positives:
        assert site.study_accession in {"PXD006140", "PXD024061"}
        assert site.evidence_level == "site_ms"
        assert site.source_sha256
        assert site.protein_accession
        assert site.cys_position_in_protein > 0


def test_unlabeled_sites_have_no_experimental_evidence_claim(
    real_benchmark: BenchmarkLabels,
) -> None:
    labels = real_benchmark

    unlabeled = [(i, s) for i, s in enumerate(labels) if s.label == "unlabeled"]
    assert unlabeled
    for _, site in unlabeled:
        # An unlabeled site must not carry an evidence claim; it is simply
        # a cysteine whose persulfidation status is unknown.
        assert site.label == "unlabeled"
        assert site.study_accession == ""


def test_benchmark_is_deterministic(
    cached_data: tuple[BenchmarkLabels, Path],
) -> None:
    first, site_root = cached_data
    tmp = Path(tempfile.mkdtemp())
    try:
        from plantpersulf.benchmark.labels import build_benchmark_labels

        second = build_benchmark_labels(
            site_output_root=site_root,
            proteome_path=Path("data/raw/references/arabidopsis_ref_proteome_v1.fasta"),
            output_directory=tmp / "benchmark_v2",
        )
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    assert len(first) == len(second)
    assert first.positive_count == second.positive_count
    assert first.unlabeled_count == second.unlabeled_count
    assert first.negative_count == 0
    assert second.negative_count == 0


def test_benchmark_refuses_missing_proteome(
    cached_data: tuple[BenchmarkLabels, Path],
) -> None:
    site_root = cached_data[1]
    with pytest.raises((RuntimeError, OSError, FileNotFoundError)):
        build_benchmark_labels(
            site_output_root=site_root,
            proteome_path=Path("data/raw/references/nonexistent.fasta"),
            output_directory=Path(tempfile.mkdtemp()) / "benchmark_v1",
        )
