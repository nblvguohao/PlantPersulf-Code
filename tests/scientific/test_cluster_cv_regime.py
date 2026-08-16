"""Repeated homology-cluster K-fold CV — the development-stability track.

User decision 2026-08-11: the Arabidopsis benchmark's development-stability
yardstick is a repeated MMseqs2-cluster K-fold CV (a cluster never crosses
train/val/test). This track is NOT the NC main evidence (that is the future
tomato lockbox cohort) and NOT Gate 2 evidence — its fold tags
(``cluster_cv_r<rep>f<fold>|...``) are structurally invisible to
``scripts/validate_external.py`` (see
tests/release/test_gate2_ignores_within_dataset_split_metrics.py).
"""

from __future__ import annotations

from scripts.run_experiment import _cluster_cv_chunks


def test_chunks_are_exhaustive_and_disjoint() -> None:
    cluster_ids = [f"C{i:03d}" for i in range(100)]
    chunks = _cluster_cv_chunks(cluster_ids, n_folds=5, seed=12345)
    assert len(chunks) == 5
    flat = [c for chunk in chunks for c in chunk]
    assert sorted(flat) == sorted(cluster_ids)
    for i in range(5):
        for j in range(i + 1, 5):
            assert set(chunks[i]).isdisjoint(chunks[j])


def test_chunks_are_deterministic_per_seed() -> None:
    cluster_ids = [f"C{i:03d}" for i in range(50)]
    assert _cluster_cv_chunks(cluster_ids, 5, 7) == _cluster_cv_chunks(
        cluster_ids, 5, 7
    )
    assert _cluster_cv_chunks(cluster_ids, 5, 7) != _cluster_cv_chunks(
        cluster_ids, 5, 8
    )


def test_chunks_roughly_balanced() -> None:
    cluster_ids = [f"C{i:03d}" for i in range(101)]
    chunks = _cluster_cv_chunks(cluster_ids, n_folds=10, seed=1)
    sizes = [len(c) for c in chunks]
    assert max(sizes) - min(sizes) <= 1
