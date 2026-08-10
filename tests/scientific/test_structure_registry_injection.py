"""RED (Phase S0, Task 1): the structure branch reads only the registry it was
explicitly handed.

Scientific integrity check. The structure-coverage experiment scores the *same*
benchmark rows twice — once against the 7-structure release and once against the
2,006-structure release — so the only way the comparison can be trusted is if
the registry is an explicit, per-call input rather than a module-level default
that silently follows whatever ``data/registry/alphafold_structures.tsv``
happens to contain today.

These tests monkeypatch **only** the registry auditor (a provenance boundary),
never a biological feature value, and the stub returns an *empty* audited
registry so no structure is fabricated: every row must then come back with
``has_structure=False``, which is the honest "no registered structure" answer.

The real-data test at the bottom audits both frozen release snapshots against
the registry data root. Registry ``local_path`` values are recorded relative to
the directory of the mutable registry (``data/registry``); a relocated snapshot
under ``data/registry/releases`` therefore needs that root supplied explicitly
rather than guessed, and guessing is exactly what a provenance auditor must
never do.

Expected RED: ``_structure_feature_vectors_with_mask()`` rejects
``structure_registry_path``.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from scripts import run_experiment

BENCHMARK = Path("data/processed/benchmark_v1/sites.tsv")
PROTEOME = Path("data/raw/references/arabidopsis_ref_proteome_v1.fasta")
RELEASES = Path("data/registry/releases")
REGISTRY_DATA_ROOT = Path("data/registry")
DEFAULT_REGISTRY = Path("data/registry/alphafold_structures.tsv")


def _real_benchmark_rows(limit: int) -> list[dict[str, str]]:
    """Read the first ``limit`` rows of the frozen benchmark verbatim."""
    rows: list[dict[str, str]] = []
    with BENCHMARK.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows.append(dict(row))
            if len(rows) == limit:
                break
    return rows


def test_feature_assembly_uses_only_explicit_registry(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    requested: list[Path] = []

    def record_registry(
        registry_path: Path = DEFAULT_REGISTRY,
        *,
        base_directory: Path | None = None,
    ) -> tuple[object, ...]:
        requested.append(registry_path)
        return ()

    monkeypatch.setattr(run_experiment, "audit_alphafold_structures", record_registry)
    explicit = tmp_path / "software-policy-registry.tsv"

    vectors, mask = run_experiment._structure_feature_vectors_with_mask(
        _real_benchmark_rows(1),
        PROTEOME,
        tmp_path,
        "injection",
        structure_registry_path=explicit,
    )

    assert requested == [explicit]
    # An empty audited registry means "no registered structure", never a
    # fabricated one: the mask must say so and the values must stay zero.
    assert mask == [False]
    assert vectors == [[0.0, 0.0]]


def test_feature_assembly_defaults_to_the_mutable_registry(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Existing v1 callers keep their behaviour: the default is unchanged."""
    requested: list[Path] = []

    def record_registry(
        registry_path: Path = DEFAULT_REGISTRY,
        *,
        base_directory: Path | None = None,
    ) -> tuple[object, ...]:
        requested.append(registry_path)
        return ()

    monkeypatch.setattr(run_experiment, "audit_alphafold_structures", record_registry)
    run_experiment._structure_feature_vectors_with_mask(
        _real_benchmark_rows(1), PROTEOME, tmp_path, "default"
    )
    assert requested == [DEFAULT_REGISTRY]


def test_branch_feature_assembly_propagates_the_explicit_registry(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    requested: list[Path] = []

    def record_registry(
        registry_path: Path = DEFAULT_REGISTRY,
        *,
        base_directory: Path | None = None,
    ) -> tuple[object, ...]:
        requested.append(registry_path)
        return ()

    monkeypatch.setattr(run_experiment, "audit_alphafold_structures", record_registry)
    explicit = RELEASES / "alphafold_structures_release_v1.tsv"

    branches = run_experiment._build_branch_features(
        _real_benchmark_rows(2),
        PROTEOME,
        tmp_path,
        "injection",
        need_esm=False,
        structure_registry_path=explicit,
        structure_registry_base=REGISTRY_DATA_ROOT,
    )

    assert requested == [explicit]
    assert branches.structure_mask == [False, False]


def test_frozen_release_snapshots_audit_against_the_registry_data_root() -> None:
    """Both relocated snapshots verify every registered local file's SHA256."""
    from plantpersulf.download.alphafold import audit_alphafold_structures

    v1 = audit_alphafold_structures(
        RELEASES / "alphafold_structures_release_v1.tsv",
        base_directory=REGISTRY_DATA_ROOT,
    )
    v2 = audit_alphafold_structures(
        RELEASES / "alphafold_structures_release_v2.tsv",
        base_directory=REGISTRY_DATA_ROOT,
    )
    assert len(v1) == 7
    assert len(v2) == 2006
    assert {source.accession for source in v1} < {source.accession for source in v2}
