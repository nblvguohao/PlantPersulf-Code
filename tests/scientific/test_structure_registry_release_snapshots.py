"""RED (Phase S0, Task 1): the two frozen AlphaFold registry releases exist,
carry the exact contracted bytes, and stand in a strict subset relationship.

Scientific integrity check. The whole structure-coverage experiment compares a
7-structure registry (the coverage that existed when Gate 2 was frozen) against
the 2,006-structure registry that later cross-species work produced. That
comparison is only meaningful if both registries are *immutable byte snapshots*
of real, previously registered downloads:

* release v1 must be the exact blob at ``c827277^:data/registry/
  alphafold_structures.tsv`` (extracted with Git plumbing, never reserialized
  from parsed rows);
* release v2 must be a byte-for-byte copy of the currently audited mutable
  registry;
* every release-v1 record must reappear unchanged in release v2, so the
  expansion is provably additive and no earlier structure was re-downloaded,
  re-versioned, or silently edited.

The SHA256 values below are the frozen scientific contract
(``docs/superpowers/plans/2026-07-27-pu-ranker-v2-structure-coverage.md``). A
mismatch is a blocking provenance failure, never something to repair by
rewriting the snapshot.

Expected RED: the release snapshots do not exist yet.
"""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path

RELEASES = Path("data/registry/releases")

V1_SHA256 = "E09D18D334A3CF41233E59EB6CF66C9611EF9B13CD407C30770F608674128326"
V2_SHA256 = "BEE2D30ED28D754A7162283BB3C6080928DBC6A4CFA563958C48BD0146188075"

V1_RECORDS = 7
V2_RECORDS = 2006

# The registry schema is owned by plantpersulf.download.alphafold.STRUCTURE_FIELDS;
# it is restated here so a silent column change is caught by this test too.
EXPECTED_COLUMNS = (
    "accession",
    "model_version",
    "source_url",
    "local_path",
    "retrieved_at",
    "size_bytes",
    "sha256",
    "downloader_version",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest().upper()


def _read_registry_records(path: Path) -> list[dict[str, str]]:
    """Read a frozen release snapshot without any repair or normalisation."""
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert tuple(reader.fieldnames or ()) == EXPECTED_COLUMNS
        return [dict(row) for row in reader]


def _canonical_record(row: dict[str, str]) -> tuple[tuple[str, str], ...]:
    return tuple(sorted((key, row[key]) for key in EXPECTED_COLUMNS))


def test_frozen_structure_registries_have_expected_identity() -> None:
    v1 = RELEASES / "alphafold_structures_release_v1.tsv"
    v2 = RELEASES / "alphafold_structures_release_v2.tsv"

    assert _sha256(v1) == V1_SHA256
    assert _sha256(v2) == V2_SHA256

    rows_v1 = _read_registry_records(v1)
    rows_v2 = _read_registry_records(v2)

    accessions_v1 = {row["accession"] for row in rows_v1}
    accessions_v2 = {row["accession"] for row in rows_v2}
    assert len(rows_v1) == len(accessions_v1) == V1_RECORDS
    assert len(rows_v2) == len(accessions_v2) == V2_RECORDS

    assert {_canonical_record(row) for row in rows_v1} <= {
        _canonical_record(row) for row in rows_v2
    }


def test_release_v1_is_a_strict_subset_of_release_v2() -> None:
    """The expansion must be additive: v1 accessions all survive into v2, and
    v2 adds strictly more (otherwise there is nothing to compare)."""
    rows_v1 = _read_registry_records(RELEASES / "alphafold_structures_release_v1.tsv")
    rows_v2 = _read_registry_records(RELEASES / "alphafold_structures_release_v2.tsv")

    accessions_v1 = {row["accession"] for row in rows_v1}
    accessions_v2 = {row["accession"] for row in rows_v2}
    assert accessions_v1 < accessions_v2


def test_release_v2_records_survive_in_the_mutable_registry() -> None:
    """Release v2 is a copy, not a re-serialisation, of the audited registry
    at freeze time; the live registry may legally GROW afterwards (new
    independent accessions, e.g. the tomato-local track appended 76 tomato
    structures on 2026-08-10) but no frozen record may be altered or
    removed — a frozen experiment input stays byte-identical for
    reproducibility, while additive growth is legitimate."""
    live = Path("data/registry/alphafold_structures.tsv")
    frozen = RELEASES / "alphafold_structures_release_v2.tsv"
    frozen_records = {_canonical_record(r) for r in _read_registry_records(frozen)}
    live_records = {_canonical_record(r) for r in _read_registry_records(live)}
    # every frozen record still present, byte-for-byte field values
    assert frozen_records <= live_records
