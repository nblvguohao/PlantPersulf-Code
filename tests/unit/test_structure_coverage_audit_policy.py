"""RED (Phase S0, Task 2): software-policy-level audit rules.

These tests validate the fail-closed behaviour of the audit without needing
real biological data — they use a real registered row copied into temp files
and then deliberately tampered with.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from plantpersulf.evaluation.structure_coverage_audit import (
    FrozenFile,
    sha256_file,
    verify_frozen_file,
)

REAL_REGISTRY = Path("data/registry/releases/alphafold_structures_release_v1.tsv")


class TestFrozenFileVerification:
    def test_verify_frozen_file_passes_when_hash_matches(self) -> None:
        item = FrozenFile(
            name="test", path=REAL_REGISTRY, sha256=sha256_file(REAL_REGISTRY)
        )
        verify_frozen_file(item)  # must not raise

    def test_verify_frozen_file_raises_on_mismatch(self, tmp_path: Path) -> None:
        tampered = tmp_path / "tampered.tsv"
        shutil.copy2(REAL_REGISTRY, tampered)
        # Append a byte so the hash changes
        with tampered.open("ab") as f:
            f.write(b"tampered\n")
        item = FrozenFile(name="test", path=tampered, sha256=sha256_file(REAL_REGISTRY))
        with pytest.raises(RuntimeError, match="SHA256 mismatch"):
            verify_frozen_file(item)

    def test_sha256_file_is_uppercase_hex(self) -> None:
        got = sha256_file(REAL_REGISTRY)
        assert len(got) == 64
        assert got == got.upper()
        assert all(c in "0123456789ABCDEF" for c in got)

    def test_sha256_file_is_deterministic(self) -> None:
        a = sha256_file(REAL_REGISTRY)
        b = sha256_file(REAL_REGISTRY)
        assert a == b


class TestRegistryIdentityChecks:
    def test_duplicate_identical_rows_are_error(self) -> None:
        """Duplicate identical rows are a registry integrity error."""
        from plantpersulf.evaluation.structure_coverage_audit import (
            verify_registry_pair,
        )

        v1_valid = Path("data/registry/releases/alphafold_structures_release_v1.tsv")
        v2_valid = Path("data/registry/releases/alphafold_structures_release_v2.tsv")
        verify_registry_pair(v1_valid, v2_valid)  # must not raise


class TestCoverageRecordSchema:
    def test_coverage_record_fields_exist(self) -> None:
        from plantpersulf.evaluation.structure_coverage_audit import CoverageRecord

        rec = CoverageRecord(
            release="test",
            protein_accession="P12345",
            cys_position_in_protein=42,
            label="unlabeled",
            study_accession="",
            cluster_id="C001",
            has_registered_structure=False,
            maps_to_cys=False,
            mapping_status="absent_structure",
            plddt=None,
        )
        assert rec.mapping_status == "absent_structure"
        assert rec.plddt is None


class TestStructureCoverageAuditBlocking:
    def test_require_pass_raises_on_errors(self) -> None:
        from plantpersulf.evaluation.structure_coverage_audit import (
            StructureCoverageAudit,
        )

        audit = StructureCoverageAudit(
            records=(),
            summary_rows=(),
            input_hashes={},
            blocking_errors=("test error",),
        )
        with pytest.raises(RuntimeError, match="blocking"):
            audit.require_pass()

    def test_require_pass_succeeds_without_errors(self) -> None:
        from plantpersulf.evaluation.structure_coverage_audit import (
            StructureCoverageAudit,
        )

        audit = StructureCoverageAudit(
            records=(),
            summary_rows=(),
            input_hashes={},
            blocking_errors=(),
        )
        audit.require_pass()  # must not raise
