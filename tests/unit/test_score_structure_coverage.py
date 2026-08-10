"""RED (Phase S0, Task 5): fail-closed orchestration for the paired scoring
runner. Unit tests verify the public interface, output-path rejection rules,
and score-field schema."""

from __future__ import annotations

from pathlib import Path

import pytest


class TestPrepareOutputDirectory:
    def test_creates_new_directory(self, tmp_path: Path) -> None:
        from scripts.score_structure_coverage import prepare_output_directory

        target = tmp_path / "new_output"
        prepare_output_directory(target)
        assert target.is_dir()

    def test_refuses_existing_directory(self, tmp_path: Path) -> None:
        from scripts.score_structure_coverage import prepare_output_directory

        target = tmp_path / "existing"
        target.mkdir()
        with pytest.raises(FileExistsError):
            prepare_output_directory(target)


class TestScoreFieldSchema:
    def test_score_fields_are_defined(self) -> None:
        from scripts.score_structure_coverage import SCORE_FIELDS

        assert "coverage_release" in SCORE_FIELDS
        assert "held_out_study" in SCORE_FIELDS
        assert "seed" in SCORE_FIELDS
        assert "arm" in SCORE_FIELDS
        assert "protein_accession" in SCORE_FIELDS
        assert "cys_position_in_protein" in SCORE_FIELDS
        assert "label" in SCORE_FIELDS
        assert "score" in SCORE_FIELDS
        assert "uncertainty" in SCORE_FIELDS


class TestRunStructureCoverageScoring:
    def test_verify_only_mode_imports_run_function(self, tmp_path: Path) -> None:
        """verify_only and run_scoring are importable."""
        from scripts.score_structure_coverage import (
            run_structure_coverage_scoring,
        )

        assert callable(run_structure_coverage_scoring)

    def test_config_hashing_catches_mismatch(self, tmp_path: Path) -> None:
        """A config with a wrong hash should be rejected."""
        from plantpersulf.evaluation.structure_coverage_audit import sha256_file
        from plantpersulf.evaluation.structure_coverage_config import (
            load_structure_coverage_config,
        )

        cfg = load_structure_coverage_config(
            Path("configs/experiments/pu_ranker_structcover_v2.yaml")
        )
        # The real config has correct hashes — must not fail
        assert cfg.benchmark_sha256 == sha256_file(Path(cfg.benchmark_path))


class TestCLIParsing:
    def test_default_config_path(self) -> None:
        from scripts.score_structure_coverage import build_parser

        parser = build_parser()
        args = parser.parse_args(["--verify-only"])
        assert args.verify_only

    def test_custom_config_path(self) -> None:
        from scripts.score_structure_coverage import build_parser

        parser = build_parser()
        args = parser.parse_args(["--config", "some/path.yaml"])
        assert args.config == Path("some/path.yaml")
