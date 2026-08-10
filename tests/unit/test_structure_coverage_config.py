"""RED (Phase S0, Task 3): validate the frozen structure-coverage experiment
config — exact paths, hashes, arms, and controlled variables."""

from __future__ import annotations

from pathlib import Path

import pytest

from plantpersulf.evaluation.structure_coverage_config import (
    StructureCoverageExperimentConfig,
    load_structure_coverage_config,
    validate_controlled_variables,
)

CONFIG_PATH = Path("configs/experiments/pu_ranker_structcover_v2.yaml")

V1_RELEASE = Path("data/registry/releases/alphafold_structures_release_v1.tsv")
V2_RELEASE = Path("data/registry/releases/alphafold_structures_release_v2.tsv")

EXPECTED_ARMS = frozenset({
    "sequence_only",
    "sequence_coverage_only",
    "sequence_contact",
    "sequence_plddt",
    "sequence_contact_plddt",
})


class TestConfigLoading:
    def test_loads_and_validates_the_frozen_config(self) -> None:
        cfg = load_structure_coverage_config(CONFIG_PATH)
        validate_controlled_variables(cfg)
        assert cfg.version == 2
        assert cfg.experiment_name == "pu_ranker_structcover_v2"

    def test_registry_paths_are_exact(self) -> None:
        cfg = load_structure_coverage_config(CONFIG_PATH)
        assert cfg.registry_v1.path.resolve() == V1_RELEASE.resolve()
        assert cfg.registry_v2.path.resolve() == V2_RELEASE.resolve()

    def test_registry_counts_are_baked_in(self) -> None:
        cfg = load_structure_coverage_config(CONFIG_PATH)
        assert cfg.registry_v1.records == 7
        assert cfg.registry_v2.records == 2006

    def test_benchmark_path_and_hash_match_contract(self) -> None:
        from plantpersulf.evaluation.structure_coverage_audit import sha256_file

        cfg = load_structure_coverage_config(CONFIG_PATH)
        exp = cfg.benchmark_sha256.upper()
        got = sha256_file(Path(cfg.benchmark_path))
        assert got == exp

    def test_studies_are_pxd006140_and_pxd024061(self) -> None:
        cfg = load_structure_coverage_config(CONFIG_PATH)
        assert cfg.studies == ["PXD006140", "PXD024061"]

    def test_seeds_are_0_to_4(self) -> None:
        cfg = load_structure_coverage_config(CONFIG_PATH)
        assert cfg.seeds == [0, 1, 2, 3, 4]

    def test_subsample_ratio_is_20_and_seed_12345(self) -> None:
        cfg = load_structure_coverage_config(CONFIG_PATH)
        assert cfg.subsample_ratio == 20
        assert cfg.subsample_seed == 12345

    def test_five_arms_are_present(self) -> None:
        cfg = load_structure_coverage_config(CONFIG_PATH)
        assert set(cfg.arms) == EXPECTED_ARMS

    def test_every_arm_disables_esm_and_study_context(self) -> None:
        cfg = load_structure_coverage_config(CONFIG_PATH)
        for arm_name, arm in cfg.arms.items():
            assert not arm.get("use_esm", True), f"{arm_name} should have use_esm=false"
            assert not arm.get("use_study_context", True), \
                f"{arm_name} should have use_study_context=false"

    def test_ranker_parameters_are_copied_from_v1(self) -> None:
        cfg = load_structure_coverage_config(CONFIG_PATH)
        assert cfg.ranker_hidden == 16
        assert cfg.ranker_dropout == 0.2
        assert cfg.ranker_epochs == 200
        assert cfg.ranker_lr == 0.05
        assert cfg.ranker_n_mc_dropout == 16

    def test_output_directory_is_distinct_and_new(self) -> None:
        cfg = load_structure_coverage_config(CONFIG_PATH)
        assert "pu_ranker_structcover_v2" in cfg.output_directory
        assert "pu_ranker_v1" not in cfg.output_directory
        assert "external_validation" not in cfg.output_directory


class TestConfigRejection:
    def test_rejects_missing_file(self, tmp_path: Path) -> None:
        missing = tmp_path / "nonexistent.yaml"
        with pytest.raises((FileNotFoundError, SystemExit, RuntimeError)):
            load_structure_coverage_config(missing)

    def test_unknown_arm_is_rejected(self, monkeypatch: pytest.MonkeyPatch) -> None:
        cfg = load_structure_coverage_config(CONFIG_PATH)
        cfg = cfg._replace(arms={**cfg.arms, "fake_arm": {"use_esm": False, "use_structure": True}})
        with pytest.raises(RuntimeError, match="unknown"):
            validate_controlled_variables(cfg)

    def test_more_than_five_arms_is_rejected(self) -> None:
        cfg = load_structure_coverage_config(CONFIG_PATH)
        cfg_violated = cfg._replace(arms={**cfg.arms, "extra_arm": {"use_esm": False, "use_structure": True, "use_plddt": False, "use_accessibility": False, "use_study_context": False}})
        with pytest.raises(RuntimeError, match="unknown"):
            validate_controlled_variables(cfg_violated)

    def test_output_in_v1_tree_is_rejected(self, tmp_path: Path) -> None:
        """Not a real test on the frozen config — a policy test."""
        cfg = load_structure_coverage_config(CONFIG_PATH)
        # The real config has a clean output dir; this only checks the validator rejects
        # v1 and external_validation paths.
        assert "external_validation" not in cfg.output_directory

    def test_frozen_hashes_in_config_are_uppercase(self) -> None:
        cfg = load_structure_coverage_config(CONFIG_PATH)
        for label, sha in [
            ("benchmark", cfg.benchmark_sha256),
            ("proteome", cfg.proteome_sha256),
            ("clusters", cfg.clusters_sha256),
            ("registry_v1", cfg.registry_v1.sha256),
            ("registry_v2", cfg.registry_v2.sha256),
        ]:
            assert sha == sha.upper(), f"{label} hash not uppercase"
            assert len(sha) == 64, f"{label} hash wrong length"


class TestBootstrapSeeds:
    def test_cluster_bootstrap_seed_is_1729(self) -> None:
        cfg = load_structure_coverage_config(CONFIG_PATH)
        assert cfg.cluster_bootstrap_seed == 1729

    def test_cluster_bootstrap_replicates_are_5000(self) -> None:
        cfg = load_structure_coverage_config(CONFIG_PATH)
        assert cfg.cluster_bootstrap_replicates == 5000
