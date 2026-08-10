"""Phase S0, Task 3 — load and validate the frozen structure-coverage config.

Every controlled variable is read once and validated against the frozen
contract before any scoring begins. Unknown arms, missing fields, hash
mismatches, an output inside v1/external-validation, or any deviation from
the exact five-arm ESM-free design is a ``RuntimeError``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml

ALLOWED_ARMS = frozenset(
    {
        "sequence_only",
        "sequence_coverage_only",
        "sequence_contact",
        "sequence_plddt",
        "sequence_contact_plddt",
    }
)

FORBIDDEN_OUTPUT_ROOTS = frozenset(
    {
        "pu_ranker_v1",
        "external_validation",
    }
)


@dataclass(frozen=True)
class RegistryReleaseConfig:
    release_name: str
    path: Path
    sha256: str
    records: int


@dataclass(frozen=True)
class StructureCoverageExperimentConfig:
    """Every field is read-only after construction. Use ``_replace()`` only in tests."""

    version: int
    experiment_name: str
    benchmark_path: str
    benchmark_sha256: str
    proteome_path: str
    proteome_sha256: str
    clusters_path: str
    clusters_sha256: str
    legacy_config_path: str
    legacy_config_sha256: str
    registry_v1: RegistryReleaseConfig
    registry_v2: RegistryReleaseConfig
    studies: list[str]
    seeds: list[int]
    cluster_bootstrap_seed: int
    cluster_bootstrap_replicates: int
    permutation_seed: int
    permutation_replicates: int
    subsample_ratio: int
    subsample_seed: int
    ranker_hidden: int
    ranker_dropout: float
    ranker_epochs: int
    ranker_lr: float
    ranker_n_mc_dropout: int
    arms: dict[str, dict[str, bool]]
    output_directory: str

    def _replace(self, **kwargs: object) -> StructureCoverageExperimentConfig:
        return StructureCoverageExperimentConfig(
            **{
                **self.__dict__,
                **kwargs,
            }
        )


def load_structure_coverage_config(
    path: Path,
) -> StructureCoverageExperimentConfig:
    """Load and return a validated frozen config.

    Raises ``RuntimeError`` on any structural or frozen-hash mismatch."""
    if not path.is_file():
        raise FileNotFoundError(f"config not found: {path}")

    try:
        loaded: object = yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise RuntimeError(f"invalid YAML in config: {path}") from exc

    if not isinstance(loaded, dict):
        raise RuntimeError(f"config is not a mapping: {path}")
    if loaded.get("version") != 2:
        raise RuntimeError(f"config version must be 2: {path}")

    exp = loaded.get("experiment", {})
    reg = loaded.get("registries", {})
    splits = loaded.get("splits", {})
    ev = loaded.get("evaluation", {})
    sub = loaded.get("subsample", {})
    rnk = loaded.get("ranker", {})
    arms_raw = loaded.get("arms", {})
    out = loaded.get("output", {})

    registry_v1 = RegistryReleaseConfig(
        release_name="structcover_v1",
        path=Path(str(reg["structcover_v1"]["path"])),
        sha256=str(reg["structcover_v1"]["sha256"]).upper(),
        records=int(reg["structcover_v1"]["records"]),
    )
    registry_v2 = RegistryReleaseConfig(
        release_name="structcover_v2",
        path=Path(str(reg["structcover_v2"]["path"])),
        sha256=str(reg["structcover_v2"]["sha256"]).upper(),
        records=int(reg["structcover_v2"]["records"]),
    )

    return StructureCoverageExperimentConfig(
        version=2,
        experiment_name=str(exp.get("name", "")),
        benchmark_path=str(exp["benchmark_labels"]["path"]),
        benchmark_sha256=str(exp["benchmark_labels"]["sha256"]).upper(),
        proteome_path=str(exp["reference_proteome"]["path"]),
        proteome_sha256=str(exp["reference_proteome"]["sha256"]).upper(),
        clusters_path=str(exp["protein_clusters"]["path"]),
        clusters_sha256=str(exp["protein_clusters"]["sha256"]).upper(),
        legacy_config_path=str(exp["legacy_config"]["path"]),
        legacy_config_sha256=str(exp["legacy_config"]["sha256"]).upper(),
        registry_v1=registry_v1,
        registry_v2=registry_v2,
        studies=[str(s) for s in splits.get("studies", [])],
        seeds=[int(s) for s in ev.get("seeds", [])],
        cluster_bootstrap_seed=int(ev.get("cluster_bootstrap_seed", 1729)),
        cluster_bootstrap_replicates=int(ev.get("cluster_bootstrap_replicates", 5000)),
        permutation_seed=int(ev.get("permutation_seed", 2718)),
        permutation_replicates=int(ev.get("permutation_replicates", 5000)),
        subsample_ratio=int(sub.get("unlabeled_per_positive", 20)),
        subsample_seed=int(sub.get("seed", 12345)),
        ranker_hidden=int(rnk.get("hidden", 16)),
        ranker_dropout=float(rnk.get("dropout", 0.2)),
        ranker_epochs=int(rnk.get("epochs", 200)),
        ranker_lr=float(rnk.get("lr", 0.05)),
        ranker_n_mc_dropout=int(rnk.get("n_mc_dropout", 16)),
        arms={
            str(k): {str(ak): bool(av) for ak, av in dict(v).items()}
            for k, v in arms_raw.items()
        },
        output_directory=str(out.get("directory", "")),
    )


def validate_controlled_variables(
    cfg: StructureCoverageExperimentConfig,
) -> None:
    """Fail-closed validation of every variable that must be identical across
    the two registry releases.

    Raises ``RuntimeError`` on the first violation.
    """
    if set(cfg.arms) != ALLOWED_ARMS:
        extra = set(cfg.arms) - ALLOWED_ARMS
        missing = ALLOWED_ARMS - set(cfg.arms)
        parts = []
        if extra:
            parts.append(f"unknown arms: {sorted(extra)}")
        if missing:
            parts.append(f"missing arms: {sorted(missing)}")
        raise RuntimeError("arm mismatch: " + "; ".join(parts))

    if cfg.studies != ["PXD006140", "PXD024061"]:
        raise RuntimeError(f"studies must be [PXD006140, PXD024061], got {cfg.studies}")

    if cfg.seeds != [0, 1, 2, 3, 4]:
        raise RuntimeError(f"seeds must be [0,1,2,3,4], got {cfg.seeds}")

    if len(set(cfg.seeds)) != len(cfg.seeds):
        raise RuntimeError("seeds must not contain duplicates")

    if cfg.subsample_ratio != 20:
        raise RuntimeError(f"subsample ratio must be 20, got {cfg.subsample_ratio}")

    if cfg.subsample_seed != 12345:
        raise RuntimeError(f"subsample seed must be 12345, got {cfg.subsample_seed}")

    if cfg.ranker_hidden != 16:
        raise RuntimeError(f"ranker hidden must be 16, got {cfg.ranker_hidden}")

    if cfg.ranker_epochs != 200:
        raise RuntimeError(f"ranker epochs must be 200, got {cfg.ranker_epochs}")

    if cfg.ranker_lr != 0.05:
        raise RuntimeError(f"ranker lr must be 0.05, got {cfg.ranker_lr}")

    if cfg.ranker_n_mc_dropout != 16:
        raise RuntimeError(f"n_mc_dropout must be 16, got {cfg.ranker_n_mc_dropout}")

    for arm_name, arm in cfg.arms.items():
        if arm.get("use_esm", True):
            raise RuntimeError(f"arm {arm_name} must have use_esm=false")
        if arm.get("use_study_context", True):
            raise RuntimeError(f"arm {arm_name} must have use_study_context=false")

    for forbidden in FORBIDDEN_OUTPUT_ROOTS:
        if forbidden in cfg.output_directory.lower():
            raise RuntimeError(
                f"output directory must not overlap with {forbidden}: "
                f"{cfg.output_directory}"
            )

    if cfg.registry_v1.records != 7:
        raise RuntimeError(f"v1 records must be 7, got {cfg.registry_v1.records}")
    if cfg.registry_v2.records != 2006:
        raise RuntimeError(f"v2 records must be 2006, got {cfg.registry_v2.records}")

    if cfg.cluster_bootstrap_seed != 1729:
        raise RuntimeError("cluster bootstrap seed must be 1729")
    if cfg.cluster_bootstrap_replicates != 5000:
        raise RuntimeError("cluster bootstrap replicates must be 5000")

    if cfg.permutation_replicates <= 0:
        raise RuntimeError("permutation replicates must be positive")
    if cfg.cluster_bootstrap_replicates <= 0:
        raise RuntimeError("cluster bootstrap replicates must be positive")
