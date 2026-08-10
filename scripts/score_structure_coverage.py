#!/usr/bin/env python
"""Phase S0, Task 5 — paired structure-coverage scoring runner.

Scores the frozen benchmark rows under two AlphaFold structure-registry
releases (7 vs 2,006 accessions), five ESM-free ablation arms each, and
identical folds/seeds/hyperparameters. The registry is selected explicitly
per release so the comparison cannot drift when the mutable registry grows.

Usage::

    python scripts/score_structure_coverage.py --verify-only
    python scripts/score_structure_coverage.py \
        --config configs/experiments/pu_ranker_structcover_v2.yaml
"""

from __future__ import annotations

import argparse
import csv
import json
import random
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Any

from plantpersulf.evaluation.structure_coverage_audit import (
    FrozenFile,
    audit_structure_coverage,
    sha256_file,
    verify_frozen_file,
    verify_registry_pair,
)
from plantpersulf.evaluation.structure_coverage_config import (
    StructureCoverageExperimentConfig,
    load_structure_coverage_config,
    validate_controlled_variables,
)
from plantpersulf.models.structure_ranker import AblationConfig
from scripts import run_experiment

SCORE_FIELDS = (
    "coverage_release",
    "held_out_study",
    "seed",
    "arm",
    "protein_accession",
    "cys_position_in_protein",
    "label",
    "study_accession",
    "cluster_id",
    "has_registered_structure",
    "mapping_status",
    "plddt_bin",
    "score",
    "uncertainty",
)


def prepare_output_directory(path: Path) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {path}")
    path.mkdir(parents=True)


def _ablation_config_for_arm(arm_cfg: dict[str, bool]) -> AblationConfig:
    return AblationConfig(
        use_esm=bool(arm_cfg.get("use_esm", True)),
        use_structure=bool(arm_cfg.get("use_structure", True)),
        use_plddt=bool(arm_cfg.get("use_plddt", True)),
        use_accessibility=bool(arm_cfg.get("use_accessibility", True)),
        use_study_context=bool(arm_cfg.get("use_study_context", True)),
    )


def _plddt_bin(plddt: float | None) -> str:
    if plddt is None:
        return "none"
    if plddt < 50.0:
        return "lt_50"
    if plddt < 70.0:
        return "50_70"
    if plddt < 90.0:
        return "70_90"
    return "ge_90"


def run_structure_coverage_scoring(
    config_path: Path,
    *,
    verify_only: bool = False,
) -> Path:
    """Run the paired scoring experiment once. Returns the output directory.

    Required order:
    1. load and validate frozen config
    2. verify all frozen hashes
    3. verify the v1/v2 registry relationship
    4. audit both releases and call require_pass()
    5. in verify_only, emit no directory and stop
    6. ensure output is absent and create it once
    7. build frozen folds/subsamples once
    8. loop release -> fold -> seed -> arm
    9. write sorted per-release TSVs
    """
    # --- Step 1: load and validate ---
    cfg = load_structure_coverage_config(config_path)
    validate_controlled_variables(cfg)

    # --- Step 2: verify frozen hashes ---
    _verify_frozen_inputs(cfg, audits=None)

    # --- Step 3-4: verify registry pair and audit both releases ---
    audits: list[Any] = []
    _verify_frozen_inputs(cfg, audits=audits)
    for a in audits:
        a.require_pass()

    # --- Step 5: verify-only mode ---
    out_dir = Path(cfg.output_directory)
    if verify_only:
        print(
            f"verify_only: all frozen files audit OK, "
            f"training_started=false, output_dir={out_dir}"
        )
        return out_dir

    # --- Step 6: ensure output absent and create ---
    prepare_output_directory(out_dir)

    # Temporary sibling directory for atomic writes
    tmp_dir = out_dir.with_name(out_dir.name + ".tmp")
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)
    tmp_dir.mkdir(parents=True)

    try:
        _run_paired_scoring(cfg, audits, tmp_dir)

        # Atomic rename on success
        for child in tmp_dir.iterdir():
            shutil.move(str(child), str(out_dir / child.name))
        tmp_dir.rmdir()

        # Write manifest last
        _write_manifest(cfg, out_dir)
    except Exception:
        # Retain failure record outside the result directory
        fail_path = out_dir.with_name(out_dir.name + ".FAILED")
        fail_path.write_text(
            f"Scoring run failed. Partial output at {out_dir} (if it was created).\n",
            encoding="utf-8",
        )
        raise

    return out_dir


def _verify_frozen_inputs(
    cfg: StructureCoverageExperimentConfig,
    *,
    audits: list[Any] | None,
) -> None:
    """Verify every input hash and audit both releases if audits list provided."""
    verify_frozen_file(
        FrozenFile(
            name="benchmark",
            path=Path(cfg.benchmark_path),
            sha256=cfg.benchmark_sha256,
        )
    )
    verify_frozen_file(
        FrozenFile(
            name="proteome",
            path=Path(cfg.proteome_path),
            sha256=cfg.proteome_sha256,
        )
    )
    verify_frozen_file(
        FrozenFile(
            name="clusters",
            path=Path(cfg.clusters_path),
            sha256=cfg.clusters_sha256,
        )
    )
    verify_frozen_file(
        FrozenFile(
            name="registry_v1",
            path=cfg.registry_v1.path,
            sha256=cfg.registry_v1.sha256,
        )
    )
    verify_frozen_file(
        FrozenFile(
            name="registry_v2",
            path=cfg.registry_v2.path,
            sha256=cfg.registry_v2.sha256,
        )
    )
    verify_frozen_file(
        FrozenFile(
            name="legacy_config",
            path=Path(cfg.legacy_config_path),
            sha256=cfg.legacy_config_sha256,
        )
    )

    verify_registry_pair(cfg.registry_v1.path, cfg.registry_v2.path)

    REGISTRY_BASE = Path("data/registry")

    if audits is not None:
        for release_name, reg, _reg_records in [
            ("structcover_v1", cfg.registry_v1, cfg.registry_v1.records),
            ("structcover_v2", cfg.registry_v2, cfg.registry_v2.records),
        ]:
            audit = audit_structure_coverage(
                release=release_name,
                benchmark_sha256=cfg.benchmark_sha256,
                proteome_sha256=cfg.proteome_sha256,
                clusters_sha256=cfg.clusters_sha256,
                registry_sha256=reg.sha256,
                registry_path=reg.path,
                registry_base=REGISTRY_BASE,
                benchmark_path=Path(cfg.benchmark_path),
                proteome_path=Path(cfg.proteome_path),
                clusters_path=Path(cfg.clusters_path),
            )
            audits.append(audit)
            print(
                f"audit {release_name}: "
                f"{len(audit.records)} rows, "
                f"mapped={sum(1 for r in audit.records if r.mapping_status == 'mapped_cys')}, "
                f"blocking={len(audit.blocking_errors)}"
            )


def _run_paired_scoring(
    cfg: StructureCoverageExperimentConfig,
    audits: list[Any],
    tmp_dir: Path,
) -> None:
    """Run paired scoring for both releases with identical rows/splits."""
    # Read benchmark once
    all_rows = _read_benchmark(Path(cfg.benchmark_path))
    proteome_path = Path(cfg.proteome_path)
    clusters_path = Path(cfg.clusters_path)
    ratio = cfg.subsample_ratio
    sub_seed = cfg.subsample_seed
    seeds = cfg.seeds
    ranker_params = {
        "hidden": cfg.ranker_hidden,
        "dropout": cfg.ranker_dropout,
        "epochs": cfg.ranker_epochs,
        "lr": cfg.ranker_lr,
        "n_mc_dropout": cfg.ranker_n_mc_dropout,
    }
    need_esm = False  # All five arms are ESM-free

    # Cluster mapping for output rows
    cluster_map = _read_clusters(clusters_path)

    # Audit records lookup: per-release (accession, position) -> audit record
    audit_lookup: dict[str, dict[tuple[str, int], Any]] = {}
    for audit in audits:
        lookup: dict[tuple[str, int], Any] = {}
        for rec in audit.records:
            lookup[(rec.protein_accession, rec.cys_position_in_protein)] = rec
        audit_lookup[audit.records[0].release if audit.records else ""] = lookup

    # Two releases in order
    releases = [
        ("structcover_v1", cfg.registry_v1),
        ("structcover_v2", cfg.registry_v2),
    ]
    REGISTRY_BASE = Path("data/registry")

    for rel_name, rel in releases:
        print(f"\n===== Release: {rel_name} =====")
        scored_rows: list[dict[str, str]] = []

        for study in cfg.studies:
            print(f"  Fold: leave_{study}_out")
            train_rows, test_rows = run_experiment._study_fold_rows(all_rows, study)
            train_rows = run_experiment._subsample_unlabeled(
                train_rows, ratio, sub_seed
            )
            test_rows = run_experiment._subsample_unlabeled(
                test_rows, ratio, sub_seed + 1
            )

            # Identical train/val split as run_experiment.py
            rng = random.Random(sub_seed)
            train_dedup = train_rows[:]
            rng.shuffle(train_dedup)
            n_val = max(1, int(len(train_dedup) * 0.2))
            train_rows_fold = train_dedup[n_val:]
            train_y = [r["label"] for r in train_rows_fold]

            scratch = Path(tempfile.mkdtemp(prefix="sc_structcover_"))
            try:
                branch_train = run_experiment._build_branch_features(
                    train_rows_fold,
                    proteome_path,
                    scratch,
                    f"{rel_name}_{study}_train",
                    need_esm,
                    structure_registry_path=rel.path,
                    structure_registry_base=REGISTRY_BASE,
                )
                branch_test = run_experiment._build_branch_features(
                    test_rows,
                    proteome_path,
                    scratch,
                    f"{rel_name}_{study}_test",
                    need_esm,
                    structure_registry_path=rel.path,
                    structure_registry_base=REGISTRY_BASE,
                )
            finally:
                shutil.rmtree(scratch, ignore_errors=True)

            for arm_name, arm_cfg in cfg.arms.items():
                ab = _ablation_config_for_arm(arm_cfg)
                for seed in seeds:
                    from plantpersulf.models.structure_ranker import (
                        structure_ranker_scores,
                    )

                    out = structure_ranker_scores(
                        branch_train,
                        train_y,
                        branch_test,
                        seed=seed,
                        ablation=ab,
                        hidden=int(ranker_params["hidden"]),
                        dropout=float(ranker_params["dropout"]),
                        epochs=int(ranker_params["epochs"]),
                        lr=float(ranker_params["lr"]),
                        n_mc_dropout=int(ranker_params["n_mc_dropout"]),
                    )

                    lookup = audit_lookup.get(rel_name, {})

                    for row, score, unc in zip(
                        test_rows, out.scores, out.uncertainty, strict=True
                    ):
                        acc = row["protein_accession"]
                        pos = int(row["cys_position_in_protein"])
                        audit_rec = lookup.get((acc, pos))
                        cluster = cluster_map.get(acc, f"__singleton__{acc}")

                        has_struct = "0"
                        mapping_status = "absent_structure"
                        plddt_val = None
                        if audit_rec is not None:
                            has_struct = (
                                "1" if audit_rec.has_registered_structure else "0"
                            )
                            mapping_status = audit_rec.mapping_status
                            plddt_val = audit_rec.plddt

                        scored_rows.append(
                            {
                                "coverage_release": rel_name,
                                "held_out_study": study,
                                "seed": str(seed),
                                "arm": arm_name,
                                "protein_accession": acc,
                                "cys_position_in_protein": str(pos),
                                "label": row["label"],
                                "study_accession": row.get("study_accession", ""),
                                "cluster_id": cluster,
                                "has_registered_structure": has_struct,
                                "mapping_status": mapping_status,
                                "plddt_bin": _plddt_bin(plddt_val),
                                "score": f"{score:.10g}",
                                "uncertainty": f"{unc:.10g}",
                            }
                        )

                    from plantpersulf.evaluation.metrics import average_precision

                    ap = average_precision(
                        list(
                            zip(
                                out.scores, [r["label"] for r in test_rows], strict=True
                            )
                        )
                    )
                    print(
                        f"    {rel_name} {study} {arm_name} "
                        f"seed={seed} test_ap={ap:.4f}"
                    )

        # Write per-release TSV
        scored_rows.sort(
            key=lambda r: (
                r["held_out_study"],
                r["seed"],
                r["arm"],
                r["protein_accession"],
                r["cys_position_in_protein"],
            )
        )
        _write_scored(tmp_dir / f"scores_{rel_name}.tsv", scored_rows)
        print(f"  wrote {len(scored_rows)} rows -> scores_{rel_name}.tsv")


def _read_benchmark(path: Path) -> list[dict[str, str]]:
    BENCHMARK_FIELDS = (
        "protein_accession",
        "cys_position_in_protein",
        "label",
        "study_accession",
        "evidence_level",
        "source_sha256",
    )
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != BENCHMARK_FIELDS:
            raise RuntimeError(f"benchmark labels have invalid columns: {path}")
        return [dict(r) for r in reader]


def _read_clusters(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            mapping[row["protein_accession"]] = row["cluster_id"]
    return mapping


def _write_scored(
    path: Path,
    rows: list[dict[str, str]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(SCORE_FIELDS),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def _write_manifest(
    cfg: StructureCoverageExperimentConfig,
    out_dir: Path,
) -> None:
    """Write a manifest recording every output file's SHA256."""
    manifest: dict[str, object] = {
        "experiment": cfg.experiment_name,
        "config_sha256": sha256_file(
            Path("configs/experiments/pu_ranker_structcover_v2.yaml")
        ),
        "output_files": {},
    }
    for child in sorted(out_dir.iterdir()):
        if child.is_file():
            manifest["output_files"] = {
                **manifest.get("output_files", {}),
                child.name: sha256_file(child),
            }
    (out_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(f"manifest written to {out_dir / 'manifest.json'}")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="paired structure-coverage scoring runner")
    p.add_argument(
        "--config",
        type=Path,
        default=Path("configs/experiments/pu_ranker_structcover_v2.yaml"),
    )
    p.add_argument("--verify-only", action="store_true")
    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    out = run_structure_coverage_scoring(args.config, verify_only=args.verify_only)
    if not args.verify_only:
        print(f"\nDone. Output at {out}")
    sys.exit(0)
