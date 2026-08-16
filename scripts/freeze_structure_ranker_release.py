#!/usr/bin/env python
"""Gate 0 — freeze the v11 structure_ranker candidate release bundle.

Fits the frozen ``structure_ranker`` configuration of
``configs/experiments/multispecies_v2_global_clusters_v11.yaml`` on the
complete v11 literature-random-track development panel — the union of all ten
seed panels (train + validation + test partitions, one row per cysteine site)
— with the release seed **20260813** (chosen at freeze time, never selected
by performance), and persists the ``StructureRankerBundle`` into the
candidate release package's ``model_weights/`` directory together with a
content-addressed fit manifest.

Guards (all enforced or recorded here):

- config SHA256 must match the value pinned in the release manifest;
- rows come exclusively from the development split (the random-track builder
  raises if a frozen test row enters the panel);
- the fit uses the frozen hyperparameters and ablation policy verbatim;
- the resulting bundle SHA256, panel SHA256 and code revision are recorded in
  ``fit_manifest.json`` so blind-time verification (Gate 3) can confirm the
  frozen artifact is byte-identical to this fit.

After Gate 0, blind-time scoring uses ONLY the saved bundle via
``inference_script.py``; no refit ever runs on blind data.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

# Allow ``python scripts/freeze_structure_ranker_release.py`` from the repo
# root: ``scripts`` is a package whose sibling modules we reuse.
_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import yaml

from plantpersulf.benchmark.literature_random_track import (
    build_literature_random_track,
)
from plantpersulf.evaluation.comparable_track import (
    build_registered_structure_features,
)
from plantpersulf.models.structure_ranker import (
    AblationConfig,
    BranchFeatures,
    fit_structure_ranker,
)
from plantpersulf.proteomics.multispecies_v2_sources import sequence_feature_map
from scripts.train_multispecies_v2 import (
    _code_revision,
    _dirty_paths,
    _development_rows,
)

RELEASE_ID = "multispecies-v2-candidate-release-v1"
RELEASE_SEED = 20260813


def _sha256_bytes(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _site_key(row: object) -> tuple[str, int]:
    return row.global_protein_id, row.cys_position  # type: ignore[attr-defined]


def _panel_sha256(rows: list) -> str:
    lines = []
    for row in rows:
        lines.append(
            "\t".join(
                (
                    row.species,  # type: ignore[attr-defined]
                    row.protein_accession,  # type: ignore[attr-defined]
                    str(row.cys_position),  # type: ignore[attr-defined]
                    row.label,  # type: ignore[attr-defined]
                    ";".join(row.study_accessions),  # type: ignore[attr-defined]
                )
            )
        )
    return hashlib.sha256("\n".join(lines).encode("utf-8")).hexdigest()


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("configs/experiments/multispecies_v2_global_clusters_v11.yaml"),
    )
    parser.add_argument(
        "--release-dir",
        type=Path,
        default=Path("results/candidates/multispecies_v2_candidate_release_v1"),
    )
    parser.add_argument("--seed", type=int, default=RELEASE_SEED)
    parser.add_argument(
        "--device", type=str, default=None, help="override PLANTPERSULF_DEVICE"
    )
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_argument_parser().parse_args(argv)

    release_dir = args.release_dir
    manifest_path = release_dir / "manifest.json"
    if not manifest_path.is_file():
        raise RuntimeError(f"release manifest not found: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if not isinstance(manifest, dict):
        raise RuntimeError("release manifest is corrupt")
    pinned_config_sha256 = manifest.get("source_config_sha256")
    if not isinstance(pinned_config_sha256, str):
        raise RuntimeError("release manifest is missing source_config_sha256")

    config_sha256 = _sha256_bytes(args.config)
    if config_sha256 != pinned_config_sha256:
        raise RuntimeError(
            "config SHA256 does not match the frozen release manifest: "
            f"{config_sha256} != {pinned_config_sha256}"
        )
    print(f"config {args.config}: sha256={config_sha256} (matches frozen manifest)")

    # --- panel: union of all ten literature-track seed panels ----------------
    rows, _frozen, proteomes = _development_rows(args.config, sample_unlabeled=False)
    cfg = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    if not isinstance(cfg, dict):
        raise RuntimeError("multispecies v2 config is required")
    random_cfg = cfg["literature_random_protein"]
    seeds = tuple(random_cfg["seeds"])
    if seeds != tuple(range(10)):
        raise RuntimeError("literature random track requires seeds 0-9")
    runs = build_literature_random_track(
        rows,
        seeds=seeds,
        test_fraction=random_cfg["test_fraction"],
        validation_fraction_of_remaining=random_cfg["validation_fraction_of_remaining"],
        unlabeled_per_positive=cfg["unlabeled_panel"]["per_positive"],
        sampling_seed=cfg["unlabeled_panel"]["seed"],
    )
    panel: dict[tuple[str, int], object] = {}
    for run in runs:
        for partition in run.partitions.values():
            for row in partition:
                panel.setdefault(_site_key(row), row)
    panel_rows = [panel[key] for key in sorted(panel)]
    n_positives = sum(1 for row in panel_rows if row.label == "positive")  # type: ignore[attr-defined]
    panel_sha256 = _panel_sha256(panel_rows)
    print(
        f"release fit panel: n_rows={len(panel_rows)} n_positives={n_positives} "
        f"sha256={panel_sha256}"
    )

    # --- features: sequence + structure, ESM disabled (use_esm=0) ------------
    raw_features = sequence_feature_map(panel_rows, proteomes)
    if any(any(value is None for value in values) for values in raw_features.values()):
        raise RuntimeError("sequence comparison features contain missing values")
    sequence_features = {
        key: tuple(float(value) for value in values if value is not None)
        for key, values in raw_features.items()
    }
    inputs_cfg = cfg.get("comparison_inputs", {})
    if not isinstance(inputs_cfg, dict):
        raise RuntimeError("comparison_inputs must be a mapping")
    structure_registry = Path(str(inputs_cfg["structure_features"]))
    build_cfg = cfg.get("comparison_feature_build")
    if not isinstance(build_cfg, dict):
        raise RuntimeError("comparison_feature_build configuration is required")
    structure_features = build_registered_structure_features(
        list(sequence_features),
        registry_path=structure_registry,
        registry_base=Path(str(build_cfg["structure_registry_base"])),
    )
    missing_structure = sum(
        1 for _, available in structure_features.values() if not available
    )
    print(
        f"structure registry: sites={len(structure_features)} "
        f"missing_structures={missing_structure}"
    )

    train = BranchFeatures(
        sequence=[list(sequence_features[key]) for key in sorted(sequence_features)],
        esm=[[0.0] for _ in panel_rows],
        structure=[
            list(values) for values, _ in (
                structure_features[key] for key in sorted(sequence_features)
            )
        ],
        structure_mask=[
            available
            for _, available in (
                structure_features[key] for key in sorted(sequence_features)
            )
        ],
        study_ids=None,
    )
    train_y = [str(row.label) for row in panel_rows]  # type: ignore[attr-defined]

    # --- frozen fit -----------------------------------------------------------
    params = cfg["literature_baseline_parameters"]["structure_ranker"]
    ablation = AblationConfig(use_esm=False, use_study_context=False)
    print(
        f"fit: seed={args.seed} epochs={params['epochs']} hidden={params['hidden']} "
        f"dropout={params['dropout']} lr={params['lr']} "
        f"holdout_fraction={params['holdout_fraction']} "
        f"n_mc_dropout={params['n_mc_dropout']} device={args.device or 'cpu'}"
    )
    bundle = fit_structure_ranker(
        train,
        train_y,
        seed=args.seed,
        ablation=ablation,
        hidden=int(params["hidden"]),
        dropout=float(params["dropout"]),
        epochs=int(params["epochs"]),
        lr=float(params["lr"]),
        holdout_fraction=float(params["holdout_fraction"]),
        n_mc_dropout=int(params["n_mc_dropout"]),
        device_name=args.device,
        batch_size=None,
    )

    # --- persist bundle + fit manifest ---------------------------------------
    weights_dir = release_dir / "model_weights"
    bundle_path = weights_dir / "structure_ranker_bundle.pt"
    bundle.save(bundle_path)
    bundle_sha256 = _sha256_bytes(bundle_path)
    dirty = _dirty_paths()
    fit_manifest = {
        "release_id": manifest.get("release_id"),
        "artifact": "model_weights/structure_ranker_bundle.pt",
        "model": "structure_ranker",
        "source_config_path": args.config.as_posix(),
        "source_config_sha256": config_sha256,
        "training_panel": {
            "definition": (
                "union of all ten v11 literature-random-track seed panels "
                "(train+validation+test partitions), one row per cysteine site; "
                "development split only"
            ),
            "n_rows": len(panel_rows),
            "n_positives": n_positives,
            "sha256": panel_sha256,
        },
        "seed": args.seed,
        "ablation": {"use_esm": False, "use_study_context": False},
        "hyperparameters": {
            "epochs": int(params["epochs"]),
            "hidden": int(params["hidden"]),
            "dropout": float(params["dropout"]),
            "lr": float(params["lr"]),
            "holdout_fraction": float(params["holdout_fraction"]),
            "n_mc_dropout": int(params["n_mc_dropout"]),
            "use_esm": 0,
        },
        "device": args.device or "cpu",
        "code_revision": _code_revision(),
        "dirty_paths": list(dirty),
        "bundle_sha256": bundle_sha256,
        "created": datetime.now(timezone.utc).isoformat(),
    }
    fit_manifest_path = release_dir / "model_weights" / "fit_manifest.json"
    fit_manifest_path.write_text(
        json.dumps(fit_manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(f"bundle: {bundle_path} sha256={bundle_sha256}")
    print(f"fit manifest: {fit_manifest_path}")


if __name__ == "__main__":
    main()
