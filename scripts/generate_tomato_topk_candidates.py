#!/usr/bin/env python
"""Gate 0 — score the tomato proteome with the frozen release bundle.

Scores every cysteine site of the SHA-pinned tomato reference proteome
with the frozen ``structure_ranker`` release bundle, EXCLUDING every site
that belongs to the release training panel (the multispecies v2
development sites the model was trained on — scoring them would be
training-contaminated, not a candidate ranking). No blind data is involved:
the tomato reference proteome is public, and the lockbox cohort (Gate 3)
has not been generated.

Deterministic pipeline (all hashes recorded in the generation manifest):

1. verify the v11 config against the frozen release manifest SHA256;
2. rebuild the release training panel exactly as the freeze fit did and
   verify its panel SHA256 against ``fit_manifest.json``;
3. build every tomato Cys site from the SHA-pinned proteome; drop panel
   sites;
4. extract sequence features (same extractor as the v11 track) and
   registered AlphaFold structure features with explicit missingness
   mask;
5. score with the saved bundle via ``score_structure_ranker_bundle``
   (CPU by default; no refit — only the frozen state dict is used);
6. write the ranked table + matched controls + a content-addressed
   generation manifest.

Matched-controls rule (frozen here; mirrored in the SAP protocol): for
each candidate site, controls are drawn from non-panel tomato sites in the
same structure-availability stratum (has_structure true/false), ranked by
Euclidean distance in z-scored sequence-feature space (z-scores from the
candidate pool itself), with a caliper of 0.25 standard deviations; at
most 5 controls per candidate, no site used twice as a control.

Outputs:

- ``top_k_candidates.tsv`` — full ranked table (site_key, cys_position,
  score, uncertainty), descending score; the SAP protocol fixes K.
- ``matched_controls.tsv`` — candidate/control pairs.
- ``generation_manifest.json`` — input and bundle hashes, counts.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

# Allow running directly from the repo root.
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
from plantpersulf.features.sequence import _load_proteome
from plantpersulf.models.structure_ranker import (
    BranchFeatures,
    StructureRankerBundle,
    score_structure_ranker_bundle,
)
from plantpersulf.proteomics.multispecies_v2_dataset import MultispeciesV2SiteRow
from plantpersulf.proteomics.multispecies_v2_sources import sequence_feature_map
from scripts.train_multispecies_v2 import _development_rows

RELEASE_ID = "multispecies-v2-candidate-release-v1"
CONTROLS_PER_CANDIDATE = 5
CONTROL_CALIPER_SD = 0.25
CANDIDATE_RANK_HEAD = 200  # controls are matched for the top CANDIDATE_RANK_HEAD


def _sha256_bytes(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _panel_sha256(rows: list) -> str:
    lines = []
    for row in rows:
        lines.append(
            "\t".join(
                (
                    row.species,
                    row.protein_accession,
                    str(row.cys_position),
                    row.label,
                    ";".join(row.study_accessions),
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
    parser.add_argument(
        "--species",
        choices=("tomato", "arabidopsis", "rice"),
        default="tomato",
        help="scan species; non-tomato runs are diagnostics (output to --output-dir)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="default: release dir for tomato; results/diagnostics/candidates_<species> otherwise",
    )
    parser.add_argument(
        "--device", type=str, default=None, help="default CPU; PLANTPERSULF_DEVICE honored"
    )
    parser.add_argument("--batch-size", type=int, default=16384)
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

    config_sha256 = _sha256_bytes(args.config)
    if config_sha256 != manifest.get("source_config_sha256"):
        raise RuntimeError(
            "config SHA256 does not match the frozen release manifest: "
            f"{config_sha256} != {manifest.get('source_config_sha256')}"
        )
    cfg = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    if not isinstance(cfg, dict):
        raise RuntimeError("multispecies v2 config is required")

    # --- release training panel (verify byte-identical to the freeze fit) ----
    rows, _frozen, _proteomes = _development_rows(args.config, sample_unlabeled=False)
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
                panel.setdefault((row.global_protein_id, row.cys_position), row)
    panel_rows = [panel[key] for key in sorted(panel)]
    panel_sha256 = _panel_sha256(panel_rows)
    fit_manifest = json.loads(
        (release_dir / "model_weights" / "fit_manifest.json").read_text(encoding="utf-8")
    )
    pinned_panel_sha256 = fit_manifest["training_panel"]["sha256"]
    if panel_sha256 != pinned_panel_sha256:
        raise RuntimeError(
            f"rebuilt panel does not match the frozen fit panel: "
            f"{panel_sha256} != {pinned_panel_sha256}"
        )
    print(f"release training panel verified: {len(panel_rows)} rows, sha256={panel_sha256}")
    species_panel_keys = {
        (row.global_protein_id, row.cys_position)
        for row in panel_rows
        if row.species == args.species
    }
    print(
        f"{args.species} panel sites excluded from candidates: "
        f"{len(species_panel_keys)}"
    )

    # --- species proteome + all Cys sites --------------------------------------
    species_ref = next(
        item for item in cfg["reference_proteomes"] if item["species"] == args.species
    )
    species_path = Path(str(species_ref["path"]))
    if _sha256_bytes(species_path) != species_ref["sha256"]:
        raise RuntimeError(f"{args.species} reference proteome SHA256 mismatch")
    proteome = _load_proteome(species_path)
    scan_rows: list[MultispeciesV2SiteRow] = []
    for accession, sequence in sorted(proteome.items()):
        global_id = f"{args.species}|{accession}"
        if not sequence:
            continue
        for position, residue in enumerate(sequence, start=1):
            if residue != "C":
                continue
            key = (global_id, position)
            if key in species_panel_keys:
                continue
            scan_rows.append(
                MultispeciesV2SiteRow(
                    species=args.species,
                    protein_accession=accession,
                    cys_position=position,
                    label="unlabeled",
                    study_accessions=(),
                    global_protein_id=global_id,
                    cluster_id="",
                    split="development",
                    development_fold=None,
                )
            )
    print(f"{args.species} scan sites (panel-excluded): {len(scan_rows)}")
    if not scan_rows:
        raise RuntimeError(f"{args.species} scan produced no candidate sites")

    # --- features -------------------------------------------------------------
    raw_features = sequence_feature_map(scan_rows, {args.species: proteome})
    if any(any(value is None for value in values) for values in raw_features.values()):
        raise RuntimeError("sequence comparison features contain missing values")
    sequence_features = {
        key: tuple(float(value) for value in values if value is not None)
        for key, values in raw_features.items()
    }
    structure_features = build_registered_structure_features(
        set(sequence_features),
        registry_path=Path(str(cfg["comparison_inputs"]["structure_features"])),
        registry_base=Path(str(cfg["comparison_feature_build"]["structure_registry_base"])),
    )
    n_with_structure = sum(1 for _, available in structure_features.values() if available)
    print(
        f"structure coverage: {n_with_structure}/{len(scan_rows)} "
        f"({n_with_structure / len(scan_rows):.2%})"
    )

    ordered_keys = sorted(sequence_features)
    train = BranchFeatures(
        sequence=[list(sequence_features[key]) for key in ordered_keys],
        esm=[[0.0] for _ in ordered_keys],
        structure=[list(structure_features[key][0]) for key in ordered_keys],
        structure_mask=[structure_features[key][1] for key in ordered_keys],
        study_ids=None,
    )

    # --- score with the frozen bundle -----------------------------------------
    bundle_path = release_dir / "model_weights" / "structure_ranker_bundle.pt"
    bundle = StructureRankerBundle.load(bundle_path)
    actual_bundle_sha = _sha256_bytes(bundle_path)
    if actual_bundle_sha != fit_manifest["bundle_sha256"]:
        raise RuntimeError(
            "bundle SHA256 mismatch: bundle changed since the frozen fit"
        )
    print(
        f"scoring {len(ordered_keys)} sites (batch {args.batch_size}, "
        f"device {args.device or 'cpu'})..."
    )
    output = score_structure_ranker_bundle(
        bundle, train, device_name=args.device, batch_size=args.batch_size
    )
    ranked = sorted(
        zip(ordered_keys, output.scores, output.uncertainty, strict=True),
        key=lambda item: item[1],
        reverse=True,
    )
    print(
        f"score range: [{min(s for _, s, _ in ranked):.6f}, "
        f"{max(s for _, s, _ in ranked):.6f}]"
    )

    # --- outputs --------------------------------------------------------------
    if args.species == "tomato" and args.output_dir is None:
        output_dir = release_dir  # frozen package artifact (unchanged behaviour)
    else:
        output_dir = args.output_dir or (
            _REPO_ROOT / "results" / "diagnostics" / f"candidates_{args.species}"
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"output directory: {output_dir}")
    candidates_path = output_dir / "top_k_candidates.tsv"
    with candidates_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("site_key", "cys_position", "score", "uncertainty"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for (global_id, position), score, uncertainty in ranked:
            writer.writerow(
                {
                    "site_key": global_id,
                    "cys_position": position,
                    "score": f"{score:.10g}",
                    "uncertainty": f"{uncertainty:.10g}",
                }
            )

    # --- matched controls (rule documented in the module docstring) -----------
    candidate_head = ranked[:CANDIDATE_RANK_HEAD]
    head_keys = {key for key, _, _ in candidate_head}
    control_pool = [key for key in ordered_keys if key not in head_keys]
    z = _z_score_features(sequence_features, control_pool)
    stratum: dict[bool, list[tuple[str, int]]] = defaultdict(list)
    for key in control_pool:
        stratum[structure_features[key][1]].append(key)
    pairs: list[tuple[tuple[str, int], tuple[str, int], float]] = []
    used_controls: set[tuple[str, int]] = set()
    for key, _score, _unc in candidate_head:
        features = z[key]
        candidates_in_stratum = [
            control
            for control in stratum[structure_features[key][1]]
            if control not in used_controls
        ]
        distances = sorted(
            (
                (control, sum((a - b) ** 2 for a, b in zip(features, z[control], strict=True)) ** 0.5)
                for control in candidates_in_stratum
            ),
            key=lambda item: item[1],
        )
        for control, distance in distances[:CONTROLS_PER_CANDIDATE]:
            if distance > CONTROL_CALIPER_SD:
                continue
            pairs.append((key, control, distance))
            used_controls.add(control)
    controls_path = output_dir / "matched_controls.tsv"
    with controls_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "candidate_site_key",
                "candidate_cys_position",
                "control_site_key",
                "control_cys_position",
                "sequence_distance_sd",
            ),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for (cand, cand_pos), (ctrl, ctrl_pos), distance in pairs:
            writer.writerow(
                {
                    "candidate_site_key": cand,
                    "candidate_cys_position": cand_pos,
                    "control_site_key": ctrl,
                    "control_cys_position": ctrl_pos,
                    "sequence_distance_sd": f"{distance:.6g}",
                }
            )
    print(f"matched controls: {len(pairs)} pairs for {len(candidate_head)} candidates")

    generation_manifest = {
        "release_id": RELEASE_ID,
        "species": args.species,
        "artifact": {
            "candidates": "top_k_candidates.tsv",
            "candidates_sha256": _sha256_bytes(candidates_path),
            "matched_controls": "matched_controls.tsv",
            "matched_controls_sha256": _sha256_bytes(controls_path),
        },
        "reference_proteome": {
            "species": args.species,
            "path": species_ref["path"],
            "sha256": species_ref["sha256"],
        },
        "training_panel": {
            "n_rows": len(panel_rows),
            f"{args.species}_sites_excluded": len(species_panel_keys),
            "sha256": panel_sha256,
        },
        "scan": {
            "n_candidate_sites": len(ordered_keys),
            "n_with_structure": n_with_structure,
            "score_range": [
                min(s for _, s, _ in ranked),
                max(s for _, s, _ in ranked),
            ],
            "candidate_head_for_controls": CANDIDATE_RANK_HEAD,
            "controls_per_candidate": CONTROLS_PER_CANDIDATE,
            "control_caliper_sd": CONTROL_CALIPER_SD,
        },
        "bundle": {
            "path": "model_weights/structure_ranker_bundle.pt",
            "sha256": actual_bundle_sha,
        },
        "model": {
            "seed": bundle.seed,
            "hyperparameters": bundle.hyperparameters,
            "ablation": {
                "use_esm": bundle.ablation.use_esm,
                "use_study_context": bundle.ablation.use_study_context,
            },
        },
        "created": datetime.now(timezone.utc).isoformat(),
    }
    generation_path = output_dir / "top_k_generation_manifest.json"
    generation_path.write_text(
        json.dumps(generation_manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(f"wrote {candidates_path}")
    print(f"wrote {controls_path}")
    print(f"wrote {generation_path}")


def _z_score_features(
    features: dict[tuple[str, int], tuple[float, ...]],
    pool: list[tuple[str, int]],
) -> dict[tuple[str, int], tuple[float, ...]]:
    """z-score each of the 3 sequence features over ``pool`` (control pool)."""
    dims = len(next(iter(features.values())))
    means = [0.0] * dims
    for key in pool:
        for d, value in enumerate(features[key]):
            means[d] += value
    n = max(1, len(pool))
    means = [value / n for value in means]
    variances = [0.0] * dims
    for key in pool:
        for d, value in enumerate(features[key]):
            variances[d] += (value - means[d]) ** 2
    stds = [(variance / n) ** 0.5 for variance in variances]
    stds = [std if std > 0 else 1.0 for std in stds]
    z: dict[tuple[str, int], tuple[float, ...]] = {}
    for key, values in features.items():
        z[key] = tuple((value - mean) / std for value, mean, std in zip(values, means, stds, strict=True))
    return z


if __name__ == "__main__":
    main()
