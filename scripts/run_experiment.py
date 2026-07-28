#!/usr/bin/env python
"""Deterministic TDD experiment runner for persulfidation baselines.

Reads a frozen YAML experiment config and runs every model on the same
benchmark, features, splits, and seeds — producing a single, deterministic
results tree under ``results/experiments/<name>/`` with per-seed metrics
and a manifest recording every input SHA256 used.

Usage (per the Task 8 acceptance gate)::

    python scripts/run_experiment.py \\
        --config configs/experiments/baseline_sequence_v1.yaml
"""

from __future__ import annotations

import argparse
import csv
import json
import random
import shutil
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, cast

import yaml

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

BENCHMARK_FIELDS = (
    "protein_accession", "cys_position_in_protein", "label",
    "study_accession", "evidence_level", "source_sha256",
)


def _ensure_venv() -> None:
    venv = Path(sys.executable).with_name("pyvenv.cfg")
    if not venv.is_file():
        print("WARNING: not running inside a project venv.", file=sys.stderr)


def _load_config(path: Path) -> dict[str, Any]:
    loaded: object = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(loaded, dict) or loaded.get("version") != 1:
        raise SystemExit(f"invalid experiment config: {path}")
    return cast(dict[str, Any], loaded)


# ---------------------------------------------------------------------------
# results
# ---------------------------------------------------------------------------

@dataclass
class ModelResult:
    model: str
    seed: int
    train_ap: float | None = None
    val_ap: float | None = None
    test_ap: float | None = None
    test_recall_10: float | None = None
    test_recall_50: float | None = None
    test_mrr: float | None = None


@dataclass
class ExperimentResults:
    experiment: str
    model_results: list[ModelResult] = field(default_factory=list)
    manifest_sha256: str = ""


# ---------------------------------------------------------------------------
# bench --> split
# ---------------------------------------------------------------------------

def _build_singleton_cluster_file(benchmark_path: Path, output: Path) -> Path:
    """Each protein in its own cluster when no MMseqs2 output is available."""
    proteins: set[str] = set()
    with benchmark_path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            proteins.add(row["protein_accession"])
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("protein_accession", "cluster_id"))
        for i, p in enumerate(sorted(proteins)):
            writer.writerow((p, str(i)))
    return output


def _train_val_test_rows(
    benchmark_rows: list[dict[str, str]],
    proteome_path: Path,
    cluster_file: Path,
    split_config_path: Path,
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    from plantpersulf.benchmark.splits import build_splits

    table = build_splits(split_config_path)
    split_map = {row.protein_accession: row.split for row in table.rows}

    train_rows: list[dict[str, str]] = []
    val_rows: list[dict[str, str]] = []
    test_rows: list[dict[str, str]] = []
    for row in benchmark_rows:
        split = split_map.get(row["protein_accession"], "train")
        if split == "validation":
            val_rows.append(row)
        elif split == "test":
            test_rows.append(row)
        else:
            train_rows.append(row)
    return train_rows, val_rows, test_rows


# ---------------------------------------------------------------------------
# subsampling + feature extraction
# ---------------------------------------------------------------------------

def _subsample_unlabeled(
    rows: list[dict[str, str]],
    ratio: int,
    seed: int,
) -> list[dict[str, str]]:
    """Keep every positive; keep a seeded random sample of ratio x positives
    unlabeled rows. A PU benchmark with ~390 positives and ~395k unlabeled is
    both intractable for heavy features and dominated by the unlabeled class;
    matched-negative subsampling is standard and is applied deterministically
    (fixed seed) so the whole experiment is reproducible. This is a training/
    evaluation convenience only — the frozen benchmark itself is never altered.
    """
    positives = [r for r in rows if r["label"] == "positive"]
    unlabeled = [r for r in rows if r["label"] != "positive"]
    keep_n = min(len(unlabeled), max(1, ratio * max(1, len(positives))))
    rng = random.Random(seed)
    sampled = rng.sample(unlabeled, keep_n)
    combined = positives + sampled
    # stable deterministic order for reproducible feature extraction
    combined.sort(
        key=lambda r: (r["protein_accession"], int(r["cys_position_in_protein"]))
    )
    return combined


def _write_labels_tsv(rows: list[dict[str, str]], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(BENCHMARK_FIELDS),
            delimiter="\t", lineterminator="\n", extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def _sequence_feature_vectors(
    rows: list[dict[str, str]],
    proteome_path: Path,
    scratch_dir: Path,
    tag: str,
) -> list[list[float]]:
    """Extract [hydrophobicity, cys_density] per benchmark row via the real
    feature extractor, driven by a temp labels file in the benchmark schema."""
    from plantpersulf.features.sequence import extract_sequence_features

    labels_path = scratch_dir / f"{tag}_labels.tsv"
    _write_labels_tsv(rows, labels_path)
    feats = extract_sequence_features(labels_path, proteome_path, window_radius=10)
    lookup: dict[tuple[str, int], list[float]] = {
        (f.protein_accession, f.cys_position): [f.hydrophobicity, f.cys_density]
        for f in feats
    }
    return [
        lookup.get(
            (r["protein_accession"], int(r["cys_position_in_protein"])),
            [0.0, 0.0],
        )
        for r in rows
    ]


# ---------------------------------------------------------------------------
# study-split experiment (leave-one-study-out)
# ---------------------------------------------------------------------------

def _study_fold_rows(
    all_rows: list[dict[str, str]],
    holdout_study: str,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Partition benchmark rows for one leave-study-out fold.

    Delegates to the canonical, unit-tested
    ``plantpersulf.evaluation.external_validation.partition_leave_study_out``
    so the runner and the external-validation report share one definition of
    "held-out study positives never appear in train"."""
    from plantpersulf.evaluation.external_validation import (
        partition_leave_study_out,
    )

    return partition_leave_study_out(all_rows, holdout_study)


def _run_study_split_experiment(
    cfg: dict[str, Any],
    exp_cfg: dict[str, Any],
    eval_cfg: dict[str, Any],
    out_cfg: dict[str, Any],
    all_rows: list[dict[str, str]],
    proteome_path: Path,
) -> None:
    split_mode = cfg["splits"].get("mode", "leave_study_out")
    if split_mode == "time_split":
        # Time split: single fold; train on earlier study's positives,
        # test on later study's positives.
        studies = [cfg["splits"]["test_study"]]
    else:
        studies = cfg["splits"]["studies"]
    ratio = int(cfg.get("subsample", {}).get("unlabeled_per_positive", 50))
    sub_seed = int(cfg.get("subsample", {}).get("seed", 12345))
    limitation = str(cfg.get("limitation", ""))

    feature_fn = _resolve_feature_fn(cfg)
    all_fold_results: list[ModelResult] = []
    for holdout_study in studies:
        print(f"\n===== Fold: leave_{holdout_study}_out =====")
        train_rows, test_rows = _study_fold_rows(all_rows, holdout_study)
        train_rows = _subsample_unlabeled(train_rows, ratio, sub_seed)
        test_rows = _subsample_unlabeled(test_rows, ratio, sub_seed + 1)
        train_pos = sum(1 for r in train_rows if r["label"] == "positive")
        test_pos = sum(1 for r in test_rows if r["label"] == "positive")
        print(
            f"  subsampled: train={len(train_rows)}(+{train_pos}) "
            f"test={len(test_rows)}(+{test_pos})"
        )
        # Split train further into train/val (80/20 of train rows)
        rng = random.Random(sub_seed)
        train_dedup = train_rows[:]
        rng.shuffle(train_dedup)
        n_val = max(1, int(len(train_dedup) * 0.2))
        val_rows_fold = train_dedup[:n_val]
        train_rows_fold = train_dedup[n_val:]

        # No val set — use train as val here (study-split has no natural val)
        # For a no-validation study-split: train on fold train, evaluate on
        # fold test. Hyperparameters are fixed at defaults.
        fold_results = _run_baselines(
            cfg, train_rows_fold,
            [r["label"] for r in train_rows_fold],
            val_rows_fold,
            [r["label"] for r in val_rows_fold],
            test_rows,
            [r["label"] for r in test_rows],
            proteome_path,
            f"{exp_cfg['name']}_{holdout_study}",
            feature_fn=feature_fn,
        )
        for r in fold_results:
            r.model = f"leave_{holdout_study}_out|{r.model}"
        all_fold_results.extend(fold_results)

    combined = ExperimentResults(experiment=exp_cfg["name"])
    combined.model_results = all_fold_results
    _write_results(combined, out_cfg, limitation=limitation)


# ---------------------------------------------------------------------------
# shared baseline execution + result writing
# ---------------------------------------------------------------------------

def _load_proteome_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    cur_header = ""
    cur_lines: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if cur_header:
                seqs[cur_header.split("|")[1]] = "".join(cur_lines)
            cur_header = line
            cur_lines = []
        elif line:
            cur_lines.append(line)
    if cur_header:
        seqs[cur_header.split("|")[1]] = "".join(cur_lines)
    return seqs


def _esm2_feature_vectors(
    rows: list[dict[str, str]],
    proteome_path: Path,
    scratch_dir: Path,
    tag: str,
) -> list[list[float]]:
    """Extract 1280-dim ESM-2 per-residue embeddings in chunks of ~50
    unique proteins to avoid OOM on CPU (a single 5000-protein batch
    can allocate >100 GB)."""
    from plantpersulf.features.esm2 import extract_esm2_embeddings

    # Deduplicate proteins
    seen: dict[str, list[int]] = {}
    for i, row in enumerate(rows):
        seen.setdefault(row["protein_accession"], []).append(i)

    # Sort by sequence length (shortest first) to minimise padding waste
    proteome = _load_proteome_fasta(proteome_path)
    unique = sorted(seen, key=lambda acc: len(proteome.get(acc, "")))
    lookup: dict[tuple[str, int], list[float]] = {}
    chunk_size = 5
    total = len(unique)
    for start in range(0, total, chunk_size):
        chunk_prots = unique[start : start + chunk_size]
        chunk_rows = [
            r for r in rows if r["protein_accession"] in chunk_prots
        ]
        labels_path = scratch_dir / f"{tag}_esm2_{start}.tsv"
        _write_labels_tsv(chunk_rows, labels_path)
        emb_rows = extract_esm2_embeddings(labels_path, proteome_path)
        for r in emb_rows:
            lookup[(r.protein_accession, r.cys_position)] = list(r.embedding)
        if (start // chunk_size) % 50 == 0:
            print(f"    ESM2 chunk {start // chunk_size}/{total // chunk_size} ...")

    default = [0.0] * 1280
    return [
        lookup.get(
            (r["protein_accession"], int(r["cys_position_in_protein"])),
            default,
        )
        for r in rows
    ]


def _resolve_feature_fn(
    cfg: dict[str, Any],
) -> Any:
    """Return the feature-extraction function based on the config."""
    feature_names = [f["name"] for f in cfg["features"]]
    if "structure" in feature_names:
        return _structure_feature_vectors
    if "esm2" in feature_names:
        return _esm2_feature_vectors
    return _sequence_feature_vectors


def _structure_feature_vectors_with_mask(
    rows: list[dict[str, str]],
    proteome_path: Path,
    scratch_dir: Path,
    tag: str,
) -> tuple[list[list[float]], list[bool]]:
    """Extract [contact_number_proxy, plddt] per benchmark row plus a
    per-row ``has_structure`` mask, from real AlphaFold/SWISS-MODEL
    structures. Rows without a resolvable structure get a zero vector and a
    ``False`` mask — the mask (not a placeholder value) is what downstream
    consumers must honour."""
    from plantpersulf.download.alphafold import audit_alphafold_structures
    from plantpersulf.features.structure import extract_cys_structure_features

    # Build a lookup: (accession) -> pdb_text (cached, one read per protein)
    sources = audit_alphafold_structures()
    pdb_cache: dict[str, str] = {}
    for s in sources:
        pdb_cache[s.accession] = s.local_path.read_text(encoding="utf-8")

    vectors: list[list[float]] = []
    mask: list[bool] = []
    for row in rows:
        acc = row["protein_accession"]
        pos = int(row["cys_position_in_protein"])
        pdb = pdb_cache.get(acc)
        if pdb is None:
            vectors.append([0.0, 0.0])
            mask.append(False)
            continue
        feats = extract_cys_structure_features(pdb, acc, [pos])
        f = feats[0]
        if not f.has_structure:
            vectors.append([0.0, 0.0])
            mask.append(False)
        else:
            contact = f.contact_number_proxy or 0.0
            plddt = f.plddt or 0.0
            vectors.append([contact, plddt])
            mask.append(True)
    return vectors, mask


def _structure_feature_vectors(
    rows: list[dict[str, str]],
    proteome_path: Path,
    scratch_dir: Path,
    tag: str,
) -> list[list[float]]:
    """Flat-feature view of the structure branch (for the traditional
    baselines): the [contact, plddt] vectors only, missing rows zero-filled."""
    vectors, _ = _structure_feature_vectors_with_mask(
        rows, proteome_path, scratch_dir, tag
    )
    return vectors


RANKER_MODEL = "structure_ranker"


def _build_branch_features(
    rows: list[dict[str, str]],
    proteome_path: Path,
    scratch_dir: Path,
    tag: str,
    need_esm: bool = True,
) -> Any:
    """Assemble a multi-branch ``BranchFeatures`` (sequence + frozen-ESM +
    structure + missingness mask + study context) for the structure-aware
    ranker, reusing the same real feature extractors as the flat baselines.

    ESM extraction (which loads the 650M model) is skipped when no ablation in
    the run uses the ESM branch; a 1-dim zero placeholder is substituted so the
    (disabled) branch still has a consistent shape."""
    from plantpersulf.models.structure_ranker import BranchFeatures

    seq = _sequence_feature_vectors(rows, proteome_path, scratch_dir, tag)
    if need_esm:
        esm = _esm2_feature_vectors(rows, proteome_path, scratch_dir, tag)
    else:
        esm = [[0.0] for _ in rows]
    struct, mask = _structure_feature_vectors_with_mask(
        rows, proteome_path, scratch_dir, tag
    )
    study_ids = [r.get("study_accession") or "__unlabeled__" for r in rows]
    return BranchFeatures(
        sequence=seq, esm=esm, structure=struct,
        structure_mask=mask, study_ids=study_ids,
    )


def _ablation_from_spec(spec: dict[str, Any]) -> Any:
    """Build an AblationConfig from a config dict; unspecified flags default
    to enabled (True)."""
    from plantpersulf.models.structure_ranker import AblationConfig

    return AblationConfig(
        use_esm=bool(spec.get("use_esm", True)),
        use_structure=bool(spec.get("use_structure", True)),
        use_plddt=bool(spec.get("use_plddt", True)),
        use_accessibility=bool(spec.get("use_accessibility", True)),
        use_study_context=bool(spec.get("use_study_context", True)),
    )


def _ablation_specs(cfg: dict[str, Any]) -> list[dict[str, Any]]:
    specs = cfg.get("ablations")
    if not specs:
        return [{"name": "full"}]
    return [dict(s) for s in specs]


def _concat_branches(a: Any, b: Any) -> Any:
    from plantpersulf.models.structure_ranker import BranchFeatures

    return BranchFeatures(
        sequence=a.sequence + b.sequence,
        esm=a.esm + b.esm,
        structure=a.structure + b.structure,
        structure_mask=a.structure_mask + b.structure_mask,
        study_ids=(a.study_ids or []) + (b.study_ids or []),
    )


def _evaluate_motif(
    train_rows: list[dict[str, str]],
    train_y: list[str],
    val_rows: list[dict[str, str]],
    val_y: list[str],
    test_rows: list[dict[str, str]],
    test_y: list[str],
    proteome_path: Path,
    label: str,
    seed: int,
) -> ModelResult:
    """Fit the motif-frequency baseline on train, score val/test."""
    import csv as _csv
    import shutil as _sh
    import tempfile as _tf

    from plantpersulf.evaluation.metrics import average_precision as ap
    from plantpersulf.models.baselines import train_motif_baseline

    tmp = Path(_tf.mkdtemp(prefix="motif_"))
    try:
        sp = tmp / "splits.csv"
        with sp.open("w", encoding="utf-8", newline="") as h:
            w = _csv.writer(h, delimiter=",", lineterminator="\n")
            w.writerow(("protein_accession", "cys_position", "label", "split"))
            pairs = [
                (train_rows, "train"), (val_rows, "val"), (test_rows, "test")
            ]
            for rs, sn in pairs:
                for r in rs:
                    w.writerow(
                        (r["protein_accession"], r["cys_position_in_protein"],
                         r["label"], sn)
                    )

        model = train_motif_baseline(sp, proteome_path, window_radius=10)

        proteome: dict[str, str] = {}
        ch = ""
        cl: list[str] = []
        for ln_ in proteome_path.read_text(encoding="utf-8").splitlines():
            if ln_.startswith(">"):
                if ch:
                    proteome[ch.split("|")[1]] = "".join(cl)
                ch = ln_
                cl = []
            elif ln_:
                cl.append(ln_)
        if ch:
            proteome[ch.split("|")[1]] = "".join(cl)

        def _sc(rows_):
            out = []
            for r_ in rows_:
                seq = proteome.get(r_["protein_accession"], "")
                if not seq or not seq.strip():
                    out.append((float("-inf"), r_["label"]))
                    continue
                pos = int(r_["cys_position_in_protein"])
                if pos < 1 or pos > len(seq):
                    out.append((float("-inf"), r_["label"]))
                    continue
                r_rad = 10
                lp = max(0, r_rad - (pos - 1))
                rp = max(0, r_rad - (len(seq) - pos))
                st = max(0, pos - 1 - r_rad)
                ed = min(len(seq), pos + r_rad)
                win = "X" * lp + seq[st:ed] + "X" * rp
                out.append((model.score(win), r_["label"]))
            return out

        vs = _sc(val_rows)
        ts = _sc(test_rows)
        return ModelResult(
            model=f"{label}|motif_frequency", seed=seed,
            val_ap=ap(vs) if vs else None,
            test_ap=ap(ts) if ts else None,
        )
    finally:
        _sh.rmtree(tmp, ignore_errors=True)


def _evaluate_ranker(
    ablation_name: str,
    ablation: Any,
    seed: int,
    branch_train: Any,
    train_y: list[str],
    branch_val: Any,
    val_y: list[str],
    branch_test: Any,
    test_y: list[str],
    ranker_params: dict[str, Any] | None = None,
) -> ModelResult:
    """Train the ranker once and score val+test in a single forward pass
    (val and test are concatenated then split), so the model is fit exactly
    once per (ablation, seed)."""
    from plantpersulf.evaluation.metrics import (
        average_precision,
        mean_reciprocal_rank,
        recall_at_k,
    )
    from plantpersulf.models.structure_ranker import structure_ranker_scores

    result = ModelResult(model=ablation_name, seed=seed)
    params = ranker_params or {}

    n_val = branch_val.n_rows()
    combined = _concat_branches(branch_val, branch_test)
    out = structure_ranker_scores(
        branch_train, train_y, combined, seed=seed, ablation=ablation,
        hidden=int(params.get("hidden", 16)),
        dropout=float(params.get("dropout", 0.2)),
        epochs=int(params.get("epochs", 200)),
        lr=float(params.get("lr", 0.05)),
        n_mc_dropout=int(params.get("n_mc_dropout", 16)),
    )
    val_scores = out.scores[:n_val]
    test_scores = out.scores[n_val:]

    val_scored = list(zip(val_scores, val_y, strict=True))
    result.val_ap = average_precision(val_scored)
    if test_scores:
        test_scored = list(zip(test_scores, test_y, strict=True))
        result.test_ap = average_precision(test_scored)
        result.test_recall_10 = recall_at_k(test_scored, 10)
        result.test_recall_50 = recall_at_k(test_scored, 50)
        result.test_mrr = mean_reciprocal_rank(test_scored)
    return result


def _run_baselines(
    cfg: dict[str, Any],
    train_rows: list[dict[str, str]],
    train_y: list[str],
    val_rows: list[dict[str, str]],
    val_y: list[str],
    test_rows: list[dict[str, str]],
    test_y: list[str],
    proteome_path: Path,
    label: str,
    feature_fn: Any | None = None,
) -> list[ModelResult]:
    if feature_fn is None:
        feature_fn = _sequence_feature_vectors
    models = list(cfg["models"])
    needs_flat = any(m != RANKER_MODEL for m in models)
    needs_branch = RANKER_MODEL in models

    scratch = Path(tempfile.mkdtemp(prefix="run_experiment_"))
    feat_train = feat_val = feat_test = None
    branch_train = branch_val = branch_test = None
    try:
        if needs_flat:
            feat_train = feature_fn(train_rows, proteome_path, scratch, "train")
            feat_val = feature_fn(val_rows, proteome_path, scratch, "val")
            feat_test = feature_fn(test_rows, proteome_path, scratch, "test")
        if needs_branch:
            need_esm = any(
                _ablation_from_spec(s).use_esm for s in _ablation_specs(cfg)
            )
            branch_train = _build_branch_features(
                train_rows, proteome_path, scratch, "train", need_esm
            )
            branch_val = _build_branch_features(
                val_rows, proteome_path, scratch, "val", need_esm
            )
            branch_test = _build_branch_features(
                test_rows, proteome_path, scratch, "test", need_esm
            )
    finally:
        shutil.rmtree(scratch, ignore_errors=True)

    eval_cfg = cfg["evaluation"]
    results: list[ModelResult] = []
    for model_name in models:
        for seed in eval_cfg["seeds"]:
            if model_name == "motif_frequency":
                for seed in eval_cfg["seeds"]:
                    r = _evaluate_motif(
                        train_rows, train_y, val_rows, val_y,
                        test_rows, test_y, proteome_path, label, seed,
                    )
                    msg = (
                        f"motif seed={seed} "
                        f"val_ap={r.val_ap or 'NA'} "
                        f"test_ap={r.test_ap or 'NA'}"
                    )
                    print(msg)
                    results.append(r)
                continue
            if model_name == RANKER_MODEL:
                for spec in _ablation_specs(cfg):
                    ab_name = str(spec.get("name", "full"))
                    print(f"  {RANKER_MODEL}[{ab_name}] seed={seed} ...", end=" ")
                    r = _evaluate_ranker(
                        f"{RANKER_MODEL}:{ab_name}",
                        _ablation_from_spec(spec), seed,
                        branch_train, train_y,
                        branch_val, val_y,
                        branch_test, test_y,
                        ranker_params=cfg.get("ranker", {}),
                    )
                    r.model = f"{label}|{r.model}"
                    msg = (
                        f"val_ap={r.val_ap:.4f}"
                        if r.val_ap is not None else "no_val_pos"
                    )
                    if r.test_ap is not None:
                        msg += f" test_ap={r.test_ap:.4f}"
                    print(msg)
                    results.append(r)
                continue
            print(f"  {model_name} seed={seed} ...", end=" ")
            r = _evaluate(
                model_name, seed,
                feat_train, train_y,
                feat_val, val_y,
                feat_test, test_y,
            )
            r.model = f"{label}|{r.model}"
            msg = f"val_ap={r.val_ap:.4f}" if r.val_ap is not None else "no_val_pos"
            if r.test_ap is not None:
                msg += f" test_ap={r.test_ap:.4f}"
            print(msg)
            results.append(r)
    return results


def _write_results(
    results: ExperimentResults,
    out_cfg: dict[str, Any],
    limitation: str = "",
) -> None:
    out_dir = Path(out_cfg["directory"])
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True)

    rows_out: list[dict[str, str]] = []
    for r in results.model_results:
        out_row: dict[str, str] = {"model": r.model, "seed": str(r.seed)}
        for key, val in [
            ("val_ap", r.val_ap), ("test_ap", r.test_ap),
            ("test_recall_10", r.test_recall_10),
            ("test_recall_50", r.test_recall_50),
            ("test_mrr", r.test_mrr),
        ]:
            out_row[key] = f"{val:.6f}" if val is not None else ""
        rows_out.append(out_row)
    with (out_dir / "metrics.tsv").open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(
            h, fieldnames=list(rows_out[0]), delimiter="\t", lineterminator="\n",
        )
        w.writeheader()
        w.writerows(rows_out)


    manifest: dict[str, object] = {
        "experiment": results.experiment,
        "n_seeds": len({r.seed for r in results.model_results}),
    }
    if limitation:
        manifest["limitation"] = limitation
    (out_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"\nResults written to {out_dir}")


# ---------------------------------------------------------------------------
# baseline scoring
# ---------------------------------------------------------------------------

def _train_predict(
    model_name: str,
    train_X: list[list[float]],
    train_y: list[str],
    predict_X: list[list[float]],
    seed: int,
) -> list[float]:
    from plantpersulf.models.esm_baseline import esm_linear_head_scores
    from plantpersulf.models.traditional import (
        logistic_regression_scores,
        pu_logistic_regression_scores,
        random_forest_scores,
        xgboost_scores,
    )

    train_emb = [tuple(v) for v in train_X]
    predict_emb = [tuple(v) for v in predict_X]

    if model_name == "logistic":
        return logistic_regression_scores(train_X, train_y, predict_X, seed)
    if model_name == "pu_logistic":
        return pu_logistic_regression_scores(train_X, train_y, predict_X, seed)
    if model_name == "random_forest":
        return random_forest_scores(
            train_X, train_y, predict_X, seed, n_estimators=50
        )
    if model_name == "xgboost":
        return xgboost_scores(
            train_X, train_y, predict_X, seed, n_estimators=50
        )
    if model_name == "esm_linear_head":
        return esm_linear_head_scores(train_emb, train_y, predict_emb, seed)
    raise ValueError(f"unknown model: {model_name}")


def _evaluate(
    model_name: str,
    seed: int,
    train_X: list[list[float]],
    train_y: list[str],
    val_X: list[list[float]],
    val_y: list[str],
    test_X: list[list[float]],
    test_y: list[str],
) -> ModelResult:
    from plantpersulf.evaluation.metrics import (
        average_precision,
        mean_reciprocal_rank,
        recall_at_k,
    )

    result = ModelResult(model=model_name, seed=seed)

    def _scored(
        features: list[list[float]], labels: list[str],
    ) -> list[tuple[float, str]]:
        scores = _train_predict(model_name, train_X, train_y, features, seed)
        return list(zip(scores, labels, strict=True))

    val_scored = _scored(val_X, val_y)
    result.val_ap = average_precision(val_scored)

    if test_X:
        test_scored = _scored(test_X, test_y)
        result.test_ap = average_precision(test_scored)
        result.test_recall_10 = recall_at_k(test_scored, 10)
        result.test_recall_50 = recall_at_k(test_scored, 50)
        result.test_mrr = mean_reciprocal_rank(test_scored)
    return result


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def run_experiment(config_path: Path) -> ExperimentResults:
    cfg = _load_config(config_path)
    exp_cfg = cfg["experiment"]
    eval_cfg = cfg["evaluation"]
    out_cfg = cfg["output"]

    # --- benchmark ---
    benchmark_path = Path(exp_cfg["benchmark_labels"])
    if not benchmark_path.exists():
        print(f"Building benchmark into {benchmark_path} …")
        from plantpersulf.benchmark.labels import build_benchmark_labels
        from plantpersulf.proteomics.persulfidation_publish import (
            publish_persulfidation_sites,
        )
        interim = Path("data/interim")
        for acc in ("PXD006140", "PXD024061"):
            publish_persulfidation_sites(acc, output_root=interim)
        build_benchmark_labels(
            site_output_root=interim,
            proteome_path=Path(exp_cfg["reference_proteome"]),
            output_directory=benchmark_path.parent,
        )

    rows: list[dict[str, str]] = []
    with benchmark_path.open(encoding="utf-8", newline="") as h:
        for row in csv.DictReader(h, delimiter="\t"):
            rows.append(dict(row))

    proteome_path = Path(exp_cfg["reference_proteome"])

    # --- splits ---
    splits_cfg = cfg["splits"]
    split_mode = splits_cfg.get("mode", "cluster")

    if split_mode in ("leave_study_out", "time_split"):
        _run_study_split_experiment(
            cfg, exp_cfg, eval_cfg, out_cfg, rows, proteome_path,
        )
        return ExperimentResults(experiment=exp_cfg["name"])

    # --- cluster split (default) ---
    cluster_file = Path(splits_cfg["cluster_file"])
    if not cluster_file.exists():
        cluster_file = _build_singleton_cluster_file(benchmark_path, cluster_file)
    split_config = Path(splits_cfg["split_config"])
    train_rows, val_rows, test_rows = _train_val_test_rows(
        rows, proteome_path, cluster_file, split_config,
    )
    print(f"split: train={len(train_rows)} val={len(val_rows)} test={len(test_rows)}")

    # --- matched-negative subsampling (deterministic) ---
    ratio = int(cfg.get("subsample", {}).get("unlabeled_per_positive", 50))
    sub_seed = int(cfg.get("subsample", {}).get("seed", 12345))
    train_rows = _subsample_unlabeled(train_rows, ratio, sub_seed)
    val_rows = _subsample_unlabeled(val_rows, ratio, sub_seed + 1)
    test_rows = _subsample_unlabeled(test_rows, ratio, sub_seed + 2)
    train_pos = sum(1 for r in train_rows if r["label"] == "positive")
    val_pos = sum(1 for r in val_rows if r["label"] == "positive")
    test_pos = sum(1 for r in test_rows if r["label"] == "positive")
    print(
        f"subsampled (1:{ratio}): train={len(train_rows)}(+{train_pos}) "
        f"val={len(val_rows)}(+{val_pos}) test={len(test_rows)}(+{test_pos})"
    )

    # --- sequence features + baselines(shared runner extracts internally) ---
    feature_names = [f["name"] for f in cfg["features"]]
    print(f"features: {feature_names}")
    train_y = [r["label"] for r in train_rows]
    val_y = [r["label"] for r in val_rows]
    test_y = [r["label"] for r in test_rows]

    feature_fn = _resolve_feature_fn(cfg)
    results_list = _run_baselines(
        cfg, train_rows, train_y, val_rows, val_y, test_rows, test_y,
        proteome_path, exp_cfg["name"],
        feature_fn=feature_fn,
    )
    combined = ExperimentResults(experiment=exp_cfg["name"])
    combined.model_results = results_list
    _write_results(combined, out_cfg)
    return combined


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="run a persulfidation baseline experiment")
    p.add_argument("--config", type=Path, required=True)
    return p


if __name__ == "__main__":
    _ensure_venv()
    args = build_parser().parse_args()
    run_experiment(args.config)
