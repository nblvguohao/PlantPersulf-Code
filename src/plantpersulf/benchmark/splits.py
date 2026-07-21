"""Task 6 — anti-leakage benchmark splits.

Produces a deterministic cluster-based train/validation/test split where
every protein in the same cluster is assigned to the same fold (no leakage
via homologous proteins), and the split is reproducible given the seed and
the frozen cluster file.

The known-mechanism holdout list (Zhang lab published sites) is recorded in
the split table but enforced at benchmark/validation time, not here; the
split itself does not reference the benchmark labels.
"""

from __future__ import annotations

import random
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

import yaml

SPLIT_LABELS = ("train", "validation", "test")


@dataclass(frozen=True)
class SplitRow:
    protein_accession: str
    cluster_id: str
    split: str


@dataclass(frozen=True)
class BenchmarkSplitTable:
    rows: tuple[SplitRow, ...]
    seed: int
    known_mechanism_holdout: tuple[str, ...]

    @property
    def split_counts(self) -> dict[str, int]:
        counts: dict[str, int] = dict.fromkeys(SPLIT_LABELS, 0)
        for row in self.rows:
            counts[row.split] += 1
        return counts


def _load_clusters(path: Path) -> list[tuple[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"cluster file not found: {path}")
    import csv

    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = tuple(reader.fieldnames or ())
        if fieldnames != ("protein_accession", "cluster_id"):
            raise RuntimeError(f"cluster file has invalid columns: {path}")
        rows = [(row["protein_accession"], row["cluster_id"]) for row in reader]
    if not rows:
        raise RuntimeError("cluster file is empty")
    return rows


def _load_config(config_path: Path) -> dict[str, Any]:
    try:
        loaded: object = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise RuntimeError(f"invalid split config: {config_path}") from exc
    if not isinstance(loaded, dict) or loaded.get("version") != 1:
        raise RuntimeError("split config requires version 1")
    return cast(dict[str, Any], loaded)


def build_splits(config_path: Path) -> BenchmarkSplitTable:
    """Build the deterministic cluster-level train/validation/test splits."""
    config = _load_config(config_path)
    cs = config["cluster_split"]
    if not isinstance(cs, dict):
        raise RuntimeError("split config requires cluster_split mapping")
    cluster_file = Path(str(cs["cluster_file"]))
    train_frac = float(cs["train_fraction"])
    val_frac = float(cs["validation_fraction"])
    test_frac = float(cs["test_fraction"])
    if abs(train_frac + val_frac + test_frac - 1.0) > 1e-9:
        raise RuntimeError("cluster split fractions must sum to 1.0")
    seed = int(cs["seed"])
    holdout = config.get("known_mechanism_holdout")
    holdout = (
        tuple(str(g) for g in holdout)
        if isinstance(holdout, list)
        else ()
    )

    protein_clusters = _load_clusters(cluster_file)
    clusters: dict[str, list[str]] = {}
    for protein, cid in protein_clusters:
        clusters.setdefault(cid, []).append(protein)

    cids = sorted(clusters)
    rng = random.Random(seed)
    rng.shuffle(cids)
    n = len(cids)
    n_train = max(1, round(n * train_frac))
    n_val = max(1, round(n * val_frac))
    train_cids = set(cids[:n_train])
    val_cids = set(cids[n_train : n_train + n_val])
    _ = set(cids[n_train + n_val :])  # test (unused inline but verified by complement)

    rows: list[SplitRow] = []
    for cid in sorted(clusters):
        split = (
            "train"
            if cid in train_cids
            else "validation"
            if cid in val_cids
            else "test"
        )
        for protein in sorted(clusters[cid]):
            rows.append(
                SplitRow(
                    protein_accession=protein, cluster_id=cid, split=split
                )
            )
    return BenchmarkSplitTable(
        rows=tuple(rows), seed=seed, known_mechanism_holdout=holdout
    )
