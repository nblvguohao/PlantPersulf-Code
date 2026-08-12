"""Frozen, globally homology-blocked multispecies split construction.

This module is deliberately independent of model code.  It creates the one
immutable test boundary used by the v2 multispecies experiments and fails
closed when a protein has no registered global MMseqs2 cluster.
"""

from __future__ import annotations

import hashlib
import random
from collections import Counter, defaultdict
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

TEST_SPLIT = "test"
DEVELOPMENT_SPLIT = "development"
SPLIT_STATUS_OK = "OK"
SPLIT_STATUS_REVIEW = "DEVIATION_REQUIRES_REVIEW"
GLOBAL_CLUSTER_COLUMNS = (
    "species",
    "protein_accession",
    "global_protein_id",
    "cluster_id",
    "identity_threshold",
    "source_proteome_sha256",
)
FROZEN_SPLIT_COLUMNS = (
    "global_protein_id",
    "cluster_id",
    "split",
    "development_fold",
    "species",
    "source_sha256",
)


def global_protein_id(species: str, protein_accession: str) -> str:
    """Return the collision-safe key for one species-specific protein."""
    if not species or not protein_accession:
        raise ValueError("species and protein_accession are required")
    return f"{species}|{protein_accession}"


@dataclass(frozen=True)
class GlobalClusterRow:
    species: str
    protein_accession: str
    global_protein_id: str
    cluster_id: str
    identity_threshold: float
    source_proteome_sha256: str

    def __post_init__(self) -> None:
        expected = global_protein_id(self.species, self.protein_accession)
        if self.global_protein_id != expected:
            raise ValueError(
                "global_protein_id must be '{species}|{protein_accession}'"
            )
        if not self.cluster_id:
            raise ValueError("cluster_id is required")
        if not 0.0 < self.identity_threshold <= 1.0:
            raise ValueError("identity_threshold must be in (0, 1]")
        if len(self.source_proteome_sha256) != 64:
            raise ValueError("source_proteome_sha256 must be a SHA256")
        try:
            int(self.source_proteome_sha256, 16)
        except ValueError as exc:
            raise ValueError("source_proteome_sha256 must be hexadecimal") from exc


@dataclass(frozen=True)
class MultispeciesSiteRow:
    species: str
    protein_accession: str
    cys_position: int
    label: str
    study_accession: str

    @property
    def global_protein_id(self) -> str:
        return global_protein_id(self.species, self.protein_accession)


@dataclass(frozen=True)
class FrozenSplitRow:
    global_protein_id: str
    cluster_id: str
    split: str
    development_fold: int | None
    species: str
    source_sha256: str


@dataclass(frozen=True)
class FrozenMultispeciesSplit:
    rows: tuple[FrozenSplitRow, ...]
    seed: int
    test_fraction: float
    n_development_folds: int
    status: str
    test_cluster_fraction: float
    positive_test_fraction_by_stratum: dict[str, float]

    @property
    def sha256(self) -> str:
        payload = "\n".join(
            "\t".join(
                (
                    row.global_protein_id,
                    row.cluster_id,
                    row.split,
                    "" if row.development_fold is None else str(row.development_fold),
                    row.species,
                    row.source_sha256,
                )
            )
            for row in self.rows
        )
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()

    def require_trainable(self) -> None:
        if self.status != SPLIT_STATUS_OK:
            raise RuntimeError(
                "frozen split requires reviewer action: "
                f"status={self.status}"
            )


@dataclass(frozen=True)
class FrozenSplitAuditSummary:
    """Machine-readable result of a structural frozen-split leakage audit."""

    split_version: str
    split_sha256: str
    protein_count: int
    cluster_count: int
    test_protein_count: int
    development_protein_count: int
    development_fold_count: int
    status: str


def write_global_cluster_table(
    path: Path, rows: Iterable[GlobalClusterRow]
) -> None:
    """Write the only accepted on-disk schema for global MMseqs2 clusters."""
    import csv

    materialized = tuple(rows)
    if not materialized:
        raise RuntimeError("refusing to write an empty global cluster table")
    seen: set[str] = set()
    for row in materialized:
        if row.global_protein_id in seen:
            raise RuntimeError(f"duplicate global protein id: {row.global_protein_id}")
        seen.add(row.global_protein_id)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=GLOBAL_CLUSTER_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in sorted(materialized, key=lambda item: item.global_protein_id):
            writer.writerow(
                {
                    "species": row.species,
                    "protein_accession": row.protein_accession,
                    "global_protein_id": row.global_protein_id,
                    "cluster_id": row.cluster_id,
                    "identity_threshold": f"{row.identity_threshold:g}",
                    "source_proteome_sha256": row.source_proteome_sha256,
                }
            )


def load_global_cluster_table(path: Path) -> tuple[GlobalClusterRow, ...]:
    """Load a registered global cluster table and reject schema drift."""
    import csv

    if not path.is_file():
        raise FileNotFoundError(f"global cluster table not found: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != GLOBAL_CLUSTER_COLUMNS:
            raise RuntimeError(f"global cluster table has invalid columns: {path}")
        rows = tuple(
            GlobalClusterRow(
                species=row["species"],
                protein_accession=row["protein_accession"],
                global_protein_id=row["global_protein_id"],
                cluster_id=row["cluster_id"],
                identity_threshold=float(row["identity_threshold"]),
                source_proteome_sha256=row["source_proteome_sha256"],
            )
            for row in reader
        )
    if not rows:
        raise RuntimeError("global cluster table is empty")
    if len({row.global_protein_id for row in rows}) != len(rows):
        raise RuntimeError("global cluster table contains duplicate global protein ids")
    return rows


def write_frozen_multispecies_split(
    path: Path, split: FrozenMultispeciesSplit
) -> None:
    """Persist a frozen split once; overwriting a boundary is prohibited."""
    import csv
    import json

    if path.exists():
        raise FileExistsError(f"refusing to overwrite frozen split: {path}")
    audit_frozen_multispecies_split(split)
    metadata = {
        "seed": split.seed,
        "test_fraction": split.test_fraction,
        "n_development_folds": split.n_development_folds,
        "status": split.status,
        "test_cluster_fraction": split.test_cluster_fraction,
        "positive_test_fraction_by_stratum": split.positive_test_fraction_by_stratum,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("# metadata=" + json.dumps(metadata, sort_keys=True) + "\n")
        writer = csv.DictWriter(
            handle,
            fieldnames=FROZEN_SPLIT_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in split.rows:
            writer.writerow(
                {
                    "global_protein_id": row.global_protein_id,
                    "cluster_id": row.cluster_id,
                    "split": row.split,
                    "development_fold": ""
                    if row.development_fold is None
                    else row.development_fold,
                    "species": row.species,
                    "source_sha256": row.source_sha256,
                }
            )


def load_frozen_multispecies_split(path: Path) -> FrozenMultispeciesSplit:
    """Load the immutable split schema and verify its structural integrity."""
    import csv
    import json

    if not path.is_file():
        raise FileNotFoundError(f"frozen split not found: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        first_line = handle.readline()
        if not first_line.startswith("# metadata="):
            raise RuntimeError("frozen split metadata is missing")
        try:
            metadata = json.loads(first_line.removeprefix("# metadata="))
        except json.JSONDecodeError as exc:
            raise RuntimeError("frozen split metadata is invalid") from exc
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != FROZEN_SPLIT_COLUMNS:
            raise RuntimeError("frozen split has invalid columns")
        rows = tuple(
            FrozenSplitRow(
                global_protein_id=row["global_protein_id"],
                cluster_id=row["cluster_id"],
                split=row["split"],
                development_fold=(
                    None
                    if not row["development_fold"]
                    else int(row["development_fold"])
                ),
                species=row["species"],
                source_sha256=row["source_sha256"],
            )
            for row in reader
        )
    if not isinstance(metadata, dict):
        raise RuntimeError("frozen split metadata is invalid")
    result = FrozenMultispeciesSplit(
        rows=rows,
        seed=int(metadata["seed"]),
        test_fraction=float(metadata["test_fraction"]),
        n_development_folds=int(metadata["n_development_folds"]),
        status=str(metadata["status"]),
        test_cluster_fraction=float(metadata["test_cluster_fraction"]),
        positive_test_fraction_by_stratum={
            str(key): float(value)
            for key, value in dict(
                metadata["positive_test_fraction_by_stratum"]
            ).items()
        },
    )
    audit_frozen_multispecies_split(result)
    return result


def _stratum_counts(
    sites: Iterable[MultispeciesSiteRow],
) -> Counter[str]:
    """Count positives by species and species-study for split balancing."""
    counts: Counter[str] = Counter()
    for site in sites:
        if site.label != "positive":
            continue
        counts[f"species:{site.species}"] += 1
        if site.study_accession:
            counts[f"study:{site.species}|{site.study_accession}"] += 1
    return counts


def _cluster_site_counts(
    sites: Iterable[MultispeciesSiteRow],
    cluster_by_protein: dict[str, str],
) -> dict[str, Counter[str]]:
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    for site in sites:
        cid = cluster_by_protein.get(site.global_protein_id)
        if cid is None:
            raise RuntimeError(
                f"site protein missing from registered global clusters: "
                f"{site.global_protein_id}"
            )
        if site.label != "positive":
            continue
        counts[cid][f"species:{site.species}"] += 1
        if site.study_accession:
            counts[cid][f"study:{site.species}|{site.study_accession}"] += 1
    return counts


def _cluster_protein_counts(
    rows: Iterable[GlobalClusterRow],
) -> dict[str, Counter[str]]:
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    for row in rows:
        counts[row.cluster_id][f"protein:{row.species}"] += 1
    return counts


def _greedy_select_clusters(
    cluster_ids: list[str],
    cluster_counts: dict[str, Counter[str]],
    total_counts: Counter[str],
    n_select: int,
    fraction: float,
    seed: int,
) -> set[str]:
    """Pick exactly ``n_select`` clusters minimizing balance deviation.

    The deterministic greedy objective uses positive-site strata and protein
    strata.  It is a split-construction procedure, never a model-dependent
    optimization.
    """
    if not cluster_ids:
        raise RuntimeError("cannot split an empty cluster table")
    rng = random.Random(seed)
    tie_order = {cid: rng.random() for cid in cluster_ids}
    chosen: set[str] = set()
    observed: Counter[str] = Counter()
    ordered = sorted(cluster_ids, key=lambda cid: (tie_order[cid], cid))
    window_size = min(64, len(ordered))
    candidates = ordered[:window_size]
    next_candidate = window_size
    for _ in range(n_select):
        if not candidates:
            break

        def objective(cid: str) -> tuple[float, float, str]:
            projected = observed + cluster_counts.get(cid, Counter())
            error = 0.0
            for stratum, total in total_counts.items():
                target = total * fraction
                scale = max(target, 1.0)
                error += ((projected[stratum] - target) / scale) ** 2
            return (error, tie_order[cid], cid)

        pick = min(candidates, key=objective)
        chosen.add(pick)
        observed.update(cluster_counts.get(pick, Counter()))
        candidates.remove(pick)
        if next_candidate < len(ordered):
            candidates.append(ordered[next_candidate])
            next_candidate += 1
    return chosen


def build_frozen_multispecies_split(
    cluster_rows: Iterable[GlobalClusterRow],
    site_rows: Iterable[MultispeciesSiteRow],
    *,
    seed: int,
    test_fraction: float = 0.2,
    n_development_folds: int = 5,
    min_development_positive_sites: int = 1,
    min_development_positive_proteins: int = 1,
) -> FrozenMultispeciesSplit:
    """Create a deterministic global cluster holdout and development folds.

    All proteins, including proteins without observed positive sites, receive
    an assignment.  Positive sites are used only to balance the frozen split;
    no model metric is computed here.
    """
    if not 0.0 < test_fraction < 1.0:
        raise ValueError("test_fraction must be in (0, 1)")
    if n_development_folds < 2:
        raise ValueError("n_development_folds must be >= 2")
    if min_development_positive_sites < 1 or min_development_positive_proteins < 1:
        raise ValueError("development positive minimums must be >= 1")

    clusters = tuple(cluster_rows)
    sites = tuple(site_rows)
    if not clusters:
        raise RuntimeError("global cluster table is empty")
    by_global_id: dict[str, GlobalClusterRow] = {}
    for row in clusters:
        if row.global_protein_id in by_global_id:
            raise RuntimeError(f"duplicate global protein id: {row.global_protein_id}")
        by_global_id[row.global_protein_id] = row

    cluster_by_protein = {
        row.global_protein_id: row.cluster_id for row in clusters
    }
    positive_counts = _stratum_counts(sites)
    cluster_positive = _cluster_site_counts(sites, cluster_by_protein)
    cluster_proteins = _cluster_protein_counts(clusters)
    all_counts: Counter[str] = Counter(positive_counts)
    for per_cluster in cluster_proteins.values():
        all_counts.update(per_cluster)

    cluster_ids = sorted({row.cluster_id for row in clusters})
    n_test = max(1, round(len(cluster_ids) * test_fraction))
    test_clusters = _greedy_select_clusters(
        cluster_ids,
        {
            cid: (
                cluster_positive.get(cid, Counter())
                + cluster_proteins.get(cid, Counter())
            )
            for cid in cluster_ids
        },
        all_counts,
        n_test,
        test_fraction,
        seed,
    )
    development_clusters = sorted(set(cluster_ids) - test_clusters)
    if len(development_clusters) < n_development_folds:
        raise RuntimeError("not enough development clusters for requested folds")

    # Spread remaining clusters greedily across five development folds.
    fold_clusters: list[set[str]] = [set() for _ in range(n_development_folds)]
    fold_counts: list[Counter[str]] = [Counter() for _ in range(n_development_folds)]
    fold_positive_counts = [0 for _ in range(n_development_folds)]
    rng = random.Random(seed + 1)
    tie_order = {cid: rng.random() for cid in development_clusters}
    development_clusters.sort(
        key=lambda cid: (
            -sum(cluster_positive.get(cid, Counter()).values()),
            tie_order[cid],
            cid,
        )
    )
    positive_clusters = [
        cid
        for cid in development_clusters
        if sum(cluster_positive.get(cid, Counter()).values()) > 0
    ]
    if positive_clusters and len(positive_clusters) < n_development_folds:
        raise RuntimeError(
            "development split cannot provide one positive cluster per fold"
        )
    seeded_clusters = set(positive_clusters[:n_development_folds])
    for index, cid in enumerate(positive_clusters[:n_development_folds]):
        profile = cluster_positive.get(cid, Counter()) + cluster_proteins.get(
            cid, Counter()
        )
        fold_clusters[index].add(cid)
        fold_counts[index].update(profile)
        fold_positive_counts[index] += sum(
            cluster_positive.get(cid, Counter()).values()
        )
    for cid in development_clusters:
        if cid in seeded_clusters:
            continue
        profile = (
            cluster_positive.get(cid, Counter())
            + cluster_proteins.get(cid, Counter())
        )

        def fold_objective(
            index: int, candidate_profile: Counter[str] = profile
        ) -> tuple[float, int]:
            projected = fold_counts[index] + candidate_profile
            error = 0.0
            for stratum, total in all_counts.items():
                target = total * (1.0 - test_fraction) / n_development_folds
                error += ((projected[stratum] - target) / max(target, 1.0)) ** 2
            return (error, index)

        positive_total = sum(cluster_positive.get(cid, Counter()).values())
        index = min(
            range(n_development_folds),
            key=lambda candidate: (
                fold_positive_counts[candidate] if positive_total else 0,
                fold_objective(candidate),
            ),
        )
        fold_clusters[index].add(cid)
        fold_counts[index].update(profile)
        fold_positive_counts[index] += positive_total
    fold_by_cluster = {
        cid: fold for fold, members in enumerate(fold_clusters) for cid in members
    }

    rows: list[FrozenSplitRow] = []
    for row in sorted(clusters, key=lambda item: item.global_protein_id):
        is_test = row.cluster_id in test_clusters
        rows.append(
            FrozenSplitRow(
                global_protein_id=row.global_protein_id,
                cluster_id=row.cluster_id,
                split=TEST_SPLIT if is_test else DEVELOPMENT_SPLIT,
                development_fold=None if is_test else fold_by_cluster[row.cluster_id],
                species=row.species,
                source_sha256=row.source_proteome_sha256,
            )
        )

    test_positive_counts: Counter[str] = Counter()
    for site in sites:
        if (
            cluster_by_protein[site.global_protein_id] in test_clusters
            and site.label == "positive"
        ):
            test_positive_counts[f"species:{site.species}"] += 1
            if site.study_accession:
                key = f"study:{site.species}|{site.study_accession}"
                test_positive_counts[key] += 1
    fractions = {
        key: test_positive_counts[key] / count
        for key, count in positive_counts.items()
        if count
    }
    cluster_fraction = len(test_clusters) / len(cluster_ids)
    status = SPLIT_STATUS_OK
    if not 0.18 <= cluster_fraction <= 0.22:
        status = SPLIT_STATUS_REVIEW
    if any(not 0.15 <= fraction <= 0.25 for fraction in fractions.values()):
        status = SPLIT_STATUS_REVIEW
    positive_proteins_by_fold: list[set[str]] = [
        set() for _ in range(n_development_folds)
    ]
    positive_sites_by_fold = [0 for _ in range(n_development_folds)]
    for site in sites:
        cluster_id = cluster_by_protein[site.global_protein_id]
        fold = fold_by_cluster.get(cluster_id)
        if site.label == "positive" and fold is not None:
            positive_sites_by_fold[fold] += 1
            positive_proteins_by_fold[fold].add(site.global_protein_id)
    if any(
        positive_sites_by_fold[fold] < min_development_positive_sites
        or len(positive_proteins_by_fold[fold]) < min_development_positive_proteins
        for fold in range(n_development_folds)
    ):
        status = SPLIT_STATUS_REVIEW
    result = FrozenMultispeciesSplit(
        rows=tuple(rows),
        seed=seed,
        test_fraction=test_fraction,
        n_development_folds=n_development_folds,
        status=status,
        test_cluster_fraction=cluster_fraction,
        positive_test_fraction_by_stratum=fractions,
    )
    audit_frozen_multispecies_split(result)
    return result


def audit_frozen_multispecies_split(split: FrozenMultispeciesSplit) -> None:
    """Raise when any global protein or cluster crosses a frozen boundary."""
    if not split.rows:
        raise RuntimeError("frozen split is empty")
    if split.n_development_folds < 1:
        raise RuntimeError("frozen split requires development folds")
    by_protein: set[str] = set()
    by_cluster: dict[str, tuple[str, int | None]] = {}
    for row in split.rows:
        if row.split not in {DEVELOPMENT_SPLIT, TEST_SPLIT}:
            raise RuntimeError(
                f"invalid frozen split assignment: {row.global_protein_id}={row.split}"
            )
        expected_prefix = f"{row.species}|"
        if not row.species or not row.global_protein_id.startswith(expected_prefix):
            raise RuntimeError(
                f"global protein/species mismatch: {row.global_protein_id}"
            )
        if not row.cluster_id:
            raise RuntimeError(f"missing cluster id: {row.global_protein_id}")
        if len(row.source_sha256) != 64 or any(
            char not in "0123456789abcdef" for char in row.source_sha256.lower()
        ):
            raise RuntimeError(f"invalid source SHA256: {row.global_protein_id}")
        if row.global_protein_id in by_protein:
            raise RuntimeError(f"global protein appears twice: {row.global_protein_id}")
        by_protein.add(row.global_protein_id)
        assignment = (row.split, row.development_fold)
        previous = by_cluster.setdefault(row.cluster_id, assignment)
        if previous != assignment:
            raise RuntimeError(f"cluster crosses split boundary: {row.cluster_id}")
        if row.split == TEST_SPLIT and row.development_fold is not None:
            raise RuntimeError("test cluster must not have a development fold")
        if row.split == DEVELOPMENT_SPLIT and row.development_fold is None:
            raise RuntimeError("development cluster requires a development fold")
        if (
            row.development_fold is not None
            and not 0 <= row.development_fold < split.n_development_folds
        ):
            raise RuntimeError(
                f"development fold outside configured range: {row.global_protein_id}"
            )


def audit_frozen_split_file(
    path: Path, *, split_version: str
) -> FrozenSplitAuditSummary:
    """Load and structurally audit one explicitly versioned immutable split."""
    import re

    if re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_]*", split_version) is None:
        raise ValueError(f"invalid split version: {split_version}")
    if path.stem != split_version:
        raise RuntimeError(
            "split path/version mismatch: "
            f"expected {split_version}.tsv, got {path.name}"
        )
    split = load_frozen_multispecies_split(path)
    audit_frozen_multispecies_split(split)
    physical_sha256 = hashlib.sha256(path.read_bytes()).hexdigest()
    test_count = sum(row.split == TEST_SPLIT for row in split.rows)
    return FrozenSplitAuditSummary(
        split_version=split_version,
        split_sha256=physical_sha256,
        protein_count=len(split.rows),
        cluster_count=len({row.cluster_id for row in split.rows}),
        test_protein_count=test_count,
        development_protein_count=len(split.rows) - test_count,
        development_fold_count=split.n_development_folds,
        status=split.status,
    )


def partition_rows_by_frozen_split(
    rows: Iterable[MultispeciesSiteRow],
    split: FrozenMultispeciesSplit,
) -> dict[str, list[MultispeciesSiteRow]]:
    """Partition sites without any unclustered-protein fallback."""
    assignment = {row.global_protein_id: row.split for row in split.rows}
    result: dict[str, list[MultispeciesSiteRow]] = {
        DEVELOPMENT_SPLIT: [],
        TEST_SPLIT: [],
    }
    for row in rows:
        designated = assignment.get(row.global_protein_id)
        if designated is None:
            raise RuntimeError(
                f"site protein missing from frozen split: {row.global_protein_id}"
            )
        result[designated].append(row)
    return result
