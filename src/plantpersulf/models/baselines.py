r"""Task 8 — leakage-safe PU baselines (no-learning, motif-frequency).

The motif baseline counts amino acid frequencies at each flanking position
relative to training-positive cysteines. A query cysteine is scored by the
average frequency of its flanking residues — this is the simplest conceivable
ranking and must be dominated by any model that extracts non-trivial signal.

All statistics (frequencies, normalisation) are computed from the train split
only. Validation and test cysteine scores are predictions — the model never
sees their labels or sequences during fitting.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

from plantpersulf.features.structure import (
    CysStructureFeature,
    read_structure_features_tsv,
)
from plantpersulf.provenance.hashing import hash_file

AA = "ACDEFGHIKLMNPQRSTVWY"
PAD_CHAR = "X"


SPLIT_FIELDS = (
    "protein_accession",
    "cys_position",
    "label",
    "split",
)


@dataclass(frozen=True)
class CysteineMotifFrequency:
    """Position-weight matrix of AA frequencies around train-positive cysteines."""

    window_radius: int
    central_residue_freq: dict[str, float]
    proteome_hash: str

    def score(self, flanking_window: str) -> float:
        """Average frequency of each observed residue at its position.

        Returns -inf (no score) when the window length mismatches.
        """
        expected = 2 * self.window_radius + 1
        if len(flanking_window) != expected:
            return float("-inf")
        total = 0.0
        n = 0
        for _pos, aa in enumerate(flanking_window):
            freq = self.central_residue_freq.get(aa, 0.0)
            total += freq
            n += 1
        return total / n if n else float("-inf")


@dataclass(frozen=True)
class RankingResult:
    split_name: str
    total_sites: int
    positives_found: int | None  # None if no positives in split
    recall_at_k: dict[int, float] | None


def _read_splits(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=",")
        fieldnames = tuple(reader.fieldnames or ())
        if fieldnames != SPLIT_FIELDS:
            raise RuntimeError(f"split file has invalid columns: {path}")
        return [dict(row) for row in reader]


def _load_proteome(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    cur_header = ""
    cur_lines: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if cur_header:
                acc = cur_header.strip().split("|")[1]
                sequences[acc] = "".join(cur_lines)
            cur_header = line
            cur_lines = []
        elif line:
            cur_lines.append(line)
    if cur_header:
        acc = cur_header.strip().split("|")[1]
        sequences[acc] = "".join(cur_lines)
    return sequences


def _flanking(seq: str, pos: int, radius: int) -> str:
    left_pad = max(0, radius - (pos - 1))
    right_pad = max(0, radius - (len(seq) - pos))
    start = max(0, pos - 1 - radius)
    end = min(len(seq), pos + radius)
    return PAD_CHAR * left_pad + seq[start:end] + PAD_CHAR * right_pad


def train_motif_baseline(
    splits_path: Path,
    proteome_path: Path,
    window_radius: int = 10,
) -> CysteineMotifFrequency:
    """Fit a position-wise AA frequency model on train-positive cysteines."""
    splits = _read_splits(splits_path)
    proteome = _load_proteome(proteome_path)
    proteome_hash = hash_file(proteome_path, "sha256")

    positions = 2 * window_radius + 1
    counts: list[dict[str, int]] = [{} for _ in range(positions)]
    total = 0
    for row in splits:
        if row["split"] != "train" or row["label"] != "positive":
            continue
        seq = proteome.get(row["protein_accession"])
        if seq is None:
            continue
        pos = int(row["cys_position"])
        if pos < 1 or pos > len(seq) or seq[pos - 1] != "C":
            continue
        window = _flanking(seq, pos, window_radius)
        for i, aa in enumerate(window):
            counts[i][aa] = counts[i].get(aa, 0) + 1
        total += 1

    # Laplace-smoothed central residue frequency (worst-case: uniform)
    central_freq: dict[str, float] = {}
    for aa in AA + PAD_CHAR:
        central_freq[aa] = 0.0
    if total > 0:
        for aa in central_freq:
            pos_counts = [c.get(aa, 0) for c in counts]
            central_freq[aa] = math.fsum(pos_counts) / (total * positions)
    else:
        # No train positives → uniform prior
        for aa in central_freq:
            central_freq[aa] = 1.0 / len(central_freq)

    return CysteineMotifFrequency(
        window_radius=window_radius,
        central_residue_freq=central_freq,
        proteome_hash=proteome_hash,
    )


def evaluate_ranking(
    model: CysteineMotifFrequency,
    splits_path: Path,
    proteome_path: Path,
    split_name: str,
) -> RankingResult:
    """Score all cysteines in a split and report positive recovery."""
    splits = _read_splits(splits_path)
    proteome = _load_proteome(proteome_path)

    scored: list[tuple[float, str]] = []
    total = 0
    pos_count = 0
    for row in splits:
        if row["split"] != split_name:
            continue
        seq = proteome.get(row["protein_accession"])
        if seq is None:
            continue
        pos = int(row["cys_position"])
        if pos < 1 or pos > len(seq) or seq[pos - 1] != "C":
            continue
        window = _flanking(seq, pos, model.window_radius)
        score = model.score(window)
        scored.append((score, row["label"]))
        total += 1
        if row["label"] == "positive":
            pos_count += 1

    if not scored:
        return RankingResult(
            split_name=split_name,
            total_sites=0,
            positives_found=None,
            recall_at_k=None,
        )
    scored.sort(key=lambda x: -x[0])
    found = 0
    recall: dict[int, float] = {}
    for k in (1, 5, 10, 25, 50, 100):
        if k > len(scored):
            recall[k] = 1.0 if pos_count > 0 else 0.0
            continue
        hits = sum(1 for _, label in scored[:k] if label == "positive")
        recall[k] = hits / pos_count if pos_count else 0.0
    for _, label in scored:
        if label == "positive":
            found += 1

    return RankingResult(
        split_name=split_name,
        total_sites=total,
        positives_found=found if pos_count > 0 else None,
        recall_at_k=recall if pos_count > 0 else None,
    )


# ---------------------------------------------------------------------------
# Accessibility (contact-number-proxy) no-learning baseline.
#
# Scores a cysteine by how solvent-exposed its structural neighbourhood
# looks, using the ``contact_number_proxy`` field produced by
# ``plantpersulf.features.structure`` (a documented contact-count proxy for
# accessibility, not a rigorous SASA/RSA calculation — see that module's
# docstring). The literature this project builds on treats solvent-exposed
# cysteines as more reactive to H2S-derived persulfidation signalling, so a
# *lower* contact-number proxy (fewer nearby residues, i.e. more exposed)
# must score *higher*. Standardisation statistics are fit on train-split,
# structure-observed cysteines only — exactly the same train-fold-only
# discipline as ``train_motif_baseline`` above. Cysteines with no structure
# are never imputed with the train mean; they always receive the same fixed
# floor score, so they rank below every structure-observed cysteine instead
# of landing arbitrarily in the middle of the ranking.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class AccessibilityBaseline:
    """Z-scored contact-number-proxy ranking model, fit on train only."""

    train_mean: float
    train_std: float
    direction: float = -1.0  # lower contact proxy (more exposed) -> higher score
    missing_score: float = float("-inf")

    def score(self, feature: CysStructureFeature) -> float:
        """Score a Cys; structure-missing Cys always get the fixed floor."""
        if not feature.has_structure or feature.contact_number_proxy is None:
            return self.missing_score
        z = (feature.contact_number_proxy - self.train_mean) / self.train_std
        return self.direction * z


def _index_structure_features(
    structure_features_path: Path,
) -> dict[tuple[str, int], CysStructureFeature]:
    rows = read_structure_features_tsv(structure_features_path)
    return {(row.protein_accession, row.cys_position): row for row in rows}


def train_accessibility_baseline(
    splits_path: Path,
    structure_features_path: Path,
) -> AccessibilityBaseline:
    """Fit contact-number-proxy mean/std on train, structure-observed Cys."""
    splits = _read_splits(splits_path)
    features_by_key = _index_structure_features(structure_features_path)

    train_values: list[float] = []
    for row in splits:
        if row["split"] != "train":
            continue
        key = (row["protein_accession"], int(row["cys_position"]))
        feature = features_by_key.get(key)
        if (
            feature is None
            or not feature.has_structure
            or feature.contact_number_proxy is None
        ):
            continue
        train_values.append(feature.contact_number_proxy)

    if not train_values:
        # No train-fold structure observations: a neutral, non-informative
        # prior (mean 0, unit spread) rather than a division by zero.
        return AccessibilityBaseline(train_mean=0.0, train_std=1.0)

    mean = math.fsum(train_values) / len(train_values)
    variance = math.fsum((v - mean) ** 2 for v in train_values) / len(train_values)
    std = math.sqrt(variance) if variance > 0 else 1.0
    return AccessibilityBaseline(train_mean=mean, train_std=std)


def evaluate_accessibility_ranking(
    model: AccessibilityBaseline,
    splits_path: Path,
    structure_features_path: Path,
    split_name: str,
) -> RankingResult:
    """Score all cysteines in a split by accessibility and report recovery."""
    splits = _read_splits(splits_path)
    features_by_key = _index_structure_features(structure_features_path)

    scored: list[tuple[float, str]] = []
    total = 0
    pos_count = 0
    for row in splits:
        if row["split"] != split_name:
            continue
        key = (row["protein_accession"], int(row["cys_position"]))
        feature = features_by_key.get(key)
        if feature is None:
            continue
        score = model.score(feature)
        scored.append((score, row["label"]))
        total += 1
        if row["label"] == "positive":
            pos_count += 1

    if not scored:
        return RankingResult(
            split_name=split_name,
            total_sites=0,
            positives_found=None,
            recall_at_k=None,
        )
    scored.sort(key=lambda x: -x[0])
    found = sum(1 for _, label in scored if label == "positive")
    recall: dict[int, float] = {}
    for k in (1, 5, 10, 25, 50, 100):
        if k > len(scored):
            recall[k] = 1.0 if pos_count > 0 else 0.0
            continue
        hits = sum(1 for _, label in scored[:k] if label == "positive")
        recall[k] = hits / pos_count if pos_count else 0.0

    return RankingResult(
        split_name=split_name,
        total_sites=total,
        positives_found=found if pos_count > 0 else None,
        recall_at_k=recall if pos_count > 0 else None,
    )
