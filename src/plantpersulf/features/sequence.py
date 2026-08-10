"""Task 7 — deterministic sequence-based cysteine features.

Extracts per-cysteine flanking windows, amino acid composition, and
physicochemical properties from the SHA256-pinned reference proteome,
keyed to the benchmark labels. No external model or network call.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

POSITIVE_CHARGE_RESIDUES: frozenset[str] = frozenset({"K", "R"})
"""Lys/Arg only — fully protonated/positive at physiological pH. His (pKa
~6) is deliberately excluded: it is only partially protonated at pH 7 and
its inclusion would weaken rather than sharpen the thiolate-stabilization
signal this feature targets (COPLBI-D-26-00068 review, Figure 1C:
persulfidation requires nucleophilic attack by the anionic thiolate form of
Cys, which nearby positive charge favours by lowering the local thiol
pKa)."""

KYTE_DOOLITTLE: dict[str, float] = {
    "A": 1.8,
    "C": 2.5,
    "D": -3.5,
    "E": -3.5,
    "F": 2.8,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "K": -3.9,
    "L": 3.8,
    "M": 1.9,
    "N": -3.5,
    "P": -1.6,
    "Q": -3.5,
    "R": -4.5,
    "S": -0.8,
    "T": -0.7,
    "V": 4.2,
    "W": -0.9,
    "Y": -1.3,
    "X": 0.0,
}

BENCHMARK_FIELDS = (
    "protein_accession",
    "cys_position_in_protein",
    "label",
    "study_accession",
    "evidence_level",
    "source_sha256",
)


@dataclass(frozen=True)
class SequenceFeatureRow:
    protein_accession: str
    cys_position: int
    label: str
    flanking_window: str
    hydrophobicity: float
    cys_density: float
    protein_length: int
    local_positive_charge_density: float


def _load_labels(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != BENCHMARK_FIELDS:
            raise RuntimeError(f"benchmark labels have invalid columns: {path}")
        return [dict(row) for row in reader]


def _header_accession(header: str) -> str:
    """Accession from a fasta header: UniProt-style ``>sp|ACC|...`` takes the
    second pipe field; anything else (e.g. EnsemblFungi ``>MGG_07573T0 pep
    chromosome:...``) takes the first whitespace-separated token."""
    parts = header.strip().split("|")
    if len(parts) >= 2 and parts[1]:
        return parts[1]
    return header.strip().lstrip(">").split()[0]


def _load_proteome(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    cur_header = ""
    cur_lines: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if cur_header:
                sequences[_header_accession(cur_header)] = "".join(cur_lines)
            cur_header = line
            cur_lines = []
        elif line:
            cur_lines.append(line)
    if cur_header:
        sequences[_header_accession(cur_header)] = "".join(cur_lines)
    return sequences


def _flanking_window(
    sequence: str,
    position: int,
    radius: int,
) -> str:
    """Extract padded flanking window centred on the modified residue."""
    left_pad = max(0, radius - (position - 1))
    right_pad = max(0, radius - (len(sequence) - position))
    start = max(0, position - 1 - radius)
    end = min(len(sequence), position + radius)
    return "X" * left_pad + sequence[start:end] + "X" * right_pad


def _hydrophobicity(window: str) -> float:
    """Average Kyte-Doolittle hydrophobicity over the window."""
    values = [KYTE_DOOLITTLE.get(aa, 0.0) for aa in window]
    return sum(values) / len(values) if values else 0.0


def _local_positive_charge_density(window: str) -> float:
    """Fraction of the flanking window that is Lys/Arg — a thiolate-
    stabilization proxy (see ``POSITIVE_CHARGE_RESIDUES`` docstring).
    Padding 'X' residues count toward the denominator (matching how
    ``_hydrophobicity`` treats them via ``KYTE_DOOLITTLE['X'] = 0.0``) so a
    site near a sequence terminus is not artificially inflated."""
    if not window:
        return 0.0
    charged = sum(1 for aa in window if aa in POSITIVE_CHARGE_RESIDUES)
    return charged / len(window)


def extract_sequence_features(
    labels_path: Path,
    proteome_path: Path,
    window_radius: int = 20,
) -> tuple[SequenceFeatureRow, ...]:
    """Extract per-cysteine sequence features from the benchmark and proteome."""
    labels = _load_labels(labels_path)
    proteome = _load_proteome(proteome_path)

    rows: list[SequenceFeatureRow] = []
    for row in labels:
        protein = row["protein_accession"]
        cys_pos = int(row["cys_position_in_protein"])
        seq = proteome.get(protein)
        if seq is None:
            continue
        if cys_pos < 1 or cys_pos > len(seq) or seq[cys_pos - 1] != "C":
            continue
        flank = _flanking_window(seq, cys_pos, window_radius)
        hydro = _hydrophobicity(flank)
        cys_count = seq.count("C")
        cys_density = cys_count / len(seq) if seq else 0.0
        charge_density = _local_positive_charge_density(flank)
        rows.append(
            SequenceFeatureRow(
                protein_accession=protein,
                cys_position=cys_pos,
                label=row["label"],
                flanking_window=flank,
                hydrophobicity=hydro,
                cys_density=cys_density,
                protein_length=len(seq),
                local_positive_charge_density=charge_density,
            )
        )
    return tuple(rows)
