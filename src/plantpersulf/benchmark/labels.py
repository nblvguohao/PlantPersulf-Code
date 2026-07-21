"""Task 5 — positive–unlabeled persulfidation benchmark v1.

Consumes the deterministic persulfidation site outputs from Phase A and a
SHA256-pinned Arabidopsis reference proteome. Every cysteine receives exactly
one label: ``positive`` if the (protein, Cys-position) pair appears in any
published site-level persulfidation evidence, ``unlabeled`` otherwise. The
label ``negative`` is never emitted.

Output is written atomically under ``data/processed/benchmark_v1/`` with a
manifest recording every input SHA256 and output file hash.
"""

from __future__ import annotations

import json
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path

from plantpersulf.proteomics.persulfidation_publish import (
    PERSULF_SITE_FIELDS,
    audit_persulfidation_sites,
)
from plantpersulf.provenance.hashing import hash_file

BENCHMARK_VERSION = 1
BENCHMARK_SITE_FIELDS = (
    "protein_accession",
    "cys_position_in_protein",
    "label",
    "study_accession",
    "evidence_level",
    "source_sha256",
)
OUTPUT_NAMES = (
    "sites.tsv",
    "manifest.json",
)
STUDIES = ("PXD006140", "PXD024061")


@dataclass(frozen=True)
class LabeledSite:
    protein_accession: str
    cys_position_in_protein: int
    label: str
    study_accession: str
    evidence_level: str
    source_sha256: str


class BenchmarkLabels:
    """Iterable, immutable collection of labeled cysteine sites."""

    def __init__(
        self,
        sites: tuple[LabeledSite, ...],
        positive_count: int,
        unlabeled_count: int,
    ) -> None:
        self.sites = sites
        self.positive_count = positive_count
        self.unlabeled_count = unlabeled_count
        self.negative_count = 0

    def __len__(self) -> int:
        return len(self.sites)

    def __iter__(self):  # type: ignore[no-untyped-def]
        return iter(self.sites)

    def __getitem__(self, index: int) -> LabeledSite:
        return self.sites[index]


def _read_tsv(path: Path, fields: tuple[str, ...]) -> list[dict[str, str]]:
    import csv

    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != fields:
            raise RuntimeError(f"benchmark source has invalid columns: {path}")
        rows = [dict(row) for row in reader]
    if any(
        None in row or any(value is None for value in row.values()) for row in rows
    ):
        raise RuntimeError(f"benchmark source has malformed row: {path}")
    return rows


def _write_tsv(
    path: Path,
    records: list[dict[str, str]],
) -> None:
    import csv

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=BENCHMARK_SITE_FIELDS,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(records)


def _tree_hashes(root: Path) -> dict[str, str]:
    return {
        path.relative_to(root).as_posix(): hash_file(path, "sha256")
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


def _parse_fasta(path: Path) -> dict[str, str]:
    """Return {accession: sequence} from a multi-record FASTA."""
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


def build_benchmark_labels(
    site_output_root: Path,
    proteome_path: Path,
    output_directory: Path,
) -> BenchmarkLabels:
    """Build the PU benchmark from published site outputs and reference proteome."""
    proteome_sha256 = hash_file(proteome_path, "sha256")

    # Collect positives from both studies
    positives: dict[tuple[str, int], dict[str, str]] = {}
    for accession in STUDIES:
        site_dir = site_output_root / accession / "persulfidation_sites_v1"
        audit_persulfidation_sites(accession, site_dir)
        site_rows = _read_tsv(site_dir / "sites.tsv", PERSULF_SITE_FIELDS)
        for row in site_rows:
            if row["evidence_level"] != "site_ms":
                continue
            key = (row["protein_accession_raw"], int(row["cys_position_in_protein"]))
            positives[key] = {
                "protein_accession": row["protein_accession_raw"],
                "cys_position_in_protein": row["cys_position_in_protein"],
                "label": "positive",
                "study_accession": accession,
                "evidence_level": row["evidence_level"],
                "source_sha256": row["source_sha256"],
            }

    # Build unlabeled set: every Cys in the proteome not already a positive
    sequences = _parse_fasta(proteome_path)
    records: list[dict[str, str]] = []
    unlabeled = 0
    seen_unlabeled: set[tuple[str, int]] = set()
    for acc in sorted(sequences):
        seq = sequences[acc]
        for pos, aa in enumerate(seq, start=1):
            if aa != "C":
                continue
            key = (acc, pos)
            if key in positives:
                records.append(positives[key])
            elif key not in seen_unlabeled:
                records.append(
                    {
                        "protein_accession": acc,
                        "cys_position_in_protein": str(pos),
                        "label": "unlabeled",
                        "study_accession": "",
                        "evidence_level": "",
                        "source_sha256": "",
                    }
                )
                seen_unlabeled.add(key)
                unlabeled += 1

    # Handle any positive whose protein is not in the reference proteome
    for key, row in positives.items():
        existing = {
            (r["protein_accession"], int(r["cys_position_in_protein"]))
            for r in records
        }
        if key not in existing:
            records.append(row)

    labels = BenchmarkLabels(
        sites=tuple(
            LabeledSite(
                protein_accession=r["protein_accession"],
                cys_position_in_protein=int(r["cys_position_in_protein"]),
                label=r["label"],
                study_accession=r["study_accession"],
                evidence_level=r["evidence_level"],
                source_sha256=r["source_sha256"],
            )
            for r in records
        ),
        positive_count=len(positives),
        unlabeled_count=unlabeled,
    )

    # Write atomic, deterministic output
    output_directory.parent.mkdir(parents=True, exist_ok=True)
    temporary_root = Path(
        tempfile.mkdtemp(
            prefix=f".{output_directory.name}-",
            dir=output_directory.parent,
        )
    )
    staged = temporary_root / output_directory.name
    staged.mkdir()
    try:
        _write_tsv(staged / "sites.tsv", records)
        manifest: dict[str, object] = {
            "schema_version": BENCHMARK_VERSION,
            "benchmark_version": BENCHMARK_VERSION,
            "positive_count": labels.positive_count,
            "unlabeled_count": labels.unlabeled_count,
            "negative_count": labels.negative_count,
            "labels_created": True,
            "labels_are_positive_or_unlabeled_only": True,
            "nondetection_labeled_negative": False,
            "biological_values_modified": False,
            "proteome_sha256": proteome_sha256,
            "inputs": [
                {
                    "name": "arabidopsis_reference_proteome",
                    "path": proteome_path.resolve().as_posix(),
                    "sha256": proteome_sha256,
                },
            ],
            "outputs": [
                {
                    "file": "sites.tsv",
                    "sha256": hash_file(staged / "sites.tsv", "sha256"),
                },
            ],
        }
        (staged / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        if output_directory.exists():
            if _tree_hashes(output_directory) != _tree_hashes(staged):
                raise RuntimeError(
                    f"existing benchmark output differs: {output_directory}"
                )
            return labels
        staged.replace(output_directory)
        return labels
    finally:
        shutil.rmtree(temporary_root, ignore_errors=True)
