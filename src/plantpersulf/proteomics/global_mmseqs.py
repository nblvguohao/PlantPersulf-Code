"""Preparation and execution helpers for traceable global MMseqs2 clusters."""

from __future__ import annotations

import hashlib
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from plantpersulf.benchmark.multispecies_splits import (
    GlobalClusterRow,
    global_protein_id,
)


@dataclass(frozen=True)
class ReferenceProteomeInput:
    species: str
    path: Path
    sha256: str


@dataclass(frozen=True)
class NamespacedProtein:
    species: str
    protein_accession: str
    global_protein_id: str
    source_proteome_sha256: str


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _parse_fasta(path: Path) -> Iterable[tuple[str, str]]:
    header: str | None = None
    sequence: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)
                header = line[1:]
                sequence = []
            else:
                if header is None:
                    raise RuntimeError(f"FASTA sequence without header: {path}")
                sequence.append(line)
    if header is not None:
        yield header, "".join(sequence)


def _accession(header: str) -> str:
    fields = header.split("|")
    if len(fields) > 1 and fields[1].strip():
        return fields[1].strip()
    return header.split()[0]


def prepare_namespaced_fasta(
    proteomes: Iterable[ReferenceProteomeInput], output_path: Path
) -> tuple[NamespacedProtein, ...]:
    """Create a header-only derived FASTA after exact source-hash checks.

    Sequence values are copied byte-for-byte at the residue level; only FASTA
    identifiers are rewritten to make cross-species MMseqs2 output reversible.
    """
    inputs = tuple(proteomes)
    if not inputs:
        raise RuntimeError("at least one registered reference proteome is required")
    if output_path.exists():
        raise FileExistsError(f"refusing to overwrite namespaced FASTA: {output_path}")
    rows: list[NamespacedProtein] = []
    records: list[tuple[str, str]] = []
    seen_species: set[str] = set()
    seen_global_ids: set[str] = set()
    for source in inputs:
        if not source.species:
            raise ValueError("reference proteome species is required")
        if source.species in seen_species:
            raise RuntimeError(
                f"duplicate reference proteome species: {source.species}"
            )
        seen_species.add(source.species)
        if not source.path.is_file():
            raise FileNotFoundError(f"reference proteome not found: {source.path}")
        observed = _sha256(source.path)
        if observed.lower() != source.sha256.lower():
            raise RuntimeError(
                f"reference proteome SHA256 mismatch for {source.species}: "
                f"expected {source.sha256}, observed {observed}"
            )
        found = 0
        for header, sequence in _parse_fasta(source.path):
            accession = _accession(header)
            key = global_protein_id(source.species, accession)
            if key in seen_global_ids:
                raise RuntimeError(f"duplicate global protein id: {key}")
            if not sequence:
                raise RuntimeError(f"empty protein sequence: {key}")
            seen_global_ids.add(key)
            rows.append(
                NamespacedProtein(
                    species=source.species,
                    protein_accession=accession,
                    global_protein_id=key,
                    source_proteome_sha256=observed,
                )
            )
            records.append((key, sequence))
            found += 1
        if not found:
            raise RuntimeError(f"reference proteome contains no records: {source.path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        for key, sequence in records:
            handle.write(f">{key}\n{sequence}\n")
    return tuple(rows)


def global_cluster_rows_from_mmseqs_pairs(
    proteins: Iterable[NamespacedProtein],
    pairs_path: Path,
    identity_threshold: float,
) -> tuple[GlobalClusterRow, ...]:
    """Translate MMseqs2 representative/member pairs into the frozen schema."""
    if not pairs_path.is_file():
        raise FileNotFoundError(f"MMseqs2 cluster pairs not found: {pairs_path}")
    materialized = tuple(proteins)
    known = {row.global_protein_id: row for row in materialized}
    if len(known) != len(materialized):
        raise RuntimeError("namespaced protein table contains duplicate identifiers")
    representative_by_member: dict[str, str] = {}
    with pairs_path.open(encoding="utf-8", newline="") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) != 2:
                raise RuntimeError(
                    f"invalid MMseqs2 pair at line {line_number}: expected two columns"
                )
            representative, member = fields
            if representative not in known or member not in known:
                raise RuntimeError(
                    f"MMseqs2 pair references unknown protein at line {line_number}"
                )
            if member in representative_by_member:
                raise RuntimeError(f"duplicate MMseqs2 member: {member}")
            representative_by_member[member] = representative
    missing = sorted(set(known) - set(representative_by_member))
    if missing:
        raise RuntimeError(f"MMseqs2 cluster pairs missing proteins: {missing[:5]}")
    return tuple(
        GlobalClusterRow(
            species=protein.species,
            protein_accession=protein.protein_accession,
            global_protein_id=protein.global_protein_id,
            cluster_id=representative_by_member[protein.global_protein_id],
            identity_threshold=identity_threshold,
            source_proteome_sha256=protein.source_proteome_sha256,
        )
        for protein in sorted(materialized, key=lambda item: item.global_protein_id)
    )


def easy_cluster_command(
    mmseqs_binary: str,
    input_fasta: Path,
    output_prefix: Path,
    temporary_directory: Path,
    *,
    identity_threshold: float,
    threads: int = 8,
) -> list[str]:
    """Return the fully pinned MMseqs2 easy-cluster invocation."""
    if identity_threshold not in {0.2, 0.3, 0.4}:
        raise ValueError("identity_threshold must be one of 0.2, 0.3, 0.4")
    if threads < 1:
        raise ValueError("threads must be >= 1")
    return [
        mmseqs_binary,
        "easy-cluster",
        str(input_fasta),
        str(output_prefix),
        str(temporary_directory),
        "--min-seq-id",
        f"{identity_threshold:g}",
        "-c",
        "0.5",
        "--cov-mode",
        "0",
        "--threads",
        str(threads),
    ]
