"""Traceable materialization of registered reference-proteome revisions."""

from __future__ import annotations

import hashlib
import json
import tempfile
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class RegisteredFastaInput:
    path: Path
    sha256: str


@dataclass(frozen=True)
class MaterializedReferenceProteome:
    output_path: Path
    output_sha256: str
    protein_count: int


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _records(path: Path) -> Iterable[tuple[str, str, str]]:
    header: str | None = None
    sequence: list[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        if raw_line.startswith(">"):
            if header is not None:
                yield header, header.split("|")[1], "".join(sequence)
            header = raw_line[1:]
            sequence = []
        elif raw_line:
            if header is None:
                raise RuntimeError(f"FASTA sequence without header: {path}")
            sequence.append(raw_line)
    if header is not None:
        yield header, header.split("|")[1], "".join(sequence)


def _verify_input(source: RegisteredFastaInput) -> None:
    if not source.path.is_file():
        raise FileNotFoundError(source.path)
    if len(source.sha256) != 64 or _sha256(source.path) != source.sha256:
        raise RuntimeError(f"registered source SHA256 mismatch: {source.path}")


def materialize_augmented_reference_proteome(
    base: RegisteredFastaInput,
    additions: tuple[RegisteredFastaInput, ...],
    output_path: Path,
    manifest_path: Path,
) -> MaterializedReferenceProteome:
    """Combine hash-verified FASTA sources without replacing an existing revision."""
    if not additions:
        raise ValueError("at least one registered addition is required")
    if output_path.exists() or manifest_path.exists():
        raise FileExistsError("refusing to overwrite reference-proteome revision")
    inputs = (base, *additions)
    for source in inputs:
        _verify_input(source)
    records: list[tuple[str, str]] = []
    seen: set[str] = set()
    for source in inputs:
        for header, accession, sequence in _records(source.path):
            if not accession or not sequence:
                raise RuntimeError(f"invalid FASTA record: {source.path}")
            if accession in seen:
                raise RuntimeError(f"duplicate protein accession: {accession}")
            seen.add(accession)
            records.append((header, sequence))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", newline="\n", dir=output_path.parent, delete=False
    ) as handle:
        temporary = Path(handle.name)
        for header, sequence in records:
            handle.write(f">{header}\n{sequence}\n")
    temporary.replace(output_path)
    output_sha256 = _sha256(output_path)
    manifest = {
        "input_sha256": {str(source.path): source.sha256 for source in inputs},
        "output_path": str(output_path),
        "output_sha256": output_sha256,
        "protein_count": len(records),
    }
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", newline="\n", dir=manifest_path.parent, delete=False
    ) as handle:
        temporary = Path(handle.name)
        handle.write(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    temporary.replace(manifest_path)
    return MaterializedReferenceProteome(output_path, output_sha256, len(records))
