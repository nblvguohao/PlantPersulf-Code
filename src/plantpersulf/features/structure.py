r"""Task 7 — structure features from real AlphaFold DB PDB files.

Parses the standard, publicly documented fixed-column PDB ``ATOM`` record
format (wwPDB format v3.3, section 9): columns 7-11 atom serial, 13-16 atom
name, 18-20 residue name, 22 chain ID, 23-26 residue sequence number,
31-38/39-46/47-54 x/y/z coordinates, 61-66 temperature factor. AlphaFold DB
stores per-residue pLDDT confidence in the temperature-factor (B-factor)
column, replicated across every atom of a residue; this module reads it
from each residue's C-alpha atom.

Solvent accessibility: this module does **not** implement a rigorous
Shrake-Rupley SASA calculation. It instead computes a documented,
explicitly-labelled proxy — ``contact_number_proxy`` — the count of other
residues' C-alpha atoms within a fixed radius of a residue's own C-alpha
atom. A *low* contact number indicates a residue with few nearby residues,
which is a common (imperfect) correlate of higher solvent exposure; a
*high* contact number correlates with burial. This is a simplification and
must never be reported or consumed as ``relative_solvent_accessibility``.

Missing structures (no AlphaFold model, a UniProt isoform accession, or a
Cys position that fails to map onto the parsed residue chain) are never
imputed with an average value. They are represented as an explicit
``has_structure=False`` mask with every structural field set to ``None``,
per the project's structure-missingness rule
(docs/PlantPersulf_Code_TDD_Codex.md, Task 7).
"""

from __future__ import annotations

import csv
import math
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

LOW_PLDDT_THRESHOLD = 50.0
DEFAULT_CONTACT_RADIUS_ANGSTROM = 12.0

STRUCTURE_FEATURE_FIELDS = (
    "protein_accession",
    "cys_position",
    "has_structure",
    "plddt",
    "low_plddt",
    "contact_number_proxy",
)


@dataclass(frozen=True)
class ResidueStructure:
    chain_id: str
    res_seq: int
    res_name: str
    plddt: float
    x: float
    y: float
    z: float


@dataclass(frozen=True)
class CysStructureFeature:
    protein_accession: str
    cys_position: int
    has_structure: bool
    plddt: float | None
    low_plddt: bool | None
    contact_number_proxy: float | None


def parse_alphafold_pdb(pdb_text: str) -> tuple[ResidueStructure, ...]:
    """Parse C-alpha atoms from standard fixed-column PDB ATOM records.

    Non-ATOM lines and non-CA atoms are ignored. Malformed ATOM lines
    (wrong length, non-numeric coordinate/B-factor fields) are skipped
    rather than crashing the whole file, mirroring the project's
    per-record fault tolerance elsewhere (e.g. peptide parsers).
    """
    residues: list[ResidueStructure] = []
    for line in pdb_text.splitlines():
        if len(line) < 66 or line[0:6] != "ATOM  ":
            continue
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue
        res_name = line[17:20].strip()
        chain_id = line[21]
        try:
            res_seq = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            plddt = float(line[60:66])
        except ValueError:
            continue
        residues.append(
            ResidueStructure(
                chain_id=chain_id,
                res_seq=res_seq,
                res_name=res_name,
                plddt=plddt,
                x=x,
                y=y,
                z=z,
            )
        )
    return tuple(residues)


def missing_structure_feature(
    protein_accession: str,
    cys_position: int,
) -> CysStructureFeature:
    """Return the explicit missing-structure mask (no field is imputed)."""
    return CysStructureFeature(
        protein_accession=protein_accession,
        cys_position=cys_position,
        has_structure=False,
        plddt=None,
        low_plddt=None,
        contact_number_proxy=None,
    )


def _contact_number_proxy(
    residues: Sequence[ResidueStructure],
    center_index: int,
    radius_angstrom: float,
) -> float:
    center = residues[center_index]
    count = 0
    for i, residue in enumerate(residues):
        if i == center_index:
            continue
        distance = math.sqrt(
            (residue.x - center.x) ** 2
            + (residue.y - center.y) ** 2
            + (residue.z - center.z) ** 2
        )
        if distance <= radius_angstrom:
            count += 1
    return float(count)


def extract_cys_structure_features(
    pdb_text: str,
    protein_accession: str,
    cys_positions: Sequence[int],
    contact_radius_angstrom: float = DEFAULT_CONTACT_RADIUS_ANGSTROM,
) -> list[CysStructureFeature]:
    """Extract pLDDT and the contact-number accessibility proxy per Cys.

    A requested position that is absent from the parsed chain, or whose
    residue name is not CYS (a coordinate-mapping failure), is returned as
    an explicit missing-structure record — it is dropped from the
    structure model but the caller is expected to retain it in any
    sequence-only model, per the Task 7 spec.
    """
    residues = parse_alphafold_pdb(pdb_text)
    by_seq = {residue.res_seq: (i, residue) for i, residue in enumerate(residues)}
    features: list[CysStructureFeature] = []
    for position in cys_positions:
        match = by_seq.get(position)
        if match is None or match[1].res_name != "CYS":
            features.append(missing_structure_feature(protein_accession, position))
            continue
        index, residue = match
        contact = _contact_number_proxy(residues, index, contact_radius_angstrom)
        features.append(
            CysStructureFeature(
                protein_accession=protein_accession,
                cys_position=position,
                has_structure=True,
                plddt=residue.plddt,
                low_plddt=residue.plddt < LOW_PLDDT_THRESHOLD,
                contact_number_proxy=contact,
            )
        )
    return features


def _format_optional(value: float | bool | None) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    return repr(value)


def write_structure_features_tsv(
    rows: Sequence[CysStructureFeature],
    path: Path,
) -> None:
    """Write structure features with an explicit missingness mask column."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = path.with_suffix(".tsv.tmp")
    with temporary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(STRUCTURE_FEATURE_FIELDS)
        for row in rows:
            writer.writerow(
                [
                    row.protein_accession,
                    str(row.cys_position),
                    _format_optional(row.has_structure),
                    _format_optional(row.plddt),
                    _format_optional(row.low_plddt),
                    _format_optional(row.contact_number_proxy),
                ]
            )
    temporary_path.replace(path)


def _parse_optional_float(value: str) -> float | None:
    return float(value) if value else None


def _parse_optional_bool(value: str) -> bool | None:
    if value == "":
        return None
    if value not in ("true", "false"):
        raise RuntimeError(f"structure feature has invalid boolean: {value!r}")
    return value == "true"


def read_structure_features_tsv(path: Path) -> list[CysStructureFeature]:
    """Read structure features back, preserving the missingness mask."""
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != STRUCTURE_FEATURE_FIELDS:
            raise RuntimeError(f"structure feature file has invalid columns: {path}")
        rows: list[CysStructureFeature] = []
        for record in reader:
            has_structure = _parse_optional_bool(record["has_structure"])
            if has_structure is None:
                raise RuntimeError(
                    f"structure feature row has no has_structure flag: {path}"
                )
            rows.append(
                CysStructureFeature(
                    protein_accession=record["protein_accession"],
                    cys_position=int(record["cys_position"]),
                    has_structure=has_structure,
                    plddt=_parse_optional_float(record["plddt"]),
                    low_plddt=_parse_optional_bool(record["low_plddt"]),
                    contact_number_proxy=_parse_optional_float(
                        record["contact_number_proxy"]
                    ),
                )
            )
    return rows
