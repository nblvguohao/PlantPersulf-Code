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


def parse_mmcif(text: str) -> tuple[ResidueStructure, ...]:
    """Parse C-alpha atoms from an mmCIF ``_atom_site`` loop.

    AlphaFold Server produces ModelCIF/mmCIF (``.cif``) which stores
    per-residue pLDDT in the ``B_iso_or_equiv`` column, same semantics as
    the PDB B-factor column. Only CA atoms are retained.
    """
    lines = text.splitlines()
    col_start: int | None = None
    col_end: int | None = None
    for i, line in enumerate(lines):
        if line.startswith("_atom_site.") and col_start is None:
            col_start = i
        if _end_of_atom_site(col_start, i, line):
            col_end = i
            break
    if col_start is None or col_end is None:
        return ()
    cols = [
        ln.split(".")[1].strip()
        for ln in lines[col_start:col_end]
        if ln.startswith("_atom_site.")
    ]

    def _ci(name: str) -> int | None:
        for j, c in enumerate(cols):
            if c == name:
                return j
        return None

    i_seq = _ci("label_seq_id")
    i_comp = _ci("label_comp_id")
    i_atom = _ci("label_atom_id")
    i_chain = _ci("label_asym_id")
    i_B = _ci("B_iso_or_equiv")
    i_x = _ci("Cartn_x")
    i_y = _ci("Cartn_y")
    i_z = _ci("Cartn_z")
    if None in (i_seq, i_comp, i_atom, i_B, i_x, i_y, i_z):
        return ()
    assert i_seq is not None and i_comp is not None and i_atom is not None
    assert i_B is not None and i_x is not None and i_y is not None and i_z is not None

    row_start = col_end
    while row_start < len(lines) and (
        lines[row_start].startswith("#") or lines[row_start].strip() == ""
    ):
        row_start += 1
    row_end = row_start
    while (
        row_end < len(lines)
        and not lines[row_end].startswith("#")
        and lines[row_end].strip()
    ):
        row_end += 1

    residues: list[ResidueStructure] = []
    for line in lines[row_start:row_end]:
        parts = line.split()
        try:
            if parts[i_atom] != "CA":
                continue
            chain = parts[i_chain] if i_chain is not None else "A"
            residues.append(
                ResidueStructure(
                    chain_id=chain,
                    res_seq=int(parts[i_seq]),
                    res_name=parts[i_comp],
                    plddt=float(parts[i_B]),
                    x=float(parts[i_x]),
                    y=float(parts[i_y]),
                    z=float(parts[i_z]),
                )
            )
        except (ValueError, IndexError):
            continue
    return tuple(residues)


def _end_of_atom_site(
    col_start: int | None,
    i: int,
    line: str,
) -> bool:
    return (
        col_start is not None
        and i > col_start
        and bool(line.strip())
        and not line.startswith(("_atom_site.", "#"))
    )


def _parse_structure_text(
    text: str,
) -> tuple[ResidueStructure, ...]:
    """Auto-detect PDB vs mmCIF and parse. Prefers PDB if both signatures
    are present; falls back to mmCIF."""
    residues = parse_alphafold_pdb(text)
    if residues:
        return residues
    return parse_mmcif(text)


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
    residues = _parse_structure_text(pdb_text)
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
