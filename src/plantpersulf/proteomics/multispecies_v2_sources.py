"""Registered, development-restricted sources for multispecies v2."""

from __future__ import annotations

import csv
from collections.abc import Iterable
from pathlib import Path

from plantpersulf.benchmark.multispecies_splits import MultispeciesSiteRow
from plantpersulf.proteomics.multispecies_v2_dataset import MultispeciesV2SiteRow

FROZEN_DEVELOPMENT_POSITIVE_COLUMNS = (
    "species",
    "protein_accession",
    "cys_position",
    "study_accessions",
    "source_sha256",
)


def _valid_sha256(value: str) -> bool:
    return len(value) == 64 and all(
        char in "0123456789abcdef" for char in value.lower()
    )


def write_frozen_development_positives(
    path: Path,
    rows: Iterable[MultispeciesSiteRow],
    *,
    source_sha256_by_study: dict[str, str],
) -> None:
    """Persist a development-only positive manifest once, merging evidence."""
    if path.exists():
        raise FileExistsError(f"refusing to overwrite frozen manifest: {path}")
    grouped: dict[tuple[str, str, int], set[str]] = {}
    for row in rows:
        if row.label != "positive":
            raise ValueError("frozen development manifest accepts positives only")
        if row.study_accession not in source_sha256_by_study:
            raise RuntimeError(
                f"missing source SHA256 for study: {row.study_accession}"
            )
        grouped.setdefault(
            (row.species, row.protein_accession, row.cys_position), set()
        ).add(row.study_accession)
    if not grouped:
        raise RuntimeError("refusing to write empty development-positive manifest")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=FROZEN_DEVELOPMENT_POSITIVE_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for (species, accession, position), studies_set in sorted(grouped.items()):
            studies = tuple(sorted(studies_set))
            hashes = tuple(
                sorted(
                    {
                        source_hash
                        for study in studies
                        for source_hash in source_sha256_by_study[study].split(";")
                    }
                )
            )
            if any(not _valid_sha256(value) for value in hashes):
                raise RuntimeError("development-positive source SHA256 is invalid")
            writer.writerow(
                {
                    "species": species,
                    "protein_accession": accession,
                    "cys_position": position,
                    "study_accessions": ";".join(studies),
                    "source_sha256": ";".join(hashes),
                }
            )


def load_frozen_development_positives(
    path: Path,
) -> tuple[MultispeciesSiteRow, ...]:
    """Load the only positive artifact accessible to an unlocked dev process."""
    rows: list[MultispeciesSiteRow] = []
    seen: set[tuple[str, str, int, str]] = set()
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != FROZEN_DEVELOPMENT_POSITIVE_COLUMNS:
            raise RuntimeError("development-positive manifest has invalid columns")
        for record in reader:
            source_hashes = record["source_sha256"].split(";")
            if not source_hashes or any(
                not _valid_sha256(value) for value in source_hashes
            ):
                raise RuntimeError("development-positive source SHA256 is invalid")
            studies = tuple(
                study for study in record["study_accessions"].split(";") if study
            )
            if not studies:
                raise RuntimeError("development-positive study accession is missing")
            for study in studies:
                key = (
                    record["species"],
                    record["protein_accession"],
                    int(record["cys_position"]),
                    study,
                )
                if key in seen:
                    raise RuntimeError(
                        f"duplicate development-positive evidence: {key}"
                    )
                seen.add(key)
                rows.append(
                    MultispeciesSiteRow(
                        species=key[0],
                        protein_accession=key[1],
                        cys_position=key[2],
                        label="positive",
                        study_accession=key[3],
                    )
                )
    if not rows:
        raise RuntimeError("development-positive manifest is empty")
    return tuple(rows)


def development_positive_rows(
    *,
    benchmark: Path,
    tomato_xlsx: Path,
    rice_sd01: Path,
    rice_sd04: Path,
    rice_ss_all: Path,
    magnaporthe_tsv: Path,
    proteomes: dict[str, dict[str, str]],
    allowed_proteins: dict[str, set[str]],
) -> tuple[MultispeciesSiteRow, ...]:
    """Parse only positive rows whose proteins are in development scope."""
    rows: list[MultispeciesSiteRow] = []
    with benchmark.open(encoding="utf-8", newline="") as handle:
        for record in csv.DictReader(handle, delimiter="\t"):
            accession = record.get("protein_accession", "")
            if (
                record.get("label") == "positive"
                and accession in allowed_proteins["arabidopsis"]
            ):
                rows.append(
                    MultispeciesSiteRow(
                        "arabidopsis",
                        accession,
                        int(record["cys_position_in_protein"]),
                        "positive",
                        record.get("study_accession", ""),
                    )
                )
    from plantpersulf.proteomics.kiae271_sites import parse_kiae271_sites
    from plantpersulf.proteomics.pxd063170_sites import parse_pxd063170_sites
    from plantpersulf.proteomics.pxd072089_sites import parse_pxd072089_sites

    tomato = parse_kiae271_sites(
        tomato_xlsx, proteomes["tomato"], allowed_accessions=allowed_proteins["tomato"]
    )
    rice = parse_pxd072089_sites(
        rice_sd01,
        rice_sd04,
        proteomes["rice"],
        rice_ss_all,
        allowed_accessions=allowed_proteins["rice"],
    )
    fungus = parse_pxd063170_sites(
        magnaporthe_tsv,
        proteomes["magnaporthe"],
        allowed_accessions=allowed_proteins["magnaporthe"],
    )
    rows.extend(
        MultispeciesSiteRow(
            "tomato",
            s.protein_accession,
            s.cys_position,
            "positive",
            "KIAE271_SUPPL",
        )
        for s in tomato.sites
    )
    rows.extend(
        MultispeciesSiteRow(
            "rice",
            s.protein_accession,
            s.cys_position,
            "positive",
            "PXD072089",
        )
        for s in rice.sites
    )
    rows.extend(
        MultispeciesSiteRow(
            "magnaporthe",
            s.protein_accession,
            s.cys_position,
            "positive",
            "PXD063170",
        )
        for s in fungus.sites
    )
    return tuple(rows)


def reference_cysteines(
    proteomes: dict[str, dict[str, str]], allowed_proteins: dict[str, set[str]]
) -> Iterable[tuple[str, str, int]]:
    for species, proteome in proteomes.items():
        for accession in allowed_proteins[species]:
            sequence = proteome.get(accession)
            if sequence is None:
                raise RuntimeError(
                    f"allowed protein missing from reference: {species}|{accession}"
                )
            for position, residue in enumerate(sequence, start=1):
                if residue == "C":
                    yield species, accession, position


def sequence_feature_map(
    rows: Iterable[MultispeciesV2SiteRow], proteomes: dict[str, dict[str, str]]
) -> dict[tuple[str, int], tuple[float | None, ...]]:
    """Return real reference-sequence features for v2 fit/validation rows."""
    from plantpersulf.features.sequence import (
        _flanking_window,
        _hydrophobicity,
        _local_positive_charge_density,
    )

    result: dict[tuple[str, int], tuple[float | None, ...]] = {}
    for row in rows:
        sequence = proteomes[row.species].get(row.protein_accession)
        if sequence is None or sequence[row.cys_position - 1] != "C":
            raise RuntimeError(f"unverified v2 Cys coordinate: {row.global_protein_id}")
        flank = _flanking_window(sequence, row.cys_position, 10)
        result[(row.global_protein_id, row.cys_position)] = (
            _hydrophobicity(flank),
            sequence.count("C") / len(sequence),
            _local_positive_charge_density(flank),
        )
    return result
