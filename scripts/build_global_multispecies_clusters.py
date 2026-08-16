#!/usr/bin/env python
"""Build registered global 20/30/40% MMseqs2 cluster tables for Task 9.

This command intentionally fails before producing a scientific output when
MMseqs2 is absent, a reference-proteome hash changes, or a cluster export is
incomplete. It never creates singleton fallback clusters.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import tempfile
from pathlib import Path
from typing import Any

import yaml

from plantpersulf.benchmark.multispecies_splits import write_global_cluster_table
from plantpersulf.proteomics.global_mmseqs import (
    ReferenceProteomeInput,
    easy_cluster_command,
    global_cluster_rows_from_mmseqs_pairs,
    prepare_namespaced_fasta,
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_accession_mapping(path: Path, proteins: tuple[Any, ...]) -> None:
    import csv

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "species",
                "protein_accession",
                "global_protein_id",
                "source_proteome_sha256",
            ),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for protein in sorted(proteins, key=lambda row: row.global_protein_id):
            writer.writerow(
                {
                    "species": protein.species,
                    "protein_accession": protein.protein_accession,
                    "global_protein_id": protein.global_protein_id,
                    "source_proteome_sha256": protein.source_proteome_sha256,
                }
            )


def _output_path(primary: Path, threshold: float) -> Path:
    percent = int(round(threshold * 100))
    return primary.with_name(primary.name.replace("30", f"{percent:02d}", 1))


def _provenance_path(primary: Path, suffix: str) -> Path:
    stem = primary.stem.replace("_30_", "_", 1)
    return primary.with_name(f"{stem}.{suffix}")


def build(config_path: Path, mmseqs_binary: str) -> list[Path]:
    cfg: dict[str, Any] = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if cfg.get("version") != 2:
        raise RuntimeError("multispecies v2 config is required")
    mmseqs = cfg["global_mmseqs2"]
    proteomes = tuple(
        ReferenceProteomeInput(
            species=str(row["species"]),
            path=Path(str(row["path"])),
            sha256=str(row["sha256"]),
        )
        for row in cfg["reference_proteomes"]
    )
    primary_output = Path(str(mmseqs["cluster_table"]))
    thresholds = [float(mmseqs["primary_identity_threshold"])] + [
        float(value) for value in mmseqs["sensitivity_identity_thresholds"]
    ]
    if len(set(thresholds)) != 3 or set(thresholds) != {0.2, 0.3, 0.4}:
        raise RuntimeError("global MMseqs2 config must request 20/30/40% clusters")
    outputs = [_output_path(primary_output, threshold) for threshold in thresholds]
    if any(path.exists() for path in outputs):
        existing = next(path for path in outputs if path.exists())
        raise FileExistsError(
            f"refusing to overwrite global cluster output: {existing}"
        )
    mapping_path = _provenance_path(primary_output, "accession_mapping.tsv")
    if mapping_path.exists():
        raise FileExistsError(
            f"refusing to overwrite global accession mapping: {mapping_path}"
        )

    primary_output.parent.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    commands: list[list[str]] = []
    version_result = subprocess.run(
        [mmseqs_binary],
        check=True,
        capture_output=True,
        text=True,
    )
    mmseqs_version = version_result.stdout.strip()
    if not mmseqs_version:
        raise RuntimeError("MMseqs2 version output is empty")
    with tempfile.TemporaryDirectory(prefix="plantpersulf_global_mmseqs_") as temp:
        root = Path(temp)
        namespaced_fasta = root / "global_multispecies.fasta"
        proteins = prepare_namespaced_fasta(proteomes, namespaced_fasta)
        _write_accession_mapping(mapping_path, proteins)
        for threshold, output in zip(thresholds, outputs, strict=True):
            prefix = root / f"global_{int(threshold * 100):02d}"
            command = easy_cluster_command(
                mmseqs_binary,
                namespaced_fasta,
                prefix,
                root / f"tmp_{int(threshold * 100):02d}",
                identity_threshold=threshold,
                threads=int(mmseqs["threads"]),
            )
            commands.append(command)
            subprocess.run(command, check=True)
            pairs = Path(str(prefix) + "_cluster.tsv")
            rows = global_cluster_rows_from_mmseqs_pairs(proteins, pairs, threshold)
            write_global_cluster_table(output, rows)
            written.append(output)

    manifest = {
        "config_path": str(config_path),
        "config_sha256": _sha256(config_path),
        "mmseqs_binary": mmseqs_binary,
        "mmseqs_version": mmseqs_version,
        "commands": commands,
        "thresholds": thresholds,
        "inputs": [
            {"species": item.species, "path": str(item.path), "sha256": item.sha256}
            for item in proteomes
        ],
        "accession_mapping": {
            "path": str(mapping_path),
            "sha256": _sha256(mapping_path),
        },
        "outputs": {str(path): _sha256(path) for path in written},
        "biological_values_modified": False,
    }
    manifest_path = _provenance_path(primary_output, "manifest.json")
    if manifest_path.exists():
        raise FileExistsError(
            f"refusing to overwrite global MMseqs2 manifest: {manifest_path}"
        )
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return written


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path, default=Path("configs/experiments/multispecies_v2.yaml")
    )
    parser.add_argument("--mmseqs", default="mmseqs")
    args = parser.parse_args()
    for path in build(args.config, args.mmseqs):
        print(path)


if __name__ == "__main__":
    main()
