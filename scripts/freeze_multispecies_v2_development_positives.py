#!/usr/bin/env python
"""Administrative one-time freezing of the v2 development-positive manifest."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import yaml

from plantpersulf.benchmark.multispecies_splits import (
    DEVELOPMENT_SPLIT,
    load_frozen_multispecies_split,
)
from plantpersulf.proteomics.multispecies_v2_sources import (
    development_positive_rows,
    write_frozen_development_positives,
)
from plantpersulf.provenance.audit import assert_registered_input


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser


def _registered(path: Path) -> None:
    for registry in (
        Path("data/registry/cross_crop_target_label_free_inputs_v1.tsv"),
        Path("data/registry/supplementary_sources.tsv"),
        Path("data/registry/model_inputs.tsv"),
    ):
        try:
            assert_registered_input(path, registry)
            return
        except RuntimeError:
            continue
    raise RuntimeError(f"input is absent from every approved registry: {path}")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(argv: list[str] | None = None) -> None:
    arguments = build_argument_parser().parse_args(argv)
    cfg = yaml.safe_load(arguments.config.read_text(encoding="utf-8"))
    evidence = cfg["positive_evidence"]
    references = {
        item["species"]: Path(item["path"]) for item in cfg["reference_proteomes"]
    }
    evidence_paths = {
        key: Path(value) for key, value in evidence.items() if key != "registry"
    }
    split_path = Path(cfg["strict_cluster_holdout"]["split_path"])
    for path in (*references.values(), *evidence_paths.values(), split_path):
        _registered(path)

    from plantpersulf.features.sequence import _load_proteome
    from plantpersulf.proteomics.pxd063170_sites import load_ensembl_fungi_proteome

    proteomes = {
        "arabidopsis": _load_proteome(references["arabidopsis"]),
        "tomato": _load_proteome(references["tomato"]),
        "rice": _load_proteome(references["rice"]),
        "magnaporthe": load_ensembl_fungi_proteome(references["magnaporthe"]),
    }
    frozen = load_frozen_multispecies_split(split_path)
    allowed = {
        species: {
            row.global_protein_id.split("|", 1)[1]
            for row in frozen.rows
            if row.species == species and row.split == DEVELOPMENT_SPLIT
        }
        for species in proteomes
    }
    positives = development_positive_rows(
        benchmark=evidence_paths["arabidopsis_benchmark"],
        tomato_xlsx=evidence_paths["tomato_kiae271"],
        rice_sd01=evidence_paths["rice_sd01"],
        rice_sd04=evidence_paths["rice_sd04"],
        rice_ss_all=evidence_paths["rice_ss_all"],
        magnaporthe_tsv=evidence_paths["magnaporthe_sites"],
        proteomes=proteomes,
        allowed_proteins=allowed,
    )
    benchmark_hash = _sha256(evidence_paths["arabidopsis_benchmark"])
    rice_hashes = ";".join(
        _sha256(evidence_paths[key])
        for key in ("rice_sd01", "rice_sd04", "rice_ss_all")
    )
    source_hashes = {
        "PXD006140": benchmark_hash,
        "PXD024061": benchmark_hash,
        "KIAE271_SUPPL": _sha256(evidence_paths["tomato_kiae271"]),
        "PXD072089": rice_hashes,
        "PXD063170": _sha256(evidence_paths["magnaporthe_sites"]),
    }
    write_frozen_development_positives(
        arguments.output,
        positives,
        source_sha256_by_study=source_hashes,
    )
    print(f"froze {len(positives)} development evidence rows: {arguments.output}")


if __name__ == "__main__":
    main()
