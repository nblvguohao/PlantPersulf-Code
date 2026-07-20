"""PlantPersulf-Code audit command line interface."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

from plantpersulf.download.registered import audit_downloaded_files
from plantpersulf.evidence.preflight import (
    audit_metadata_preflight,
    build_metadata_preflight,
)
from plantpersulf.proteomics.peptide_parser import parse_proteomics_accession
from plantpersulf.proteomics.site_normalizer import audit_site_output
from plantpersulf.provenance.registry import audit_registry


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    audit = subparsers.add_parser(
        "audit-registry",
        help="verify configured datasets and cached metadata provenance",
    )
    audit.add_argument(
        "--config",
        type=Path,
        default=Path("configs/data_sources.yaml"),
    )
    audit.add_argument(
        "--registry-dir",
        type=Path,
        default=Path("data/registry"),
    )
    audit_files = subparsers.add_parser(
        "audit-files",
        help="verify downloaded files against official and local checksums",
    )
    audit_files.add_argument("--accession", required=True)
    audit_files.add_argument(
        "--registry-dir",
        type=Path,
        default=Path("data/registry"),
    )
    audit_files.add_argument(
        "--selection-config",
        type=Path,
        default=Path("configs/download_selection.yaml"),
    )
    parse_proteomics = subparsers.add_parser(
        "parse-proteomics",
        help="parse registered proteomics results without assigning labels",
    )
    parse_proteomics.add_argument("--accession", required=True)
    parse_proteomics.add_argument(
        "--output-root",
        type=Path,
        default=Path("data/interim"),
    )
    parse_proteomics.add_argument(
        "--registry-dir",
        type=Path,
        default=Path("data/registry"),
    )
    parse_proteomics.add_argument(
        "--selection-config",
        type=Path,
        default=Path("configs/download_selection.yaml"),
    )
    audit_sites = subparsers.add_parser(
        "audit-sites",
        help="audit sequence-verified proteomics site coordinates",
    )
    audit_sites.add_argument("--accession", required=True)
    audit_sites.add_argument(
        "--output-root",
        type=Path,
        default=Path("data/interim"),
    )
    audit_sites.add_argument(
        "--registry-dir",
        type=Path,
        default=Path("data/registry"),
    )
    build_preflight = subparsers.add_parser(
        "build-evidence-preflight",
        help="build metadata-only evidence inventory without assigning labels",
    )
    build_preflight.add_argument("--version", choices=("v1",), required=True)
    build_preflight.add_argument(
        "--output-root",
        type=Path,
        default=Path("data/interim"),
    )
    build_preflight.add_argument(
        "--registry-dir",
        type=Path,
        default=Path("data/registry"),
    )
    build_preflight.add_argument(
        "--policy",
        type=Path,
        default=Path("configs/evidence_preflight_v1.yaml"),
    )
    build_preflight.add_argument(
        "--selection-config",
        type=Path,
        default=Path("configs/download_selection.yaml"),
    )
    audit_preflight = subparsers.add_parser(
        "audit-evidence-preflight",
        help="audit the metadata-only evidence inventory",
    )
    audit_preflight.add_argument("--version", choices=("v1",), required=True)
    audit_preflight.add_argument(
        "--output-root",
        type=Path,
        default=Path("data/interim"),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = build_parser().parse_args(argv)
    if arguments.command == "audit-registry":
        registry_summary = audit_registry(
            registry_dir=arguments.registry_dir,
            config_path=arguments.config,
        )
        print(json.dumps(asdict(registry_summary), sort_keys=True))
        return 0
    if arguments.command == "audit-files":
        file_summary = audit_downloaded_files(
            accession=arguments.accession,
            selection_path=arguments.selection_config,
            files_registry_path=arguments.registry_dir / "files.tsv",
            datasets_registry_path=arguments.registry_dir / "datasets.tsv",
            downloads_registry_path=arguments.registry_dir / "downloads.tsv",
        )
        print(json.dumps(asdict(file_summary), sort_keys=True))
        return 0
    if arguments.command == "parse-proteomics":
        parse_summary = parse_proteomics_accession(
            accession=arguments.accession,
            output_root=arguments.output_root,
            registry_dir=arguments.registry_dir,
            selection_path=arguments.selection_config,
        )
        print(json.dumps(asdict(parse_summary), sort_keys=True))
        return 0
    if arguments.command == "audit-sites":
        site_summary = audit_site_output(
            accession=arguments.accession,
            output_root=arguments.output_root,
            registry_dir=arguments.registry_dir,
        )
        print(json.dumps(asdict(site_summary), sort_keys=True))
        return 0
    if arguments.command == "build-evidence-preflight":
        preflight_summary = build_metadata_preflight(
            policy_path=arguments.policy,
            selection_path=arguments.selection_config,
            registry_dir=arguments.registry_dir,
            output_directory=(
                arguments.output_root / "evidence_preflight_v1"
            ),
        )
        print(json.dumps(asdict(preflight_summary), sort_keys=True))
        return 0
    if arguments.command == "audit-evidence-preflight":
        audited_preflight = audit_metadata_preflight(
            arguments.output_root / "evidence_preflight_v1"
        )
        print(json.dumps(asdict(audited_preflight), sort_keys=True))
        return 0
    raise RuntimeError(f"unsupported command: {arguments.command}")


if __name__ == "__main__":
    raise SystemExit(main())
