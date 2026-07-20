"""PlantPersulf-Code audit command line interface."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

from plantpersulf.download.registered import audit_downloaded_files
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
    raise RuntimeError(f"unsupported command: {arguments.command}")


if __name__ == "__main__":
    raise SystemExit(main())
