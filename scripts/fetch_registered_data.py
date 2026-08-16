"""Fetch configured official metadata and rebuild provenance registries."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

from plantpersulf.download.registered import download_registered_files
from plantpersulf.provenance.registry import fetch_registered_metadata


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--metadata-only",
        action="store_true",
        help="fetch official metadata without downloading scientific source files",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("configs/data_sources.yaml"),
    )
    parser.add_argument(
        "--registry-dir",
        type=Path,
        default=Path("data/registry"),
    )
    parser.add_argument("--accession")
    parser.add_argument("--file-class")
    parser.add_argument(
        "--selection-config",
        type=Path,
        default=Path("configs/download_selection.yaml"),
    )
    parser.add_argument(
        "--raw-dir",
        type=Path,
        default=Path("data/raw"),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    arguments = parser.parse_args(argv)
    if arguments.metadata_only:
        if arguments.accession or arguments.file_class:
            parser.error("--metadata-only cannot be combined with download mode")
        metadata_summary = fetch_registered_metadata(
            config_path=arguments.config,
            registry_dir=arguments.registry_dir,
        )
        print(json.dumps(asdict(metadata_summary), sort_keys=True))
    else:
        if not arguments.accession or not arguments.file_class:
            parser.error("download mode requires --accession and --file-class")
        file_classes = tuple(
            value.strip() for value in arguments.file_class.split(",") if value.strip()
        )
        download_summary = download_registered_files(
            accession=arguments.accession,
            file_classes=file_classes,
            selection_path=arguments.selection_config,
            files_registry_path=arguments.registry_dir / "files.tsv",
            datasets_registry_path=arguments.registry_dir / "datasets.tsv",
            downloads_registry_path=arguments.registry_dir / "downloads.tsv",
            raw_dir=arguments.raw_dir,
        )
        print(json.dumps(asdict(download_summary), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
