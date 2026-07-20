"""Fetch configured official metadata and rebuild provenance registries."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

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
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    arguments = parser.parse_args(argv)
    if not arguments.metadata_only:
        parser.error("Task 1 only permits --metadata-only")
    summary = fetch_registered_metadata(
        config_path=arguments.config,
        registry_dir=arguments.registry_dir,
    )
    print(json.dumps(asdict(summary), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
