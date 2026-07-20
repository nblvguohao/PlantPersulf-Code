"""PlantPersulf-Code audit command line interface."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

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
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = build_parser().parse_args(argv)
    if arguments.command == "audit-registry":
        summary = audit_registry(
            registry_dir=arguments.registry_dir,
            config_path=arguments.config,
        )
        print(json.dumps(asdict(summary), sort_keys=True))
        return 0
    raise RuntimeError(f"unsupported command: {arguments.command}")


if __name__ == "__main__":
    raise SystemExit(main())
