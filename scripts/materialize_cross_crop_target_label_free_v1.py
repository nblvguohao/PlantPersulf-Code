#!/usr/bin/env python
"""Materialize the frozen Task 9A source batches and tomato candidate request."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from plantpersulf.workflows.cross_crop_materialization import (
    materialize_cross_crop_request,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("configs/experiments/cross_crop_materialization_v1.yaml"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/processed/cross_crop_target_label_free_v1"),
    )
    parser.add_argument("--request-registry", type=Path)
    args = parser.parse_args()
    manifest = materialize_cross_crop_request(
        args.config, args.output_dir, args.request_registry
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
