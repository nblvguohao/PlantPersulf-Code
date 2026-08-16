"""Run the separately approved frozen tomato-v2 workflow."""

from __future__ import annotations

import argparse
from pathlib import Path

from plantpersulf.workflows.tomato_ranker_v2 import run_tomato_ranker_v2


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    run_tomato_ranker_v2(args.config, args.output_dir)


if __name__ == "__main__":
    main()
