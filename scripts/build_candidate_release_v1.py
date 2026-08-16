"""Build an immutable blind release only from explicitly frozen inputs."""

from __future__ import annotations

import argparse
from pathlib import Path

from plantpersulf.reporting.candidate_release import (
    ReleaseInputs,
    build_candidate_release,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tomato-model-manifest", type=Path, required=True)
    parser.add_argument("--cross-crop-model-manifest", type=Path)
    parser.add_argument("--candidate-table", type=Path, required=True)
    parser.add_argument("--analysis-plan", type=Path, required=True)
    parser.add_argument("--total-noncontrol-k", type=int, required=True)
    parser.add_argument("--cross-crop-k", type=int, required=True)
    parser.add_argument("--matched-unlabeled-k", type=int, required=True)
    parser.add_argument("--process-control-k", type=int, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    build_candidate_release(
        ReleaseInputs(
            tomato_model_manifest=args.tomato_model_manifest,
            cross_crop_model_manifest=args.cross_crop_model_manifest,
            candidate_table=args.candidate_table,
            analysis_plan=args.analysis_plan,
            total_noncontrol_k=args.total_noncontrol_k,
            cross_crop_k=args.cross_crop_k,
            matched_unlabeled_k=args.matched_unlabeled_k,
            process_control_k=args.process_control_k,
        ),
        args.output_dir,
    )


if __name__ == "__main__":
    main()
