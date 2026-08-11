from pathlib import Path

import pytest

from plantpersulf.reporting.candidate_release import (
    ReleaseInputs,
    build_candidate_release,
)


def test_release_requires_frozen_analysis_plan(tmp_path: Path) -> None:
    inputs = ReleaseInputs(
        tmp_path / "model.json", None, tmp_path / "scores.tsv", None, 20, 0, 20, 2
    )
    with pytest.raises(RuntimeError, match="analysis plan"):
        build_candidate_release(inputs, tmp_path / "release")
