from pathlib import Path

import pytest

from plantpersulf.reporting.candidate_release import ReleaseInputs


def test_cross_crop_quota_cannot_exceed_one_quarter() -> None:
    with pytest.raises(RuntimeError, match="25%"):
        ReleaseInputs(
            Path("model"), Path("cross"), Path("scores"), Path("plan"), 20, 6, 20, 2
        ).validate()
