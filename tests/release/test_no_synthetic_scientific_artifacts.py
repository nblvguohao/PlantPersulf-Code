from pathlib import Path

import pytest


def test_release_rejects_synthetic_scientific_artifact(tmp_path: Path) -> None:
    from plantpersulf.provenance.audit import (
        assert_no_forbidden_scientific_artifacts,
    )

    release_path = tmp_path / "release"
    release_path.mkdir()
    (release_path / "synthetic_results.tsv").write_text(
        "policy-test-only\n",
        encoding="utf-8",
    )

    with pytest.raises(RuntimeError, match="forbidden scientific artifact"):
        assert_no_forbidden_scientific_artifacts(
            release_path,
            config_path=Path("configs/scientific_integrity.yaml"),
        )


def test_release_without_forbidden_artifact_names_is_accepted(tmp_path: Path) -> None:
    from plantpersulf.provenance.audit import (
        assert_no_forbidden_scientific_artifacts,
    )

    release_path = tmp_path / "release"
    release_path.mkdir()
    (release_path / "registered_results.tsv").write_text(
        "policy-test-only\n",
        encoding="utf-8",
    )

    assert_no_forbidden_scientific_artifacts(
        release_path,
        config_path=Path("configs/scientific_integrity.yaml"),
    )
