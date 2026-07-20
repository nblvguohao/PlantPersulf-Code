import hashlib
from pathlib import Path

import pytest


@pytest.mark.network
def test_pxd006140_metadata_matches_expected_project(tmp_path: Path) -> None:
    from plantpersulf.download.pride import PrideClient

    project = PrideClient(cache_dir=tmp_path).get_project("PXD006140")

    assert project.accession == "PXD006140"
    assert "persulfidation" in project.title.lower()
    assert project.publication_date == "2018-10-24"
    assert any("Arabidopsis thaliana" in organism for organism in project.organisms)
    assert project.files

    cache_path = tmp_path / "PXD006140.json"
    assert cache_path.is_file()
    assert hashlib.sha256(cache_path.read_bytes()).hexdigest() == project.cache_sha256
