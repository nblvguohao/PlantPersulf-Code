import hashlib
from pathlib import Path

import pytest


@pytest.mark.network
def test_pxd051570_uses_its_official_iprox_host(tmp_path: Path) -> None:
    from plantpersulf.download.iprox import IproxClient

    dataset = IproxClient(cache_dir=tmp_path).get_dataset("PXD051570")

    assert dataset.accession == "PXD051570"
    assert "ripening" in dataset.title.lower()
    assert dataset.publication_date == "2024-04-19"
    assert dataset.organisms == ("Solanum lycopersicum",)
    assert dataset.files

    cache_path = tmp_path / "PXD051570.json"
    assert cache_path.is_file()
    assert hashlib.sha256(cache_path.read_bytes()).hexdigest() == dataset.cache_sha256
