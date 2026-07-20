import hashlib
from pathlib import Path

import pytest


@pytest.mark.network
def test_gse163745_metadata_is_tomato_ripening_study(tmp_path: Path) -> None:
    from plantpersulf.download.geo import GeoClient

    series = GeoClient(cache_dir=tmp_path).get_series("GSE163745")

    assert series.accession == "GSE163745"
    assert series.organism == "Solanum lycopersicum"
    assert series.publication_date == "2020-12-24"
    assert len(series.samples) == 8
    assert all(sample.accession.startswith("GSM") for sample in series.samples)

    cache_path = tmp_path / "GSE163745.json"
    assert cache_path.is_file()
    assert hashlib.sha256(cache_path.read_bytes()).hexdigest() == series.cache_sha256
