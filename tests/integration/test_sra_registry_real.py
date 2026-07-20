from pathlib import Path

import pytest


@pytest.mark.network
def test_srp299054_run_metadata_matches_geo_samples(tmp_path: Path) -> None:
    from plantpersulf.download.sra import SraClient

    study = SraClient(cache_dir=tmp_path).get_runs("SRP299054")

    assert study.query_accession == "SRP299054"
    assert study.bioproject_accession == "PRJNA687418"
    assert study.sra_study_accession == "SRP299054"
    assert len(study.runs) == 8
    assert {run.run_accession for run in study.runs} == {
        "SRR13292593",
        "SRR13292594",
        "SRR13292595",
        "SRR13292596",
        "SRR13292597",
        "SRR13292598",
        "SRR13292599",
        "SRR13292600",
    }
    assert {run.experiment_accession for run in study.runs} == {
        "SRX9721582",
        "SRX9721583",
        "SRX9721584",
        "SRX9721585",
        "SRX9721586",
        "SRX9721587",
        "SRX9721588",
        "SRX9721589",
    }
    assert study.cache_path.is_file()
    assert len(study.cache_sha256) == 64
