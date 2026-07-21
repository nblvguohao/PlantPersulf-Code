import csv
import json
from pathlib import Path


def test_stop_release_is_label_free_and_fail_closed(
    tmp_path: Path,
    capsys,  # type: ignore[no-untyped-def]
) -> None:
    from plantpersulf.cli import main

    exit_code = main(
        [
            "build-benchmark-readiness",
            "--version",
            "v1",
            "--output-root",
            str(tmp_path / "interim"),
        ]
    )
    result = json.loads(capsys.readouterr().out)
    readiness = tmp_path / "interim/benchmark_readiness_v1"

    assert exit_code == 0
    assert result["decision"] == "STOP"
    assert result["eligible_site_count"] == 0
    assert result["eligible_study_count"] == 0
    manifest = json.loads(
        (readiness / "manifest.json").read_text(encoding="utf-8")
    )
    assert manifest["benchmark_created"] is False
    assert manifest["labels_created"] is False
    assert manifest["nondetection_labeled_negative"] is False
    for path in readiness.glob("*.tsv"):
        with path.open(encoding="utf-8", newline="") as handle:
            fields = tuple(csv.DictReader(handle, delimiter="\t").fieldnames or ())
        assert not any("label" in field.lower() for field in fields)
    assert not (tmp_path / "processed/benchmark_v1").exists()
