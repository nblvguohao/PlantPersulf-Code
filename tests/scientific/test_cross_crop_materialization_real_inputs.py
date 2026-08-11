"""Task 9A materialization must stay traceable and tomato-label-free."""

import json
from pathlib import Path

from plantpersulf.provenance.audit import assert_registered_input
from plantpersulf.workflows.cross_crop_materialization import (
    materialize_cross_crop_request,
)


def test_real_four_domain_request_is_pu_only_and_has_no_tomato_labels(
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "cross_crop_request"
    request_registry = tmp_path / "request_registry.tsv"
    manifest = materialize_cross_crop_request(
        Path("configs/experiments/cross_crop_materialization_v1.yaml"),
        output_dir,
        request_registry,
    )

    request_path = output_dir / "request.json"
    assert_registered_input(request_path, request_registry)
    request = json.loads(request_path.read_text(encoding="utf-8"))
    assert manifest["source_domain_count"] == 4
    assert manifest["source_species"] == [
        "Arabidopsis thaliana",
        "Magnaporthe oryzae",
        "Oryza sativa",
    ]
    assert "labels" not in request["target_candidates"]
    assert len(request["target_candidates"]["site_keys"]) == 251787
    assert all(
        set(source["labels"]) == {"positive", "unlabeled"}
        for source in request["sources"]
    )
    pxd006140 = next(
        source
        for source in request["sources"]
        if source["study_accession"] == "PXD006140"
    )
    assert "P42737-2" in pxd006140["protein_ids"]
    assert pxd006140["labels"][pxd006140["protein_ids"].index("P42737-2")] == "positive"
