from pathlib import Path

import pytest

from plantpersulf.reporting.candidate_release import (
    ReleaseInputs,
    build_candidate_release,
    release_id_from_inputs,
)


def test_existing_release_directory_is_never_overwritten(tmp_path: Path) -> None:
    output = tmp_path / "release"
    output.mkdir()
    with pytest.raises(FileExistsError):
        build_candidate_release(
            ReleaseInputs(
                tmp_path / "model",
                None,
                tmp_path / "table",
                tmp_path / "plan",
                1,
                0,
                0,
                0,
            ),
            output,
        )


def test_changed_candidate_order_gets_a_new_release_id(tmp_path: Path) -> None:
    model, plan, table = (
        tmp_path / "model.json",
        tmp_path / "plan.json",
        tmp_path / "candidates.tsv",
    )
    model.write_text(
        '{"complete": true, "selected_model": "pu_logistic", '
        '"model_release_created": true}\n',
        encoding="utf-8",
    )
    plan.write_text(
        '{"status":"frozen","power_passed":true,"collaborator_protocol_frozen":true,"total_noncontrol_k":1,"cross_crop_k":0,"matched_unlabeled_k":0,"process_control_k":0}\n',
        encoding="utf-8",
    )
    header = "site_key\tprotein_accession\tarm\tselection_order\tmatched_stratum\n"
    table.write_text(header + "site-a\tprotein-a\tT\t1\tstratum-a\n", encoding="utf-8")
    first = ReleaseInputs(model, None, table, plan, 1, 0, 0, 0)
    first_id = release_id_from_inputs(first)
    table.write_text(header + "site-a\tprotein-a\tT\t2\tstratum-a\n", encoding="utf-8")
    assert (
        release_id_from_inputs(ReleaseInputs(model, None, table, plan, 1, 0, 0, 0))
        != first_id
    )
