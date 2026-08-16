import inspect

from plantpersulf.workflows.cross_crop_target_label_free import (
    run_cross_crop_target_label_free,
)


def test_cross_crop_runner_has_no_target_label_parameter() -> None:
    names = set(inspect.signature(run_cross_crop_target_label_free).parameters)
    assert not any("target_label" in name for name in names)
