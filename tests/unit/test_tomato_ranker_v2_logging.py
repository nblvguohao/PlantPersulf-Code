"""The real-data workflow exposes durable, non-biological progress logs."""

from pathlib import Path

from plantpersulf.workflows.tomato_ranker_v2 import configure_run_logger


def test_run_logger_writes_a_phase_marker(tmp_path: Path) -> None:
    log_path = tmp_path / "run.log"
    logger = configure_run_logger(log_path)
    logger.info("phase=unit_test")
    for handler in logger.handlers:
        handler.flush()
        handler.close()
    assert "phase=unit_test" in log_path.read_text(encoding="utf-8")
