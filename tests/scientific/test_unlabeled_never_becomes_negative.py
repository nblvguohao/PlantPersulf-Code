"""Scientific policy guard for the PU-loss implementation."""

from pathlib import Path


def test_new_model_source_has_no_hard_negative_label_literal() -> None:
    source = Path("src/plantpersulf/models/ranking_loss.py").read_text(
        encoding="utf-8"
    )
    assert '"negative"' not in source
    assert "'negative'" not in source
