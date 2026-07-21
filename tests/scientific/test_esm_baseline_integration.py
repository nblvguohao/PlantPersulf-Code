"""RED (Task 8): ESM embedding + linear head baseline, end-to-end.

Integration test using the real frozen ESM-2 model on tiny synthetic
sequences (same pattern as tests/scientific/test_esm2_features.py) to
confirm esm_linear_head_scores works on genuine 1280-dim embeddings
produced by extract_esm2_embeddings, not just on toy numeric fixtures.

Expected RED: ``plantpersulf.models.esm_baseline`` does not exist yet.
"""

from __future__ import annotations

from pathlib import Path

from plantpersulf.features.esm2 import extract_esm2_embeddings
from plantpersulf.models.esm_baseline import (  # RED: module missing
    esm_linear_head_scores,
)


def _write_fixtures(tmp_path: Path) -> tuple[Path, Path]:
    proteome = tmp_path / "mini.fasta"
    proteome.write_text(
        ">sp|P1\nMVCGK\n>sp|P2\nACDEF\n>sp|P3\nGCCCC\n", encoding="utf-8"
    )
    labels = tmp_path / "sites.tsv"
    labels.write_text(
        "protein_accession\tcys_position_in_protein\tlabel\t"
        "study_accession\tevidence_level\tsource_sha256\n"
        "P1\t3\tpositive\tPXD006140\tsite_ms\taaa\n"
        "P2\t2\tunlabeled\t\t\t\n"
        "P3\t2\tunlabeled\t\t\t\n"
        "P3\t3\tpositive\tPXD006140\tsite_ms\tbbb\n",
        encoding="utf-8",
    )
    return labels, proteome


def test_linear_head_scores_real_esm2_embeddings(tmp_path: Path) -> None:
    labels_path, proteome_path = _write_fixtures(tmp_path)
    rows = extract_esm2_embeddings(labels_path, proteome_path)
    assert len(rows) == 4

    train_rows = [r for r in rows if r.protein_accession != "P3"]
    predict_rows = [r for r in rows if r.protein_accession == "P3"]

    scores = esm_linear_head_scores(
        train_embeddings=[r.embedding for r in train_rows],
        train_y=[r.label for r in train_rows],
        predict_embeddings=[r.embedding for r in predict_rows],
        seed=0,
    )

    assert len(scores) == len(predict_rows)
    for score in scores:
        assert 0.0 <= score <= 1.0
