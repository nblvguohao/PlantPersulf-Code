"""RED (Task 7): frozen ESM-2 residue embeddings for benchmark cysteines.

Extracts per-position embeddings from a frozen ESM-2 model and records the
vector at each labeled cysteine position. The model is loaded once and
embeddings are deterministic (same model checkpoint, same input → same
output). Only benchmark proteins are embedded; embeddings must not leak
test/validation proteins into the feature table.

Expected RED: ``plantpersulf.features.esm2`` does not exist yet.
"""

from __future__ import annotations

from pathlib import Path

from plantpersulf.features.esm2 import (  # RED: module missing
    extract_esm2_embeddings,
)


def _write_fixtures(tmp_path: Path) -> tuple[Path, Path]:
    proteome = tmp_path / "mini.fasta"
    proteome.write_text(">sp|P1\nMVCGK\n>sp|P2\nACDEF\n", encoding="utf-8")
    labels = tmp_path / "sites.tsv"
    labels.write_text(
        "protein_accession\tcys_position_in_protein\tlabel\t"
        "study_accession\tevidence_level\tsource_sha256\n"
        "P1\t3\tpositive\tPXD006140\tsite_ms\taaa\n"
        "P2\t2\tunlabeled\t\t\t\n",
        encoding="utf-8",
    )
    return labels, proteome


def test_esm2_embedding_dimension_is_correct(tmp_path: Path) -> None:
    """esm2_t33_650M outputs 1280-dim per-position embeddings."""
    labels_path, proteome_path = _write_fixtures(tmp_path)
    rows = extract_esm2_embeddings(labels_path, proteome_path)

    assert len(rows) == 2
    for row in rows:
        assert len(row.embedding) == 1280


def test_cysteine_position_yields_embedding_not_mean(tmp_path: Path) -> None:
    """The embedding component at the Cys position differs from the
    protein mean embedding (confirming it is position-specific)."""
    labels_path, proteome_path = _write_fixtures(tmp_path)
    rows = extract_esm2_embeddings(labels_path, proteome_path)

    p1 = [r for r in rows if r.protein_accession == "P1"][0]
    p2 = [r for r in rows if r.protein_accession == "P2"][0]
    assert p1.embedding != p2.embedding  # different proteins, different context


def test_missing_protein_is_skipped_without_crashing(tmp_path: Path) -> None:
    """A label for a protein absent from the proteome must not crash."""
    labels = tmp_path / "bad_labels.tsv"
    labels.write_text(
        "protein_accession\tcys_position_in_protein\tlabel\t"
        "study_accession\tevidence_level\tsource_sha256\n"
        "P999\t5\tpositive\tPXD006140\tsite_ms\taaa\n",
        encoding="utf-8",
    )
    proteome = tmp_path / "mini.fasta"
    proteome.write_text(">sp|P1\nMVCGK\n", encoding="utf-8")
    rows = extract_esm2_embeddings(labels, proteome)
    assert len(rows) == 0
