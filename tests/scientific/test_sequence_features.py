"""RED (Task 7): sequence-based cysteine features from the reference proteome.

Extracts per-cysteine features — flanking windows, amino acid composition,
physicochemical properties (hydrophobicity, charge, volume), protein length,
and cysteine density — driven by the benchmark labels and the SHA256-pinned
reference proteome. Features are deterministic and do not depend on any
external model.

Expected RED: ``plantpersulf.features.sequence`` does not exist yet.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from plantpersulf.features.sequence import (  # RED: module missing
    extract_sequence_features,
)


def _write_fixtures(tmp_path: Path) -> tuple[Path, Path]:
    """Write minimal proteome + benchmark labels with all labels needed for tests."""
    proteome = tmp_path / "mini.fasta"
    proteome.write_text(
        ">sp|P1\nMVCGK\n"
        ">sp|P2\nACDEF\n"
        ">sp|P3\nGCCCC\n",
        encoding="utf-8",
    )
    labels = tmp_path / "sites.tsv"
    labels.write_text(
        "protein_accession\tcys_position_in_protein\tlabel\t"
        "study_accession\tevidence_level\tsource_sha256\n"
        "P1\t3\tpositive\tPXD006140\tsite_ms\taaa\n"
        "P2\t2\tunlabeled\t\t\t\n"
        "P3\t2\tpositive\tPXD024061\tsite_ms\tbbb\n"
        "P3\t3\tunlabeled\t\t\t\n"
        "P3\t4\tunlabeled\t\t\t\n"
        "P3\t5\tunlabeled\t\t\t\n",
        encoding="utf-8",
    )
    return labels, proteome


def test_flanking_window_is_centered_on_cysteine(tmp_path: Path) -> None:
    labels_path, proteome_path = _write_fixtures(tmp_path)
    rows = extract_sequence_features(labels_path, proteome_path, window_radius=2)

    by_prot = {r.protein_accession: r for r in rows}
    # P1: MVCGK, C at pos 3, window ±2 = MVCGK (full protein, no padding needed)
    assert by_prot["P1"].flanking_window == "MVCGK"
    # P3: GCCCC, C at pos 4, window ±2 = CCCCC
    p3_rows = [r for r in rows if r.protein_accession == "P3"]
    p3_by_pos = {r.cys_position: r.flanking_window for r in p3_rows}
    assert len(p3_by_pos[2]) == 5
    assert p3_by_pos[2][2] == "C"  # centre of the window is the cysteine
    assert len(p3_by_pos[4]) == 5
    assert p3_by_pos[5][-1] == "X"  # C-term padded


def test_every_labeled_cysteine_has_a_feature_row(tmp_path: Path) -> None:
    labels_path, proteome_path = _write_fixtures(tmp_path)
    rows = extract_sequence_features(labels_path, proteome_path)

    # 6 labeled cysteines → 6 feature rows
    assert len(rows) == 6
    keys = {(r.protein_accession, r.cys_position) for r in rows}
    assert ("P1", 3) in keys
    assert ("P2", 2) in keys


def test_feature_values_are_deterministic(tmp_path: Path) -> None:
    labels_path, proteome_path = _write_fixtures(tmp_path)
    first = extract_sequence_features(labels_path, proteome_path)
    second = extract_sequence_features(labels_path, proteome_path)

    assert len(first) == len(second)
    for a, b in zip(first, second, strict=True):
        assert a.hydrophobicity == b.hydrophobicity
        assert a.cys_density == b.cys_density
        assert a.protein_length == b.protein_length


def test_cysteine_density_is_correct(tmp_path: Path) -> None:
    labels_path, proteome_path = _write_fixtures(tmp_path)
    rows = extract_sequence_features(labels_path, proteome_path)

    by_prot = {}
    for r in rows:
        by_prot.setdefault(r.protein_accession, []).append(r)
    # P3: GCCCC → 4 Cys / 5 residues = 0.8
    for r in by_prot["P3"]:
        assert r.cys_density == pytest.approx(4 / 5)


def test_hydrophobicity_is_frame_centered_on_cys(tmp_path: Path) -> None:
    """Hydrophobicity is computed from the flanking window center ± window_radius."""
    labels_path, proteome_path = _write_fixtures(tmp_path)
    rows = extract_sequence_features(labels_path, proteome_path, window_radius=1)

    # P2: ACDEF, C at pos 2, window ±1 = ACD (radius=1, center C)
    p2 = [r for r in rows if r.protein_accession == "P2"][0]
    assert len(p2.flanking_window) == 3
    assert p2.flanking_window[1] == "C"


def test_proteome_loader_accepts_non_uniprot_headers(tmp_path: Path) -> None:
    """EnsemblFungi-style headers (first whitespace token = accession, no
    pipe fields) must load too — the cross-species track scores Magnaporthe
    sites against the MG8 proteome."""
    proteome = tmp_path / "ensembl.fasta"
    proteome.write_text(
        ">MGG_00001T0 pep chromosome:MG8:1:100:200:1 gene:MGG_00001\n"
        "MACDEFGHIK\n"
        ">MGG_00002T0 pep chromosome:MG8:1:300:400:-1\n"
        "MADEFGCCLM\n",
        encoding="utf-8",
    )
    labels = tmp_path / "sites.tsv"
    labels.write_text(
        "protein_accession\tcys_position_in_protein\tlabel\t"
        "study_accession\tevidence_level\tsource_sha256\n"
        "MGG_00001T0\t3\tpositive\tPXD063170\tsite_ms\tccc\n",
        encoding="utf-8",
    )
    rows = extract_sequence_features(labels, proteome)
    assert [(r.protein_accession, r.cys_position) for r in rows] == [
        ("MGG_00001T0", 3)
    ]
