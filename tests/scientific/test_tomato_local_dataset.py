"""RED: within-tomato PU dataset assembly and leakage-safe grouped folds.

The kiae271 panel (99 coordinate-verified Solanum lycopersicum Cys
persulfidation sites over 88 proteins) is used here as *target-species
supervision*, not as a cross-species holdout. That requires three things the
existing machinery does not provide:

1. A PU background drawn from the tomato proteome, with every ambiguous
   Dataset S1 row (the 20 that failed coordinate/localization verification)
   removed from the unlabeled pool — such a row is neither a positive nor a
   negative and must not be silently recycled as background.
2. Two explicitly separated evaluation arenas: ``proteome`` (all-other-Cys in
   the whole tomato proteome, subsampled) and ``panel`` (only the other Cys of
   the 88 detected proteins). The second controls for MS protein-level
   detection/abundance bias, which the first cannot.
3. Grouped k-fold assignment whose unit is a homology cluster, so no protein
   (and no homolog of it) has sites in both the train and test side of a fold.

Expected RED: ``plantpersulf.proteomics.tomato_local_dataset`` does not exist.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from plantpersulf.proteomics.tomato_local_dataset import (  # RED: module missing
    ARENA_PANEL,
    ARENA_PROTEOME,
    assign_grouped_folds,
    build_tomato_pu_rows,
    kiae271_excluded_keys,
)

# Two panel proteins (one with 3 Cys, one with 2) plus two background proteins.
MINI_PROTEOME = {
    "PANEL1": "MACDEFGHIKCLMNPQRSTCVWY",  # Cys at 3, 11, 20
    "PANEL2": "MCADEFGHIKLMNPQRSTVWYC",  # Cys at 2, 22
    "BACK1": "MGCDEFGHIKCLMNPQRSTVWY",  # Cys at 3, 11
    "BACK2": "MGGCDEFGHIKLMNPQRSTVWY",  # Cys at 4
}

HEADER = (
    "Proteins",
    "Positions within proteins",
    "Leading proteins",
    "Protein",
    "Fasta headers",
    "Localization prob",
    "PEP",
    "Score",
    "Score for localization",
    "Number of S(C)",
    "Amino acid",
    "Sequence window",
    "S(C) Probabilities",
    "Position in peptide",
    "Charge",
    "Intensity LCD_1",
    "Intensity LCD_2",
    "Intensity WT_1",
    "Intensity WT_2",
    "Intensity LCD",
    "Intensity WT",
)


def _write_xlsx(path: Path, rows: list[tuple]) -> None:
    import openpyxl

    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Supplementary Dataset S1"
    ws.append(("Title row (ignored)",))
    ws.append(HEADER)
    for row in rows:
        ws.append(row)
    wb.save(path)


def _row(
    proteins: str,
    positions: str,
    leading: str,
    amino: str = "C",
    loc_prob: float = 1.0,
    intensity_lcd: float = 0.0,
    intensity_wt: float = 0.0,
) -> tuple:
    return (
        proteins,
        positions,
        leading,
        leading,
        f"{leading} description",
        loc_prob,
        0.001,
        100,
        100,
        1,
        amino,
        "SEQWINDOW",
        "1",
        "PEP",
        2,
        intensity_lcd,
        0,
        intensity_wt,
        0,
        intensity_lcd,
        intensity_wt,
    )


@pytest.fixture()
def mini_xlsx(tmp_path: Path) -> Path:
    """Two verified positives (PANEL1 C3, PANEL2 C2) and one row that fails
    localization (PANEL1 C11) — the failing row is the ambiguity that must be
    kept out of the unlabeled pool."""
    xlsx = tmp_path / "DSs.xlsx"
    _write_xlsx(
        xlsx,
        [
            _row("tr|PANEL1|P1_SOLLC", "3", "tr|PANEL1|P1_SOLLC", intensity_lcd=1000.0),
            _row("tr|PANEL2|P2_SOLLC", "2", "tr|PANEL2|P2_SOLLC", intensity_wt=800.0),
            _row(
                "tr|PANEL1|P1_SOLLC",
                "11",
                "tr|PANEL1|P1_SOLLC",
                loc_prob=0.4,
                intensity_lcd=500.0,
            ),
        ],
    )
    return xlsx


def test_excluded_keys_cover_every_unverified_dataset_s1_site(
    mini_xlsx: Path,
) -> None:
    excluded = kiae271_excluded_keys(mini_xlsx, MINI_PROTEOME)
    assert ("PANEL1", 11) in excluded
    # Verified positives are NOT "excluded keys" — they are positives.
    assert ("PANEL1", 3) not in excluded
    assert ("PANEL2", 2) not in excluded


def test_ambiguous_row_is_neither_positive_nor_unlabeled(mini_xlsx: Path) -> None:
    rows = build_tomato_pu_rows(
        mini_xlsx,
        MINI_PROTEOME,
        arena=ARENA_PANEL,
        ratio=20,
        seed=12345,
    )
    keys = {(r.protein_accession, r.cys_position, r.label) for r in rows}
    assert ("PANEL1", 11, "positive") not in keys
    assert ("PANEL1", 11, "unlabeled") not in keys


def test_panel_arena_background_is_only_panel_protein_cysteines(
    mini_xlsx: Path,
) -> None:
    rows = build_tomato_pu_rows(
        mini_xlsx,
        MINI_PROTEOME,
        arena=ARENA_PANEL,
        ratio=20,
        seed=12345,
    )
    assert {r.protein_accession for r in rows} == {"PANEL1", "PANEL2"}
    unlabeled = {
        (r.protein_accession, r.cys_position) for r in rows if r.label == "unlabeled"
    }
    # PANEL1 C20 and PANEL2 C22 remain; C3/C2 are positives, C11 is excluded.
    assert unlabeled == {("PANEL1", 20), ("PANEL2", 22)}


def test_proteome_arena_background_is_subsampled_at_the_fixed_ratio(
    mini_xlsx: Path,
) -> None:
    rows = build_tomato_pu_rows(
        mini_xlsx,
        MINI_PROTEOME,
        arena=ARENA_PROTEOME,
        ratio=2,
        seed=12345,
    )
    n_pos = sum(1 for r in rows if r.label == "positive")
    n_unl = sum(1 for r in rows if r.label == "unlabeled")
    assert n_pos == 2
    assert n_unl == 4  # ratio 2 x 2 positives


def test_proteome_arena_subsample_is_deterministic(mini_xlsx: Path) -> None:
    def _keys(seed: int) -> list[tuple[str, int]]:
        return [
            (r.protein_accession, r.cys_position)
            for r in build_tomato_pu_rows(
                mini_xlsx,
                MINI_PROTEOME,
                arena=ARENA_PROTEOME,
                ratio=2,
                seed=seed,
            )
        ]

    assert _keys(12345) == _keys(12345)


def test_positive_is_never_also_in_the_unlabeled_pool(mini_xlsx: Path) -> None:
    rows = build_tomato_pu_rows(
        mini_xlsx,
        MINI_PROTEOME,
        arena=ARENA_PROTEOME,
        ratio=2,
        seed=12345,
    )
    positives = {
        (r.protein_accession, r.cys_position) for r in rows if r.label == "positive"
    }
    unlabeled = {
        (r.protein_accession, r.cys_position) for r in rows if r.label == "unlabeled"
    }
    assert positives & unlabeled == set()


def test_grouped_folds_never_split_a_homology_cluster(mini_xlsx: Path) -> None:
    """The leakage-safety property: PANEL1 and BACK1 are homologs (same
    cluster), so they must land in the same fold no matter the seed."""
    rows = build_tomato_pu_rows(
        mini_xlsx,
        MINI_PROTEOME,
        arena=ARENA_PROTEOME,
        ratio=2,
        seed=12345,
    )
    clusters = {"PANEL1": "cl_A", "BACK1": "cl_A", "PANEL2": "cl_B", "BACK2": "cl_C"}
    folded = assign_grouped_folds(rows, clusters, n_folds=2, seed=20260810)
    by_cluster: dict[str, set[int]] = {}
    for row in folded:
        by_cluster.setdefault(row.cluster_id, set()).add(row.fold)
    assert all(len(folds) == 1 for folds in by_cluster.values())


def test_unclustered_protein_gets_a_singleton_cluster(mini_xlsx: Path) -> None:
    rows = build_tomato_pu_rows(
        mini_xlsx,
        MINI_PROTEOME,
        arena=ARENA_PANEL,
        ratio=20,
        seed=12345,
    )
    folded = assign_grouped_folds(rows, {}, n_folds=2, seed=20260810)
    assert all(r.cluster_id == f"__singleton__{r.protein_accession}" for r in folded)


def test_grouped_folds_are_deterministic_and_cover_every_row(
    mini_xlsx: Path,
) -> None:
    rows = build_tomato_pu_rows(
        mini_xlsx,
        MINI_PROTEOME,
        arena=ARENA_PROTEOME,
        ratio=2,
        seed=12345,
    )
    clusters = {"PANEL1": "cl_A", "BACK1": "cl_A", "PANEL2": "cl_B", "BACK2": "cl_C"}
    first = assign_grouped_folds(rows, clusters, n_folds=2, seed=20260810)
    second = assign_grouped_folds(rows, clusters, n_folds=2, seed=20260810)
    assert [(r.protein_accession, r.cys_position, r.fold) for r in first] == [
        (r.protein_accession, r.cys_position, r.fold) for r in second
    ]
    assert len(first) == len(rows)
    assert {r.fold for r in first} == {0, 1}


def test_every_fold_holds_at_least_one_positive(mini_xlsx: Path) -> None:
    """Average precision is undefined on a fold with no positive, so the
    dealer must spread positive-bearing clusters across folds."""
    rows = build_tomato_pu_rows(
        mini_xlsx,
        MINI_PROTEOME,
        arena=ARENA_PROTEOME,
        ratio=2,
        seed=12345,
    )
    clusters = {"PANEL1": "cl_A", "BACK1": "cl_A", "PANEL2": "cl_B", "BACK2": "cl_C"}
    folded = assign_grouped_folds(rows, clusters, n_folds=2, seed=20260810)
    for fold in (0, 1):
        assert any(r.label == "positive" and r.fold == fold for r in folded), (
            f"fold {fold} has no positive"
        )


def test_unknown_arena_is_rejected(mini_xlsx: Path) -> None:
    with pytest.raises(ValueError):
        build_tomato_pu_rows(
            mini_xlsx,
            MINI_PROTEOME,
            arena="whole_genome",
            ratio=2,
            seed=1,
        )
