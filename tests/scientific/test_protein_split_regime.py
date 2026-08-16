"""Random protein-level split — the literature-comparable regime.

Sul-BertGRU (Bioinformatics 2025, btaf078) evaluates with a random
protein-level 80/20 split and 10 repetitions; pCysMod-style tools use a
within-integrated-dataset split. ``partition_random_protein`` mirrors that
regime so our literature-comparable number uses the identical split
geometry: ``test_ratio`` of *proteins* held out entirely (all their Cys
rows, positive and unlabeled), NO homology control (deliberately — matching
the published regime). It is NOT a cross-study split and must never enter
Gate 2 (see tests/release/test_gate2_ignores_within_dataset_split_metrics.py).
"""

from __future__ import annotations

from plantpersulf.evaluation.external_validation import partition_random_protein


def _rows(n_proteins: int = 10) -> list[dict[str, str]]:
    rows = []
    for p in range(n_proteins):
        protein = f"P{p:02d}"
        rows.append(
            {
                "protein_accession": protein,
                "cys_position_in_protein": "10",
                "label": "positive",
                "study_accession": "PXD_A",
            }
        )
        rows.append(
            {
                "protein_accession": protein,
                "cys_position_in_protein": "20",
                "label": "unlabeled",
                "study_accession": "",
            }
        )
    return rows


def test_all_rows_of_a_protein_stay_together() -> None:
    train, test = partition_random_protein(_rows(), seed=0)
    train_proteins = {r["protein_accession"] for r in train}
    test_proteins = {r["protein_accession"] for r in test}
    assert train_proteins.isdisjoint(test_proteins)
    # both rows of every protein ended up in exactly one fold
    for p in range(10):
        assert f"P{p:02d}" in train_proteins or f"P{p:02d}" in test_proteins


def test_approximately_20_percent_of_proteins_held_out() -> None:
    rows = _rows(n_proteins=100)
    train, test = partition_random_protein(rows, seed=1, test_ratio=0.2)
    n_test = len({r["protein_accession"] for r in test})
    assert n_test == 20


def test_partition_is_deterministic_per_seed() -> None:
    rows = _rows(n_proteins=50)
    t1a, t1b = partition_random_protein(rows, seed=7)
    t2a, t2b = partition_random_protein(rows, seed=7)
    assert {r["protein_accession"] for r in t1a} == {
        r["protein_accession"] for r in t2a
    }
    assert {r["protein_accession"] for r in t1b} == {
        r["protein_accession"] for r in t2b
    }


def test_different_seeds_give_different_splits() -> None:
    rows = _rows(n_proteins=100)
    s0_train, _ = partition_random_protein(rows, seed=0)
    s1_train, _ = partition_random_protein(rows, seed=1)
    assert {r["protein_accession"] for r in s0_train} != {
        r["protein_accession"] for r in s1_train
    }
