"""Hash registration for the frozen tomato model-input snapshot."""

from pathlib import Path

from plantpersulf.provenance.audit import assert_registered_input


def test_tomato_snapshot_and_clusters_are_hash_registered() -> None:
    registry = Path("data/registry/model_inputs.tsv")

    assert_registered_input(
        Path("data/raw/references/tomato_ref_proteome_v1.fasta"), registry
    )
    assert_registered_input(
        Path("data/processed/clusters/tomato_proteome_clusters_v1.tsv"), registry
    )
