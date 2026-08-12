"""Fail-closed materialization of an augmented registered proteome."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from plantpersulf.proteomics.reference_proteome import (
    RegisteredFastaInput,
    materialize_augmented_reference_proteome,
)
from plantpersulf.provenance.audit import assert_registered_input


def test_materialization_records_hash_verified_base_and_isoform_inputs(
    tmp_path: Path,
) -> None:
    """A missing source hash would make the v2 proteome non-traceable."""
    base = tmp_path / "base.fasta"
    isoform = tmp_path / "isoform.fasta"
    base.write_text(">sp|MARKER1|base\nACDE\n", encoding="utf-8")
    isoform.write_text(">sp|MARKER2-2|isoform\nCDEA\n", encoding="utf-8")
    output = tmp_path / "reference_v2.fasta"
    manifest = tmp_path / "reference_v2.manifest.json"

    result = materialize_augmented_reference_proteome(
        RegisteredFastaInput(base, hashlib.sha256(base.read_bytes()).hexdigest()),
        (
            RegisteredFastaInput(
                isoform, hashlib.sha256(isoform.read_bytes()).hexdigest()
            ),
        ),
        output,
        manifest,
    )

    assert result.output_sha256 == hashlib.sha256(output.read_bytes()).hexdigest()
    assert output.read_text(encoding="utf-8") == (
        ">sp|MARKER1|base\nACDE\n>sp|MARKER2-2|isoform\nCDEA\n"
    )
    recorded = json.loads(manifest.read_text(encoding="utf-8"))
    assert recorded["input_sha256"] == {
        str(base): hashlib.sha256(base.read_bytes()).hexdigest(),
        str(isoform): hashlib.sha256(isoform.read_bytes()).hexdigest(),
    }
    assert recorded["output_sha256"] == result.output_sha256


def test_materialized_arabidopsis_v2_proteome_is_hash_registered() -> None:
    """The v2 global clustering input must resolve in the input registry."""
    assert_registered_input(
        Path("data/raw/references/arabidopsis_ref_proteome_v2.fasta"),
        Path("data/registry/model_inputs.tsv"),
    )
