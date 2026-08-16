"""Scientific-input tests for site-biology feature extraction."""

import math
from pathlib import Path

from plantpersulf.features.sequence import (
    SequenceFeatureRow,
    _flanking_window,
    _hydrophobicity,
    _local_positive_charge_density,
)
from plantpersulf.features.site_biology import (
    BIOLOGY_FEATURE_NAMES,
    build_site_biology_vector,
    build_site_biology_vector_from_sequence,
)
from plantpersulf.provenance.audit import assert_registered_input
from plantpersulf.provenance.hashing import hash_file

REFERENCE_PATH = Path("data/registry/cache/uniprot/Q9ZW96.fasta")
REFERENCE_REGISTRY = Path("data/registry/reference_sequences.tsv")
REFERENCE_SHA256 = "f4a9ba53775e0433e3efe92599d023a5f54f8cf8995806cba42db27a3bb78b19"


def test_registered_q9zw96_yields_six_finite_core_features() -> None:
    """Feature extraction is supported by a registered public sequence."""
    assert_registered_input(REFERENCE_PATH, REFERENCE_REGISTRY)
    assert hash_file(REFERENCE_PATH, "sha256") == REFERENCE_SHA256

    sequence = "".join(REFERENCE_PATH.read_text(encoding="utf-8").splitlines()[1:])
    cys_position = 5  # PXD006140 coordinate-validated Q9ZW96 site.
    radius = 10
    window = _flanking_window(sequence, cys_position, radius)
    row = SequenceFeatureRow(
        protein_accession="Q9ZW96",
        cys_position=cys_position,
        label="psm_coordinate_only",
        flanking_window=window,
        hydrophobicity=_hydrophobicity(window),
        cys_density=sequence.count("C") / len(sequence),
        protein_length=len(sequence),
        local_positive_charge_density=_local_positive_charge_density(window),
    )

    from_row = build_site_biology_vector(row)
    from_sequence = build_site_biology_vector_from_sequence(
        "Q9ZW96", cys_position, sequence, radius=radius
    )

    assert from_row.names == BIOLOGY_FEATURE_NAMES
    assert from_row == from_sequence
    assert len(from_row.values) == 6
    assert all(math.isfinite(value) for value in from_row.values)
