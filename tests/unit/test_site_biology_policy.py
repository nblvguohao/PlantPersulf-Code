"""Policy tests for the shortcut-free site-biology feature contract."""

from plantpersulf.features.site_biology import BIOLOGY_FEATURE_NAMES


def test_biology_feature_schema_exposes_only_the_six_core_features() -> None:
    """The public schema cannot grow shortcut features without review."""
    assert BIOLOGY_FEATURE_NAMES == (
        "hydrophobicity",
        "protein_cys_density",
        "local_positive_charge_density",
        "local_negative_charge_density",
        "local_cys_density",
        "local_sequence_entropy",
    )
    assert {
        "protein_length",
        "plddt",
        "has_structure",
        "study_id",
    }.isdisjoint(BIOLOGY_FEATURE_NAMES)
