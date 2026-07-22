"""Registered known-mechanism persulfidation controls (master TDD sect. 3.3).

Each control is a published, independently validated persulfidation site whose
rank in the trained model is reported transparently — including failed and
unmappable controls. Two scope classes exist:

* **cross-species** (tomato, ``in_benchmark_as="absent"``): training positives
  are exclusively Arabidopsis, so recovery cannot be a training artefact; the
  control is scored from the tomato reference proteome.
* **same-species** (Arabidopsis, ``in_benchmark_as="unlabeled"``): the control
  protein IS in the benchmark, with the control site present as an unlabeled
  PU-pool row — never a training positive. Before scoring, the control row
  must be excluded from the scorer's training sample (``filter_out_keys``);
  the control is scored from the Arabidopsis reference proteome. This is
  declared here so no consumer can mistake it for a cross-species control.
"""

from __future__ import annotations

from collections.abc import Collection, Iterable
from dataclasses import dataclass


@dataclass(frozen=True)
class KnownControl:
    mechanism_lineage_id: str
    gene: str
    uniprot_accession: str
    cys_position: int
    doi: str
    control_type: str
    status: str  # mapped | unmappable | position_shift
    control_species: str
    in_benchmark_as: str  # absent | unlabeled
    provenance: str


REGISTERED_CONTROLS: tuple[KnownControl, ...] = (
    KnownControl(
        mechanism_lineage_id="SLWRKY6_H2S_PHOSPHORYLATION",
        gene="SlWRKY6",
        uniprot_accession="A0A3Q7F586",
        cys_position=396,
        doi="10.1093/plphys/kiae271",
        control_type="strong_single_site_control",
        status="mapped",
        control_species="Solanum lycopersicum",
        in_benchmark_as="absent",
        provenance=(
            "UniProt tomato reference proteome v1 — peptide position verified as Cys"
        ),
    ),
    KnownControl(
        mechanism_lineage_id="SLERFD2_H2S_ETHYLENE",
        gene="SlERF.D2",
        uniprot_accession="A0A3Q7JX06",
        cys_position=35,
        doi="10.1111/tpj.70000",
        control_type="strong_single_site_control",
        status="mapped",
        control_species="Solanum lycopersicum",
        in_benchmark_as="absent",
        provenance=(
            "UniProt tomato reference proteome v1 — peptide position verified as Cys"
        ),
    ),
    KnownControl(
        mechanism_lineage_id="BRG3_H2S_UBIQUITINATION",
        gene="BRG3",
        uniprot_accession="",
        cys_position=0,
        doi="10.1093/plphys/kiad070",
        control_type="conditional_site_group_control",
        status="unmappable",
        control_species="Solanum lycopersicum",
        in_benchmark_as="absent",
        provenance=(
            "UniProt tomato reference proteome v1 (36988 records): "
            "no entry matching gene name BRG3 found at the expected "
            "protein length (Cys206/Cys212 positions exceed the only "
            "matched entry K4BRG3, which is 143aa)."
        ),
    ),
    KnownControl(
        mechanism_lineage_id="ERFD3_H2S_CONTEXT",
        gene="ERF.D3",
        uniprot_accession="A0A3Q7ESP9",
        cys_position=128,
        doi="10.1093/plphys/kiae560",
        control_type="conditional_site_group_control",
        status="position_shift",
        control_species="Solanum lycopersicum",
        in_benchmark_as="absent",
        provenance=(
            "UniProt tomato reference proteome v1: mapped to "
            "A0A3Q7ESP9 (AP2/ERF domain-containing protein, GN=ERF-D3, "
            "316aa, 4 Cys at [128,131,136,139]). The paper reports "
            "Cys115/Cys118; the 13-residue offset may reflect a "
            "signal peptide, an isoform numbering difference, or a "
            "post-translational cleavage — MUST be confirmed with the "
            "authors before any recovery claim."
        ),
    ),
    KnownControl(
        mechanism_lineage_id="ATG6PD6_H2S_G6PD_SALT",
        gene="AtG6PD6",
        uniprot_accession="Q9FJI5",
        cys_position=159,
        doi="10.1111/nph.19188",
        control_type="strong_single_site_control",
        status="mapped",
        control_species="Arabidopsis thaliana",
        in_benchmark_as="unlabeled",
        provenance=(
            "Arabidopsis reference proteome v1 — Q9FJI5 (G6PD6_ARATH, "
            "515aa) residue 159 verified as Cys; benchmark row "
            "(Q9FJI5,159) is unlabeled, not a training positive. "
            "Li et al. 2023 (New Phytologist, NWAFU; iProX PXD043969): "
            "MS + mutagenesis-validated persulfidation at AtG6PD6 Cys159 "
            "under salt stress. Same-species, independent-lab control — "
            "scored against the Arabidopsis unlabeled reference with the "
            "control row excluded from the scorer's training sample."
        ),
    ),
)


def scores_against_benchmark_proteome(control: KnownControl) -> bool:
    """Same-species (Arabidopsis) controls are scored from the benchmark
    reference proteome; cross-species controls from their own."""
    return control.control_species == "Arabidopsis thaliana"


def filter_out_keys(
    rows: Iterable[dict[str, str]],
    keys: Collection[tuple[str, int]],
) -> list[dict[str, str]]:
    """Drop rows whose (protein_accession, cys_position) is in ``keys``.

    Used to remove same-species control rows from a PU scorer's training
    unlabeled sample, so the control is scored as a genuinely held-out
    (never-trained-on) unlabeled row. Order is preserved.
    """
    return [
        row
        for row in rows
        if (row["protein_accession"], int(row["cys_position_in_protein"])) not in keys
    ]
