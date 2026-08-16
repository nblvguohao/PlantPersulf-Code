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
        uniprot_accession="A0A3Q7EW23",
        cys_position=206,
        doi="10.1093/plphys/kiad070",
        control_type="conditional_site_group_control",
        status="mapped",
        control_species="Solanum lycopersicum",
        in_benchmark_as="absent",
        provenance=(
            "Resolved (2026-07-28): the prior gene-name search for "
            "'BRG3' matched only the unrelated 143aa K4BRG3 entry. The "
            "paper's own Accession numbers section gives NCBI GeneID "
            "LOC101267168 ('probable BOI-related E3 ubiquitin-protein "
            "ligase 3', Solanum lycopersicum chr1). That GeneID maps to "
            "UniProt A0A3Q7EW23 (RING-type domain-containing protein, "
            "GN=LOC101267168, 243aa) in the tomato reference proteome. "
            "Both Cys206 and Cys212 are verified as true Cys residues "
            "at the expected positions in this sequence "
            "(...KSCNSRSSCMICLPCRH... and ...SSCMICLPCRHLSSCKT..., "
            "1-based); the paper's BRG3-Cys-Ala mutant targets both "
            "sites jointly (conditional_site_group_control). Cys206 is "
            "registered as the representative scored position, "
            "consistent with the single-position-per-lineage "
            "convention used for ERF.D3's paired Cys115/Cys118 site "
            "group below; Cys212 is the second site in the same group "
            "and is not scored separately."
        ),
    ),
    KnownControl(
        mechanism_lineage_id="WRKY71_H2S_UBIQUITINATION",
        gene="SlWRKY71",
        uniprot_accession="A0A3Q7FNU4",
        cys_position=35,
        doi="10.1093/plphys/kiad070",
        control_type="conditional_site_group_control",
        status="position_shift",
        control_species="Solanum lycopersicum",
        in_benchmark_as="absent",
        provenance=(
            "Resolved (2026-08-10): the paper's own Accession numbers "
            "section gives NCBI GeneID LOC101264783 ('WRKY transcription "
            "factor 71-like', Solanum lycopersicum chr2) — prior public "
            "gene-name search for 'WRKY71' had returned zero hits because "
            "it never checked this section (same root cause as the BRG3 "
            "fix above). NCBI RefSeq XP_004233015.1 (317aa) is the linked "
            "full-length protein; its sequence contains the paper's "
            "antibody antigen peptide 'CQVKKRVERSYQDP' verbatim at "
            "residues 199-211, and its Cys193/Cys198 match the paper's "
            "WRKY71-Cys193Ala-Cys198Ala double-mutant sites exactly with "
            "zero offset. UniProt's own GeneID/RefSeq cross-reference "
            "tables do not yet index this locus (UniProt ID mapping and "
            "REST gene-name/xref search both return zero hits), but the "
            "tomato reference proteome v1 fasta already used for scoring "
            "contains A0A3Q7FNU4 (159aa, 'WRKY domain-containing "
            "protein'), whose sequence is a 100%-identical, contiguous "
            "match to RefSeq residues 159-317 (an N-terminally truncated "
            "automatic gene model of the same locus, not a different "
            "gene). In A0A3Q7FNU4's own numbering this is Cys35/Cys40 "
            "(offset -158 from the RefSeq/paper numbering: 193-158=35, "
            "198-158=40). Cys35 is registered as the representative "
            "scored position, consistent with the BRG3/ERF.D3 "
            "single-position-per-lineage convention; Cys40 is the second "
            "site in the same paired-mutant group and is not scored "
            "separately. Unlike ERF.D3's unexplained 13-residue shift, "
            "this offset is fully explained (N-terminal truncation of "
            "the same locus) and independently corroborated via NCBI "
            "RefSeq — status is 'position_shift' rather than 'mapped' "
            "only because the scored UniProt entry's numbering differs "
            "from the paper's, not because of residual doubt about "
            "identity."
        ),
    ),
    KnownControl(
        mechanism_lineage_id="RNF144B_H2S_UBIQUITINATION",
        gene="RNF144b",
        uniprot_accession="A0A3Q7GXU6",
        cys_position=122,
        doi="10.1093/plphys/kiad070",
        control_type="ms_detected_no_functional_validation",
        status="mapped",
        control_species="Solanum lycopersicum",
        in_benchmark_as="absent",
        provenance=(
            "Resolved (2026-08-10) from the paper's own Supplemental "
            "Table S3 ('Identification of the persulfidation peptide by "
            "LC-MS/MS'), not from the main text (which names RNF-144b "
            "only as a second E3 ligase detected alongside BRG3, with no "
            "position given there). Table S3 row: peptide "
            "'FYCPYKDCSAMLVNDSDEIVR', Protein Group Accession "
            "XP_004242195.1, Modifications 'C3(S); M11(Oxidation)' "
            "(1-based within-peptide numbering; position 11 lands "
            "exactly on the peptide's only Met, cross-checking the "
            "numbering convention). NCBI RefSeq XP_004242195.1 ('E3 "
            "ubiquitin-protein ligase RSL1-like', 233aa) confirmed via "
            "GeneID LOC101265447 esearch; the peptide is an exact, "
            "unique substring at protein residues 120-141, placing the "
            "modified Cys (in-peptide position 3) at protein position "
            "122. The tomato reference proteome v1 fasta contains a "
            "longer isoform of the same locus, A0A3Q7GXU6 (318aa, 'RBR-"
            "type E3 ubiquitin transferase'), whose first 233 residues "
            "are identical to the RefSeq entry — Cys122 is confirmed "
            "identical in both numbering systems (zero offset), so "
            "status is 'mapped', not 'position_shift'. Unlike BRG3 and "
            "WRKY71 in this same paper, RNF-144b's persulfidation was "
            "NOT functionally validated by an Ala-substitution mutant — "
            "only detected by LC-MS/MS on the recombinant protein and "
            "tested (unsuccessfully, weaker signal than BRG3) for "
            "protein-protein interaction with WRKY71 by luciferase "
            "complementation. control_type reflects this lower "
            "evidentiary tier honestly."
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
        mechanism_lineage_id="CAT1_H2S_CUONP_OXIDATIVE",
        gene="CAT1",
        uniprot_accession="P30264",
        cys_position=234,
        doi="10.1016/j.plaphy.2020.09.020",
        control_type="strong_single_site_control",
        status="mapped",
        control_species="Solanum lycopersicum",
        in_benchmark_as="absent",
        provenance=(
            "Resolved (2026-08-10) via Corpas et al. COPLBI-D-26-00068 "
            "(in-review Current Opinion in Plant Biology review, Table 1) "
            "citing Li et al. 2020 (Plant Physiol Biochem 156:257-266) — "
            "a third independent lab (Li/Shi/Wang/Liao), distinct from "
            "both the Seville Romero/Gotor/Aroca training network and "
            "Zhang Hua's HFUT lab. UniProt P30264 (CATA1_SOLLC, reviewed "
            "Swiss-Prot entry, 492aa) gene name 'CAT1' resolves cleanly "
            "via public gene-name search (unlike the tomato TrEMBL "
            "automatic-annotation entries this project has repeatedly "
            "had to resolve via GeneID cross-reference instead); residue "
            "234 verified as Cys with zero offset from the paper's "
            "reported position. Ten-day-old tomato cv. Liger seedlings, "
            "CuO-nanoparticle oxidative stress ± NaHS: persulfidation of "
            "CAT1 Cys234 DECREASES catalase activity (contrast with "
            "APX1/POD5 from the same paper, whose persulfidation "
            "increases their activity — not registered here because "
            "their UniProt accessions could not be independently "
            "resolved: public gene-name search for tomato 'APX1' returns "
            "candidates A0A3Q7GLU9/Q3I5C4, but position 168 in both is "
            "Ala, not Cys — a real mismatch, not merely an unindexed "
            "GeneID as with BRG3/WRKY71/RNF-144b; 'POD5' returns zero "
            "Solanum lycopersicum hits. Both need the paper's own "
            "Accession numbers section, which is not available from this "
            "review's citation alone.)."
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
    KnownControl(
        mechanism_lineage_id="PAD3_H2S_HCN_OSMOTIC",
        gene="PAD3",
        uniprot_accession="Q9LW27",
        cys_position=440,
        doi="10.1111/pce.70593",
        control_type="strong_single_site_control",
        status="mapped",
        control_species="Arabidopsis thaliana",
        in_benchmark_as="unlabeled",
        provenance=(
            "Arabidopsis reference proteome v1 — Q9LW27 (C71BF_ARATH, "
            "camalexin synthase CYP71B15/PAD3, 490aa) residue 440 "
            "verified as Cys; benchmark row (Q9LW27,440) is unlabeled, "
            "not a training positive. Zhang et al. 2026 (Plant, Cell & "
            "Environment, doi:10.1111/pce.70593; Pei/Jin lab — no "
            "author overlap with the Seville Romero/Gotor/Aroca network "
            "across their >20 independent H2S-signaling papers since "
            "2013): biotin-switch-assay-confirmed persulfidation of "
            "PAD3 at Cys440 enhances its HCN-synthase activity under "
            "osmotic stress. Second same-species, independent-lab "
            "control — scored against the Arabidopsis unlabeled "
            "reference with the control row excluded from the scorer's "
            "training sample."
        ),
    ),
    KnownControl(
        mechanism_lineage_id="DES1_H2S_SELF_RBOHD_ABA",
        gene="DES1",
        uniprot_accession="F4K5T2",
        cys_position=44,
        doi="10.1105/tpc.19.00826",
        control_type="strong_single_site_control",
        status="mapped",
        control_species="Arabidopsis thaliana",
        in_benchmark_as="unlabeled",
        provenance=(
            "Resolved (2026-08-10) via Corpas et al. COPLBI-D-26-00068 "
            "review Table 1, citing Shen et al. 2020 (Plant Cell "
            "32(4):1000-1017, doi:10.1105/tpc.19.00826). UniProt F4K5T2 "
            "(CGL_ARATH, 'L-cysteine desulfhydrase 1 (DES1)', 323aa) "
            "resolves cleanly via public gene-name search; residue 44 "
            "verified as Cys, zero offset from the paper's reported "
            "position; benchmark row (F4K5T2,44) is unlabeled, not a "
            "training positive. DES1 is the H2S-generating enzyme itself "
            "in this guard-cell ABA-signaling story — under ABA, DES1 "
            "persulfidates itself at Cys44/Cys205 (Cys205 is the second "
            "site in the same paired report and is not scored "
            "separately, matching the BRG3/WRKY71 single-position-per-"
            "lineage convention) as well as RBOHD (registered "
            "separately below, same paper). Same-species, "
            "training-independent control — scored against the "
            "Arabidopsis unlabeled reference with the control row "
            "excluded from the scorer's training sample."
        ),
    ),
    KnownControl(
        mechanism_lineage_id="RBOHD_H2S_ROS_ABA",
        gene="RBOHD",
        uniprot_accession="Q9FIJ0",
        cys_position=825,
        doi="10.1105/tpc.19.00826",
        control_type="strong_single_site_control",
        status="mapped",
        control_species="Arabidopsis thaliana",
        in_benchmark_as="unlabeled",
        provenance=(
            "Resolved (2026-08-10), same source/paper as DES1 above "
            "(Shen et al. 2020, Plant Cell, doi:10.1105/tpc.19.00826). "
            "UniProt Q9FIJ0 (RBOHD_ARATH, 'Respiratory burst oxidase "
            "homolog protein D', 921aa) resolves cleanly via public "
            "gene-name search; residue 825 verified as Cys, zero offset. "
            "DES1 persulfidates RBOHD at Cys825/Cys890 (Cys890 is the "
            "second site in the paired report, not scored separately); "
            "under high ROS, RBOHD/DES1 become persulfide-oxidised "
            "(-SSOnH), desensitising ABA signalling — a feedback loop "
            "reducible by thioredoxin. Same-species, training-"
            "independent control, same exclusion treatment as DES1."
        ),
    ),
    KnownControl(
        mechanism_lineage_id="SNRK26_H2S_PHOSPHORYLATION",
        gene="SnRK2.6",
        uniprot_accession="Q940H6",
        cys_position=131,
        doi="10.1016/j.molp.2021.07.002",
        control_type="strong_single_site_control",
        status="mapped",
        control_species="Arabidopsis thaliana",
        in_benchmark_as="unlabeled",
        provenance=(
            "Resolved (2026-08-10) via Corpas et al. COPLBI-D-26-00068 "
            "review Table 1, citing Chen et al. 2020 (Mol Plant "
            "13:732-744, doi:10.1016/j.molp.2020.01.004 — first "
            "discovery) and Chen et al. 2021 (Mol Plant 14(11):1814-1830, "
            "doi:10.1016/j.molp.2021.07.002 — mechanistic follow-up, "
            "used here as the primary DOI). UniProt Q940H6 (SRK2E_ARATH, "
            "'SNF1-related kinase 2.6 / OST1', 362aa) resolves cleanly "
            "via public gene-name search; residue 131 verified as Cys, "
            "zero offset. **This is the single most important literature "
            "precedent for the project's PTM-crosstalk-grammar phase "
            "(docs/superpowers/plans/2026-07-28-ptm-crosstalk-grammar.md): "
            "unlike SlWRKY6/SlERF.D2 (functional antagonism between "
            "persulfidation and phosphorylation at DIFFERENT residues), "
            "SnRK2.6 is a directly-demonstrated case of persulfidation "
            "(Cys131/Cys137) and phosphorylation (Ser175/Ser267) on the "
            "SAME protein mutually influencing each other via an "
            "intramolecular structural change (Chen et al. 2021's title "
            "claim) — Cys137 persulfidation facilitates Ser175 "
            "phosphorylation, while Ser267 phosphorylation positively "
            "regulates Cys137 persulfidation. Cys137 is the second site "
            "in the paired report and is not scored separately.** Same-"
            "species, training-independent control, same exclusion "
            "treatment as DES1/RBOHD above."
        ),
    ),
    KnownControl(
        mechanism_lineage_id="ABI4_H2S_MAPKKK18_ABA",
        gene="ABI4",
        uniprot_accession="A0MES8",
        cys_position=250,
        doi="10.1016/j.molp.2021.03.007",
        control_type="strong_single_site_control",
        status="mapped",
        control_species="Arabidopsis thaliana",
        in_benchmark_as="unlabeled",
        provenance=(
            "Resolved (2026-08-10) via Corpas et al. COPLBI-D-26-00068 "
            "review Table 1, citing Zhou et al. 2021 (Mol Plant "
            "14(6):921-936, doi:10.1016/j.molp.2021.03.007). UniProt "
            "A0MES8 (ABI4_ARATH, 'Ethylene-responsive transcription "
            "factor ABI4 / Protein ABSCISIC ACID INSENSITIVE 4', 328aa) "
            "resolves cleanly via public gene-name search; residue 250 "
            "verified as Cys, zero offset. DES1-dependent persulfidation "
            "of ABI4 at Cys250 enhances its DNA-binding to the "
            "MAPKKK18 promoter, amplifying ABA-MAPK-cascade signalling — "
            "a second independent DES1-target report from a different "
            "author group/paper than the DES1/RBOHD entry above (Zhou "
            "et al., not Shen et al.), though within the same Seville "
            "Gotor/Romero network. Same-species, training-independent "
            "control, same exclusion treatment as DES1/RBOHD/SnRK2.6."
        ),
    ),
    KnownControl(
        mechanism_lineage_id="ATG4A_H2S_AUTOPHAGY_ABA",
        gene="ATG4a",
        uniprot_accession="Q8S929",
        cys_position=170,
        doi="10.1105/tpc.20.00766",
        control_type="strong_single_site_control",
        status="mapped",
        control_species="Arabidopsis thaliana",
        in_benchmark_as="unlabeled",
        provenance=(
            "Resolved (2026-08-10) via Corpas et al. COPLBI-D-26-00068 "
            "review Table 1, citing Laureano-Marin et al. 2020 (Plant "
            "Cell 32(12):3902-3920, doi:10.1105/tpc.20.00766). The "
            "paper's title says only 'the Cys Protease ATG4' — Arabidopsis "
            "has two paralogs, ATG4a (Q8S929, 467aa) and ATG4b (Q9M1Y0, "
            "477aa); disambiguated by checking residue 170 in both: "
            "ATG4a position 170 is Cys (context "
            "'...SDVNWGC[170]MIRSSQ...'), ATG4b position 170 is Asn — "
            "only ATG4a is consistent with the paper's reported Cys170, "
            "confirming the isoform. Endogenous H2S negatively regulates "
            "autophagy by maintaining ATG4a persulfidated at Cys170 "
            "(inhibiting its protease activity on ATG8); ABA lowers "
            "persulfidation, activating ATG4a and triggering "
            "autophagosome formation. Same-species, training-independent "
            "control, same exclusion treatment as the other four "
            "Corpas-review Arabidopsis entries above."
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
