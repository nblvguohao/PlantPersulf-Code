# Research canon

## Traceable project facts

1. The committed tomato v2 configuration uses the panel arena as primary, repeated five-fold homology-cluster CV with five repetitions, model seeds 0–4, and `development_k=50` (`configs/experiments/tomato_ranker_v2.yaml`). Every real run must bind the exact execution-time config hash and code revision; an approved acceleration-only revision is a new executable version, not an invisible change to an active run.
2. Tomato v2 compares an additive PU candidate with a PU-logistic baseline and falls back to the baseline unless all three paired metrics pass the conservative admission rule (`src/plantpersulf/evaluation/model_admission.py`).
3. The current v2 claim class is `development_candidate_ranking_not_gate2`; its outputs cannot change Gate 2 (`configs/experiments/tomato_ranker_v2.yaml`).
4. The tomato candidate universe is built from a registered frozen proteome and MMseqs2 cluster table with SHA256 records (`data/registry/model_inputs.tsv`).
5. KIAE271 contributes 99 coordinate-verified tomato sites over 88 proteins, while unverified rows are excluded rather than recycled as negatives (`docs/tomato_local_model_v1_2026-08-10.md`).
6. The repository already defines a frozen ESM-2 650M baseline interface, but model weights and any derived embedding cache must be independently provenance-locked before a new scientific run (`configs/experiments/baseline_esm2_v1.yaml`).
7. Sul-BertGRU uses iCysMod positives and putative negatives with confident-learning cleanup; it is a literature comparator, not a source of experimental negatives for this project (Wei et al., 2025, DOI: `10.1093/bioinformatics/btaf078`).
8. PU learning for PTM prediction predates this project, so the project must not claim the first generic PTM-PU method (He et al., 2017, PMID: `28872627`).

## Definitions

- `B_v2`: the model selected by the frozen tomato v2 admission policy, either additive PU or its PU-logistic fallback.
- `P1`: frozen PLM representation plus a regularized linear nnPU/pairwise head.
- `P2`: interpretable biological score plus a gated PLM residual trained under the same PU objective.
- `panel arena`: other cysteines in the KIAE271-detected proteins; the primary development arena.
- `proteome arena`: sampled cysteines from the whole tomato proteome; an audit and deployment-pressure arena, not a model-selection arena.
- `blind release`: a later immutable, hash-addressed candidate package whose rankings cannot be changed after wet-lab results.

## Forbidden claims

- The current or proposed model is a general cross-crop predictor.
- A high model score proves persulfidation, function, mechanism, or causality.
- PLM, AlphaFold, docking, or MD output is experimental evidence.
- Unobserved cysteines are biological negatives.
- A v3 model is superior before paired homology-blocked evaluation passes.
- Nature Communications readiness follows from algorithmic performance alone.

## Unresolved items

- The outcome and admission decision of the ongoing v2 run.
- Whether the exact ESM-2 weights used locally have a complete accession/version/license/SHA256 record.
- Whether a PLM residual adds stable panel-arena value over `B_v2`.
- Wet-lab capacity, expected control hit rate, target enrichment, and therefore the final candidate count.
- Whether a target-label-free cross-crop arm will pass its separate computational and power gates.
