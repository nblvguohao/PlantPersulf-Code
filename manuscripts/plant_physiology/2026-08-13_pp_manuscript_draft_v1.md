<!--
====================================================================
DRAFT STATUS — REMOVE THIS BLOCK BEFORE SUBMISSION
====================================================================
Draft v1, 2026-08-13. Target journal: Plant Physiology (Research Article).
Strategy context: this is the pre-registered downgrade-path manuscript of
docs/superpowers/plans/2026-08-13-nc-entry-gate-roadmap.md. The same draft
carries the Nature Communications path if Gates 3-4 pass; the Plant
Physiology version stands on Gates 0-2 + Gate 3 regardless of direction.

PENDING before submission (in roadmap order):
  1. ~~Co-signature of sap_protocol.json / analysis_plan.json~~ DONE 2026-08-13
     (both parties confirmed; audit --expect-signed PASS; git tag
     gate0-release-v1-20260813). Remaining: preregistration submission -> DOI
     [insert in R5 and Data Availability]
  2. Gate 1: strict-track frozen test (one-shot) -> insert numbers in R2
  3. Gate 3 blind cohort -> complete R6/R7, Table 5, Figure 6, Abstract
     outcome sentence (two pre-written variants kept below)
  4. Gate 6: public DOI for results/data_archive/nc_entry_gates/
  5. Author list, affiliations, ORCIDs, funding, acknowledgements
  6. All [VERIFY ...] reference entries
  7. Figures 1-5 rendered from the scripts listed in each legend

WORDING GATE (binding): every outward text must pass
verify_predictive_claims -> [] against the STOP decision
(plantpersulf.evaluation.conclusion_gate). This file was written under
that discipline; re-run the gate after ANY edit. Negations of forbidden
phrases also trip the matcher — reword, never negate.
====================================================================
-->

# A frozen, preregistered framework for proteome-wide ranking of candidate cysteine persulfidation sites in tomato (*Solanum lycopersicum*)

**Running title:** Preregistered ranking of persulfidation candidates

**One-sentence summary:** We freeze a multispecies-trained site-ranking model, release a fully hash-registered tomato candidate table with matched controls, calibrate expectations against published sites, and preregister a blind experimental test whose outcome will be reported regardless of direction.

**Authors:** Lu Guohao^1,\*, [Zhang Hua]^2,\*, [additional authors TBD]

^1 PlantPersulf project, [affiliation TBD]
^2 [South China Agricultural University, Guangzhou, China — VERIFY full affiliation]

\* Corresponding authors: [emails TBD]

ORCIDs: [TBD]

---

## Abstract

Protein persulfidation, the hydrogen sulfide (H2S)-driven conversion of cysteine thiols to persulfides, is an emerging redox regulatory mechanism in plants, but site-level experimental mapping remains costly and low-throughput. Computational ranking of candidate sites could focus experimental effort, yet published tools are scarce, their training data are difficult to obtain or reuse, and performance claims rest almost exclusively on within-dataset splits. Here we present a transparent, preregistered alternative. We integrated a provenance-tracked, four-species dataset of 389,609 cysteine sites (2,334 experimentally supported positives) and benchmarked six models under literature-comparable random-protein splits; a gated-fusion positive-unlabeled (PU) ranker achieved the highest within-dataset macro average precision (0.386; next best 0.047). We then froze this model — weights, features, hyperparameters, candidate table, matched controls, and analysis plan — under SHA256 registration before any blind data existed, scored all 179,736 non-panel cysteine sites of the tomato reference proteome, and calibrated honest expectations using published tomato sites: previously unseen confirmed sites rank moderately above the list median but do not concentrate at the list top. A blind experimental cohort (Top-200 candidates plus 400 matched controls) with a preregistered one-sided Fisher exact test was designed with full power tabulation. **[OUTCOME SENTENCE — Gate 3; two pre-written variants in Results R6.]** All artifacts, negative results, and exclusion ledgers are openly archived. This work establishes a reproducible, preregistered pipeline for candidate ranking in plant redox proteomics and reports its blind test regardless of outcome.

*(Word count target ≤250; verify at finalization.)*

---

## Introduction

Hydrogen sulfide (H2S) has moved from toxic byproduct to recognized signaling molecule in plants, where it participates in stomatal movement, stress responses, autophagy, and fruit ripening [VERIFY refs: Aroca/Gotor reviews; Zhang et al. 2021]. Its best-characterized molecular mechanism is persulfidation (S-sulfhydration): the covalent conversion of reactive cysteine thiols (-SH) to persulfides (-SSH), which alters the activity, localization, or interactions of target proteins. Site-resolved mapping of persulfidation relies on mass-spectrometric chemoproteomics (e.g. tag-switch derivatives), and each candidate site then requires individual experimental confirmation — a costly bottleneck that limits the field to a few large datasets in Arabidopsis [PXD006140, PXD024061], rice [PXD072089], the fungal pathogen *Magnaporthe oryzae* [PXD063170], and, notably for this work, tomato, where a persulfidome linked to fruit ripening was reported by Zhang et al. (2024, *Plant Physiology*, doi:10.1093/plphys/kiae271; hereafter "kiae271").

Computational prioritization of candidate sites is an attractive complement, but the available tooling is thin. Sul-BertGRU [VERIFY ref] is the only published persulfidation-specific deep-learning tool with retrievable training data, and our audit of that data found substantial label problems (1,003 of 2,705 entries have a non-cysteine center residue; 78% human, 21% Arabidopsis sequences) [see companion audit: docs/competitor_data_audit_2026-08-10.md]; pCysMod's training data could not be obtained [VERIFY availability]; iCysMod is a database rather than a site-ranking tool [VERIFY ref]. No new persulfidation-specific predictors appeared in the 2025 literature as of our audit date. More fundamentally, performance claims across the cysteine-PTM tool literature rest on within-dataset random splits, and — as in much of computational biology — blinded, preregistered experimental evaluation is essentially absent.

We take the position that the correct response to this situation is not a larger claim but a stricter protocol. In prior work on Arabidopsis we showed that a sequence+structure PU ranker carries a real but modest within-dataset signal (leave-study-out pooled AP 0.102, 2.1× base rate; both studies from one laboratory), and we froze a binding limitation statement that governs all outward text from this project:

> Current public data are insufficient to demonstrate cross-study predictive ability; the model is used only for candidate organisation and hypothesis generation.

Accordingly, everything in the present paper is organized around a candidate-organisation tool whose value is decided by a preregistered blind experiment, not by retrospective metrics.

Here we report: (i) a provenance-tracked, four-species site-level dataset and a literature-comparable within-dataset benchmark in which a gated-fusion PU ranker leads six models; (ii) a fully frozen, hash-registered candidate release for the tomato proteome — model weights, feature schema, inference script, ranked candidate table, and matched controls; (iii) a pre-blind calibration against published tomato and Arabidopsis sites that sets honest expectations for the blind test; (iv) the preregistered blind-cohort design, power analysis, and statistical analysis plan, co-signed by the modeling and wet-laboratory partners on 2026-08-13; preregistration DOI **[TBD at submission]**; and (v) **[Gate 3: the blind-cohort result, reported regardless of outcome]**. All negative results, exclusions, and intermediate artifacts are preserved in a public data archive.

---

## Results

### R1. A provenance-tracked, four-species site-level persulfidation dataset

We integrated site-level persulfidation evidence from four public studies into a single development panel in which each row is one cysteine site with full study provenance (accession, coordinate, label, source study): PXD006140 and PXD024061 (Arabidopsis; same laboratory and chemistry), PXD072089 (rice; independent laboratory), PXD063170 (*M. oryzae*), and the kiae271 tomato sites. The release training panel — the union of the ten seed panels used for the literature-comparable benchmark (development split only; the strict-track frozen test partition was never touched) — comprises **389,609 sites with 2,334 positives** (SHA256 `30e432fb…`; Table 1 gives the per-species candidate-scan accounting; per-species positive counts in Supplemental Table S1 **[TODO: from dataset manifest]**). The panel is dominated by unlabeled sites (positive-unlabeled learning setting): experimentally supported positives are certain, but "unlabeled" pools true negatives with sites never assayed.

Two provenance facts shape everything downstream and are stated up front. First, the Arabidopsis positives come from one laboratory with one chemistry; no amount of modeling converts them into cross-study evidence. Second, all 79 tomato positives in the panel are the kiae271 published sites (verified by direct key overlap; see R4) — the model has seen every published tomato site there is.

### R2. Within-dataset benchmark under literature-comparable splits

We benchmarked six models on the integrated dataset under the split protocol customary in the cysteine-PTM tool literature (random protein-level holdout, 10 seeds): the gated-fusion PU ranker (`structure_ranker`), an ESM-2 linear head, XGBoost, a PU logistic baseline, random forest, and a retrained Sul-BertGRU. The reported metric is macro average precision (AP) averaged over the three crop species (Arabidopsis, rice, tomato). This is a **within-dataset** evaluation; per the project's claim gate it supports exactly one claim class, quoted verbatim:

> Literature-comparable within-dataset performance improvement.

`structure_ranker` achieved mean macro AP **0.386 ± 0.173** (10 seeds), 7.9× the mean of the best-per-seed baseline and leading on every seed (Table 2). Per-species APs were 0.493 ± 0.267 (Arabidopsis), 0.647 ± 0.257 (rice), and **0.0185 ± 0.006 (tomato)** — tomato is by far the hardest species within the dataset, a point we return to below. A parallel homology-cluster-split development track (5 folds) is maintained for model selection; its 20% frozen test partition remains locked and will be unlocked exactly once under audit before any blind-cohort unblinding **[Gate 1 — insert strict-track frozen-test numbers here when unlocked; one sentence plus table row]**.

None of these numbers enters any external-evidence claim. Their role in this paper is to justify the *choice* of the model that was frozen for the blind test.

### R3. A frozen, hash-registered candidate release for the tomato proteome

On 2026-08-13 we froze the complete candidate release `multispecies-v2-candidate-release-v1` (Table 3): the trained model bundle (gated fusion of three sequence features — local hydrophobicity, cysteine density, local positive charge density — plus an AlphaFold structure branch with explicit missingness masking; ESM and study-context branches disabled; hidden width 16, dropout 0.2, 200 epochs, lr 0.05, 16 MC-dropout passes; training seed fixed at the freeze date, not performance-selected), the feature schema, a self-contained inference script without training code, the ranked tomato candidate table, and matched controls.

The candidate table scores **179,736 tomato cysteine sites** — every cysteine of the SHA-pinned reference proteome except the 72,051 sites in the training panel. Scores span 0.0053–0.305 with an extremely flat top (the top 200 sites occupy 0.295–0.305), a direct consequence of positives being ~0.6% of training data. A structural finding of the release: the frozen AlphaFold registry contains **zero tomato accessions**, so every tomato site is scored with the structure branch masked (coverage 0/179,736). Diagnostic scans with the same frozen model quantify the contrast: Arabidopsis 284,387 candidates with 1.07% structure coverage, rice 469,748 with 0.75% (Table 1). Tomato ranking in this release is therefore driven entirely by the three sequence features — a feature-starvation state that the pre-blind calibration (R4) suggests is consequential, and that any future release may address only under a new release number and new preregistration.

For each of the top-200 candidates we pre-computed up to five matched controls (966 pairs total): non-panel tomato sites from the same structure-availability stratum, nearest by Euclidean distance in z-scored sequence-feature space, caliper 0.25 SD, no control reused. The blind analysis uses exactly the two nearest controls per candidate (rule locked at signing; the artifact itself is unchanged).

### R4. Pre-blind calibration: what published sites say to expect

Before any blind data existed, we scored the 99 coordinate-verified kiae271 sites with the frozen bundle and ranked them within the frozen candidate table — a "mock blind" diagnostic registered in advance as expectation-setting only, with results reported regardless of direction. Three facts emerged:

1. **Contamination fact.** 79 of the 99 published sites are in the release training panel (100% of the panel's tomato positives). They cannot inform expectations and are excluded from interpretation.
2. **Clean subset.** The 20 sites the model never saw rank moderately high: 13/20 (65%) above the candidate-table median (mean percentile 60.5; two-sided binomial vs 50% p = 0.26).
3. **The decisive observation.** **0 of 99** published sites appear in the candidate table's top 2,000.

The same pattern repeats on seven Arabidopsis control sites from laboratories independent of all training data (DES1, RBOHD, SnRK2.6/OST1, ABI4, ATG4a, AtG6PD6, PAD3; all coordinate-verified against the v2 proteome): 5/7 above the median (mean percentile 56.7, range 4.9–85.8; two-sided p = 0.45).

Our reading, stated before unblinding: the frozen model carries a genuine but **mild** mid-list ranking signal, and true sites do **not** concentrate at the list top. The blind Top-K enrichment we consider plausible a priori is modest (odds ratio on the order of 1.5–2.5), which directly motivates the large-K design of R5. This calibration is registered with full per-site data in the project archive (`results/known_controls/`, `results/data_archive/nc_entry_gates/negative_results/`); it changes nothing in the frozen artifacts — it only prevents the field (and us) from over-reading the blind outcome in either direction.

### R5. Preregistered blind-cohort design and statistical analysis plan

The blind validation is fixed in two co-signed documents — the Site-Aware Prioritization (SAP) protocol and the analysis plan (`sap_protocol.json`, `analysis_plan.json` in the release package; preregistration DOI **[TBD at submission]**). Core elements:

- **Cohort.** Top-K candidates from the frozen table plus matched controls, assayed blind by the wet-laboratory partner (Zhang Hua laboratory); labels remain with the wet laboratory until the preregistered unblinding step. Sites on proteins carrying published kiae271 persulfidation are excluded from the primary analysis and retained only as a preregistered sensitivity stratum; process controls are defined by the wet laboratory before measurement.
- **Primary endpoint.** Enrichment of experimentally confirmed persulfidation sites among Top-K candidates relative to matched controls; one-sided Fisher exact test (α = 0.05, preregistered direction) with odds ratio and exact 95% CI. **Success criteria:** p < 0.05 **and** OR > 1 **and** ≥ 2 confirmed candidate sites.
- **Design and power.** Because the effect size is unknowable before blind data exist, power was tabulated over the full grid K ∈ {50,100,200,300} × controls-per-candidate ∈ {1,2,5} × background confirmation rate ∈ {1%,2%,5%} × assumed OR ∈ {1.5,2,3,5} (10,000 simulations per cell, seed 20260813; complete 144-cell table in `power_table.json`). Structurally, at fixed assay budget a larger K with fewer controls dominates. The recommended design is **K = 200 with 2 controls per candidate (600 assayed sites)**: at 2% background rate, power 0.305 / 0.677 / 0.978 at OR 2 / 3 / 5 (Table 4). If throughput allows ≥ 900 sites, K = 300 (900 sites) raises OR = 3 power to 0.850. The minimum-viable 400-site design (K = 200, 1:1) is an honest coin flip at OR = 3 (power 0.524) and is flagged as such. K and the control ratio are fixed before unblinding as the largest assayable design; the full power table is reported, not filtered to favorable rows.
- **Unblinding procedure.** Two independent persons verify the bundle, fit-manifest, candidate-table, control-table, and inference-script hashes; the blind-cohort scoring output is written and hash-locked **before** labels are released; the preregistered analysis then runs verbatim. Forbidden after unblinding: any model/feature/Top-K/threshold change, adding samples to rescue a failing endpoint, or re-ordering the candidate table.
- **Reporting.** The primary analysis is reported regardless of outcome; all negative results and exclusions are preserved with provenance; any deviation requires a dated, hashed, co-signed amendment.

### R6. Blind-cohort primary endpoint **[TO BE COMPLETED AT GATE 3 — DO NOT EDIT ANYTHING ABOVE THIS LINE AFTER UNBLINDING]**

**[Insert verbatim: the 2×2 table (Table 5), Fisher exact p, OR with exact 95% CI, per-endpoint numbers, process-control pass rate. Pre-written abstract outcome variants:**

- *Variant A (success):* "The preregistered blind cohort met its success criteria (odds ratio X.XX, exact 95% CI […], one-sided p = …, N confirmed candidates of K = …), demonstrating enrichment of confirmed persulfidation sites among frozen top-ranked candidates under fully blinded conditions."
- *Variant B (no success):* "The preregistered blind cohort did not meet its success criteria (odds ratio X.XX, exact 95% CI […], one-sided p = …); we report the complete result, the pre-registered sensitivity analyses, and the resource artifacts, which remain usable for candidate organisation under the stated limitations."

**]**

### R7. Sensitivity and secondary analyses **[TO BE COMPLETED AT GATE 3]**

**[Preregistered list, to be filled verbatim: (i) exclusion of any site on a training-panel protein; (ii) exclusion of kiae271-protein sites; (iii) leave-one-out over candidate proteins and over homology clusters; (iv) rank-percentile calibration of confirmed sites (KS / percentile bootstrap); (v) MC-dropout uncertainty stratification; (vi) per-structure-stratum analysis — noting that no tomato structure stratum exists in this release (R3), so (vi) is vacuous here and is retained for release-lineage completeness.]**

### R8. Mechanism follow-up **[OPTIONAL — Gate 4; only if pursued]**

**[If ≥1 confirmed blind hit is carried to a mechanism chain (site confirmation → Cys→Ser/Ala substitution → functional effect → plant phenotype → rescue), report here. Absence of this section does not affect the primary claim; its presence upgrades the journal path per the preregistered decision matrix.]**

---

## Discussion

This paper does three things that, together, we believe the cysteine-PTM field currently lacks: it freezes everything before the experiment, it calibrates expectations honestly before unblinding, and it commits to publishing the blind outcome either way.

**The within-dataset signal is real but its meaning is bounded.** Under literature-comparable splits the gated-fusion PU ranker leads five alternatives by a wide margin (macro AP 0.386 vs ≤ 0.047), and its per-species pattern (Arabidopsis 0.49, rice 0.65, tomato 0.019) is itself informative: performance tracks the availability of training positives and structure features per species, not any intrinsic "learnability" of the modification. Tomato combines the fewest positives (79, all from one paper) with zero structure coverage — and yields near-baseline within-dataset AP. The model nonetheless ranks previously unseen confirmed tomato sites moderately above the list median (R4), consistent with a weak but real sequence-context signal. We deliberately make no claim beyond these within-dataset statements; the mandated limitation quoted in the Introduction applies to every number in this paper.

**Calibration before unblinding is the contribution other fields call registered reporting.** The mock-blind diagnostic told us, before any blind data existed, that the enrichment to expect is modest (prior OR ~ 1.5–2.5) and that success criteria built on list-top concentration would be self-deceiving. The consequence was not to tune the model — that path is closed by the freeze — but to size the experiment for a modest effect (large K, few controls per candidate; R5) and to pre-commit to success criteria that a modest effect can actually meet or honestly fail. We suggest this two-step (calibrate on published sites, then size the blind design) is reusable for any candidate-ranking tool in low-data proteomics.

**Limitations.** (i) All training positives trace to four studies, two of them same-laboratory Arabidopsis; the dataset cannot support cross-study statements, and none are made. (ii) Tomato structure coverage is 0%; the structure branch, which contributes within-dataset value for Arabidopsis and rice, is inert for tomato — the ranking rests on three sequence features. (iii) The flat score distribution (top-200 within 0.01 score units) means rank order at the very top is fragile to small perturbations; we mitigate by reporting enrichment over K = 200–300 rather than any small "top-N". (iv) MC-dropout uncertainty is deterministic per device but not bit-identical across CPU and GPU (different RNG streams); point scores are bit-identical. (v) The matched-control design controls sequence-feature strata, not expression, localization, or assay detectability; residual confounding by detectability cannot be excluded and is why process controls are part of the primary report. (vi) A blind failure would not falsify persulfidation biology — it would bound the usefulness of *this* feature set for *this* species, which is itself the published answer this design guarantees.

**What happens next is already written down.** Per the preregistered decision matrix: blind success plus a complete mechanism chain (R8) upgrades the submission target; blind success without mechanism stands as a Plant Physiology-grade resource and validation report; blind failure with intact resources stands as the same, with the negative result as a first-class finding. No model development is permitted to "rescue" a failed endpoint; any improved model (e.g. with tomato structure coverage, recalibrated scores, or expanded training positives — all registered in the project's deferred-improvement register) requires a new release number and a new preregistration.

---

## Materials and Methods

### Dataset assembly and provenance

Site-level positives and unlabeled cysteine pools were integrated from PXD006140, PXD024061 (Arabidopsis), PXD072089 (rice), PXD063170 (*M. oryzae*), and kiae271 (tomato) into the multispecies v2 dataset with per-row study provenance; coordinates verified against SHA-pinned reference proteomes (tomato reference proteome SHA256 `1a46efe6…`). The release training panel (389,609 rows, 2,334 positives; SHA256 `30e432fb…`) is the union of the ten literature-random-track seed panels, development split only. Exclusion ledger (8 entries) and full input inventory (243 files with SHA256) are in the data archive.

### Features

Per-site input features: three sequence-window features (local hydrophobicity, cysteine density, local positive charge density; exact definitions frozen in `feature_schema.json`) and two AlphaFold-derived structure features (contact proxy, pLDDT) with an explicit availability mask — missing structures are masked, never mean-imputed. ESM-2 embeddings and study-context features are implemented but disabled in this release (`use_esm = 0`, `use_study_context = false`).

### Model and training

`structure_ranker` is a gated-fusion positive-unlabeled ranker (branch encoders with gated combination; PU objective). Frozen hyperparameters: hidden width 16, dropout 0.2, 200 epochs, learning rate 0.05, holdout fraction 0.2, 16 MC-dropout passes for uncertainty. Trained once on the full panel with seed 20260813 (freeze date; not performance-selected) on CUDA; serialized as `StructureRankerBundle` (PyTorch, pickle protocol 4, atomic write; SHA256 `ab8a0353…`). Point scores are bit-identical across CPU and GPU; MC-dropout uncertainty is deterministic per device. Code revision `d0c740c…` [full hash in fit_manifest.json].

### Benchmark protocols

Literature-comparable track: random protein-level holdout, 10 seeds, 6 models; metric = macro AP across the three crop species (scripts and manifests in `results/experiments/multispecies_v2_global_clusters_v11/`). Strict track: 5-fold homology-cluster (MMseqs2, 30% identity) development splits; 20% frozen test partition locked **[Gate 1 unlock procedure: scripts/train_multispecies_v2.py --score-test --test-unlock, under manifest audit]**.

### Candidate table and matched controls

`scripts/generate_tomato_topk_candidates.py`: verifies the config against the frozen manifest, rebuilds and hash-verifies the training panel, scans every cysteine of the tomato reference proteome, drops panel sites, extracts features, scores with the frozen bundle, and writes the ranked table (179,736 rows; SHA256 `36963c71…`) plus matched controls (top-200 head; ≤5 controls per candidate; same structure-availability stratum; z-scored sequence-feature Euclidean distance; caliper 0.25 SD; no reuse; 966 pairs; SHA256 `c1752b07…`). Arabidopsis and rice diagnostic scans used the identical script with species-specific proteomes (outputs under `results/diagnostics/candidates_<species>/`).

### Pre-blind calibration

kiae271 sites (99 coordinate-verified) and 7 registered Arabidopsis independent-laboratory controls were scored with the frozen bundle and ranked by percentile within the respective frozen candidate tables (scripts `evaluate_kiae271_with_release_bundle.py`, `evaluate_arabidopsis_controls_with_release_bundle.py`); two-sided binomial tests against the 50th percentile; panel overlap computed by exact key match. Per-site results: `results/known_controls/kiae271_release_bundle_recovery_v1.rows.tsv`, `arabidopsis_release_bundle_recovery_v1.rows.tsv`.

### Blind-cohort statistics and power analysis

Primary test: one-sided Fisher exact test on the 2×2 (candidate/control × confirmed/not) table; odds ratio with exact 95% CI; α = 0.05. Power: Monte-Carlo simulation (10,000 replicates per cell, seed 20260813) over the 144-cell grid described in R5 (`scripts/compute_release_power_table.py`; `power_table.json`). Secondary: KS / percentile-bootstrap rank calibration; MC-dropout uncertainty stratification (exploratory); leave-one-out sensitivity at protein and homology-cluster level.

### Reproducibility and wording governance

All frozen artifacts carry SHA256 sums (`SHA256SUMS`, `manifest.json`); `scripts/audit_candidate_release.py` re-verifies hashes, bundle determinism, document signatures, and wording-gate compliance in one pass. Environment lock (`pip freeze`) and reproduction command index are in the data archive. All outward text of this project, this manuscript included, must pass the project's automated wording gate (`verify_predictive_claims` → no matches) before release.

---

## Tables

**Table 1. Proteome-wide candidate scans with the frozen model.**
| species | role | panel sites excluded | candidate sites | structure-covered | max score |
|---|---|---|---|---|---|
| tomato (*S. lycopersicum*) | frozen release | 72,051 | 179,736 | 0 (0%) | 0.3053 |
| Arabidopsis | diagnostic | 111,496 | 284,387 | 3,050 (1.07%) | 0.3615 |
| rice (*O. sativa*) | diagnostic | 185,379 | 469,748 | 3,539 (0.75%) | 0.3614 |

**Table 2. Within-dataset literature-comparable benchmark (macro AP, 10 seeds).**
| model | mean | std | min | max |
|---|---|---|---|---|
| structure_ranker | 0.386 | 0.173 | 0.166 | 0.554 |
| esm_linear_head | 0.047 | 0.007 | 0.037 | 0.055 |
| xgboost | 0.036 | 0.008 | 0.030 | 0.059 |
| pu_logistic | 0.035 | 0.007 | 0.028 | 0.050 |
| random_forest | 0.030 | 0.004 | 0.026 | 0.038 |
| sul_bertgru (retrained) | 0.021 | 0.002 | 0.018 | 0.024 |

*Claim class for this table (verbatim, binding): "Literature-comparable within-dataset performance improvement." These numbers are within-dataset evidence only.*

**Table 3. Frozen artifact registry (release `multispecies-v2-candidate-release-v1`, 2026-08-13).**
| artifact | SHA256 |
|---|---|
| model bundle (`structure_ranker_bundle.pt`) | `ab8a03532898726b765d3da18fe34736e7920cfed82134ec06195dc21b156cd3` |
| fit manifest | `caf3c7f4d82b7dc88836bc2c16a48f5209ef95ac600bc0fbfc7a8722d07f8dbb` |
| training panel (389,609 rows) | `30e432fb3c558c69e92362fb36d26abe1dee0f25bd23df417ad40886449156cd` |
| tomato candidate table (179,736 rows) | `36963c7144792752a6d6b6e197da6c24fc0c8844d9540d8b63cf861dc7b7e949` |
| matched controls (966 pairs) | `c1752b071e4413b17e41e4206d013c3c32759a6d7825da6cc299520abaefb404` |
| frozen experiment config (v11) | `ed8b76bef0b201bb669322ff5610e84c42ba422aa8cf5b9eb274427325f7dbaf` |
| tomato reference proteome | `1a46efe6461b239ec80b915a4c23ea39bbe39e37d43f9cc7c38b648523efca0f` |
| training seed | 20260813 |
| code revision | `d0c740cc24a0ce81157d23921f7875ea5f862606` |

**Table 4. Preregistered design options and power (p0 = 2% background confirmation rate; one-sided Fisher, α = 0.05; 10,000 simulations, seed 20260813).**
| design | K | controls/candidate | total sites | power OR = 2 | power OR = 3 | power OR = 5 |
|---|---|---|---|---|---|---|
| recommended | 200 | 2 | 600 | 0.305 | 0.677 | 0.978 |
| escalation (if throughput allows) | 300 | 2 | 900 | — | 0.850 | — |
| minimum viable | 200 | 1 | 400 | — | 0.524 | — |

*Full 144-cell grid: `power_table.json` (release package). The full table is reported unfiltered per the analysis plan.*

**Table 5. [GATE 3 PLACEHOLDER] Primary 2×2 table: candidates vs matched controls × confirmed/not-confirmed, with Fisher exact p and OR (exact 95% CI).**

---

## Figure legends

- **Figure 1. Study design and gate pipeline.** [TODO: schematic — dataset → benchmark → freeze (hashes) → calibration → preregistration → blind cohort → decision matrix.]
- **Figure 2. Within-dataset benchmark.** Per-seed macro AP for six models (10 seeds; data from Table 2). [Script: TODO `scripts/figures/fig2_benchmark.py` from literature manifest.]
- **Figure 3. Frozen tomato candidate release.** Score distribution of 179,736 candidates with top-200 head inset (flat-top detail); matched-control distance distribution. [Data: `top_k_candidates.tsv`, `matched_controls.tsv`.]
- **Figure 4. Pre-blind calibration.** Percentile ranks of 20 never-trained kiae271 sites and 7 Arabidopsis independent-laboratory controls within their frozen candidate tables; panel-overlap annotation (79/99). [Data: `results/known_controls/*.rows.tsv`.]
- **Figure 5. Power landscape.** Heatmaps of simulated power over K × assumed OR at p0 = 2% for control ratios 1/2/5; recommended design marked. [Data: `power_table.json`.]
- **Figure 6. [GATE 3 PLACEHOLDER] Blind-cohort result:** confirmation rate in candidates vs controls with exact CIs; rank-percentile distribution of confirmed sites.

## Supplemental data

- S1: per-species dataset composition and provenance [TODO: from dataset manifest].
- S2: full power table (144 cells) — ships as `power_table.json`.
- S3: per-site calibration results (kiae271, Arabidopsis controls).
- S4: exclusion ledger and negative-results register (verbatim from data archive).
- S5: SAP protocol and analysis plan (co-signed versions) + preregistration record [DOI TBD].
- S6: feature schema and inference script (frozen copies).

## Accession numbers

Proteomics data: PXD006140, PXD024061, PXD072089, PXD063170 (ProteomeXchange); kiae271 (doi:10.1093/plphys/kiae271). Reference proteomes: [UniProt proteome IDs TBD]. All other accessions in Supplemental Table S1.

## Data and code availability

Frozen release package: `results/candidates/multispecies_v2_candidate_release_v1/` (see Table 3 hashes). Data archive (inventories, exclusion ledger, negative-results register, environment lock): `results/data_archive/nc_entry_gates/`. Code: [repository URL TBD] at revision `d0c740c…`; permanent DOI **[TBD at submission]**. Preregistration DOI **[TBD]**.

## Author contributions

[CRediT TBD — L.G.: conceptualization, methodology, software, formal analysis, writing — original draft; Zhang Hua lab: investigation (blind assays), resources, writing — review; ...]

## Funding

[TBD]

## Conflict of interest

The authors declare no conflict of interest. [VERIFY]

## Acknowledgements

[TBD — including the laboratories whose public datasets made this work possible.]

---

## References [skeleton — VERIFY each entry before submission]

1. Zhang et al. (2024) Tomato persulfidome and fruit ripening. *Plant Physiology* — doi:10.1093/plphys/kiae271.
2. Aroca Á, Gotor C, Romero LC, et al. — Arabidopsis persulfidome (PXD006140) [VERIFY full citation].
3. Jurado-Flores A, et al. — Arabidopsis persulfidation study (PXD024061) [VERIFY].
4. Xie et al. (2026) Rice persulfidome (PXD072089). *PNAS* [VERIFY].
5. *Magnaporthe oryzae* persulfidation dataset (PXD063170) [VERIFY authors/journal].
6. Zhang et al. (2026) PAD3 persulfidation. *Plant, Cell & Environment* — doi:10.1111/pce.70593.
7. Sul-BertGRU [VERIFY citation]; companion data audit: this paper, Methods.
8. iCysMod database [VERIFY]; pCysMod [VERIFY].
9. Jumper J et al. (2021) Highly accurate protein structure prediction with AlphaFold. *Nature* 596:583–589.
10. Steinegger M, Söding J (2017) MMseqs2. *Nature Biotechnology* 35:1026–1028.
11. Lin Z et al. (2023) ESM-2. *Science* 379:1123–1130. [VERIFY]
12. Elkan C, Noto K (2008) Learning classifiers from only positive and unlabeled data. *KDD* 213–220. [VERIFY relevance]
13. Nosek BA et al. (2018) The preregistration revolution. *PNAS* 115:2600–2606. [VERIFY]
14. Chambers CD (2013) Registered reports. *Cortex* 49:609–610. [VERIFY]
