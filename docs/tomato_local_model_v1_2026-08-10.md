# Tomato-local model v1: within-species training on kiae271

> Date: 2026-08-10 | Branch: phase-a-persulfidation-site-evidence
> Predecessor: `docs/gate2_condition3_kiae271_expansion_2026-08-10.md` (cross-species negative result)
> Config: `configs/experiments/tomato_local_v1.yaml`

## 1. Motivation

The kiae271 panel (99 coordinate-verified tomato Cys persulfidation sites, Zhang et al. 2024, doi:10.1093/plphys/kiae271) was previously used as an **external** cross-species control panel. An Arabidopsis-trained 2-feature pu_logistic model scored those 99 sites: 48/99 recovered (48.5%, binomial p=0.84 vs 50% chance) -- statistically indistinguishable from random.

The question now: if we train AND evaluate within tomato, using the same 99 kiae271 sites as target-species supervision, can we build a useful tomato-specific candidate-prioritisation model? This is a different, legitimate scientific question from "does an Arabidopsis-trained model transfer to tomato" (answered: no).

## 2. Design decisions (frozen before any result was seen)

### Cross-validation

**5-fold homology-cluster-grouped.** The full tomato reference proteome (36,988 proteins) was clustered with MMseqs2 easy-cluster (30% identity, >=50% coverage, cov-mode 0) into 17,758 homology clusters. Positive-bearing clusters are dealt round-robin into 5 folds first; background-only clusters follow. Every member of a homology cluster lands in exactly one fold -- no protein family crosses a fold boundary. Seed: 20260810.

The 88 kiae271-detected proteins belong to 85 distinct homology clusters (3 pairs of homologs share clusters). With 5 folds, this gives ~17 positives per fold.

### Two evaluation arenas

| Arena | Background | Controls for |
|-------|-----------|-------------|
| `proteome` | All other Cys in the tomato proteome (subsampled 1:20, seed=12345) | MS protein-level detection bias |
| `panel` | Only other Cys of the 88 kiae271-detected proteins (subsampled 1:10, seed=12346) | Protein-abundance confound |

The `panel` arena holds protein identity constant -- all rows are from the same 88 proteins, so the model must learn which Cys within already-detected proteins carries persulfidation, not which proteins are abundant enough to appear in the screen.

### Features

- `seq_2`: [hydrophobicity, cys_density] -- replicates the 2-feature cross-species baseline exactly
- `seq_3`: [hydrophobicity, cys_density, local_positive_charge_density] -- adds the thiolate-stabilisation proxy (K/R fraction in flanking window; mechanistically motivated by Corpas et al. COPLBI-D-26-00068)
- `seq_structure`: seq_3 + [contact_number_proxy, plddt] from AlphaFold DB predicted structures

### Models and evaluation

pu_logistic (Elkan-Noto) across all three arms. 5 model seeds (0-4); ensemble = mean score per site across seeds. 5-fold grouped cross-validation; per-fold AP plus 1000-iteration cluster-level bootstrap CI.

### The 20 unverified Dataset S1 rows

The 20 kiae271 supplementary rows that failed coordinate verification, position alignment, or the localization-probability bar are **never recycled as background**. The `kiae271_excluded_keys` function enumerates every (accession, Cys position) the supplementary table touches, and the builder removes all of them that are not verified positives from the unlabeled pool. For rows where position alignment is ambiguous (parallel `Proteins` and `Positions` lists of different length), the cartesian product is excluded -- deliberately over-inclusive, since dropping a handful of Cys from a 252,000-cysteine background is harmless while recycling one as background would corrupt the negative class.

### Ambiguous-study input rows

Rows where the two conditions disagree (a peptide was detected in both LCD1-OE and WT, but with different MS intensities -- a quantitative "up" or "down") are included as positives irrespective of the direction. The task is *persulfidation-site detection* (any condition), not differential regulation between conditions. Each site's `regulation` label (lcd_gain / wt_only / both) is recorded but not used as a training target.

## 3. AlphaFold structure coverage

Panel proteins: 82/88 (93.2%) have AlphaFold DB structures (all v4), downloaded via the existing `fetch_alphafold_structure` pipeline and registered in `data/registry/alphafold_structures.tsv`. Six accessions genuinely absent from AlphaFold DB:

`A0A3Q7G7S9, A0A3Q7HD77, A0A3Q7HRP0, A0A3Q7I6P2, A0A3Q7IBY6, A0A3Q7IQY4`

These are true "no predicted structure" facts (HTTP 404 from both v4 and v6 URLs), not download failures.

Background proteins: in the `proteome` arena, 2 out of 1,880 unique background proteins have AlphaFold structures (both happen to be other UniProt accessions of Arabidopsis orthologs also registered for the Arabidopsis benchmark). The remaining 1,878 background proteins have no structure. **This 82:2 asymmetry between positive and background structure coverage is the critical confound flagged in Section 4.3 below.**

## 4. Results

### 4.1 Proteome arena (background = all-tomato-proteome Cys, base_rate = 4.8%)

| Arm | Mean fold AP | Ensemble AP | Bootstrap 95% CI | vs base_rate |
|-----|-------------|-------------|-------------------|-------------|
| seq_2 | 0.050 | 0.046 | [0.036, 0.062] | 1.0x (random) |
| seq_3 | 0.115 | 0.095 | [0.071, 0.134] | **2.4x** |
| seq_structure | 0.943 | 0.885 | [0.809, 0.954] | see Section 4.3 |

Per-fold AP (seq_3): [0.155, 0.079, 0.115, 0.110, 0.118]

The 2-feature baseline (seq_2) is indistinguishable from the base rate in every fold. Adding `local_positive_charge_density` roughly doubles the AP in every fold -- the improvement is consistent and statistically significant (bootstrap 95% CI of seq_3 AP [0.071, 0.134] excludes the base rate of 0.048).

### 4.2 Panel arena (background = only other Cys of the same 88 proteins, base_rate = 9.1%)

| Arm | Mean fold AP | Ensemble AP | Bootstrap 95% CI | vs base_rate |
|-----|-------------|-------------|-------------------|-------------|
| seq_2 | 0.110 | 0.094 | [0.076, 0.119] | 1.0x (random) |
| seq_3 | 0.203 | 0.163 | [0.128, 0.215] | **2.2x** |
| seq_structure | 0.192 | 0.156 | [0.121, 0.216] | 1.7x |

Per-fold AP (seq_3): [0.239, 0.131, 0.248, 0.195, 0.201]

The panel arena is the harder question: holding protein identity constant, can the model distinguish which Cys within the same protein is persulfidated? The answer: yes, but modestly. seq_2 is again indistinguishable from the base rate; seq_3 achieves ~2.2x enrichment, with the bootstrap CI excluding the base rate. Crucially, **structure features add nothing in this arena** -- seq_structure AP (0.192) is slightly *below* seq_3 AP (0.203). Within the same protein, pLDDT is identical for every Cys, and the contact-number proxy does not provide per-residue discrimination above what the flanking-window charge density already provides.

### 4.3 Structure leakage diagnosis: proteome arena seq_structure

The proteome arena seq_structure AP of 0.943 is not a genuine persulfidation signal. The structure features act as a near-perfect proxy for "which protein is detected in the kiae271 screen":

- 82/88 (93%) positive-bearing proteins have AlphaFold structures
- 2/1880 (0.1%) background proteins have AlphaFold structures

The model learns in effect: "if structure features are non-zero, rank high." This is memorisation of the protein-detection status, not learning of persulfidation biochemistry. The fact that `seq_structure` adds nothing in the panel arena (where protein identity is held constant) confirms this interpretation -- when every row comes from a protein with a structure, the structure features provide no additional discrimination.

**This is an honest result about coverage limits, not an error.** It demonstrates why structure-based modelling requires near-complete structure coverage of the evaluation background to be meaningful. For this specific dataset, the coverage gap is order-of-magnitude: 93% for positives vs 0.1% for background. A future iteration that downloads AlphaFold structures for all 1,880 background proteins would close this gap, but at present the structure arm cannot be evaluated fairly in the proteome arena.

## 5. What `local_positive_charge_density` specifically contributes

The ablation is clean and consistent across both arenas:

| Arena | seq_2 AP | seq_3 AP | Charge density gain |
|-------|---------|---------|-------------------|
| proteome | 0.050 | 0.115 | +0.065 (+130%) |
| panel | 0.110 | 0.203 | +0.093 (+85%) |

In both arenas, seq_2 is indistinguishable from the base rate. **The entire signal above random comes from `local_positive_charge_density`.** Without it, the model is random; with it, the model achieves ~2x enrichment across both arenas.

This is mechanistically plausible. Persulfidation requires nucleophilic attack by the anionic thiolate (deprotonated Cys-S^-), and nearby positive charge (Lys/Arg) lowers the local thiol pKa, favouring the deprotonated/reactiive form. The feature is computed as the fraction of K/R residues in the flanking window, and it is fully sequence-based -- it works for tomato without any AlphaFold structure requirement.

## 6. Comparison to the cross-species negative result

| Experiment | Training data | Evaluation | Result |
|-----------|--------------|-----------|--------|
| Cross-species (kiae271 panel) | Arabidopsis (PXD006140+PXD024061, 390 positives) | 99 tomato kiae271 sites, held out | 48/99 recovered (48.5%, p=0.84) -- random |
| Within-tomato seq_2 | Tomato kiae271 (99 positives) | 5-fold CV, homology-cluster-grouped | AP 0.050, indistinguishable from base_rate 0.048 |
| Within-tomato seq_3 | Tomato kiae271 (99 positives) + charge_density | 5-fold CV, homology-cluster-grouped | AP 0.115, 2.4x above base_rate |

The cross-species failure was not because persulfidation is inherently unpredictable. It was because a model trained on Arabidopsis-persulfidation-correlated features does not transfer to tomato. Within tomato, with a mechanistically motivated tomato-local feature (charge density), there IS a real signal -- modest (2.4x enrichment at best) but statistically significant and consistent across folds.

## 7. Limitations (honest)

1. **Single-study, single-lab, single-chemistry, single-species.** The 99 positives come from exactly one experiment (Zhang lab, 4D label-free + iodoTMT, LCD1-OE vs WT tomato leaf). Cross-validation is method-level only; it does not demonstrate chemical-, laboratory-, or species-independent predictive value.
2. **Small positive set.** 99 positives over 88 proteins is a very small training set by any standard. The between-fold variance is substantial (AP ranges from 0.08 to 0.24 across folds in the panel arena). Individual-fold bootstrap CIs would be wide.
3. **Structure confound in the proteome arena.** The seq_structure arm is not meaningfully evaluable without downloading AlphaFold structures for all background proteins.
4. **The model is a prioritisation tool, not a predictor.** With a base_rate of 5-9% and an AP of 11-20%, the model enriches persulfidation candidates 2-2.4x above random guessing -- useful for organising a ranked candidate list, but not strong enough to be treated as a reliable binary classifier.
5. **No held-out external tomato validation set exists** beyond the tiny CAT1 control (P30264 Cys234, a different lab's single site, present in the proteome arena background as an unlabeled row -- not used as a positive).

## 8. Products

- `src/plantpersulf/proteomics/tomato_local_dataset.py` -- PU dataset assembly, excluded-keys enumeration, grouped folds
- `tests/scientific/test_tomato_local_dataset.py` -- 11 unit tests (11/11 GREEN)
- `scripts/run_tomato_local_v1.py` -- main experiment script
- `configs/experiments/tomato_local_v1.yaml` -- frozen experiment config
- `data/processed/clusters/tomato_proteome_clusters_v1.tsv` -- MMseqs2 clusters for the tomato proteome (36,988 proteins, 17,758 clusters, sha256: b129e6c9...)
- `results/tomato_local_v1/summary.json` -- full results
- `results/tomato_local_v1/scores_{arena}_{arm}.tsv` -- per-site ensembled scores (4 files)
- `data/registry/alphafold_structures.tsv` -- 82 new tomato panel entries appended (total now 2,088)

## 9. Next steps (uncommitted, for future decision)

1. **Download AlphaFold structures for the full proteome-arena background** (~1,880 proteins), then re-run the seq_structure arm so the coverage asymmetry is eliminated and the arm can be evaluated fairly. This would close the one honest gap in the current analysis.
2. **Explore per-residue additional features**: B-factors beyond pLDDT, SASA via FreeSASA (more rigorous than the contact-number proxy), or secondary-structure type at the Cys position -- any of these might add value in the panel arena where protein-level features are held constant.
3. **Stratify by regulation**: the 99 positives include 43 `lcd_gain` (H2S-elevated), 40 `wt_only`, and 16 `both`. A model trained specifically on `lcd_gain` might capture a different signal.
4. **CAT1 cross-lab holdout**: P30264 Cys234 (doi:10.1016/j.plaphy.2020.09.020) is an independent tomato persulfidation site from a third lab. It is present as an unlabeled row in the proteome arena background (not excluded); checking its rank under the seq_3 model would provide a tiny, fully-external sanity check.
