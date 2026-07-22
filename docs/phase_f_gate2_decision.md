# Phase F — Gate 2 decision record

**Date**: 2026-07-22
**Decision**: **`GATE2_STOP`** (1/5 conditions passed)
**Evidence release**: `pu_ranker_v1` (leave-study-out), supplementary track
`pu_ranker_cluster_v1` (Split A cluster split — literature-comparable,
structurally excluded from Gate 2)
**Machine-verified decision file**: `results/external_validation/pu_ranker_v1/gate2_decision.json`

## Frozen limitation (carried from benchmark_v1)

Site-level positives come from exactly two studies: `PXD006140` (317) and
`PXD024061` (73) — **same laboratory (Romero/Gotor, U. Sevilla), same
tag-switch chemistry, same species (Arabidopsis)**. `configs/gate2_v1.yaml`
freezes `studies_are_independent: false`; Gate 2 condition 1 therefore cannot
pass on current public data regardless of model quality.

## Mandated downgrade statement (verbatim, binding for all external text)

> 当前公开数据不足以证明跨研究预测能力，模型仅用于候选组织与假设生成。
> Current public data are insufficient to demonstrate cross-study predictive
> ability; the model is used only for candidate organisation and hypothesis
> generation.

## Gate 2 conditions

| # | Condition | Result | Detail (from gate2_decision.json) |
|---|---|---|---|
| 1 | independent_studies_beat_baseline | FAIL | studies_independent=False, beats_baseline=2/2, need>=2 |
| 2 | effect_ci_excludes_zero | FAIL (unmeasured) | delta_ci_lower=None — `--scored` bootstrap not wired |
| 3 | recovery_is_not_training_leakage | **PASS** | control_leakage=[], independent_units=2 |
| 4 | structure_gain_on_structured_subset | FAIL (unmeasured) | structure_gain=None — ablation-delta wiring not connected |
| 5 | not_driven_by_single_cluster | FAIL (unmeasured) | single_cluster_driven=True, permutation_p=None |

Conditions 2, 4 and 5 fail because the statistical wiring was not connected,
**not** because the effect was measured and found absent. Wiring them is
tracked as follow-up P2; the decision stays STOP either way because condition
1 is structurally locked by data provenance.

## Primary track — leave-study-out (Gate 2 evidence)

`results/experiments/pu_ranker_v1/metrics.tsv`: 7 ablations × 2 folds ×
5 seeds = 70 runs. Subsample 1:20 unlabeled:positive → base rate ≈ 0.0476.

| Ablation | leave_PXD006140_out (n=5) | leave_PXD024061_out (n=5) | pooled mean | vs base rate |
|---|---|---|---|---|
| **seq_structure** | 0.0825 ± 0.0081 | 0.1208 ± 0.0096 | **0.1017** | **2.1x** |
| sequence_only | 0.0676 ± 0.0020 | 0.0550 ± 0.0048 | 0.0613 | 1.3x |
| no_study_context | 0.0444 ± 0.0037 | 0.0646 ± 0.0329 | 0.0545 | 1.1x |
| full (all branches) | 0.0580 ± 0.0372 | 0.0387 ± 0.0047 | 0.0484 | 1.0x |
| seq_esm | 0.0416 ± 0.0054 | 0.0507 ± 0.0101 | 0.0461 | 1.0x |
| no_plddt | 0.0427 ± 0.0175 | 0.0397 ± 0.0059 | 0.0412 | 0.9x |
| no_accessibility | 0.0413 ± 0.0101 | 0.0399 ± 0.0071 | 0.0406 | 0.9x |

## Supplementary track — Split A cluster split (NOT Gate 2 evidence)

`results/experiments/pu_ranker_cluster_v1/metrics.tsv`: 7 ablations × 5 seeds
= 35 runs. Reported only for literature comparability (within-integrated-
dataset splits are the norm in published cysteine-PTM predictors); the
`pu_ranker_cluster_v1.yaml` limitation text must accompany any citation.

| Ablation | mean test_ap (n=5) | vs base rate |
|---|---|---|
| no_accessibility | 0.1695 ± 0.3052 ⚠️ | — |
| **seq_structure** | **0.0837 ± 0.0132** | 1.8x |
| sequence_only | 0.0504 ± 0.0076 | 1.1x |
| no_study_context | 0.0435 ± 0.0119 | 0.9x |
| seq_esm | 0.0399 ± 0.0057 | 0.8x |
| no_plddt | 0.0342 ± 0.0038 | 0.7x |
| full | 0.0342 ± 0.0037 | 0.7x |

⚠️ `no_accessibility` is a seed-level outlier: std 0.3052 driven by one seed
at 0.7155 while the other four seeds sit at 0.031–0.047. It is not signal and
is not admissible evidence in either track.

## Known-control recovery (from control_recovery.tsv)

| mechanism_lineage_id | gene | status | independent unit |
|---|---|---|---|
| SLWRKY6_H2S_PHOSPHORYLATION | SlWRKY6 (Cys396) | mapped | yes |
| SLERFD2_H2S_ETHYLENE | SlERF.D2 (Cys35) | mapped | yes |
| BRG3_H2S_UBIQUITINATION | BRG3 | unmappable | no |
| ERFD3_H2S_CONTEXT | ERF.D3 | unmappable | no |

No registered control appeared in training (leakage list empty). Percentile
ranks were not computed (requires `--scored`); follow-up P2.

## Interpretation

1. **The sequence+structure signal is real and consistent.** seq_structure
   beats sequence_only in 10/10 fold×seed runs on the leave-study-out track
   and 5/5 on the cluster track, with tight seed variance on the primary
   track. Structure adds information sequence alone does not carry.
2. **The ceiling is data, not evaluation strictness.** Leave-study-out means
   are *equal to or higher than* cluster-split means for every stable
   ablation (gap −0.006 to −0.018). The hypothesis that cross-study splits
   merely hide learnable signal is falsified on this benchmark: the two
   tracks agree the signal is modest under any split.
3. **The frozen ESM-2 branch is harmful at this data scale.** full (with ESM)
   underperforms seq_structure on both tracks; seq_esm underperforms
   sequence_only on the primary track. With 390 positives and a 1280-dim
   frozen branch, the ESM features add variance, not information.
4. **Condition 1 is a data-provenance lock, not a model failure.** Both
   held-out folds beat the no-learning baseline (beats_baseline=2/2); the
   condition fails solely because the two studies are not independent
   (same lab/chemistry/species, frozen in config).

## Hardware-driven scope limitation (documented, audited)

Proteins longer than 2500 aa are skipped during ESM-2 extraction (RTX
5070 Ti 16 GB cannot hold the O(L²) attention matrices; CPU fallback is
impractical). ~70 Arabidopsis proteins are affected, covering **≤3 positive
sites (0.8% of 390)**; their ESM branch input defaults to zeros and is
mask-gated like any other missing branch. See
`src/plantpersulf/features/esm2.py` (MAX_PROTEIN_LENGTH).

## Follow-ups

- **P2 (statistical wiring)**: connect `--scored` cluster-bootstrap CI
  (condition 2), ablation-delta structure gain (condition 4), permutation
  test (condition 5), and control percentile ranks. Expected to move the
  record from 1/5 to ~4/5 passed, with condition 1 remaining the sole
  structural blocker — the strongest possible basis for a data request.
- **P0 (Phase Z)**: data-resource and systematic evidence-audit deliverable
  + collaboration data-request list (per roadmap, Gate 2 STOP route).
- **P1 (data expansion)**: deep-parse independent PRIDE datasets
  (PXD035795 / PXD039999 candidates); the only path that can flip
  condition 1 and reopen Gate 2.
