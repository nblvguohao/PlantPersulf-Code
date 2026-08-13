# Phase F — Gate 2 decision record

**Date**: 2026-07-22 (updated 2026-07-23: PXD072089 expanded to 929 sites
via SS-all-peptides.tsv; UniParc recovery checked and found not to add
positives; PXD072300 same-paper recombinant controls registered; updated
2026-07-24: real MMseqs2 clustering — see "Cluster-granularity remediation"
below, conditions 2/5 re-verified on real homology clusters, decision
unchanged)
**Decision**: **`GATE2_STOP`** (3/5 conditions passed)
**Evidence release**: `pu_ranker_v1` (leave-study-out), supplementary track
`pu_ranker_cluster_v1` (Split A cluster split — literature-comparable,
structurally excluded from Gate 2)
**Cross-species tracks**: PXD063170 (Magnaporthe) + PXD072089 (rice, v2 =
929 sites) — both weak ∼1.3x base rate, both seed-unstable; wording locked
to "transfer only, not Gate 2 evidence" (see `docs/phase_z_evidence_audit.md`
§5.2–5.4)
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

## Gate 2 conditions (after P2 statistical wiring)

Release arm: `structure_ranker:seq_structure`; baseline arm: `pu_logistic`;
structure-ablated arm: `structure_ranker:sequence_only`. Per-site scores
regenerated deterministically on CPU by `scripts/score_release.py` under the
frozen config (same partition/subsample/split as the release run); paired
statistics in `src/plantpersulf/evaluation/effect_size.py`.

| # | Condition | Result | Detail (from gate2_decision.json / external_validation.json) |
|---|---|---|---|
| 1 | independent_studies_beat_baseline | FAIL (structural) | studies_independent=False, beats_baseline=2/2, need>=2 |
| 2 | effect_ci_excludes_zero | FAIL (**measured**) | delta_ci_lower=−0.0032 — paired cluster-bootstrap delta vs `pu_logistic`, min over folds: PXD006140 +0.0395 [−0.0032, +0.0947] crosses zero; PXD024061 +0.0766 [+0.0115, +0.1538] excludes zero. Re-verified on real MMseqs2 clusters (below): same fold fails, delta_ci_lower=−0.0049 |
| 3 | recovery_is_not_training_leakage | **PASS** | control_leakage=[], independent_units=4 |
| 4 | structure_gain_on_structured_subset | **PASS** | structure_gain=+0.0239 (paired per-fold×seed deltas seq_structure−sequence_only, mean +0.0441, 95% CI [+0.0239, +0.0637]) |
| 5 | not_driven_by_single_cluster | **PASS** | single_cluster_driven=False (top cluster = 7 rows, removal retains 77.2% of AP), permutation_p=0.001. Re-verified on real MMseqs2 clusters (below): top cluster = 197 rows (real paralog family), removal *increases* AP (retention 101.7%), permutation_p=0.001 |

Reading of the 3/5: the two failing conditions are different in kind.
Condition 1 is a **data-provenance lock** (frozen config, same lab) that only
new independent studies can open. Condition 2 is a **measured fragility**:
the ranker's margin over the PU baseline is positive in 10/10 runs but does
not survive cluster-level resampling in the PXD006140 fold — the margin
depends on cluster composition, not only on the model.

Condition-4 nuance (reported, not gated): on the structure-covered subset
itself (16 test rows, 14 positives — only 7 AlphaFold structures cover the
benchmark) the paired delta is +0.009 with CI [−0.087, +0.056] —
inconclusive at that n. The run-level gain is real, but it cannot be
attributed to the 14 covered positives alone; the masked structure branch
also acts as a regulariser on the 96% of rows without structures.

## Per-site score regeneration (P2)

`scripts/score_release.py` (CPU, deterministic) re-trains all three arms per
fold×seed and writes aligned per-site scores to
`results/external_validation/pu_ranker_v1/scored/{model,ablated,baseline}.tsv`.
Sanity: regenerated mean test APs match `metrics.tsv` within device tolerance
(e.g. PXD024061 fold 0.1205 CPU vs 0.1208 GPU). Seed-ensembled pooled AP
0.097, cluster bootstrap 95% CI [0.049, 0.108].

**Cluster-granularity remediation (2026-07-24, DONE)**: `protein_clusters_v1.tsv`
(`data/processed/clusters/`) was a **singleton placeholder** (one cluster per
protein, `_build_singleton_cluster_file`); real MMseqs2 clustering required a
Linux host, which this machine did not have. Resolved via a two-hop SSH
session to a lab A100 server (no local Linux/WSL needed — see
`docs/ops/remote_a100_via_lab_jump_host.md` for the reusable procedure):
`mmseqs easy-cluster` on the full Arabidopsis reference proteome (54,646
sequences, 30% identity, ≥50% coverage, `--cov-mode 0`, matching the TDD
Codex §4.3 default) produced **13,967 real clusters** in 6.9 seconds, written
to `data/processed/clusters/protein_clusters_v2.tsv` (SHA256
`e13d16ae…967197`). `protein_clusters_v1.tsv` stays frozen and unused going
forward; v2 is now the cluster file for anything that needs real homology
grouping. One benchmark accession (`P42737-2`, an isoform suffix not present
in the canonical reference-proteome FASTA) has no v2 cluster entry and
defaults to the `train` split via the existing `.get(acc, "train")` fallback
in `_train_val_test_rows` — a pre-existing benchmark data quirk, not a
clustering defect; it affects one positive row.

Conditions 2 and 5 were re-verified end-to-end on real clusters
(`scripts/score_release.py --clusters protein_clusters_v2.tsv` →
`scripts/validate_external.py`, both rerun with the AlphaFold structure
registry pinned to the frozen 7-structure release state so only the cluster
variable changed — see `results/external_validation/pu_ranker_v1_clusterv2_check/`):

| Condition | Singleton (v1, official) | Real MMseqs2 clusters (v2, verification) | Changed? |
|---|---|---|---|
| 2 — effect CI | PXD006140 +0.0395 [−0.0032, +0.0947] (crosses 0); PXD024061 +0.0766 [+0.0115, +0.1538] (excludes 0) | PXD006140 +0.0395 [**−0.0049**, +0.0931] (crosses 0); PXD024061 +0.0766 [+0.0128, +0.1495] (excludes 0) | No — same fold fails, marginally wider on the low end |
| 5 — single-cluster dominance | top "cluster" = 7 rows (one protein, O03042); removal retains 77.2% of AP; permutation p=0.001 | top cluster = **197 rows, a real 30%-identity paralog family (`Q3E937`)**; removal *increases* pooled AP (0.0974→0.0990, retention 101.7%); permutation p=0.001 | No — passes, and now on a test that can actually detect homology-family concentration |

Point estimates for condition 2 (+0.0395 / +0.0766) and condition 4's
structure gain (unaffected — it is a paired per-fold×seed comparison, not a
cluster-grouped bootstrap) are identical between v1 and v2 by construction:
only the resampling/grouping key changed, not the underlying model scores.
The practical upshot: the singleton placeholder was not hiding a
homology-driven optimism bias on this benchmark — real clustering reproduces
the same STOP-relevant conclusions, and condition 5 is now backed by a
genuinely meaningful family-level test (the previous "cluster" was never
capable of grouping more than one protein's own Cys rows together).

**Consequences for the supplementary Split A track**: rerunning
`pu_ranker_cluster_v1` under real clusters (`pu_ranker_cluster_v2.yaml` /
`split_config_v2.yaml`) produced **noisier, not cleaner** numbers — see
`docs/phase_z_evidence_audit.md` §3.3 for the full table and the mechanism
(real clusters correctly keep large paralog families intact as one
train/val/test unit, e.g. a 179-member family with 5 distinct positive
paralogs landing entirely in the v2 test partition, versus the singleton file
scattering such paralogs independently across all three splits). This is
evidence the remediation was worth doing, not a regression to explain away.

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

## Known-control recovery (from control_recovery.tsv + recovery_v2.json)

Percentile ranks use the repaired shared-PU-scorer procedure in
`scripts/evaluate_known_controls.py` (the original single-class logistic
stand-in was degenerate and could not run; the registry and integrity rules
are unchanged). 2026-07-22: the registry moved to
`src/plantpersulf/evaluation/known_controls.py` and gained the first
same-species (Arabidopsis) independent-lab control; same-species controls
are excluded from the scorer's training unlabeled sample before scoring,
which shifted the tomato percentiles by ≤0.2 points (reference resample).
**2026-07-23**: a second same-species independent-lab control was
registered — PAD3 Cys440 (Zhang et al. 2026, Plant Cell Environ,
doi:10.1111/pce.70593; Pei/Jin lab, Shanxi — a >20-paper independent
Chinese H2S-signaling lineage with zero author overlap with the Seville
Romero/Gotor/Aroca network across a decade of publications, surfaced during
the P1 same-species dataset search below).

| mechanism_lineage_id | gene | status | percentile | independent unit |
|---|---|---|---|---|
| SLWRKY6_H2S_PHOSPHORYLATION | SlWRKY6 (Cys396) | mapped | 18.9% (not recovered) | yes |
| SLERFD2_H2S_ETHYLENE | SlERF.D2 (Cys35) | mapped | 96.5% (recovered) | yes |
| ATG6PD6_H2S_G6PD_SALT | AtG6PD6 (Cys159) | mapped (same-species, NWAFU) | 28.6% (not recovered) | yes |
| PAD3_H2S_HCN_OSMOTIC | PAD3 (Cys440) | mapped (same-species, Pei/Jin lab) | 9.8% (not recovered) | yes |
| BRG3_H2S_UBIQUITINATION | BRG3 | unmappable | — | no |
| ERFD3_H2S_CONTEXT | ERF.D3 | position_shift (unconfirmed) | — | no |

No registered control appeared in training (leakage list empty; AtG6PD6
Cys159 and PAD3 Cys440 are both unlabeled PU-pool rows in benchmark_v1,
verified before registration). Recovery is mixed — 1/4 mapped controls
above the 50th percentile of the unlabeled reference — and is reported
as-is; the sequence-only PU scorer carries no structure or
species-specific features, so low recovery of SlWRKY6, AtG6PD6, and PAD3
is informative about the feature set, not evidence against any mechanism.
PAD3's especially low percentile (9.8%, the lowest of all four mapped
controls) is consistent with this pattern, not an outlier requiring
separate explanation.

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

- **P2 (statistical wiring) — DONE (this update).** Cluster-bootstrap effect
  CI, run-level structure gain, top-cluster dominance, permutation test and
  control percentile ranks are all wired and measured; the record moved from
  1/5 (four unmeasured) to 3/5 (one measured fragility + one structural
  lock).
- **P0 (Phase Z)**: data-resource and systematic evidence-audit deliverable
  + collaboration data-request list (per roadmap, Gate 2 STOP route). The
  measured condition-2 fragility and the condition-4 subset caveat are
  first-class inputs to the audit.
- **P1 (data expansion)**: search for genuinely independent (non-Seville)
  site-level persulfidation datasets — the only path that can flip
  condition 1 and reopen Gate 2. **2026-07-22 update**: PXD072089 (Xie et al.
  2026, PNAS) — rice persulfidome, independent lab/chemistry/species, 897
  coordinate-verified sites on 646 proteins — has been parsed and evaluated
  as a cross-species transfer track (§5.3 of the Phase Z audit). The
  cross-species AP is 0.0618 (1.30x base rate, 3/5 seeds below baseline,
  recall@50=0), consistent with PXD063170's 0.0635 (1.33x, same instability).
  **2026-07-23 update**: joining the PRIDE-deposited `SS-all-peptides.tsv`
  (submitter's full MaxQuant peptide export) as a third source raised the
  coordinate-verified site count 897→929; the rerun is materially unchanged
  (AP 0.0618, delta 0.0166, CI excludes zero). Separately, 455/483 deleted
  UniProt accessions were recovered via UniParc, but re-verification found
  **zero additional positive sites** — the deleted accessions' peptides
  were never called as persulfidation sites in the first place, so sequence
  recovery cannot recover positions the source tables never assigned.
  929 is the ceiling reachable from public+deposited data; closing the gap
  to the paper's claimed 1,691 requires the authors' original site-level
  MaxQuant output, not further sequence-database work. A future
  `pu_ranker_v2` 3-study leave-out config (PXD006140/PXD024061/PXD072089
  mapped to Arabidopsis ortholog space via SD02) was scoped and found
  impractical: only 192/929 rice sites sit on a 1:1 Arabidopsis ortholog,
  and a simple window alignment cleanly maps only ~25 of those to a
  verified Arabidopsis Cys — too few for a standalone fold. PXD072300
  (10 recombinant rice proteins, same paper/lab, in-vitro persulfidation
  validated by PEAKS) was parsed as a same-paper orthogonal control set
  (23 sites / 8 proteins; TKT unmappable) — explicitly NOT an independent
  unit for conditions 1 or 3, used only as an internal-consistency probe
  for the rice track. Note: PXD039999 is the same laboratory
  (Jurado-Flores/Aroca/Romero/Gotor 2023, protein-level public material);
  PXD035795's method record is hosted at
  IDUS Seville (same group ecosystem) and its site-level route requires a
  reviewed DCP/NBF→site mapping rule.
- **P3 (cluster-granularity remediation) — DONE, 2026-07-24.** Real MMseqs2
  clustering (see "Cluster-granularity remediation" above) closes the gap
  flagged in the original P2 update and in `docs/phase_z_evidence_audit.md`
  §3.3. Conditions 2 and 5 re-verified on real clusters; decision unchanged
  (still STOP 3/5). `docs/ops/remote_a100_via_lab_jump_host.md` documents the
  no-local-Linux workflow used (two-hop SSH to a lab A100 server) for reuse
  on future remediation items that need Linux tooling.
- **P4 (structure-coverage expansion, NOT YET VALIDATED — flagged, not
  resolved)**: while re-scoring for P3, `data/registry/alphafold_structures.tsv`
  was found to have grown from 7 to 2,006 registered structures (98.0%
  coverage, uncommitted at time of writing — see
  `docs/phase_z_evidence_audit.md` §5.5.2 for the batch-download work that
  produced it) since the official `pu_ranker_v1` release was scored. A quick,
  informal check (re-running `score_release.py` against the *current*,
  uncommitted structure registry instead of the frozen 7-structure release
  state) showed the `seq_structure` arm's leave-study-out test AP jumping to
  0.05–0.54 across seeds on the PXD006140 fold — an order of magnitude above
  the frozen release's 0.0825±0.0081, and highly seed-unstable. **This number
  is not validated and must not be cited anywhere** — it was produced by an
  informal single-variable-isolation check, not a designed experiment, and
  was never intended to answer "does more structure coverage help." Whether
  this reflects a genuine capability improvement, overfitting to the newly
  available structural signal, or a confound in the feature pipeline is
  unknown. If pursued, this needs its own frozen `pu_ranker_v2` release under
  the full TDD Codex process (frozen split, fixed seeds, no cherry-picking),
  not an ad hoc rerun. Filed here so it is not lost, not acted on.

---

## 2026-08-11 附录：评估架构变更声明 + 蛋白切分展示轨数字

### A. 评估架构变更（用户决策，本文档为决策记录）

| 场景 | 评估方式 | 角色 |
|---|---|---|
| 番茄本地模型 | 5 折同源簇分组（MMseqs2 cluster CV） | 主评估（番茄侧） |
| 拟南芥（现有基准） | 补充重复 5/10 折同源簇分组（`pu_ranker_cluster_cv_v1`） | 开发稳定性（模型选择/消融） |
| **NC 主结论** | **未来新番茄队列 = lockbox**（隐藏验证，标签冻结前不开放） | **主证据——文章核心主张的唯一锚点** |
| leave-one-study-out | **有第二个独立研究后**再做 | 补充/未来 |
| 跨物种（PXD063170/PXD072089） | 额外挑战集 | robustness 展示，不替代外部验证 |

含义：
- 现有 `pu_ranker_v1` leave-study-out 冻结结果**保留为参考轨**（本附录下方表格仍可引用），不再承担 NC 主证据角色；`gate2_decision.json` 未重判、字节不变（2026-08-11 验证）。
- NC 稿件"外部验证/泛化"主张只能引用 lockbox 结果；Level 2 在跨研究证据出现前不可声称（codex §12）。
- lockbox 协议进入阶段 2 冻结文档（见 `docs/superpowers/plans/2026-08-11-nc-submission-execution-manual.md`）。

### B. 蛋白切分展示轨数字（Sul-BertGRU 同口径，2026-08-11 A100 完成）

配置 `configs/experiments/pu_ranker_protein_split_v1.yaml`：随机 20% 蛋白留出 × 10 seeds（10 次重复）、无同源控制、subsample 1:20。结果 `results/experiments/pu_ranker_protein_split_v1/{metrics.tsv,manifest.json,summary.json,summary.md}`。

| ablation | 蛋白切分 test_ap (10 seeds) | LSO 参考（冻结 5 seeds） | Split A cluster 参考 |
|---|---|---|---|
| sequence_only | 0.0763 ± 0.0202 | 0.0613 | 0.0504 |
| seq_esm | 0.0493 ± 0.0145 | 0.0461 | 0.0399 |
| **seq_structure** | **0.1135 ± 0.0332** | **0.1017** | **0.0837** |
| full | 0.0358 ± 0.0063 | 0.0484 | 0.0342 |
| no_plddt | 0.0382 ± 0.0114 | 0.0412 | 0.0342 |
| no_accessibility | 0.0356 ± 0.0048 | 0.0406 | 0.1695 ⚠️(单 seed 离群) |
| no_study_context | 0.0478 ± 0.0094 | 0.0545 | 0.0435 |

- 结构增益（seq_structure − sequence_only，10 paired deltas → `paired_delta_ci`）：**+0.0372 [95% CI 0.0167, 0.0606]** —— CI 排除零
- 效应 delta（seq_structure − pu_logistic，Route A 基线逐 seed 配对）：**+0.0594 [95% CI 0.0400, 0.0785]** —— CI 排除零；基线配置 `pu_ranker_protein_split_v1_baseline.yaml`（切分/抽样与展示轨逐 seed 一致，本地 CPU 10 seeds）
- base rate 0.0476；seq_structure = 2.4x base rate

**解读（如实）**：蛋白切分（无同源控制）数字高于 LSO/cluster 参考轨，且两个关键 CI 排除零——这是"数据集内可学习信号"在同行口径下的量化，同时量化了**同行口径的乐观上限**（同源蛋白跨划分的泄漏贡献）。**不构成跨研究证据**：本轨 structurally 排除于 Gate 2（测试锁定），条件 1 数据来源锁不变（`studies_are_independent: false`），Gate 2 仍为 STOP 3/5。对外展示必须附 limitation verbatim（见 `summary.json` / 两轨政策文档）。

### C. Route 记录
- 基线缺口（展示轨 `models: [structure_ranker]` 未含基线）→ **Route A 执行**：冻结 `pu_ranker_protein_split_v1_baseline.yaml`，本地 CPU 10 seeds 完成，逐 seed 配对成立（切分/抽样模型无关性，`_run_protein_split_experiment`）。Route B 措辞未启用。
- 开发稳定性轨 `pu_ranker_cluster_cv_v1`（5 折 × 5 重复同源簇分组）运行记录见 `docs/superpowers/plans/2026-08-11-nc-submission-execution-manual.md` 阶段 1。

### D. 措辞锁定（本附录数字可伴随的唯一表述）
> 蛋白切分（Sul-BertGRU 同口径）数字显示数据集内存在可学习信号（结构增益与效应 delta 的 95% CI 均排除零），但该轨与 LSO/cluster 轨同属数据集内评估，不构成跨研究、跨实验室或跨物种证据；NC 主结论锚定于未来番茄 lockbox 队列。

---

## 2026-08-13 附录：multispecies v2 双轨结果与 Gate 2 关系

`configs/experiments/multispecies_v2_global_clusters_v11.yaml` 已完成：

- **严格同源簇 development 轨**：5 folds（结构：`strict_cluster_holdout|...`）。
- **文献同口径随机蛋白轨**：6 models × 10 seeds，含 Sul-BertGRU baseline（结构：`literature_random_protein|...`）。

两轨均**结构性地排除于 Gate 2**：

1. 数据仍来自同一实验室/同一化学体系的拟南芥 + 水稻 + 番茄整合集（`studies_are_independent: false` 未变）。
2. 模型名前缀 `strict_cluster_holdout` 与 `literature_random_protein` 均无法通过 `leave_<study>_out` 前缀检查（`tests/release/test_gate2_ignores_within_dataset_split_metrics.py` 已锁定）。
3. 配置中显式声明 `limitation: INTERNAL MULTISPECIES EVALUATION ONLY; must never enter Gate 2` 与 `comparator_policy.gate2_eligible: false`。

因此，multispecies v2 结果可作为**数据集内可学习信号与文献同口径改进**的证据，并支撑番茄盲法 lockbox 候选工具的选择，但**不能**用于声称跨研究、跨实验室或跨物种预测能力。NC 主结论仍唯一锚定于未来独立番茄 lockbox 队列的结果。

