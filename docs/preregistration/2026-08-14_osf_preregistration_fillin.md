# OSF 预注册填写指南 — Tomato 盲法队列(Gate 0 → Gate 3)

生成日期:2026-08-14
用途:把 Gate 0 冻结的 SAP + 分析计划在 OSF Registries 注册,拿回 DOI,写进 manifest.json 与 PP 手稿。
内容来源:全部英文填写内容抽取自已签署冻结的 `sap_protocol.json` / `analysis_plan.json` / `manifest.json`(v1.0, signed_2026-08-13,git tag `gate0-release-v1-20260813`)。

---

## 0. 分工边界

我能做的:全部表单文字(本文件)、上传文件清单与打包。
你必须自己做的(账号身份绑定):注册 OSF 账号、登录、点提交。全程复制粘贴,预计 20–30 分钟。

## 1. 点击路径

1. 打开 https://osf.io ,注册账号(建议用学校邮箱;姓名用你论文署名的英文拼写)。
2. 顶部 `Create new project`,项目名建议:
   `PlantPersulf — Preregistered tomato blind cohort (multispecies-v2-candidate-release-v1)`
3. 进入项目,左侧 `Files` → `OSF Storage`,上传第 3 节的全部文件。
4. 顶部 `Registrations` 标签 → `New registration` → 选择模板 **`OSF Preregistration`**(完整版,不是 AsPredicted 简版)。
5. 按第 4 节逐栏粘贴。
6. 提交前设置 **embargo(延期公开)**——见第 2 节,这一步关乎盲法完整性,不要跳过。
7. License 选 **CC-BY 4.0**;Subjects 加 `Life sciences / Plant sciences`;Tags 建议:`persulfidation`, `proteomics`, `tomato`, `preregistration`, `blind cohort`。
8. 可选:Add contributor 填张华的邮箱(需要对方有 OSF 账号,可后置)。
9. 提交。DOI 即时生成,形如 `10.17605/OSF.IO/XXXXX`。

## 2. ⚠️ Embargo 设置(保护盲法,必读)

上传的文件里有候选清单和 blind_id 钥匙。**如果注册立即公开,张华实验室理论上能看到候选/对照标签,盲法就破了。**

- 提交时选择 **Embargo**,结束日期设为 **2027-06-30**(或揭盲完成日,取较晚者;OSF 最长支持 4 年)。
- Embargo 期间:时间戳已锁定(预注册效力完整),内容对外不可见;到期自动公开,届时揭盲早已完成,公开反而成为"标签分配从未改动"的证据。
- 如果你希望部分内容立即公开(如 SAP 正文):也可以分两个 project——公开 project 放 SAP/analysis_plan/manifest,embargo project 放 assay_list 和 key。简单起见建议整体 embargo。

## 3. 上传文件清单(先传文件,再建注册)

| 文件 | 来源路径(本仓库) | 说明 |
|---|---|---|
| sap_protocol.json | results/candidates/multispecies_v2_candidate_release_v1/ | 已签署 SAP v1.0 |
| analysis_plan.json | 同上 | 已签署分析计划 v1.0 |
| co_signature_package.md | 同上 | 联合签署包(含签署记录) |
| manifest.json + SHA256SUMS | 同上 | 冻结清单与全部哈希 |
| feature_schema.json | 同上 | 特征定义 |
| inference_script.py | 同上 | 冻结推理脚本 |
| power_table.json | 同上 | 144 格功效表 |
| top_k_candidates.tsv | 同上 | 候选表(embargo 保护) |
| matched_controls.tsv | 同上 | 对照表(embargo 保护) |
| assay_list_blind.tsv | results/handoff/zhang_lab_blind_cohort_v1/ | 盲法测定清单(embargo 保护) |
| blind_id_key.DO_NOT_SEND.tsv | results/handoff/zhang_lab_blind_cohort_v1/ | 标签钥匙(embargo 保护;到期公开即防篡改证据) |
| handoff_manifest.json | results/handoff/zhang_lab_blind_cohort_v1/ | 交接清单(含钥匙 SHA256) |
| README.md | results/handoff/zhang_lab_blind_cohort_v1/ | 湿实验室操作说明 |
| structure_ranker_bundle.pt + fit_manifest.json | results/candidates/.../model_weights/ | 冻结模型本体(很小) |

我已把以上文件复制到 `results/handoff/osf_registration_bundle_v1/`,直接整个文件夹拖进 OSF Storage 即可。

## 4. 表单逐栏内容(英文,直接粘贴)

> 注:OSF 表单字段名如与下列略有出入,按含义对应即可。以下全部文字已通过本仓库的两道措辞门(verify_predictive_claims / find_forbidden_external_claims)。

### Title

```
Frozen, hash-verified prioritization of tomato persulfidation candidate sites: a preregistered independent blind-cohort enrichment test
```

### Description

```
This registration locks the design, cohort, analysis plan, and frozen computational artifacts for an independent blind-cohort test of a persulfidation-site prioritization release (multispecies-v2-candidate-release-v1). The release ranks cysteine sites of the tomato (Solanum lycopersicum) reference proteome using a frozen ranking model (structure_ranker bundle, SHA256-pinned) fit on previously published, fully documented persulfidation datasets. The ranking model, features, candidate table, matched-control table, and analysis plan were frozen and co-signed by both participating parties on 2026-08-13 (git tag gate0-release-v1-20260813) before any blind assay data existed. The wet-lab partner (Zhang Hua lab, South China Agricultural University) will assay 593 sites (200 candidates + 393 matched controls) under blind IDs. The primary endpoint is enrichment of experimentally confirmed sites among candidates relative to matched controls. All results, including null and negative outcomes, will be reported.
```

### Hypotheses

```
H1 (primary, directional): within the assayed cohort, the fraction of experimentally confirmed persulfidation sites is higher in the candidate arm (Top-200 ranked sites) than in the matched-control arm. H1 is supported only if all three pre-registered success criteria hold simultaneously: (i) one-sided Fisher's exact test p < 0.05; (ii) odds ratio point estimate > 1; (iii) at least 2 confirmed candidate sites.
H2 (secondary): confirmed sites concentrate in the upper rank percentiles of the frozen candidate ranking.
H3 (secondary, robustness): the enrichment is not driven by any single protein or any single homology cluster (leave-one-out sensitivity).
We explicitly register the possibility that H1 is not supported; the study will be reported regardless of outcome.
```

### Design Plan — Study type

```
Other — a frozen computational ranking followed by an independent wet-lab confirmation assay of a fixed site cohort.
```

### Design Plan — Blinding

```
Three-way separation. (1) The wet-lab partner receives the cohort only as blind IDs (deterministic shuffle, fixed seed 20260813) and generates assay data without knowledge of candidate/control status. (2) The modeling team holds the label key (its SHA256 is recorded in the handoff manifest) and has no access to assay data until the pre-registered unblinding step. (3) Unblinding requires two independent persons to verify all artifact hashes and inference-script determinism before labels are released; the scoring output for the blind cohort is written and hash-recorded before labels are received.
```

### Design Plan — Study design

```
Fixed two-arm site cohort: 200 candidates / 393 matched controls, 593 assayed sites total. Candidates are the top 200 entries of a frozen ranked table covering all tomato reference-proteome cysteine sites, excluding every site in the release training panel. Controls are matched per candidate on structure-availability stratum and Euclidean distance <= 0.25 SD in z-scored sequence-feature space (at most 5 per candidate in the control table, no reuse); the analysis uses exactly the 2 nearest controls per candidate, fixed at signing. 4 of the 200 candidates have fewer than 2 controls within the frozen caliper (3 with 0, 1 with 1); they remain in the candidate arm, the control arm is simply smaller, and Fisher's exact test handles the unequal arm sizes. Pre-registered exclusions from the primary analysis: any site in the release training panel; any site on a protein whose persulfidation status is published in kiae271 (Zhang et al. 2024) — the latter are analyzed only as a pre-registered sensitivity stratum. Assay-specific positive and negative process controls are defined by the wet-lab partner and recorded before any measurement.
```

### Design Plan — Randomization

```
None. Arm assignment is deterministic (frozen ranking plus frozen matching rule). Blind IDs are assigned by a deterministic shuffle with fixed seed 20260813; the ID-to-label mapping file is hash-recorded and held under embargo with this registration.
```

### Sampling Plan — Existing data

```
Registration prior to creation of data.

Explanation: the outcome data of this study (assay confirmation labels for the 593 sites) do not exist yet; they will be created by the wet-lab partner after this registration. The ranking model was fit on previously published datasets; those training data are fully documented in the uploaded manifest and are not outcome data of this study.
```

### Sampling Plan — Data collection procedures

```
The wet-lab partner (Zhang Hua lab, South China Agricultural University) assays all 593 sites with its established persulfidome workflow, working solely from the blind-ID assay list. The judgment rule is fixed before measurement: a site counts as confirmed only if it meets the laboratory's pre-specified detection criteria; non-detection in a valid assay is counted as not confirmed (non-detection does not by itself establish absence of the modification).
```

### Sampling Plan — Sample size

```
593 protein cysteine sites: 200 candidates and 393 matched controls, assayed in a single round.
```

### Sampling Plan — Sample size rationale

```
A simulation grid over candidate count K in {50, 100, 200, 300}, controls per candidate r in {1, 2, 5}, background confirmation rate p0 in {0.01, 0.02, 0.05}, and assumed odds ratios in {1.5, 2, 3, 5} (one-sided Fisher's exact test, alpha = 0.05, 10,000 simulations per cell, seed 20260813; the full 144-cell table is uploaded as power_table.json). At a fixed assay budget, a larger K with fewer controls dominates a smaller K with more controls. The selected design (K = 200, r = 2, nominal 600 sites) reaches power 0.68 at p0 = 0.02 and OR = 3, and 0.98 at OR = 5, with 11.5 expected candidate hits at OR = 3 against the >= 2 confirmed-site criterion. Because 4 of 200 candidates have fewer than 2 matched controls within the frozen 0.25 SD caliper, the actual cohort is 593 sites; the power impact is negligible.
```

### Sampling Plan — Stopping rule

```
The cohort is fixed at 593 sites in a single assay round. No post-hoc sample additions are permitted after unblinding.
```

### Variables — Manipulated variables

```
None. The grouping variable (candidate vs matched control) is assigned by the frozen ranking and the frozen matching rule; the investigators manipulate no experimental variable.
```

### Variables — Measured variables

```
Per site: persulfidation confirmation status (binary: confirmed / not confirmed) from the wet-lab assay; assay process-control pass/fail.
```

### Variables — Indices

```
The primary 2x2 table (arm x confirmation status); the odds ratio with exact 95% confidence interval; the rank percentile of each confirmed site within the frozen candidate ranking.
```

### Analysis Plan — Statistical models

```
Primary: one-sided Fisher's exact test on the 2x2 table (directional alternative: enrichment in the candidate arm), plus the odds ratio with exact 95% CI.
Pre-registered sensitivity analyses: (a) exclusion of sites on proteins present in the release training panel (protein-level leakage guard); (b) exclusion of sites on proteins with published persulfidation evidence in kiae271 (published-positive contamination guard); (c) leave-one-out over candidate proteins and over homology clusters (not single-protein / single-cluster driven); (d) per structure-availability stratum — the frozen structure registry contains no tomato accessions, so every tomato candidate is scored with the structure branch masked and the structure branch cannot contribute to the ranking in this release.
Secondary: rank-percentile calibration of confirmed sites (Kolmogorov-Smirnov test and percentile bootstrap); process-control pass rate reported alongside the primary result.
```

### Analysis Plan — Transformations

```
None.
```

### Analysis Plan — Inference criteria

```
One-sided alpha = 0.05. The primary endpoint is supported only if all of the following hold: p < 0.05; odds ratio point estimate > 1; at least 2 confirmed candidate sites. Effect sizes are reported with exact 95% confidence intervals. The study is reported regardless of outcome.
```

### Analysis Plan — Data exclusion

```
Pre-registered exclusions from the primary analysis: sites in the release training panel; sites on proteins with persulfidation published in kiae271 (analyzed only as a sensitivity stratum); sites whose assays fail the wet-lab partner's pre-specified quality criteria (reported with provenance, never replaced).
```

### Analysis Plan — Missing data

```
Non-detection in a valid assay is counted as not confirmed. Sites with failed assays are excluded from the 2x2 table and reported with provenance; they are not replaced and no new sites are added.
```

### Analysis Plan — Exploratory analysis

```
One pre-registered exploratory check: whether MC-dropout uncertainty stratifies confirmed vs unconfirmed candidates within the Top-200. Any analysis not listed in this registration will be labeled exploratory in the report.
```

### Other — Anything else you would like to pre-register?

```
Frozen artifacts and verification: the model bundle, fit manifest, candidate table, control table, feature schema, and inference script are SHA256-pinned in the uploaded manifest.json (status: frozen; git tag gate0-release-v1-20260813).
Unblinding procedure: two independent persons verify all artifact hashes and inference-script determinism; the scoring output for the blind cohort is written and hash-recorded before labels are released; the pre-registered analysis then runs verbatim.
Forbidden after unblinding: any model, feature, threshold, Top-K, or control-set change; adding samples to a failing endpoint; re-ordering the candidate table.
Co-signature: the protocol (sap_protocol.json) and analysis plan (analysis_plan.json) were co-signed by both parties on 2026-08-13; any amendment requires a dated, hashed re-signature by both parties.
Claim discipline: this study tests within-cohort enrichment of a candidate-organisation tool; benchmark metrics from model development are within-dataset evidence and do not enter the claims of this study.
Reporting: all outcomes, including null and negative results and all exclusions, will be reported with full provenance.
```

## 5. 提交后

把 OSF DOI(`10.17605/OSF.IO/XXXXX`)发给我,我会:
1. 写入 `manifest.json` 的 `pending_actions`(替换唯一的 pending 项);
2. 写入 PP 手稿 Methods 的 preregistration 段;
3. 更新 roadmap 日志与记忆文件。
