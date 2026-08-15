# 分物种结构标定 — 新发布候选变更预注册（草案，待联合签署）

> **状态**：DRAFT（建模方起草，待联合签署后定稿）。
> **发布号**（占位，签署时定稿）：`multispecies-v2-candidate-release-v2`
> **变更类别**：结构分支**输入投影**（特征缩放规则），**非新架构**。评估层诊断 B 已
> 证明其方向（claim_class `diagnostic_only`）；本文件把正式发布评估的**统计口径写死**，
> 供联合签署后按纪律执行。
> **执行门槛**：register §1——新一代模型立项唯一路径为 **Gate 3 盲法验证完成之后**；
> 本文件在此前仅锁定协议，不触发任何训练/候选表改动。
> 关联：`2026-08-13-model-improvement-register.md`（§2.2 W4、§3）、
> `2026-08-15-regime-routing-moe-design.md`（§4B/§5）、
> `2026-08-15-model-accuracy-campaign.md`（W1 + DECISION_SUMMARY）、
> `results/experiments/accuracy_campaign/w4_per_species_scaling/`。

---

## 1. 背景与证据链

| 阶段 | 证据 | 出处 |
|---|---|---|
| W1（阴性） | v3 番茄结构覆盖 0→90%，文献轨 macro 0.3863→0.2117、严格轨 0.5383→0.3257；tomato 无增益（严格轨 −48%） | W1 report；register §2.2 L1 |
| 机制线索 #1 | 全局 `_BranchScalers.struct` 在全部 present-structure 行拟合，v3 后番茄占 ~78%，arab/rice 结构输入被番茄统计重缩放 | W1 report §3 结论 |
| 诊断 B（本 session，阳性） | 分物种结构标定：3-crop macro 0.2117→**0.3799**（v2 基线 0.3863 的 98%），配对 Wilcoxon **p=0.0098**（9/10 seeds）；arab +0.265（9/10 up）、rice +0.242（8/10 up）、tomato 近地板持平、magnaporthe −0.034 非显著 | `w4_per_species_scaling/comparison.json` |

诊断 B 是**已预注册的诊断**（MoE 设计文档 §4B，成功标准在运行前写死），claim_class
`diagnostic_only`，未触冻结包/候选表/盲法队列。本文件在其之上把**正式发布评估**的
协议定稿。

## 2. 变更定义（精确，签署后不可改）

**范围**：`structure_ranker` 结构分支的**输入缩放规则**。

- **现状（冻结 v11）**：结构特征 → 单一全局 `TrainOnlyScaler`（全部 present-structure
  行 pooled 拟合）→ 分支输入。
- **变更**：结构特征 → 先按**物种** z 标定（每物种在 train present 行拟合独立
  `TrainOnlyScaler`）→ 再进入模型（rank 内部全局二次 fit 因每物种 mean0/var1 近似
  恒等，不重新引入跨物种拉偏）。
- **不变**：模型架构（gated-fusion PU ranker，hidden-16）、全部超参（v11 冻结）、
  特征集、缺失结构掩蔽语义（无结构行结构分支强制为零，绝不均值插补）、训练面板
  构建、ESM/study-context 分支关闭。
- **fallback 规则**：物种无 present-structure train 行（或行 masked）→ 恒等返回
  （无伪造统计量），与诊断 B 实现一致。

**与 L8 分物种 isotonic 校准的关系**：**正交**。L8 是 Top-K **池化后处理**（排序后
单调校准），本变更是**模型输入投影**（训练侧）。两者可独立预注册，也可组合（组合
时各自协议分别锁定，不互相干扰）。

## 3. 统计口径（写死；K 与 claim gate 同口径）

三轨全部在**冻结内部轨**上评估，不触碰盲法数据：

### 3.1 文献轨（主轨，10-seed 配对方差）

- 面板：v11 `literature_random_protein` 10-seed 面板，`panel_sha256` 与冻结 summary
  逐 seed 一致（gate）。
- 指标：per-species **AP@K / recall@K**，K ∈ {50, 200}（claim gate 口径），test 分区。
- 对比：per_species 臂 vs global 臂（global 臂须逐 seed 复现 W1 v3，reproduction
  gate 容差 1e-6），每 seed 配对。
- **主要检验**：per-species AP@200 的 10-seed 配对 **Wilcoxon signed-rank**（单侧，
  显著为正），并报告 per-species 逐物种 delta 的方向一致性（n_up/n）。

### 3.2 严格轨（5-fold dev，折叠配对）

- 面板：严格轨 5-fold dev CV（v11 同构）。
- 指标：validation AP，3-crop macro + per-species。
- 对比：per_species 臂 vs global 臂，**折叠配对**差异。

### 3.3 已知对照（n=12，hit@2 不降）

- 12 个 mapped controls 全部 Cys 冻结架构蛋白内排序。
- 指标：hit@2（真位点进入蛋白内 top-2 的比例）、总负担（相对随机 58.5）。

## 4. 成功 / 失败定义（写死）

**成功（全部满足）**：

1. **文献轨**：per-species AP@200 的 10-seed 配对 Wilcoxon **显著为正（p<0.05）**，
   且 **arabidopsis 或 rice 至少其一**逐物种 AP 配对方差显著改善；
2. **tomato 不灾难性恶化**：tomato per-species AP 配对 delta 不显著为负（近地板
   ±0.005 视为持平）；
3. **严格轨**：5-fold 折叠配对 macro 不显著恶化；
4. **已知对照**：hit@2 不降（≥ 现状）。

**失败（任一不满足）**：登记为阴性结论（register 保留，不删除），不发布、不改候选表。

## 5. 数据、产物与发布

- 结构注册表：`ALPHAFOLD_STRUCTURES_RELEASE_V3`（SHA256 `60aa3af0…`，已注册
  `model_inputs.tsv`），与诊断 B / W1 同一版本。
- 训练面板：v11 冻结内部轨（panel gate）。
- **新发布产物**（签署后、Gate 3 后执行）：
  - 新发布号（占位 `multispecies-v2-candidate-release-v2`）；
  - 新 `fit_manifest`（面板/种子/超参/代码修订/SHA256）；
  - 新候选打分表 + 匹配对照 + 特征 schema + 盲法推理脚本（镜像 release v1 结构）；
  - 新 SAP / 分析计划 JSON + 联合签署 + OSF 注册（embargo 至揭盲）。

## 6. 纪律承诺（本文件即约束）

1. **不触碰冻结 v1**（bundle `ab8a0353…`、候选表、盲法队列、已消费 frozen test）。
2. **不在盲法数据上做任何反向开发**；Gate 3 完成前不重跑候选表。
3. **先发布预注册（联合签署），再训练**；训练产物全新发布号。
4. 文献轨是 1:20 PU 面板，池化 Top-K 数字不可直接外推全蛋白组候选表——候选表层级
   验证在 Gate 3 后单独做（与 L8 同局限，如实记录）。
5. 结果无论正负均写入 register §4；阴性不删除。

## 7. 状态与下一步

- [ ] 建模方起草（本文件，2026-08-15）
- [ ] 建模方 PI 审阅定稿
- [ ] 张华实验室确认（对照比 / K / 功效口径，若涉及候选表变更）
- [ ] 联合签署 → 锁定统计口径
- [ ] Gate 3 盲法验证完成后执行：训练 → 三轨评估 → 新候选表 → SAP/OSF

---

**附录：诊断 B 关键数字（预注册草案的证据来源）**

- global 臂 3-crop macro = 0.2117（=W1 v3 精确复现，reproduction gate 10/10）。
- per_species 臂 3-crop macro = 0.3799；配对 Wilcoxon p=0.0098（9/10 seeds 正向，
  mean delta +0.1682）；4-species p=0.0195（8/10）。
- 逐物种（10 seeds 均值）：arab 0.2302→0.4957（9/10 up）、rice 0.3844→0.6265
  （8/10 up）、tomato 0.0204→0.0175（6/10 up，近地板持平）、magnaporthe
  0.4536→0.4200（4/10 up，0 结构行 spillover，非显著）。
- 产物：`results/experiments/accuracy_campaign/w4_per_species_scaling/`
  （comparison.json + 每 seed 每臂 test/validation TSV）；claim_class
  `diagnostic_only`。
