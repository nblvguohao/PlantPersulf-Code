# 模型改进登记册 — 延迟开发清单（Gate 0 冻结期）

> **状态**：登记不实施。本文件记录 `structure_ranker` 已知的、技术上可行的改进方向及其证据，
> 但**明确禁止在当前冻结期内实施其中任何一项**。它是信息存档，不是行动指令。
>
> **与路线图的关系**：承接 `docs/superpowers/plans/2026-08-13-nc-entry-gate-roadmap.md`
> 全局约束 1（禁止反向开发）与约束 3（不得追加无边界模型开发"抢救"）。
> 本登记册的目的正是把"以后可能的改进"与"现在必须冻结的模型"在纸面上分开。

## 1. 为什么现在不做（不可违反）

- 内部基准（文献轨 macro AP 0.386）**不是目标**；盲法番茄队列富集才是。针对内部基准继续调模型，
  本质是对基准的过拟合，且不会转化为盲法证据。
- 任何改进 = 特征/模型/超参变化 → 必须新建发布号 + 新预注册 + 重新联合签署（见
  `results/candidates/multispecies_v2_candidate_release_v1/` 冻结纪律）。
- 解锁条件（唯一路径）：**Gate 3 盲法验证完成之后**，作为"新一代模型"立项，配完整的新预注册流程。

## 2. 改进杠杆清单

### L1. 番茄 AlphaFold 结构覆盖（最大杠杆）
- **现状**：`data/registry/releases/alphafold_structures_release_v2.tsv` 中番茄 accession 为 **0 条**；
  `top_k_generation_manifest.json` 记录番茄候选结构覆盖 0/179,736。番茄排序 100% 由 3 维序列特征驱动，
  结构分支对番茄完全掩蔽。
- **预期收益**：结构分支（contact proxy + pLDDT，缺失掩蔽机制已就绪）对番茄恢复输入；
  kiae271 分析显示张华团队从未对自身位点做过结构分析（记忆记录），是差异化角度。
- **代价/风险**：属于特征变化；需新注册表发布、重跑候选表、新预注册。开发轨上结构缺失率
  386,009/389,609，结构分支在训练中贡献有限，收益不保证。
- **触发条件**：Gate 3 之后；或作为独立方法学论文（bib 级）的候选主题。

### L2. ESM 分支（use_esm: 0 → 1）
- **现状**：v11 冻结为 `use_esm=0`（`literature_baseline_parameters.structure_ranker`）；
  同一基准上 esm_linear_head 基线 mean AP 0.0466（弱）。
- **预期收益**：gated 融合中的 ESM 分支可能捕获序列语义；但 v11 消融已选 0，证据不支持短期翻盘。
- **代价**：推理成本、特征管线依赖 esm2_t33 窗口产物；需新预注册。

### L3. Study-context 分支（use_study_context: false → true）
- **现状**：冻结关闭；盲法队列无 study 标签，此分支对番茄盲法无直接收益。
- **预期收益**：仅在开发轨上有 study 归一化收益，属"基准内增益"，不是盲法杠杆。

### L4. 训练数据扩充
- **现状**：2,334 阳性（4 物种）。已知可用的公开位点：
  - **kiae271 119 肽表**（`docs/zhanghua/_extract_kiae271/dataset_s1.tsv`，99 验证位点，
    与当前 multispecies v2 训练面板无重叠或仅部分重叠，需核验）；
  - COPLBI 综述（Corpas 等，审稿中）~35 个多物种位点表；
  - 张华实验室的原始 MS 搜索参数（外联请求中，`docs/zhang_outreach_draft.md`）。
- **预期收益**：番茄阳性样本扩充 10 倍级，是"同物种训练数据"缺口的直接补丁——
  记忆结论：跨物种迁移信号 ≈ 0（leave-one-species-out 富集 0.90–1.00×），同物种数据是唯一实证有效方向。
- **触发条件**：需先解决位点级 provenance 核验（对齐 119 表与参考蛋白组版本）；
  属新模型代际。

### L5. 序列特征扩充（机制驱动）
- **现状**：3 维（疏水性、Cys 密度、局部正电荷密度）。已实现但未入冻结轨道的
  `local_positive_charge_density`（硫醇 pKa 代理）是同一思路的延续。
- **候选**：跨物种 Cys 位置保守性（MSA）、硫醇 pKa 的经验代理、窗口内芳香族/极性组成。
- **风险**：单特征增益小（v1 经验：特征增量收益 < 0.01 AP 量级）；易过拟合开发轨。

### L6. 架构调整
- **候选**：基于图的位点上下文、1D-CNN 窗口编码、Sul-BertGRU 风格适配器。
- **现状**：Codex 约束"最小实现"；v11 竞品轨显示复杂基线（sul_bertgru）在同基准上最弱（0.0213），
  架构复杂度与收益无正相关证据。

### L7. PU/排序方法学
- **候选**：标签频率 c 的稳健估计（多折交叉）、分数校准、Top-K 不确定度引导重排。
- **风险**：属"阈值/重排"范畴，冻结纪律对 Top-K 重排尤其严格；任何改动 = 新预注册。

## 2.1 诊断校准（2026-08-13 模拟盲测后）——杠杆优先级更新

kiae271 模拟盲测（`scripts/evaluate_kiae271_with_release_bundle.py`，
99 个已发表番茄位点 × 冻结 bundle，百分位在候选表内）给出三点新事实：

1. **训练数据事实**：v2 面板的 79 个番茄阳性 = kiae271 位点（100% 重叠）。
2. **干净子集（20 个未训练位点）**：65% 排在中位数以上（平均百分位 60.5，
   双侧 p=0.26，n=20）——信号存在但弱。
3. **关键校准**：**0/99 位点进入候选表 Top-2000**。模型把真实位点提到
   列表中部，但不堆到顶端。原始分数分布极扁平（0.005–0.305 覆盖 18 万
   位点，Top-200 挤在 0.295–0.305 之间），顶部被"特征极端但未必真实"
   的位点占据。

**由此重排的杠杆优先级（仍是登记不实施）**：

- **L1 升级（番茄结构覆盖）**：诊断直接支持"番茄特征饥饿"假说——结构分支
  对番茄贡献为 0，而它对其他三个物种的排序贡献显著。补 AlphaFold 番茄结构
  是当前最有针对性的杠杆，但属特征变化 → 新发布 + 新预注册。
- **新增 L8（分数校准/选择策略）**：信号在中段而非顶端，提示问题可能部分
  在"选择规则"而非"模型"：候选按原始分数取 Top-N 可能不是最优选择。
  候选策略：MC-dropout 不确定度引导（避开高不确定的极端位点）、
  先验校准后验（positives 仅占 ~0.6%，分数被压扁）、按分数带分层抽样。
  **任何选择规则改动 = Top-K 变化 = 新预注册**。
- **L4 更新（训练数据）**：99 个 kiae271 位点中 20 个从未进入训练——
  新一代模型可纳入为训练阳性，但代价是失去这批"干净校准位点"。
- **L2/L5/L6 不变**：诊断未提供针对 ESM/序列特征/架构的额外证据。

## 2.2 评估战役 2026-08-15（进行中）——结构与结论更新

用户（建模方 PI）于 2026-08-15 提出"针对数据结构、模型、已有类似位点预测工具提升
准确度"，确认三线并进（W1+W2+W3）。计划文档：
`docs/superpowers/plans/2026-08-15-model-accuracy-campaign.md`。
全部评估在冻结内部轨上运行，产物在 `results/experiments/accuracy_campaign/`，
**不触碰冻结包/候选表/盲法队列，不再解锁 frozen test**。

### 基础设施：文献轨逐位点重跑管线（已建成并验证）

- `scripts/accuracy_campaign/run_literature_ranker_scores.py`：逐 seed 重建 v11 面板
  （同生产函数），`panel_sha256` 与冻结 summary **10/10 一致**；structure_ranker
  test 逐物种 AP 与 v11 **逐位一致**（reproduction_mismatches 0/10）→ 管线位级等价，
  此后任何臂（use_esm、新结构注册表）的分数均可与 v11 同口径直接对比，且顺带补上
  生产轨丢弃的**逐位点分数 + MC-dropout uncertainty**（此前只有聚合指标落盘）。

### W2 已有工具集成：Sul-BertGRU 堆叠 —— 阴性，不立项

- 等权秩聚合（sr × sul）：macro AP 0.3863 → **0.1181**（被弱模型拖累）。
- 验证集网格最优凸组合：10/10 seeds 最优 alpha ≈ 0.05（=忽略 sul），组合分数与
  sr_only 排序一致。
- 根因：sul 分数在该面板近乎退化——发布头部 sigmoid→softmax
  （`competitors/sul_bertgru.py:91`）使多数行输出同一常数值（0.26894…=σ(−1)），
  排名信息极弱。standalone AP 0.0213 即此症状。
- 结论：**Sul-BertGRU（当前移植形态）对 structure_ranker 无增量价值**；登记不实施。
  与 L6 既有判断一致（架构复杂度与收益无正相关证据）。

### W3-L8 校准/选择策略 —— 一阴一阳

- **阴性：MC-dropout σ 门控不成立**。真阳性系统性集中在**高不确定度区**：仅丢弃
  top-10% 高 σ 行即丢掉约 60% 真阳性（133/223），hits@50 从 24.7 掉到 11.3。
  即"顶部假阳性并非高不确定度位点"——与 kiae271 诊断的直观假说相反。
  σ 门控方向不立项。
- **阳性：分物种 isotonic 校准（验证集拟合、单调后处理）改善跨物种池化 Top-K**：
  hits@50 24.7→29.8（+21%），recall_full@200 0.497→0.581（+17%），不丢行。
  增益集中在 raw 最弱 seeds（0–5→9–19），强 seeds 持平，无一 seed 灾难性恶化。
- 含义：候选表为跨物种池化 Top-K；"分物种 isotonic 校准"列为**新发布候选变更**
  （仍须新发布号+新预注册，Top-K 变更纪律不变）。局限：文献轨 PU 面板（1:20）的
  池化数字不可直接外推全蛋白组候选表。

### W3-L2 ESM 分支（use_esm=1）——阴性，不立项

- 同口径重跑（仅 use_esm 0→1，其余 v11 超参冻结）：macro AP **0.1988±0.2006** vs
  v11 冻结 0.3863±0.1727；tomato 0.0129 vs 0.0185。10/10 seeds 劣于或接近 v11，
  且呈双峰（部分 seed 0.42–0.52、部分 0.03–0.05）——ESM 分支增加不稳定。
- 判定：**L2 不立项**。印证 v11 消融选择 use_esm=0；与登记册 L2 既有判断一致。
  与 L1/L4/L8 无关：L1（番茄结构）与 L8（分物种校准）的正面/待评估证据不受影响。

### W1 番茄 AlphaFold 结构覆盖（L1）——阴性（按现有管线形态），不立项

- 数据：29,767 番茄 accession → 27,542 下载 / 1,105 not_found / 0 错误；
  release v3 快照 30,681 行全量审计通过（SHA256 `60aa3af0…`），已注册
  `model_inputs.tsv`。番茄结构覆盖 0 → 90%（严格面板 7,803/8,676 行解掩蔽）。
- 文献轨（v11 冻结超参，仅换注册表）：macro **0.3863→0.2117**；arabidopsis
  0.4934→0.2302、rice 0.6470→0.3844、tomato 0.0185→0.0204（名义 +0.002，方差
  加倍，无增益证据）。
- 严格轨 5-fold：macro **0.5383→0.3257**；tomato **0.0413→0.0215（−48%）**。
- 机制线索（登记不实施）：① 结构 scaler 被番茄行主导（v3 后番茄占有结构行 ~78%），
  arabidopsis/rice 结构输入被重缩放——与两作物大幅退化一致；② 79 番茄阳性 × 7,803
  新激活结构行在 hidden-16 门控融合中弱监督过拟合。附带：magnaporthe（0 结构行）
  在 v3 臂反而上升，证实共享隐层被番茄结构行间接改变。
- 可检验下一步（新代设计候选）：分物种结构 scaler；番茄"仅覆盖率"臂；
  结构分支预训练/冻结。
- 注：v11 冻结模型的结构分支本就只对 ~3% arabidopsis/rice 行激活——"结构贡献显著"
  的既有印象需按此修正；L1 证据现为阴性。

### W4 分物种结构标定（诊断 B，MoE 设计轨）——W1 机制线索 #1 阳性，新代候选

- 复测 W1 机制线索①（"全局结构 scaler 被番茄行主导"）：文献轨 v11 冻结超参、
  同一 v3 番茄结构注册表、同一 10 seeds 面板（panel_sha256 gate 10/10），仅把
  结构特征从"全局标定"改为"分物种 z 标定"（scaler 只在 train present 行按物种
  拟合；未知/无结构物种恒等 fallback；后续全局二次 fit 因每物种 mean0/var1 近似
  恒等，不重新引入跨物种拉偏）。
- global 臂逐物种 test AP 与 W1 v3 逐 seed 一致（reproduction gate 10/10，
  容差 1e-6）：3-crop macro **0.2117（=W1 v3 精确复现）**；per_species 臂
  3-crop macro **0.3799**——恢复到 v2 基线 0.3863 的 98%。4-species macro
  0.2722→0.3899。
- 配对 Wilcoxon（10 seeds）：3-crop macro **p=0.0098**（9/10 seeds 正向，mean
  delta +0.1682）；4-species p=0.0195（8/10，+0.1177）。逐物种 test AP（10 seeds
  均值）：arab 0.2302→0.4957（9/10 up）、rice 0.3844→0.6265（8/10 up）、
  tomato 0.0204→0.0175（近地板，无灾难性恶化）、magnaporthe 0.4536→0.4200
  （0 结构行，训练端 spillover，非显著）。
- 结论：**W1 机制线索①实证成立**——全局结构 scaler 确被番茄主导；分物种标定把
  arabi/rice 恢复到 v2 水平且 tomato 不塌。命中 MoE 设计文档 §4B/§5 成功定义。
  **分物种结构标定列为新发布候选变更**（结构输入缩放=特征投影，非新架构）——
  实施仍需新发布号 + 新预注册 + 联合签署，本阶段不落地。
- 产物：`results/experiments/accuracy_campaign/w4_per_species_scaling/`
  （comparison.json + 每 seed 每臂 test/validation TSV）；claim_class
  `diagnostic_only`。

### 战役结论一览（2026-08-15）

- 已评估：L1（阴）、L2（阴）、L8-σ 门控（阴）、L8-分物种校准（**阳**）、
  W2 Sul-BertGRU 堆叠（阴）。
- 战役唯一正向杠杆：**分物种 isotonic 校准**（纯后处理，Top-K 池化规则），列为新发布
  候选变更；其余杠杆不立项。
- 战役后（MoE 设计轨）**W4 诊断 B 阳性**：分物种结构标定恢复 arabi/rice 到 v2 水平
  （3-crop p=0.0098），W1 机制线索 #1 实证成立，同样列为新发布候选变更（详见上文 W4 节）。
- 全部证据与报告：`results/experiments/accuracy_campaign/`（含 `DECISION_SUMMARY.md`）。


## 3. 评估纪律

- 任何杠杆评估必须：先在**冻结内部轨**（文献轨同口径）上跑通并记录数字，再决定是否立项；
  禁止直接改候选表或盲法流程。
- 评估结果无论正负都写入本登记册更新日志；阴性结果不删除。
- 新一代模型立项必须：新发布号 + 新 `fit_manifest` + 新 SAP/分析计划 + 双方重新签署。

## 4. 更新日志

- 2026-08-13：建立登记册。背景：Gate 0 冻结包完成（bundle `ab8a0353…`、候选表
  `36963c71…`、SAP/分析计划草案待签署）；用户询问"模型还能提升吗"，本文件作为
  "登记不实施"的正式答复存档。
- 2026-08-13：kiae271 模拟盲测（`evaluate_kiae271_with_release_bundle.py`）后更新
  第 2.1 节：0/99 位点进入 Top-2000 + 20 个干净位点 65% 中位数以上 → 新增 L8
  （分数校准/选择策略），L1（番茄结构覆盖）升为最高优先级；仍登记不实施。
- 2026-08-15：**共肽阴性注册 + 冻结 bundle 诊断**（增量建议 §1，最高价值项）——
  新增 `data/registry/copeptide_negatives_v1.tsv`（5 行，2 蛋白：BRG3 C206/C209/C212、
  RNF144b C122/C127），为项目注册**首个显式阴性证据**（AGENTS.md 允许注册的类别）：
  源自 kiad070 Fig.9B 肽段 SSCMICLPCR 与 Table S3 肽段 FYCPYKDCSAMLVNDSDEIVR，
  逐残基核对注册番茄参考蛋白组（唯一命中）。新增 `evidence/copeptide_negatives.py`
  三态分类（positive/negative/undetermined，"未定位≠未修饰"纪律：无
  site-determining-ion 覆盖时未修饰 Cys 标 undetermined 而非 negative）。
  **诊断结果（冻结 bundle 蛋白内排名，结构分支掩蔽）**：BRG3 阳性 C206(rank 8)/
  C212(rank 6) **未排到**阴性 C209(rank 5) 之上；RNF144b C122(rank 17) 未排到
  C127(rank 5) 之上。**阴性结论但信息量高**：±10 窗口序列特征在原理上无法分辨
  相距 3 残基的共肽 Cys，实证证实——支持"区分必须由结构承载"的 v2 方向；
  同时共肽阴性是唯一不受 Gate 2 条件 1（实验室独立性）限制的负样本轴。
  产物：`results/diagnostics/copeptide_negatives_v1.json`。
- 2026-08-15：**蛋白内排序诊断**（增量建议 §2）——训练用了 within-protein
  pairwise 损失但从未报告。12 个 mapped controls 全 Cys 冻结 bundle 打分：
  **SlWRKY6 C396 rank 1/7**（突变负担 1 vs 随机 4.0）、**PAD3 C440 rank 1/8**
  （1 vs 4.5）、5/12 进入 top-2；总负担 54 vs 随机 58.5。旗舰控制 SlWRKY6 在
  蛋白内排第 1 的同时不在全局 Top-2000 —— 任务口径匹配：模型训练目标就是蛋白内
  排序，不是蛋白组盲扫。**叙事修正**：工具适用于"已知蛋白的位点定位"；蛋白组盲扫
  结论需按此降级表述。WRKY71/ERF.D3 因 position_shift 状态未纳入（release 惯例）。
  产物：`results/known_controls/within_protein_ranking_v1.json`。
- 2026-08-15：**cys_density 体制假说诊断**（增量建议 §4）——关键模型事实：冻结
  模型唯一密度特征是**蛋白级** cys_density（蛋白内常数），无法驱动蛋白内排序；
  假设中的"同一特征需符号相反权重"在现有模型上不可实现（新特征需新发布）。
  诊断局部密度（±10 窗口 C 计数，不进模型）：SlWRKY6 真位点 0.048 < 干扰项 0.095
  （IDR 体制，方向符合假说）；BRG3 真位点 0.238 ≥ 蛋白均值 0.158（金属簇体制，
  方向符合假说）。两典型体制数据方向均与假说一致 → **体制路由（MoE）列入 v2
  架构候选证据 +1**；局部窗口密度作为新特征候选（新发布预注册）。
  产物：`results/diagnostics/regime_hypothesis_v1.json`。
- 2026-08-15：**共肽阴性系统化扩产（PXD024061，增量建议 §1 下一步）**——
  新增 `evidence/maxquant_sites_negatives.py`（MaxQuant sites 表解析：从
  `Sulfide(C) Probabilities` 列解析肽段序列+候选定位概率，按 mod-peptide ID
  分组合并；保守规则：候选概率 ≥0.75 → positive、<0.75 弱候选 → undetermined
  （"未定位≠未修饰"）、非候选 Cys 仅当修饰总数确定时 → negative；肽段在注册
  蛋白组唯一命中重新定位——MaxQuant 搜索坐标空间与注册蛋白组不一致）。
  `classify_peptide_cys` 增加 weak-candidate 维度（无行为回归）。扫描脚本
  `scripts/scan_copeptide_negatives_pxd024061.py`（幂等，site_id 去重）：
  PXD024061（Aroca et al. 2021, Antiox 10:508）Sulfide(C)Sites.txt（82 行）+
  CianoBiotin(C)Sites.txt（6 行）→ 74+6 肽段中 **7 个肽段注册、16 行新增**
  （5→21 行，6 蛋白 7 肽段，拟南芥；如 Q9FYD1 三 C 肽段 KPCFICGSLEHGAKQCSK
  = C191 修饰 + C194/C204 未修饰）。全部行通过拟南芥参考蛋白组 v2 逐残基
  核验（残基=C、肽段唯一命中、in-peptide 坐标一致）。**扩展诊断（9 肽段，
  冻结 bundle 蛋白内排名）**：8/9 未分离；唯一分离者 Q944L8 C204(rank 2) >
  C209(rank 3) 为噪声水平单例（9 肽段中 1 例在随机排序下非罕见），不改变
  "序列模型无共肽分离能力"结论；顺带 BRG3/RNF144b 原结论复现。
  产物：`results/diagnostics/copeptide_negatives_v1.json`（重跑覆盖）。
  证据轴价值：共肽阴性从 2 蛋白 2 肽段扩至 8 蛋白 9 肽段（跨番茄/拟南芥），
  v2 显式阴性池雏形。
- 2026-08-15：**P2 早期判定点：结构特征分离共肽位点（方法设计 §8 P2，首个
  阳性判定）**——新增 `evaluation/structure_features.py`（纯 numpy：PDB 解析、
  Shrake-Rupley SASA、Sγ 几何（8Å Cys 数/最近 Sγ 距离，AFDB 无金属故为配位
  代理）、6Å 正电残基计数、Sγ 库仑静电势近似、接触数；无 DSSP/APBS 依赖）。
  `scripts/evaluate_p2_structure_separation.py`：SlBRG3 + SlWRKY6 用注册 AFDB
  v6 结构（PDB/蛋白组逐残基一致校验）；PyMYB10 无 AFDB 文件（API 条目存在但
  模型文件全 404，折叠另行决策，如实记录未虚构）。
  **结果（决定性）**：① 逐特征方向——BRG3 上 3/7 特征在真位点 {206,212} vs
  金标准阴性 {209} 有干净方向（RSA 真位点更暴露、接触数更低、静电势更低）；
  WRKY6 上 RSA 与 Sγ 最近距离把 IDR 真位点 C396 排蛋白内第 2（负担 2 vs
  随机 4）。② 体制内判定——**RING 域（C197-C231）内带符号合成：C206 排
  rank 1/9、C209 排 rank 9/9（最后），负担 1 vs 随机 3.33**；全蛋白同符号
  合成反而 rank 9/10（N 端无序区 C24/C184/C185 暴露度更高，抢走方向）——
  **特征方向是体制内的，跨体制等权合成必失败**：命题三"路由而非平均"的首个
  结构侧实证。③ C209 是簇内埋藏最深（rsa 0.007、接触 26）、静电势最高
  （0.352）的 Cys = 最像配位 Cys；C206 是簇内最暴露（rsa 0.238）、静电势
  最低（0.159）的 Cys —— 结构侧支持设计文档 §1.3 RING 假设（9 Cys 必有
  游离 Cys，游离 Cys 才是 HS⁻ 靶点）。④ 监督树模型不可估（labeled 仅 4 个，
  21 Cys）——诚实记录；结构证据需体制内手工作合/域内 gate，树模型待 MIL
  级标签（L0）。产物：`results/diagnostics/p2_structure_separation_v1.json`。
  **gate 结论**：序列模型无法分离的共肽位点，结构特征在体制内可以
  （C206 rank 1 vs C209 rank 9）；P2 判定通过，结构承载方向获得首个正证据。
- 2026-08-15：**已知对照蛋白结构体制内检验（方法设计 §8 P2 的系统化扩展，
  n=12）**——新增 `evaluation/structure_regime.py`（pLDDT 体制分桶、蛋白内 z、
  按特征聚合排序/负担、体制局部带符号合成）+ `scripts/evaluate_structure_regime
  _separation.py` + 10 个单元测试。**新下载并注册 7 个拟南芥对照的 AFDB v6
  结构**（Q9FJI5/Q9LW27/F4K5T2/Q9FIJ0/Q940H6/A0MES8/Q8S929；PDB/蛋白组逐
  残基一致、SHA256 审计通过），使 12/12 mapped controls 全部有结构覆盖。
  产物：`results/diagnostics/structure_regime_separation_v1.json`。
  **关键结果（诊断级证据，n=12）**：① **aggregate 方向翻转 naive 暴露假说**
  ——折叠体制真位点系统性**更埋藏**：`contact_number_10a` 蛋白内 z +0.72
  （hit@2 6/9），全 n=12 上 mean_z +0.41、8/12 高于蛋白中位、负担 36 vs 随机
  58.5（38% 削减）；而 `rsa_relative` 负担 69 > 随机 58.5（暴露方向更差）。
  ② **方向随体制翻转**：disordered 体制（WRKY6/SlERF.D2/ABI4，n=3）真位点
  RSA 与最近 Sγ 距离 hit@2 全 3/3（暴露/孤立），folded 体制真位点是接触数
  高（埋藏/堆积）→ **单一全局符号复合无法同时服务两种体制**，命题三"路由
  而非平均"在 n=12 上的结构性实证。③ **体制局部排序无符号即全特征有效**：
  regime-local 排序（仅在真位点自身 pLDDT 桶内）对所有 7 特征负担都下降
  （最近 Sγ 59→45、RSA 69→53、接触 36→30）。④ 带符号复合上限：全局符号
  蛋白内负担 37、体制局部 32（随机 58.5/47.5），体制局部 top1 5/12、
  rank≤2 9/12。⑤ **BRG3 RING 连续性**：P2 level-5 方法（RING 局部 z + 暴露
  符号）复现 C206 rank 1、C209 rank 9、负担 1 vs 3.33；**显式记录张力**——
  同一 RING 在全局埋藏符号下 C209 升到 rank 7，说明 P2 的"游离暴露 Cys"
  解读是体制内的，不推广为一般金属簇规则（cluster-like 真位点 RNF144b C122
  完全埋藏、RSA 0.000、接触 z +1.60）。局限：contact_number 为主导是 7 特征
  集内的探索性选择（非预注册端点）；监督树仍不可估（labeled 12 个）。
  **v2 设计输入**：结构分支必须体制路由——folded 体制用埋藏/堆积（接触数），
  disordered 体制用暴露/孤立（RSA、最近 Sγ 距离）；RING/金属簇内另需簇内
  暴露 gate。
- 2026-08-15：**LOO + 排列置换检验（对 C 复合 ceiling 的诚实修正）**——新增
  `structure_regime.py` 的 `signs_from_aggregate`/`loo_composite_burden`/
  `permutation_null`（+5 单元测试，共 15）+ `scripts/evaluate_structure_regime
  _loo.py` → `results/diagnostics/structure_regime_loo_v1.json`（B=999、
  固定种子、确定性）。**修正结论**：① **contact_number_10a 是唯一显著特征**
  ——蛋白内负担 36 vs 排列零分布中位 59，**p=0.018**（regime-local 30，
  p=0.032）；埋藏/堆积是真信号。② 其余 6 特征全不显著（RSA p=0.84/0.85，
  暴露假说在已知对照上**无超出随机信号**）。③ 7 特征等权带符号复合**不能留一
  泛化**（LOO 负担 39/38 vs 随机 58.5，p≈0.3）——§9.2 的复合 ceiling 与
  "方向随体制翻转"需降级；disordered 侧暴露信号（n=3）排列下不显著。
  **v2 设计输入修正**：结构分支以 contact_number 为首要特征，其余特征待 MIL
  级标签（L0）学权重/特征选择；体制路由方向保留但暴露分支降为假设。
- 2026-08-15：**共肽结构分离检验（把 C 的已验证特征显式应用到共肽阴性池）**
  ——新增 `copeptide_structure.py`（肽段内正/负对比 + 保组成排列零分布，
  +12 单元测试）+ `scripts/evaluate_copeptide_structure_separation.py` →
  `results/diagnostics/copeptide_structure_separation_v1.json`。9 共肽组
  （2 番茄金标准 + 7 拟南芥，全部 AFDB v6 结构，PDB/蛋白组一致审计通过）。
  **零结果**：7 特征无一在任一方向分离共肽正/负（contact_number p=0.759/
  0.234，其余全 p>0.3）；方向混杂——BRG3 负 C209 最埋藏（26 vs 22/16），
  RNF144B 反而正 C122 更埋藏（24 vs 14）。**解读**：共肽阴性是结构与序列共同
  的"不可分辨区"（序列 8/9 败、结构 0/9 显著），验证了共肽阴性设计价值
  （检测匹配时连结构也不能假称分离）；不推翻跨蛋白埋藏信号（contact
  p=0.018 是"相对同蛋白其他 Cys"，簇内分辨是另一层面）；P2"游离暴露 Cys"
  是簇内 case-level 读数、不推广。**v2 输入**：结构分支声称范围限定在蛋白内
  排序，不扩展到共肽内正/负分辨；共肽池保持为有效性约束与标定资源。
- 2026-08-15：**PXD072089 水稻共肽阴性扩产（增量建议 §1 下一步，第 4 物种）**
  ——`evidence/pnas_site_negatives.py`（纯逻辑 + TypedDict，坐标偏移安全的
  秩映射）+ `scripts/scan_copeptide_negatives_pxd072089.py`（pandas I/O +
  水稻蛋白组唯一重定位 + registry 写入）+ 19 unit + 2 scientific 测试。
  5 个部分修饰肽段 → 8 蛋白-肽段组 → **19 行**（10 阳性 + 9 阴性），
  `copeptide_negatives_v1.tsv` 21 → 40 行。PNAS 2025 独立实验室（Gate 2
  条件 1 实验室独立性在阴性维度的补充）；坐标安全映射处理 Q5ZCB1 的 +1 偏移。
  这批行此前由一次性脚本生成、未提交；本次可复现代码**精确复现**（幂等
  `rows_added=0`）。逐残基核验水稻参考蛋白组 v1。
- 2026-08-15：**冻结 bundle 共肽诊断重跑（40 行 / 17 蛋白 / 3 物种）+ 置换检验**
  ——`evaluate_copeptide_negatives.py` 加 rice 蛋白组并新增保组成置换零分布
  （B=999）。**17 蛋白-肽段实例中 3/17 分离，随机期望 7.6（零分布中位 8、
  p5=4/p95=11）；p(分离能力超随机)=0.997、观测位于左尾**。**跨番茄/拟南芥/
  水稻三物种、17 实例证实冻结序列模型对共肽正/负零分离能力**——"序列无法区分
  共肽位点"从 2 蛋白 2 肽段升级为 17 实例 3 物种 + 置换检验口径。
- 2026-08-15：**水稻共肽结构重测（结构零结果跨物种稳健，9 → 14 组）**——
  `evaluate_copeptide_structure_separation.py` 加 rice 蛋白组、按**不同肽段**
  去重（IIPTPNC 4 旁系同源实例只留首个，避免 4 倍加权）、缺 AFDB 结构跳过并
  诚实计数（本次 0 dropped）。P0C361/Q5ZCB1 用批量下载器补下（2/2 成功，
  PDB/蛋白组长度+残基校验通过）。**14 组（2 番茄 + 7 拟南芥 + 5 水稻）7 特征
  无一在任一方向分离共肽正/负（全 p>0.05，最小 coulomb p=0.39）；方向依旧
  混杂**——水稻 5 组 P0C361/A0A0P0Y2A9/A0A0P0WV74 pos_higher、A0A0P0Y2K3
  pos_lower、Q5ZCB1 5-Cys 肽段 mixed，与番茄两金标准相反方向一致。**结论升级**：
  共肽"结构与序列共同不可分辨区"的零结果从 2 物种 9 肽段稳健扩展到 **3 物种
  14 肽段**（独立实验室/独立富集方案）；v2 声称范围（结构分支限蛋白内排序、
  不扩共肽内分辨）跨物种证据更足。
- 2026-08-15：**体制路由 MoE 诊断（L6 架构候选评估，诊断 A 阴性）**——
  `evaluation/regime_moe.py`（专家方向写死为结构先验、非 labeled 学习，无 LOO
  学习偏差；+6 单元测试）+ `scripts/evaluate_regime_moe.py` →
  `results/diagnostics/regime_moe_v1.json`（n=12 已知对照，B=999 置换）。
  设计文档 `2026-08-15-regime-routing-moe-design.md`（命题三 → 预注册草案 +
  三轨诊断）。**结果：MoE 路由不优于 contact 单特征**——routed regime 31 vs
  contact regime 30（p 0.044→0.057，不再显著）；routed protein 38 vs 36。
  folded 专家 contact 方向承载全部信号（9 folded 对照 6 个 rank≤2）；
  disordered 专家 RSA/Sγ 方向无增益且稀释（A0MES8 rank 3/3，n=3）——与 §9.4
  "RSA 排列下不显著"一致。**结论：命题三方向正确，但当前特征分辨率下
  disordered 专家无物可路由；MoE 架构在 labeled n=12 + 7 特征下不立项。**
  分物种/分体制结构标定（W1 scaler 线索）与路由正交，列新代候选（Gate 3 后
  再评估）。
- 2026-08-15：**诊断 B（分物种结构标定）——W1 机制线索 #1 阳性，新代候选**——
  `evaluation/species_structure_scaling.py`（纯函数：fit_species_struct_scalers/
  transform_species_struct，每物种在 train present 行 fit TrainOnlyScaler、masked
  与无 scaler 物种恒等 fallback；+5 单元测试）+ `scripts/accuracy_campaign/
  eval_per_species_structure_scaling.py`（复用 W1 runner 面板构建与逐位点路径，
  两臂仅结构输入投影不同：global=冻结全局标定 / per_species=分物种 z 标定；
  +4 脚本 helper 单元测试）→ `results/experiments/accuracy_campaign/
  w4_per_species_scaling/comparison.json`（claim_class `diagnostic_only`）。
  **口径**：文献轨 v11 冻结超参、同一 v3 番茄结构注册表、同一 10 seeds 面板
  （panel_sha256 gate 10/10），test 分区 AP@200。global 臂逐物种 test AP 与 W1 v3
  逐 seed 一致（reproduction gate 10/10，容差 1e-6）→ 3-crop macro 0.2117 精确复现
  W1 v3。per_species 臂 3-crop macro **0.3799**（v2 基线 0.3863 的 98%）。
  **结果**：配对 Wilcoxon 3-crop macro **p=0.0098**（9/10 seeds 正向，mean delta
  +0.1682）；4-species p=0.0195（8/10）。逐物种 test AP（10 seeds 均值）：
  arab 0.2302→0.4957（+0.265，9/10 up）、rice 0.3844→0.6265（+0.242，8/10 up）、
  tomato 0.0204→0.0175（近地板，无灾难性恶化）、magnaporthe 0.4536→0.4200
  （0 结构行，训练端 spillover，非显著）。**结论：W1 机制线索①"全局结构 scaler 被
  番茄行主导"实证成立**；分物种标定把 arabi/rice 恢复到 v2 水平且 tomato 不塌，
  命中 MoE 设计文档 §4B/§5 成功定义。**分物种结构标定列为新发布候选变更**（特征
  投影，非新架构）——实施须新发布号 + 新预注册 + 联合签署，本阶段不落地。
