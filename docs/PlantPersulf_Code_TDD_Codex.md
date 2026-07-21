# PlantPersulf-Code 技术设计文档与 TDD 实施计划

> **用途**：本文件可直接交给 Codex 执行。  
> **项目定位**：利用真实、可追查的公共蛋白质组、翻译后修饰组、转录组和结构数据，开发一个结构与生物学上下文感知的植物蛋白硫巯基化候选发现框架，并形成可交付给张华老师团队的番茄候选验证包。  
> **首阶段目标**：做出可信的计算证据与可复现 benchmark，不提前声称发现了新的生物学机制。  
> **开发范式**：严格 Test-Driven Development（RED → GREEN → REFACTOR），不允许先写生产代码再补测试。

---

## 0. Codex 首条执行指令

将下面整段作为 Codex 新会话的第一条指令：

```text
你现在负责执行 PlantPersulf-Code 项目。

请完整阅读仓库根目录的：
1. docs/PlantPersulf_Code_TDD.md
2. AGENTS.md
3. configs/data_sources.yaml
4. configs/scientific_integrity.yaml

这是一个真实科研项目，必须严格遵守以下不可变约束：

A. 所有用于科学分析、模型训练、验证、图表和结论的数据，都必须来自真实、公开、可追查的数据集，或者后续合作方正式提供且带有样本表和实验元数据的数据。
B. 禁止使用合成数据、随机生成的生物学数值、模拟标签、伪造样本、占位结果或手工编造的模型指标。
C. 单元测试不得发明生物学记录。需要小样本时，只能从真实数据中提取“真实微型夹具”，并记录来源 accession、原始文件、行号/谱图/肽段标识和 SHA256。
D. 临时目录、空文件、非法配置等纯软件错误测试可以人工构造，但不得被用于任何科研结果或性能结论。
E. 所有下载文件必须进入 provenance 清单，记录 accession、官方来源、下载时间、URL、文件大小、SHA256、许可证/使用说明、解析脚本版本。
F. “未检测到硫巯基化”不等于负样本。主任务必须按正样本—未标记（PU learning）设计，除非存在明确实验阴性证据。
G. 任何图表和结论必须能够由一条命令从已登记的真实数据重新生成。
H. 每个功能必须先写失败测试，确认失败原因正确，再写最小实现使其通过。
I. 每完成一个任务，运行目标测试、全量快速测试和科学完整性检查，然后提交一次独立 commit。
J. 禁止在 README、论文草稿或报告中使用“显著优于”“SOTA”“发现新机制”等表述，除非对应统计检验、外部验证和来源表已经通过 release gate。

不要一次性实现整个项目。先执行 Task 0，只完成 Task 0 并报告：
- 创建/修改的文件
- RED 测试输出
- GREEN 测试输出
- commit SHA
- 尚未满足的条件

Task 0 通过审查后再执行 Task 1。
```

---

# 1. 项目目标与边界

## 1.1 科学目标

构建一个可解释、可外部验证的计算框架，回答：

1. 哪些植物蛋白半胱氨酸位点具有较高的实验性硫巯基化支持？
2. 哪些序列、结构和生物学上下文特征与已检测的硫巯基化位点相关？
3. 模型能否在**未见蛋白家族、未见研究、未见物种或较晚发表的数据集**中恢复真实实验阳性位点？
4. 在番茄成熟过程中，哪些候选位点同时具备：
   - 跨研究模型支持；
   - 蛋白表达或转录动态支持；
   - 与磷酸化等 PTM 串扰的真实数据支持；
   - H₂S/乙烯/成熟机制背景支持；
   - 明确的实验可操作性？
5. 能否形成一份冻结、可追踪、不过度解释的候选验证包，供张华老师团队开展盲法湿实验？

## 1.2 第一阶段明确不做

首阶段不得：

- 宣称模型证明某个位点具有因果作用；
- 把未检测到的半胱氨酸直接标成生物学阴性；
- 使用随机生成的氨基酸序列、随机标签或模拟蛋白组作为科研验证；
- 把 AlphaFold 预测结构当作实验结构；
- 把模型注意力直接解释为因果机制；
- 在没有真实实验验证时使用“regulatory code 已被证实”的表述；
- 为了提升指标改变预先冻结的数据划分；
- 在测试集上调参；
- 将张华老师论文中的已知阳性位点同时放入训练集和独立验证集；
- 把公开论文补充表手工抄写后当作未经审计的真值；
- 使用任何无法回溯到 accession、论文补充材料或正式合作数据清单的生物学输入。

## 1.3 第一阶段交付物

1. 真实数据注册表和下载审计系统；
2. 植物硫巯基化正样本—未标记 benchmark v1；
3. 防同源泄漏、leave-study-out、leave-species-out 和时间切分；
4. 序列基线、结构增强基线、PU 学习基线；
5. 严格外部验证报告；
6. 番茄成熟上下文数据整合；
7. 候选位点证据卡；
8. `candidate_release_v1` 冻结包；
9. 面向张华老师的两页合作摘要和湿实验优先级矩阵；
10. 完整代码、测试、配置、日志、来源表和可重现命令。

---

## 1.4 合作科学定位与项目增量

基于已公开论文形成的项目工作判断如下：

- 张华老师团队已经具备 LC-MS/MS 位点鉴定、硫巯基化检测、Cys 定点突变、CRISPR/过表达、蛋白互作、转录调控、PTM 串扰和番茄表型验证能力；
- 本项目不以代做 RNA-seq、常规差异分析、富集分析或通用生信服务作为合作价值；
- 本项目拟补充可追踪位点级证据整合、防泄漏跨研究预测、训练前候选冻结、匹配对照和前瞻性盲法验证；
- 单独增加 AI 模型、随机交叉验证或候选排名表，不构成 Nature-family 级生物学创新；
- 联合研究要从逐个单蛋白机制上升到可推广、可预测并经前瞻实验检验的植物硫巯基化规律。

上述判断只用于项目设计，不评价团队未公开能力，也不保证任何期刊结果。

---

# 2. 科学完整性合同

## 2.1 数据等级

每条输入记录必须标记以下等级之一：

| 等级 | 定义 | 是否可用于训练 | 是否可用于主要结论 |
|---|---|---:|---:|
| A | 官方公共数据库中的原始数据，具 accession 和文件校验值 | 是 | 是 |
| B | 官方公共数据库中的处理后数据，具 accession 和方法说明 | 是 | 是，但需说明处理来源 |
| C | 同行评议论文的正式补充数据，具 DOI、文件名和校验值 | 可用于阳性控制/外部验证 | 可以，需标明非原始数据 |
| D | 合作方提供的真实实验数据，具样本表、实验方案、原始文件和数据使用确认 | 是 | 是 |
| E | 文献正文中单个已验证位点或机制 | 不作为大规模训练主数据 | 仅作机制卡与阳性控制 |
| X | 无 accession、无文件、无法追踪来源或来源冲突 | 否 | 否 |

## 2.2 真实微型夹具

自动化测试需要快速运行时，只能使用从真实数据截取的微型夹具：

```text
tests/fixtures/real/
├── PXD006140/
│   ├── source_manifest.json
│   ├── peptides_subset.tsv
│   └── SHA256SUMS
├── PXD051570/
│   ├── source_manifest.json
│   ├── phosphosites_subset.tsv
│   └── SHA256SUMS
└── GSE163745/
    ├── source_manifest.json
    ├── expression_subset.tsv
    └── SHA256SUMS
```

`source_manifest.json` 必须包括：

```json
{
  "accession": "PXD006140",
  "source_file": "官方文件名",
  "source_sha256": "原始文件校验值",
  "extraction_command": "生成微型夹具的精确命令",
  "row_selector": "精确筛选条件",
  "fixture_sha256": "微型夹具校验值",
  "scientific_use": "tests_only",
  "biological_values_modified": false
}
```

禁止人为修改真实微型夹具中的蛋白、肽段、位点、表达量或标签。

## 2.3 结论来源链

每个最终候选必须可追踪到：

```text
candidate_id
→ model_release_id
→ split_id
→ feature_version
→ protein sequence accession/version
→ structure accession/version
→ source study accession
→ source file SHA256
→ parser commit
→ statistical result file
→ report figure/table
```

缺少任一关键环节，候选不得进入冻结发布包。

---

# 3. 真实公共数据注册表 v1

下面的数据是启动注册表。Codex 必须通过官方 API/记录重新获取元数据，不得把本表文字当作下载结果。

## 3.1 植物硫巯基化核心数据

| Accession | 物种/场景 | 用途 | 首阶段角色 |
|---|---|---|---|
| PXD006140 | Arabidopsis；WT 与 L-CYSTEINE DESULFHYDRASE 1 相关材料；tag-switch 硫巯基化蛋白质组 | 构建实验阳性集合 | 主要发现数据 |
| PXD024061 | Arabidopsis；氮饥饿根与碳饥饿叶；定量硫巯基化蛋白质组 | 条件特异性与外部研究验证 | 外部/leave-study-out |
| PXD035795 | Arabidopsis；非光呼吸与空气转换条件下硫化物调控 | 条件迁移验证 | 外部研究验证 |
| PXD039999 | Arabidopsis；干旱及 H₂S 预处理；硫巯基化研究 | 胁迫场景验证 | 外部研究验证 |

要求：

- 优先下载原始文件、搜索结果文件、样本设计文件和 SDRF；
- 若某项目没有足够的位点级结果，必须标记为 `protein_level_only`，不得伪造位点；
- 不同实验使用的富集、化学标记和搜索方法必须作为 study context 特征记录；
- 同一蛋白在不同研究中重复出现必须保留 study 维度，不能无条件去重。

## 3.2 番茄成熟真实上下文数据

| Accession | 数据 | 用途 |
|---|---|---|
| PXD051570 | 番茄五个成熟阶段的蛋白质组和磷酸化组，包含大规模蛋白与磷酸化位点定量 | 构建番茄成熟/PTM 上下文 |
| GSE163745 | 番茄 WT、单突变和双突变成熟调控转录组 | 成熟状态和经典调控因子验证 |
| GSE142713 | SlJMJ6 过表达与 WT 的 RNA-seq | 成熟相关表达网络 |
| GSE142712 | SlJMJ6 过表达与 WT 的 H3K27me3 ChIP-seq | 表观调控上下文 |
| GSE267238 | 番茄 WT 与成熟突变体在 mature green/breaker 阶段的 RNA-seq | 独立成熟状态验证 |
| GSE49289 | WT、FUL1/FUL2 抑制材料和 rin 突变体 RNA-seq | 经典成熟网络补充 |
| GSE40257 | RIN 靶基因 ChIP 数据 | 直接调控证据 |
| GSE125306 | 番茄果实 m6A-seq，包含 WT 阶段与 Cnr 突变体 | 可选 PTM/表观转录组上下文，不进入 MVP 必选模块 |

要求：

- MVP 必选：PXD051570、GSE163745；
- 第二批必选：GSE142713、GSE142712、GSE267238；
- 其余数据只有在完成前述数据质量审计后才能启用；
- 每个 GEO/SRA 数据集必须记录 BioProject、BioSample、Run accession 和样本条件；
- 若使用作者提供的 processed matrix，必须同时保存其处理说明；
- 不同 ITAG/SL 基因组版本必须通过版本化映射表统一，禁止静默转换。

## 3.3 张华老师团队已发表机制作为训练外阳性控制

已发表机制只能作为 `known_positive_mechanism_card`，不能用于声称新发现，也不能用于解除训练 benchmark 的数据不足。

| 控制 | DOI | 控制类型 | 使用限制 |
|---|---|---|---|
| SlWRKY6 Cys396 | `10.1093/plphys/kiae271` | `strong_single_site_control` | 训练外单点回顾性恢复 |
| SlERF.D2 Cys35 | `10.1111/tpj.70000` | `strong_single_site_control` | 训练外单点回顾性恢复 |
| BRG3 Cys206/Cys212 | `10.1093/plphys/kiad070` | `conditional_site_group_control` | 无单点拆分证据时不得解释为两个独立功能阳性 |
| ERF.D3 Cys115/Cys118 | `10.1093/plphys/kiae560` | `conditional_site_group_control` | 无单点拆分证据时不得解释为两个独立功能阳性 |

磷酸化、泛素化和转录调控位点只能作为 PTM 串扰解释证据。2026 年 SlWRKY6–SlGRF1–SlGIF2 研究（`10.1093/plphys/kiag512`）复用同一 SlWRKY6/H₂S 机制链，不得重复计为新的 `independent_validation_unit`。

处理规则：

1. 登记正文或补充材料、DOI、官方 URL、文件名和 SHA256；
2. 登记 canonical protein accession/version、位点映射、证据层级和冲突；
3. 人工复核并记录复核状态；
4. 登记 `control_type`、`mechanism_lineage_id` 和 `independent_validation_unit`；
5. 所有控制标记 `positive_control_only`、`retrospective_recovery_test` 和 `excluded_from_training`；
6. 控制身份、独立性分组和评价规则在训练前冻结；
7. 控制排名不得用于特征、模型、阈值、超参数或停止条件选择；
8. 同一蛋白、位点或机制谱系的后续论文不得重复计数。

## 3.4 参考序列和结构

- Arabidopsis 与 Solanum lycopersicum 的 UniProt reference proteome；
- 对应版本的 NCBI/Ensembl Plants/Phytozome 基因注释，用于 ID 映射；
- AlphaFold Protein Structure Database 中对应模型；
- 每个序列和结构文件记录版本、下载时间和 SHA256；
- 无结构或映射冲突的蛋白不得被静默删除，必须进入缺失报告。

---

# 4. 研究设计

## 4.1 任务定义

主任务是**正样本—未标记位点排序**，而不是普通二分类。

### 实验阳性

满足以下条件之一：

1. 官方数据处理结果中明确定位到含目标 Cys 的 persulfidated peptide；
2. 论文正式补充材料中通过 LC-MS/MS、位点突变或化学方法确认；
3. 合作方后续真实实验明确验证。

### 未标记

- 同一真实蛋白中的其他 Cys；
- 同一物种蛋白组中的其他 Cys；
- 未在当前实验检测到的位点。

未标记不等于阴性。

## 4.2 分析层级

分别报告：

1. 肽段层级；
2. Cys 位点层级；
3. 蛋白层级；
4. 研究/条件层级；
5. 物种层级。

禁止将蛋白级结果包装为位点级精度。

## 4.3 数据切分

必须预先生成并冻结以下 split：

### Split A：蛋白家族聚类切分

- 使用 MMseqs2 对蛋白序列聚类；
- 默认阈值：30% sequence identity、至少 50% coverage；
- 同一 cluster 不得跨 train/validation/test；
- 阈值敏感性：20%、30%、40%。

### Split B：leave-study-out

每次完全留出一个 PXD 研究，训练过程不得读取其标签和归一化统计量。

### Split C：leave-species-out

训练以 Arabidopsis 数据为主，番茄作为迁移/候选场景；如果后续加入其他植物物种，必须执行物种留一。

### Split D：时间切分

- 从官方论文元数据获取发布日期；
- 在代码中固定 cutoff；
- 早期研究用于训练；
- 较晚公开研究用于外部验证；
- cutoff 必须在任何模型训练前写入 `configs/splits/time_split_v1.yaml` 并提交。

### Split E：已知机制排除切分

张华老师团队已发表机制位点全部从训练集剔除，专门用于 retrospective recovery。

## 4.4 模型层次

### Baseline 0：无学习基线

- Cys 局部氨基酸频率；
- 溶剂可及性排序；
- 蛋白丰度排序；
- 随机但固定的 matched ranking 仅用于统计参考。

注意：随机参考不得作为科研数据；只可用于计算 enrichment 的零假设，必须固定种子并明确标记为统计随机化，不可生成生物学记录。

### Baseline 1：传统机器学习

- Logistic Regression / PU Logistic；
- XGBoost；
- Random Forest；
- 特征包括局部序列、理化性质、保守性、结构可及性。

### Baseline 2：蛋白语言模型

- 冻结 ESM-2 residue embedding；
- 不在小样本阶段端到端微调大型模型；
- 训练轻量分类/排序头；
- 按 protein cluster 严格划分。

### Model 1：Structure-aware PU ranker

输入：

- Cys 周围序列窗口；
- residue-level PLM embedding；
- pLDDT；
- relative solvent accessibility；
- 二级结构；
- disorder；
- 蛋白长度和 Cys 密度；
- study context。

输出：

- 位点得分；
- 不确定性；
- applicability domain；
- 证据贡献。

### Model 2：Context router

动态路由以下证据：

1. sequence；
2. structure；
3. protein abundance；
4. transcript dynamics；
5. phosphosite proximity/dynamics；
6. pathway/mechanism context。

MVP 只有 Model 1 达到 release gate 后才允许实现 Model 2。

## 4.5 评价指标

由于是 PU 问题，主要指标不能只写普通 Accuracy。

必报：

- held-out positive Recall@K；
- Mean Reciprocal Rank；
- Average Precision（明确 positive-vs-unlabeled 假设）；
- matched-cysteine enrichment；
- cluster-level bootstrap 95% CI；
- leave-study-out recovery；
- retrospective known-mechanism recovery；
- calibration / reliability；
- applicability-domain 分层结果；
- 每个数据集单独结果，不只报告 pooled 结果。

统计比较：

- paired cluster bootstrap；
- permutation test；
- 多重比较校正；
- 效应量与置信区间；
- 不以单次 seed 最优结果作为主结果。

---

# 5. 仓库架构

```text
PlantPersulf-Code/
├── AGENTS.md
├── README.md
├── pyproject.toml
├── environment.yml
├── Makefile
├── configs/
│   ├── scientific_integrity.yaml
│   ├── data_sources.yaml
│   ├── id_mapping.yaml
│   ├── features/
│   │   ├── sequence_v1.yaml
│   │   └── structure_v1.yaml
│   ├── splits/
│   │   ├── cluster_split_v1.yaml
│   │   ├── study_split_v1.yaml
│   │   ├── species_split_v1.yaml
│   │   ├── time_split_v1.yaml
│   │   └── known_mechanism_holdout_v1.yaml
│   └── experiments/
│       ├── baseline_sequence_v1.yaml
│       ├── baseline_structure_v1.yaml
│       └── pu_ranker_v1.yaml
├── data/
│   ├── registry/
│   │   ├── datasets.tsv
│   │   ├── files.tsv
│   │   ├── samples.tsv
│   │   ├── publications.tsv
│   │   └── SHA256SUMS
│   ├── raw/                 # gitignored
│   ├── external/            # gitignored
│   ├── interim/             # gitignored
│   └── processed/           # gitignored
├── src/plantpersulf/
│   ├── provenance/
│   │   ├── schema.py
│   │   ├── registry.py
│   │   ├── hashing.py
│   │   └── audit.py
│   ├── download/
│   │   ├── pride.py
│   │   ├── geo.py
│   │   ├── sra.py
│   │   ├── uniprot.py
│   │   └── alphafold.py
│   ├── proteomics/
│   │   ├── metadata.py
│   │   ├── peptide_parser.py
│   │   ├── site_normalizer.py
│   │   └── study_context.py
│   ├── transcriptomics/
│   │   ├── metadata.py
│   │   ├── expression.py
│   │   └── gene_id_mapping.py
│   ├── benchmark/
│   │   ├── labels.py
│   │   ├── unlabeled.py
│   │   ├── cluster_split.py
│   │   ├── study_split.py
│   │   ├── time_split.py
│   │   └── leakage_audit.py
│   ├── features/
│   │   ├── sequence.py
│   │   ├── esm2.py
│   │   ├── structure.py
│   │   ├── conservation.py
│   │   └── context.py
│   ├── models/
│   │   ├── baselines.py
│   │   ├── pu_risk.py
│   │   ├── structure_ranker.py
│   │   ├── router.py
│   │   └── calibration.py
│   ├── evaluation/
│   │   ├── metrics.py
│   │   ├── bootstrap.py
│   │   ├── permutation.py
│   │   └── external_validation.py
│   ├── tomato/
│   │   ├── ripening_context.py
│   │   ├── phospho_context.py
│   │   ├── mechanism_cards.py
│   │   ├── candidate_score.py
│   │   └── evidence_cards.py
│   ├── reporting/
│   │   ├── tables.py
│   │   ├── figures.py
│   │   ├── model_card.py
│   │   └── candidate_release.py
│   └── cli.py
├── tests/
│   ├── fixtures/real/
│   ├── unit/
│   ├── integration/
│   ├── scientific/
│   └── release/
├── scripts/
│   ├── fetch_registered_data.py
│   ├── build_real_fixtures.py
│   ├── build_benchmark.py
│   ├── extract_features.py
│   ├── run_experiment.py
│   ├── validate_external.py
│   ├── rank_tomato_candidates.py
│   └── build_candidate_release.py
├── docs/
│   ├── PlantPersulf_Code_TDD.md
│   ├── data_dictionary.md
│   ├── scientific_claims_policy.md
│   ├── zhang_collaboration_brief.md
│   └── wetlab_validation_matrix.md
├── results/                 # gitignored, release snapshot除外
├── releases/
│   └── candidate_release_v1/
└── manuscripts/
    └── computational_pilot/
```

---

# 6. 全局工程约束

- Python 3.10 或 3.11，环境版本锁定；
- 使用 `pytest`、`ruff`、`mypy`；
- 配置驱动，不在代码中硬编码 accession 和路径；
- 所有随机操作要求显式 seed；
- 所有科学输出包含 config hash、git commit 和数据 registry hash；
- 原始大文件不提交 Git；
- release 中只保存可公开的来源表、结果表、图、配置和下载脚本；
- 所有外部 API 请求带重试、超时和响应缓存；
- 下载失败必须失败退出，不得自动生成替代文件；
- 解析不到字段必须输出缺失报告，不得猜测；
- ID 映射一对多时保留冲突表，不得任意选择第一项；
- 每个实验至少 5 个固定 seed，主结果报告均值、置信区间和全部 seed；
- 不允许把 validation set 用于最终外部测试结论；
- 不允许根据测试结果修改 split；
- 不允许在代码中存在 `synthetic`, `fake`, `dummy`, `mock_biology` 科研路径；
- CI 中扫描上述禁用词及未登记数据路径；
- 纯软件 mock 只能用于 API 异常和文件系统错误，不可用于科学测试。

---

# 7. TDD 实施任务

## Task 0：仓库初始化与科研诚信闸门

**目标**：任何科学模块开始前，先建立禁止伪造数据和来源不明数据的自动闸门。

**创建文件**

- `AGENTS.md`
- `pyproject.toml`
- `configs/scientific_integrity.yaml`
- `src/plantpersulf/provenance/schema.py`
- `src/plantpersulf/provenance/audit.py`
- `tests/unit/test_integrity_config.py`
- `tests/scientific/test_no_unregistered_scientific_inputs.py`
- `tests/release/test_no_synthetic_scientific_artifacts.py`

### RED 1：缺少来源字段必须失败

```python
def test_dataset_record_requires_traceable_source():
    from plantpersulf.provenance.schema import DatasetRecord

    with pytest.raises(ValueError, match="accession"):
        DatasetRecord(
            accession="",
            repository="PRIDE",
            source_url="",
            scientific_role="training",
        )
```

运行：

```bash
pytest tests/unit/test_integrity_config.py::test_dataset_record_requires_traceable_source -v
```

预期：因模块或类不存在而 FAIL。

### GREEN 1

实现最小 `DatasetRecord`，强制验证：

- accession；
- repository；
- source_url；
- scientific_role；
- metadata_retrieved_at；
- metadata_sha256。

### RED 2：未登记文件不能进入分析

```python
def test_unregistered_file_is_rejected(tmp_path):
    from plantpersulf.provenance.audit import assert_registered_input

    unknown = tmp_path / "unknown.tsv"
    unknown.write_text("protein\tvalue\nP1\t1\n", encoding="utf-8")

    with pytest.raises(RuntimeError, match="not registered"):
        assert_registered_input(unknown, registry_path=tmp_path / "files.tsv")
```

### GREEN 2

实现严格 registry 检查。registry 不存在、文件未登记或 SHA256 不一致都必须失败。

### RED 3：release 中出现合成科研产物必须失败

测试创建文件名为 `synthetic_results.tsv`，release audit 必须拒绝。

注意：这是纯软件策略测试，不包含任何合成生物学记录。

### 验收

```bash
pytest tests/unit/test_integrity_config.py \
       tests/scientific/test_no_unregistered_scientific_inputs.py \
       tests/release/test_no_synthetic_scientific_artifacts.py -v
ruff check .
mypy src/plantpersulf
```

全部通过后提交：

```bash
git add .
git commit -m "chore: establish scientific integrity and provenance gates"
```

---

## Task 1：真实数据注册表与官方元数据解析

**目标**：建立 accession 驱动的数据源配置，验证项目存在、物种、标题、发布日期和文件清单。

**创建文件**

- `configs/data_sources.yaml`
- `src/plantpersulf/download/pride.py`
- `src/plantpersulf/download/geo.py`
- `src/plantpersulf/provenance/registry.py`
- `tests/integration/test_pride_registry_real.py`
- `tests/integration/test_geo_registry_real.py`

`configs/data_sources.yaml` 至少登记：

```yaml
pride:
  - accession: PXD006140
    role: discovery_persulfidation
  - accession: PXD024061
    role: external_persulfidation
  - accession: PXD035795
    role: external_persulfidation
  - accession: PXD039999
    role: external_persulfidation
  - accession: PXD051570
    role: tomato_ripening_proteome_phosphoproteome

geo:
  - accession: GSE163745
    role: tomato_ripening_transcriptome
  - accession: GSE142713
    role: tomato_ripening_transcriptome
  - accession: GSE142712
    role: tomato_ripening_chipseq
  - accession: GSE267238
    role: tomato_ripening_transcriptome
```

### RED

真实 API 集成测试：

```python
@pytest.mark.network
def test_pxd006140_metadata_matches_expected_project():
    project = PrideClient().get_project("PXD006140")
    assert project.accession == "PXD006140"
    assert "persulfidation" in project.title.lower()
    assert project.files
```

```python
@pytest.mark.network
def test_gse163745_metadata_is_tomato_ripening_study():
    series = GeoClient().get_series("GSE163745")
    assert series.accession == "GSE163745"
    assert series.organism == "Solanum lycopersicum"
    assert len(series.samples) > 0
```

### GREEN

实现 API 客户端和本地 JSON 缓存。缓存文件本身也要登记 SHA256。

### 验收

```bash
pytest -m network tests/integration/test_pride_registry_real.py \
                  tests/integration/test_geo_registry_real.py -v
python scripts/fetch_registered_data.py --metadata-only
python -m plantpersulf.cli audit-registry
```

输出：

- `data/registry/datasets.tsv`
- `data/registry/files.tsv`
- `data/registry/samples.tsv`
- `data/registry/publications.tsv`

提交：

```bash
git add configs src tests data/registry
git commit -m "feat: add official PRIDE and GEO dataset registry"
```

---

## Task 2：真实文件下载、校验与失败保护

**目标**：下载真实文件并确保完整性，禁止静默替代。

**创建文件**

- `src/plantpersulf/download/base.py`
- `src/plantpersulf/provenance/hashing.py`
- `tests/integration/test_real_file_download.py`
- `tests/unit/test_download_failure_is_fatal.py`

### RED

1. 下载一个官方小型 metadata/SDRF 文件；
2. 校验大小和 SHA256；
3. 模拟 HTTP 失败时必须抛出异常，且不得生成空“成功文件”。

### 验收

```bash
pytest tests/unit/test_download_failure_is_fatal.py -v
pytest -m network tests/integration/test_real_file_download.py -v
python scripts/fetch_registered_data.py \
  --accession PXD006140 \
  --file-class metadata,results
python -m plantpersulf.cli audit-files --accession PXD006140
```

提交：

```bash
git add src tests scripts data/registry
git commit -m "feat: add checksummed real-data download pipeline"
```

---

## Task 3：从真实数据构建测试微型夹具

**目标**：保证快速测试仍来自真实数据。

**创建文件**

- `scripts/build_real_fixtures.py`
- `tests/scientific/test_fixture_provenance.py`
- `tests/fixtures/real/...`

### 规则

- 脚本只能从已经校验的真实文件截取；
- 不能修改生物学值；
- 每个 fixture 记录 extraction command；
- fixture 仅用于测试，不进入训练和论文结果。

### RED

```python
def test_every_real_fixture_has_complete_manifest():
    manifests = list(Path("tests/fixtures/real").glob("**/source_manifest.json"))
    assert manifests
    for path in manifests:
        manifest = json.loads(path.read_text())
        assert manifest["accession"]
        assert manifest["source_sha256"]
        assert manifest["fixture_sha256"]
        assert manifest["biological_values_modified"] is False
```

### 验收

```bash
python scripts/build_real_fixtures.py --accessions PXD006140,PXD051570,GSE163745
pytest tests/scientific/test_fixture_provenance.py -v
```

提交：

```bash
git add scripts tests/fixtures tests/scientific
git commit -m "test: add provenance-locked real-data fixtures"
```

---

## Task 4：蛋白质组元数据和肽段解析

**目标**：统一真实搜索结果，保留研究、样本、肽段、修饰和蛋白映射。

**创建文件**

- `src/plantpersulf/proteomics/metadata.py`
- `src/plantpersulf/proteomics/peptide_parser.py`
- `src/plantpersulf/proteomics/site_normalizer.py`
- `tests/unit/test_peptide_parser_real_fixture.py`
- `tests/scientific/test_site_coordinates_match_sequence.py`

### 核心输出 schema

```text
study_accession
sample_id
source_file
spectrum_id
peptide_sequence
modified_sequence
protein_accession_raw
protein_accession_canonical
cys_position_in_peptide
cys_position_in_protein
modification_name_raw
evidence_level
quant_value
quant_unit
parser_version
source_sha256
```

### 必须处理

- 一条肽段映射多个蛋白；
- 同一蛋白多个 isoform；
- 位点位置冲突；
- 没有位点级信息的 protein-level 项目；
- 不同搜索引擎修饰命名；
- decoy 和 contaminant；
- 缺失样本条件。

### RED

使用 PXD006140 的真实 fixture 检查：

- 解析后肽段序列不变；
- Cys 位置能在真实 UniProt 序列中匹配；
- 不一致记录进入冲突表，不得自动纠正；
- protein-level-only 记录不得生成虚构位点。

### 验收

```bash
pytest tests/unit/test_peptide_parser_real_fixture.py \
       tests/scientific/test_site_coordinates_match_sequence.py -v
python -m plantpersulf.cli parse-proteomics --accession PXD006140
python -m plantpersulf.cli audit-sites --accession PXD006140
```

提交：

```bash
git add src tests
git commit -m "feat: parse traceable persulfidation peptide and site evidence"
```

---

## Task 5：正样本—未标记 benchmark v1

**目标**：建立不把“未检测到”误当阴性的真实 benchmark。

**创建文件**

- `src/plantpersulf/benchmark/labels.py`
- `src/plantpersulf/benchmark/unlabeled.py`
- `tests/scientific/test_no_unobserved_site_is_labeled_negative.py`
- `tests/scientific/test_positive_sites_have_experimental_evidence.py`

### RED

```python
def test_unobserved_cysteines_are_unlabeled_not_negative(real_benchmark):
    assert set(real_benchmark["label"].unique()) <= {"positive", "unlabeled"}
    assert "negative" not in set(real_benchmark["label"])
```

```python
def test_each_positive_has_traceable_experimental_evidence(real_benchmark):
    positives = real_benchmark.query("label == 'positive'")
    assert positives["study_accession"].notna().all()
    assert positives["source_sha256"].notna().all()
    assert positives["evidence_level"].isin(
        ["site_ms", "site_mutagenesis", "site_biochemical"]
    ).all()
```

### 输出

- `data/processed/benchmark_v1/sites.parquet`
- `data/processed/benchmark_v1/proteins.parquet`
- `data/processed/benchmark_v1/studies.parquet`
- `data/processed/benchmark_v1/conflicts.tsv`
- `data/processed/benchmark_v1/manifest.json`

### 验收

```bash
python scripts/build_benchmark.py --version v1
pytest tests/scientific/test_no_unobserved_site_is_labeled_negative.py \
       tests/scientific/test_positive_sites_have_experimental_evidence.py -v
python -m plantpersulf.cli audit-benchmark --version v1
```

提交：

```bash
git add src tests configs
git commit -m "feat: build positive-unlabeled persulfidation benchmark"
```

---

## Task 6：防泄漏数据切分

**目标**：防止同源蛋白、相同研究和已知机制泄漏。

**创建文件**

- `src/plantpersulf/benchmark/cluster_split.py`
- `src/plantpersulf/benchmark/study_split.py`
- `src/plantpersulf/benchmark/time_split.py`
- `src/plantpersulf/benchmark/leakage_audit.py`
- `tests/scientific/test_no_cluster_leakage.py`
- `tests/scientific/test_no_study_leakage.py`
- `tests/scientific/test_known_mechanisms_are_held_out.py`

### RED

```python
def test_no_protein_cluster_crosses_split(split_table):
    counts = split_table.groupby("cluster_id")["split"].nunique()
    assert counts.max() == 1
```

```python
def test_known_zhang_lab_sites_are_not_in_training(split_table):
    held_out = split_table.query("known_mechanism_holdout == True")
    assert not (held_out["split"] == "train").any()
```

### 验收

```bash
python -m plantpersulf.cli make-splits --benchmark v1 --split-version v1
python -m plantpersulf.cli audit-leakage --split-version v1
pytest tests/scientific/test_no_cluster_leakage.py \
       tests/scientific/test_no_study_leakage.py \
       tests/scientific/test_known_mechanisms_are_held_out.py -v
```

所有 split 文件在首次模型训练前 commit，后续不得覆盖，只能创建 v2。

提交：

```bash
git add src tests configs/splits
git commit -m "feat: freeze leakage-controlled benchmark splits"
```

---

## Task 7：真实序列和结构特征

**目标**：从版本化 UniProt 和 AlphaFold 数据提取特征。

**创建文件**

- `src/plantpersulf/download/uniprot.py`
- `src/plantpersulf/download/alphafold.py`
- `src/plantpersulf/features/sequence.py`
- `src/plantpersulf/features/esm2.py`
- `src/plantpersulf/features/structure.py`
- `tests/scientific/test_sequence_version_and_checksum.py`
- `tests/scientific/test_structure_residue_mapping.py`

### 特征

序列：

- Cys 邻域 ±5、±10、±20；
- 氨基酸组成；
- 电荷、疏水性和体积；
- 蛋白长度；
- Cys 密度；
- 冻结 ESM-2 residue embedding。

结构：

- pLDDT；
- relative solvent accessibility；
- 二级结构；
- 邻域原子/残基密度；
- 是否处于低置信区；
- 与已知 phosphosite 的序列/结构距离。

### 科学约束

- 结构缺失不能填充为“平均真实结构”；
- 缺失使用显式 mask；
- 低 pLDDT 单独分层报告；
- 坐标映射失败的记录不进入结构模型，但保留在序列模型；
- 模型结果必须分别报告有结构和无结构样本。

### 验收

```bash
python scripts/extract_features.py --benchmark v1 --feature-set sequence_v1
python scripts/extract_features.py --benchmark v1 --feature-set structure_v1
pytest tests/scientific/test_sequence_version_and_checksum.py \
       tests/scientific/test_structure_residue_mapping.py -v
```

提交：

```bash
git add src tests configs/features
git commit -m "feat: add versioned sequence and structure features"
```

---

## Task 8：无学习基线和传统 PU 基线

**目标**：先建立可信基线，禁止直接跳到复杂深度模型。

**创建文件**

- `src/plantpersulf/models/baselines.py`
- `src/plantpersulf/models/pu_risk.py`
- `src/plantpersulf/evaluation/metrics.py`
- `tests/unit/test_metrics.py`
- `tests/scientific/test_train_only_normalization.py`
- `tests/scientific/test_test_labels_never_used_for_model_selection.py`

### 模型

- motif/frequency baseline；
- accessibility baseline；
- Logistic/PU Logistic；
- Random Forest；
- XGBoost；
- ESM embedding + linear head。

### TDD 重点

- 标准化参数只能由训练折估计；
- 超参数只能由 validation 选择；
- test 只运行一次；
- 结果文件包含 seed、split、config hash；
- 每个模型使用相同 split。

### 验收

```bash
python scripts/run_experiment.py \
  --config configs/experiments/baseline_sequence_v1.yaml
pytest tests/unit/test_metrics.py \
       tests/scientific/test_train_only_normalization.py \
       tests/scientific/test_test_labels_never_used_for_model_selection.py -v
```

提交：

```bash
git add src tests configs/experiments
git commit -m "feat: add leakage-safe PU baselines"
```

---

## Task 9：Structure-aware PU ranker

**启动门槛**

只有 Task 8 满足以下条件才启动：

- 所有防泄漏测试通过；
- 至少两个 leave-study-out fold 可执行；
- ESM + linear baseline 能输出稳定结果；
- 结果可由命令重现；
- 没有未登记输入。

**创建文件**

- `src/plantpersulf/models/structure_ranker.py`
- `src/plantpersulf/models/calibration.py`
- `tests/unit/test_structure_ranker_shapes.py`
- `tests/scientific/test_missing_structure_mask_is_respected.py`
- `tests/scientific/test_model_release_is_reproducible.py`

### 架构

最小实现：

```text
sequence branch
+ frozen residue embedding branch
+ structure feature branch
+ missingness masks
→ gated fusion
→ PU ranking head
→ uncertainty/calibration head
```

禁止在 MVP 添加：

- 大型异质图；
- 多任务生成模型；
- 复杂自监督预训练；
- 端到端微调 650M 模型；
- 无法做消融的模块。

### 必须消融

- sequence only；
- sequence + ESM；
- sequence + structure；
- full；
- no pLDDT；
- no accessibility；
- no study context。

### 验收

```bash
python scripts/run_experiment.py \
  --config configs/experiments/pu_ranker_v1.yaml
pytest tests/unit/test_structure_ranker_shapes.py \
       tests/scientific/test_missing_structure_mask_is_respected.py \
       tests/scientific/test_model_release_is_reproducible.py -v
```

提交：

```bash
git add src tests configs
git commit -m "feat: add structure-aware persulfidation PU ranker"
```

---

## Task 10：严格外部验证与已知机制回顾性恢复

**目标**：验证跨研究能力，不能只做随机 CV。

**创建文件**

- `src/plantpersulf/evaluation/bootstrap.py`
- `src/plantpersulf/evaluation/permutation.py`
- `src/plantpersulf/evaluation/external_validation.py`
- `tests/scientific/test_external_study_not_seen_during_training.py`
- `tests/release/test_claims_have_supporting_tables.py`

### 必须输出

- 每个 PXD 的 leave-study-out 结果；
- protein cluster bootstrap；
- 时间切分结果；
- 已知张华团队机制位点的排名；
- 失败位点和不可映射位点；
- 模型适用域；
- 与所有 baseline 的效应量和 CI；
- 所有 seed 的结果，不只最佳 seed；
- 每个控制的 percentile rank、applicability-domain 状态和不确定性；
- 每个控制与无学习、传统模型和 ESM 基线的并列结果；
- `mechanism_lineage_id` 和独立验证单元状态；
- 所有未恢复、不可映射和证据不足的控制，不得只展示成功案例。

### 结论闸门

只有满足下列条件才允许在合作摘要中写“具有预测价值”：

1. 完整模型在至少两个独立 held-out study 上优于无学习基线；
2. 95% CI 不完全跨越零效应；
3. 已知机制恢复不是训练泄漏；
4. 结构增强至少在有结构、足够 pLDDT 的子集上产生可解释增益；
5. 结果不是由一个高同源 protein cluster 驱动。

否则必须写：

> 当前公开数据不足以证明跨研究预测能力，模型仅用于候选组织与假设生成。

### 未来 RED 规格

1. `test_known_control_cannot_be_training_and_recovery`：控制同时进入训练并被报告为独立恢复时失败；
2. `test_duplicate_mechanism_lineage_is_not_independent`：同一机制谱系被计为多个独立验证单元时失败；
3. `test_control_rank_is_not_used_for_model_selection`：控制排名参与模型或超参数选择时失败；
4. `test_predictive_claim_requires_gate2_go`：Gate 2 为 `STOP` 时声称通用预测价值时失败；
5. `test_failed_and_unmappable_controls_are_reported`：失败或不可映射控制被遗漏时失败。

测试只能使用已登记真实控制记录或纯软件策略文本，不得创建虚构位点。

### 验收

```bash
python scripts/validate_external.py --model-release pu_ranker_v1
pytest tests/scientific/test_external_study_not_seen_during_training.py \
       tests/release/test_claims_have_supporting_tables.py -v
```

提交：

```bash
git add src tests manuscripts
git commit -m "feat: add leave-study-out and mechanism holdout validation"
```

---

## Task 11：番茄成熟上下文整合

**目标**：使用真实番茄成熟数据给候选添加上下文，不能把上下文关联解释为硫巯基化阳性。

**创建文件**

- `src/plantpersulf/transcriptomics/expression.py`
- `src/plantpersulf/transcriptomics/gene_id_mapping.py`
- `src/plantpersulf/tomato/ripening_context.py`
- `src/plantpersulf/tomato/phospho_context.py`
- `tests/scientific/test_tomato_samples_match_official_metadata.py`
- `tests/scientific/test_gene_id_mapping_conflicts_are_not_silenced.py`

### 处理

PXD051570：

- 五成熟阶段蛋白定量；
- 磷酸化位点定量；
- stage-specific dynamics；
- 蛋白 abundance 与 phosphosite 解耦；
- phosphosite 与候选 Cys 的序列/结构距离。

GSE163745 等：

- 样本条件来自 GEO metadata；
- 使用训练集/对比定义生成差异表达；
- 不跨版本静默映射；
- 每个数据集单独分析后再做方向一致性整合。

### 输出

- `tomato_protein_stage_context.parquet`
- `tomato_phosphosite_stage_context.parquet`
- `tomato_expression_context.parquet`
- `tomato_id_mapping_conflicts.tsv`
- `tomato_context_manifest.json`

### 验收

```bash
python -m plantpersulf.cli build-tomato-context \
  --proteomics PXD051570 \
  --transcriptomics GSE163745,GSE142713,GSE267238
pytest tests/scientific/test_tomato_samples_match_official_metadata.py \
       tests/scientific/test_gene_id_mapping_conflicts_are_not_silenced.py -v
```

提交：

```bash
git add src tests
git commit -m "feat: integrate real tomato ripening and phosphoproteomic context"
```

---

## Task 12：机制卡与去循环候选评分

**目标**：将已知 H₂S/成熟机制转为可审计证据，不让已知基因标签直接决定候选排名。

**创建文件**

- `configs/mechanism_cards/h2s_wrky6_phosphorylation.yaml`
- `configs/mechanism_cards/h2s_brg3_ubiquitination.yaml`
- `configs/mechanism_cards/h2s_erfd2_ethylene.yaml`
- `src/plantpersulf/tomato/mechanism_cards.py`
- `src/plantpersulf/tomato/candidate_score.py`
- `tests/scientific/test_candidate_score_is_decircularized.py`

### 评分组成

候选主评分：

- 模型外部验证得分；
- 模型不确定性；
- 番茄成熟蛋白动态；
- 转录动态一致性；
- phosphosite 串扰支持；
- 结构可操作性；
- 新颖性；
- 实验可行性。

已知 H₂S/成熟基因集只能：

- 作为解释证据；
- 做富集；
- 做阳性控制；
- 生成机制卡邻域。

不得直接作为训练标签或主评分中的硬编码加分项。

### RED

将一个仅因出现在已知机制卡而得高分的候选输入评分器，测试必须保证它不会在没有独立数据支持时自动进入 Top 候选。

### 验收

```bash
python scripts/rank_tomato_candidates.py \
  --model-release pu_ranker_v1 \
  --context-version tomato_context_v1 \
  --score-version candidate_score_v1
pytest tests/scientific/test_candidate_score_is_decircularized.py -v
```

提交：

```bash
git add src tests configs/mechanism_cards
git commit -m "feat: add decircularized tomato candidate prioritization"
```

---

## Task 13：证据卡和冻结候选发布包

**目标**：形成可以和张华老师讨论、且不会随意改变的真实候选包。

**创建文件**

- `src/plantpersulf/tomato/evidence_cards.py`
- `src/plantpersulf/reporting/model_card.py`
- `src/plantpersulf/reporting/candidate_release.py`
- `docs/zhang_collaboration_brief.md`
- `docs/wetlab_validation_matrix.md`
- `tests/release/test_candidate_release_traceability.py`
- `tests/release/test_release_rebuild_is_identical.py`

## 13.1 每个候选证据卡必须包含

- tomato protein/gene ID；
- canonical protein sequence accession/version；
- Cys 位点；
- 序列上下文；
- AlphaFold model/version 和 pLDDT；
- solvent accessibility；
- 模型分数和置信区间；
- applicability-domain 状态；
- 外部研究恢复证据；
- PXD051570 蛋白动态；
- PXD051570 phosphosite 邻近/结构邻近证据；
- GEO 表达动态；
- 机制卡邻域；
- 是否为已知机制；
- 文献新颖性审计；
- 实验建议；
- 所有来源 accession；
- 所有关键文件 SHA256；
- 生成代码 commit。

## 13.2 发布等级

- Tier A：多源真实证据、模型适用域内、低不确定性、可实验；
- Tier B：模型和部分上下文支持，但存在缺失模态；
- Tier C：探索性候选，不进入首轮湿实验。

首轮建议：

- Tier A：5–10 个；
- Tier B：10–20 个；
- 低排名 matched controls：5–10 个；
- 已知阳性 controls：3–5 个。

具体数量依据真实结果决定，禁止为了凑数量降低门槛。

## 13.3 湿实验盲法建议包

`docs/wetlab_validation_matrix.md` 应提出：

1. 盲法编号，不向实验人员显示模型排名；
2. Top 候选、低排名 matched controls、已知阳性 controls；
3. 统一 tag-switch/biotin-switch 或团队成熟方法；
4. 靶向 LC-MS/MS；
5. Cys→Ser/Ala 突变；
6. H₂S donor、LCD1 相关材料或团队现有体系；
7. 不成功结果完整返回；
8. 冻结候选后不得因实验结果重新排序并冒充前瞻预测；
9. 记录批次、操作者、抗体、仪器、原始文件和排除标准；
10. 湿实验前冻结 `frozen_analysis_plan`，声明终点、命中定义、排除标准、匹配变量、统计比较和失败处理；
11. `blind_id` 是实验人员可见的唯一候选标识，排名和 Tier 保持隐藏；
12. 高排名候选、matched controls 和已知阳性使用一致检测与排除规则；
13. 系统主张依赖多个独立新位点或蛋白的整体富集，不得只挑一个成功案例；
14. 深入机制对象可从真实命中中选择，但不得回写或重排 `candidate_release_v1`；
15. 数据返回后新增 `prospective_validation_v1`，不覆盖原发布包。

### 未来 RED 规格

1. `test_nature_readiness_requires_prospective_blind_validation`：无冻结候选、匹配对照和盲法结果时，不得标记 Nature-family ready；
2. `test_candidate_release_cannot_be_reranked_after_wetlab`：湿实验后修改顺序或覆盖原发布包时失败；
3. `test_unsuccessful_candidates_are_retained`：失败候选被删除时失败；
4. `test_prediction_is_not_causal_mechanism`：模型排名、位点修饰和因果表型混为同一证据层级时失败。

### 验收

```bash
python scripts/build_candidate_release.py \
  --candidate-score candidate_score_v1 \
  --output releases/candidate_release_v1
pytest tests/release/test_candidate_release_traceability.py \
       tests/release/test_release_rebuild_is_identical.py -v
sha256sum -c releases/candidate_release_v1/SHA256SUMS
```

提交和标签：

```bash
git add releases docs src tests
git commit -m "release: freeze tomato persulfidation candidate package v1"
git tag -a candidate-release-v1 -m "Frozen pre-wetlab candidate release v1"
```

---

# 8. 全量测试命令

## 快速测试

```bash
pytest tests/unit -q
```

## 真实微型夹具科学测试

```bash
pytest tests/scientific -q
```

## 官方网络集成测试

```bash
pytest -m network tests/integration -q
```

## release gate

```bash
pytest tests/release -q
python -m plantpersulf.cli audit-registry
python -m plantpersulf.cli audit-files
python -m plantpersulf.cli audit-benchmark --version v1
python -m plantpersulf.cli audit-leakage --split-version v1
ruff check .
mypy src/plantpersulf
```

## 一键复现

```bash
make metadata
make download-registered
make benchmark-v1
make features-v1
make baselines-v1
make model-v1
make external-validation-v1
make tomato-context-v1
make candidate-release-v1
make verify-release-v1
```

一键复现不得：

- 调用合成模式；
- 自动跳过下载失败；
- 用旧结果冒充新运行；
- 在缺失输入时继续生成空图；
- 修改冻结 split；
- 自动覆盖发布包。

---

# 9. Stop/Go 科学决策门槛

## Gate 1：数据可用性

**GO**

- 至少两个硫巯基化研究可解析到可靠蛋白/位点证据；
- 来源表和样本表完整；
- ID 映射冲突率可量化；
- 真实 benchmark 可建立。

**STOP/调整**

- 只有蛋白级列表且无法获得位点；
- 关键文件损坏或无法合法获取；
- 研究间修饰定义不可比；
- 阳性位点数过少，不足以训练。

若 STOP，项目转为：

> 植物硫巯基化数据资源与系统性证据审计

不得强行训练深度模型。

## Gate 2：跨研究预测

**GO**

- 至少两个 held-out study 显示稳定富集；
- 结果不由单一同源 cluster 驱动；
- 置信区间和统计检验支持；
- 模型校准可接受。

**STOP/调整**

- 只在随机 split 有效；
- leave-study-out 失效；
- ESM/结构模型不优于简单基线；
- 指标高度依赖研究批次。

若 STOP，候选排名必须降级为：

> evidence integration，不声称通用 prediction。

## Gate 3：番茄候选价值

**GO**

- 至少一批候选位于模型适用域；
- 有真实成熟蛋白/转录/PTM 上下文；
- 不是全部已知经典成熟基因；
- 能提出可证伪实验。

**STOP/调整**

- 番茄跨物种迁移不稳定；
- 候选主要由缺失值或 ID 映射驱动；
- 候选全部是已知机制复述。

## Gate 4A：非正式联系

即使 Task 5 为 `STOP`，仍可沟通研究方向、科学问题、已审计数据缺口、能力互补、数据请求和盲法实验设想。不得展示不存在的模型分数或候选排名，不得把计划写成结果，不得声称预测价值、新机制或确定的期刊结果。

## Gate 4B：带结果正式汇报

带模型结果和候选证据卡正式汇报至少要求：

1. registry 和可训练 benchmark 通过审计；
2. 防同源、研究、物种、时间和已知机制泄漏的 split 冻结；
3. 完成无学习、传统模型和 ESM 基线；
4. Gate 2 状态和允许主张明确；
5. 所有控制均有训练外恢复结果，包括失败和不可映射项；
6. 有 5–20 个初步番茄候选证据卡，不得为凑数量降低门槛；
7. 报告适用域、不确定性、基线和限制；
8. 准备好数据、实验能力和盲法验证请求清单。

若 Gate 2 为 `STOP`，只能展示证据审计和 `evidence_integration`，不得称为通用 predictor。

## Gate 5：Nature-family 故事就绪度

本 Gate 只评价联合研究的故事证据结构，不预测或保证期刊结果。

只有同时满足以下条件，才允许内部标记 `nature_family_story_ready`：

1. 主张超越另一个单蛋白—单个位点机制；
2. Gate 2 已 `GO`；
3. 候选、matched controls、终点和统计方案在湿实验前冻结；
4. 完成前瞻性盲法验证并保留全部结果；
5. 高排名组显示多个独立新位点或蛋白的整体富集；
6. 至少一至两个中心命中完成深入生化或遗传验证；
7. 提出可检验的 PTM 串扰或其他一般原则；
8. 与主张相称地覆盖条件、材料、品种或物种；
9. 连接生长、产量、成熟、采后品质或更广泛植物问题；
10. 数据、模型、协议、来源链和失败结果在协议允许范围内可复现。

缺失核心条件时，降级为证据审计、benchmark/候选组织、回顾性预测、无一般机制主张的前瞻富集，或单候选机制研究。模型、随机交叉验证、单个新位点或单次成功实验均不足以通过 Gate 5。

---

# 10. 与张华老师沟通时的数据请求

优先询问是否存在以下真实、未发表或可合作产生的数据：

1. SlLCD1 overexpression、sllcd1、WT 的原始 persulfidome；
2. 不同成熟阶段的 H₂S 处理与对照 persulfidome；
3. 对应 total proteome；
4. 对应 phosphoproteome 或 ubiquitinome；
5. 原始质谱文件和搜索参数；
6. sample sheet、实验批次、重复、处理浓度和时间；
7. RNA-seq 原始数据；
8. 已有但尚未深入验证的候选列表；
9. 可用于盲法批量验证的样品和实验方法；
10. 可做遗传救援的材料储备。

合作数据进入项目必须新增：

```text
data/registry/collaboration_datasets.tsv
data/registry/data_use_confirmation/
configs/data_sources_private.yaml
```

私有数据不得上传公开仓库。结果发布遵循双方书面约定。

---

# 11. 首次汇报材料结构

## 两页概念书

### 第 1 页

- 科学问题；
- 公开数据来源；
- 防泄漏设计；
- 当前模型结果；
- 已知机制回顾性恢复；
- 主要限制；
- 当前最高通过 Gate；
- 当前类型：`evidence_audit`、`evidence_integration`、`retrospective_prediction` 或 `prospective_validation`；
- 排除训练的已知控制；
- 缺失关键证据；
- `STOP` Gate 对外主张的精确降级文字。

### 第 2 页

- 番茄候选证据卡示例；
- 候选/阴性/阳性盲法实验设计；
- 需要张华老师团队提供的数据；
- 双方分工；
- 可能的文章层级和降级路线。

## 十页 PPT

1. H₂S—硫巯基化—番茄成熟研究背景；
2. 张华团队已完成的机制与尚未解决的系统问题；
3. 真实数据地图；
4. 数据完整性与防泄漏协议；
5. PU benchmark；
6. 模型及外部验证；
7. 已知机制恢复；
8. 番茄候选证据卡；
9. 前瞻性盲法验证设计；
10. 合作分工与数据需求。

---

# 12. 预期论文路线

## 计算预研不足以支持高水平机制结论时

可形成：

- 数据资源/benchmark；
- 植物 PTM 位点预测方法；
- 番茄候选假设生成。

## 获得张华团队系统数据和批量验证后

文章主张可以逐步升级：

## Level 0：数据与证据审计
> 建立可追踪数据资源，报告证据、方法差异、冲突和缺口。

## Level 1：benchmark 与候选组织
> 建立真实 PU benchmark、基线和候选组织框架，不声称通用预测。

## Level 2：跨研究预测与回顾性恢复
> 在防泄漏 held-out studies 上证明跨研究价值，并恢复排除训练的机制。

## Level 3：前瞻性盲法富集
> 冻结番茄候选在盲法实验中相对 matched controls 富集真实位点。

## Level 4：系统规律、深入机制与植物意义
> 多个独立新命中支持可推广原则，一至两个中心命中完成机制，并连接重要植物或农业性状。

只有 Level 4 且 Gate 5 为 `READY`，才允许内部按 Nature-family 故事组织材料；不保证期刊结果。

---

# 13. Codex 每个任务的汇报模板

```text
Task:
Commit:

1. RED
- 测试命令：
- 失败输出摘要：
- 确认失败是因为目标功能缺失，而不是测试错误：是/否

2. GREEN
- 最小实现：
- 通过命令：
- 测试输出摘要：

3. REAL DATA
- 使用 accession：
- 使用文件：
- SHA256：
- 是否修改生物学值：否
- 是否存在未登记输入：否

4. SCIENTIFIC INTEGRITY
- 是否使用合成/模拟生物学数据：否
- 是否把未检测位点标成阴性：否
- 是否读取测试标签调参：否
- 是否覆盖冻结 split：否

5. FILES
- 新增：
- 修改：
- 删除：

6. LIMITATIONS
- 当前无法完成的内容：
- 原因：
- 对后续科学结论的影响：

7. NEXT
- 下一任务：
- 开始条件是否满足：

8. CLAIM SAFETY
- 当前最高通过 Gate：
- 当前允许 claim level：
- 已知控制是否全部排除训练：是/否/不适用
- 是否存在机制谱系重复计数：否/是（若是，任务失败）
- 是否产生前瞻性盲法结果：否/是
```

---

# 14. 最终不可违反的红线

1. 没有真实数据，不生成科研结果。
2. 下载失败，不生成替代结果。
3. 解析失败，不猜测字段。
4. 未检测，不等于阴性。
5. 同源泄漏，不得作为泛化。
6. 已知机制进入训练，不得再作为独立验证。
7. 随机划分结果，不得代替跨研究结果。
8. 单次最优 seed，不得作为主结果。
9. 模型注意力，不等于因果机制。
10. AlphaFold 结构，不等于实验结构。
11. 公开数据关联，不等于湿实验验证。
12. 失败候选必须保留，不能只展示命中者。
13. 冻结候选后不得根据实验结果回填排名。
14. 所有结论必须有表、图、配置、来源和 commit 支撑。
15. 达不到门槛时主动降级结论，不得包装升级。
16. 同一蛋白、位点或机制谱系不得重复计算为多个独立验证。
17. Gate 2 为 `STOP` 时不得声称通用预测价值。
18. 没有冻结候选和前瞻性盲法验证，不得标记 Nature-family ready。
19. 单个成功候选不得替代多个独立命中的整体富集证据。
20. 预测、位点修饰、生化功能和因果表型必须分层表述。
21. Nature-family story readiness 是内部闸门，不是期刊结果保证。

---

## 当前项目状态

Phase A 证据补齐已完成：从两个已登记、SHA256 审计的来源解析出约 393 个坐标验证过的位点级 persulfidation 证据——PXD006140 Dataset S3（作者自定义 Sulfide/CN-Biotin-Sulfide 修饰，320 个 site_ms）与 PXD024061 MaxQuant Sulfide(C)/CianoBiotin(C) 位点表（73 个 class-I）。fail-closed 的 readiness v2 gate 据此判定 **Gate 1（数据可用性）= GO**，并已构建 PU benchmark v1、防泄漏 split、序列/ESM-2 特征与 motif 基线。

关键限制：两个研究同属 Seville（Romero/Gotor）实验室、同一 tag-switch 化学、同一物种（拟南芥）。因此当前 GO **仅为数据可用性层面**，**不等于** Gate 2 的跨研究预测价值，也**不满足** Gate 5。在 Gate 2 以防泄漏 held-out study 证明预测价值之前，对外只能作 `evidence_integration`，不得声称通用 predictor，不得标记 Nature-family ready。该状态已超过 Gate 4A（非正式联系），但**尚未**满足 Gate 4B（需 Gate 2 状态明确且所有控制有训练外恢复结果）。真正独立化学/实验室/物种的位点级 persulfidome 仍是对张华团队的优先数据请求。

---

# 15. 立即执行顺序

现在只做以下顺序：

1. Task 0：科研诚信闸门；
2. Task 1：官方元数据注册表；
3. Task 2：真实文件下载和 SHA256；
4. Task 3：真实测试微型夹具；
5. Task 4：PXD006140 解析；
6. Task 5：PU benchmark；
7. Task 6：防泄漏 split；
8. Task 7–10：模型和外部验证；
9. Task 11–13：番茄上下文和候选发布。

在 Task 6 之前，不实现复杂模型。  
在 Task 10 之前，不生成对外候选清单。  
在 Task 13 之前，不使用“前瞻性候选发布”表述。  
在湿实验真实结果返回前，不声称发现新机制。
