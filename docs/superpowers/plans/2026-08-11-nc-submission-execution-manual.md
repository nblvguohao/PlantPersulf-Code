# NC 投稿操作执行手册 — PlantPersulf 拟南芥→番茄验证路线

> **用途**：本文是从当前状态到 Nature Communications 投稿的**操作执行手册**（living document）。
> 它不替代 `docs/PlantPersulf_Code_TDD_Codex.md`（权威规格）、`docs/phase_f_gate2_decision.md`
> （Gate 2 决策记录）、`docs/superpowers/specs/2026-07-21-phase-f-pu-ranker-design.md`
> （two-track 政策）或 `docs/zhang_collaboration_workpackages_2026-08-10.md`（合作工作包）——
> 本文提供它们都没有的：**顺序执行骨架、每阶段验收标准、时间线与验证协议**。
> 每完成一个阶段，在 §0 状态表追加带日期条目。

## 0. 状态总览（living anchor）

| 阶段 | 状态 | 负责人 | 关键产物 |
|---|---|---|---|
| 0 版本清理与可复现 | 未开始 | 我们 | 干净 clone 全绿 |
| 1 展示轨数字 + 开发稳定性轨 + 决策记录 | **进行中** | 我们 | 蛋白切分 summary、cluster CV 轨、决策附录 |
| 2 实验设计冻结（含 lockbox 协议） | 未开始 | 我们起草 / 张华确认 | 冻结协议包 + candidate_release_v1 |
| 3 番茄发现队列 | 未开始 | 张华提供 / 我们接收 | tomato_benchmark_v1 |
| 4 lockbox 外部验证 + 公平比较 | 未开始 | 双方 | lockbox 解盲分析（**NC 主结论**） |
| 5 前瞻盲法候选实验 | 未开始 | 张华执行 | 盲法富集结果 |
| 6 机制闭环 | 未开始 | 双方 | 1-2 命中证据链 |
| 7 投稿包 | 未开始 | 我们 / 张华协议 | NC 投稿材料 |

## 1. 为什么存在这份手册

### 1.1 要讲的故事
植物 persulfidation 可泛化规律 + 番茄成熟验证：从"可追踪位点级证据整合与防泄漏建模"到
"冻结模型与候选 → 盲法富集验证 → 1-2 个因果机制 → 独立复现与完整开放发布"（Level 4，
codex §12）。

### 1.2 诚实的当前位置
- Gate 2 = STOP 3/5；条件 1 是数据来源锁（2 研究同实验室，`configs/gate2_v1.yaml` 冻结
  `studies_are_independent: false`）。
- 按 2026-08-11 评估架构决策：**NC 主结论由未来番茄 lockbox 队列承载**；当前 LSO 结果
  保留为参考轨；Level 2 在跨研究证据出现前不可声称。
- 所有对外文本必须通过 `verify_predictive_claims`（`conclusion_gate.py:154`）→ `[]`。

### 1.3 来源索引（引用不复制）
- Gate 1-5 定义：codex §9；Level 0-4：codex §12
- Gate 2 决策：`docs/phase_f_gate2_decision.md` + `results/external_validation/pu_ranker_v1/gate2_decision.json`
- 评估政策：`docs/superpowers/specs/2026-07-21-phase-f-pu-ranker-design.md`（two-track + 2026-08-10 报告顺序决策）
- 合作工作包：`docs/zhang_collaboration_workpackages_2026-08-10.md`
- 竞争工具审计：`docs/competitor_data_audit_2026-08-10.md`
- 数据审计：`docs/phase_z_evidence_audit.md`（§5.2-5.4 跨物种措辞锁定）

## 2. 评估架构（2026-08-11 决策——验证主线）

| 场景 | 评估方式 | 角色 |
|---|---|---|
| 番茄本地模型 | **5 折同源簇分组**（MMseqs2 cluster CV） | 主评估（番茄侧） |
| 拟南芥（现有基准） | **补充重复 5/10 折同源簇分组**（`pu_ranker_cluster_cv_v1`） | 开发稳定性（模型选择/消融） |
| **NC 主结论** | **未来新番茄队列 = lockbox**（隐藏验证） | **主证据——唯一锚点** |
| leave-one-study-out | **有第二个独立研究后**再做 | 补充/未来 |
| 跨物种（PXD063170/PXD072089） | 额外挑战集 | robustness 展示，不替代外部验证 |

硬约束：NC 稿件"外部验证/泛化"主张只能引用 lockbox 结果；LSO、跨物种、蛋白切分数字各归
其位；lockbox 协议进入阶段 2 冻结文档；现有 LSO 冻结结果保留为参考轨（不删除、不重判）。

## 3. 阶段

### 阶段 0 版本清理与可复现（我们；张华确认可公开范围）
- [ ] 修复 Ruff/mypy/pytest；环境锁定（依赖 pin）
- [ ] 清理在审稿件/版权 PDF（`docs/compete/btaf078.pdf` 等）公开展示风险；`git ls-files` 审计
- [ ] 修正番茄 cluster 数表述不一致；刷新 `docs/data_reproduction.md`
**PASS**：干净 clone 全绿；`git ls-files` 无版权材料；复现文档 1:1。**FAIL**：任何红 → 阶段 7 DOI/代码发布被阻塞。

### 阶段 1 展示轨数字 + 开发稳定性轨 + 决策记录（我们；进行中）
1. [ ] 下载 A100 蛋白切分结果（metrics.tsv + manifest.json）→ 校验 `experiment.name` 身份 + SHA256
2. [ ] `scripts/analyze_protein_split.py`（已落位 + 测试通过）：10 次重复均值±SD、结构增益
       `paired_delta_ci`、效应 delta（基线缺口 Route A/B）
3. [ ] Route A 基线（`pu_ranker_protein_split_v1_baseline.yaml`，已落位，本机已跑完 10 seeds）
4. [ ] cluster CV 开发稳定性轨（`pu_ranker_cluster_cv_v1.yaml`，runner + 测试已落位，运行中）
5. [ ] 决策记录附录（`docs/phase_f_gate2_decision.md`）：评估架构声明 + 展示数字 + 三轨对照表 + Route 记录
6. [ ] 回归：`validate_external.py` → `gate2_decision.json` 字节不变；tests/release + tests/scientific 全绿
7. [ ] 附录 A 措辞模板逐条过 `verify_predictive_claims` → `[]`
**PASS**：见上全部；**FAIL**：不得对外展示任何数字，进入 DP1。
**里程碑 M1**。

### 阶段 2 实验设计冻结（我们起草；张华确认）——含 lockbox 协议
- [ ] `results/candidates/candidate_release_v1/`（blind_id、证据卡、matched controls、已知阳性、哈希 manifest；分数对实验者不可见）
- [ ] 6 份协议：data_dictionary / sample_numbering / split_boundary / **lockbox 协议**（隐藏标签机制、批次、样本量功效分析、解盲条件、标签在冻结前不开放）/ statistical_analysis_plan（命中定义、终点、预注册富集检验、事后不加样本）/ 盲法协议
- [ ] 番茄本地模型 5 折同源簇分组配置冻结；全部 SHA256 + git tag
**PASS**：新验证数据到达前全部哈希+tag；冻结后零重排；张华确认实验条件与 lockbox 机制。
**里程碑 M2**。Gate 3 证据 + Gate 5 条件 3。

### 阶段 3 番茄发现队列（张华提供；我们接收）
- [ ] kiae271 可复现包接收（原始 MS、搜库参数、修饰定义、FDR/定位概率、sample sheet、ERF.D3 CDS）
- [ ] 校验→解析→坐标核验（复用 PXD006140 parser 模式）→ 番茄 PU benchmark + 5 折同源簇分组 + 泄漏审计 → 冻结 `tomato_benchmark_v1`
- [ ] 数据到达前先建 dry-run 接收管线（与阶段 0 并行）
**PASS**：完整 provenance 链；泄漏审计干净；训练前 split 冻结。**DP2**：数据不可得 → 降级。

### 阶段 4 lockbox 外部验证 + 公平比较（双方）——NC 主结论锚点
- [ ] lockbox 解盲（按冻结协议）→ 冻结模型/特征/阈值哈希核验 → 预注册分析 → **NC 主结论**
- [ ] 与 pCysMod/Sul-BertGRU 同队列公平比较（按 `docs/competitor_data_audit_2026-08-10.md` 前置：Sul-BertGRU 复现 2,705 集/中心 C 偏差/物种分层/环境隔离；pCysMod 仅 web server 定性）
- [ ] 跨物种挑战集展示（PXD063170/PXD072089，措辞锁定）；LSO 重做（第二独立研究后）
**PASS**：lockbox 按冻结协议完成；哈希一致；零事后调参；各证据各归其位。**里程碑 M4**。

### 阶段 5 前瞻盲法候选实验（张华执行）
- [ ] 冻结候选包 → 实验者只见 blind_id → 冻结 SAP 预注册分析 → 全结果（含阴性）返回后解盲
**PASS**：按冻结 SAP 显著富集；100% 返回。**DP3**：不显著 → 回顾性论文。

### 阶段 6 机制闭环（双方）
- [ ] 1-2 中心命中完整链：位点确认 → Cys 突变 → 功能 → 表型 → 救援 →（若主张串扰）直接测量 persulfidation×phosphorylation 互赖（对齐 `docs/superpowers/plans/2026-07-28-ptm-crosstalk-grammar.md`）→ 连接成熟/胁迫表型
**PASS**：≥1 命中完整链；串扰必须有直接测量。**里程碑 M6**。Level 4；Gate 5 条件 6/7/9。

### 阶段 7 投稿包（我们；张华确认作者/数据协议）
- [ ] 图表 + Source Data + Nature ML checklist + Reporting Summary + 代码/数据 DOI + Data Availability + 作者治理 + 模拟审稿
- [ ] 合规扫描：每个展示数字附 limitation verbatim；`verify_predictive_claims` 全扫；`find_unsupported_claims` → `[]`
**PASS**：清单完整；主张↔表映射完整；DOI 已铸；作者协议已签。**里程碑 M7** = Gate 5 READY → 投稿；否则 DP1 降级。

## 4. 决策点与降级路线

| 决策点 | 触发 | 走向 |
|---|---|---|
| DP1 | 阶段 1 后：评估架构落地；或 lockbox 路径不可行 | 论文上限 Level 0/1（证据审计 / benchmark 方法学）；NC 主结论在 lockbox 完成前不存在 |
| DP2 | 阶段 3/4：张华数据不可得 | 跨物种迁移 + 恢复叙事（措辞锁定 phase_z §5.2-5.4） |
| DP3 | 阶段 5：富集不显著 | 仅回顾性论文（Level 2/1 主张，无前瞻声明） |
| DP4 | 持续 | Gate 5 十条件追踪（下表） |

**Gate 5 十条件追踪**：1 超单蛋白机制→阶段 5/6；2 Gate 2 GO→**阶段 4 lockbox**（LSO 待第二独立研究）；3 冻结前置→阶段 2；4-5 盲法完整/独立命中富集→阶段 5；6-7 深验证/可检验原则→阶段 6；8 覆盖度→阶段 3-4；9 植物意义→阶段 6；10 可复现→阶段 0/7。

## 5. 时间线与依赖

- **Track A（内部，立即并行）**：阶段 0 ∥ 阶段 1 ∥ 阶段 2 起草 ∥ 张华外联（草稿已存在）∥ 阶段 3 接收管线 dry-run
- **关键路径（NC 主证据）**：2 → 3 → 4（lockbox）→ 5 → 6 → 7；里程碑 M1-M7；张华门槛弧：2 确认、3、4、5、6、7 协议
- **支线**：LSO 重做（第二独立研究后）、跨物种挑战集
- **关键路径上无模型开发**（"不要继续堆叠更大模型"；模型开发只在阶段 3 番茄 5 折框架内）

## 6. 验证与手册更新协议

- 机械化检查：pytest 三套件（锁定测试：`test_gate2_ignores_within_dataset_split_metrics.py`、`test_protein_split_regime.py`、`test_protein_split_summary.py`、`test_cluster_cv_regime.py`）；`validate_external.py` → `gate2_decision.json` 字节比对；SHA256 manifest 于每次接收；措辞扫描；ruff/mypy
- 手册更新：§0 状态表为 living anchor；每阶段追加带日期条目（日期、负责人、产物、PASS 证据、偏差）；冻结产物修订一律带日期+哈希，绝不静默编辑

## 附录 A：措辞锁定

（阶段 1 完成后填充，逐条过 `verify_predictive_claims`）
- 降级声明（Gate 2 STOP 期间）：现有冻结文本 verbatim
- 展示轨 limitation（蛋白切分 / cluster CV）：对应 config limitation verbatim
- 跨物种转移措辞：per phase_z §5.2-5.4

## 附录 B：产物清单

（随阶段推进更新：路径、哈希来源、冻结状态）
- `results/experiments/pu_ranker_protein_split_v1/{metrics.tsv,manifest.json,summary.json,summary.md}`
- `results/experiments/pu_ranker_protein_split_v1_baseline/`
- `results/experiments/pu_ranker_cluster_cv_v1/`
- `docs/phase_f_gate2_decision.md`（2026-08-11 附录）
