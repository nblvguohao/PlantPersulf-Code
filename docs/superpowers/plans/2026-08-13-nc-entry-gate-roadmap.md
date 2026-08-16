# Nature Communications 入场门路线图 — PlantPersulf 多物种 v2 后续计划

> **目标**：在不使用冻结测试或盲法结果反向开发模型的前提下，构建达到 *Nature Communications*（NC）入场门的完整证据链；若预注册门槛未满足，则直接按对应证据等级正式投稿 *Plant Physiology*。
>
> **与现有文档的关系**：本计划是对 `docs/superpowers/plans/2026-08-11-nc-submission-execution-manual.md` 的**策略更新与门控细化**，并承接 `docs/superpowers/plans/2026-08-13-multispecies-v2-task9-6-claim-gate.md` 已完成的统计闸门基础设施。

## 0. 当前状态

- `configs/experiments/multispecies_v2_global_clusters_v11.yaml` 双轨实验已完成：
  - 严格同源簇 development 轨：5 folds
  - 文献同口径随机蛋白轨：6 models × 10 seeds（含 Sul-BertGRU）
- 文献轨道结果：`structure_ranker` 10-seed 平均 macro AP = 0.386，显著高于每轮最优 baseline（0.049），`literature_track_leads = True`。
- Task 9.6 claim gate 当前输出：`literature_comparable_within_dataset_improvement_only`。
- 严格轨道尚未解锁 frozen test；番茄盲法队列尚未启动。

## 1. 全局约束

1. **禁止反向开发**：任何盲法数据查看之后，不得再调整模型、特征、Top-K 或阈值。模型冻结必须在 Gate 3 之前完成。
2. **预注册优先**：盲法主终点、SAP（Site-Aware Prioritization）协议、成功标准必须在数据产生/解盲前写成不可静默编辑的冻结文档。
3. **证据等级如实对应期刊**：NC 仅接受完整证据链；任何核心门槛失败即触发预注册降级，不得追加无边界模型开发来“抢救”。
4. **数据透明**：全部阴性结果、排除项、原始数据与中间产物保留并公开。
5. **Gate 2 措辞一致**：所有对外文本必须通过 `verify_predictive_claims` → `[]`；multispecies v2 轨道不得进入 Gate 2 证据。

## 2. 强制入场门与降级路线

| 门 | 内容 | PASS → | FAIL → |
|---|---|---|---|
| Gate 0 | 策略与预注册冻结：最终模型/特征/Top-K/对照/SAP/代码哈希冻结 | Gate 1 | 暂停，先完成冻结 |
| Gate 1 | Task 9 严格轨道一次性 frozen test + 完整审计 | Gate 2 | 只有内部 benchmark → 不以 NC 为目标 |
| Gate 2 | Gate 2 措辞与真实证据一致；NC 主结论仅锚定未来 lockbox | Gate 3 | 修正主张或降级 |
| Gate 3 | 番茄盲法主终点达到预注册标准；非单蛋白/单簇驱动 | Gate 4 | 盲法失败但资源完整 → Plant Physiology / Communications Biology |
| Gate 4 | ≥1 独立命中完成机制链：位点确认 → Cys 突变 → 功能效应 → 植物表型 → 救援 | Gate 5 | 盲法成功但机制不足 → Plant Physiology |
| Gate 5 | 若主张 PTM 串扰，必须有直接测量修饰间依赖关系 | Gate 6 | 删除串扰主张或补实验 |
| Gate 6 | 全部阴性结果、排除项与原始数据保留并公开 | Gate 7 | 补全数据档案 |
| Gate 7 | Nature-style 模拟审稿无未解决致命问题 | Gate 8 | 修正或降级 |
| Gate 8 | Gate 5 达到 READY；Gate 2 状态与措辞最终一致 | **投稿 NC** | 按 FAIL 门降级 |

### 降级期刊矩阵

| 已满足证据 | 目标期刊 |
|---|---|
| 全部 Gate 0–8 | Nature Communications |
| Gate 3 盲法富集成功，但 Gate 4 机制链不足 | Plant Physiology |
| Gate 3 盲法未成功，但 Gate 0–1 + 资源/方法完整 | Plant Physiology 或 Communications Biology（诚实方法/资源主张） |
| 仅 Gate 0–1（内部 benchmark） | 不以 NC 为目标；按实际证据投稿 |

## 3. 任务拆解

### Task 1：Gate 0 — 投稿策略与预注册冻结包

**负责人**：我们起草，张华确认。

- [x] 写定 `docs/superpowers/plans/2026-08-13-nc-entry-gate-roadmap.md`（本文件）并锁定版本。
- [x] 明确最终候选模型选择规则：
  - 从 v11/v12/v13 结果中综合严格轨道与文献轨道选出；
  - 选择标准写入冻结文档，不得事后调整。
- [ ] 生成 `results/candidates/candidate_release_v1/`：
  - [ ] 冻结模型权重与推理脚本哈希；
  - [ ] 特征模式与 Top-K 候选位点列表；
  - [ ] matched controls；
  - [ ] SAP 协议；
  - [x] `manifest.json` 记录当前已冻结的元数据与 SHA256（见 `results/candidates/multispecies_v2_candidate_release_v1/`）。
- [ ] 预注册番茄盲法队列：
  - 主终点定义（如 Top-K 富集倍数、FDR、跨蛋白/跨簇稳健性）；
  - 样本量与功效分析；
  - 解盲条件与标签隐藏机制；
  - 统计分案计划（命中定义、预注册检验、事后不加样本）。
- [ ] 全部冻结文档 git tag + SHA256 记录；冻结后零重排。

**PASS**：所有盲法数据产生前，模型、特征、Top-K、SAP、代码哈希已冻结并审计通过。

### Task 2：Gate 1 — Task 9 严格轨道 frozen test

- [x] 准备 test unlock 文件并按授权流程解锁 frozen test。
- [x] 运行 `scripts/score_v2_frozen_test.py --test-unlock <path>`（实际打分入口；`train_multispecies_v2.py` 承担面板准备）。
- [x] 生成 strict-track per-species / pooled AP、物种级 delta、`strict_track_bootstrap_delta`（10,000 replicate，实际 seed 20260813，与冻结发布一致）。
- [x] 通过 `audit_v2_run_manifest` 完整审计 root + literature manifests。
- [x] 更新 Task 9.6 claim gate：strict CI lower > 0 但 tomato delta 为负 → 升级条件不满足，维持 `literature_comparable_within_dataset_improvement_only`。

**PASS**：frozen test 一次性完成，manifest 审计通过，统计报告与 claim 明确。

### Task 3：Gate 2 — Gate 2 措辞与证据一致性审查

- [x] 运行 `verify_predictive_claims` 扫描全部展示文本，确保 `[]`。
- [x] 确认 `tests/release/test_gate2_ignores_within_dataset_split_metrics.py` 仍然全绿。
- [x] 更新 `docs/phase_f_gate2_decision.md`：明确 NC 主结论只能引用未来 lockbox；v11 内部轨道仅作数据集内证据。

**PASS**：没有任何 within-dataset 数字被误标为跨研究/跨物种证据。

### Task 4：Gate 3 — 番茄盲法主终点

**负责人**：张华执行实验；我们按冻结模型打分。

- [ ] 接收独立番茄盲法样本，记录 provenance 与批次。
- [ ] 使用 Gate 0 冻结的模型/特征/阈值进行打分；打分脚本哈希核验。
- [ ] 按预注册 SAP 进行解盲与统计检验。
- [ ] 敏感性分析：排除单一蛋白或单一同源簇驱动。

**PASS**：主终点达到预注册成功标准且稳健；**FAIL**：触发降级矩阵。

### Task 5：Gate 4 — 机制证据链

**负责人**：双方；张华实验室主导湿实验。

- [ ] 从盲法阳性命中中挑选 ≥1 独立命中。
- [ ] 完成证据链：
  1. 位点确认（质谱或抗体）；
  2. Cys 突变（Cys→Ser/Ala）；
  3. 功能效应（酶活、结合等）；
  4. 植物表型（成熟、胁迫等）；
  5. 救援实验（回补野生型 Cys）。

**PASS**：至少一条完整机制链；**FAIL**：盲法成功但机制不足 → Plant Physiology。

### Task 6：Gate 5 — PTM 串扰直接测量（若主张）

- [ ] 若稿件中主张 persulfidation 与 phosphorylation / 其他 PTM 串扰：
  - 设计直接测量实验（如修饰位点互赖、竞争/促进关系）；
  - 参考 `docs/superpowers/plans/2026-07-28-ptm-crosstalk-grammar.md` 的措辞规则。
- [ ] 若不主张串扰，则删除相关措辞，避免 overclaim。

**PASS**：串扰主张有直接测量支持，或无串扰主张。

### Task 7：Gate 6 — 数据透明与阴性结果档案

- [ ] 建立 `results/data_archive/nc_entry_gates/`：
  - 全部原始数据、中间产物、排除样本清单；
  - 阴性结果与失败命中记录；
  - 可复现脚本与 environment lock。
- [ ] 生成公共 DOI 或数据仓库链接。

**PASS**：任何审稿人可在 24 小时内定位到任意数字的原始来源。

### Task 8：Gate 7 — Nature-style 模拟审稿

- [ ] 按 NC 标准组装投稿包：图表、Source Data、Nature ML checklist、Reporting Summary。
- [ ] 组织模拟审稿：识别致命缺陷（样本量、机制、统计、可复现性等）。
- [ ] 对致命缺陷要么解决，要么诚实降级。

**PASS**：模拟审稿无未解决致命问题。

### Task 9：Gate 8 — 最终路径决策

- [ ] 复核 Gate 5 十条件全部 READY。
- [ ] 复核 Gate 2 状态与最终措辞一致。
- [ ] 决策：投稿 NC 或按降级矩阵执行。

## 4. 关键路径与并行工作

- **立即并行**：Task 1（冻结包起草）∥ Task 3（Gate 2 措辞审查）∥ Task 7（数据档案骨架）。
- **关键路径**：Task 1 → Task 2 → Task 4 → Task 5 → Task 8 → Task 9。
- **关键路径上无模型开发**：模型选择/冻结在 Task 1 完成；后续只使用冻结模型。

## 5. 验证协议

- 每次冻结产物更新后运行：`pytest tests/unit tests/scientific tests/release -q`。
- 每次 manifest 生成后运行：`python -m plantpersulf.cli audit-registry` / `audit-files` / `audit-leakage`。
- 每次对外文本过 `verify_predictive_claims` → `[]`。
- 盲法数据解盲前，由两人独立核验模型/特征/阈值哈希与冻结 manifest 一致。

## 6. 版本控制

- 本文件为 living document；每完成一个 Task 后在本节追加带日期、负责人、PASS 证据的条目。
- 冻结产物修订必须带日期+哈希；禁止静默编辑。

### 更新日志

- 2026-08-13：根据当前 v11 双轨实验结果与张华合作策略，建立 NC 入场门计划；锁定 Plant Physiology 预注册降级路径。
- 2026-08-13（Task 1，Gate 0 冻结包）：
  - `structure_ranker` 序列化 API 完成（`fit_structure_ranker` → `StructureRankerBundle` → `score_structure_ranker_bundle`；save/load 原子化，往返位级一致；655 tests 全绿）。
  - 冻结模型已拟合并保存：389,609 开发位点（10 种子面板并集，面板 SHA256 `30e432fb…`）、种子 20260813、`model_weights/structure_ranker_bundle.pt`（SHA256 `ab8a0353…`）。
  - 番茄蛋白组 Top-K 候选已生成：179,736 位点（排除训练面板 72,051），分数区间 [0.0053, 0.3053]；966 对照对；生成清单 `top_k_generation_manifest.json`。
  - 特征模式 `feature_schema.json`、盲法打分脚本 `inference_script.py`（无训练代码，bundle 哈希核验）已入包。
  - `sap_protocol.json` / `analysis_plan.json` 草案已起草（含 K 功效表、解盲条件、排除规则、敏感性分析）；**待张华实验室共同签署**。
  - 措辞纪律：发布包全部文本 `verify_predictive_claims → []`。
  - 待办：联合签署 → 预注册提交 → git tag + SHA256 全包冻结。
  - 新增 `2026-08-13-model-improvement-register.md`：模型改进杠杆登记册（番茄结构覆盖/ESM/数据扩充等 7 项），**登记不实施**，冻结期内禁止改动；解锁条件 = Gate 3 盲法完成后的新一代模型 + 新预注册。
- 2026-08-13（Task 7，Gate 6 数据透明骨架）：`scripts/build_nc_data_archive.py` 生成 `results/data_archive/nc_entry_gates/`——注册输入清单（243 文件+SHA256）、实验产物哈希、发布包清单、排除账本（8 项）、阴性结果登记（7 项，永不删除）、pip freeze 环境锁与复现命令索引；全档案文本 `verify_predictive_claims → []`。待办：公共 DOI/仓库链接（随预注册一并提交）。
- 2026-08-13（投稿稿件骨架，用户决策"先按 Plant Physiology 写"）：`manuscripts/plant_physiology/2026-08-13_pp_manuscript_draft_v1.md` 建立——PP Research Article 格式完整骨架：R1–R5（数据集/文献同口径基准/冻结发布包/盲法前校准/预注册设计）按现有证据写就，R6–R8 为 Gate 3/4 占位块（含两个预写摘要结局变体，成功/不成功各一）；含 4 表 6 图例、方法学、数据可用性与参考骨架。全文通过 `verify_predictive_claims → []` 及第二措辞族（external generalization 等）零命中。填写依赖：联合签署+预注册 DOI → Gate 1 数字 → Gate 3 盲法结果 → Gate 6 公共 DOI。
- 2026-08-13（**Gate 0 冻结完成**）：双方联合签署确认（由建模方 PI 记录）→ `sap_protocol.json` / `analysis_plan.json` 升 1.0（`signed_2026-08-13`），`co_signature_package.md` 增 §6 签署记录 → 全包文件 SHA256 写入 `manifest.json`（status=`frozen`，pending 仅剩预注册 DOI）→ `SHA256SUMS` 刷新 → `audit_candidate_release.py --expect-signed` **PASS**（bundle 确定性、措辞门 `[]`、JSON 解析、签署状态全绿）→ git tag `gate0-release-v1-20260813`。**冻结后零重排**：任何变更须新发布号 + 新预注册。剩余唯一 pending：OSF 预注册提交（用户人工步骤）。
- 2026-08-13（盲法设计定稿 + 张华交接包）：用户决策按 600 位点方案执行（K=200 + 每候选 2 对照），`analysis_plan.json` 标记 `decided_2026-08-13_pending_lab_co_signature`。新增 `scripts/build_zhang_lab_handoff.py` 生成 `results/handoff/zhang_lab_blind_cohort_v1/`：盲态位点清单（593 行，blind_id + accession/基因/位置/±10aa 窗口，种子 20260813 确定性打散）、分组密钥（`DO_NOT_SEND`，建模方留存，SHA256 已入 `handoff_manifest.json` 供解盲核验）、中文 README（判定规则：confirmed 定义由实验室开测前书面确认；未检出 ≠ 未被修饰，一律计 not confirmed）。**队列实为 593 而非 600**：4 个候选在 0.25 SD 卡尺内不足 2 个对照（3 个 0 对照、1 个 1 对照）——特征极端位点无近邻，如实记录在 `analysis_plan.json` 的 `cohort_actual`，统计检验与功效不受影响（<1%）。交接纪律：湿实验方只收盲态清单 + 两份待签署文件 + 中文说明；完整候选表/分数/密钥不交予实验方。交接包与修订后分析计划均过 `verify_predictive_claims → []`，`audit_candidate_release.py` 全绿。
- 2026-08-14（**Task 2 完成，Gate 1 frozen test 一次性审计通过**）：unlock `v2_frozen_test_unlock_20260813.json`（code_revision `476faa9d…`，4 次重签均为文档/审计变更）→ `scripts/score_v2_frozen_test.py --test-unlock`（beflb50s2，exit 0）→ **manifest audited**（root manifest 全字段 + 双臂指纹/checkpoint/training 事件全覆盖）。结果：269,801 位点、566 阳性（arabidopsis 63 / rice 190 / tomato 20 / magnaporthe 293）；pooled bootstrap delta **+0.1402 [0.1183, 0.1699]**（10,000 聚类 bootstrap，seed 20260813）；per-species 模型 AP：arabidopsis 0.0996（基线 0.0015，125× base）、rice 0.1978（0.0023，133×）、**tomato 0.00032（0.00074，0.76×，低于 base rate）**、magnaporthe 0.0158（0.0264，0.82×）。三份分数文件与 run-2 记录 SHA256 逐字节一致（三次打分完全复现）。claim gate：tomato delta −0.0004 → **维持** `literature_comparable_within_dataset_improvement_only`（升级需三作物 delta 均为正），`gate2_eligible=False`；数字写入 `task9_6_claim_gate_report.md` §2/§4/§5，稿件 R2 以此定稿。
