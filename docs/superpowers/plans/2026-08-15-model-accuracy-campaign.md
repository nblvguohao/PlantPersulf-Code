# 模型准确度提升评估战役 — 数据结构 / 模型 / 已有工具三线计划

> **目标**：在**不触碰冻结包、候选表、盲法队列与已消费的 frozen test** 的前提下，
> 针对三个方向——①数据结构（番茄 AlphaFold 结构覆盖）、②模型（校准/选择、ESM 分支）、
> ③已有类似位点预测工具（Sul-BertGRU 堆叠、pCysMod 定性）——在**冻结内部轨**
> （文献同口径随机蛋白轨 + 严格同源簇 dev folds）上跑通评估并记录数字，
> 为 Gate 3 解锁后的新一代模型立项提供决策证据。
>
> **与现有文档的关系**：承接 `2026-08-13-model-improvement-register.md` 的
> "评估纪律"（§3：任何杠杆评估必须先在校内轨跑通并记录数字，再决定是否立项），
> 以及 `2026-08-13-nc-entry-gate-roadmap.md` 全局约束 1（禁止反向开发）与约束 3
> （不得追加无边界模型开发"抢救"）。本战役是**评估**，不是**实施**：
> 结果只写入登记册，任何模型/特征/Top-K 变更仍然走"Gate 3 解锁后新发布号 +
> 新预注册 + 联合签署"的唯一路径。

## 1. 动机与证据

- Gate 1 frozen test：arabidopsis 0.0996（125× base）、rice 0.1978（133×）、
  **tomato 0.00032（0.76×，低于 base rate，n=20）**、magnaporthe 0.0158（0.82×）。
  claim gate 因 tomato delta −0.0004 无法升级（升级需三作物 delta 全正）。
- kiae271 模拟盲测诊断（`scripts/evaluate_kiae271_with_release_bundle.py`）：
  0/99 位点进入候选表 Top-2000；20 个干净位点 65% 排在中位数以上（信号存在但弱）；
  原始分数分布极扁平，顶部被"特征极端但未必真实"的位点占据。
- 登记册优先级（`2026-08-13-model-improvement-register.md` §2.1）：
  L1（番茄结构覆盖）为最高杠杆；L8（分数校准/选择策略）次之；L2（ESM）/L5/L6 无新增证据。
- 用户决策（2026-08-15）：三线并进（W1+W2+W3），启动本战役。

## 2. 纪律边界（不可违反）

1. **不触碰冻结产物**：`results/candidates/multispecies_v2_candidate_release_v1/`
   （bundle `ab8a0353…`、候选表、matched controls、SAP/分析计划）、
   Gate 0 tag `gate0-release-v1-20260813` 均保持原样。
2. **不触碰盲法数据**：张华实验室 593 位点盲法队列按冻结模型继续，本战役的任何
   产物不得进入该流程（`results/handoff/zhang_lab_blind_cohort_v1/` 只读）。
3. **不再解锁 frozen test**：`v2_frozen_test_unlock_20260813.json` 为一次性授权，
   已消费。本战役全部评估在**文献轨（10 seeds）与严格轨 5-fold dev** 上完成。
4. **不做反向开发**：本战役不看盲法数据、不用 frozen test 调任何东西；
   评估数字只是登记册证据。
5. **阴性结果保留**：所有评估结果（含阴性）写入登记册更新日志，不删除。

## 3. 评估协议（三线共用）

- **内部轨口径**：与 v11 完全一致 ——
  - 文献轨：`literature_random_protein`（test 0.2 / val 0.2-of-remaining，seeds 0–9），
    主指标 `macro_species_average_precision`（物种内 AP 的宏平均）。
  - 严格轨：`multispecies_strict_v4.tsv` 5-fold dev，per-species AP，
    折叠间 bootstrap（seed 与 v11 一致）。
- **对照基线**：v11 已记录数字（structure_ranker 文献轨 0.3863±0.1727；
  严格轨 dev fold val AP 见 `task9_6_claim_gate_report.md` §2）作为对照。
  评估任何杠杆时**模型超参冻结为 v11 值**（hidden 16 / dropout 0.2 / epochs 200 /
  lr 0.05 / use_esm 0 / use_study_context 按 v11），只变化被测变量。
- **番茄为约束指标**：三作物 per-species AP 均报告；tomato 的提升是战役核心目标
  （它是 Gate 2 升级与盲法富集的绑定约束）。
- **统计**：文献轨 10 seeds 报 mean±std 与 per-seed 对比（配对方差）；
  严格轨沿用 v11 的 bootstrap 口径。

## 4. W1 — 数据结构：番茄 AlphaFold 结构覆盖（L1，最高杠杆）

- **假说**：结构分支对番茄完全掩蔽（registry 0 条番茄结构）是番茄排序退化的
  主因之一；补齐番茄结构特征后，tomato per-species AP 应向 arabidopsis/rice 方向恢复。
- **方法**：
  1. 枚举内部轨面板（文献轨 10-seed 面板 + 严格轨 dev folds）中的番茄蛋白 accession；
  2. 用 `scripts/download_alphafold_structures_bulk.py`（fail-closed，404=真实缺失
     记录）下载 AFDB 结构（AFDB 含 Solanum lycopersicum, UP000004994）；
  3. 生成 registry release v3 tsv（`download_alphafold_structures_bulk.py` 输出格式，
     SHA256 齐全）并过 `audit-registry`；
  4. 重建 contact_number_proxy + pLDDT 结构特征行（只重建番茄行，其余物种行沿用 v11 缓存）；
  5. 相同冻结超参重跑文献轨 structure_ranker（10 seeds）与严格轨 dev folds；
  6. 记录 per-species AP delta（尤其 tomato）与结构分支 gate 权重变化。
- **产出**：`results/experiments/accuracy_campaign/w1_tomato_structures/`
  （registry tsv、特征 manifest、双轨分数、对比报告）。
- **风险**：番茄阳性少（79 dev / 20 frozen），增益量级可能有限；AFDB 可能对部分
  番茄 accession 无模型（404 记为真实缺失，不算失败）。

## 5. W2 — 已有工具：Sul-BertGRU 堆叠（近零成本）

- **现状**：Sul-BertGRU 已移植（`competitors/sul_bertgru.py`），文献轨 10 seeds 分数
  已落盘（`literature_random_protein/sul_bertgru_seed{0..9}/scores_batch128.tsv`，
  partition 含 test/validation；`shared_panel.tsv` 含 train/validation/test）。
  单独跑弱（macro AP 0.0213），但其分数作为集成输入从未评估。
- **方法**：
  - **方案 A（秩聚合，训练自由）**：structure_ranker 与 sul_bertgru 的 per-seed
    test 分数做中位秩/Borda 聚合（复用 `evaluation/rank_ensemble.py` 口径），
    与单模型对比 10-seed macro AP。零新训练。
  - **方案 B（堆叠特征）**：sul 分数作为 structure_ranker 的额外输入。
    前提：train 分区行也要有 OOF sul 分数（当前只有 test/validation 有分）。
    可行性：`.bert_windows.npy` 已缓存 BERT embedding，只需对 train 分区交叉拟合
    GRU 头（`fit_sul_bertgru` 已支持任意 train_indices/score_indices）。
    若成本允许则执行；否则方案 A 结论先行。
- **pCysMod**：维持 `docs/competitor_data_audit_2026-08-10.md` 结论——训练数据
  不可获取、仅 web server → **定量互跑不可行**，本战役只做定性记录（v11
  `comparator_policy.pcysmod` 已定义姿态），不安排 web server 批量提交。
- **产出**：`results/experiments/accuracy_campaign/w2_tool_stacking/`。

## 6. W3 — 模型：校准/选择（L8）+ ESM 分支（L2）

- **L8 校准/选择策略评估**（针对 kiae271 诊断"信号在中段非顶端"）：
  - 在文献轨 structure_ranker 分数上评估：MC-dropout 不确定度引导选择
    （排除高不确定的极端位点）、分数带分层抽样、校准后验重排；
  - 指标：与原始分数对比 per-species AP@K 与召回@K（K=50/200，同 claim gate 口径）；
  - **只评估选择规则，不改任何冻结产物**（任何 Top-K 变更 = 新预注册，登记册 §L8）。
- **L2 ESM 分支**：ESM 特征已缓存（`data/processed/features/multispecies_v2_esm2_t33_windows_v1/`），
  use_esm=1 重跑文献轨 structure_ranker（其余超参冻结），与 v11 use_esm=0 对比。
- **产出**：`results/experiments/accuracy_campaign/w3_calibration_esm/`。

## 7. 决策规则与交付

- 每条杠杆的"通过"定义：文献轨 10-seed macro AP 或 tomato per-species AP 出现
  **配对方差显著且方向为正**的改善（预注册在对应报告内写死统计口径），
  且严格轨 dev folds 不恶化。
- 全部数字 → `2026-08-13-model-improvement-register.md` 更新日志
  （`## 2.2 评估战役 2026-08-15`），含阴性结果。
- 立项决策仍守：Gate 3 盲法验证完成后 → 新发布号 + 新预注册 + 联合签署。
  本战役结束后生成《立项决策摘要》，列出各杠杆的证据强度、预期收益、代价。
- 验收：新增/改动代码走 TDD（tests/unit 全绿）；全部对外文本
  `verify_predictive_claims → []`；registry/特征产物过 `audit-registry`/`audit-files`。

## 8. 版本控制

- 本文件为 living document；每条工作流完成后在更新日志追加日期、负责人、数字、结论。
- 评估产物目录 `results/experiments/accuracy_campaign/` 不入冻结发布包，标注
  `INTERNAL EVALUATION ONLY`。

### 更新日志

- 2026-08-15：战役立项。背景：用户（建模方 PI）提出"针对数据结构、模型、已有类似
  位点预测工具提升预测准确度"；经 `AskUserQuestion` 确认三线并进（W1+W2+W3）。
  登记册纪律（评估先行、立项守 Gate 3）不变。
- 2026-08-15（基础设施）：`scripts/accuracy_campaign/` 建成并验证——文献轨逐位点重跑
  runner（panel_sha256 与冻结 summary 10/10 一致、test AP 逐位复现）、W2 堆叠评估、
  W3 校准评估、严格轨 fold 评估脚本；unit tests 全绿（test_stacking /
  test_accuracy_campaign_helpers）。
- 2026-08-15（W2 完成，阴性）：Sul-BertGRU 堆叠无增益——等权秩聚合 macro AP
  0.3863→0.1181；验证集最优凸组合=忽略 sul（10/10 seeds alpha≈0.05，sul 分数因
  sigmoid→softmax 头部近乎退化）。详见 `results/experiments/accuracy_campaign/w2_tool_stacking/report.md`。
- 2026-08-15（W3-L8 完成，一阴一阳）：σ 门控不成立（丢 top10% 高 σ 行即丢 ~60%
  阳性，hits@50 24.7→11.3）；分物种 isotonic 校准正面（hits@50 +21%、recall_full@200
  +17%，增益集中在弱 seeds，无灾难性恶化）——列为新发布候选变更（仍守 Top-K 预注册纪律）。
  详见 `results/experiments/accuracy_campaign/w3_calibration_esm/report.md`。
- 2026-08-15（W3 完成，L2 阴性）：use_esm=1 重跑 macro AP 0.1988±0.2006 vs v11
  0.3863±0.1727，10/10 seeds 劣于或接近 v11，呈双峰不稳定 → L2 不立项，印证 v11 消融
  选择 use_esm=0。W3 全部结论见 `results/experiments/accuracy_campaign/w3_calibration_esm/report.md`。
- 2026-08-15（W1 数据就绪）：29,767 番茄 accession 下载完成——**27,542 下载 /
  1,105 not_found / 0 错误**（6.6 GB）；修复两处生产缺陷（`--workers` 并行批次、
  API 404→not_found）；release v3 快照完成（`alphafold_structures_release_v3.tsv`，
  30,681 行全量 SHA256 审计通过，SHA256 `60aa3af0…`）并注册 `model_inputs.tsv`；
  严格轨 5-fold 基线臂完成（tomato 平均 val AP 0.0413，macro 0.5383）。
- 2026-08-15（**W1 完成，L1 阴性**）：文献轨 macro 0.3863→0.2117（tomato 0.0185→0.0204
  无增益）；严格轨 macro 0.5383→0.3257（**tomato 0.0413→0.0215，−48%**）。机制线索：
  结构 scaler 被番茄行主导 + 弱监督过拟合；分物种 scaler 等改进留作新代设计候选。
  详见 `results/experiments/accuracy_campaign/w1_tomato_structures/report.md`。
- 2026-08-15（**战役完成**）：全部杠杆评估完毕（一阳四阴）。唯一正向杠杆 = 分物种
  isotonic 校准（新发布候选变更）；`DECISION_SUMMARY.md` 定稿；登记册 §2.2 同步。
