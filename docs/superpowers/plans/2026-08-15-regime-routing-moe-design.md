# 体制路由 MoE v2 设计文档（预注册草案 + 评估层诊断）

> 状态：**设计 + 评估层诊断阶段**（Gate 3 盲法验证完成前，不训练、不改冻结包）。
> 关联：`2026-08-13-model-improvement-register.md`（§2.2 / §L6）、
> `2026-08-15-upgrade-iteration-diagnostics.md`（§3.2 / §9.2-9.6）、
> `植物巯基化位点筛选_方法设计.md`（命题三，§3 / P6）、
> `2026-08-15-model-accuracy-campaign.md`（W1 结构 scaler 线索）。
> 本文件为 living document；诊断结论与预注册定稿在更新日志追加。

## 1. 目标与纪律

**目标**：把 v2 结构分支从"全局加性门控融合"改为**体制路由**——由 pLDDT
（AF 无序度代理）等结构上下文驱动 gate，把位点路由到体制专家（IDR / 折叠 /
金属簇），每专家用自己的特征方向评分。命题三（方法设计 §3）认为**单一
分类器必然在两个体制间取折中，两边都做不好**；v2 已收集的证据支持这一方向。

**硬纪律（本阶段不可违反）**：
1. **不触碰冻结包**（bundle `ab8a0353…`）、候选表、盲法队列
   （`results/handoff/zhang_lab_blind_cohort_v1/` 只读）、已消费 frozen test。
2. **任何新架构训练 = 新发布号 + 新预注册 + 联合签署**。本阶段只做
   （a）评估层诊断（冻结内部轨同口径，记录数字再定立项）、（b）预注册草案。
3. 评估先用**文献轨 10-seed 同口径**（`scripts/accuracy_campaign/run_literature_ranker_scores.py`，
   `panel_sha256` 与冻结 summary 一致）跑通并记录，再决定是否立项。
4. **不在盲法数据上做任何反向开发**；不在冻结轨外另起训练。

## 2. 证据基础（为什么是体制路由）

| 来源 | 证据 | 方向 |
|---|---|---|
| §3.2 体制假说 | SlWRKY6（IDR）真位点局部密度 0.048 < 干扰 0.095；BRG3（金属簇）真位点 0.238 ≥ 均值 0.158 | 两典型体制局部密度方向**均符合** → 路由而非平均 |
| §9.2-9.3 n=12 体制内检验 | folded 体制真位点系统性更埋藏（contact hit@2 7/12）；方向随体制翻转（folded 用埋藏，disordered 用暴露） | 体制内带符号排序有效（ceiling） |
| §9.4 LOO 修正 | **contact_number_10a 唯一显著**（protein p=0.018、regime-local p=0.032）；7 特征等权带符号复合不泛化（p≈0.3）；disordered 暴露子主张不显著（n=3 过小） | 复合被 6 个噪声特征稀释；contact 是首要特征 |
| W1 番茄结构（L1 阴性） | 结构 scaler 被番茄行主导（v3 后番茄占结构行 ~78%），arabidopsis/rice 结构输入被重缩放 → 两作物大幅退化 | **按体制/按物种标定结构输入**是 W1 之后的第一候选修复 |
| 共肽 9.5-9.6 零结果 | 共肽内正/负结构不可分辨（3 物种 14 肽段，全 p>0.05） | 路由只服务**蛋白内排序**，不承诺共肽内分辨 |

**诚实局限**（不可逾越）：
- 监督树在 labeled n=12 下**不可估**（多次记录）；体制分桶后每桶 labeled 更少
  （disordered 仅 ~3）。任何"从 labeled 学到的方向"在 LOO 下都不显著
  （§9.4 已证）。
- 因此本设计把**路由的方向先验写死为结构假说**（不是从 labeled 学），用
  置换检验评估；真正的端到端 MoE 训练依赖 L0/MIL 级标签（待建）。

## 3. 架构设计（v2 候选，预注册端点）

```
        pLDDT / RSA / contact 上下文
                    │
              gate（软路由）
              │    │    │
        disordered  folded  RING/金属簇
        专家       专家      专家
        （暴露/孤立）（埋藏/堆积）（簇内暴露 gate）
              │    │    │
              └─ soft 加权融合 ──→ 结构分支分数
```

- **gate 语义**：pLDDT 分桶硬路由（`regime_bucket`，folded≥70 / 50≤linker<70 /
  disordered<50）作 base；软路由（RSA/contact 连续上下文）为候选扩展。
  gate 的**输入只含结构上下文，永远不含位点标签**（结构性先验，非 learned）。
- **专家方向**（写死为结构假说，非 labeled 学习）：
  - folded：埋藏/堆积方向（contact_number 高 = 更像功能位点，§9.4 唯一显著特征）；
  - disordered：暴露/孤立方向（RSA 高、最近 Sγ 远）——**假设**（§9.4 排列下
    n=3 不显著，标记为待更数据）；
  - RING/金属簇：folded 的亚专家，额外簇内暴露 gate（方法设计 §1.3 游离暴露
    Cys 读数，P2 实证）。
- **与冻结模型关系**：不替换冻结 bundle。v2 结构分支若立项，作为**新发布
  模型**的一部分（`build_candidate_release_v2`），冻结 v1 保持现状。
- **缺失结构掩蔽**：沿用冻结模型语义——缺 PDB 行结构分支输出强制为零，
  绝不均值插补（不可变承诺）。

## 4. 评估层诊断（本阶段可做，不训练）

三轨并行，全部在冻结内部轨 / 已知对照上：

**A. 路由 vs 单一 contact（已知对照 n=12，LOO + 置换）**
在 `structure_regime_loo_v1.json` 基础上，直接比较
"contact_number 单特征蛋白内排序" vs "按 pLDDT 体制路由后各专家自评"
的总负担与排列 p 值。**成功标准（诊断级）**：路由负担 ≤ 单特征负担，
且路由的置换 p ≤ 单特征 p；反之记录为"路由在当前数据上不优于 contact
单特征"（阴性也是结论，登记不立项）。

**B. 体制感知结构标定（文献轨，W1 机制线索）**
针对 W1 "结构 scaler 被番茄行主导"的机制线索：在文献轨同口径上，对比
"全局结构标定" vs "按体制（pLDDT 桶）分桶标定"的逐位点分数差异与 per-species
AP@K/召回@K（K=50/200，同 claim gate 口径）。**成功标准**：分桶标定使
arabidopsis/rice per-species AP 改善且 tomato 不灾难性恶化（配对方差显著为正）。
仅评估标定规则，不改任何冻结产物（结构输入缩放属于特征投影，仍记为新发布候选）。

**C. 体制分布审计（全蛋白组，数据面）**
统计冻结候选表 / 全蛋白组 Cys 在 three-regime 桶的分布、contact 特征的中位
分位数——确认路由在**推断侧**可用（gate 输入在候选表上全覆盖，无 missing
结构导致的路由饥饿）。

## 5. 预注册端点（若诊断 A/B 正向，正式立项前定稿）

- **统计口径写死**：文献轨 10-seed 配对方差（Wilcoxon signed-rank）显著为正的
  改善 + 严格轨 5-fold dev 不恶化 + 已知对照 hit@2 不降；K=50/200 同 claim gate。
- **成功定义**：路由 MoE ≥ 当前 v1 结构分支在全部三轨不劣，且至少在
  arabidopsis/rice 之一配对方差显著改善。
- **承诺**：立项后先发布预注册（联合签署），再训练；训练产物全新发布号；
  冻结 v1 候选表/盲法流程不受影响。

## 6. 实施清单（按顺序）

1. `src/plantpersulf/evaluation/regime_moe.py`（纯函数，TDD）：
   `moe_route_score`（gate 路由 + 专家自评 + soft 加权）、
   `routed_burden`（路由负担 vs 随机）、`loo_routed_burden` + `permutation_null`。
   → 诊断 A。
2. `scripts/evaluate_regime_moe.py` → `results/diagnostics/regime_moe_v1.json`。
3. 诊断 B：文献轨体制标定对比（复用 accuracy_campaign runner）。
4. 诊断 C：体制分布审计脚本。
5. 全量测试 + register §4 更新日志 + memory 同步 + commit。
6. 诊断结果（正/负）进 register 与《立项决策摘要》；若正向，起草正式预注册。

## 7. 更新日志

- 2026-08-15：设计文档建立。背景：用户指令"继续推进体制路由 MoE 实现"；
  梳理 register §2.2/L6、诊断 §3.2/§9.2-9.6、方法设计命题三、W1 scaler 线索后
  定稿纪律边界与三轨诊断方案。本阶段不训练、不改冻结包。
- 2026-08-15：**诊断 A 完成（阴性）**。新增 `evaluation/regime_moe.py`
  （+6 单元测试，专家方向写死为结构先验、非 labeled 学习，故无 LOO 学习偏差）
  + `scripts/evaluate_regime_moe.py` → `results/diagnostics/regime_moe_v1.json`
  （n=12 已知对照，B=999 置换，seed 20260815）：

  | 变体 | observed | random | nullMed | p |
  |---|---|---|---|---|
  | contact protein | 36 | 58.5 | 59 | 0.016 |
  | contact regime | 30 | 58.5 | 44 | 0.044 |
  | routed protein | 38 | 58.5 | 58 | 0.030 |
  | routed regime | 31 | 58.5 | 44 | **0.057** |

  **结论：MoE 路由（当前先验方向）不优于 contact 单特征。** regime 模式
  几乎持平但 p 从 0.044 恶化到 0.057（不再显著）；protein 模式略差
  （38 vs 36）。机制：folded 专家的 contact 方向承载了全部信号（9 个
  folded 对照 6 个 rank≤2）；disordered 专家的 RSA/Sγ 方向（n=3，含
  A0MES8 rank 3/3）无增益且轻微稀释——与 §9.4 的 "RSA 方向排列下不显著"
  一致。**命题三的路由方向正确，但当前特征分辨率下 disordered 专家没有
  可路由的显著特征**：路由无物可路由。MoE 架构在 labeled n=12 + 7 特征下
  不立项。
- 2026-08-15：**诊断 B 状态调整**。W1 的"全局结构 scaler 被番茄主导"是
  **分物种/分体制标定**问题，与 MoE 路由正交；诊断 A 表明体制方向先验无法
  带来增益。分桶标定需要完整文献轨重训，列为新代设计候选（Gate 3 后，
  若标定轨有独立证据再立项），本阶段不重训。
- 2026-08-15：**诊断 C（体制分布审计）**。已知对照 n=12 体制分布：
  folded 9 / disordered 3 / linker 0（本诊断确认 gate 输入在对照集上无
  linker 饿死风险；全蛋白组分布审计待候选表快照延伸）。
