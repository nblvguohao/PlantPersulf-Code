# Phase P — 过硫化×磷酸化串扰语法（PTM Crosstalk Grammar）

**日期:** 2026-07-28
**状态:** 待项目评审人批准；批准后按任务逐个执行，每个任务单独评审
**父目标:** 在 `GATE2_STOP` 结构性锁定的前提下，产出一个不依赖跨研究预测能力的、
可发表的机制性发现

---

## 1. 为什么做这个

Phase Z 已经确立了两层跨物种关联性证据（§5.5.1 家族级 PANTHER 富集、§5.5.2
结构级包埋趋同），但 v10 的新颖性核查诚实地把可主张的新颖性收窄为
"**跨物种统计趋同性**"——因为"氧化还原/中心碳代谢富集"这一定性格局在
Aroca et al. 2017 单物种层面已被报道。

也就是说：**目前的发现是"更严格地确认了别人已经定性观察到的东西"，
不是一个新机制。** 这不足以支撑 Plant Communications 以上量级。

本 Phase 提出一个**尚未被任何人系统检验过的机制假说**，且该假说：

- 不需要跨研究预测能力（绕开 `GATE2_STOP` 的结构性阻塞）；
- 使用的统计设计与 §5.5.2 已验证的"同蛋白内对照"完全同构（复用已测试代码）；
- 有两个来自张华团队的**突变体级锚定案例**，不是凭空假设。

### 假说的来源（真实文献观察，非推测）

| 蛋白 | 过硫化位点 | 磷酸化位点 | 激酶 | 已报道的相互关系 |
|---|---|---|---|---|
| SlWRKY6 | Cys396 | Ser33 | SlMAPK4 | 过硫化**削弱** MAPK4 介导的磷酸化（doi:10.1093/plphys/kiae271）|
| SlERF.D2 | Cys35 | Ser42 | SlMAPK4 | 过硫化削弱转录活性，磷酸化增强（doi:10.1111/tpj.70000）|

两个独立蛋白、两篇独立论文、同一个激酶、同向的拮抗关系。SlERF.D2 的
Cys35 与 Ser42 在一级序列上**相距 7 个残基**。

> **假说 H:** 过硫化的半胱氨酸在序列和/或三维空间上系统性地邻近磷酸化位点，
> 即 H₂S 信号通过在 Cys 上"占位"来干扰邻近 Ser/Thr/Tyr 的激酶可及性。

### 必须先声明的反向预期（防止确认偏误）

§5.5.2 已确立：过硫化 Cys **系统性更包埋**（三物种方向一致，2/3 显著）。
而磷酸化位点通常位于**表面/无序区**。

因此结构层面（H2）的先验方向其实是**相反的**——包埋的 Cys 应该离表面磷酸化位点
更远。这意味着：

- **H2 若成立，是一个逆先验的强结果**；
- **H2 若为零结果，高度可能只是包埋效应的机械后果，不能解读为"假说被证伪"**。

这一非对称性必须在设计阶段写死，不得在看到结果后再补充解释。

---

## 2. 科学边界（硬约束）

本 Phase **不能**：

- 改变 `configs/gate2_v1.yaml` 或 `studies_are_independent: false`；
- 把 `GATE2_STOP` 改为任何其他判定；
- 用番茄的四个已登记控制位点做**统计输入**（番茄无 persulfidome，它们只能作为
  案例叙述，见 Task X6）；
- 把"未观测到过硫化的 Cys"当作实验阴性——全流程为正例-未标注（PU），
  输出字段与文档措辞必须用 `unlabeled`，禁止出现 `negative`；
- 在看到检验结果后调整零模型、分层变量或显著性阈值。

沿用的降级措辞不变：

> 当前公开数据不足以证明跨研究预测能力，模型仅用于候选组织与假设生成。

---

## 3. 前置任务 X0：本地数据可复现性修复（**阻塞性，必须先做**）

### 已发现的严重缺口

复核发现**注册表记录的资产与本地磁盘实际内容严重不一致**：

| 资产 | 注册表声明 | 本地磁盘实际 | 状态 |
|---|---|---|---|
| AlphaFold 结构 | 2,006 条（`alphafold_structures.tsv` 2007 行）| `data/raw/alphafold/` **18 个文件** | ❌ 缺失 |
| 同源聚类 v2（MMseqs2 真实簇）| §3.3 声明 13,967 簇，SHA256 `e13d16ae…` | `data/processed/clusters/` **只有 v1 singleton** | ❌ 缺失 |
| 水稻位点表 PXD072089 | §5.3 声明 929 个坐标核验位点 | 磁盘上**无任何 072089 文件** | ❌ 缺失 |
| 稻瘟菌位点表 PXD063170 | §5.2 已解析 | 磁盘上**无文件** | ❌ 缺失 |
| 拟南芥 benchmark_v1 | 395,879 行 | ✅ 存在 | OK |
| Figure 1 | 三面板已发布 | ✅ svg/pdf/png 存在 | OK |

**这意味着 `docs/phase_z_evidence_audit.md` §5.5 的全部数字目前在本机无法重跑。**
图已经画出来了，但支撑它的中间产物不在。

这是**谈合作之前必须堵上的信誉风险**：如果张华老师团队要求复现，或者审稿人要求
提供数据，现在无法响应。

### X0 交付要求

按 `docs/data_reproduction.md` 的清单与 `data/registry/reproduction_downloads.tsv`
重新获取，并逐项 SHA256 比对注册表：

1. 2,006 个 AlphaFold v6 结构；
2. PXD072089 六份 PNAS 补充表 + `SS-all-peptides.tsv`，重跑解析器，
   **断言位点数 == 929**；
3. PXD063170 位点表，重跑解析器；
4. `protein_clusters_v2.tsv`（若 A100 不可用，按
   `docs/ops/remote_a100_via_lab_jump_host.md` 重新走两跳 SSH；
   **禁止**用 v1 singleton 替代并声称是 v2）；
5. 重跑 §5.5.1 / §5.5.2 全部检验，**逐位比对已发布的 p 值**
   （拟南芥×水稻、水稻×稻瘟菌、拟南芥×稻瘟菌；结构层三物种）。

**验收:** 任何一个数字对不上，立即停止本 Phase 并开 issue。
不得"就近取整"或"方向一致即可"。

---

## 4. 数据获取任务

### Task X1 — 磷酸化位点数据源预检与登记

沿用 `src/plantpersulf/evidence/preflight.py` 的既有模式：**先审计，再使用**。

**候选源（按可溯源性排序）:**

| 物种 | 主源候选 | 交叉校验源 |
|---|---|---|
| 拟南芥 | PhosPhAt 4.0 批量下载 | UniProt `ft_mod_res`（已有 REST 基础设施）|
| 水稻 | 已发表的水稻磷酸化蛋白组 PRIDE 数据集 | UniProt `ft_mod_res` |

**RED:** `tests/scientific/test_phosphosite_source_registry.py`

- 断言 `data/registry/phosphosite_sources.tsv` 存在，且每行有
  `source_accession / species / file_name / sha256 / retrieved_at / data_level / scientific_use`；
- 断言任何未登记来源被 `test_no_unregistered_scientific_inputs.py` 的既有机制拒绝
  （扩展该测试，而非新写一套）。

**GREEN:** 新增 `configs/phosphosite_sources_v1.yaml` + 复用
`src/plantpersulf/download/registered.py` 下载路径。

**预检产物:** `src/plantpersulf/ptm/phosphosite_content.py` —— 输出实际观测到的
列名、行数、物种、位点数、accession 命名空间。**覆盖率不足则本 Phase 降级**
（判定标准见 Task X5 的 `insufficient_coverage`）。

---

### Task X2 — 磷酸化位点解析为规范表 + 坐标核验

**规范 schema**（与 `BENCHMARK_SITE_FIELDS` 风格一致）:

```
protein_accession, residue_position, residue_aa, source_accession,
evidence_level, source_sha256
```

**RED:** `tests/scientific/test_phosphosite_coordinates_match_sequence.py`

直接复用 `test_site_coordinates_match_sequence.py` 的既有断言模式：
每个 `residue_position` 在参考蛋白组 FASTA 中的实际残基必须等于
`residue_aa` ∈ {S, T, Y}。不匹配的行进入 `conflicts` 报告并被排除，
**不得静默丢弃、不得猜测偏移**。

fixture 必须是真实文件切片（受 `test_fixture_provenance.py` 约束）。

---

## 5. 分析任务

### Task X3 — 一级序列串扰（H1）

**统计设计（同蛋白内对照，与 §5.5.2 同构）:**

对每个满足以下条件的蛋白：至少 1 个过硫化 Cys、至少 1 个磷酸化位点、
至少 2 个 Cys——

1. 对该蛋白的**每个** Cys，计算到最近磷酸化位点的序列距离 `min|pos_cys − pos_phos|`；
2. 计算过硫化 Cys 在该蛋白全部 Cys 的距离排序中的**归一化秩** ∈ [0,1]；
3. 零假设下期望均值 = 0.5。

**为什么用同蛋白内对照:** 自动控制蛋白丰度、检测偏倚、蛋白长度、
磷酸化位点数量、蛋白级无序度——这些正是 PTM 共定位分析最容易翻车的混杂。

**检验:** 复用 `src/plantpersulf/evaluation/permutation.py`；
效应量复用 `effect_size.py`；置信区间复用 `bootstrap.py`。
**不引入 scipy**（与仓库既有从零实现的统计保持一致）。

**RED:**
- `tests/unit/test_ptm_crosstalk_rank.py` —— 归一化秩的纯数学性质
  （允许非生物学构造数据，受 AGENTS.md "纯软件策略测试" 条款许可）；
- `tests/scientific/test_ptm_crosstalk_pu_discipline.py` —— 断言输出字段与报告中
  **不出现 `negative`**，未标注 Cys 一律为 `unlabeled`。

**GREEN:** `src/plantpersulf/evaluation/ptm_crosstalk.py`

**必须同时输出的分层结果**（防止单一汇总数掩盖机制）:
- 按物种分层（拟南芥 / 水稻）；
- 按 Cys 包埋程度分层（用 §5.5.2 已有的真实 SASA 分箱）——
  **这是区分"真串扰"与"包埋效应机械后果"的关键**；
- 按蛋白无序度分层（若可得）。

---

### Task X4 — 三维空间串扰（H2）

**依赖 X0 完成**（需要 2,006 个结构）。

- 距离定义：Cys **SG 原子** 到磷酸化残基功能原子（Ser OG / Thr OG1 / Tyr OH）
  的最小欧氏距离；
- 复用 `features/structure.py` 的定宽 PDB ATOM 解析与
  `features/sasa.py` 的原子级基础设施；
- 同蛋白内归一化秩，与 X3 完全同构；
- **必须携带 pLDDT 掩码**：低 pLDDT 区域的坐标不可信，
  按 `features/structure.py` 既有的 `low_plddt` 语义处理，不得均值填补。

⚠️ **按 §2 声明的反向预期解读结果**：H2 的零结果不能算作假说被证伪。

---

### Task X5 — 判定门 `CROSSTALK_GATE`

沿用 Gate 2 的 fail-closed 风格，**条件在看到结果前冻结于
`configs/crosstalk_gate_v1.yaml`**：

| 条件 | 要求 |
|---|---|
| 1 覆盖率 | 两物种各自参与检验的蛋白数 ≥ 预设下限，否则 `insufficient_coverage` |
| 2 方向一致 | 拟南芥与水稻效应方向相同 |
| 3 统计显著 | 至少一个物种在 Bonferroni **与** BH 下都显著（复用 `bonferroni_and_bh_correction`）|
| 4 非包埋假象 | 效应在 SASA 分层内仍然存在 |
| 5 非单簇驱动 | 移除最大同源簇（`protein_clusters_v2.tsv`）后效应保持——**与 Gate 2 条件 5 同构** |

任一不满足 → `CROSSTALK_STOP`，结果**作为零结果如实发表**，不得重跑换设计。

**RED:** `tests/scientific/test_crosstalk_gate_is_fail_closed.py` ——
断言缺任一条件时判定为 STOP，且判定 JSON 可机器核验
（复用 `conclusion_gate.py` 模式）。

---

### Task X6 — 图与案例叠加

- 图：`scripts/make_crosstalk_figure.py` → `figures/fig2_ptm_crosstalk.{svg,pdf,png}`，
  遵循 `nature-figure` 技能规范，与 Figure 1 视觉体系一致；
- **张华团队案例叠加**：SlWRKY6 Cys396/Ser33、SlERF.D2 Cys35/Ser42 作为
  **注释性案例**画在图上，明确标注 `case_study_not_statistical_input`；
- 更新 `docs/phase_z_evidence_audit.md` 新增 §5.6，沿用既有的
  "变更记录 + 约束措辞" 格式。

---

## 6. 完成门（每个任务都要）

按 AGENTS.md：任务测试、全部快速测试、科学完整性测试、发布测试、Ruff、mypy，
并报告 RED/GREEN 命令与输出、数据溯源、完整性检查、局限性、commit SHA。

---

## 7. 预期结局与对应动作

| 结局 | 概率判断 | 动作 |
|---|---|---|
| H1 显著、H2 显著 | 低 | 逆先验的强结果，直接冲 NC 级；立刻找张华团队做湿实验验证 |
| H1 显著、H2 零 | **中** | 序列邻近串扰成立；主张收窄为"局部序列层面"，Plant Communications 稳 |
| H1 零、H2 零 | 中 | 如实发表零结果；论文回落到 Phase Z 的保守性发现，需张华团队数据补强 |
| 覆盖率不足 | 中高 | `insufficient_coverage`，本 Phase 降级——此时**番茄 persulfidome 成为唯一出路**，谈判优先级提到最高 |

**注意最后一行**：即使本 Phase 失败，它也产出了一个明确的、可以拿给张华老师看的
"为什么我们需要你们的数据"的技术论证——这比空口索要数据强得多。
