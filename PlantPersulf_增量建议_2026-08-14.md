# PlantPersulf-Code 阅读后的增量建议

> 阅读范围：`AGENTS.md`、`docs/phase_f_gate2_decision.md`、`docs/cross_species_v2_findings_2026-08-14.md`、
> `manuscripts/plant_physiology/2026-08-13_pp_manuscript_draft_v1.md`、
> `src/plantpersulf/{features,models,evaluation,evidence,benchmark}/`、`results/known_controls/`、git log 25 条
> 日期：2026-08-14

---

## 0. 先说结论

我上一份文档里**大部分内容你们已经做了，而且做得比我写的更严格**。具体地：

| 我提的 | 你们的状态 |
|---|---|
| 阴性样本可检测性偏倚 | ✅ `observation_propensity.py`；AGENTS.md 明确 PU 框架 |
| MIL / 标签稀缺 | ✅ 选了更合适的 PU 而非 MIL |
| 蛋白内排序而非全局 AUC | ⚠️ **训练用了**（`ranking_loss.within_protein_pairs`），**但从未报告** ← 见 §2 |
| 跨物种迁移评估 | ✅ leave-study-out + 4 物种 + MMseqs2 真聚类 + 置换检验 |
| 与 Sul-BertGRU 对比 | ✅ 而且做了训练数据审计（2705 条里 1003 条中心残基不是 Cys —— 这个发现本身就值得写） |
| 结构层价值 | ✅ 且你们已经用 pLDDT 对照**推翻了自己的 v1 结论**，比我写的更进一步 |
| 共形预测 / 冻结 / 预注册 | ✅ Gate 0-6 + SAP 共签 + OSF 预注册，远超我的设想 |

所以下面只写**我确认代码库里没有、且能直接接上你们现有 gate 纪律**的四件事。按价值排序。

---

## 1. 【最高价值】共肽阴性：BRG3 的 Cys209

### 事实

`grep -rIn -iE "co-?peptide|same.?peptide|observed.?unmodified|explicit.?negative|negative.?evidence"` 在 `src/` 和 `docs/` 里返回 **零结果**。

而 AGENTS.md 写的是：

> Do not treat an unobserved persulfidation site as an experimental negative.
> The primary task is **positive-unlabeled unless explicit negative evidence is registered**.

**kiad070 提供了一个可注册的显式阴性，你们已经注册了它的两个阳性，但漏了这个阴性。**

kiad070 Fig.9B 的质谱肽段（我已对 A0A3Q7EW23 序列逐残基核验）：

```
蛋白位置    204 205 206 207 208 209 210 211 212 213
残基         S   S   C   M   I   C   L   P   C   R
修饰状态     -   -   ★   -   -   ✗   -   -   ★   -
                  C206         C209         C212
                  已注册      未注册       组内第二位点
```

**C209 的"未被修饰"是一次观测，不是一次缺席。** 同一条肽段、同一张谱图、同一次胰酶消化、同一次富集、同一个蛋白 —— 可检测性完全匹配。这正是你们 AGENTS.md 里说的 "explicit negative evidence"。

### 为什么这对你们**当前的具体困境**特别关键

`cross_species_v2_findings` §3 的结论是：结构层的 SASA 差异在 pLDDT≥70 之后全部消失，因为
`pLDDT ↔ SASA ↔ 质谱可检出性` 三者耦合。你们的处理是**过滤**（≥70），代价是番茄丢掉 40% 阳性（86→52）。

**共肽阴性从构造上消除这个混杂**：同一条肽段内的两个 Cys，不可能有可检测性差异。这是唯一一类
不需要 pLDDT 校正、不需要倾向性加权就干净的对照。

而且它**在原理上判定序列窗口模型**：C206/C209/C212 两两相距 3 个残基，±10/±15 窗口输入几乎相同。
任何纯序列模型在这三个位点上的输出必然接近，无法同时把 C206、C212 排在 C209 之上。
你们的冻结模型三个序列特征（hydrophobicity / cys_density / positive_charge_density）在这三个位点上
数值几乎不可区分 —— **这是一个可以今天就在冻结产物上跑、且结论是二值的诊断**。

### 可扩展成一条新的证据轴

这不是一个孤例。**任何持硫化质谱研究里，含 ≥2 个 Cys 而只有部分被定位的肽段，都产出共肽阴性。**
你们手上有 PXD006140 / PXD024061 / PXD072089 / PXD063170 的处理数据和 `evidence/mzidentml.py`，
以及 kiae271 DSs.xlsx 的 `Sequence window` + `Localization prob` 列。

粗估：Cys 在植物蛋白组约占 1.5–2%，胰酶肽平均 ~14 aa，含 ≥2 Cys 的肽段比例不低；
在数千个已定位位点的规模上，**共肽阴性可能有几百到上千条**。这是四物种通用的、
全新的、且完全符合你们数据纪律的负样本来源。领域内无人做过。

### 必须诚实处理的方法学细节（这部分本身就可发表）

共肽阴性有一个真实陷阱：**质谱"未定位"≠"未修饰"**。需要两层条件：

1. 该肽段上**至少一个** Cys 的 localization probability ≥ 阈值（你们已用 0.75）；
2. 阴性位置有**足够的 site-determining ions 覆盖** —— 即谱图确实有能力区分该位置。
   缺这一条时应标为 `undetermined` 而非 `negative`，进入第三类而不是阴性池。

建议在 `evidence/` 下新增 `copeptide_negatives.py`，输出三态标签 `{positive, negative, undetermined}`，
并像现有 controls 一样带完整 provenance（谱图 ID、localization prob、site-determining ion 覆盖判据）。

### 落地形态（不破坏冻结）

- 新增 registry：`data/registry/copeptide_negatives_v1.tsv`（带 SHA256）
- 新增诊断：`scripts/evaluate_copeptide_negatives.py`，claim_class = `diagnostic_only`
- **Gate 2 条件 1 的结构性死锁**（正样本全来自同一实验室）**在阴性维度上不适用** ——
  阴性来自不同论文、不同实验室、不同物种。这可能是唯一一条不需要等新研究就能打开的路。

---

## 2. 蛋白内排序：你们训练用了，但从未报告

### 事实

`src/plantpersulf/models/ranking_loss.py::within_protein_pairs` 存在，且被 `additive_pu_ranker.py`
用作训练损失。但 `results/known_controls/*.json` 和稿件 R4 报告的全是
**全候选表百分位**（`percentiles within the frozen candidate table`），核心陈述是：

> **0 of 99** published sites appear in the candidate table's top 2,000.

### 问题

**这两个是不同的任务，而 0/99 那个数字对应的任务没有人真的需要。**

kiad068/070/100/271 四篇的实际工作流全部是：先有目标蛋白（转录组/互作/表型挑出来的），
再问"这个蛋白的哪个 Cys"。没有任何一篇是从 179,736 个位点里盲扫。

- **全蛋白组 Top-K**：179,736 选 200，先验 0.1%
- **蛋白内排序**：SlWRKY6 是 7 选 1，PyMYB10 是 5 选 1，BRG3 是 14 选 2

在全蛋白组尺度上排不进前 2,000（前 1.1%），与在 SlWRKY6 的 7 个 Cys 里排第 1，
**在数学上完全可以同时成立** —— 因为跨蛋白的分数可比性本来就不是这个模型训练目标
（你们训的是蛋白内 pairwise）。

### 建议

新增一个诊断，claim class 与 R4 的 mock-blind 完全同级（`pre_blind_expectation_setting_only`，
诊断非门控，不改模型/特征/Top-K/SAP）：

对每个已注册 control，报告**其在自身蛋白全部 Cys 中的排名**：

| 蛋白 | Cys 总数 | 真位点 | 冻结模型蛋白内排名 | 随机期望 |
|---|---|---|---|---|
| SlWRKY6 | 7 | C396 | ? | 4.0 |
| PyMYB10 | 5 | C194/C218 | ? | 3.0（以功能位点 C218 计） |
| SlBRG3 | 14 | C206/C212 | ? | 5.0 |
| SlWRKY71 | ? | C35(/C40) | ? | — |
| RNF144b | ? | C122 | ? | — |
| CAT1 | ? | C234 | ? | — |
| + 7 个拟南芥 controls | | | | |

指标：Top-1 / Hit@2 / MRR / **突变负担**（= 命中真位点的期望 C→A 构建数，随机基线 `(n+1)/(k+1)`）。

**这是一个半天的活，跑在已冻结的 bundle 上，不动任何冻结产物。** 而且它可能显著改变稿件的
叙事重心 —— 如果蛋白内排名好而全局百分位平，那结论是"该工具适用于已知蛋白的位点定位，
不适用于蛋白组盲扫"，这是一个**有用的正面结论**，而不是现在这个近乎全阴的结论。

⚠️ 前提是诚实：如果蛋白内排名同样接近随机，那就照实写，并且它比"0/99 in top 2000"更能
说明问题出在哪。两个方向都值得跑。

---

## 3. 番茄结构覆盖 = 0，是 v2 最大的单一杠杆

稿件 R3：

> the frozen AlphaFold registry contains **zero tomato accessions**, so every tomato site is
> scored with the structure branch masked (coverage 0/179,736)

而 Gate 2 条件 4（唯一通过的结构条件）测的是拟南芥上的 structure gain（+0.0441）。
番茄 per-species AP = **0.0185**（拟南芥 0.493、水稻 0.647）。

**也就是说：在唯一要做盲测的物种上，跑的是一个被砍掉结构分支的残缺模型。**

而且 `cross_species_v2_findings` §0.4 显示你们**后来已经下了 76 个番茄 AlphaFold v6 结构**
用于保守性分析 —— 说明这不是数据不可得，只是冻结发生在下载之前。

v2 发布的第一顺位就是这个。AFDB 有完整的 *S. lycopersicum* 蛋白组；缺口部分可用 Boltz-2
本地补（MIT 协议、开源、A100 上很快）。当然按你们的纪律，这必须走**新 release 号 + 新预注册**。

---

## 4. 特征集的机理缺口，和一个可以立刻验的假说

### 现有特征全集

```
序列：hydrophobicity, protein_cys_density, local_positive_charge_density,
      local_negative_charge_density, local_cys_density, local_sequence_entropy
结构：plddt, low_plddt, contact_number_proxy
ESM：（冻结发布中已禁用）
```

**没有**：硫醇 pKa、金属/Zn 配位、二硫键状态、Sγ 处静电势、二级结构/无序度（pLDDT 之外）、
保守性（只做了事后分析，未进模型）、磺烯化前体（H2S 不与还原态 R-SH 反应，需先有 -SOH）。

而你们自己的 v2 结论已经指出：冻结模型的结构分支 = `{plddt, low_plddt, contact_number_proxy}`，
**恰好就是那个被证明在编码模型置信度而非结构生物学的三元组**。

### 一个可以今天验证的假说：`cys_density` 在两种体制间符号相反

把三个已知案例按局部 Cys 密度排一下（我已核验序列）：

| 位点 | 上下文 | local_cys_density | 体制 |
|---|---|---|---|
| SlWRKY6 **C396** ★ | `RAMLP·C·SSNMA`（C 端 IDR） | **低** | IDR 型 |
| SlWRKY6 C313/C320/C326/C335 | WRKY 锌指簇 | **高** | 干扰项 |
| PyMYB10 **C218** ★ | `LSARS·C·ANFPE`（C 端 IDR） | **低** | IDR 型 |
| PyMYB10 C47/C51 | `AGLNR·C·RKSCR`（R2R3 DBD） | **高** | 干扰项 |
| SlBRG3 **C206/C212** ★ | RING 域，35 aa 内 9 个 Cys | **很高** | 金属簇型 |

**同一个特征，在 WRKY6/MYB10 上要求负权重，在 BRG3 上要求正权重。**
一个 hidden width 16 的加性门控模型必然取折中，两边都做不好 —— 这与番茄 AP=0.0185
和 "0/99 in top 2000" 完全自洽，而且给出了一个**可证伪的机制解释**，不只是"数据不够"。

**诊断做法**（冻结产物上，零成本）：对这 3 个蛋白，把每个 Cys 的 6 维序列特征和最终分数导出来，
看真位点与干扰项在 `local_cys_density` 上的分离方向是否如上表所料，以及模型分数是否被它主导。

**如果成立**，v2 的架构方向就明确了：**按体制路由（IDR / 结构域 / 金属簇）而非全局加性融合**
—— gate 由 pLDDT + contact proxy + 结构域注释驱动。这条正好接你 CANOPY-Router 的 MoE 门控经验，
实现风险低。

---

## 5. 如果只做一件事

**做 §1 的 C209。**

理由：它同时命中你们三个最硬的约束 ——
(a) Gate 2 条件 1 的实验室独立性死锁，在阴性维度上不适用；
(b) `cross_species_v2` 揭示的 pLDDT/可检出性混杂，共肽阴性从构造上免疫；
(c) AGENTS.md 允许注册显式阴性，而当前一条都没有。

而且它有一个干净的二值结局：**冻结模型能否把 C206 和 C212 排在 C209 之上。**
排得上，结构分支有真实内容；排不上，序列特征主导的诊断被证实。两个结果都值得写进稿件。

---

## 附：已核验的序列事实

```
SlWRKY6  A0A3Q7F586  550 aa   7 Cys: 82, 298, 313, 320, 326, 335, 396
         C396 上下文 RAMLPCSSNMA —— 与 kiae271 报告的肽段 AMLPC396SSNM 一致
         C298-C335 五个位于 WRKY 结构域/锌指区

PyMYB10  A0A0U2QCQ1  244 aa   5 Cys: 26, 47, 51, 194, 218
         C26/C47/C51 位于 R2R3 DNA 结合域；C194/C218 位于 C 端 IDR
         kiad100：两者均被修饰，仅 C218 有功能表型

SlBRG3   A0A3Q7EW23  243 aa  14 Cys: 24, 30, 173, 184, 185,
                                     197, 200, 206, 209, 212, 218, 221, 228, 231
         肽段 SSCMICLPCR 起始于 204 → C206(★) / C209(✗) / C212(★)
         RING 区 35 aa 内 9 个 Cys；C3HC4 型 RING 仅需 7 Cys + 1 His 作配体
         → 很可能存在非配体游离 Cys（可被结构方法解出的明确预测）
```
