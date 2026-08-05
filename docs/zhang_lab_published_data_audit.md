# 张华团队(合肥工业大学)已发表数据审计

> 日期:2026-07-28 | 分支:claude/zhanghua-team-literature-review-9f48b1(基于 phase-a-persulfidation-site-evidence @ fcbd903)
> 检索工具:PubMed MCP、Europe PMC REST、PRIDE Archive v2/v3 REST、GEO esearch、UniProtKB REST、NCBI Gene esummary
> 前置说明:本次审计**先于**上游 Phase Z 工作展开,同步分支后发现上游已独立完成大量重叠工作(AtG6PD6/PAD3 同物种独立对照、PXD072089/PXD072300 水稻跨物种轨道、真实 MMseqs2 聚类)。本文档已核对并合并上游现状,只保留仍属**增量**的发现,避免与 `docs/phase_z_evidence_audit.md` 重复。

---

## 1. 裁定摘要

**张华团队已发表的过硫化(persulfidation)论文中,没有可直接用于构建位点级 benchmark 的公开数据集**——6 篇相关论文全部非开放获取、全部不在 PMC 开放全文库、且 PRIDE/GEO 投递记录均为零。位点信息只存在于正文/补充材料中(Data 等级 C),需订阅访问或直接向作者索取。

本次审计的增量价值有三点:

1. **BRG3 的 `unmappable` 状态已解决并已修复入代码**——问题不在于数据缺失,而在于此前按字面基因名"BRG3"检索,而论文的 Accession numbers 段实际给的是 NCBI GeneID `LOC101267168`。已解出 UniProt `A0A3Q7EW23`(243aa),Cys206/Cys212 均为真实 Cys,序列验证通过。见 §3。
2. **补齐了 registry 中缺失的出版物登记**——`known_controls.py` 里引用的多个 DOI(kiad070、kiae271、kiae560、tpj.70000、PNAS 水稻文等)此前从未出现在 `data/registry/publications.tsv` 中;已补齐可验证的部分。见 §5。
3. **发现两个尚未登记、且暂无法独立验证 accession 的候选对照**(SlWRKY71 C193/C198、PyMYB10 C194/C218)——按项目"无 accession 不得使用"的原则,本次**不**写入 `REGISTERED_CONTROLS`,只记录待办。见 §4。

同时更正一处此前(本次调研早期)对水稻数据的误判:PNAS 水稻 persulfidome(PXD072089)**不能**直接翻转 Gate 2——上游已经跑过跨物种迁移轨道,结论是"独立轴(实验室/化学/物种)均满足,但 Gate 2 条件 1 要求的是同一 benchmark 内的留一研究独立性,跨物种迁移证据仍在该框架之外"(`docs/phase_z_evidence_audit.md` v4/v5)。Gate 2 判定维持 STOP。

---

## 2. 张华团队论文逐篇裁定

| DOI | 年 | 期刊 | 对象 | 位点 | 是否已登记为已知对照 | 备注 |
|---|---|---|---|---|---|---|
| `10.1093/plphys/kiae271` | 2024 | Plant Physiol | SlWRKY6 | Cys396(+Ser33 磷酸化)| ✅ 已登记(`SLWRKY6_H2S_PHOSPHORYLATION`,mapped)| 摘要明确提到"SlWRKY6 在 SlLCD1-OE 叶片中呈现**差异 persulfidation**"——暗示存在未公开的番茄 persulfidome 全谱筛选,见 §6 优先请求 |
| `10.1093/plphys/kiad070` | 2023 | Plant Physiol | BRG3 / SlWRKY71 | BRG3 Cys206+Cys212;**SlWRKY71 Cys193+Cys198** | BRG3 ✅ 本次修复为 mapped;SlWRKY71 ❌ 未登记(accession 无法独立验证,见 §4) | 同一篇论文报告了两个持久化蛋白 |
| `10.1093/plphys/kiad100` | 2023 | Plant Physiol | **PyMYB10(红皮梨)** | **Cys194+Cys218,附 LC-MS/MS 肽段** `RAAC194PSIELEEELFTTFWFDDRL`、`RSC218ANFPEEGQSRS` | ❌ 未登记(accession 无法独立验证,见 §4) | 张华团队**唯一非番茄**、且**唯一给出位点级 LC-MS/MS 证据**的过硫化论文;跨物种对照价值最高 |
| `10.1111/tpj.70000` | 2025 | Plant J | SlERF.D2 | Cys35(+Ser42 磷酸化)| ✅ 已登记(`SLERFD2_H2S_ETHYLENE`,mapped)| — |
| `10.1093/plphys/kiae560` | 2025 | Plant Physiol | ERF.D3 | Cys115+Cys118 | ✅ 已登记(`ERFD3_H2S_CONTEXT`,`position_shift`)| A0A3Q7ESP9(316aa,以 Y 非 M 开头)N 端比论文模型多 13 残基;128−13=115、131−13=118 数值吻合,但仍需作者确认所用 CDS/Solyc 编号 |
| `10.1093/plphys/kiag512` | 2026 | Plant Physiol | SlWRKY6-SlGRF1-SlGIF2 | 无新位点(复用 kiae271 的 Cys396)| 不适用 | 已在上游文档中明确标注"不得重复计数为新独立验证单元" |
| `10.1111/nph.20431`(SlWRKY71/SlDCD1)、`10.1021/acs.jafc.4c09530`(SlIMPA3)、`10.1093/hr/uhad014`(SlDCD2)、`10.3390/ijms232012239`(SlMS1)、`10.1038/s41438-020-00439-1`(SlLCD1)等 | 2020–2025 | 多刊 | 番茄 H₂S 生成/信号通路机制 | **无 persulfidation 位点** | 不适用 | 仅提供机制背景,不含任何过硫化位点数据 |

---

## 3. BRG3 修复(已落地代码)

`src/plantpersulf/evaluation/known_controls.py` 中 BRG3 此前状态为 `unmappable`,理由是"按基因名 BRG3 检索只匹配到 143aa 的 K4BRG3,长度不足以容纳 Cys206/Cys212"。

**根因**:检索用了论文里的功能别名"BRG3",而不是论文 Accession numbers 段给出的正式标识符。

**解析链路**(2026-07-28,全部通过公开 REST API 验证,可复现):

| 步骤 | 结果 |
|---|---|
| 论文 Accession numbers 段 | NCBI GeneID `LOC101267168` |
| NCBI Gene esummary | "probable BOI-related E3 ubiquitin-protein ligase 3", *Solanum lycopersicum*, 1 号染色体 |
| UniProt REST(按 GeneID 反查)| `A0A3Q7EW23`(`A0A3Q7EW23_SOLLC`,RING-type domain-containing protein,GN=LOC101267168,**243aa**)|
| 序列验证(UniProt FASTA 直接下载)| Cys 全表:24,30,173,184,185,197,200,**206**,209,**212**,218,221,228,231 — 206/212 均为真 Cys,上下文 `...KSCNSRSSC[206]MICLPC[212]RH...` |

**已应用的代码修改**:`BRG3_H2S_UBIQUITINATION` 的 `uniprot_accession` 由空字符串改为 `A0A3Q7EW23`,`cys_position` 由 `0` 改为 `206`(代表性位点,沿用 ERF.D3 的"一个 lineage 一行"惯例——项目现有 `test_same_species_control.py::test_new_control_is_a_distinct_lineage` 强制每个 `mechanism_lineage_id` 恰好一行,因此 Cys212 未单独建行,而是记录在同一行的 provenance 说明中),`status` 由 `unmappable` 改为 `mapped`。

**测试验证**:`tests/scientific/test_same_species_control.py`、`test_known_control_recovery_integrity.py`、`test_split_leakage.py` 共 20 个测试全部通过(4 个因缺参考蛋白组 fixture 跳过,符合预期);全量 `pytest tests/ -k "not network"` 前后对比 267→268 passed、69→68 failed,失败集合完全一致(均为本 worktree 缺失 `data/raw/...` 原始文件导致,与本次改动无关)。

**同步更新**:`configs/splits/split_config_v1.yaml` 与 `split_config_v2.yaml` 的 `known_mechanism_holdout` 列表补齐为 6 项(原来只有 `SlWRKY6`/`SlERF.D2`/`BRG3` 三项,遗漏了已登记的 `ERF.D3`/`AtG6PD6`/`PAD3`)。该字段目前只是随 split 结果透传的文档性元数据,未被其他逻辑按名称消费,修改不影响任何断言。

---

## 4. 两个新对照线索——因无法独立验证 accession,暂不登记

调研发现张华团队论文中还有两处可作为额外训练外对照的持久化位点,但**均无法仅凭公开数据库独立验证其 UniProt/NCBI accession**(论文全文在付费墙后,PMC 无开放全文;NCBI/UniProt 按基因名搜索均为零命中):

| 候选 | 来源 | 报告位点 | 检索尝试 | 结论 |
|---|---|---|---|---|
| SlWRKY71 | kiad070(与 BRG3 同一篇)| Cys193 + Cys198(论文构建了 `WRKY71Cys193AlaCys198Ala` 双突变体,biotin-switch 信号消失)| UniProt `gene:WRKY71 AND organism_id:4081` 零命中;NCBI Gene `WRKY71[Gene Name] AND Solanum lycopersicum[Organism]` 零命中 | **不登记**,待作者提供 accession |
| PyMYB10 | kiad100(红皮梨)| Cys194 + Cys218,LC-MS/MS 肽段已知 | 论文正文提及 `Pbr016663.1`,但 NCBI/UniProt 均无法独立检索到该标识符对应的记录 | **不登记**,`Pbr016663.1` 需向作者/梨基因组数据库直接核实后才能使用 |

按 `AGENTS.md` 的原则("不注册 accession + 源文件 + SHA256 就不得使用"),本次**没有**把这两个对照写入 `REGISTERED_CONTROLS`——宁可留白,不猜号。这两项已列入 §6 的数据请求清单第一优先级。

---

## 5. 投递状态核查(可复现查询记录)

| 查询 | 工具/端点 | 日期 | 结果 |
|---|---|---|---|
| `persulfidation` 关键词全项目检索 | PRIDE Archive v2 `search/projects` | 2026-07-28 | 12 个项目命中,**无番茄、无梨**;唯一中国投递为南京林业大学谢彦杰组水稻项目(PXD072300/072089/072035)|
| `tomato AND (hydrogen sulfide OR LCD1 OR persulfidation)` | NCBI GEO esearch(`db=gds`)| 2026-07-28 | **0 命中** |
| 6 篇张华团队 persulfidation 论文逐篇 | Europe PMC REST(`isOpenAccess`/`inEPMC`)| 2026-07-28 | 全部 `isOpenAccess=N`、`inEPMC=N`,无 PMCID |

**已修复的 registry 缺口**(`data/registry/publications.tsv`,原文件只有 4 个拟南芥 PRIDE 研究 + 1 个番茄 iProX 研究 + 4 个 GEO,遗漏了上游本次同步带入的三个跨物种数据集的出版物登记):

- 新增 `PXD072089`、`PXD072300` → PMID 42479832,`10.1073/pnas.2608150123`(Lin, Zhou et al., PNAS 2026,水稻 persulfidome + 重组蛋白正交验证)
- 新增 `PXD063170` → PMID 40617858,`10.1038/s41467-025-61582-8`(Chen et al., Nat Commun 2025,稻瘟菌 CSE_OE/WT persulfidome)

**未处理、留待后续的同类缺口**:`PXD038309`(Nat Chem Biol 2023,人 MPST 研究,已被上游降级为"无发表位点表")的确切 DOI 本次未独立核实,未登记,避免猜测式录入。

**发现但未登记的第三个 PRIDE 兄弟数据集**:`PXD072035`——与已登记的 PXD072089(persulfidome)、PXD072300(重组蛋白正交对照)同属 Lin/Zhou et al. PNAS 2026 投递,但内容是**水稻叶片总蛋白质组**(DIA/Spectronaut 定量,150mM NaCl 盐胁迫 0/1/3/6/12h,非过硫化位点数据)。当前项目管线未对任何数据集做丰度归一化,故不构成阻塞;仅记录以备将来水稻轨道扩展时使用其作为 total-proteome 参考。

---

## 6. 收敛后的张华团队数据请求清单(供 `docs/zhang_collaboration_brief.md` 引用)

按本次审计结果重新排定优先级:

1. **【最高优先】kiae271 补充材料中 SlLCD1-OE vs WT 的差异 persulfidation 全表**——论文摘要已证实该筛选存在("SlWRKY6 undergoes differential protein persulfidation in SlLCD1-overexpressing leaves"),这是获得番茄 persulfidome 的最短路径,且团队自己已经做过、只是未公开数据表。
2. **SlWRKY71(kiad070)与 PyMYB10(kiad100)的准确 accession 确认**——两者均有位点级证据(甚至 PyMYB10 有 LC-MS/MS 肽段),但按公开数据库无法独立核实标识符,需要作者直接提供 UniProt/NCBI/梨基因组数据库的正式编号。
3. **ERF.D3 的 13 残基位置偏移确认**——需要作者提供其克隆/测序时实际使用的 CDS 或 Solyc 编号,以确定 A0A3Q7ESP9 是否为正确的参考序列版本。
4. **kiad070(BRG3/WRKY71)、kiad100(PyMYB10)的原始 LC-MS/MS 文件与搜索参数**——用于独立审计,复用项目已有的 PRIDE/UniProt 解析管线。

补充说明(用于协作信函措辞):独立化学/物种的验证集**已有公开替代**——PNAS 水稻 persulfidome(PXD072089,南京林业大学谢彦杰组,与本项目训练数据无实验室重叠)已被本项目用于跨物种迁移验证轨道,但按 Gate 2 的严格定义,该轨道**不能**替代"同一 benchmark 内独立研究留一验证"这一条件。因此张华团队的番茄位点级数据仍然是唯一能真正翻转 Gate 2 条件 1 的数据来源,而不再是"唯一出路"——这有助于把协作请求定位为"补齐关键环节"而非"救急"。

---

## 7. 与上游 `docs/phase_z_evidence_audit.md` 的关系

本文档不重复上游已经系统性完成的工作,仅记录增量:BRG3 修复、发现的两个未验证对照线索、registry 出版物登记缺口、以及对"水稻数据能否翻转 Gate 2"这一问题的确认性复核(结论:不能,维持 STOP)。完整的 Gate 2 判定记录、跨物种迁移轨道、结构覆盖率/聚类工作请参阅 `docs/phase_z_evidence_audit.md` 与 `docs/phase_f_gate2_decision.md`。
