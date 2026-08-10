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
| `10.1093/plphys/kiad070` | 2023 | Plant Physiol | BRG3 / SlWRKY71 / RNF-144b | BRG3 Cys206+Cys212;SlWRKY71 Cys193+Cys198;RNF-144b(位点未公开,见 Supplemental Table S3)| BRG3 ✅ mapped;**SlWRKY71 ✅ 2026-08-10 解出并登记(`WRKY71_H2S_UBIQUITINATION`,position_shift,见 §8.1)**;RNF-144b ❌ 未登记(缺具体 Cys 位点编号,见 §8.2)| 同一篇论文 LC-MS/MS 报告了三个持久化蛋白/位点组 |
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
| ~~SlWRKY71~~ | kiad070(与 BRG3 同一篇)| Cys193 + Cys198(论文构建了 `WRKY71Cys193AlaCys198Ala` 双突变体,biotin-switch 信号消失)| UniProt `gene:WRKY71 AND organism_id:4081` 零命中;NCBI Gene `WRKY71[Gene Name] AND Solanum lycopersicum[Organism]` 零命中 | **已于 2026-08-10 解出并登记,见 §8** |
| PyMYB10 | kiad100(红皮梨)| Cys194 + Cys218,LC-MS/MS 肽段已知 | 论文正文提及 `Pbr016663.1`,但 NCBI/UniProt 均无法独立检索到该标识符对应的记录 | **不登记**,`Pbr016663.1` 需向作者/梨基因组数据库直接核实后才能使用(2026-08-10 复核未变) |

按 `AGENTS.md` 的原则("不注册 accession + 源文件 + SHA256 就不得使用"),本次**没有**把 PyMYB10 写入 `REGISTERED_CONTROLS`——宁可留白,不猜号。这一项已列入 §6 的数据请求清单第一优先级。SlWRKY71 的解决过程和一个新发现的候选(RNF-144b)见 §8。

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

> **2026-08-10 更新**:原第 1、2 项(kiae271 差异全表、kiad070 Table S3 的 RNF-144b 位点)均已从公开补充材料中直接解出,不再需要请求,见 §8.2、§8.3。清单重排如下。

1. **【最高优先】PyMYB10(kiad100)的准确 accession 确认**——已有位点级 LC-MS/MS 肽段证据(Cys194/Cys218),但正文给出的 `Pbr016663.1` 按公开数据库无法独立核实,需要作者直接提供 UniProt/NCBI/梨基因组数据库的正式编号。`kiad100_supplementary_data.pdf` 尚未精读,下次应优先检查其中是否已含答案(参照 RNF-144b 在 kiad070 Table S3 里被找到的先例)。
2. **ERF.D3 的 13 残基位置偏移确认**——需要作者提供其克隆/测序时实际使用的 CDS 或 Solyc 编号,以确定 A0A3Q7ESP9 是否为正确的参考序列版本。
3. **kiad070(BRG3/WRKY71/RNF-144b)、kiad100(PyMYB10)的原始 LC-MS/MS 搜索参数(检索引擎/数据库/FDR 阈值)**——三篇论文正文和目前读到的补充材料均未报告这些参数,独立审计仍需要原始质谱文件或至少完整的搜索软件配置。
4. **番茄成熟阶段更完整的组学数据**(见 docs/zhang_collaboration_brief.md 第二阶段)——119 位点表虽已到手,但只是叶片、单一 SlLCD1-OE vs WT 对比;番茄果实各成熟阶段的位点级数据、对应总蛋白组/磷酸化蛋白组仍是从"拟南芥迁移+单一对比公开表"升级到"番茄体系内系统建模"的关键缺口。

补充说明(用于协作信函措辞):独立化学/物种的验证集**已有公开替代**——PNAS 水稻 persulfidome(PXD072089,南京林业大学谢彦杰组,与本项目训练数据无实验室重叠)已被本项目用于跨物种迁移验证轨道,但按 Gate 2 的严格定义,该轨道**不能**替代"同一 benchmark 内独立研究留一验证"这一条件。因此张华团队的番茄位点级数据仍然是唯一能真正翻转 Gate 2 条件 1 的数据来源,而不再是"唯一出路"——这有助于把协作请求定位为"补齐关键环节"而非"救急"。

---

## 7. 与上游 `docs/phase_z_evidence_audit.md` 的关系

本文档不重复上游已经系统性完成的工作,仅记录增量:BRG3 修复、发现的两个未验证对照线索、registry 出版物登记缺口、以及对"水稻数据能否翻转 Gate 2"这一问题的确认性复核(结论:不能,维持 STOP)。完整的 Gate 2 判定记录、跨物种迁移轨道、结构覆盖率/聚类工作请参阅 `docs/phase_z_evidence_audit.md` 与 `docs/phase_f_gate2_decision.md`。

---

## 8. 2026-08-10 更新:SlWRKY71、RNF-144b 均解出并登记;kiae271 补充材料发现完整 119 位点表;Corpas 综述交叉核对

本次是对 `docs/zhanghua/` 下 kiad068/kiad070/kiad100/kiae271 四篇论文全文(而非摘要/正文片段)的逐篇精读,加上一篇在审的 Corpas 等 *Current Opinion in Plant Biology* 综述(`COPLBI-D-26-00068`,尚未发表)。

### 8.1 SlWRKY71:accession 就在论文自己的 Accession numbers 段里

kiad070 正文末尾"Accession numbers"段原话列出 `WRKY71 (LOC101264783)`——此前(§4)按基因名"WRKY71"检索零命中的原因,和 BRG3 当初一样:检索用了基因别名,没去查论文自己给的 GeneID 段。解析链路(全部公开 REST API 可复现):

| 步骤 | 结果 |
|---|---|
| 论文 Accession numbers 段 | NCBI GeneID `LOC101264783` |
| NCBI Gene esummary | "WRKY transcription factor 71-like", *Solanum lycopersicum*, 2 号染色体 |
| NCBI protein esearch(`LOC101264783[gene]`) | RefSeq `XP_004233015.1`(317aa)|
| 序列核对(RefSeq FASTA) | 论文抗体抗原肽 `CQVKKRVERSYQDP` 原样出现在残基 199-211;Cys 全表为 [193, 198] —— 与论文 `WRKY71Cys193AlaCys198Ala` 双突变体位点**零偏移完全吻合** |
| UniProt REST/ID mapping(按 GeneID、按 RefSeq accession 反查)| **零命中**——UniProt 尚未收录此 locus 的 GeneID/RefSeq 交叉引用 |
| 本地 `tomato_ref_proteome_v1.fasta` 全文比对 | 命中 `A0A3Q7FNU4`(159aa,"WRKY domain-containing protein"),与 RefSeq 残基 159-317 完全一致(100% 恒等、无空位)——是同一 locus 的 N 端截断版自动基因模型,不是不同基因 |
| A0A3Q7FNU4 自身编号中的 Cys 位点 | [35, 40](= RefSeq/论文编号 193/198,偏移 -158)|

**已登记**:`src/plantpersulf/evaluation/known_controls.py` 新增 `WRKY71_H2S_UBIQUITINATION`(`uniprot_accession="A0A3Q7FNU4"`, `cys_position=35`, `status="position_shift"`,provenance 中记录 Cys40 为同组第二位点)。`configs/splits/split_config_v1.yaml`/`split_config_v2.yaml` 的 `known_mechanism_holdout` 补充 `SlWRKY71`。测试验证:`tests/scientific/test_same_species_control.py`、`test_known_control_recovery_integrity.py`、`test_split_leakage.py` 共 24 个测试全部通过;全量 `pytest tests/ -k "not network"` 327 passed / 25 failed(失败集合与改动前一致,均为本 worktree 缺失原始数据文件导致,与本次改动无关)。

注:虽然 identity 证据比 ERF.D3 更充分(GeneID 直接匹配 + 抗原肽验证 + Cys 位点零偏移),但因为**实际用于打分的 UniProt 条目编号与论文编号不同**(-158 偏移,源于该 UniProt 自动条目缺失 N 端 158 个残基),按项目惯例仍标记 `status="position_shift"` 而非 `"mapped"`,provenance 中已写明这只是编号系统差异,不是身份存疑。

### 8.2 RNF-144b(LOC101265447):补充材料 Table S3 给出了具体位点,已解出并登记

kiad070 正文(L487)提到 LC-MS/MS 在两个番茄 E3 连接酶中检出硫巯基化信号:BRG3 和 **RNF-144b(LOC101265447)**,但正文本身没给位点编号。**kiad070 的补充材料 PDF(`kiad070_supplementary_data.pdf`)里的 Supplemental Table S3("Identification of the persulfidation peptide by LC-MS/MS")给出了这个位点**:

| 字段 | 值 |
|---|---|
| 肽段(Table S3 原文,大小写标记修饰残基) | `FYcPYKDCSAmLVNDSDEIVR` |
| Modifications 列 | `C3(S); M11(Oxidation)` |
| Protein Group Accession | `XP_004242195.1` |

解析:NCBI RefSeq `XP_004242195.1`("E3 ubiquitin-protein ligase RSL1-like", 233aa,GeneID `LOC101265447` 反查确认)全长序列里,肽段 `FYCPYKDCSAMLVNDSDEIVR` 精确唯一匹配到残基 120-141(字符串检索,非估算),Modifications 列的 `M11(Oxidation)` 落在该肽段唯一的 Met 上,交叉验证了 1-based 肽内编号规则;由此 `C3(S)` 对应**蛋白第 122 位 Cys**。本地 `tomato_ref_proteome_v1.fasta` 中的 `A0A3Q7GXU6`("RBR-type E3 ubiquitin transferase",318aa,同一 locus 的更长注释)前 233 个残基与 RefSeq 逐字符一致,Cys122 编号零偏移。

**已登记**:`known_controls.py` 新增 `RNF144B_H2S_UBIQUITINATION`(`uniprot_accession="A0A3Q7GXU6"`, `cys_position=122`, `status="mapped"`)。`control_type` 特意标为 `ms_detected_no_functional_validation`(而非其他对照用的 `strong_single_site_control`),如实反映证据强度差异:BRG3、SlWRKY71 都有 Cys→Ala 双突变体的功能验证,RNF-144b 只有 LC-MS/MS 检出 + 与 WRKY71 的荧光互补实验(信号弱于 BRG3,论文未继续深入),没有针对硫巯基化本身的功能验证。`split_config_v1.yaml`/`split_config_v2.yaml` 的 `known_mechanism_holdout` 已补充 `RNF144b`。测试验证:同 §8.1,24 个相关测试全部通过。

### 8.3 kiae271 补充材料:Dataset S1 其实是完整的 119 位点差异硫巯基化数据表,不只是 SlWRKY6

这是本次最重要的发现。`kiae271_supplementary_data.zip` 里的 `DSs.xlsx` 含 5 个数据集,其中 **Supplementary Dataset S1**(标题写的是"Persulfidation site analysis data of SlWRKY6 by LC-MS/MS",但实际内容是**完整的 119 行差异硫巯基化位点表**,SlLCD1-OE vs WT 番茄叶片,MaxQuant 风格的 site 表(真实 UniProt tr|...|..._SOLLC accession、Localization prob、PEP、Score、Sequence window、每个生物学重复的定量 Intensity)):

- 119 行、102 个不同蛋白;108/119(91%)定位概率 ≥0.75
- 49 个位点仅在 LCD1-OE 组检出(H2S 诱导获得)、49 个仅在 WT 组检出(H2S 抑制/组成性)、18 个两组都有、3 个两组都无(边缘/噪声)
- **SlWRKY6 Cys396 确实在表里**(`A0A3Q7F586`,位置 396,定位概率 1,仅 LCD1-OE 组检出,强度 221591 vs WT 组 0)——与正文叙述完全一致,证实这张表就是正文所说"119 peptides"筛选的原始数据,不是只截取了 SlWRKY6 一行
- 表里还有一个**未登记的 WRKY 家族蛋白**(`A0A3Q7GA72`,位置 230/240)和一个**未登记的 RBR 型 E3 连接酶**(`A0A3Q7HGP6`,位置 417/422,与 RNF-144b 是不同的蛋白/accession,但同为"H2S 靶向 E3 连接酶"这一功能类别的又一例)

**这意味着此前协作信里"请求 kiae271 差异硫巯基化全表"这一最高优先级数据请求,其实已经在公开补充材料里满足了**——该论文在 *Plant Physiology*(OUP/ASPB)发表,期刊政策为 Open Access,补充材料本身没有访问限制。这份数据尚未被写入 `known_controls.py` 或任何 benchmark/registry 文件,是否要以及如何使用(单独作为对照集扩展,还是作为番茄本地训练/测试数据的种子)是下一步的策略决定,不在本次"解析两个 accession"的范围内,留待 §9(或专门讨论)处理。原始表已导出至本地 `docs/zhanghua/_extract_kiae271/dataset_s1.tsv`(工作文件,未纳入版本库正式产物)。

### 8.4 kiad070 补充材料:phylogeny/CRISPR/primer 附图附表,无新增位点线索

`kiad070_supplementary_data.pdf` 的 Figure S1-S7、Table S1(转录组质控)、S2(表达相关性)、S4-S6(qPCR/克隆引物)均为背景性方法学材料,不含额外未登记的持久化位点(RNF-144b 的位点来自 Table S3,已见 §8.2)。`kiad100_supplementary_data.pdf`(13 页)本次未及精读,留待下次(可能含 PyMYB10 accession 的进一步线索,见 §6)。

### 8.3 kiae271:差异硫巯基化筛选的具体规模首次确认

此前(§6 第 1 项)只能引用摘要的定性描述"SlWRKY6 undergoes differential protein persulfidation"。精读原文后确认:该筛选是 SlLCD1-OE vs WT **番茄叶片**(非果实)4D label-free 定量蛋白组学,共找到 **119 个发生差异硫巯基化的肽段**,其中只有 1 个(SlWRKY6 上的)被深入表征;公开的 Supplementary Data Set 1 也只覆盖 SlWRKY6 单蛋白的位点谱图,119 肽段全表始终未公开。这个具体数字应当替换协作信草稿(`docs/zhang_outreach_draft.md`)中"请问这个筛选是否为全蛋白质组范围"的模糊问法,直接引用"119 个肽段"。

### 8.4 与 Corpas 综述(在审)交叉核对

`COPLBI-D-26-00068`(Corpas/Taboada/Molina-Escobar/Palma,CSIC 西班牙,*Current Opinion in Plant Biology* 在审)是目前 NO/H2S 衍生 PTM 领域最新的综述性汇总,给出跨物种(拟南芥/番茄/水稻/黄瓜/白菜/西瓜/小麦)约 35 条 persulfidation/S-nitrosation 位点表。**该综述参考文献列表中没有引用张华团队任何一篇番茄/梨硫巯基化论文**(kiad070/kiad100/kiae271 均缺席,只引了 SlERF.D2 与 SlDCD2 两篇)——这是当前领域文献综述里一个可以被填补的引用空白。该综述也证实磷酸化-硫巯基化竞争同一 Cys 已有发表先例(SnRK2.6 Cys131/Cys137 vs Ser175/Ser267,Chen et al. 2020/2021 *Mol Plant*),这对 PTM crosstalk grammar 阶段的假说是有利佐证,但也意味着"系统性规律"的立论必须做到结构/基因组层面的系统性证据,不能停留在单蛋白案例。
