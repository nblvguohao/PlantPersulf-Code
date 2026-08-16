# Kimi_Agent「PlantOxiPTM」方案差距与可行性审计（2026-08-16）

> 审计对象：`Kimi_Agent/` 下三份外部方案稿
> （`plan.md`、`nature_family_idea.md`、`idea1_executable_plan.md`），
> 代号 PlantOxiPTM——植物气体信号氧化还原开关的计算解码。
> 目的：与本仓库根目录已进行的工作比对差距，评估方案可行性与新颖性主张。
> 方法：逐条对照本仓库注册表（`data/registry/`）、已注册诊断结论、
> 竞争工具实测（Sul-BertGRU 阴性）与 PRIDE API 实时代理验证。

## 0. 结论摘要

1. **差距**：M1（数据库）与我们注册表 + 证据管线重叠约 60–70%，且方案对
   工作量乐观了一个量级；M2（预测器）高度重叠且其"首个预测器"新颖性主张
   已被实测反驳；M3（番茄 GRN + 结构/MD）是**真正空白**；M4（验证）已覆盖。
2. **可行性**：方案的核心科学主张（"过硫化位点可学"）在我们自己的诊断中
   **已实测为负**（oxiPTM discrimination p=0.976，方向相反）。这是标签定义
   问题，不是工程问题——在该标签下训练预测器为时过早。
3. **可吸收动作 3 项**：SNO 数据线（PXD055278，已动手完成注册 + 预检）、
   PXD005168 补注册、番茄 GRN 数据已注册但未使用的空白。
4. **建议**：不整体采纳该方案；吸收 3 项动作，GRN 工作流留待 Gate 3 后
   作为新预注册提出，其余以本仓库既有纪律推进。

## 1. 方案结构（四层金字塔）

```
Layer 4  故事层：把计算发现收束到"成熟/着色调控"的具体可检验预测
Layer 3  机制层：候选开关分子的结构模拟 + 调控网络定位
Layer 2  规律层：首个植物过硫化深度学习预测器 + 多 PTM 串扰图谱
Layer 1  资源层：PlantOxiPTM 数据库（PRIDE 重搜库 + 文献位点整合）
```

四条创新点主张：
- **N1**：首个植物专用、覆盖过硫化–亚硝基化串扰的 oxiPTM 整合数据库；
- **N2**：首个专门的过硫化深度学习预测器 + 同一位点多 PTM 竞争的全蛋白质组图谱；
- **N3**：气体信号 PTM 层首次系统接入果实成熟调控网络。

## 2. 逐层差距分解

| 层 | 方案内容 | 本仓库现状 | 差距判定 |
|---|---|---|---|
| M1 数据库 | 15+ 套 PRIDE 原始质谱重搜库（FragPipe），统一 FDR，整合文献位点，网页数据库 | 注册表已有 14 个数据集（PRIDE 9 / iProX 1 / GEO 4），证据预检管线 + 位点注册（`oxiptm_sites.py`），候选排序与结构特征 | **60–70% 重叠**。方案"重搜原始质谱"（自估 50–100 CPU 核周）在注册优先（register-first）路线下**不必要**；其"1–2 万位点"目标对工作量乐观一个量级 |
| M2 预测器 | 首个过硫化 deep predictor（ESM2/ProtT5 + 分类头），基准 pCysMod；过硫化 × SNO × 磷酸化 × 泛素化串扰图谱 | 竞争工具 Sul-BertGRU（*Bioinformatics* 2025 btaf078）已实测并登记为**阴性**；L11 oxiPTM 判别诊断实测 p=0.976 方向相反；串扰 grammar 为已注册主题 | **高度重叠 + 新颖性被反驳**。N2 的"首个预测器"不成立（Sul-BertGRU 已存在且我们的互跑为负）；"串扰图谱"与 L11 判别重叠 |
| M3 成熟 GRN + 结构/MD | fruitENCODE + TomExpress GRN 重建，候选枢纽结构 + MD 模拟 | 番茄多组学数据已注册（PXD051570、GSE163745/142713/142712/267238）；结构特征已覆盖候选枢纽（SlWRKY6、BRG3、PyMYB10…）；**GRN 重建本身与 MD 不在本仓库工作流** | **真正空白**。数据就绪但未建 GRN；MD 模拟超出本仓库范围（依赖外部计算线） |
| M4 验证 | 独立数据集回测 + PRM 靶向/点突变验证清单 | 湿实验验证矩阵 + Zhang lab Cys→Ala 功能验证工作包 | **已覆盖**。且方案的"回测用数据集"（PXD056815、PXD039999）中 PXD039999 已在我们注册表内 |

## 3. 新颖性主张逐条裁定（以实测证据为准）

| 主张 | 裁定 | 证据 |
|---|---|---|
| N1 首个植物 oxiPTM 数据库 | ⚠️ **部分成立**（资源缺口真实，但表述过度） | Plant PTM Viewer 无过硫化、dbPTM 无 persulfidation 独立类型已核实；但"重搜原始谱才可建库"的路径不成立，注册优先已覆盖大半；位点量级 1–2 万为乐观估计 |
| N2 首个过硫化 deep predictor | ❌ **不成立** | Sul-BertGRU（btaf078）已存在且经我们互跑**阴性**（2026-08-10 竞争工具审计）；本仓库 oxiPTM 判别诊断 p=0.976（方向相反）说明**当前标签下过硫化位点与未修饰 Cys 在特征空间不可分**——标签定义问题，非工程问题 |
| N3 PTM 层接入成熟网络 | ✅ **部分成立**（集成本身无先例） | 我们的候选体系已做"气体信号 PTM → 成熟/着色 TF"单点映射（SlWRKY6、PyMYB10 等）；"系统重建成熟 GRN 再投影"这一体量是新增量 |

## 4. 可行性总评

| 维度 | 判定 | 依据 |
|---|---|---|
| 数据可得性 | ✅ 清单大体真实 | 逐条代理验证：PXD005168（2015 过硫化蛋白，tag-switch）、PXD056815（番茄 SNO 定量，Plant Cell 2025）、PXD055278（拟南芥 SNO×磷酸化串扰）均真实存在；但方案**低估本仓库注册表**——其多数数据集我们已注册 |
| 方法可行性 | ❌ 核心主张已实测为负 | 预测器方向被 oxiPTM 判别 p=0.976 反驳；在该位点标签不可分的现状下训练预测器不具科学基础 |
| 算力/工期 | ⚠️ 乐观 | 重搜库 50–100 CPU 核周为真实成本但可规避；12 个月完成四层不现实，M1 重搜一项即是一个完整项目 |
| 数据线价值 | ✅ 有可吸收增量 | SNO 数据线可喂 L9 两步架构 stage-1；番茄 GRN 数据已注册未使用 |

## 5. 可吸收动作（3 项）

### 动作 1：SNO 数据线 —— PXD055278（已完成注册 + 预检）
- 数据集：*Nitric Oxide Regulates Stomatal Development and Stress Responses by
  S-Nitrosylation-Mediated Inhibiting Phosphorylation of MPK6*（PRIDE，2025-02-21，CC0）。
- 已入注册表：`data_sources.yaml` 记 `role: external_snitrosylation`；
  `download_selection.yaml` 选取 SEARCH xlsx + checksum.txt；metadata 已抓取
  （dataset_count 13→14）；小文件已下载并 SHA256 登记
  （xlsx `9bbe83da…`，10,764 B；checksum.txt `a810ec63…`）。
- 内容预检（见 §6）：SEARCH xlsx 为 Sequest/MaxQuant 风格肽段表，
  含 2 个 biotin-M（SNO 标签）Cys 位点；**MPK6 Cys201（AT2G43790.1 = Q39026）
  biotin-M 与注册表 `oxiptm_sites.py` 中 MPK6 C201=Q39026 SNO 对照一致**——
  独立数据集交叉确认已注册位点。
- 用途判定：数据可得性 + 位点交叉验证动作；**任何模型级使用须按新发布
  预注册纪律**（`evidence_preflight_v1.yaml` + 预注册），本轮不做。

### 动作 2：PXD005168 补注册（方案的主打过硫化源，2017，2015 蛋白 tag-switch）
- 当前未在我们注册表内；建议按标准流程补 metadata 注册，作为 Seville 体系
  之外的独立过硫化源交叉检查。

### 动作 3：番茄 GRN 数据已注册未使用（M3 的真正空白）
- PXD051570（成熟蛋白/磷酸化组）、GSE163745/142713/142712/267238 已注册；
  GRN 重建未启动。作为 Gate 3 后新预注册候选，而非随本方案仓促立项。

## 6. PXD055278 预检记录（内容层）

| 项 | 结果 |
|---|---|
| 文件 | `data/raw/PXD055278/ZJR_WDF_GSNO_1_20181120.xlsx`（10,764 B，CC0，SHA256 `9bbe83da…`） |
| 结构 | Sheet1：35 行 × 17 列；蛋白块结构（蛋白头行 + 肽段子头行 + 肽段行） |
| 蛋白 | AT2G43790.1 MAP kinase 6（15 唯一肽段）、ATCG00490.1 RuBisCO、AT1G14110.1 岩藻糖基转移酶 9、AT3G10550.1 Myotubularin-like |
| 修饰列 | 空、M8(Oxidation)、C7(Carbamidomethyl)、**C13(biotin-M)**、C3(biotin-M) |
| SNO 位点 | MPK6 内 `DLKPSNLLLNANCDLK`（C13 biotin-M）→ Cys201；AT3G10550.1 `ARcRLPVITWCQPGSGAVIAR`（C3 biotin-M） |
| RAW 文件 | 463 MB（`ZJR_WDF_GSNO_1_20181120.raw`）登记为 `remote_only`，未下载（非本轮需要） |
| 交叉验证 | MPK6 Cys201 与注册 SNO 位点一致 → 独立数据源确认（无新增位点待注册） |

## 7. 建议

1. **不整体立项**：方案四层中三层与既有工作重叠或已被反驳，唯一空白（GRN）
   应在 Gate 3 后以新预注册提出，而非随外部方案推进。
2. **吸收动作 1（已完成）与动作 2、3（待定）**：PXD005168 补注册与番茄 GRN
   工作流分别作为数据线与新预注册候选登记。
3. **守住标签纪律**：在 oxiPTM 判别为负的现状下，不因外部方案而重开
   预测器训练；任何预测器主张须先解决位点标签定义问题。

## 附：相关已注册文档
- 竞争工具审计：`docs/competitor_data_audit_2026-08-10.md`（Sul-BertGRU 阴性）
- 诊断登记：`docs/superpowers/plans/2026-08-13-model-improvement-register.md`
  （L11 oxiPTM 判别 p=0.976）
- 番茄多组学注册：`configs/data_sources.yaml`
