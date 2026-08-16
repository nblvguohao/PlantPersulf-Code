# 2026-08-15 升级迭代：共肽阴性 + 蛋白内排序 + 体制假说诊断

> 依据：《植物巯基化位点筛选_方法设计.md》与《PlantPersulf_增量建议_2026-08-14.md》，
> 在准确度战役（W1/W2/W3，`2026-08-15-model-accuracy-campaign.md`）定稿之后，
> 执行增量建议中**尚未落地**的三项：§1 共肽阴性（最高价值）、§2 蛋白内排序、
> §4 cys_density 体制假说。§3（番茄结构覆盖）已由战役 W1 完成（阴性）。
>
> 全部产物为 **diagnostic_only / pre_blind_expectation_setting_only** 级别：
> 不触碰冻结包（bundle `ab8a0353…`）、候选表、盲法队列与 SAP；无新预注册需求。
> 数字已同步 `2026-08-13-model-improvement-register.md` §4 更新日志。

---

## 1. 共肽阴性（增量建议 §1）：首个显式阴性证据 + 二值诊断

### 1.1 注册（`data/registry/copeptide_negatives_v1.tsv`，5 行 / 2 蛋白）

| 蛋白 | 肽段（注册蛋白组唯一命中） | 阳性 | **阴性** |
|---|---|---|---|
| BRG3 (A0A3Q7EW23) | `SSCMICLPCR` @204–213（kiad070 Fig.9B） | C206, C212 | **C209** |
| RNF144b (A0A3Q7GXU6) | `FYCPYKDCSAMLVNDSDEIVR` @120–140（kiad070 Table S3，`C3(S)`） | C122 | **C127** |

- 共肽阴性 = 同肽段、同谱图、同胰酶消化、同富集的未修饰 Cys —— 可检测性
  按构造匹配，是 AGENTS.md 允许注册的显式阴性类别；项目此前一条都没有。
- 纪律细节（`evidence/copeptide_negatives.py` 实现）：三态标签
  {positive, negative, undetermined}；无 site-determining-ion 覆盖时未修饰
  Cys 必须标 `undetermined`（未定位 ≠ 未修饰）；注册行逐残基核对注册
  番茄参考蛋白组，谱图级原始数据不可得已在 provenance 如实记录。

### 1.2 冻结 bundle 诊断（`results/diagnostics/copeptide_negatives_v1.json`）

| 蛋白 | 阳性排名 | 阴性排名 | 结论 |
|---|---|---|---|
| BRG3 (n=14) | C206 → **rank 8**, C212 → **rank 6** | C209 → **rank 5** | **未分离（阴性）** |
| RNF144b (n=24) | C122 → **rank 17** | C127 → **rank 5** | **未分离（阴性）** |

- 冻结模型无法把共肽阳性排到共肽阴性之上 —— 与信息论论证一致：三个 Cys 相距
  各 3 残基，±10 窗口输入几乎相同，序列特征（疏水性/蛋白级 Cys 密度/正电荷密度）
  在原理上无法区分。**实证证实"区分必须由结构承载"**（v2 方向证据 +1）。
- 两个结局都值得写：排得上 → 结构分支有真实内容；排不上 → 序列特征主导被证实。
  本次为后者。
- 方法论价值：这是唯一不受 Gate 2 条件 1（正样本实验室独立性死锁）限制的负样本轴
  —— 阴性来自论文层面即可注册，无需等待新研究。

## 2. 蛋白内排序（增量建议 §2）：训练目标 vs 报告口径的错位修复

### 2.1 12 个 mapped controls 全 Cys 冻结打分（`results/known_controls/within_protein_ranking_v1.json`）

| 蛋白 (n Cys) | 真位点 | 蛋白内排名 | 突变负担 (随机基线) | Top-1 / Hit@2 |
|---|---|---|---|---|
| **SlWRKY6** (7) | C396 | **1** | **1 (4.0)** | ✓ / ✓ |
| SlERF.D2 (2) | C35 | 2 | 2 (1.5) | ✗ / ✓ |
| BRG3 (14) | C206 | 8 | 8 (7.5) | ✗ / ✗ |
| RNF144b (24) | C122 | 17 | 17 (12.5) | ✗ / ✗ |
| CAT1 (10) | C234 | 8 | 8 (5.5) | ✗ / ✗ |
| AtG6PD6 (6) | C159 | 2 | 2 (3.5) | ✗ / ✓ |
| **PAD3** (8) | C440 | **1** | **1 (4.5)** | ✓ / ✓ |
| DES1 (3) | C44 | 3 | 3 (2.0) | ✗ / ✗ |
| RBOHD (10) | C825 | 2 | 2 (5.5) | ✗ / ✓ |
| SnRK2.6 (6) | C131 | 3 | 3 (3.5) | ✗ / ✗ |
| ABI4 (3) | C250 | 3 | 3 (2.0) | ✗ / ✗ |
| ATG4a (12) | C170 | 4 | 4 (6.5) | ✗ / ✗ |

汇总：Top-1 2/12、Hit@2 5/12、MRR 中位 ~0.5、总负担 54 vs 随机 58.5（k=1 口径）。

### 2.2 叙事修正

- 旗舰控制 **SlWRKY6 C396 在自身蛋白 7 个 Cys 中排第 1**，同时不在全局候选表
  Top-2000 —— 两个陈述可以同时成立，因为模型训练目标就是蛋白内 pairwise 排序，
  不是跨蛋白分数可比性。
- **结论**：该工具适用于"已知蛋白的位点定位"（四篇已发表工作流的真实任务），
  不适用于蛋白组盲扫；稿件叙事应以此为正面结论，取代近乎全阴的全局百分位口径。
- 局限：WRKY71 / ERF.D3 因 position_shift 状态未纳入（release 单位置惯例一致）；
  结构分支全程掩蔽（release 番茄打分语义）。

## 3. 体制假说（增量建议 §4）：模型事实 + 数据方向双重确认

### 3.1 关键模型事实（诊断前必须修正的前提）

冻结模型的序列特征为 3 维：疏水性、**蛋白级** cys_density（count(C)/蛋白长度）、
局部正电荷密度。**蛋白级密度在蛋白内是常数** —— 假说中的"同一特征需符号相反
权重"在现有模型上**不可实现**（局部窗口密度不是模型特征；加入需新发布预注册）。

### 3.2 局部密度分离方向（±10 窗口 C 计数，`results/diagnostics/regime_hypothesis_v1.json`）

| 蛋白 | 体制 | 真位点局部密度 | 干扰项均值 | 方向 |
|---|---|---|---|---|
| SlWRKY6 | IDR | C396: 0.048 | 0.095 | **真 < 干扰（符合）** |
| BRG3 | 金属簇 | C206: 0.238 | 0.158 | **真 ≥ 均值（符合）** |

两个典型体制的数据方向均与假说一致 → **体制路由（IDR / 结构域 / 金属簇专家，
MoE）列入 v2 架构候选证据 +1**；局部窗口密度列为新特征候选（新发布预注册）。

## 4. 文件清单

新增：
- `src/plantpersulf/evidence/copeptide_negatives.py`（三态分类 + registry 读写校验）
- `src/plantpersulf/evaluation/release_scoring.py`（冻结 bundle 共享打分路径）
- `src/plantpersulf/evaluation/within_protein_ranking.py`（蛋白内排序指标）
- `data/registry/copeptide_negatives_v1.tsv`（首个显式阴性注册）
- `scripts/evaluate_copeptide_negatives.py` / `scripts/diagnose_within_protein_ranking.py`
  / `scripts/diagnose_regime_hypothesis.py`
- `tests/unit/{test_copeptide_negatives,test_release_scoring,test_within_protein_ranking}.py`
- `tests/scientific/test_copeptide_negatives_registry.py`
- 产物：`results/diagnostics/{copeptide_negatives_v1,regime_hypothesis_v1}.json(+rows.tsv)`、
  `results/known_controls/within_protein_ranking_v1.json(+rows.tsv)`

## 5. 验收与纪律

- TDD：RED（模块缺失）→ GREEN（30 个新测试）→ 重构；全量 `tests/unit + tests/scientific`
  **706 passed**；ruff / mypy 全过；`audit-registry` / `audit-files` PASS。
- 未触碰冻结产物：bundle、候选表、matched controls、SAP/分析计划、盲法队列均只读。
- 无合成数据、无未登记输入；注册行全部对照注册蛋白组逐残基核验。

## 6. 追加：共肽阴性系统化扩产（PXD024061，2026-08-15 同日完成）

### 6.1 数据源与解析（`evidence/maxquant_sites_negatives.py`）

PXD024061（Aroca et al. 2021, Antioxidants 10:508，已注册 supplementary source）
的 MaxQuant sites 表——`Sulfide(C)Sites.txt`（82 行）+ `CianoBiotin(C)Sites.txt`
（6 行）——每行是一个候选定位；`Sulfide(C) Probabilities` 列携带完整肽段序列
与每个候选的定位概率（如 `VPSPTC(0.5)WC(0.5)SK`）。按 mod-peptide ID 分组后，
肽段的修饰状态可完全确定（或可证明不确定）：

- 候选概率 ≥ 0.75 → **positive**；< 0.75（弱候选/替代定位）→ **undetermined**
  （"未定位 ≠ 未修饰"纪律）；**非候选 Cys 仅当修饰总数确定**（Number 列或全
  强候选推导）时 → **negative**（确定的修饰总数封死"被修饰但没看见"的可能）。
- 坐标重新定位：MaxQuant 搜索坐标空间与注册蛋白组不一致（实测 Q9C9U5 pos313
  在蛋白组是 S）→ 肽段在注册拟南芥蛋白组 **v2 唯一命中**才注册。
- `classify_peptide_cys` 增加 weak-candidate 维度（默认空，无行为回归）。

### 6.2 注册结果（`copeptide_negatives_v1.tsv` 5 → 21 行，幂等追加）

| 肽段（蛋白） | 结构 | 结论 |
|---|---|---|
| 7 个新肽段 / 6 蛋白（拟南芥） | 全部同肽段正/负对比 | 11 条新增 negative + 5 positive |
| 代表：KPCFICGSLEHGAKQCSK (Q9FYD1) | C191 修饰，C194/C204 未修饰 | 3-Cys 三态注册 |
| 代表：SDEVKACIVTCGGLCPGINTVIR (Q9FKG3, CianoBiotin 探针) | C146 修饰，C150/C154 未修饰 | 2 条 negative |

跳过 73 肽段（74+6 中）：绝大多数为单 Cys 肽段（无共肽对比）、定位模糊
（如 0.5/0.5）、非唯一命中或 N_mod 不确定 —— 全部 fail-closed 处理。

### 6.3 扩展诊断（9 肽段冻结 bundle 蛋白内排名，`copeptide_negatives_v1.json` 重跑）

- **8/9 未分离**（positive 未排到 negative 之上）；BRG3/RNF144b 原结论复现。
- **唯一分离：Q944L8 C204(rank 2) > C209(rank 3)**——噪声水平单例：9 肽段中
  1 例分离在随机排序下并不罕见（2-Cys 肽段随机分离概率 0.5），且实际分离数
  低于随机期望（~4/9）→ 不构成"序列模型具备共肽分离能力"的证据，结论不变：
  **共肽区分必须由结构承载**。
- 统计口径修正：此前的"2 肽段全败"升级为"9 肽段 8 败 + 1 噪声例"。

### 6.4 验收

- 38 个针对性测试（新增 unit 20 + scientific 扩展 4 + 既有）全过；ruff/mypy 干净；
  `audit-registry` / `audit-files` PASS；重跑脚本幂等（0 新增）。
- 未触碰冻结产物；全部新行逐残基核验拟南芥参考蛋白组 v2（残基=C、肽段唯一、
  坐标一致、三态互斥）。

### 6.5 追加：PXD072089 水稻共肽阴性扩产（21 → 40 行，2026-08-15 同日完成）

第四物种（Oryza sativa）共肽阴性——PNAS 2025（`10.1073/pnas.2608150123`）
独立实验室，正是 Gate 2 条件 1 实验室独立性在阴性维度需要的补充。

- **数据源**：PNAS 补充 Dataset S4（论文确认的 -SSH 位点，每行一个修饰 Cys）
  与 Dataset S1（MaxQuant 肽段库存——每个被检持硫化肽段的全部 Cys 位点）。
  两表联合：检测肽段中论文确认修饰的 Cys 少于肽段总 Cys 时，其余 Cys 为
  共肽阴性候选。
- **坐标安全映射**（`evidence/pnas_site_negatives.py`）：S4/S1 的蛋白坐标在论文
  搜索空间——4/5 部分修饰肽段与注册水稻蛋白组一致，但 Q5ZCB1 差 1 个残基。
  通过肽段自身 Cys **秩**映射（修饰位点在库存中的排名 → 肽段内位置，坐标偏移
  不变），再在注册水稻蛋白组 v1 中唯一重定位取规范坐标。
- **注册结果**：5 个部分修饰肽段 → 8 个蛋白-肽段组 → **19 行**
  （10 阳性 + 9 阴性，`site_determining_ion_coverage=True` 论文报告型、同
  kiad070 番茄基础；raw spectra 不可得记为 limitation）。IIPTPNCALSSLGLPLRPGEPICTFYSR
  在同肽段内 4 个蛋白同源位点出现（A0A0P0Y2A9/A0A0P0Y253/U5KNJ1/Q2R4J4），
  Q5ZCB1 5-Cys 肽段（SLPPICHCADEVASCAAACKECDMVNSSSEPPR）3 修饰 / 2 阴性。
- **可复现代码**：这批行此前由一次性脚本生成、代码丢失且未提交；本次新增
  `evidence/pnas_site_negatives.py`（纯逻辑 + TypedDict，19 unit 测试）+

  `scripts/scan_copeptide_negatives_pxd072089.py`（pandas I/O + 重定位 +
  registry 写入，幂等）**精确复现**既有 19 行（`rows_added=0`），并新增
  scientific 测试（rice 专属 + Q5ZCB1 坐标偏移回归）。
- **验收**：unit 19 + scientific 10 全过；ruff/mypy 干净；逐残基核验水稻蛋白组 v1。

### 6.6 重跑冻结 bundle 诊断（40 行 / 17 蛋白 / 3 物种 + 置换检验）

`scripts/evaluate_copeptide_negatives.py` 加 rice 蛋白组后重跑（`copeptide_negatives_v1.json`
重写）。**17 个蛋白-肽段实例**（7 拟南芥 + 8 水稻 + 2 番茄），全部 Cys 冻结 bundle
打分（结构分支掩蔽）、蛋白内排名，判定"全部阳性是否排在全部阴性之上"。

- **观测**：**3/17 分离**（Q944L8——既知的噪声单例；A0A0P0Y2A9 与 A0A0P0Y2K3
  两个水稻 1v1 单例）。
- **置换零分布**（B=999，保每组 k_pos/k_neg、只在肽段自身 Cys 间洗标签）：
  **随机期望 7.6，零分布中位 8（p5=4 / p95=11），观测 3 位于左尾**。
  `p(分离能力超随机) = 0.997`；`p(≤随机) = 0.015`。
- **结论（比 9 组更强）**：冻结序列模型跨 **番茄 / 拟南芥 / 水稻 3 物种、17 个
  蛋白-肽段实例对共肽正/负零分离能力**——观测分离数甚至低于随机（近并列分数
  下相对顺序被噪声支配）。3 个"分离"单例全部在随机解释内。**"序列无法区分
  共肽位点"从 2 蛋白 2 肽段 → 17 实例 3 物种，统计口径升级为置换检验。**

## 7. P2 早期判定点（方法设计 §8）：结构特征分离共肽位点 —— 首个阳性判定

### 7.1 问题与工具

方法设计 P2 判定点："只用手工结构特征（RSA/pKa/金属配位/Sγ 静电势）+ 树模型，
看 SlBRG3 的 C206/C212 能否与 C209 分开。这一步决定整个项目值不值得做。"
此前序列侧已证伪（冻结模型 9 肽段 8 败）；**结构侧从未检验过**——本项补上。

- `evaluation/structure_features.py`：纯 numpy 特征管线（PDB 解析、Shrake-Rupley
  SASA、Sγ 几何、正电残基计数、库仑静电势、接触数），无 DSSP/FreeSASA/APBS 依赖
- 结构来源：注册表已有 AFDB v6（SlBRG3/SlWRKY6，PDB 与蛋白组逐残基一致校验）
- PyMYB10 无 AFDB 文件（API 有条目但模型文件全 404）→ 折叠另行决策，如实记录
- **监督树模型不可估**：labeled 仅 4 个（3 真 + 1 金标准阴性），21 Cys ——
  诚实记录为 gate 本身的样本量发现；改用逐特征方向 + 域内带符号合成

### 7.2 结果

| 检验 | 结果 |
|---|---|
| 逐特征方向（BRG3 真 {206,212} vs 阴 {209}） | **3/7 特征干净分离**：RSA 真更高、接触数真更低、库仑静电势真更低 |
| 单特征排序（WRKY6） | **RSA 与 Sγ 最近距离把 C396 排蛋白内第 2**（负担 2 vs 随机 4.0） |
| 全蛋白带符号合成（BRG3） | rank 9/10 —— **N 端无序区 Cys（C24/C184/C185）暴露度更高，抢走方向** |
| **RING 域内带符号合成（BRG3）** | **C206 rank 1/9、C209 rank 9/9（最后），负担 1 vs 随机 3.33** ✓ |

关键生物学读数：**C209 是 RING 簇内埋藏最深（rsa 0.007）、接触最多（26）、
静电势最高（0.352）的 Cys = 最像配位 Cys**；**C206 是簇内最暴露（0.238）、
静电势最低（0.159）的 Cys = 最像游离（HS⁻ 可攻击）Cys** —— 结构侧支持
方法设计 §1.3 的 RING 假设（9 个 Cys 必有非配位游离 Cys）。

### 7.3 判定结论

**P2 通过**：序列模型无法分离的共肽位点，结构特征在体制内可以分离
（C206 rank 1 vs C209 rank 9）——"区分必须由结构承载"从信息论论证升级为
**首个结构侧正证据**。同时产出两个直接输入：

1. **特征方向是体制内的**：同一符号在 RING 域内有效、跨体制失效 → 命题三
   "路由而非平均"的结构侧实证（MoE gate 必要性 +1）
2. **体制先验可操作化**：RING/金属簇域注释（或 pLDDT+域注释 gate）作为
   v2 结构分支的体制路由信号

局限：AFDB 无金属离子（配位为几何代理）；库仑势为介电简化的相对值（APBS
PB 为严格后续）；n=2 蛋白的案例证据，统计功效待 MIL 级标签（L0）扩充。

产物：`results/diagnostics/p2_structure_separation_v1.json`。

## 8. 下一步（新代设计输入，仍守 Gate 3 纪律）

1. **共肽阴性扩产**：用 PXD024061（Sulfide/CianoBiotin sites.txt 含
   "Number of Sulfide(C)"、localization prob）与 PXD072089 肽段表做系统性
   共肽阴性扫描（需肽段级链接，粗估数百条）→ 成为 v2 训练的显式阴性池。
2. **体制路由架构**（MoE：IDR / 结构域 / 金属簇专家，gate 由 pLDDT + 结构域
   注释驱动）：体制数据方向已验证（见 §9——**方向随体制翻转**，n=12 实证）；
   实现需新发布号 + 新预注册。
3. **结构承载区分**：AF3/Boltz-2 结构特征（金属配位、Sγ 静电势、RSA）加入后
   重测共肽分离（BRG3 C206/C212 vs C209）—— 这是"结构桥跨物种"命题的直接
   检验点，也是方法设计文档 P2 早期判定的对应物。**已部分推进**：§9 的 12
   对照体制内检验给出结构特征系统性排序证据；金属配位精确化（AF3+Zn /
   Metal3D）仍待做。
4. **稿件叙事**：蛋白内排序结果（SlWRKY6 rank 1/7）与共肽阴性诊断
   （序列模型原理性失效）写入稿件 R4 的已知控制章节。

## 9. 已知对照蛋白结构体制内检验（方法设计 §8 P2 的系统化扩展，n=12）

### 9.1 结构与工具

P2（§7）在 2 个蛋白上证明"结构特征能在体制内分离共肽位点"。本项把它系统化
到**全部 12 个 mapped controls**（5 番茄 + 7 拟南芥）。新增 7 个拟南芥对照的
AFDB v6 结构下载并注册（Q9FJI5/Q9LW27/F4K5T2/Q9FIJ0/Q940H6/A0MES8/Q8S929，
PDB/蛋白组逐残基一致、SHA256 审计通过）——**12/12 结构覆盖达成**。

- `evaluation/structure_regime.py`：pLDDT 体制分桶（folded ≥70 / linker 50-70 /
  disordered <50）、蛋白内 z、按特征聚合方向与负担、体制局部带符号合成
- `scripts/evaluate_structure_regime_separation.py` → `structure_regime_separation_v1.json`
- 10 个单元测试；ruff/mypy 干净；audit-registry / audit-files PASS

### 9.2 结果

**① Aggregate 方向翻转 naive 暴露假说（最重要）**：折叠体制真位点系统性
**更埋藏**——`contact_number_10a` 蛋白内 z 均值 +0.41、8/12 高于蛋白中位、
hit@2 7/12、总首中负担 **36 vs 随机 58.5（削减 38%）**；`rsa_relative` 反而
负担 69 > 随机 58.5。已验证功能位点（催化 Cys ATG4a C170、激酶 SnRK2.6
C131、RBOHD C825 等）多为埋藏/堆积残基，不是表面暴露 Cys。

**② 方向随体制翻转（命题三 n=12 实证）**：

| 体制 | n | 主导特征 | 证据 |
|---|---|---|---|
| folded | 9 | contact_number_10a | mean_z +0.72、hit@2 6/9 |
| disordered | 3 | rsa_relative / nearest_sg_distance | hit@2 全 3/3（暴露/孤立） |

单一全局符号复合无法同时服务两种体制（对 disordered 少数群体是近随机）→
**路由是必需的，不是可选优化**。

**③ 体制局部排序无符号即全特征有效**：仅在真位点自身 pLDDT 桶内排序，
所有 7 特征的负担都下降（最近 Sγ 59→45、RSA 69→53、接触 36→30、库仑
63→47）。

**④ 带符号复合上限**（全局符号，标注 ceiling）：蛋白内负担 37、体制局部
**32**（随机 58.5/47.5）；体制局部 top1 **5/12**、rank≤2 **9/12**。

**⑤ BRG3 RING 连续性 + 张力显式化**：P2 level-5 方法（RING 局部 z + 暴露
符号）复现 **C206 rank 1、C209 rank 9、负担 1 vs 3.33**；同一 RING 在全局
埋藏符号下 C209 升到 rank 7——**P2 的"游离暴露 Cys"是体制内解读**。且
cluster-like 真位点 RNF144b C122 完全埋藏（RSA 0.000、接触 z +1.60），
"金属簇内游离 Cys = 暴露"不推广为一般金属簇规则。

### 9.3 判定与 v2 设计输入

- 结构特征在 n=12 对照集上**系统性把真位点排前**（体制局部复合 rank≤2
  9/12），且方向**依赖体制**。
- **v2 结构分支必须体制路由**：folded 体制用埋藏/堆积（接触数），disordered
  体制用暴露/孤立（RSA、最近 Sγ 距离）；RING/金属簇内另需簇内暴露 gate。
- 局限（诚实记录）：contact_number 为主导是 7 特征集内的探索性选择，非预注册
  端点；监督树仍不可估（labeled 12）；AFDB 无金属（配位为几何代理）；MIL
  级标签（L0）是监督方向的先决条件。

### 9.4 LOO + 排列置换修正（2026-08-15 同日，重要）

§9.2-9.3 的带符号复合是**标注的 ceiling**（符号从 n=12 全量学）。`scripts/
evaluate_structure_regime_loo.py`（`structure_regime_loo_v1.json`，B=999
排列）把 ceiling 换成**留一估计 + 交换性零分布**，结果修正如下：

| 检验 | 结果 | 解读 |
|---|---|---|
| **contact_number_10a 单特征** | protein-wide 负担 36 vs 零分布中位 59，**p=0.018**；regime-local 30，**p=0.032** | **唯一显著特征**——埋藏/堆积是真信号 |
| 其余 6 特征单特征 | 全部 p>0.5（RSA protein p=0.84、regime p=0.85） | **全部不显著**，暴露/静电/簇几何无超出随机信号 |
| LOO 带符号复合（global/regime_grouped） | 39/38 vs 随机 58.5，p≈0.3 | **复合不能留一泛化**——7 特征等权带符号被 6 个噪声特征稀释 |
| disordered 暴露子主张（RSA hit@2 3/3） | 排列下不显著（n=3 过小） | **不成立为稳健证据** |

**修正后的结论（比 §9.3 更诚实）**：
1. **结构侧唯一稳健信号 = 埋藏/堆积（contact_number）**，p=0.018——12 个功能
   验证对照的真位点比同蛋白其他 Cys 显著更堆积。这是"功能位点结构化/埋藏"
   的稳健证据。
2. **§9.2 的"方向随体制翻转"与复合 ceiling 需要降级**：disordered 侧暴露信号
   （RSA）在排列零分布下不显著（n=3）；等权复合不泛化。
3. **v2 设计输入修正**：结构分支不直接用 7 特征等权复合；以 contact_number
   （埋藏）为首要特征，其余特征经 MIL 级标签（L0）学权重或做特征选择后再定；
   "暴露/孤立用于 disordered 体制"降为待更数据的假设，不是已证方向。

### 9.5 共肽结构分离检验（2026-08-15 同日完成）——零结果

把已验证的结构特征（含 contact_number）应用到**显式共肽阴性池**：9 个同肽段
组内（2 番茄金标准 + 7 拟南芥 PXD024061，全部有 AFDB v6 结构），对每肽段比较
MS 确认修饰（POS）vs 未修饰（NEG）Cys 的 7 个结构特征。`copeptide_structure.py`
（单元测试 12）+ `scripts/evaluate_copeptide_structure_separation.py` →
`copeptide_structure_separation_v1.json`（B=999 排列，保每组 k_pos/k_neg、只在
该肽段自身 Cys 间洗标签）。

| 特征 | posHi | posLo | mixed | p(posHi) | p(posLo) |
|---|---|---|---|---|---|
| plddt | 2 | 4 | 3 | 0.929 | 0.482 |
| rsa_relative | 4 | 4 | 1 | 0.601 | 0.556 |
| cys_count_8a | 0 | 1 | 8 | 1.0 | 0.464 |
| nearest_sg_distance | 3 | 1 | 5 | 0.338 | 0.947 |
| positive_residue_count_6a | 2 | 1 | 6 | 0.801 | 0.870 |
| coulomb_potential_sg | 3 | 5 | 1 | 0.839 | 0.370 |
| **contact_number_10a** | 3 | 5 | 1 | 0.759 | 0.234 |

**结论：7 特征无一在任一方向分离共肽正/负（全 p>0.05）。** 方向本身是混杂的：
两个金标准番茄肽段方向**相反**——BRG3（SSCMICLPCR）负 C209 最埋藏
（contact 26 vs 正 206/212 的 22/16，P2 RING 解读），RNF144B（FYCPY…）反而
正 C122 更埋藏（24 vs 负 127 的 14）。

**解读（与 §9.4 兼容，不矛盾）**：
1. **共肽阴性是结构与序列共同的"不可分辨区"**。序列模型 8/9 失败（§1.2），
   结构 7 特征 0/9 显著——**同一肽段内修改/未修饰 Cys 的分辨低于当前特征
   分辨率**（±10 序列窗口与 10Å 几何都够不到）。这不是结构假说的失败：它
   恰好验证共肽阴性的设计价值——检测匹配时**连结构也不能假称分离**，任何
   v2 分支对共肽内分离的声称都需要超出当前特征集的证据。
2. **不推翻跨蛋白埋藏信号**。§9.4 的 contact p=0.018 是"功能位点相对**同蛋白
   其他 Cys**更埋藏"（跨蛋白层面，LOO 证明）；共肽检验是"同一簇内正/负可
   分辨"（簇内层面）。两个层面方向可以不同且同时为真——P2 的"游离暴露 Cys"
   是簇内 case-level 读数，已验证**不推广**（§9.2 ⑤ + 本检验的方向混杂）。
3. **v2 输入**：结构分支的声称范围限定在**蛋白内排序**（功能位点相对同蛋白
   背景更埋藏），不扩展到**共肽内正/负分辨**；共肽池继续作为有效性约束
   （fail-closed 校验）与标定资源，而非当前特征分辨率下的可训练判别目标。
   （2026-08-16 §9.7 追加：Sγ/6 Å 换基重算后零结果不变，此结论与声称范围不受
   影响；注意 Sγ 与 Cα 接触数同属堆积密度一类，实为一次观测，化学占据类单独
   检验见 §9.9。）

### 9.6 水稻跨物种重测（2026-08-15 同日完成）——零结果稳健（14 组 / 3 物种）

用 PXD072089 水稻共肽池（§6.5）把 §9.5 的结构分离检验**跨物种复现**：
按**不同肽段**去重（IIPTPNC 肽段在 4 个旁系同源蛋白重复出现，只保留首个
实例，避免 4 倍加权），9 组 → 14 组（2 番茄 + 7 拟南芥 + 5 水稻）。水稻
5 组需要 AFDB 结构：3 个已注册（A0A0P0WV74/A0A0P0Y2A9/A0A0P0Y2K3），
P0C361/Q5ZCB1 用 `download_alphafold_structures_bulk.py` 补下（2/2 成功，
PDB/蛋白组长度+残基逐一校验通过）。脚本改为**缺结构跳过并计数**（诚实记录
dropped，本次 0 dropped）。`scripts/evaluate_copeptide_structure_separation.py`
→ 重写 `copeptide_structure_separation_v1.json`（B=999 排列，seed 20260815）。

| 特征 | posHi | posLo | mixed | p(posHi) | p(posLo) |
|---|---|---|---|---|---|
| plddt | 5 | 5 | 4 | 0.729 | 0.746 |
| rsa_relative | 6 | 6 | 2 | 0.635 | 0.552 |
| cys_count_8a | 2 | 2 | 10 | 0.822 | 0.701 |
| nearest_sg_distance | 3 | 5 | 6 | 0.884 | 0.433 |
| positive_residue_count_6a | 2 | 3 | 9 | 0.952 | 0.547 |
| coulomb_potential_sg | 7 | 5 | 2 | 0.390 | 0.833 |
| **contact_number_10a** | 6 | 6 | 2 | 0.502 | 0.538 |

**结果：水稻 5 组加入后，7 特征仍无一在任一方向分离共肽正/负（全 p>0.05），
方向依旧混杂。** 水稻 5 组 contact 方向：P0C361 pos_higher（19 vs 14）、
A0A0P0Y2A9 pos_higher（19 vs 15）、A0A0P0WV74 pos_higher（21 vs 19）、
A0A0P0Y2K3 pos_lower（10 vs 12）、Q5ZCB1 5-Cys 肽段 mixed
（正 145/152/156 = 24/22/18，负 143/159 = 21/18）——与番茄（BRG3 neg C209
最埋藏 vs RNF144B pos C122 更埋藏）的方向混杂一致。

**结论升级（§9.5 → §9.6）**：共肽"结构与序列共同不可分辨区"的零结果从
2 物种 9 肽段稳健扩展到 **3 物种 14 肽段**（番茄金标准 + 拟南芥 PXD024061 +
水稻 PXD072089，独立实验室/独立富集方案）。§9.5 的三条解读与 v2 声称范围
（结构分支限蛋白内排序、不扩共肽内分辨）**不变**——跨物种证据更足。
**Sγ/6 Å 换基重算后零结果依然稳健（§9.7）——但该零结果只覆盖几何/堆积密度
这一类（Cα 与 Sγ 接触数同属一类、实为一次观测）；化学占据类单独检验（§9.9）
亦打空。两类均无共肽内判别子。**

### 9.7 Sγ 层接触数重算（2026-08-16）——换原子基与半径，零结果稳健

把 §9.5/9.6 的接触数特征**从 Cα 计数换成 Sγ 层重原子计数**重算：新特征
`contact_number_sg_6a` = 该 Cys **Sγ 原子 6 Å 内其他残基重原子（非氢）数**
（排除同残基自身原子，测的是硫原子周围的局部堆积），替换原 Cα 10 Å 计数
（`contact_number_10a`）。同一 14 组 / 3 物种池、同一 B=999 组内置换零分布
（seed 20260815）。`scripts/evaluate_copeptide_structure_separation.py` →
`copeptide_structure_separation_sg6a_v1.json`（**新 track**，Cα 版 v1 JSON 原样保留）。

| 特征 | posHi | posLo | mixed | p(posHi) | p(posLo) |
|---|---|---|---|---|---|
| plddt | 5 | 5 | 4 | 0.729 | 0.746 |
| rsa_relative | 6 | 6 | 2 | 0.635 | 0.552 |
| cys_count_8a | 2 | 2 | 10 | 0.822 | 0.701 |
| nearest_sg_distance | 3 | 5 | 6 | 0.884 | 0.433 |
| positive_residue_count_6a | 2 | 3 | 9 | 0.952 | 0.547 |
| coulomb_potential_sg | 7 | 5 | 2 | 0.390 | 0.833 |
| **contact_number_sg_6a** | 6 | 5 | 3 | **0.505** | **0.738** |

**结果：换原子基（Cα → Sγ）与半径（10 Å → 6 Å）后，7 特征仍无一在任一方向
分离共肽正/负（全 p>0.05），聚合零分布位置不动（Cα 版 contact 6/6/2、p 0.502/
0.538 vs Sγ 版 6/5/3、p 0.505/0.738）。** 特征本身确实变了——14 组中 6 组的逐肽
方向翻转（BRG3 从 pos_lower → mixed：C206 的 Sγ 接触 **31** 为簇内最低、C212 46
最高，C209 42；P0C361 从 pos_higher → pos_lower 等）——但聚合不可分辨不变。

**判定修正（2026-08-16，重要口径更正）**：`contact_number_sg_6a` 与
`contact_number_10a` 是**同一个量的两个几何版本，不是两类假设**——两者都在回答
"这个 Cys 周围挤不挤"（堆积密度）；换 Cα/10 Å → Sγ/6 Å 改的是**怎么量**，不是
**量什么**。证据：p 值几乎重合（0.505/0.738 vs 0.502/0.538）、零结果位置几乎重合
（6/5/3 vs 6/6/2）。**因此它们同时为零是一次观测，不是两次——不能报成"在两个
结构分辨率上独立验证"（稿件措辞必须遵守，见下）。** §9.5/9.6 的零结果在
**几何/堆积密度这一类**成立且对测量方式稳健，但**未覆盖化学占据类**：配位 Cys
与游离 Cys 可以有相同堆积密度（金属位点所有配体挤在同一口袋），contact number
对"Sγ 是否已被占用"结构性失明——这是问错了问题，不是分辨率不够。化学占据类
的检验见 §9.9。

**v2 方向（不因本结果关闭，声称范围不变）**：Sγ 层下游特征（Sγ 特异 SASA、
Sγ 埋藏深度、最近 Sγ–Sγ 距离/二硫态、Sγ 3 Å 内 His/Cys/Asp/Glu 的 N/O/S 计数
= 金属配位代理、χ1 取向）仍可直接在既有 AFDB 结构上计算，是候选特征方向；但
**共肽内分离 gate 未被 Sγ/6 Å 接触换基点亮**。BRG3 仍是内置判据：35 aa 内 9
Cys 的 C3HC4 RING 只需 7 Cys + 1 His 作配体，必有 ≥2 非配体 Cys——金属配位代理
能否挑出 C206/C212（修饰）vs C209（未修饰）是下一个针对性的二值判别子检验
（Sγ 接触本身不能：C209 的 Sγ 接触 42 位于 C206 的 31 与 C212 的 46 之间）。

### 9.8 BRG3 金属配位代理（2026-08-16，最后一枪）——3 Å 特征退化，扫描亦非判别子

按 §9.7 末列的下一个检验执行：`metal_coordination_sg_3a`（该 Cys Sγ 3 Å 内
His/Cys/Asp/Glu 残基的 N/O/S 原子数，金属配位代理；AFDB 无金属离子故以配体原子
密度代替 Zn 配位）加入结构特征管线 + 单元测试 1；聚焦脚本 `scripts/evaluate_
brg3_metal_coordination.py` → `brg3_metal_coordination_v1.json`（同 14 组池、同
B=999 组内置换）。

**结果（3 Å，提交口径）：特征完全退化——14 组全部 mixed（p=1.000，零分布中位
0/0），BRG3 全部 14 个 Cys 的 metal 值全 0。** 原因：apo 态 AFDB 模型无 Zn 拉住
配位簇，最近配体原子（其他 Cys SG、His 咪唑 N 等）通常 >3 Å，低于 3 Å 噪声底。

**探索性半径扫描（未改提交特征，仅解释退化是否只是截断）**：

| BRG3 Cys | 4 Å | 5 Å | 6 Å |
|---|---|---|---|
| C206（POS，修饰） | 0 | 2 | 3 |
| C209（NEG，未修饰） | 0 | 4 | 5 |
| C212（POS，修饰） | **5** | 5 | **11** |

**判定：不是截断问题，而是模式与预言相反——金属配位代理不是二值判别子。**
预测"修饰=非配体=低、未修饰=配体=高"不成立：两个修饰 Cys 落在 Sγ 配体密度两个
**极端**（C206 最低：4 Å 0；C212 最高：4 Å 5），未修饰 C209 与 C206 在 4 Å 持平
（均 0）、5-6 Å 仅略高。任何半径下共肽三连体都是 mixed，C209 从未整体高于修饰集。
**局限（诚实记录）**：这是 apo 模型的局部堆积读数，AFDB 无金属离子，真实全酶 Zn
配位需 Metal3D/AF3-with-Zn 才可测；本结果为"apo 模型下该代理不构成判别子"，
不扩展为"金属配位不影响修饰状态"。

**v2 含义**：Sγ 层两个针对性候选（6 Å 接触换基、3 Å 金属配位代理）在共肽分离上
均为零结果；结构分支声称仍限蛋白内排序，共肽池继续作为 fail-closed 有效性约束。

### 9.9 化学占据类最后一枪（2026-08-16，收工判据）——几何聚类读出 Zn 位点，打空

堆积密度（contact/SASA/depth）对化学占据结构性失明，故最后一枪直接测占据类：
BRG3 RING 区（197–231）的 9 个 Cys Sγ 按四面体互距聚类读出两个 Zn 位点（AF 无
显式 Zn，但 RING 围着 Zn 折叠、AF 学到 holo 构象；RING 区 9 Cys 的 pLDDT
89.2–95.3，配体几何可信）。脚本 `scripts/evaluate_sg_chemical_occupancy.py` →
`sg_chemical_occupancy_v1.json`。方向不定死（"配位 Cys 被锁住没空被修饰" vs
"配位 Cys 是低 pKa 硫醇盐、亲核性最强，持硫化打掉锌指是文献机制"），任一方向
只要配位状态分开 C206/C212 vs C209 即判别子。**不用聚合置换检验**：占据特征在
大多数组内零方差（无金属位点/二硫键），聚合置换结构性测不出只活在少数占据组的
信号——改为方差筛选后逐组报告，BRG3 为事先声明的单案例。

**BRG3 RING Zn 聚类（Sγ 互距 cutoff 4.5 Å）**：

| 位点 | Cys |
|---|---|
| Zn 位点 A（4 Cys） | 197, 200, 218, 221 |
| Zn 位点 B（3 Cys） | **212**, 228, 231 |
| 游离（非配体） | **206, 209** |

7/9 Cys 落入两个 Zn 位点（4+3 = C3HC4 的 7 个 Cys 配体，正好）；两个非配体是
**C206（修饰）与 C209（未修饰）**。

**共肽三连体判定：打空。** 配位状态不能分开 C206/C212 vs C209：
C206（POS 修饰）游离 nearestSG=6.18 Å；C209（NEG 未修饰）游离 nearestSG=4.68 Å；
**C212（POS 修饰）在 Zn 位点 B，nearestSG=3.50 Å，是配体**。两个修饰 Cys 被劈开
（206 游离、212 配体），未修饰 C209 与 C206 同侧游离——与"修饰=非配体"（C212
反例）和"未修饰=配体"（C209 反例）都冲突。

**池内逐组（方差筛选后 7/14 组，逐组报告不合并）**：方向不一致——pos_closer 3
（Q944L8、A0A0P0Y2A9、A0A0P0WV74）、pos_farther 1（Q9FYD1）、mixed 3（BRG3、
RNF144B、Q5ZCB1）。注：A0A0P0Y2A9（C3H1 锌指）里**修饰 C80 是配体**（nearestSG
3.59 Å，占据）而**未修饰 C64 游离**（9.31 Å）——支持"配位 Cys 是低 pKa 硫醇盐
靶点"的第二种逻辑，但该方向在池内不系统（Q9FYD1 反之）。

**收工判定（停止规则"打空"分支命中）**：零结果同时覆盖**几何/堆积密度**（§9.7，
测两遍实为一次）与**化学占据**（本枪）两类。共肽内修饰/未修饰分辨在当前 apo
特征集下无判别子；结构分支声称限蛋白内排序。**特征搜索到此为止，全力转 Gate 3
盲测。**

### 9.10 机制重释指针（2026-08-16，文献驱动）

§9.5–9.9 的共肽零结果按 Corpas 综述（COPLBI-D-26-00068，审稿中）Fig.1C 的
两步持硫化约束获得**任务误设重释**：P(持硫化)=P(可氧化|结构)×P(持硫化|已氧化)，
共肽阴性的负类是 gate-1/gate-2 失败潜伏混合 → 单阶段标签下零分离是原理性上限，
不是结构假说失败。§9.9 的"收工"指**当前 apo 特征集下的判别子搜索停止**，不关闭
机制问题——两步框架使 §9.9 的 C212（配体且修饰）从反例翻转为例内支持（配位=
低 pKa 硫醇盐激活），并给 contact p=0.018 埋藏信号一个机制自洽解释。完整记录见
`2026-08-13-model-improvement-register.md` §2.3；稿件叙事草稿见
`manuscripts/plant_physiology/2026-08-16_twostep_mechanism_discussion_draft.md`；
Table 1 oxiPTM 位点核验审计见 `results/diagnostics/oxiptm_sites_v1.json`。
