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
