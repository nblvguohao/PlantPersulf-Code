# 植物硫巯基化数据资源与系统性证据审计(Phase Z)

> 版本: v11 | 日期: 2026-07-24 | 分支: phase-a-persulfidation-site-evidence
> 路线: Gate 2 = **GATE2_STOP** → 按 roadmap 降级路线产出的数据资源交付物
> 机器可核验判定: `results/external_validation/pu_ranker_v1/gate2_decision.json`
> 完整判定记录: `docs/phase_f_gate2_decision.md`
> v2 变更: 新增 §5.2 跨物种迁移轨道首跑结果(PXD063170);候选 B(PXD038309)
> 经核查降级(无发表位点表)。
> v3 变更: AtG6PD6 Cys159(PXD043969,西北农林)登记为同物种独立实验室
> 控制位点(§3.6);控制登记表迁入 src 并新增同物种语义;Gate 2 条件 3
> 独立验证单元 2 → 3(判定不变,STOP)。
> v4 变更: PXD072089(水稻 persulfidome,NM-biotin+DTT,谢彦杰/南京林业大学,
> PNAS 2026)解析建成并首跑跨物种迁移轨道(§5.3)。897 个坐标核验位点/646
> 蛋白从 6 份 PNAS 补充表中提取(SD01 肽段级 + SD04 位点级联合)。
> 跨物种迁移结果与 PXD063170 一致——弱不稳定(~1.3x 基准率,3/5 seed 低于基线,
> recall@50=0)。独立轴均满足(实验室/化学/物种),为未来 3 研究留一法验证
> 提供了数据基础;但 Gate 2 条件 1 的逐留一研究跨物种效应仍在框架之外。
> 措辞约束扩展至水稻轨道(§7)。
> v5 变更(2026-07-23): PRIDE 原始沉积核查——(a) PXD072089 的
> `SS-all-peptides.tsv`(提交者存放的完整 MaxQuant 肽段表)并入解析器,
> 坐标核验位点 897→**929**(§5.3 更新,轨道结果无实质变化);(b) UniProt
> 删除的 483 个 accession 全部通过 UniParc 找回序列,但核实后**未恢复任何
> 阳性位点**——发表/沉积表在生成时已排除这些位点,序列恢复≠位点恢复
> (§3.7);(c) PXD072300(10 个重组蛋白体外验证)解析为**同论文同实验室
> 正交验证控制位点**(23 个位点/8 个蛋白,§3.6),明确不计入 Gate 2 条件 1/3
> 的独立单元。
> v6 变更(2026-07-23): P1 同物种独立数据集检索第二轮(PRIDE 3 个关键词
> 变体 41 条项目 + Europe PMC 200 篇文献,2020–2026)——**未发现新的独立
> 蛋白质组级拟南芥/近缘种位点级数据集**,确认条件 1 的结构性阻塞维持原判。
> 但发现裴瑶/靳老师课题组(山西大学系,与 Seville Romero/Gotor/Aroca 网络
> 十余年零作者重叠)持续发表拟南芥单蛋白过硫化位点研究,登记其 **PAD3
> Cys440**(doi:10.1111/pce.70593)为第二个同物种独立实验室控制位点
> (§3.6)。控制恢复评估机器重跑:PAD3 9.8 百分位(未恢复,4 个已映射控制
> 中最低);Gate 2 条件 3 独立验证单元 3 → **4**(`gate2_decision.json`
> 机器重生成,判定不变,仍 STOP 3/5)。
> v7 变更(2026-07-23): 框架重构——从"预测"转向"保守性/趋同性"(§5.5)。
> 两层发现,均使用被 Gate 2 拒绝为预测证据的同一批三物种数据,作为
> 关联性证据(不改变 GATE2_STOP):(a) 家族级——PANTHER 同源子家族富集,
> 三对物种全部显著(p=7.98e-6 / 1.75e-2 / 1.31e-12),三方交集 8 个家族里
> 7 个是核心氧化还原/中心碳代谢酶;(b) 结构级——AlphaFold 结构覆盖率从
> 3.6%(旧 v4 模型链接已全面失效)提升到 98.0%(改用 v6,2,006/2,048 个
> 结构批量下载注册),持久化 Cys 相对同蛋白其他 Cys 系统性更包埋,**方向
> 三物种一致,但显著性仅 2/3**(水稻 p=0.009、稻瘟菌 p=0.001;拟南芥
> p=0.104,方向相同但样本量/方差下未达常规阈值——批量下载中途多次取用
> 中间进度数字曾短暂显示拟南芥 p=0.034,完整下载后已更新为 0.104,
> 早前"3/3 显著"表述已撤回,以本次为准)。两层发现互相印证,支持"趋同
> 发生在结构/通路层面而非局部序列层面"这一假说,但结构层证据强度
> 为中等(2/3 显著),如实标注,不夸大。
> v8 变更(2026-07-23): 自我审稿(模拟 Redox Biology/New Phytologist 审稿
> 视角)发现 §5.5.1 家族级富集三个检验未做多重检验校正即报告"全部显著"
> ——统计上不严谨。补做 Bonferroni + Benjamini-Hochberg 校正
> (`bonferroni_and_bh_correction`,不依赖 scipy)后:拟南芥×水稻、
> 水稻×稻瘟菌两对在两种方法下都稳健显著;**拟南芥×稻瘟菌这一对
> 方法依赖**——Bonferroni 下不显著(0.0525>0.05),BH 下仍显著
> (q=0.0175)。§5.5.1/§5.5.3 措辞已改为"2/3 物种对方法无关地显著",
> 撤回"三对全部显著"的表述。判定不变(GATE2_STOP,§5.5 全部数字仍为
> 关联性证据,不改变 Gate 2)。
> v9 变更(2026-07-23): 继续按自我审稿意见修正——§5.5.2 结构级发现此前
> 只依赖 `contact_number_proxy`(明确标注的非严格可及性代理),审稿人
> 指出应补做真实 SASA 验证。新增 `src/plantpersulf/features/sasa.py`——
> 纯 Python Shrake-Rupley 实现(Bondi 范德华半径 + 确定性 Fibonacci
> 球面采样,不依赖 numpy/scipy/Biopython),对已下载的 2,006 个结构算出
> 全残基与 SG 原子专属 SASA(全数据集耗时 <1 分钟)。**结果:真实 SASA
> 与代理指标给出完全一致的定性格局**(三物种方向一致,水稻/稻瘟菌显著、
> 拟南芥不显著)——结构级发现现在有两条独立方法交叉验证,不再是单一代理
> 指标的孤证。副产物:发现并修复了一个可复现性缺陷——`build_structural_
> context_rows`(以及新增的 SASA 版本)此前按未排序的 Python `set` 迭代
> 蛋白列表,而字符串哈希按进程随机化,导致同一 `--seed` 在不同进程运行间
> 给出略有差异的 p 值(拟南芥 p 值曾在 0.034–0.104 之间波动);已改为
> 按 accession 排序处理,两次独立进程运行现在逐位精确一致。判定不变
> (GATE2_STOP)。
> v10 变更(2026-07-23): 三项并行交付——(a) 发布图表
> `figures/fig1_cross_species_conservation.{svg,pdf,png}`(`scripts/
> make_conservation_figures.py`,遵循 nature-figure 技能规范),三面板
> 汇总 §5.5.1 家族级富集(含三重校正对比)与 §5.5.2 结构级双路径趋同性;
> (b) 新增 §5.5.4 归因偏差检验,回应自我审稿 Reviewer 2 关切——检验
> UniProt accession 删除是否按氧化还原功能选择性丢弃蛋白(会混淆
> §5.5.1/§5.5.2 的氧化还原富集发现)。用 SD01/SD04 中即使已删除
> accession 仍保留的"Protein description"字段做 kept/dropped 两组
> 氧化还原关键词构成的双侧 Fisher 精确检验:**p=0.360,未发现选择性
> 偏差的证据**,§5.5 的氧化还原富集发现不能用差异性归因混杂来解释;
> (c) 新增 §5.5.5 真实文献新颖性核查(非模拟检索)——逐一核查三方交集
> 8 个家族是否已被单独报道为过硫化靶点:甘氨酸脱氢酶复合体、乙醛酸
> 转氨酶网络**已是 Seville 团队自己 2023 年综述中明确列出的已知靶点**,
> FBPase 是已知氧化还原敏感位点(经由其他 PTM);另外 5 个家族(NADH-
> 泛醌氧化还原酶 75kDa、山梨醇脱氢酶、L-苏糖酸脱氢酶、PPIase、RCC1)
> 检索未发现已发表的单蛋白过硫化报道。**关键澄清**:本项目自己拟南芥
> 阳性集合的源论文(Aroca et al. 2017)本身已报告"中心碳代谢富集"这一
> 定性格局,因此该定性观察在单物种层面并不新颖;如实可主张的新颖性
> 收窄为**跨物种统计趋同性证据**(§5.5.1/§5.5.2 的富集检验与结构方向
> 一致性),而非"发现全新靶点家族"。判定不变(GATE2_STOP)。
> v11 变更(2026-07-24): **§3.3 的同源聚类缺口已解决**——本机无
> Linux/WSL,经两跳 SSH(本机→实验室 Windows 主机→实验室 A100 服务器,
> 流程记录于 `docs/ops/remote_a100_via_lab_jump_host.md`,工具通用化为
> `scripts/remote/lab_a100_hop.py`,可跨项目复用)在 A100 上装
> MMseqs2、对拟南芥全参考蛋白组(54,646 条序列)跑
> `easy-cluster --min-seq-id 0.3 -c 0.5 --cov-mode 0`(6.9 秒),产出
> **13,967 个真实同源簇**(`data/processed/clusters/protein_clusters_v2.tsv`,
> SHA256 `e13d16ae…967197`);`protein_clusters_v1.tsv`(singleton 占位)
> 冻结保留,不覆盖。用真实簇重新验证 Gate 2 条件 2(效应 CI)与条件 5
> (单簇支配)——两者结论均不变(条件 2 仍在同一折跨零、条件 5 仍
> 不受单簇驱动),但条件 5 现在建立在一个真正有意义的检验上(移除的
> "最大簇"从 v1 的单蛋白 7 行变成 v2 的真实 30% 同源家族 197 行,
> 移除后 AP 不降反升)。补充轨道 `pu_ranker_cluster_v2`(用真实簇重跑
> Split A)数字比 v1 更不稳定而非更干净——机制已查明并记录(§3.3):
> 真实聚类正确地把大型旁系同源家族作为整体分入单一切分,这是
> singleton 占位文件在结构上无法做到的,方差增大是"检验变严格"的
> 代价而非缺陷。判定不变(GATE2_STOP)。**副产品(已记录未处理)**:
> 复核过程中发现 AlphaFold 结构注册表已从 7 个扩到 2,006 个(见
> §5.5.2 的批量下载工作,当时未提交),若不隔离该变量会与聚类变量
> 混杂;已定位并隔离(临时还原注册表到冻结的 7-结构状态完成本次
> 验证,随后完整恢复 2,006-结构版本,SHA256 校验一致)。一次非正式
> 的(未隔离前的)复核显示 `seq_structure` 臂在全结构覆盖下 leave-
> study-out AP 可能大幅跳升且逐 seed 极不稳定——**该数字未经验证,
> 不得引用**,是否代表真实能力提升、过拟合还是特征管线混淆尚不清楚,
> 需要独立的、走完整 TDD 流程的 `pu_ranker_v2` 冻结实验来回答,本次
> 不处理,详见 `docs/phase_f_gate2_decision.md` "P4" follow-up。

**约束措辞(对所有引用本文的外部文本具有约束力)**:

> 当前公开数据不足以证明跨研究预测能力,模型仅用于候选组织与假设生成。
> Current public data are insufficient to demonstrate cross-study predictive
> ability; the model is used only for candidate organisation and hypothesis
> generation.

---

## 1. 摘要

本审计回答一个问题:**当前公开数据能把植物硫巯基化(persulfidation)位点优先级排序支撑到哪一步?** 答案是:

1. **数据资源已建成且可复现**:`benchmark_v1`(390 个经逐残基坐标验证的位点级阳性、395,488 未标记、0 阴性、纯 PU 语义),全部输入带 SHA256 溯源。
2. **方法信号真实存在但有限**:序列+结构分支的 PU ranker 在留一研究验证下达到基准率 2.1 倍,结构分支增益经配对检验稳健(CI [+0.024, +0.064]);但相对 PU 基线的优势在其中一折上经蛋白级重采样后 CI 跨零(实测脆弱性,非假设)。
3. **天花板在数据不在方法**:留一研究与数据集内切分结果基本一致;限制因素是公开数据的来源集中度(两研究、同实验室、同化学、同物种)与结构覆盖率(3.6%),而非评估口径或模型容量。
4. **结论闸门 Gate 2 = STOP(3/5 条件通过)**:唯一结构性阻塞是研究独立性;唯一的实测弱项是基线优势的蛋白级脆弱性。翻转 Gate 2 需要新的独立数据,而非更多模型工作。

This audit records what the public Arabidopsis persulfidation data can and
cannot support: a reproducible PU benchmark and a modest, honestly-bounded
ranking signal, blocked from any cross-study claim by study non-independence
and by a measured cluster-composition fragility — not by modelling choices.

---

## 2. 数据资源:`benchmark_v1` 数据卡

| 项 | 值 |
|---|---|
| 行数 | 395,878(每个蛋白每个 Cys 一行) |
| 阳性 | 390(全部 `evidence_level = site_ms`,位点级质谱证据) |
| 未标记 | 395,488(PU 语义:未检出 ≠ 阴性,从不作硬阴性训练) |
| 阴性 | 0(设计如此) |
| 覆盖蛋白 | 50,892(拟南芥参考蛋白组) |
| 含阳性蛋白 | 350 |
| 单蛋白最多阳性位点 | 7(O03042) |
| sites.tsv SHA256 | `38617833…d5904c`(见 `data/processed/benchmark_v1/manifest.json`) |
| 蛋白组 SHA256 | `51559016…d1033bf` |

**阳性来源构成(决定性事实)**:

| 研究 | 阳性位点 | 阳性蛋白 | 实验室 | 化学 | 物种 |
|---|---|---|---|---|---|
| PXD006140(Aroca 2017) | 317 | 284 | Romero/Gotor (Seville) | tag-switch | 拟南芥 |
| PXD024061(Jurado-Flores 2021) | 73 | 70 | 同上 | 同上 | 同上 |

**配套资源**:番茄参考蛋白组 v1(控制位点提取用)、7 个 AlphaFold/SWISS-MODEL
结构文件(官方 `pu_ranker_v1` 发布冻结状态;结构注册表此后另行扩至
2,006 个,见 §3.2/§5.5.2,尚未用于任何官方重新验证)、
`protein_clusters_v1.tsv`(singleton 占位,冻结保留)与
`protein_clusters_v2.tsv`(真实 MMseqs2 聚类,13,967 簇,见 §3.3)。

---

## 3. 系统性证据审计

### 3.1 证据集中度

- 全部 390 个阳性来自**同一实验室、同一 tag-switch 化学、同一物种**;无独立实验室、无第二化学方法(如 BTD/CyMPL)、无第二位点级物种。
- 证据等级单一:全部为 `site_ms`(质谱位点);无逐位点生化验证的独立阳性集可充当同物种金标准。
- 张华团队已发表的番茄机制位点(4 个登记控制)是唯一跨物种、跨实验室证据,其中 2 个可映射、1 个不可映射、1 个位置待确认(§3.6)。

### 3.2 结构覆盖率

仅 7 个结构文件覆盖基准:14/390 阳性(3.6%,6 个蛋白)+ 113 行未标记
——这是官方 `pu_ranker_v1` 发布时冻结的结构注册表状态,条件 4 的
结构增益数字(§3.5)即在此状态下计算,保持有效。**登记(2026-07-24)**:
`data/registry/alphafold_structures.tsv` 此后在 §5.5.2 的批量下载工作中
扩到了 2,006 个结构(98.0% 覆盖三物种持久化蛋白集合),但**尚未提交、
尚未用于任何正式重新验证的官方实验**——一次非正式复核显示全结构覆盖下
`seq_structure` 臂的 leave-study-out AP 可能大幅且不稳定地跳升,该数字
明确标注为不可引用,详见 `docs/phase_f_gate2_decision.md` "P4" follow-up。
结构分支的实测增益(§3.5 条件 4,基于冻结的 7-结构状态)因此**不能归因于
结构覆盖的阳性本身**;结构覆盖子集(16 行)上的配对差值 CI 跨零
(+0.009 [−0.087, +0.056]),该事实如实记录,不被运行级增益掩盖。

### 3.3 同源性与聚类控制(2026-07-24 已解决)

**已解决**:`protein_clusters_v1.tsv` 曾是**每蛋白一簇的 singleton 占位
文件**(`_build_singleton_cluster_file` 生成);真 MMseqs2 聚类需要
Linux 主机,本机无 mmseqs/WSL/conda。本机无法安装 WSL(需管理员权限
启用系统功能,当前会话权限不足),改走**两跳 SSH**:本机 → 实验室
Windows 主机(`100.66.9.81`,已有对 A100 的免密钥登录)→ 实验室
A100 服务器(`100.112.165.109`,双卡 A100-80GB,128 核,251GB 内存)。
完整流程、坑点、可跨项目复用的通用工具见
`docs/ops/remote_a100_via_lab_jump_host.md` 与
`scripts/remote/lab_a100_hop.py`。

在 A100 用户空间(无 sudo)下载 MMseqs2 官方静态二进制
(`mmseqs-linux-avx2.tar.gz`),对拟南芥全参考蛋白组(54,646 条序列,
28MB FASTA,上传前后 SHA256 校验一致)运行:

```bash
mmseqs easy-cluster arabidopsis_ref_proteome_v1.fasta clusterRes/arabidopsis_v2 tmp \
    --min-seq-id 0.3 -c 0.5 --cov-mode 0 --threads 32
```

参数对应 TDD Codex §4.3 默认阈值(30% 一致性、≥50% 覆盖度)。**6.9 秒**
跑完,产出 **13,967 个真实同源簇**,下载回本地并转换为项目格式
(`protein_accession\tcluster_id`,取代表序列 accession 作为 cluster_id),
写入 `data/processed/clusters/protein_clusters_v2.tsv`
(SHA256 `e13d16aee81e0b98c9687bc98f8565122fbc027249d9e020a336a4f509671978`)。
`protein_clusters_v1.tsv` 冻结保留、不覆盖,不再用于任何新分析。

一个基准位点(`P42737-2`,同工型后缀 accession)在参考蛋白组 FASTA 中不存在
(该 FASTA 只含 canonical accession `P42737`,无同工型条目),因此在 v2
聚类文件中没有条目;下游 `_train_val_test_rows` 的
`split_map.get(acc, "train")` 回退逻辑已妥善处理(默认分入 train,不报错),
这是基准数据本身的既有特征(同工型 accession 在阳性表中出现,但参考蛋白组
只登记 canonical 序列),不是本次聚类引入的缺陷,影响 1 个阳性位点。

**条件 2/5 用真实簇重新验证**(`scripts/score_release.py --clusters
protein_clusters_v2.tsv` → `scripts/validate_external.py`,重跑时把
AlphaFold 结构注册表临时钉回冻结的 7-结构官方状态,确保只有聚类这一个
变量改变,跑完立即恢复 2,006-结构工作区状态,SHA256 校验一致):

| 条件 | v1(singleton 占位,官方) | v2(真实 MMseqs2 聚类,验证) | 结论是否改变 |
|---|---|---|---|
| 2 效应 CI | PXD006140 +0.0395 [−0.0032, +0.0947](跨零);PXD024061 +0.0766 [+0.0115, +0.1538](不含零) | PXD006140 +0.0395 [**−0.0049**, +0.0931](跨零);PXD024061 +0.0766 [+0.0128, +0.1495](不含零) | 否——同一折失败,下界略微更负 |
| 5 单簇支配 | 最大"簇" = 7 行(单蛋白 O03042);移除后保留 77.2% AP;置换 p=0.001 | 最大簇 = **197 行,真实 30% 同源家族(`Q3E937`)**;移除后 AP **不降反升**(0.0974→0.0990,保留 101.7%);置换 p=0.001 | 否——通过,且现在是一个真正能检测同源家族级别问题的检验 |

条件 2 的点估计(+0.0395/+0.0766)与条件 4 的结构增益(不受影响——
它是配对的逐 fold×seed 比较,不是簇分组 bootstrap)在 v1/v2 之间完全
一致,这是应然的:只有重采样/分组的 key 变了,底层模型打分没变。
**如实结论**:singleton 占位文件并没有掩盖同源驱动的乐观偏差——真实
聚类复现了相同的 STOP 相关结论,条件 5 现在建立在一个真正有意义的
家族级检验上(此前的"簇"永远不可能把超过一个蛋白的 Cys 行分到一起)。

**补充轨道的代价(诚实记录,非缺陷)**:用真实簇重跑 `pu_ranker_cluster_v1`
(→ `pu_ranker_cluster_v2.yaml`/`split_config_v2.yaml`)得到的数字**比 v1
更不稳定,不是更干净**:

| 消融 | v1(singleton,n=5) | v2(真实簇,n=5) |
|---|---|---|
| seq_structure | 0.0837 ± 0.0132(1.8x) | 0.1472 ± 0.0896(3.1x,极不稳定) |
| no_study_context | 0.0435 ± 0.0119(0.9x) | 0.1071 ± 0.1511(2.3x,单 seed 离群主导) |
| sequence_only | 0.0504 ± 0.0076(1.1x) | 0.0559 ± 0.0073(1.2x) |
| seq_esm | 0.0399 ± 0.0057(0.8x) | 0.0494 ± 0.0115(1.0x) |
| no_accessibility | 0.1695 ± 0.3052 ⚠️(单 seed 离群,非信号) | 0.0355 ± 0.0030(0.7x) |
| full | 0.0342 ± 0.0037(0.7x) | 0.0341 ± 0.0001(0.7x) |
| no_plddt | 0.0342 ± 0.0038(0.7x) | 0.0341 ± 0.0001(0.7x) |

机制已查明(非猜测):v2 的测试切分里,最大的真实同源家族(cluster
`Q3E937`,179 个成员落入 test 分区)完整包含 **5 个不同的阳性旁系同源
蛋白**(Q9LZ91、Q9LJJ0、F4JI14、Q84X54、O48583);另有 `Q9M9C5`
(155 成员,3 个阳性旁系同源)、`A0A1P8BD58`(105 成员,2 个)、
`Q9M9P4`(33 成员,3 个)。这些同源家族在 v1 singleton 文件下会被
**随机独立打散**到 train/val/test 三个切分(每个蛋白自成一簇);在 v2
真实聚类下,整个家族**作为不可分割单元**整体分进同一切分——这正是
`_build_singleton_cluster_file` 在结构上无法检测的"同源家族级别的
乐观偏差"风险(见下方旧版缺口描述)。390 个阳性分布在相对较少的大型
旁系同源家族中,一旦某个大家族整体落入 test,该切分的结果就由"模型
对这一个家族学得好不好"主导,逐 seed 训练随机性因此被放大成表面上的
"不稳定"——这是**评估变严格的代价**,不是聚类流程的缺陷,也不应被
解读为"真实聚类更差"。

**遗留登记(供历史参照,已解决)**:此前(v1 阶段)条件 2/5 的
"cluster bootstrap / 簇支配"测量只能在**蛋白粒度**上进行(同一蛋白
的所有 Cys 行同进同出),**不能**检测同源家族级别的乐观偏差,Split A
`pu_ranker_cluster_v1` 的"簇切分"也实为随机蛋白切分——上述 v2 重跑
已解决此项(§3.3 上文)。留一研究主轨道全程不受影响(按研究切分,
不依赖簇文件)。

### 3.4 双轨评估结果(详表见判定记录)

| 轨道 | 切分 | 最优消融 | mean test AP | 基准率倍数 |
|---|---|---|---|---|
| 主轨道(Gate 2 证据) | 留一研究 | seq_structure | 0.102 | 2.1x |
| 补充轨道 v1(文献可比,**禁入 Gate 2**) | 蛋白切分(singleton 占位簇) | seq_structure | 0.084 | 1.8x |
| 补充轨道 v2(文献可比,**禁入 Gate 2**) | 蛋白切分(真实 MMseqs2 簇) | seq_structure | 0.147 ± 0.090(极不稳定) | 3.1x |

主轨道与两条补充轨道对全部稳定消融方向一致:**"严格跨研究切分掩盖了
信号"这一假设被本基准否定**。补充轨道 v1 中 `no_accessibility` 的
0.170(std 0.305)与补充轨道 v2 中 `no_study_context` 的 0.107
(std 0.151)均为单 seed 离群(机制见 §3.3:真实簇下,单个大型同源
家族整体落入某一 seed 的评估路径可主导该 seed 的 AP),非信号,不得
作为消融排名依据。

### 3.5 Gate 2 完整判定(3/5)

| # | 条件 | 结果 | 要点 |
|---|---|---|---|
| 1 | ≥2 独立研究优于基线 | **FAIL(结构性)** | 两折均优于基线(2/2),但研究不独立(配置冻结) |
| 2 | 效应 CI 不含零 | **FAIL(实测)** | 相对 pu_logistic 的配对簇 bootstrap 差值:PXD006140 折 +0.0395 [−0.0032, +0.0947] 跨零;PXD024061 折 +0.0766 [+0.0115, +0.1538] 不含零 |
| 3 | 机制恢复非训练泄漏 | PASS | 0 泄漏,3 个独立验证单元(新增同物种 AtG6PD6) |
| 4 | 结构增益 | PASS | 配对运行级差值均值 +0.044,CI [+0.024, +0.064];但子集层面不可归因(§3.2) |
| 5 | 非单簇驱动 | PASS | 最大"簇"=7 行(单蛋白),移除后保留 77% AP;置换 p=0.001;已用真实 MMseqs2 簇重新验证(§3.3),最大真实同源家族(197 行)移除后 AP 不降反升,结论更稳 |

两个失败条件性质不同:条件 1 是**数据来源锁**(只有新的独立研究能解开);
条件 2 是**实测脆弱性**(模型对 PU 基线的优势在其中一折依赖蛋白组成,
不是纯粹的模型效应)。

### 3.6 已知机制控制恢复(当前登记程序,取代旧表)

恢复百分位由修复后的共享 PU 打分器产生
(`scripts/evaluate_known_controls.py`:一个 pu_logistic 同时给未标记参照
与控制位点打分;旧的单类别 logistic 程序在数学上退化、无法运行,其在
`docs/zhang_collaboration_brief.md` 第 1 页留下的百分位表**作废**,以本
表为准)。控制登记表位于 `src/plantpersulf/evaluation/known_controls.py`
(2026-07-22 自脚本迁入 src,登记语义新增 `control_species` 与
`in_benchmark_as` 两列)。

| 控制 | 状态 | 百分位 | 判定 |
|---|---|---|---|
| SlWRKY6 Cys396 | mapped(跨物种) | 18.9% | 未恢复 |
| SlERF.D2 Cys35 | mapped(跨物种) | 96.5% | 恢复 |
| AtG6PD6 Cys159 | mapped(**同物种**,西北农林独立实验室,doi:10.1111/nph.19188) | 28.6% | 未恢复 |
| PAD3 Cys440 | mapped(**同物种**,裴瑶/靳老师独立实验室,doi:10.1111/pce.70593) | 9.8% | 未恢复 |
| BRG3 | unmappable(番茄蛋白组无匹配) | — | — |
| ERF.D3 | position_shift(编号偏移 13 位,待作者确认) | — | — |

零训练泄漏(AtG6PD6 Cys159 与 PAD3 Cys440 在 benchmark_v1 中均为
unlabeled 行,登记前已核验;同物种控制在打分前从打分器训练样本中显式
剔除)。1/4 恢复率为如实结果:该打分器只有 2 维序列特征,SlWRKY6、
AtG6PD6、PAD3 低恢复反映的是特征集信息量的边界,不构成对机制的否定。
PAD3 9.8% 是四个已映射控制中最低值,与该模式一致,非异常。旧表中
SlWRKY6 81.9% 的数字出自无法复现的退化程序,不得再引用。

**PAD3 Cys440 发现背景(2026-07-23,P1 同物种数据集检索的副产物)**:
本次检索目标是寻找新的独立蛋白质组级拟南芥位点级数据集(见 §5,结论:
未找到)。检索过程中系统排查了 PRIDE 全部 persulfidation/sulfhydration
关键词命中项目(41 条)与 Europe PMC 200 篇 2020–2026 文献,发现裴瑶/
靳老师课题组(山西大学系)自 2013 年起独立发表了 20 余篇植物 H2S 信号
论文,与 Seville Romero/Gotor/Aroca 网络**零作者重叠**。该课题组多篇
拟南芥/十字花科论文报道了具体过硫化位点(BraATO2 Cys416、peach SnRK1α
Cys419/430/505、PAD3 Cys440 等),其中仅 PAD3(拟南芥、CYP71B15/Q9LW27、
490aa)有明确 biotin-switch 验证的单一位点且物种与基准一致,故登记为
控制;其余为番茄/桃/白菜等跨物种候选,或未给出单一位点编号(如 MPK4——
8 个 Cys 逐一突变后过硫化信号均未消失,无法定位单一位点),暂不登记。

**附:PXD072300 同论文正交验证控制(2026-07-23,不计入上表)**

PXD072089 论文同时沉积了 PXD072300(10 个纯化重组水稻蛋白的体外
persulfidation 验证,PEAKS 搜索输出)。解析器
`src/plantpersulf/proteomics/pxd072300_recombinant_sites.py` 逐蛋白定位
目标条目(污染物常常排名更高,不能按行序取,需按 AC 匹配)、修正
`Positions` 列的前置残基偏移(该列给出的是肽段**前一个**残基的位置,
肽段本身起始于该值 + 1——按字面值读会整体错位 1 个残基)、只接受带
`Carbamidomethyl[C]` 标注的持久化位点,再用水稻参考蛋白组坐标核验。

结果:**23 个核验位点、8 个蛋白**(MDHAR3/4/5、FBA1/3、PFP、RPI、TAL 各
2–5 个位点);ALDP 该实验中无 Cys 证据;TKT/转酮醇酶因重组构建体序列与
现有 UniProt 条目均不能 100% 匹配(最佳候选仅覆盖 60/84 条肽段)判定
unmappable,不强行归入。

**这批位点明确不是独立单元**:与 PXD072089 同一篇论文、同一实验室,是
正交验证(单蛋白体外 persulfidation 反应,而非全蛋白质组 MS 发现)而非
独立复现。用途仅限于水稻跨物种迁移轨道(§5.3)的论文内部一致性检验——
检查 Arabidopsis 训练模型对这些同论文生化验证位点的打分是否合理,不得
作为 Gate 2 条件 1/3 的证据。

### 3.7 统计稳健性发现汇总

1. 结构分支增益在运行级稳健,但不能归因到结构覆盖阳性(§3.2)——结构分支
   对 96% 无结构行同时充当正则化先验。
2. 冻结 ESM-2 分支在此数据规模下**有害**(full 0.048 < seq_structure 0.102):
   1280 维冻结特征在 390 阳性下只增加方差。
3. 种子集成(seed-ensemble)汇总 AP 0.097(簇 bootstrap CI [0.049, 0.108])
   高于单 seed 均值——集成可稳定 ranker。
4. 全部结论在 5 个固定 seed 上报告,无最佳 seed 挑选。

---

## 4. 翻转 Gate 2 的差距清单

| 条件 | 当前状态 | 需要什么 |
|---|---|---|
| 1 独立性 | 结构性 FAIL | ≥1 个**不同实验室**的位点级研究进入基准。注意:已登记的两个候选 PRIDE 数据集**都不能**解锁此条件——PXD039999 作者即 Romero/Gotor 同实验室(纯蛋白级公开材料);PXD035795 的方法 PDF 托管于 IDUS Seville(同课题组体系),且其 DCP/NBF 标签需经审阅的"标签→persulfidation 位点"映射规则(Route M)。真正的独立来源=张华团队番茄数据,或经新一轮文献/PRIDE 检索发现的非 Seville 数据集 |
| 2 效应 CI | 一折跨零(真实同源簇下重测确认,§3.3) | 更多阳性;独立研究加入后效应口径重新计算 |
| 3 恢复 | PASS | 维持:新控制位点继续登记、继续排除训练 |
| 4 结构增益 | PASS(运行级,基于冻结的 7-结构状态) | AlphaFold 全蛋白组批量下载已完成(98% 覆盖,§5.5.2),但**尚未**用于任何官方重新验证——需独立、走完整 TDD 流程的 `pu_ranker_v2` 冻结实验,当前未处理(`docs/phase_f_gate2_decision.md` "P4") |
| 5 单簇 | PASS(已用真实同源家族粒度重测,§3.3) | 维持:新数据加入后随基准重测 |

---

## 5. 数据请求清单(张华团队合作,衔接 `docs/zhang_collaboration_brief.md`)

按对 Gate 2 的杠杆排序:

1. **番茄位点级硫巯基化数据**(任何化学方法):位点表(蛋白登录号+Cys 位置)
   与对应原始数据访问方式。用途:第一个跨物种、跨实验室的位点级验证单元;
   同时直接服务 Phase G 番茄迁移。
2. **ERF.D3 编号确认**:论文报告 Cys115/Cys118;UniProt A0A3Q7ESP9 中对应
   Cys128/131/136/139(偏移 13 位)。需确认是否为信号肽/异构体编号差异——
   未确认前该控制不参与任何恢复声明。
3. **BRG3 序列信息**:番茄蛋白组中无预期长度匹配条目;需要作者使用的蛋白
   序列或登录号,否则该控制永久 unmappable。
4. **任何已知的独立实验室位点级公开数据集线索**(不同化学方法优先:
   BTD、CyMPL、qPerS-Ros 等)——一个即可把条件 1 从结构锁变为可测量。
   仓库内候选的局限已查明:PXD039999 为同实验室数据(可补充同框架
   位点,不解锁独立性);PXD035795 属 Seville 体系且需 Route M 映射审阅
   (可作为第二化学方法的补充,独立性待证)。
5. **合作实验设计建议**(若产生新数据):保留时间切分可能性(先冻结训练集、
   后产出验证集),位点级质谱证据 + 逐位点生化验证子集,以便同时充当
   独立研究与金标准。

工程侧配套(本仓库内完成,不依赖合作方):MMseqs2 真聚类
(§3.3,**2026-07-24 已完成**)、AlphaFold 全蛋白组结构批量下载
(§5.5.2 已完成,基于此的官方重新验证仍未处理,见 P4)、新一轮非
Seville 位点级数据集的文献/PRIDE 检索(P1,2026-07-22 已完成一轮,
结果见 §5.1)。

### 5.1 系统文献检索结果(2026-07-22,PRIDE/iProX + 文献,截至 2026-07)

**结论:植物领域不存在非 Seville 的、公开的、全蛋白质组位点级
persulfidation 数据集。** 可用的独立数据集全部是跨物种的;植物内的非
Seville 研究均为单蛋白机制研究(各含 1–2 个验证位点)。

| 数据集 | 物种 | 实验室 | 化学 | 位点级规模 | 沉积 | 判定 |
|---|---|---|---|---|---|---|
| PXD063170(Nat Commun 2025) | 稻瘟菌 M. oryzae | Xiao-Lin Chen(独立) | IAA-PEO-biotin | 全蛋白质组位点级 | PRIDE ✓ | **候选 A → 已建成并首跑(§5.2)** |
| PXD038309(Nat Chem Biol 2023) | 人(MPST 研究) | Dick, Heidelberg(独立) | DCP-Bio1 | 蛋白级补充表(70 蛋白)+ PSM 级 .msf(2.8GB);**无发表位点表** | PRIDE ✓ | 候选 B → **降级**:位点需 Route M 式 .msf 派生(§5.2) |
| PXD072089(PNAS 2026) | 水稻 O. sativa | 谢彦杰,南京林业大学(独立) | NM-biotin + DTT 选择性洗脱 | 全蛋白质组位点级,897 位点/646 蛋白(坐标核验) | PNAS 补充数据 ✓ | **候选 A2 → 已建成并首跑(§5.3)** |
| PXD0701xx(Extremophiles 2026) | 深海古菌 T. aciditolerans | J. Yang(独立) | 化学蛋白组 | 204 位点/171 蛋白 | iProX ✓ | 跨域压力测试(生物学偏远,备选) |
| PXD043969(New Phytol 2023) | 拟南芥+番茄 G6PD | Jisheng Li, 西北农林(独立) | MS 验证 | 2 个验证位点(AtG6PD6 Cys159、SlG6PDC Cys155) | iProX ✓ | 不足作研究级证据;**AtG6PD6 Cys159 已登记为同物种、独立实验室、生化验证控制位点(§3.6,2026-07-22)** |
| 小麦 TaATG6c(Stress Biol 2026) | 小麦 | Xiaojing Wang, 西北农林(独立) | LC-MS/MS | 2 个验证位点 | **无数据库沉积** | 不可用 |
| 小鼠限食 persulfidome(Nat Commun 2021) | 小鼠多组织 | Bithi et al.(独立) | — | 全蛋白质组 | 沉积号未确认 | 候选 C(确认沉积后可并入) |

补充:Sul-BertGRU(Bioinformatics 2025, btaf078)的训练集为 2,705 阳性位点,
源自 iCysMod 数据库(icysmod.omicsbio.info)——与 pCysMod 同属
omicsbio.info 体系的文献整理库;正文未报告物种组成,应用导向为人类
疾病(心血管、神经退行)。其评估为随机蛋白级 80/20 切分 + 10 次重复,
无同源控制(2026-08-10 原文核对:docs/compete/btaf078.pdf)。再次印证
"非 Seville 植物位点级数据"在全领域都稀缺,SOTA 同样依赖混合物种
文献位点。据此新增同行可比轨道 pu_ranker_protein_split_v1(随机蛋白
切分+10 次重复,与 Sul-BertGRU 同口径;不得进入 Gate 2),见
docs/superpowers/specs/2026-07-21-phase-f-pu-ranker-design.md 评估轨道章节。

**执行状态(2026-07-22)**:候选 A 的"跨物种迁移"轨道已建成并首跑
(拟南芥 PU 训练 → 对真菌蛋白组 Cys 打分 → 回收其 persulfidation 位点;
设计、登记、结果与如实解读见 §5.2);候选 B 经核查无发表位点表,降级
(§5.2)。该轨道**不**直接解锁条件 1 的"植物跨研究"原义;它测量的是
"所学决定因素的跨物种迁移",措辞必须如此表述,Gate 2 v2 配置需逐轴
记录独立性(实验室✓ 化学✓ 物种✓)。植物跨研究声明在真正的独立植物
数据集出现前维持 STOP 与降级措辞。

### 5.1a 第二轮检索(2026-07-23):确认无新独立蛋白质组级数据集,发现新控制位点

用户要求继续沿 P1 方向找同物种独立数据集。第二轮检索方法:PRIDE 关键词
搜索 `persulfidation`/`sulfhydration`/`persulfidome` 三个变体去重后共
**41 个项目**(全部人工核对物种+实验室归属);Europe PMC 检索
`persulfidation AND (plant OR Arabidopsis) AND site`,2020–2026,**200
篇文献**(标题+作者逐一扫描,重点核查作者列表与 Seville 团队
Romero/Gotor/Aroca 的重叠情况)。

**结果——4 个 PRIDE 候选全部排除**:
- PXD061767(TRXo1 去过硫化,2025):作者含 Aroca A、Romero LC、Gotor C——
  Seville 网络成员,非独立。
- PXD039852(菜豆根瘤,2024):labPI = Luis C. Romero(IBVF-CSIC/塞维利亚
  大学)——同一实验室。
- PXD019802(ATG4 自噬,2020):labPI 同为 Luis C. Romero——同一实验室。
- PXD005168(2017):`otherOmicsLinks` 直接指向 PXD006140,是同一篇论文
  (doi:10.1093/jxb/erx294)的另一个存档 ID,已计入基准,非新数据。

**文献检索同样未发现新的蛋白质组级位点表**——2020–2026 年全部植物
persulfidation 文献均为单蛋白/少蛋白靶向机制研究(生化验证 1–3 个具体
Cys 位点),没有一篇是全蛋白质组 MS 位点级数据集。**这确认并加固了
2026-07-22 P1 首轮的结论:公开数据库中不存在非 Seville 的植物内独立
全蛋白质组位点级 persulfidation 数据集**——条件 1 的结构性阻塞没有
新数据可以解除。

**副产物:发现裴瑶/靳老师课题组(山西大学系)是真正独立于 Seville 的
拟南芥 H2S 信号课题组**——自 2013 年起发表 20 余篇论文,与
Romero/Gotor/Aroca 网络零作者重叠。虽然都是单蛋白靶向研究(不能填补
条件 1 缺口),但其中 PAD3 Cys440(2026,Plant Cell Environ)有明确
biotin-switch 验证的单一位点,已登记为第二个同物种独立实验室控制
(§3.6)。其余候选(BraATO2 Cys416/Brassica rapa、SnRK1α
Cys419+430+505/桃、BraFLCs/白菜)因物种非拟南芥,归入跨物种控制候选池,
暂未登记(与 SlWRKY6/SlERF.D2 同类,如需要可后续补充)。

### 5.2 跨物种迁移轨道:PXD063170 首跑结果(2026-07-22)

**设计**(`scripts/validate_cross_species.py`,CPU 确定性;全部输入经
supplementary 注册表 SHA256 核验,xlsx 为登记真源、TSV 为其派生物):

- **训练**:全量拟南芥 benchmark_v1(390 阳性 + 1:20 未标记子采样,
  seed 12345;无内部留出——评估完全在外部物种上进行)。模型 = 冻结发布臂
  `structure_ranker:seq_structure`;基线 = `pu_logistic`;5 个固定 seed(0–4)。
- **评估**:PXD063170 发表的 CSE_OE/WT 位点表(Nat Commun 2025 MOESM3),
  经 MG8 参考蛋白组(EnsemblFungi release-62)逐残基坐标核验:
  **1,482 位点全部通过、0 丢弃**。背景 = MG8 蛋白组全部 71,823 个
  半胱氨酸(PU 语义,从不作硬阴性),1:20 子采样(seed 12346)→
  31,122 行,基准率 4.76%。
- **结构分支**:MG8 无登记结构,31,122/31,122 行全部 mask——本轨道测的
  是"序列表示 + 结构门控"的跨物种迁移,不涉及结构增益归因。

**结果**(`results/cross_species/pxd063170_transfer_v1/`):

| 指标 | seq_structure(发布臂) | pu_logistic(基线) |
|---|---|---|
| 逐 seed AP(seed 0–4) | 0.041 / 0.038 / 0.067 / 0.071 / 0.035 | 0.048(全部 5 seed) |
| seed 集成 AP | **0.0635**(基准率 1.33 倍;置换 p=0.001) | 0.0482(≈基准率) |
| 配对蛋白级 bootstrap 差值(集成) | **+0.0153 [+0.0115, +0.0198],不含零** | — |
| recall@50 / recall@500 | **0.0** / 0.8% | — |

**如实解读**:

1. **信号非随机但弱**:集成层面相对基线的优势经蛋白级重采样后 CI 不含零,
   置换检验 p=0.001;但绝对富集仅为基准率 1.33 倍(同框架留一研究为
   2.1 倍),且 top-50 回收为零。
2. **逐 seed 不稳定**:5 个 seed 中 3 个(0、1、4)模型低于基线——集成
   优势由少数 seed 贡献,不得表述为"模型在真菌上稳定有效"。
3. **与 Gate 2 判定一致**:同实验室、同化学、同物种框架内学到的信号,
   迁移到不同物种 + 不同化学(IAA-PEO-biotin)后大幅衰减。该结果支持
   维持 STOP 与降级措辞,不构成任何跨研究/跨物种声明的证据。
4. **轨道价值定位**:候选组织信号的负向/弱向边界测量;评估机器(解析器、
   注册、核验、打分)已就位,未来独立数据出现即可复用。

**候选 B(PXD038309)降级依据(同日核查)**:逐一检查其 Nature 附件
MOESM3–8,均为图 1–6 源数据(roGFP2 OxD 时间序列、蛋白级丰度表等),
**不含位点级 persulfidation 表**;PRIDE 档案仅 4 个 .raw + 1 个 2.8GB
.msf(Proteome Discoverer)。位点表须从 .msf 按 Route M 式"标签→位点"
映射规则自行派生并经受审阅,在规则审阅完成前不入轨。

### 5.3 跨物种迁移轨道:PXD072089 水稻 persulfidome 首跑结果(2026-07-22,v5 更新 2026-07-23)

**数据集**(Xie et al. 2026, PNAS, doi:10.1073/pnas.2608150123):

独立实验室(谢彦杰,南京林业大学)、独立化学(NM-biotin 烷基化 Cys-SH
与 Cys-SSH → DTT 选择性切割 -SSH 二硫键 → 仅释放 persulfidated
肽段进行 LC-MS/MS)、独立物种(水稻 Japonica + Indica)。论文声称 1,691
位点/1,122 蛋白;经 UniProt 2025–2026 年大规模清理(约 55% TrEMBL
条目直接删除、非合并)后,坐标可核验的位点为 **929 个**(646→671 蛋白,
见下文三源联合)。

**数据获取与解析**(`src/plantpersulf/proteomics/pxd072089_sites.py`):

三个独立来源联合(2026-07-23 新增第三源):
- **SD01**(PNAS 补充数据,肽段级 MaxQuant 鉴定表,1,573 行):单 Cys 肽段的
  Position 列即 persulfidation 位点;多 Cys 肽段(179 条)无法单独解析。
- **SD04**(PNAS 补充数据,位点级表,656 行):`-SSH Cys position` 列明确
  标注每个 persulfidation 位点——包括多 Cys 肽段的逐位点独立标注。
- **SS-all-peptides.tsv**(PRIDE PXD072089 提交者存放的完整 MaxQuant
  peptides.txt 导出,3,853 行,发现于核查 PXD072089/PXD072300/PXD072035
  三个关联沉积时):SD01 的预筛选前版本,单 Cys 肽段坐标验证后比 SD01
  的单 Cys 覆盖多出 **32 个位点**(816 个 vs SD01 的 784 个)。

联合策略(三源取并集,取自 `parse_pxd072089_sites(..., ss_all_peptides_path=...)`):

1. SD01 单 Cys 肽段 → 784 个坐标核验唯一键。
2. SD04 全部行 → 396 个坐标核验唯一键。
3. SS-all-peptides 单 Cys 肽段(`C Count == 1`)→ 816 个坐标核验唯一键。
4. 取并集:**929 个唯一 (accession, position) 键 / 671 蛋白**。
   来源分布:`sd01+ss`=501、`sd01+sd04+ss`=283、`sd04`=113、`ss`=32(仅
   SS 覆盖,SD01/SD04 均未采用——这 32 个是纯粹的净增量)。
5. 坐标核验:水稻参考蛋白组 UniProtKB `taxonomy_id:4530`(Japonica +
   Indica,102,397 蛋白、655,127 Cys)。核验规则:accession 存在、position
   在范围内、`sequence[position-1] == 'C'`。

丢弃统计:SD01 567 行 accession 缺失、0 行坐标错配、179 条多 Cys 跳过;
SD04 255 行 accession 缺失、5 行坐标错配;SS-all-peptides 619 行 accession
缺失、0 行坐标错配、192 条多 Cys 跳过。

**UniParc 序列恢复核查(2026-07-23,负结果,重要边界澄清)**:被删除的
483 个 accession 中 455 个(94%)通过 UniProt UniParc 归档(按旧 accession
反查 UPI)成功找回完整序列;但用这些恢复序列重新核验 SD01/SD04/
SS-all-peptides 后,**没有任何一个来自被删 accession 的位点通过验证**——
因为这些 accession 对应的肽段本就未被论文/沉积表判定为 persulfidation
位点(它们只是"检测到的肽段",不是"发表的位点")。结论:UniProt 序列
删除发生在位点判定的**下游**,序列恢复不能反推位点恢复。929 是三个已知
来源能给出的坐标可核验上限;论文声称的 1,691 与 929 之间的差距只能通过
向作者索取原始 `modificationSpecificPeptides.txt` 或等价的位点级 MaxQuant
输出来缩小,而非通过更多序列检索。

**设计**(`scripts/validate_cross_species_rice.py`,CPU 确定性;全部输入经
supplementary 注册表 SHA256 核验):

- **训练**:全量拟南芥 benchmark_v1(390 阳性 + 1:20 未标记子采样,
  seed 12345;无内部留出)。模型 = 冻结发布臂
  `structure_ranker:seq_structure`;基线 = `pu_logistic`;5 个固定 seed(0–4)。
- **评估**:929 阳性 + 1:20 未标记子采样(seed 12346) →
  19,509 行,基准率 4.76%。
- **结构分支**:水稻无登记结构,19,509/19,509 行全部 mask——本轨道测的
  是"序列表示 + 结构门控"的跨物种迁移,不涉及结构增益归因。

**结果**(`results/cross_species/pxd072089_transfer_v2/`;v1=897 位点结果
保留作历史记录,不覆盖):

| 指标 | seq_structure(发布臂) | pu_logistic(基线) |
|---|---|---|
| 逐 seed AP(seed 0–4) | 0.041 / 0.041 / 0.069 / 0.074 / 0.039 | 0.045(全部 5 seed) |
| seed 集成 AP | **0.0618**(基准率 1.30 倍;置换 p=0.001) | 0.0452(≈基准率) |
| 配对蛋白级 bootstrap 差值(集成) | **+0.0166 [+0.0108, +0.0299],不含零** | — |
| recall@50 / recall@500 | **0.0** / 2.6% | — |

929 位点相对 897 位点的结果**在噪声范围内完全一致**(AP 0.06183 vs
0.06184,delta 几乎相同)——32 个新增位点没有实质改变轨道结论。

**如实解读**:

1. **与 PXD063170 模式一致**:水稻跨物种迁移 AP 0.0618(1.30x 基准率)
   与稻瘟菌 0.0635(1.33x)近乎一致。差值 CI 均不含零,置换 p<0.001——
   信号非随机但弱于同一物种的留一研究(2.1x)。
2. **逐 seed 不稳定性相同**:5 seed 中 3 个(0、1、4)模型低于基线,
   集成优势由少数 seed 贡献。不得表述为"模型在水稻上稳定有效"。
3. **独立轴满足但不构成 Gate 2 条件 1 证据**:实验室 ✓ / 化学 ✓ / 物种 ✓,
   三条独立轴客观上全部独立于 Seville 基准;但本轨道仍为"跨物种迁移",
   条件 1 原义为留一植物研究的交叉验证。将 PXD072089 纳入 3 研究
   留一法(拟南芥→拟南芥→水稻)的 pu_ranker_v2 配置是未来翻转条件 1
   的路径,但该配置需要将水稻阳性位点映射到拟南芥同源空间(SD02 提供
   了拟南芥同源映射,但仅 192/929 位点落在 1:1 同源蛋白上,简单窗口
   比对能干净映射到拟南芥 Cys 坐标的只有约 25 个——统计力不足以单独
   成折)并注册到基准中——此为独立任务,当前不可行。
4. **数据价值**:929 个坐标核验位点(671 蛋白)在 655,127 Cys 背景上构成
   已知最大的跨物种、跨化学、跨实验室的植物 persulfidation 评估集。
   为未来条件 1 的 3 研究留一法验证提供了数据基础,但同源映射的统计力
   限制意味着该路径当前不可行(见上)。
5. **PXD072300 同论文控制**(§3.6 附):8 个重组蛋白、23 个体外验证位点
   可用于检验模型对这些论文内部生化验证位点的打分是否合理,但明确不是
   条件 1/3 的独立单元。

---

### 5.4 跨物种轨道对比

| 轨道 | 物种 | 化学 | 位点 | 蛋白 | 背景 Cys | AP | x 基准率 | recall@50 |
|---|---|---|---|---|---|---|---|---|
| PXD063170 | 稻瘟菌 | IAA-PEO-biotin | 1,482 | MG8 蛋白组 | 71,823 | 0.0635 | 1.33x | 0.0 |
| PXD072089(v2,+SS) | 水稻 | NM-biotin+DTT | 929 | 671 | 655,127 | 0.0618 | 1.30x | 0.0 |

两条独立轨道(不同物种、不同化学、不同实验室)的跨物种迁移 AP 几乎一致,
均表现为弱信号(1.3x 基准率) + 逐 seed 不稳定性(3/5 低于基线)。两条
**都不**构成 Gate 2 证据——措辞约束见 §7。

### 5.5 框架重构(2026-07-23):从"预测"转向"保守性/趋同性"——两层正向发现

**动机**:§5.2–5.4 的跨物种迁移轨道被 Gate 2 结构性拒绝为*预测*证据(条件 1
要求植物内跨研究,跨物种迁移不满足)。但同一批数据(拟南芥 benchmark_v1、
PXD072089 水稻、PXD063170 稻瘟菌——三个独立实验室、独立化学、独立物种)
可以回答一个不同、同样定量、但不需要预测有效性的问题:**这些独立数据集
的持久化靶点选择,是否在家族/通路层面和结构层面呈现出超出偶然的趋同性？**
这是标准的比较基因组学富集分析方法论,不是评估标准的降级。

#### 5.5.1 家族级趋同性(PANTHER 同源子家族富集)

方法:用 UniProtKB `xref_panther` 批量注释三个物种全蛋白组的 PANTHER
子家族(`PTHR#####:SF#`,精细同源组代理);稻瘟菌的 MG8 基因 ID
(`MGG_#####T0`)通过 UniProt 条目的 `Gene Names (ORF)` 字段桥接到
UniProt accession(12,656/12,755 = 99.2% 覆盖)。对每对物种,限定在
两个蛋白组共有的 PANTHER 子家族全集内,构建"该家族在物种 A 是否被持久化"
×"在物种 B 是否被持久化"的 2×2 列联表,单侧 Fisher 精确检验(超几何分布
精确计算,不依赖 scipy,与本项目现有 `permutation.py`/`effect_size.py`
风格一致)。

结果(`scripts/analyze_cross_species_conservation.py` →
`results/cross_species_conservation/conservation_v1.json`):

| 物种对 | 共享子家族 | 共同持久化 | 期望值 | 富集倍数 | OR | 原始 p 值 | Bonferroni | BH q 值 |
|---|---|---|---|---|---|---|---|---|
| 拟南芥 × 水稻 | 9,188 | 21 | 7.11 | 2.96x | 3.34 | 7.98×10⁻⁶ | **7.98×10⁻⁶ → 2.39×10⁻⁵**(显著) | 1.20×10⁻⁵(显著) |
| 拟南芥 × 稻瘟菌 | 2,204 | 20 | 12.63 | 1.58x | 1.90 | 1.75×10⁻² | **1.75×10⁻² → 5.25×10⁻²**(**不显著**) | 1.75×10⁻²(显著) |
| 水稻 × 稻瘟菌 | 2,320 | 67 | 30.22 | 2.22x | 3.60 | 1.31×10⁻¹² | **1.31×10⁻¹² → 3.93×10⁻¹²**(显著) | 3.93×10⁻¹²(显著) |

**多重检验校正(2026-07-23,自我审稿发现并修正)**:三对检验的原始
p 值未经校正就报告为"全部显著",这在统计上不严谨——三个检验构成一个
检验族,需要族错误率或错误发现率校正。补做 Bonferroni(保守,控制族
错误率,`m=3`)与 Benjamini-Hochberg(FDR,`scripts/
analyze_cross_species_conservation.py` 内 `bonferroni_and_bh_correction`,
不依赖 scipy)后:**拟南芥 × 稻瘟菌这一对在 Bonferroni 校正下跌出
0.05 阈值**(0.0525),但在 BH 校正下仍然显著(q=0.0175)——两种方法在
这一对上意见不一致,必须如实同时报告,不能只挑对结论有利的一个。**如实
结论:3 对中有 2 对在两种校正方法下都稳健显著(拟南芥×水稻、水稻×
稻瘟菌);拟南芥×稻瘟菌这一对的显著性方法依赖(method-dependent),
不得笼统表述为"三对全部显著"。**

三方交集(在三个物种都被持久化的子家族,与两两检验的显著性校正无关,
是独立于假设检验之外的描述性重叠)有 8 个,经 UniProt PANTHER
EntryName 人工核对身份,**7/8 是核心氧化还原/中心碳代谢酶**:果糖-1,6-
二磷酸酶(卡尔文循环)、甘氨酸脱氢酶与丙氨酸-乙醛酸转氨酶(光呼吸)、
NADH-泛醌氧化还原酶 75kDa 亚基(线粒体复合体 I)、山梨醇脱氢酶、
L-苏糖酸脱氢酶(醛缩酶超家族)、肽基脯氨酰顺反异构酶;唯一非典型的是
RCC1(染色体凝集调控)。这与"H2S 通过过硫化调控 NAD(P)H 结合结构域活性
半胱氨酸"这一领域内已有机制假说吻合,但作为跨三个独立物种的统计学佐证,
其强度应理解为"2/3 物种对方法无关地显著 + 1 个小规模三方交集的描述性
观察",而非无条件的"全面证实"。

#### 5.5.2 结构层趋同性(AlphaFold 可及性上下文)

**结构覆盖率的重大提升**:Gate 2 判定记录(§3.2)中 390 个基准阳性只有
14 个(3.6%)有登记结构,原因是早期登记使用的 AlphaFold 模型版本(v4)
已经全面 404——AlphaFold DB 已升级到 v6。改用 v6 后,对三个物种全部
持久化蛋白(拟南芥 350、水稻 671、稻瘟菌 1,032,经 PANTHER 同一套
UniProt 桥接,去重后 2,048 个 UniProt accession)批量查询,**命中率
接近 100%**(抽样 50/50 全部命中)。批量下载(`scripts/
download_alphafold_structures_bulk.py`,复用已有的 fail-closed
`plantpersulf.download.alphafold` 模块,SHA256 逐条登记入
`data/registry/alphafold_structures.tsv`)最终注册 **2,006 个结构
(98.0%)**
(2,048 中的绝大多数;不可用的是真实的"无预测结构"事实,如同工型
accession,而非下载失败)。

**两条独立测量路径,都必须报告(2026-07-23 自我审稿新增)**:

1. **可及性代理**(`features/structure.py::contact_number_proxy`,原有
   方法):半径内其他残基 Cα 数量,数值越高代表越包埋——明确标注**不是**
   严格 SASA 计算。
2. **真实 Shrake-Rupley SASA**(`features/sasa.py`,新增,响应审稿意见
   "应换成真正的 SASA"):纯 Python 实现,不引入 numpy/scipy/Biopython
   依赖,与本项目一贯的"从零实现统计/几何算法"风格一致。用完整重原子
   坐标(C/N/O/S,取自 PDB 元素符号字段第 77–78 列,而非从原子名猜测)+
   Bondi(1964)范德华半径 + 1.4Å 水探针半径,通过确定性 Fibonacci 球面
   采样(无随机数,给定采样点数完全可复现)判定每个采样点是否被邻近原子
   遮蔽,得到每原子 SASA,再按残基求和。空间网格做邻居剪枝,单个 Cys 残基
   计算约 1.5–2.7ms,全数据集(约 16,600 个 Cys)不到 1 分钟。同时报告
   **全残基 SASA** 与 **SG 原子(过硫化直接作用的硫原子)专属 SASA** 两个
   指标。因为组内比较全部是 Cys-vs-Cys(同一残基类型),不需要跨残基类型
   的相对 SASA 归一化,原始 Å² 数值已经直接可比。

两条路径都对每个物种取全部"至少含一个持久化 Cys 且有结构"的蛋白,提取
该蛋白**每一个** Cys 的指标值。持久化位点标 `positive`,同蛋白内其他
Cys 标 `unlabeled`(PU 语义:未检出≠阴性,与全项目一致,绝不标
`negative`)。用置换检验(复用 `evaluation/permutation.py` 引擎,双侧,
`|mean_diff|` 作为检验统计量)判定组间差异是否超出偶然。

**可复现性修复(2026-07-23)**:早期实现按 Python `set` 迭代顺序处理蛋白
列表,而 Python 字符串哈希按进程随机化(`PYTHONHASHSEED`),导致同一
`--seed` 参数在不同进程运行间给出略有差异的 p 值(实测拟南芥 p 值曾在
0.034/0.102/0.104/0.1019 之间波动)。已修复为按 accession 排序后处理,
两次独立进程运行现在给出逐位精确一致的结果(已验证:两次运行 JSON 输出
`diff` 结果为空)。

结果(`results/cross_species_conservation/structural_context_v1.json`,
2,006/2,048 = 98.0% 结构覆盖率,可复现性修复后的最终数字):

| 物种 | 蛋白数 | 阳性 Cys | 指标 | mean_diff | 方向 | p 值(双侧置换) |
|---|---|---|---|---|---|---|
| 拟南芥 | 344 | 384 | 代理(contact) | +1.023 | 更包埋 | 0.101(不显著) |
| 拟南芥 | 344 | 384 | SASA(全残基) | −1.19 Å² | 更包埋 | 0.551(不显著) |
| 拟南芥 | 344 | 384 | SASA(SG 原子) | −0.17 Å² | 更包埋 | 0.881(不显著) |
| 水稻 | 669 | 922 | 代理(contact) | +1.016 | 更包埋 | 0.011 |
| 水稻 | 669 | 922 | SASA(全残基) | −5.22 Å² | 更包埋 | **0.001** |
| 水稻 | 669 | 922 | SASA(SG 原子) | −2.14 Å² | 更包埋 | **0.002** |
| 稻瘟菌 | 984 | 1,387 | 代理(contact) | +1.318 | 更包埋 | **0.001** |
| 稻瘟菌 | 984 | 1,387 | SASA(全残基) | −4.34 Å² | 更包埋 | **0.001** |
| 稻瘟菌 | 984 | 1,387 | SASA(SG 原子) | −1.77 Å² | 更包埋 | **0.002** |

**核心发现:真实 SASA 与代理指标给出完全一致的定性结论**——三个物种
方向全部一致(更包埋),显著性格局也完全复制:水稻、稻瘟菌两个指标都
显著(SASA 甚至比代理指标更显著),拟南芥两个指标都不显著。这不是巧合
凑出来的一致,而是**两种独立实现的测量方法(一个 Cα 邻居计数代理,一个
全原子 Shrake-Rupley 精确计算)在同一批结构数据上给出相同答案**——这
把此前"依赖单一代理指标"的技术弱点转正为交叉验证证据。SG 原子专属 SASA
(过硫化直接修饰的硫原子本身)显示的模式与全残基 SASA 一致,进一步支持
"这是真实的侧链暴露度差异,不是残基整体大小或构象的伪影"。

**如实结论(不变)**:**方向 3/3 一致**(全部"更包埋"),**显著性
仍是 2/3**(水稻、稻瘟菌两个指标都显著;拟南芥两个指标都不显著)——
拟南芥效应量方向相同、量级合理,但样本量/方差下未过常规阈值,更可能是
统计力不足而非真实缺失效应,但不得声称拟南芥"显著"。这个方向本身并不
违反直觉——过硫化通常发生在酶活性位点的催化/调控半胱氨酸上,这类残基
常位于底物结合口袋内、部分包埋但化学反应活性高(pKa 被周围残基调制),
而非游离暴露在蛋白表面。这与 §5.5.1 识别出的靶点身份(NADH 脱氢酶、
山梨醇脱氢酶等 NAD(P)H 结合酶的活性位点半胱氨酸)自洽。

#### 5.5.3 两层发现的统一叙事

家族级趋同(§5.5.1,**2/3 物种对方法无关地显著**——拟南芥×水稻、
水稻×稻瘟菌在 Bonferroni 与 BH 校正下都稳健;拟南芥×稻瘟菌方法依赖,
Bonferroni 下不显著)+结构级趋同(§5.5.2,**方向 3/3 一致、显著性
2/3**——水稻、稻瘟菌显著,拟南芥方向一致但未达常规阈值;且该结构级
发现现在由**两条独立测量路径交叉验证**:可及性代理指标与真实
Shrake-Rupley SASA 给出完全一致的定性格局,不是单一指标的偶然结果)
共同支持一个机制假说:**H2S 介导的过硫化靶点选择由结构/功能上下文
(特定酶家族的活性位点半胱氨酸)决定,而非局部序列基序**。这与
§5.2–5.4 测得的弱序列迁移信号(1.3x 基准率、逐 seed 不稳定)方向一致,
但两层证据本身强度都是中等而非压倒性(家族级 2/3 经校正显著,结构级
2/3 显著、且经真实 SASA 复现)——如实报告为"支持性但非压倒性,且结构
层已通过独立方法交叉验证"的机制证据,不得表述为"证实"或笼统的
"三对/三物种全部显著"。三条证据线(序列迁移弱、家族富集中等、结构方向
一致且双方法复现)互相印证,拼出一个连贯、诚实标注局限性、且关键结构
发现已经过方法学交叉验证的生物学论断。

**性质澄清(措辞约束,§7 同步)**:§5.5 的全部数字都是**关联性证据**
(富集/趋同性),不是预测性证据,不构成 Gate 2 任何条件的证据,也不
改变 GATE2_STOP 的判定。

#### 5.5.4 归因偏差检验:UniProt 删除是否按氧化还原功能选择性丢弃蛋白?(2026-07-23,自我审稿 Reviewer 2 关切回应)

**动机**:§5.3 指出 PXD072089 原始报告的位点因 UniProt 2025–2026 清理大量
accession 被直接删除而无法坐标核验。一个具体、需要正面回应的担忧是:如果
被删除的 accession **系统性地偏向非氧化还原功能蛋白**,那么"幸存"下来
可核验的蛋白集合会被人为地富集氧化还原/脱氢酶类功能——这将成为
§5.5.1 家族级富集与 §5.5.2 结构级包埋发现的一个混杂解释,而非真实的
生物学趋同信号。

**方法**(`src/plantpersulf/evaluation/attrition_bias_check.py`,
`scripts/check_attrition_bias.py`):SD01/SD04 的原始行(核验之前)都携带
一个 UniProt 风格的自由文本"Protein description"字段,这个字段是原始
MS 检索时捕获的,即使该 accession 之后被删除也依然保留在补充表里——
这正是判断"被删除的蛋白原本是什么功能"所需要的资源。对 SD01+SD04 中
出现的每个唯一 accession(与坐标是否核验无关,只看 accession 字符串本身
是否非空),用一份固定、公开的氧化还原关键词表(dehydrogenase、oxidase、
oxidoreductase、reductase、peroxidase、thioredoxin、glutathione、
ferredoxin、cytochrome、superoxide dismutase、NAD(P)、flavin、
aminotransferase 等 15 个关键词,大小写不敏感,只匹配 `OS=`
之前的功能描述文本,不匹配物种名部分)分类是否"氧化还原相关"。按该
accession 是否能在当前水稻参考蛋白组(`uniprot_rice_v1.fasta`,与
`pxd072089_sites.py` 判定坐标核验时用的同一份蛋白组)中解析,分为
kept(存活)/dropped(已删除)两组,构建 2×2 列联表,做双侧 Fisher
精确检验(超几何分布精确枚举,`fisher_exact_two_sided`,不依赖 scipy,
测试集含教科书对照表验证正确性)。

结果(`results/cross_species_conservation/attrition_bias_v1.json`):

| 分组 | 唯一 accession 数 | 氧化还原相关 | 非氧化还原 | 氧化还原占比 |
|---|---|---|---|---|
| kept(存活于当前蛋白组) | 687 | 74 | 613 | 10.77% |
| dropped(已从 UniProt 删除) | 439 | 39 | 400 | 8.88% |

SD01+SD04 唯一 accession 共 1,126 个,蛋白级归因率 38.99%(注意:这与
§5.3 引用的"约 55% TrEMBL 条目"是不同口径的数字——§5.3 的 55% 是对
UniProt 2025–2026 年该批清理规模的一般性描述,本检验的 38.99% 是
SD01+SD04 原始行覆盖的 accession 全集中,无法在当前水稻参考蛋白组解析
的比例,两者分母不同,不应直接比较)。kept 组氧化还原占比(10.77%)
数值上略高于 dropped 组(8.88%),OR=1.24,但**双侧 Fisher 精确检验
p=0.360,远未达到 0.05 阈值——两组的氧化还原关键词构成在统计上没有
可检测的差异**。

**如实结论**:未发现 UniProt accession 删除存在氧化还原功能选择性偏差
的证据。§5.5.1/§5.5.2 中氧化还原/脱氢酶类靶点的富集,不能用"删除偏向
移除了非氧化还原蛋白从而人为抬高幸存集合氧化还原占比"这一混杂机制来
解释——该检验没有支持这个混杂假说。**局限性(如实标注,不隐藏)**:
(a) 这是一个描述性关键词分类,15 个关键词的覆盖面不是穷尽的,存在假
阴性(未被关键词捕获的氧化还原酶)与假阳性(关键词命中但功能实际不是
经典氧化还原,例如某些含"NAD(P)"字样的结合结构域蛋白)的双向噪声,
但该噪声对 kept/dropped 两组是均匀施加的(同一套关键词、同一套匹配规则),
不构成差异性偏差;(b) 该检验只覆盖 SD01+SD04,未纳入 SS-all-peptides.tsv
(§5.3 第三来源)的独立 accession 集合,如果后续需要更完整的归因宇宙,
应扩展本脚本纳入该文件;(c) 该检验是蛋白级的(以 accession 为单位),
不是位点级的,与 §5.3 位点级 929/1,691 的最终归因率是不同的分析单元,
两者不能相互替代。

#### 5.5.5 文献新颖性核查:三方交集 8 个家族此前是否已被单独报道为过硫化靶点?(2026-07-23,真实文献检索,回应自我审稿 Reviewer 2"跨物种趋同有多新颖"的关切)

**动机**:§5.5.1 描述性报告了三个物种交集的 8 个 PANTHER 子家族,其中
7/8 是氧化还原/中心碳代谢酶。模拟审稿(`/nature-reviewer`)的 Reviewer 2
指出:如果这些酶家族本来就是持久化研究领域内众所周知的"热门"过硫化
靶点,那么"三物种趋同"的新颖性主张就需要更谨慎地限定范围。本节是
**真实的**(非模拟)文献检索,逐一核查 8 个家族身份,而非笼统断言。

**逐家族检索结果**(WebSearch + WebFetch,2026-07-23):

| 家族(三方交集身份) | 是否已有单蛋白/靶向过硫化文献报道 | 关键来源 |
|---|---|---|
| 甘氨酸脱氢酶(GDC 复合体 P/H 蛋白)| **是,且是该领域内已充分表征的目标**——一篇 2023 年综述(Aroca/García-Díaz/García-Calderón 团队,即 Seville Romero/Gotor 网络)专门总结光呼吸酶的过硫化修饰,明确列出 GDC-H 蛋白(GDH1/GDH2)"across multiple proteomic studies"被持久化过硫化 | García-Díaz et al. 2023, *J. Exp. Bot.* 74(19):6023, doi:10.1093/jxb/erad291 |
| 丙氨酸-乙醛酸转氨酶(AGT,乙醛酸转氨酶家族)| **部分是**——同一篇综述报道了功能相邻的乙醛酸转氨酶亚型 GGAT1/GGAT2(谷氨酸-乙醛酸转氨酶)与 SGAT(丝氨酸-乙醛酸转氨酶)被过硫化,但未特指 AGT/AGX1 本身;三者同属乙醛酸转氨酶反应网络,身份高度相关但非同一基因 | 同上 |
| 果糖-1,6-二磷酸酶(FBPase,卡尔文循环)| **否(过硫化层面)**——该 Cys(Cys153/Cys173 区域)是已充分表征的**氧化还原调控位点**,但已发表机制是硫氧还蛋白/S-亚硝基化调控,未检索到 H2S 过硫化的直接报道;是"已知氧化还原敏感位点,新 PTM 类型"而非全新靶点 | Jacquot et al. 1997, *FEBS Lett.*;综述 doi:10.3390/antiox10101631 |
| NADH-泛醌氧化还原酶 75kDa 亚基(线粒体复合体 I)| **未检索到**该亚基或复合体 I 整体的过硫化文献报道 | — |
| 山梨醇脱氢酶 | **未检索到**该酶的过硫化文献报道(检索到功能相关的琥珀酸脱氢酶被 H2S 过硫化激活以诱导气孔关闭的报道,提示脱氢酶类整体是易感底物,但非同一酶) | Plant and Soil (2024),doi:10.1007/s11104-024-06837-x(琥珀酸脱氢酶,非山梨醇脱氢酶) |
| L-苏糖酸脱氢酶(醛缩酶超家族)| **未检索到**过硫化文献报道 | — |
| 肽基脯氨酰顺反异构酶(环孢素蛋白/PPIase)| **未检索到**该酶类过硫化的具体文献报道 | — |
| RCC1(染色体凝集调控)| **未检索到**——包括人类同源物在内均未检索到过硫化报道;§5.5.1 已将其标注为"唯一非典型"身份,文献检索确认其确实缺乏任何已知过硫化关联 | — |

**背景性发现(重要,影响新颖性主张的表述方式)**:本项目自己的拟南芥
阳性集合(`benchmark_v1`)本身来自 Aroca et al. 2017 的过硫化蛋白组学
研究(`doi:10.1093/jxb/erx175`);该论文已经报告"约 5% 的拟南芥全蛋白组
(2,015 个蛋白)在基线条件下可被过硫化",并明确指出**中心碳代谢是被
持久化富集的功能类别之一**。这意味着"氧化还原/碳代谢酶是过硫化热点"
这一定性观察,**在拟南芥单物种层面并不新颖**——它已经是 Aroca 团队
2017 年原始论文自带的结论,不是本项目的新发现。

**如实的新颖性定位**:
1. **不新颖的部分**:8 个交集家族中至少 2 个(甘氨酸脱氢酶复合体、
   乙醛酸转氨酶网络)已被 Seville 同一研究网络的后续综述明确列为
   拟南芥过硫化的已知靶点;FBPase 是已知的氧化还原敏感位点(经由其他
   PTM 机制)。"氧化还原/中心碳代谢酶富集"这一定性格局本身已经隐含在
   本项目所用的拟南芥源数据自己的原始论文结论里。
2. **检索未能证伪、但也无法证实为全新的部分**:NADH-泛醌氧化还原酶
   75kDa 亚基、山梨醇脱氢酶、L-苏糖酸脱氢酶、肽基脯氨酰顺反异构酶、
   RCC1 五个家族,未检索到任何单蛋白过硫化文献报道——但这是一次
   有限的 WebSearch/WebFetch 检索(非系统性数据库比对,如 dbGSH 或
   PerSpect 等专用过硫化位点数据库未逐一核对),不能排除遗漏,只能
   诚实地表述为"未发现已发表报道",不能表述为"首次发现"。
3. **本项目真正可以主张的新颖性,范围明确限定为**:据本次检索所及,
   **尚未发现任何已发表工作用独立采集的跨物种过硫化蛋白质组数据、
   在同源家族层级上做过统计学趋同性检验**(即 §5.5.1 的 Fisher 精确
   检验富集分析、§5.5.2 的跨物种结构包埋方向一致性)——这是方法学/
   统计学层面的新颖性主张,而不是"发现了全新的过硫化靶点家族"这一
   靶点发现层面的新颖性主张。**这个区分必须在任何对外文本(投稿、
   摘要)中保持清晰,不得把"统计趋同性证据新颖"混同为"生物学靶点
   本身新颖"。**

**局限性**:本节检索基于 WebSearch/WebFetch 工具,覆盖 PubMed/期刊
索引到的公开可及摘要与综述,未做穷尽式文献计量学检索,亦未核对专用
过硫化位点数据库;"未检索到"应理解为"检索范围内未发现",不是"经系统
排查确认不存在"。

---

## 6. 可复现性

```bash
# 主轨道 + 补充轨道(冻结配置,GPU 可选)
python scripts/run_experiment.py --config configs/experiments/pu_ranker_v1.yaml
python scripts/run_experiment.py --config configs/experiments/pu_ranker_cluster_v1.yaml
# 补充轨道 v2(真实 MMseqs2 簇,§3.3;本仓库的可编辑安装会让裸 import
# 解析到 C:/PlantPersulf-Code 姊妹仓库,运行前需 PYTHONPATH=src)
python scripts/run_experiment.py --config configs/experiments/pu_ranker_cluster_v2.yaml
# 逐位点分数再生(CPU 确定性) + 控制恢复 + Gate 2 判定
python scripts/score_release.py --config configs/experiments/pu_ranker_v1.yaml \
    --output-dir results/external_validation/pu_ranker_v1/scored
python scripts/evaluate_known_controls.py --output results/known_controls/recovery_v1.json
python scripts/validate_external.py --model-release pu_ranker_v1 \
    --model-tag "structure_ranker:seq_structure" \
    --scored    results/external_validation/pu_ranker_v1/scored/model.tsv \
    --scored-baseline results/external_validation/pu_ranker_v1/scored/baseline.tsv \
    --scored-ablated  results/external_validation/pu_ranker_v1/scored/ablated.tsv \
    --recovery results/known_controls/recovery_v1.json
# 条件 2/5 真实簇重验证(§3.3;重跑前必须把 alphafold_structures.tsv 钉回
# 官方冻结的 7-结构状态,否则会与结构覆盖率变量混杂,见 P4)
python scripts/score_release.py --config configs/experiments/pu_ranker_v1.yaml \
    --clusters data/processed/clusters/protein_clusters_v2.tsv \
    --output-dir results/external_validation/pu_ranker_v1_clusterv2_check/scored
python scripts/validate_external.py --model-release pu_ranker_v1_clusterv2_check \
    --model-tag "structure_ranker:seq_structure" \
    --clusters data/processed/clusters/protein_clusters_v2.tsv \
    --scored    results/external_validation/pu_ranker_v1_clusterv2_check/scored/model.tsv \
    --scored-baseline results/external_validation/pu_ranker_v1_clusterv2_check/scored/baseline.tsv \
    --scored-ablated  results/external_validation/pu_ranker_v1_clusterv2_check/scored/ablated.tsv \
    --recovery results/known_controls/recovery_v2.json
# 跨物种迁移轨道(PXD063170,CPU 确定性,措辞=跨物种迁移,非 Gate 2 证据)
python scripts/validate_cross_species.py \
    --output-dir results/cross_species/pxd063170_transfer_v1
# 跨物种迁移轨道(PXD072089 水稻,929 位点=SD01+SD04+SS-all-peptides 三源并集,
# CPU 确定性,措辞=跨物种迁移,非 Gate 2 证据)
python scripts/validate_cross_species_rice.py \
    --output-dir results/cross_species/pxd072089_transfer_v2
# 保守性/趋同性框架(§5.5,关联性证据,非 Gate 2 证据)
python scripts/analyze_cross_species_conservation.py \
    --output results/cross_species_conservation/conservation_v1.json
python scripts/download_alphafold_structures_bulk.py \
    --accessions data/raw/references/persulfidated_protein_accessions_v1.txt
python scripts/analyze_structural_context.py \
    --output results/cross_species_conservation/structural_context_v1.json
```

`protein_clusters_v2.tsv` 的生成(需要 Linux 主机跑 MMseqs2,本机没有)
不在上面的命令列表内——完整的两跳远程流程见
`docs/ops/remote_a100_via_lab_jump_host.md`,核心命令是在目标服务器上跑
`mmseqs easy-cluster arabidopsis_ref_proteome_v1.fasta clusterRes/arabidopsis_v2 tmp
--min-seq-id 0.3 -c 0.5 --cov-mode 0 --threads 32`。

全部输入 SHA256 记录于 `results/external_validation/pu_ranker_v1/manifest.json`;
跨物种轨道输入 SHA256 记录于
`results/cross_species/pxd063170_transfer_v1/cross_species_summary.json` 与
`results/cross_species/pxd072089_transfer_v2/cross_species_summary.json`
(`pxd072089_transfer_v1` 为 897 位点的历史记录,保留不覆盖);
硬件限制(>2500 aa 蛋白跳过 ESM 提取,影响 ≤3 阳性、0.8%)见判定记录。

## 7. 措辞约束

- Gate 2 = STOP 期间,任何外部文本(论文、合作书、汇报)禁止出现预测能力
  类措辞;判定记录页首的降级声明必须随行。
- Split A(补充轨道,v1 或 v2)数字可用于文献口径对比,但必须附带对应
  配置文件(`pu_ranker_cluster_v1.yaml` 或 `pu_ranker_cluster_v2.yaml`)的
  limitation 原文,且不得作为跨研究/跨实验室/跨物种外推的证据。v1 是
  singleton 占位聚类(§3.3),v2 是真实 MMseqs2 聚类但逐 seed 方差大
  (§3.3 已查明机制:大型同源家族整体落入单一切分),引用 v2 数字时
  必须同时给出 std,不得只报均值。
- 条件 2/5 引用时,若要点明"已用真实同源簇验证",必须同时给出 v1/v2
  两组数字(§3.3 表格),不得只引用其中一组暗示"验证了但不给数据"。
- 跨物种迁移轨道(§5.2 PXD063170 稻瘟菌、§5.3 PXD072089 水稻)数字只能以
  "cross-species transfer(跨物种迁移)"表述,并同时给出逐 seed 不稳定与
  recall@50=0 的事实;不得据此声称模型在真菌/水稻/其他物种上"有效",
  不得作为 Gate 2 任何条件的证据。
- §5.5 家族级/结构级趋同性数字只能以"关联性证据(associational evidence,
  富集/趋同性)"表述,不得表述为"预测"或"预测能力";`contact_number_proxy`
  必须原样称呼为"可及性代理指标",不得称为"溶剂可及性"或
  "solvent accessibility"(该模块文档明确其非严格 SASA 计算);不得作为
  Gate 2 任何条件的证据,且不改变 GATE2_STOP 判定。
- §5.5.1 家族级富集**不得**表述为"三对物种全部显著"或"3/3 显著"——
  经 Bonferroni/BH 双重校正后如实为"2/3 物种对方法无关地显著,1 对
  (拟南芥×稻瘟菌)方法依赖"(§5.5.1/v8 变更)。任何引用必须同时给出
  未校正 p 值、Bonferroni 校正后 p 值、BH q 值三者,不得只挑单一数字。
- §5.5.2 结构级趋同性引用时必须同时给出**两条独立路径**的数字(可及性
  代理 + 真实 Shrake-Rupley SASA,`features/sasa.py`),不得只引用其中
  一条来暗示"已验证";真实 SASA 结果可称"solvent-accessible surface
  area"或"SASA",但代理指标(`contact_number_proxy`)的称呼限制不变。
  两条路径均为"方向 3/3 一致、显著性 2/3(水稻、稻瘟菌)",不得表述为
  "3/3 显著"(v9 变更)。
- §5.5.1 三方交集 8 个家族的"新颖性"表述必须限定为**统计趋同性方法学
  新颖**,不得表述为"发现全新过硫化靶点家族"——§5.5.5 真实文献核查
  确认其中至少 2 个家族(甘氨酸脱氢酶复合体、乙醛酸转氨酶网络)已是
  Seville 团队自己 2023 年综述中列出的已知拟南芥靶点,且"氧化还原/
  碳代谢富集"这一定性观察已隐含在本项目拟南芥源数据自己的原始论文
  (Aroca et al. 2017)结论中;未检索到既往报道的家族(NADH-泛醌
  氧化还原酶 75kDa、山梨醇脱氢酶、L-苏糖酸脱氢酶、PPIase、RCC1)只能
  表述为"检索范围内未发现已发表报道",不得表述为"首次发现"或"证实
  为新靶点"(§5.5.5/v10 变更)。
- §5.5.4 归因偏差检验(kept vs. dropped 蛋白氧化还原关键词构成,
  p=0.360)可用于反驳"氧化还原富集是 UniProt 删除伪影"这一具体质疑,
  但不得反向表述为"证明富集是真实生物学信号"——该检验只能排除一种
  特定混杂机制,不能充当富集发现本身的正面证据(§5.5.4/v10 变更)。
