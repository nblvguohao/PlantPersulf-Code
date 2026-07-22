# 植物硫巯基化数据资源与系统性证据审计(Phase Z)

> 版本: v1 | 日期: 2026-07-22 | 分支: phase-a-persulfidation-site-evidence
> 路线: Gate 2 = **GATE2_STOP** → 按 roadmap 降级路线产出的数据资源交付物
> 机器可核验判定: `results/external_validation/pu_ranker_v1/gate2_decision.json`
> 完整判定记录: `docs/phase_f_gate2_decision.md`

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
结构文件、`protein_clusters_v1.tsv`(当前为 singleton 占位,见 §3.3)。

---

## 3. 系统性证据审计

### 3.1 证据集中度

- 全部 390 个阳性来自**同一实验室、同一 tag-switch 化学、同一物种**;无独立实验室、无第二化学方法(如 BTD/CyMPL)、无第二位点级物种。
- 证据等级单一:全部为 `site_ms`(质谱位点);无逐位点生化验证的独立阳性集可充当同物种金标准。
- 张华团队已发表的番茄机制位点(4 个登记控制)是唯一跨物种、跨实验室证据,其中 2 个可映射、1 个不可映射、1 个位置待确认(§3.6)。

### 3.2 结构覆盖率

仅 7 个结构文件覆盖基准:14/390 阳性(3.6%,6 个蛋白)+ 113 行未标记。
结构分支的实测增益(§3.5 条件 4)因此**不能归因于结构覆盖的阳性本身**;
结构覆盖子集(16 行)上的配对差值 CI 跨零(+0.009 [−0.087, +0.056]),
该事实如实记录,不被运行级增益掩盖。

### 3.3 同源性与聚类控制现状(登记缺口)

`protein_clusters_v1.tsv` 目前是**每蛋白一簇的 singleton 占位文件**
(`_build_singleton_cluster_file` 生成;真 MMseqs2 聚类需要 Linux 主机,
本机无 mmseqs/WSL/conda)。后果如实登记:

- 条件 2/5 的"cluster bootstrap / 簇支配"测量在**蛋白粒度**上进行(同一蛋白
  的所有 Cys 行同进同出),**不能**检测同源家族级别的乐观偏差;
- Split A `pu_ranker_cluster_v1` 的"簇切分"实为随机蛋白切分,只控制蛋白
  重叠,不控制同源性;其数字引用时必须附带此说明;
- 留一研究主轨道不受影响(按研究切分,不依赖簇文件)。
- **补救**:在 Linux 主机运行 `mmseqs easy-cluster`(30% 一致性),产出
  `protein_clusters_v2.tsv`(v1 冻结不覆盖),登记 SHA256 后重跑簇轨道与
  bootstrap。已列入 §5 工作清单。

### 3.4 双轨评估结果(详表见判定记录)

| 轨道 | 切分 | 最优消融 | mean test AP | 基准率倍数 |
|---|---|---|---|---|
| 主轨道(Gate 2 证据) | 留一研究 | seq_structure | 0.102 | 2.1x |
| 补充轨道(文献可比,**禁入 Gate 2**) | 蛋白切分(singleton 簇) | seq_structure | 0.084 | 1.8x |

两轨道对全部稳定消融一致(差距 −0.006…−0.018):**"严格跨研究切分掩盖了
信号"这一假设被本基准否定**。补充轨道中 `no_accessibility` 的 0.170 为单
seed 离群(std 0.305,其余 4 个 seed 0.031–0.047),非信号。

### 3.5 Gate 2 完整判定(3/5)

| # | 条件 | 结果 | 要点 |
|---|---|---|---|
| 1 | ≥2 独立研究优于基线 | **FAIL(结构性)** | 两折均优于基线(2/2),但研究不独立(配置冻结) |
| 2 | 效应 CI 不含零 | **FAIL(实测)** | 相对 pu_logistic 的配对簇 bootstrap 差值:PXD006140 折 +0.0395 [−0.0032, +0.0947] 跨零;PXD024061 折 +0.0766 [+0.0115, +0.1538] 不含零 |
| 3 | 机制恢复非训练泄漏 | PASS | 0 泄漏,2 个独立验证单元 |
| 4 | 结构增益 | PASS | 配对运行级差值均值 +0.044,CI [+0.024, +0.064];但子集层面不可归因(§3.2) |
| 5 | 非单簇驱动 | PASS | 最大"簇"=7 行(单蛋白),移除后保留 77% AP;置换 p=0.001;蛋白粒度警告见 §3.3 |

两个失败条件性质不同:条件 1 是**数据来源锁**(只有新的独立研究能解开);
条件 2 是**实测脆弱性**(模型对 PU 基线的优势在其中一折依赖蛋白组成,
不是纯粹的模型效应)。

### 3.6 已知机制控制恢复(当前登记程序,取代旧表)

恢复百分位由修复后的共享 PU 打分器产生
(`scripts/evaluate_known_controls.py`:一个 pu_logistic 同时给未标记参照
与控制位点打分;旧的单类别 logistic 程序在数学上退化、无法运行,其在
`docs/zhang_collaboration_brief.md` 第 1 页留下的百分位表**作废**,以本
表为准):

| 控制 | 状态 | 百分位 | 判定 |
|---|---|---|---|
| SlWRKY6 Cys396 | mapped | 18.8% | 未恢复 |
| SlERF.D2 Cys35 | mapped | 96.4% | 恢复 |
| BRG3 | unmappable(番茄蛋白组无匹配) | — | — |
| ERF.D3 | position_shift(编号偏移 13 位,待作者确认) | — | — |

零训练泄漏(训练阳性全部拟南芥)。1/2 恢复率为如实结果:该打分器只有
2 维序列特征,SlWRKY6 低恢复反映的是特征集信息量的边界,不构成对该机制
的否定。旧表中 SlWRKY6 81.9% 的数字出自无法复现的退化程序,不得再引用。

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
| 2 效应 CI | 一折跨零 | 更多阳性+真实同源簇文件下的重测;独立研究加入后效应口径重新计算 |
| 3 恢复 | PASS | 维持:新控制位点继续登记、继续排除训练 |
| 4 结构增益 | PASS(运行级) | 扩大结构覆盖(AlphaFold 全蛋白组批量下载),使子集层面可判定 |
| 5 单簇 | PASS(蛋白粒度) | MMseqs2 真聚类(§3.3)后按同源家族粒度重测 |

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

工程侧配套(本仓库内完成,不依赖合作方):MMseqs2 真聚类(§3.3)、
AlphaFold 全蛋白组结构批量下载、新一轮非 Seville 位点级数据集的
文献/PRIDE 检索(P1)。

---

## 6. 可复现性

```bash
# 主轨道 + 补充轨道(冻结配置,GPU 可选)
python scripts/run_experiment.py --config configs/experiments/pu_ranker_v1.yaml
python scripts/run_experiment.py --config configs/experiments/pu_ranker_cluster_v1.yaml
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
```

全部输入 SHA256 记录于 `results/external_validation/pu_ranker_v1/manifest.json`;
硬件限制(>2500 aa 蛋白跳过 ESM 提取,影响 ≤3 阳性、0.8%)见判定记录。

## 7. 措辞约束

- Gate 2 = STOP 期间,任何外部文本(论文、合作书、汇报)禁止出现预测能力
  类措辞;判定记录页首的降级声明必须随行。
- Split A(补充轨道)数字可用于文献口径对比,但必须附带
  `pu_ranker_cluster_v1.yaml` 的 limitation 原文与 §3.3 的 singleton 说明,
  且不得作为跨研究/跨实验室/跨物种外推的证据。
