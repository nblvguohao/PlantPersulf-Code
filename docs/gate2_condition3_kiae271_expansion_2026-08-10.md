# Gate 2 条件 3(recovery_no_leakage)扩容:kiae271 差异硫巯基化对照池

> 日期:2026-08-10 | 分支:phase-a-persulfidation-site-evidence
> 前置文档:`docs/zhang_lab_published_data_audit.md` §8.3(数据发现)、`docs/phase_f_gate2_decision.md`(冻结的 Gate 2 判定,本文档不修改它,只作为条件 3 的补充证据记录)

## 1. 做了什么

把 `known_controls.py` 里手工维护的 9 个对照(SlWRKY6、SlERF.D2、BRG3、SlWRKY71、RNF144b、ERF.D3、AtG6PD6、PAD3)之外,新增了一条系统性、可复现的对照来源:kiae271 论文自己发表的补充材料 Supplementary Dataset S1(标题写"SlWRKY6 位点分析"但实际是完整的 119 位点差异硫巯基化表)。

- 新模块 `src/plantpersulf/proteomics/kiae271_sites.py`:解析该表,做坐标验证(accession 存在、位置在范围内、该位置确为 Cys)+ 定位概率过滤(≥0.75)+ 至少一个条件强度非零,失败原因逐条计数,不做任何静默修复。
- 新脚本 `scripts/evaluate_kiae271_controls.py`:复用 `evaluate_known_controls.py` 完全相同的方法学(在拟南芥 benchmark 上训练一个共享 PU 打分器,对每个对照位点算相对未标注分布的百分位排名),保证与现有 9 控制对照的结果可直接比较。
- 数据来源已登记:`data/registry/supplementary_sources.tsv` + `configs/supplementary_sources_v1.yaml`(study_accession `KIAE271_SUPPL`,SHA256 校验通过),原始文件位于 `data/raw/supplements/KIAE271_SUPPL/kiae271_DSs.xlsx`。
- 10 个新增单元测试(`tests/scientific/test_kiae271_sites.py`)覆盖:位置对齐(Positions within proteins 与 Leading proteins 的正确对应,而非默认取第一个)、定位概率过滤、零强度剔除、坐标不匹配剔除、重复位点去重。

## 2. 结果

119 行原始数据 → **99 个坐标验证通过的位点**(88 个不同蛋白);20 行因位置对齐歧义(8)、accession 未匹配(3)、坐标不符(1)、定位概率 <0.75(8)被剔除并计数。

| 分层 | n | 恢复数(百分位>50) | 恢复率 | 平均百分位 |
|---|---|---|---|---|
| 全部 | 99 | 48 | 48.5% | 49.8 |
| lcd_gain(仅 LCD1-OE 检出) | 43 | 22 | 51.2% | 48.4 |
| wt_only(仅 WT 检出) | 40 | 19 | 47.5% | 50.9 |
| both(两组都检出) | 16 | 7 | 43.8% | 50.6 |

二项检验:48/99 vs 50% 随机基线,p = 0.84(95% CI [38.6%, 58.3%])——**与随机无统计学差异**。三个分层(H2S 诱导获得、组成性、两者皆有)之间也没有实质差异。

## 3. 结论与边界

**这不改变 Gate 2 条件 1 的判定**(跨物种数据结构上不能替代同一 benchmark 内独立研究留一验证,与水稻/真菌迁移轨道此前的结论一致)。它做的是把**条件 3 的统计功效提高了一个数量级**(样本量从 6-9 个手工挑选对照,变成 99 个系统性抽取、真实定量、坐标验证过的独立位点),结果与此前用小样本对照观察到的"≈随机"现象完全吻合,但现在有了统计显著性意义上的确认(而不是 6 个样本时那种"看起来像随机但样本太小说不清"的状态)。

这本身就是一个更有力的诚实负结果:在现有 2 维序列特征下,模型对独立物种/独立实验室/独立化学方法产出的硫巯基化位点没有超过随机水平的判别力,N=99 时binomial p=0.84 排除了"小样本恰好不显著"这一解释。

## 4. 产物清单

- `src/plantpersulf/proteomics/kiae271_sites.py`(解析器)
- `scripts/evaluate_kiae271_controls.py`(评估脚本)
- `tests/scientific/test_kiae271_sites.py`(10 个单元测试)
- `data/registry/supplementary_sources.tsv` + `configs/supplementary_sources_v1.yaml`(新增 KIAE271_SUPPL 条目)
- `data/raw/supplements/KIAE271_SUPPL/kiae271_DSs.xlsx`(登记的原始数据,SHA256 `a437f1a941a348ccf15d538bcea193d9feefe29019c5162f0427551416b1b98a`)
- `results/known_controls/kiae271_differential_recovery_v1.json` + `.rows.tsv`(评估输出)

## 5. 下一步(未做,留待决定)

- 是否要把这 99 个位点也接入结构特征(AlphaFold embedding/SASA)评估,看结构分支在这个更大样本上是否仍有微弱改善——这是回应上一轮讨论中"结构层面 PTM crosstalk"角度的自然延伸。
- 是否要把 `evaluate_known_controls.py` 和 `evaluate_kiae271_controls.py` 的结果合并进一份统一的 Gate 2 条件 3 报告(目前是两个独立产物,分别对应"手工精选高置信对照"和"系统性大样本对照"两种不同证据强度,暂不合并陈述,以保持两者可追溯的差异)。
