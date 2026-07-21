# PlantPersulf-Code 后续任务总纲（Post-Gate1 Roadmap）

> **用途**：本文件是 Task 5 readiness gate 发布 STOP 判定之后，对全部后续任务的统一执行总纲。
> 它不替代 `docs/PlantPersulf_Code_TDD_Codex.md`（那是权威规格），而是把规格里的 Task 5→13
> 重新组织为“先解锁 Gate 1、再决策分叉、最后建模/降级”的可执行路线。
>
> **For agentic workers**：REQUIRED SUB-SKILL —— 用 superpowers:executing-plans 或
> subagent-driven-development 逐 Phase 执行；每个 Step 用 `- [ ]` 复选框跟踪。
> 每个功能严格 RED → GREEN → REFACTOR，单独 commit，跑完整完成闸门（unit/scientific/release/ruff/mypy）。

---

## 0. 当前状态基线（写作时的事实）

- 已完成：Task 0（诚信闸门）→ Task 5 的 **readiness 判定**（`benchmark/readiness.py`）。
- 最新 commit：`feat: publish fail-closed Task 5 readiness decision`。
- **Gate 1 判定 = STOP**，来自 `data/interim/benchmark_readiness_v1/`：

  | 指标 | 当前值 | Gate 1 要求 |
  |---|---|---|
  | `eligible_site_count` (`site_ms`/`site_mutagenesis`/`site_biochemical`) | **0** | > 0 |
  | `eligible_study_count` | **0** | **≥ 2** |
  | `coordinate_only_count` (PXD006140，均为甲硫氨酸氧化，非 persulfidation) | 3 | — |
  | `non_site_evidence_count` | 10,787 | — |
  | `unresolved_candidate_count` | 25 | — |
  | **decision** | **STOP** | — |

- **STOP 根因**：四个硫巯基化研究目前只解析了外围文件——
  - PXD006140：仅两个 OMSSA `.mgf.txt`（原始谱图/搜索结果，无位点定位）；
  - PXD035795：仅 identification/protein 级 `peptide.csv`/`proteins.csv`/`mzid`；
  - PXD024061、PXD039999：**只下了 `checksum.txt`**，尚未触达结果文件。
  - 位点级 persulfidated-Cys 阳性真实存在，但主要在**同行评议补充表（Data 等级 C）**或尚未下载的结果文件里。

- **红线约束（不可违反）**：在拿到 ≥2 研究的真实位点证据前，禁止建 benchmark、禁止训模型、禁止把
  identification/protein 级证据“升级”为位点级。宁可 STOP，不伪造合格证据。

---

## 1. 全局架构

```
Phase A  证据补齐冲刺 ──→ 重跑 readiness gate ──→ ┬─ GO  ──→ Phase B~G（建模主线）
（解锁 Gate 1，决定性）                          └─ STOP ──→ Phase Z（降级为数据资源/证据审计）
Phase X  贯穿始终的工程基建（可与 A 并行）
```

**执行顺序硬约束**：
- Phase A 未产出 GO 之前，不启动 Phase B 及之后任何建模任务。
- Phase C 的 split 首次训练前冻结 commit，之后只能建 v2，不得覆盖。
- Phase F 第 10 节结论闸门未过，不得在对外材料写“具有预测价值”。
- 湿实验真实结果返回前，不声称发现新机制。

---

## 2. 技术栈与既有资产

- Python 3.10–3.11，pytest / ruff / mypy(strict)（已就绪，全绿）。
- 复用现有：`provenance/{schema,registry,hashing,audit}.py`、`download/{base,pride,geo,sra,iprox,registered}.py`、
  `proteomics/{metadata,peptide_parser,site_normalizer}.py`、`evidence/*`、`benchmark/readiness.py`。
- 复用现有审计器：`audit_site_output`、`audit_content_output`、`audit_benchmark_readiness`。
- 真实微型夹具机制：`tests/fixtures/real/<ACCESSION>/source_manifest.json`（禁改生物学值）。

---

# Phase A —— 位点级证据补齐冲刺（最高优先级 / 决定项目走向）

**Goal**：把四个硫巯基化研究里真实存在的 Cys 位点级证据，以完整 provenance 合法登记并解析进
`site_normalizer`，使 readiness gate 有机会翻成 GO。

**成功判据**：≥2 个研究产出真实 `site_*` 等级证据，且全部通过既有的
`test_site_coordinates_match_sequence` 与内容审计。

## Task A1：位点证据来源地图与登记方案

**Files:**
- Create: `docs/superpowers/specs/2026-07-2X-evidence-site-acquisition-design.md`
- Create（若采用补充表）: `data/registry/supplementary_sources.tsv`（模式复用 `evidence_methods.tsv`）
- Modify: `configs/data_sources.yaml`（登记新增位点证据文件类别）

**内容要求（逐研究）**：
- [ ] 对 PXD006140 / PXD024061 / PXD035795 / PXD039999 各列出：位点级 persulfidation 证据的**确切载体**
  （PRIDE 带修饰定位的 search result 文件 vs 论文补充表），Data 等级（预期多为 C），官方 URL/DOI，
  license，预期位点数量级，能否合法获取。
- [ ] 为每个可获取来源规划完整 provenance 字段：accession/DOI + 官方 URL + 下载时间 + size +
  SHA256 + license + downloader/parser 版本。
- [ ] 明确标注**不可获取**或**仅 protein-level**的研究，作为后续 STOP 风险证据（进入 Phase Z 材料）。

**验收**：设计文档评审通过；不写任何代码。

## Task A2：每研究位点解析器（严格 TDD）

**Files:**
- Create: `src/plantpersulf/proteomics/site_sources/<study>.py`（或扩展现有 `peptide_parser.py`）
- Modify: `src/plantpersulf/proteomics/site_normalizer.py`
- Create: `tests/scientific/test_<study>_site_evidence.py`
- Create: `tests/fixtures/real/<ACCESSION>/...`（位点证据的真实微型夹具 + `source_manifest.json`）

**RED（示例，逐研究）：**
```python
def test_site_positions_match_uniprot_sequence(real_site_fixture):
    records = parse_persulfidation_sites(real_site_fixture)
    for r in records:
        assert r.evidence_level in {"site_ms", "site_mutagenesis", "site_biochemical"}
        assert sequence_at(r.protein_accession, r.cys_position_in_protein) == "C"

def test_non_persulfidation_modifications_are_not_promoted(real_site_fixture):
    records = parse_persulfidation_sites(real_site_fixture)
    # 甲硫氨酸氧化 / carbamidomethyl 等不得被当成 persulfidation 位点
    assert all(r.modification_name_raw != "oxidation" for r in records)

def test_unmappable_sites_go_to_conflicts_not_silently_fixed(real_site_fixture):
    records, conflicts = parse_with_conflicts(real_site_fixture)
    assert conflicts_are_recorded(conflicts)
```

**GREEN 约束：**
- [ ] 只从已校验真实文件/补充表截取夹具；禁改蛋白/肽段/位点/定量/标签。
- [ ] Cys 位点必须在对应 **UniProt 版本化序列**上精确匹配（复用 `test_site_coordinates_match_sequence`）。
- [ ] 修饰命名正确归一到 persulfidation；甲硫氨酸氧化等**不得升级**。
- [ ] 冲突/不可映射位点进 conflicts 表，不自动纠正；protein-level-only 研究**不生成虚构位点**。

**验收：**
```bash
python -m plantpersulf.cli parse-sites --accession <ACCESSION>
pytest tests/scientific/test_<study>_site_evidence.py -v
```

## Task A3：重跑 readiness gate，复判 GO/STOP（决策点）

**验收：**
```bash
python -m plantpersulf.cli build-benchmark-readiness --version v1
python -m plantpersulf.cli audit-benchmark-readiness --version v1
```
- [ ] 检查 `eligible_study_count`。**这是本 plan 的决策分叉点。**

| 结果 | 条件 | 走向 |
|---|---|---|
| **GO** | ≥2 研究有真实 `site_*` 证据 | Phase B |
| **STOP（持续）** | 尽力后仍 <2 研究 | Phase Z |

---

# Phase B~G —— 建模主线（仅在 Gate 1 = GO 时启动）

> 全部照 TDD Codex 第 7 节推进；此处只给 Phase 级目标、关键测试、闸门与依赖。

## Phase B（Task 5 真正版）：PU benchmark v1
- **产出**：`data/processed/benchmark_v1/{sites,proteins,studies}.parquet` + `conflicts.tsv` + `manifest.json`。
- **关键测试**：无未观测位点被标 `negative`；每个 `positive` 有 `study_accession`+`source_sha256`+合格 `evidence_level`。
- **闸门**：`audit-benchmark --version v1` 通过；label 集合 ⊆ {positive, unlabeled}。

## Phase C（Task 6）：防泄漏 split（首次训练前冻结）
- **产出**：`configs/splits/{cluster,study,species,time,known_mechanism_holdout}_split_v1.yaml`；`benchmark/{cluster_split,study_split,time_split,leakage_audit}.py`。
- **关键测试**：同一 MMseqs2 cluster 不跨 split；已知张华机制位点不在 train；time cutoff 训练前写入并 commit。
- **硬约束**：split 冻结后只能建 v2，禁覆盖。

## Phase D（Task 7）：真实序列/结构特征
- **产出**：`download/{uniprot,alphafold}.py`；`features/{sequence,esm2,structure}.py`；`configs/features/{sequence,structure}_v1.yaml`。
- **科学约束**：结构缺失用显式 mask 不填均值；低 pLDDT 单独分层；坐标映射失败不进结构模型但保留在序列模型。

## Phase E（Task 8）：无学习基线 + 传统 PU 基线
- **模型**：motif/frequency、accessibility、Logistic/PU-Logistic、RF、XGBoost、ESM+linear head。
- **TDD 重点**：标准化只用训练折；超参只用 validation；test 只跑一次；≥5 固定 seed；结果含 seed/split/config hash。

## Phase F（Task 9→10）：Structure-aware PU ranker + 严格外部验证
- **启动门槛**：Phase E 全防泄漏测试通过、≥2 个 leave-study-out fold 可执行、ESM+linear 稳定可复现、无未登记输入。
- **架构（最小）**：sequence + frozen ESM + structure + missingness mask → gated fusion → PU 排序头 → 校准/不确定性头。禁止大图/端到端微调 650M/无法消融模块。
- **外部验证**：每 PXD leave-study-out、cluster bootstrap、time split、已知机制回顾恢复、全 seed。
- **第 10 节结论闸门**：五项全满足才可写“具有预测价值”，否则写降级措辞。→ 对应 **Gate 2**。

## Phase G（Task 11→13）：番茄上下文 + 候选发布
- **Task 11**：PXD051570 + GEO 番茄成熟上下文（蛋白/磷酸化/转录动态，ID 冲突不静默）。→ 对应 **Gate 3**。
- **Task 12**：机制卡 + 去循环候选评分（已知 H₂S 基因只作解释/富集/阳性控制，不作训练标签或硬编码加分）。
- **Task 13**：证据卡 + 冻结 `candidate_release_v1` + `docs/zhang_collaboration_brief.md` + `docs/wetlab_validation_matrix.md` + `candidate-release-v1` tag。

---

# Phase Z —— STOP 持续时的降级交付（诚信兜底）

**触发**：Phase A 证明公开数据确实拿不到 ≥2 研究的位点证据。
**动作**：不硬训模型，按项目自身契约转为——

> **《植物硫巯基化数据资源与系统性证据审计》**

**交付物：**
- [ ] 四研究证据等级全景表（A–X 分级）+ 位点级证据缺口报告。
- [ ] ID 映射冲突量化 + 研究间修饰定义不可比性分析。
- [ ] 可复现的“为何现有公开数据不足以支撑跨研究位点预测”论证（对应第 12 节 Level 1 data-resource 论文路线）。
- [ ] 与张华老师的数据请求清单（第 10 节 10 项真实/合作数据）。

---

# Phase X —— 贯穿始终的工程基建（可与 Phase A 并行）

补齐架构图已规划但仓库尚缺的资产：
- [ ] `README.md`（根目录当前完全缺失——对外协作硬伤）。
- [ ] `Makefile`（第 8 节承诺的 `make metadata/download-registered/benchmark-v1/...` 一键复现全缺）。
- [ ] `environment.yml`（版本锁定环境）。
- [ ] `data/registry/SHA256SUMS` 汇总校验。
- [ ] `configs/` 骨架补齐：`splits/`、`features/`、`experiments/`、`id_mapping.yaml`、`mechanism_cards/`。
- [ ] CI：禁用词扫描（`synthetic`/`fake`/`dummy`/`mock_biology` 科研路径）+ 未登记数据路径扫描。

---

# 附录：Stop/Go 检查点汇总

| Gate | 位置 | GO 条件（摘要） | STOP 走向 |
|---|---|---|---|
| **Gate 1 数据可用性** | Phase A3 | ≥2 研究有真实位点证据 | Phase Z 降级为数据资源/证据审计 |
| **Gate 2 跨研究预测** | Phase F | ≥2 held-out study 稳定富集、非单 cluster 驱动、CI/检验支持、校准可接受 | 候选降级为 evidence integration，不声称通用 prediction |
| **Gate 3 番茄候选价值** | Phase G | ≥1 批候选在适用域、有真实成熟上下文、非全已知基因、可提可证伪实验 | 缩减候选/仅作假设生成 |
| **Gate 4 谈合作** | Phase G 末 | registry+benchmark 过审、split 冻结、有基线+ESM、已知机制恢复可展示、5–20 候选证据卡、限制清楚 | 携 Phase Z 材料谈数据需求 |

---

# 立即执行建议

**只做一件事：启动 Phase A1（证据来源地图）。** 它是唯一能决定“项目能否建模”的动作，其余所有
Phase 都是它的下游。在拿到证据地图之前，Task 6–13 的排期都是空中楼阁。
