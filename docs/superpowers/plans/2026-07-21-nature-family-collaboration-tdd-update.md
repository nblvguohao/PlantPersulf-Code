# Nature-family Collaboration TDD Update Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Update the master TDD with an explicit collaboration value proposition, leakage-safe published controls, two contact gates, and a claim-safe Nature-family readiness ladder.

> **2026-07-21 更新（GO 对齐）**：本 plan 初稿写于 Task 5 仍为 `STOP` 的时点。此后 Phase A 证据补齐工作已用真实、坐标验证的位点级证据（PXD006140 Dataset S3 320 个 + PXD024061 MaxQuant 73 个 = ~393 个 `site_ms`）通过 fail-closed 的 readiness v2 gate，将 **Gate 1（数据可用性）翻为 GO**，并据此构建了 PU benchmark v1、防泄漏 split、序列/ESM-2 特征与 motif 基线（均已单独 commit）。本 plan 现以该 GO 为准做了对齐修订：涉及"当前状态=STOP、只允许 Gate 4A"的事实性表述已更新；但**全部 claim-safety 脚手架（控制分类学、Gate 4A/4B/5、Level 0-4、红线 16-21）保持不变且依然适用**——因为 Gate 1 的 GO 只是数据可用性层面，Gate 2 的跨研究预测价值尚未验证，对外主张仍受这些闸门约束。

**Architecture:** Apply three documentation-only revisions to the master TDD's conceptual boundaries. Preserve Task 0–13 and every scientific-integrity constraint. The Task 5 state is now `GO` at Gate 1 (data availability); reflect that in the current-state note while keeping every claim-safety gate intact.

**Tech Stack:** Markdown, PowerShell, `rg`, Git.

## Global Constraints

- These three **documentation tasks** modify only `docs/PlantPersulf_Code_TDD_Codex.md`; they change no code, tests, configs, data, or task order. (The Phase A implementation that produced the Gate 1 GO is committed separately and is out of scope for these doc edits.)
- Do not create candidate scores or claim mechanisms; do not claim cross-study predictive value until Gate 2 is `GO`.
- Keep all Zhang-laboratory published controls excluded from training.
- Do not imply that AI alone closes the Nature-family gap or guarantee publication.
- The current Gate 1 GO rests on two same-lab, same-chemistry (Romero/Gotor tag-switch), same-species (Arabidopsis) studies — a data-availability GO only; frame it as such everywhere.
- Commit each task separately.

---

### Task 1: Define collaboration value and known-control taxonomy

**Files:**
- Modify: `docs/PlantPersulf_Code_TDD_Codex.md:48-240`
- Reference: `docs/superpowers/specs/2026-07-21-nature-family-collaboration-tdd-update-design.md`

**Interfaces:**
- Consumes: sections 1 and 3.3.
- Produces: section 1.4 plus `control_type`, `mechanism_lineage_id`, and `independent_validation_unit` definitions.

- [ ] **Step 1: Run RED**

```powershell
rg -n -e "合作科学定位" -e "strong_single_site_control" -e "conditional_site_group_control" -e "independent_validation_unit" docs/PlantPersulf_Code_TDD_Codex.md
```

Expected: exit 1 because the definitions are absent.

- [ ] **Step 2: Add section 1.4 before section 2**

Insert this content:

```markdown
## 1.4 合作科学定位与项目增量

基于已公开论文形成的项目工作判断如下：

- 张华老师团队已经具备 LC-MS/MS 位点鉴定、硫巯基化检测、Cys 定点突变、CRISPR/过表达、蛋白互作、转录调控、PTM 串扰和番茄表型验证能力；
- 本项目不以代做 RNA-seq、常规差异分析、富集分析或通用生信服务作为合作价值；
- 本项目拟补充可追踪位点级证据整合、防泄漏跨研究预测、训练前候选冻结、匹配对照和前瞻性盲法验证；
- 单独增加 AI 模型、随机交叉验证或候选排名表，不构成 Nature-family 级生物学创新；
- 联合研究要从逐个单蛋白机制上升到可推广、可预测并经前瞻实验检验的植物硫巯基化规律。

上述判断只用于项目设计，不评价团队未公开能力，也不保证任何期刊结果。
```

- [ ] **Step 3: Replace section 3.3**

Use this exact control table and policy:

```markdown
## 3.3 张华老师团队已发表机制作为训练外阳性控制

已发表机制只能作为 `known_positive_mechanism_card`，不能用于声称新发现，也不能用于解除训练 benchmark 的数据不足。

| 控制 | DOI | 控制类型 | 使用限制 |
|---|---|---|---|
| SlWRKY6 Cys396 | `10.1093/plphys/kiae271` | `strong_single_site_control` | 训练外单点回顾性恢复 |
| SlERF.D2 Cys35 | `10.1111/tpj.70000` | `strong_single_site_control` | 训练外单点回顾性恢复 |
| BRG3 Cys206/Cys212 | `10.1093/plphys/kiad070` | `conditional_site_group_control` | 无单点拆分证据时不得解释为两个独立功能阳性 |
| ERF.D3 Cys115/Cys118 | `10.1093/plphys/kiae560` | `conditional_site_group_control` | 无单点拆分证据时不得解释为两个独立功能阳性 |

磷酸化、泛素化和转录调控位点只能作为 PTM 串扰解释证据。2026 年 SlWRKY6–SlGRF1–SlGIF2 研究（`10.1093/plphys/kiag512`）复用同一 SlWRKY6/H₂S 机制链，不得重复计为新的 `independent_validation_unit`。

处理规则：

1. 登记正文或补充材料、DOI、官方 URL、文件名和 SHA256；
2. 登记 canonical protein accession/version、位点映射、证据层级和冲突；
3. 人工复核并记录复核状态；
4. 登记 `control_type`、`mechanism_lineage_id` 和 `independent_validation_unit`；
5. 所有控制标记 `positive_control_only`、`retrospective_recovery_test` 和 `excluded_from_training`；
6. 控制身份、独立性分组和评价规则在训练前冻结；
7. 控制排名不得用于特征、模型、阈值、超参数或停止条件选择；
8. 同一蛋白、位点或机制谱系的后续论文不得重复计数。
```

- [ ] **Step 4: Run GREEN and commit**

```powershell
rg -n -e "## 1.4 合作科学定位与项目增量" -e "SlWRKY6 Cys396" -e "SlERF.D2 Cys35" -e "BRG3 Cys206/Cys212" -e "ERF.D3 Cys115/Cys118" -e "mechanism_lineage_id" -e "independent_validation_unit" docs/PlantPersulf_Code_TDD_Codex.md
rg -n -e "excluded_from_training" -e "不得用于特征、模型、阈值" -e "不能用于解除训练 benchmark" docs/PlantPersulf_Code_TDD_Codex.md
git diff --check
git status --short
```

Expected: all markers present, no whitespace errors, and only the master TDD modified.

```bash
git add -- docs/PlantPersulf_Code_TDD_Codex.md
git commit -m "docs: define collaboration value and known controls"
```

---

### Task 2: Specify retrospective and prospective validation

**Files:**
- Modify: `docs/PlantPersulf_Code_TDD_Codex.md:1160-1410`

**Interfaces:**
- Consumes: Task 1 control fields.
- Produces: complete recovery reporting, future RED specifications, and a frozen blind-validation package.

- [ ] **Step 1: Run RED**

```powershell
rg -n -e "test_duplicate_mechanism_lineage_is_not_independent" -e "test_predictive_claim_requires_gate2_go" -e "frozen_analysis_plan" -e "test_nature_readiness_requires_prospective_blind_validation" docs/PlantPersulf_Code_TDD_Codex.md
```

Expected: exit 1 because the requirement set is absent.

- [ ] **Step 2: Extend Task 10 outputs**

Add:

```markdown
- 每个控制的 percentile rank、applicability-domain 状态和不确定性；
- 每个控制与无学习、传统模型和 ESM 基线的并列结果；
- `mechanism_lineage_id` 和独立验证单元状态；
- 所有未恢复、不可映射和证据不足的控制，不得只展示成功案例。
```

- [ ] **Step 3: Add Task 10 future RED specifications**

Insert before Task 10 acceptance:

```markdown
### 未来 RED 规格

1. `test_known_control_cannot_be_training_and_recovery`：控制同时进入训练并被报告为独立恢复时失败；
2. `test_duplicate_mechanism_lineage_is_not_independent`：同一机制谱系被计为多个独立验证单元时失败；
3. `test_control_rank_is_not_used_for_model_selection`：控制排名参与模型或超参数选择时失败；
4. `test_predictive_claim_requires_gate2_go`：Gate 2 为 `STOP` 时声称通用预测价值时失败；
5. `test_failed_and_unmappable_controls_are_reported`：失败或不可映射控制被遗漏时失败。

测试只能使用已登记真实控制记录或纯软件策略文本，不得创建虚构位点。
```

- [ ] **Step 4: Strengthen Task 13 section 13.3**

Retain rules 1–9 and use these final rules:

```markdown
10. 湿实验前冻结 `frozen_analysis_plan`，声明终点、命中定义、排除标准、匹配变量、统计比较和失败处理；
11. `blind_id` 是实验人员可见的唯一候选标识，排名和 Tier 保持隐藏；
12. 高排名候选、matched controls 和已知阳性使用一致检测与排除规则；
13. 系统主张依赖多个独立新位点或蛋白的整体富集，不得只挑一个成功案例；
14. 深入机制对象可从真实命中中选择，但不得回写或重排 `candidate_release_v1`；
15. 数据返回后新增 `prospective_validation_v1`，不覆盖原发布包。
```

- [ ] **Step 5: Add Task 13 future RED specifications**

```markdown
### 未来 RED 规格

1. `test_nature_readiness_requires_prospective_blind_validation`：无冻结候选、匹配对照和盲法结果时，不得标记 Nature-family ready；
2. `test_candidate_release_cannot_be_reranked_after_wetlab`：湿实验后修改顺序或覆盖原发布包时失败；
3. `test_unsuccessful_candidates_are_retained`：失败候选被删除时失败；
4. `test_prediction_is_not_causal_mechanism`：模型排名、位点修饰和因果表型混为同一证据层级时失败。
```

- [ ] **Step 6: Run GREEN and commit**

```powershell
rg -n -e "test_known_control_cannot_be_training_and_recovery" -e "test_duplicate_mechanism_lineage_is_not_independent" -e "test_predictive_claim_requires_gate2_go" -e "frozen_analysis_plan" -e "test_nature_readiness_requires_prospective_blind_validation" -e "test_candidate_release_cannot_be_reranked_after_wetlab" docs/PlantPersulf_Code_TDD_Codex.md
git diff --check
git status --short
```

Expected: all identifiers present, no whitespace errors, only the master TDD modified.

```bash
git add -- docs/PlantPersulf_Code_TDD_Codex.md
git commit -m "docs: specify retrospective and prospective validation"
```

---

### Task 3: Add contact gates and Nature-family claim ladder

**Files:**
- Modify: `docs/PlantPersulf_Code_TDD_Codex.md:1474-1725`

**Interfaces:**
- Consumes: Gates 1–3 and Tasks 10/13.
- Produces: Gates 4A/4B/5, Level 0–4, report declarations, and red lines.

- [ ] **Step 1: Run RED**

```powershell
rg -n -e "Gate 4A：非正式联系" -e "Gate 4B：带结果正式汇报" -e "Gate 5：Nature-family 故事就绪度" -e "## Level 4" docs/PlantPersulf_Code_TDD_Codex.md
```

Expected: exit 1.

- [ ] **Step 2: Replace Gate 4**

```markdown
## Gate 4A：非正式联系

即使 Task 5 为 `STOP`，仍可沟通研究方向、科学问题、已审计数据缺口、能力互补、数据请求和盲法实验设想。不得展示不存在的模型分数或候选排名，不得把计划写成结果，不得声称预测价值、新机制或确定的期刊结果。

## Gate 4B：带结果正式汇报

带模型结果和候选证据卡正式汇报至少要求：

1. registry 和可训练 benchmark 通过审计；
2. 防同源、研究、物种、时间和已知机制泄漏的 split 冻结；
3. 完成无学习、传统模型和 ESM 基线；
4. Gate 2 状态和允许主张明确；
5. 所有控制均有训练外恢复结果，包括失败和不可映射项；
6. 有 5–20 个初步番茄候选证据卡，不得为凑数量降低门槛；
7. 报告适用域、不确定性、基线和限制；
8. 准备好数据、实验能力和盲法验证请求清单。

若 Gate 2 为 `STOP`，只能展示证据审计和 `evidence_integration`，不得称为通用 predictor。
```

- [ ] **Step 3: Add Gate 5**

```markdown
## Gate 5：Nature-family 故事就绪度

本 Gate 只评价联合研究的故事证据结构，不预测或保证期刊结果。

只有同时满足以下条件，才允许内部标记 `nature_family_story_ready`：

1. 主张超越另一个单蛋白—单个位点机制；
2. Gate 2 已 `GO`；
3. 候选、matched controls、终点和统计方案在湿实验前冻结；
4. 完成前瞻性盲法验证并保留全部结果；
5. 高排名组显示多个独立新位点或蛋白的整体富集；
6. 至少一至两个中心命中完成深入生化或遗传验证；
7. 提出可检验的 PTM 串扰或其他一般原则；
8. 与主张相称地覆盖条件、材料、品种或物种；
9. 连接生长、产量、成熟、采后品质或更广泛植物问题；
10. 数据、模型、协议、来源链和失败结果在协议允许范围内可复现。

缺失核心条件时，降级为证据审计、benchmark/候选组织、回顾性预测、无一般机制主张的前瞻富集，或单候选机制研究。模型、随机交叉验证、单个新位点或单次成功实验均不足以通过 Gate 5。
```

- [ ] **Step 4: Update section 11 and replace section 12 ladder**

Add to external report requirements:

```markdown
- 当前最高通过 Gate；
- 当前类型：`evidence_audit`、`evidence_integration`、`retrospective_prediction` 或 `prospective_validation`；
- 排除训练的已知控制；
- 缺失关键证据；
- `STOP` Gate 对外主张的精确降级文字。
```

Replace Level 1–3 with:

```markdown
## Level 0：数据与证据审计
> 建立可追踪数据资源，报告证据、方法差异、冲突和缺口。

## Level 1：benchmark 与候选组织
> 建立真实 PU benchmark、基线和候选组织框架，不声称通用预测。

## Level 2：跨研究预测与回顾性恢复
> 在防泄漏 held-out studies 上证明跨研究价值，并恢复排除训练的机制。

## Level 3：前瞻性盲法富集
> 冻结番茄候选在盲法实验中相对 matched controls 富集真实位点。

## Level 4：系统规律、深入机制与植物意义
> 多个独立新命中支持可推广原则，一至两个中心命中完成机制，并连接重要植物或农业性状。

只有 Level 4 且 Gate 5 为 `READY`，才允许内部按 Nature-family 故事组织材料；不保证期刊结果。
```

- [ ] **Step 5: Extend report template, red lines, and current state**

Add report fields:

```text
- 当前最高通过 Gate：
- 当前允许 claim level：
- 已知控制是否全部排除训练：是/否/不适用
- 是否存在机制谱系重复计数：否/是（若是，任务失败）
- 是否产生前瞻性盲法结果：否/是
```

Append red lines:

```markdown
16. 同一蛋白、位点或机制谱系不得重复计算为多个独立验证。
17. Gate 2 为 `STOP` 时不得声称通用预测价值。
18. 没有冻结候选和前瞻性盲法验证，不得标记 Nature-family ready。
19. 单个成功候选不得替代多个独立命中的整体富集证据。
20. 预测、位点修饰、生化功能和因果表型必须分层表述。
21. Nature-family story readiness 是内部闸门，不是期刊结果保证。
```

Add before section 15:

```markdown
## 当前项目状态

Phase A 证据补齐已完成：从两个已登记、SHA256 审计的来源解析出约 393 个坐标验证过的
位点级 persulfidation 证据——PXD006140 Dataset S3（作者自定义 Sulfide/CN-Biotin-Sulfide
修饰，320 个 site_ms）与 PXD024061 MaxQuant Sulfide(C)/CianoBiotin(C) 位点表（73 个
class-I）。fail-closed 的 readiness v2 gate 据此判定 **Gate 1（数据可用性）= GO**，并已构建
PU benchmark v1、防泄漏 split、序列/ESM-2 特征与 motif 基线。

关键限制：两个研究同属 Seville（Romero/Gotor）实验室、同一 tag-switch 化学、同一物种
（拟南芥）。因此当前 GO **仅为数据可用性层面**，**不等于** Gate 2 的跨研究预测价值，也
**不满足** Gate 5。在 Gate 2 以防泄漏 held-out study 证明预测价值之前，对外只能作
`evidence_integration`，不得声称通用 predictor，不得标记 Nature-family ready。该状态已超过
Gate 4A（非正式联系），但**尚未**满足 Gate 4B——因为 Gate 4B 要求"Gate 2 状态与允许主张
明确"及"所有控制均有训练外恢复结果"，二者都还未完成。真正独立化学/实验室/物种的位点级
persulfidome 仍是对张华团队的优先数据请求。
```

- [ ] **Step 6: Run semantic and Markdown checks**

```powershell
rg -n -e "Gate 4A：非正式联系" -e "Gate 4B：带结果正式汇报" -e "Gate 5：Nature-family 故事就绪度" -e "nature_family_story_ready" -e "## Level 0：" -e "## Level 4：" -e "Gate 1（数据可用性）= GO" docs/PlantPersulf_Code_TDD_Codex.md
$docPath = 'docs/PlantPersulf_Code_TDD_Codex.md'
$fencePattern = '^' + (([char]96).ToString() * 3)
$fenceCount = (Select-String -Path $docPath -Pattern $fencePattern -AllMatches).Matches.Count
if (($fenceCount % 2) -ne 0) { throw "unbalanced fenced code blocks: $fenceCount" }
if (Select-String -Path $docPath -Pattern '^## Gate 4：与张华老师谈合作$') { throw 'obsolete Gate 4 remains' }
if (Select-String -Path $docPath -Pattern '^### Level [123]$') { throw 'obsolete ladder remains' }
```

Expected: all markers present, balanced fences, and no obsolete gate or ladder.

- [ ] **Step 7: Run repository verification**

```powershell
git diff --check
git status --short
rg -n "^## Task (0|1|2|3|4|5|6|7|8|9|10|11|12|13)：" docs/PlantPersulf_Code_TDD_Codex.md
pytest tests/unit -q
pytest tests/scientific -q
pytest tests/release -q
ruff check .
mypy src/plantpersulf scripts
```

Expected: only the master TDD is modified; Task 0–13 remain ordered; every verification command exits 0. No acquisition, analysis, training, or candidate command runs.

- [ ] **Step 8: Commit**

```bash
git add -- docs/PlantPersulf_Code_TDD_Codex.md
git commit -m "docs: add collaboration and Nature-family evidence gates"
```

---

## Final verification report

Report the three master-TDD commits, RED/GREEN checks, test and static-analysis results, unchanged Task 0–13 order, the current Task 5 `GO` state at Gate 1 (data availability) with its same-lab/same-chemistry limitation, documentation-only scope of these doc edits, remaining blockers (Gate 2 predictive value not yet demonstrated), and the highest currently permitted contact gate (above Gate 4A, not yet Gate 4B).
