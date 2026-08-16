# 投稿总清单 — Plant Physiology

主稿：`../pp_manuscript.tex`（两支柱手稿）。本文件是投稿动作的唯一入口：
OSF 预注册、Zenodo 归档、格式核对、Gate 3 触发流程，全部在此汇总。

## 1. 状态总表（2026-08-14）

| # | 事项 | 状态 | 属主 |
|---|---|---|---|
| 1 | 两支柱手稿（保守性 + 冻结番茄资源） | ✅ 完成并提交 git（810dfc0、3e6cdf2） | 我 |
| 2 | 8 幅图 Nature 风格 + fig3 复刻 | ✅ 完成并提交 git（3e6cdf2） | 我 |
| 3 | 补充表 S1（生成器 + 测试 + PDF/TSV） | ✅ 完成（已入库，results 侧物化） | 我 |
| 4 | 封面信（两支柱版，措辞门 []） | ✅ 完成（`../pp_cover_letter.md`） | 我 |
| 5 | PP 格式核对（摘要/字数/alt text/行号） | ✅ 2026-08-14 核对，见第 2 节 | 我 |
| 6 | OSF 预注册 DOI | ⏳ 按 `osf_preregistration/README.md` 操作 | **你** |
| 7 | Zenodo 归档 DOI（含 LICENSE 拍板） | ⏳ 按 `zenodo_archive/README.md` 操作 | **你** |
| 8 | ORCID（全部 7 位作者） | ⏳ | **你** |
| 9 | CRediT 角色确认（全部合作者） | ⏳ | **你** |
| 10 | 基金信息核对（`作者信息和基金.txt`） | ⏳ | **你** |
| 11 | Gate 3 盲法队列数据（湿实验室） | ⏳ 硬阻塞，见第 3 节 | 张华实验室 |

条目 6–10 完成后把 DOI/信息发给我，我填入稿件 5 处「to be inserted」位置并重建 PDF。
条目 11 是唯一的外部硬阻塞：盲数据到达前，本手稿不可提交（Gate 3 占位章节不能提前填写）。

## 2. Plant Physiology 要求核对结果（2026-08-14 核对现行 General Instructions）

| 要求 | 本稿状态 |
|---|---|
| 摘要 ≤ 250 词（不得含参考文献或缩写） | ✅ 恰好 250 词（2026-08-14 逐词复核）；无参考文献引用；缩写已全部消除（2026-08-14：pLDDT 展开为 predicted local distance difference test，CI 展开为 confidence interval） |
| 正文 ≤ 7,000 词（含摘要；图注、表格、参考文献不计） | ✅ 实测 ≈ 6,250 词（粗测，排除图注/表格/书目） |
| 图+表 typical 6–10 个 | ⚠️ 本稿 7 图 + 7 表 = 14 个；属「typically」软性建议而非硬限制，记录在案，若编辑提出再合并（候选：表 2–4 并入 S 表） |
| 参考文献 typical 30–50 条 | ⚠️ 本稿 27 条；同为软性建议，不虚构引用 |
| 行号 | 说明书未作要求；已启用 `lineno`（审稿便利，提交定稿时可移除） |
| 每图 Alt text（图注正下方，以 "Alt text:" 开头） | ✅ 已为图 1–7 加入（2026-08-14）；图 8 为 Gate 3 占位，随图补 |
| 图以单独文件提交（位图 ≥350 dpi，矢量不限；多 panel 合为一图，panel 字母标左上角） | 📋 投稿系统上传时按图 1–7 单独文件提交 `figures/*.pdf`（矢量）；panel 已在图内字母标注 |
| 审稿期图注须同时出现在正文文本文件与图文件中 | 📋 正文已含图注（Figure legends 节）；图文件本身无文字图注——如需，用 `figures/figures_with_legends.pdf` 合并件提交 |
| 补充材料随主稿同时提交且须在正文引用 | ✅ S1 已在正文引用（R2「per-species positive counts in Supplemental Table S1」） |
| 质谱数据：须给检索/统计参数与肽段细节 | ✅ Methods 已含（四数据集、坐标校验、SHA 锁定参考蛋白组）；原始数据已在 ProteomeXchange（Data availability） |
| 预印本 | 允许，须标注「under review」——本稿未挂预印本 |

## 3. Gate 3 触发流程（盲数据到达后）

1. **验包**：两人独立核对 bundle、fit-manifest、候选表、对照表、推理脚本哈希（`scripts/audit_candidate_release.py`）。
2. **先锁分后开标**：盲队列打分输出写出并哈希锁定 → 释放标签 → 按 `analysis_plan.json` 原样运行主分析。
3. **填入稿件**（只动 Gate 3 占位区，占位注释以上一律不动）：
   - R6「Blind-cohort primary endpoint」：2×2 表（Table 5）、Fisher p、OR 及精确 95% CI、逐端点数字、过程对照通过率；
   - R7 灵敏度与次级分析（预注册清单已在占位中列全）；
   - Table 5、Figure 8（`fig8_blind_cohort.pdf`，替换占位注释并取消 float 注释）；
   - 摘要结局句：在摘要最后一句之前插入下述 A 或 B 变体（逐字）；
   - 更新图 8 的图注与 Alt text。
4. **措辞门**：`verify_predictive_claims -> []`、`find_forbidden_external_claims -> []` 全量重跑。
5. **禁则**（解盲后永久）：改模型/特征/Top-K/阈值；给失败端点加样；重排候选表。任何偏离 = 有日期、带哈希、共同签署的修订件。

### 摘要结局句变体（2026-08-14 从稿件源注释移入，逐字使用）

**A（成功）：** "The preregistered blind cohort met its success criteria (odds ratio X.XX, exact 95% CI [...], one-sided p = ..., N confirmed candidates of K = 200), demonstrating enrichment of confirmed persulfidation sites among frozen top-ranked candidates under fully blinded conditions."

**B（未达成功标准）：** "The preregistered blind cohort did not meet its success criteria (odds ratio X.XX, exact 95% CI [...], one-sided p = ...); we report the complete result, the pre-registered sensitivity analyses, and the resource artifacts, which remain usable for candidate organisation under the stated limitations."

插入位置：摘要末句「All artifacts, negative results, and exclusion ledgers are openly archived; the blind outcome will be reported regardless of direction.」之前。

## 4. ScholarOne 提交文件清单（全部就绪后）

1. 主稿 PDF（含行号审稿版；定稿版再去行号）
2. 图 1–7 单独 PDF（`figures/fig1_conservation.pdf` … `fig7_power_landscape.pdf`）；如编辑部要求图文件含图注，另附合并件
3. 补充材料：`../supplements/supplemental_table_s1.pdf`（+ `.tsv` + `.tex`）
4. 封面信（`../pp_cover_letter.md` 转 PDF）
5. OSF 预注册 DOI（R5、R7、Methods、Data availability 填入后）
6. Zenodo 归档 DOI（Data availability 两处填入后）
7. 全部作者 ORCID / CRediT / 基金信息（来自 `作者信息和基金.txt`）
