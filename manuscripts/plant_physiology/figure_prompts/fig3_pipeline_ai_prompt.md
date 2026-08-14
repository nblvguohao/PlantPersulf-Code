# Figure 3 — AI 绘制提示词(Plant Physiology 期刊风格)

> 用途:把下面的英文提示词粘贴到 AI 绘图工具(Midjourney / DALL·E / Gamma 等)。
> 生成结果建议作为**版面设计样稿**:AI 图像工具渲染大段文字极易产生错字、改数字,
> 而定稿图必须文字可编辑(PP 要求矢量可编辑文字)。
> 拿到满意的版式后,把样稿发回给我,我按同一版式用 matplotlib 复刻(
> 精确字符串、600 dpi、SVG/PDF/TIFF,继续纳入版本管理)。

---

## English prompt(直接粘贴)

```
Scientific journal figure for Plant Physiology. A clean editorial schematic,
flat vector style, white background, showing the six-stage study design of a
preregistered blind-validation project, laid out as a horizontal left-to-right
flow of six equal rounded-rectangle boxes connected by thin straight arrows,
with two annotation bands beneath.

Style: modern research-article figure, Plant Physiology journal convention —
Arial/Helvetica-like sans-serif typography, 174 mm double-column width,
restrained flat color palette (deep blue #0F4D92, grey #767676, dark grey
#4D4D4D), no gradients, no 3D, no drop shadows, no clip-art icons, no
photographs, generous white space, 1 pt strokes.

Layout: one row of 6 boxes; boxes 1, 2, 4, 6 are white with a blue outline;
boxes 3 ("Freeze") and 5 ("Preregistration") are filled deep blue with white
text to mark the two binding commitment points. Below the row, a full-width
light-grey band (grey #F2F2F2, rounded) carrying the binding rule; below it,
three small outlined cells in a row carrying the preregistered decision matrix
(condition on top in dark grey, outcome beneath in italic teal #42949E).

Box contents (render EXACTLY, including numbers and dates):
1. Dataset — 4 species, 5 studies; 389,609 sites; 2,334 positives; full
   provenance
2. Benchmark — 6 models, 10 seeds; within-dataset splits; literature-
   comparable claim class only
3. Freeze — weights + features + candidates + controls; SHA256-registered;
   co-signed 2026-08-13
4. Calibration — mock-blind diagnostic; 20 unseen tomato sites; 7 Arabidopsis
   controls; expectations set pre-blind
5. Preregistration — SAP + analysis plan; OSF-registered DOI; primary endpoint
   fixed; success criteria fixed
6. Blind cohort — 593 sites assayed blind; one-sided Fisher exact; reported
   regardless of outcome

Binding-rule band text (single line, centred): "Binding rule: no model,
feature, candidate-table, threshold, or control-set change after unblinding;
every artifact hash-verified; negative results and exclusions preserved and
archived."

Decision matrix (three cells, condition above, outcome below):
"enrichment met + mechanism" -> "upgraded submission target"
"enrichment met, no mechanism" -> "resource and validation report"
"enrichment not met" -> "negative result, first-class finding"

If exact text rendering is not possible, keep every number and date exactly as
written or replace a label with a clean placeholder bar — never alter or
invent any number, date, or scientific term. Output as a single wide image,
no extra panels, no caption.
```

---

## 内容保真清单(生成后逐项核对,任何一处不符即弃用该图)

| 项 | 正确值 |
|---|---|
| 站点总数 | 389,609 |
| 阳性数 | 2,334 |
| 共签日期 | 2026-08-13 |
| 番茄未训练位点 | 20 |
| 拟南芥对照 | 7 |
| 盲测队列规模 | 593 |
| 统计检验 | one-sided Fisher exact |
| 阶段数 | 6(顺序: Dataset → Benchmark → Freeze → Calibration → Preregistration → Blind cohort) |
| 决策矩阵 | 3 格(upgraded submission target / resource and validation report / negative result) |

## 提示词使用说明

1. **为什么给出两段提示词**:AI 绘图模型对长文本的还原不可靠,提示词里明确
   允许"用占位条代替无法还原的文字",是为了保护数字与日期不被篡改。
2. **定稿要求**:Plant Physiology 要求图中文字可编辑(矢量)。AI 生成的位图
   只能作为版式样稿;最终图必须由脚本/矢量工具排入精确文字。
3. **后续流程**:选定版式 → 发我样稿 → 我在 `fig3_design_pipeline.py` 中按
   同版式重绘(精确内容、600 dpi、SVG/PDF/TIFF、断言不变)→ 重新编译 PDF。
