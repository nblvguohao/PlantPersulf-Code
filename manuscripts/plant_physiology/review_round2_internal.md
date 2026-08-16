# 内部评审包 — pp_manuscript.tex(v3, 双支柱版, 2026-08-14)

## Review setup

- **Input scope**: full LaTeX manuscript after the two-pillar restructure (new title/abstract; R1 cross-kingdom conservation; R2 structural-context model-confidence control; renumbered R3–R7 resource strand; Tables 1–2 new; Figures 1–2 new; bibliography 27 entries). Round-1 fix list (7 items) verified applied.
- **Assessment boundary**: blind-cohort data do not exist yet (Gate 3); R8–R10, Table 7, Figure 8 remain preregistered placeholders and are assessed as such. The conservation and structural-context results are final, computed from registered inputs (`conservation_v2.json`, `structural_context_v2.json`).
- **Shared manuscript claim summary**: (pillar 1) persulfidation targeting is conserved at the ortholog-subfamily level beyond chance in the deeply sampled species, with tomato reported as underpowered-not-negative via an explicit detectability floor; an apparent structural signature of persulfidated cysteines is shown to be a model-confidence artefact and is corrected in-paper; (pillar 2) a frozen, hash-registered tomato candidate release with pre-blind calibration and a co-signed preregistered blind cohort.
- **Visible evidence base**: 6 pairwise Fisher tests with Bonferroni+BH and detectability floors; 4-species spectrum (universe 2,012) with Poisson-binomial null and 10,000-replicate permutation tests; structural context across 4 species × 3 metrics, pLDDT control and pLDDT≥70-restricted retest with per-species retained counts; the pillar-2 artifacts as in round 1.
- **Missing materials affecting confidence**: blind-cohort outcome; OSF preregistration DOI; Zenodo archive DOI; ORCIDs; Supplemental Table S1 (per-species panel composition, referenced in R3).

---

## Reviewer 1(technical soundness 权重最高)

- **Overall assessment**: The conservation layer is statistically careful in the ways that matter: the shared-universe restriction, the one-sided enrichment direction held consistently, both corrections reported side by side (and they are allowed to disagree on Arabidopsis–*Magnaporthe*), and — the part I have not seen done before — a detectability floor attached to every non-significant cell, which converts "not significant" into "could not have detected less than a 3.1–4.0× enrichment". That is the correct way to report an underpowered null. The spectrum layer uses a marginal-preserving permutation null rather than a naive binomial, which is right. The structural-context self-correction is executed properly: the confounder (pLDDT) is measured on the same residues, the direction match across four species is reported, and the claim is fully retracted rather than hedged.
- **Major strengths**:
  1. The detectability floor is a genuine methodological contribution to how sparse-omics conservation results should be reported; Table 1 is the template.
  2. The pLDDT confound chain (naive signal → confidence difference → confident-restricted collapse, with the tomato wt_only stratum identified as the driver) is a complete, self-contained correction — the strongest part of the paper for a methods reader.
  3. Cross-layer consistency checks pass: abstract numbers match Table 1/2; R2 numbers match `structural_context_v2.json`; permutation seeds and replicate counts match the scripts.
- **Major concerns**:
  1. **Confident-only null needs its denominators.** R2 reports that restricting to pLDDT ≥ 70 removes every SASA signal, but does not say how many sites survive the restriction. From the underlying output: Arabidopsis 311/384 (81%), rice 875/922 (95%), tomato 52/86 (61%), *Magnaporthe* 1,350/1,387 (97%). For three species the null is well-powered and the point stands; for tomato the null rests on 52 positives and is correspondingly weak. As written, a reader can misread the collapse as sample exhaustion, or conversely over-trust the tomato null. One sentence with the retained counts fixes both misreadings.
  2. **The spectrum tests are uncorrected across cells.** k≥2, k≥3, k≥4 are nested and therefore dependent; the permutation p = 1e-4 for k≥3 would survive a 3-cell Bonferroni, so the conclusion is safe, but the manuscript never says the spectrum layer is secondary to the pairwise layer. One clause in Methods would pre-empt the question.
  3. **The +0.140 in the abstract is a difference, not an AP.** "on an independent, one-shot homology-cluster frozen test (+0.140 pooled average precision, 95% CI 0.118–0.170)" reads as if the model scored AP 0.140; the number is the cluster-bootstrap difference against the strongest baseline (whose own AP is ~0.0015–0.0023 per species). The word "difference" must appear, or a careful reader accuses the abstract of inflating the strict-track result by two orders of magnitude.
- **Technical failings**: none beyond the above; the one-sided Fisher direction and the hypergeometric from-scratch implementation are consistent with the code and its tests.
- **Recommendation posture**: technically sound; the concerns are reporting completeness, not correctness.

## Reviewer 2(originality + significance 权重最高)

- **Overall assessment**: Pillar 1 is the first direct cross-kingdom test of persulfidation targeting and is, on its own, a publishable unit: the answer is "yes at family level among deeply sampled species, currently untestable for tomato, and here is exactly what tomato would have needed". The structural-context correction elevates it further — most groups would have published the naive buried/exposed story. The two-pillar architecture (associational biology + preregistered resource) is coherent and the Introduction now earns it.
- **Major strengths**: first cross-kingdom conservation test for this PTM; the in-paper self-correction with a reusable caution; the "state the test before seeing its answer" through-line that now genuinely unifies both halves.
- **Major concerns**:
  1. **The functional interpretation overreaches the family-level evidence.** The Discussion says the conserved set "points toward NAD(P)H-binding active-site cysteines in central metabolism". Two problems. (a) The data say a *family* is persulfidated, not *which* cysteine — nothing here localizes the modification to an active-site or cofactor-proximal residue. (b) Even at the enzyme-class level the label is strained: fructose-1,6-bisphosphatase is the canonical thioredoxin-regulated Calvin-cycle enzyme (a disulfide target, not an NAD(P)H enzyme), and glycine dehydrogenase/alanine-glyoxylate aminotransferase are PLP enzymes. What the list actually supports is "central carbon and energy metabolism, dominated by NAD(P)-dependent oxidoreductases and photosynthetic/photorespiratory enzymes" — which is still a good sentence.
  2. **The eight-family count is wrong in one cell.** The peptidyl-prolyl *cis-trans* isomerase D-related protein is a protein-folding catalyst, not a "core redox or central-carbon-metabolism enzyme". The correct count is six of eight, with two exceptions (the PPIase and RCC1). This appears in R1 ("Seven of these eight… and a peptidyl-prolyl cis-trans isomerase D-related protein; the eighth is RCC1"), in the abstract ("seven of them core redox or central-carbon-metabolism enzymes"), and in Figure 1c's encoding ("the one non-redox exception (RCC1)"). A plant biochemistry referee will catch this on first read; it is exactly the kind of slip that damages an otherwise careful paper.
- **Recommendation posture**: enthusiastic about the architecture; the two interpretation slips are small in word count but large in reviewer-trust cost.

## Reviewer 3(跨学科兴趣 + 可读性权重最高)

- **Overall assessment**: The restructure reads much better than v2: the paper now opens with a biological question rather than a governance argument, and the "two things about biology, three about method" Discussion opening is an effective map. The provenance vocabulary is now introduced where needed.
- **Major concerns (readability/format)**:
  1. **Abstract is ~318 words; Plant Physiology's limit is 250.** Round 1 already enforced this limit; the two-pillar abstract regressed. The trim should come from the dataset enumeration and the structural-correction clause, not from the numbers.
  2. **Internal count inconsistency**: the Introduction says the ranker "leads six models"; Results and abstract say six models total (ranker + five alternatives). One word fix.
  3. **Figure 1 legend (b) promises k = 2, 3, 4 bars; the figure shows k = 2, 3 as bars with k = 4 as an axis note.** Say so in the legend.
  4. "as of our audit date" (Introduction, tool landscape) — give the date (2026-08); it will age better.
- **Recommendation posture**: positive; all fixes are one-line except the abstract trim.

---

## Cross-review synthesis

- **Consensus strengths**: the detectability floor as a reporting device (R1); first cross-kingdom test + in-paper self-correction (R2); much improved narrative architecture (R3).
- **Consensus risks**: (1) eight-family composition misclassified (six, not seven; two exceptions) — R2/R3; (2) "NAD(P)H-binding active-site cysteines" overreaches family-level evidence — R2; (3) abstract over word limit and the +0.140 reads as absolute — R1/R3; (4) confident-only null lacks retained counts — R1.
- **Emphasis differences**: R1 wants denominators and one clause on nested spectrum cells; R2 wants the interpretation clipped to what family identity supports; R3 wants format compliance.
- **Most important issues to resolve**: (1) and (2) are scientific-accuracy fixes; (3)–(4) are precision/completeness fixes; all are small.

## Risk / unsupported claims

- "points toward NAD(P)H-binding active-site cysteines" — **not supported** at family resolution; must be softened (fix 2).
- "seven of them core redox or central-carbon-metabolism enzymes" — **incorrect count**; six, with PPIase-D and RCC1 as the two exceptions (fix 1).
- "+0.140 pooled average precision" in the abstract — **ambiguous to the point of misleading**; it is a difference over the strongest baseline (fix 4).
- Everything else checked: pairwise/spectrum numbers vs `conservation_v2.json`, structural numbers vs `structural_context_v2.json`, retained counts, permutation seeds — all consistent.

---

## 修改清单(据此修订 pp_manuscript.tex 与 fig1)

1. **八家族分类纠正(三处 + 图)**: R1 正文 "Seven of these eight…" → "Six of these eight…",PPIase D-related 移入例外;摘要 "seven of them" → "six of them";Discussion 列举保持五个酶不变;`fig1_conservation.py` 灰点编码加入 PPIase D-related、注记改为 "the two non-redox exceptions",重新生成 fig1 全部格式;Fig. 1 图注 (c) 同步说明。
2. **Discussion 功能解读软化**: "points toward NAD(P)H-binding active-site cysteines in central metabolism" → 家族层面可支撑的表述(中心碳/能量代谢为主:NAD(P) 依赖脱氢酶、复合体 I NADH 模块、Calvin 循环/光呼吸酶;不断言活性位点定位)。
3. **摘要压缩至 ≤250 词**(PP 上限;删除数据集枚举细节,保留全部数字)。
4. **摘要 +0.140 精确化**: 明确为 "pooled average-precision difference over the strongest baseline"。
5. **Intro (iii) 计数一致**: "leads six models" → 与 R4/摘要一致的 "leads five alternatives"。
6. **R2 补保留计数**: pLDDT≥70 过滤后各物种保留的 persulfidated Cys 数/比例 (311/384, 875/922, 52/86, 1350/1387),并注明番茄 null 权重相应较弱。
7. **Methods 一句话**: spectrum 各细胞嵌套、作为 pairwise 层的次级证据,其显著性在 3 细胞 Bonferroni 下仍成立。
8. **Fig. 1 图注 (b)**: 说明 k=4 以轴下注记给出(无家族满足)。
9. **Intro 工具审计句加日期**: "as of our audit date" → "as of our audit date (2026-08)"。

全部修改完成后:重跑两道 wording gate → 重新生成 fig1 → 重编译 PDF。
