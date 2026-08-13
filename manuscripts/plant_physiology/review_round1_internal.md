# 内部评审包 — pp_manuscript.tex(v2, 2026-08-14)

## Review setup

- **Input scope**: full LaTeX manuscript (title/abstract/intro/R1–R8/discussion/methods/5 tables/5 rendered figures + 1 placeholder), frozen release artifacts, Gate 1 statistical report.
- **Assessment boundary**: blind-cohort data do not exist yet (Gate 3); R6/R7, Table 5, Figure 6 are preregistered placeholders. This review assesses the pre-blind manuscript as it stands.
- **Shared manuscript claim summary**: a provenance-tracked four-species persulfidation site dataset; a within-dataset benchmark led by a gated-fusion PU ranker; a frozen, hash-registered tomato candidate release; an honest pre-blind calibration; and a co-signed, preregistered blind-cohort design whose outcome will be reported regardless of direction.
- **Visible evidence base**: benchmark summary (60 runs), Gate 1 strict frozen-test report (byte-identical across three scoring runs, audit-passed manifest pending at review time), candidate table (179,736 rows), matched controls (966 pairs), mock-blind calibration (20 clean kiae271 + 7 Arabidopsis controls), 144-cell power grid.
- **Missing materials affecting confidence**: blind-cohort outcome; OSF preregistration DOI; public archive DOI; Supplemental Table S1 (per-species composition).

---

## Reviewer 1(technical soundness 权重最高)

- **Overall assessment**: The engineering and statistical hygiene is far above field norms: frozen artifacts, SHA256 registration, a one-shot unlock with run-manifest audit, preregistered claim gate, and honest power tabulation. The technical case for the *framework* is established. The technical case for the *candidate table's usefulness in tomato* is deliberately not established yet, and the manuscript says so.
- **Who would be interested**: plant redox/proteomics groups tired of over-claimed PTM tools; methods reviewers who care about evaluation hygiene.
- **Major strengths**:
  1. The strict frozen test is genuinely leak-proof (homology-cluster split, one-shot unlock, audit) and its negative tomato result is reported, not hidden.
  2. The mock-blind calibration is reported with per-site data and used to size the design rather than to sell the model.
  3. The claim gate refused a stronger claim class on a −0.0004 tomato delta; the machinery demonstrably has teeth.
- **Major concerns**:
  1. **Flat-top fragility**: the top 200 candidates span 0.295–0.305 of the score range; rank order within the head is nearly arbitrary. Limitation (iii) acknowledges this, but the primary cohort *is* the flat head. The manuscript should state explicitly why K=200 enrichment (rather than rank position) is the right readout under a flat head — it gestures at this but should say it in R5, not only in Discussion.
  2. **Tiny positive counts in the strict test for tomato (20 positives)**: AP 0.00032 vs base rate 0.00042 on 20 positives is noise-dominated; the "below base rate" framing is fine but the manuscript should give the positive counts alongside the APs so readers can calibrate.
  3. **PU setting with ~0.6% positives**: the ranking objective treats unlabeled as negative; with kiae271 covering only part of the tomato persulfidome, some "controls" may be true sites. The detectability confound (limitation v) is acknowledged; the PU label-noise direction (controls contaminated by unassayed true sites biases enrichment toward null) should be stated — it strengthens rather than weakens the design's honesty.
- **Technical failings to address before the case is established**:
  1. The sentence "No new persulfidation-specific ranking tools appeared in the 2025 literature" contradicts the manuscript's own citation of Sul-BertGRU (Bioinformatics, 2025). Must be corrected.
  2. "In prior work on Arabidopsis we showed..." has no citation or archive pointer. Either cite a public artifact or reword as an earlier phase of this project with an archive reference.
  3. The wet-laboratory partner's identity relative to kiae271 is not disclosed in the manuscript. Since the blind-cohort laboratory is the laboratory that produced the tomato reference persulfidome, the independence/blinding handling (exclusion of kiae271-protein sites from the primary analysis) should be stated explicitly in the design section; silence here reads as concealment.
- **Assessment against Nature-style criteria**: technically sound within its self-declared bounds; the within-dataset claims are properly fenced.
- **Recommendation posture**: minor-to-moderate revision for a field journal; the failings are presentational/honesty-of-framing, not scientific.

## Reviewer 2(originality + significance 权重最高)

- **Overall assessment**: The individual components (PU ranking, candidate release, power analysis) are not new; the *combination* — freeze + hash registry + pre-blind calibration + co-signed preregistration + committed reporting — is genuinely new for the cysteine-PTM field and rare in plant proteomics generally. Significance is field-local but real: it raises the evidentiary standard.
- **Who would be interested**: PTM tool developers (as a template), plant signaling groups (as a candidate source), journal editors (as a preregistration exemplar).
- **Major strengths**: the preregistered decision matrix with pre-written outcome variants; the negative-results-first archiving policy; the two-paper strategy implicit in the structure.
- **Major concerns**:
  1. The manuscript cannot escape the question "what do we know now that we did not know before?" — the honest answer is "a dataset, a frozen candidate list, calibrated expectations, and a binding protocol", which is a resource contribution. The abstract's final sentence should own this framing rather than leaning on the future blind test.
  2. Originality vs prior preregistered/resource papers in proteomics is asserted only by the Nosek/Chambers citations; one sentence situating this among registered-resource precedents in other omics fields would strengthen the novelty case.
- **Technical failings**: none beyond Reviewer 1's list.
- **Recommendation posture**: supportive if the framing is owned as a resource + protocol paper.

## Reviewer 3(跨学科兴趣 + 可读性权重最高)

- **Overall assessment**: The story is compelling and unusual; a plant-science reader can follow the logic. The provenance/hash vocabulary (SHA256, unlock records, run manifests) is foreign to most PP readers and needs one orienting sentence each time it appears.
- **Major concerns (readability)**:
  1. The abstract is at/over the 250-word limit and tries to carry the whole pipeline; trimming the calibration sentence would buy room.
  2. Figure 1 is the right anchor; consider moving the "binding rule" text into the caption so the figure breathes.
  3. "mock-blind diagnostic" is jargon; define it at first use (it is defined, but late).
- **Recommendation posture**: positive; readability fixes are minor.

---

## Cross-review synthesis

- **Consensus strengths**: evaluation hygiene; honest negative/calibration reporting; preregistered machinery with demonstrated teeth; complete artifact archiving.
- **Consensus technical risks**: (1) the 2025-tools sentence contradicts the Sul-BertGRU citation; (2) uncited "prior work" claim; (3) wet-lab/kiae271 relationship undisclosed; (4) flat-top cohort rationale belongs in R5; (5) strict-test tomato positive count (n=20) should accompany APs.
- **Emphasis differences**: R1 wants quantitative context (positive counts, label-noise direction); R2 wants the resource framing owned in the abstract; R3 wants jargon orientation and a leaner abstract.
- **Broad-interest readout**: field-internal but exemplary; the reusable message is the two-step "calibrate, then size the blind design".
- **Most important issues to resolve**: items (1)–(3) are factual/consistency fixes; item (4)–(5) are one-sentence additions; abstract trim.

## Risk / unsupported claims

- "establishes a reproducible, preregistered pipeline" — supported as a resource claim; not yet supported as a *validated* prioritization tool (by design; Gate 3 pending).
- Strict-test tomato AP "below base rate" — supported numerically but on 20 positives; must carry the count.
- Any implication of wet-lab independence from kiae271 — not supported; must be disclosed and handled by the exclusion design.

---

## 修改清单(据此修订 pp_manuscript.tex)

1. Introduction: "No new ... appeared in the 2025 literature" → 纠正为与 Sul-BertGRU (2025) 一致的表述。
2. Introduction: "In prior work on Arabidopsis we showed..." → 改为指向项目档案的表述。
3. R5 Cohort: 明确披露湿实验室伙伴即 kiae271 数据来源实验室,并指出排除设计如何处理。
4. R5: 一句话说明平坦头部分数下 K=200 富集(而非名次)是正确读出。
5. R2 strict-test 段: 标注番茄阳性数 n=20。
6. Discussion limitations: 增加 PU 标签噪声方向(对照臂可能含未测定真位点 → 偏置趋向零假设)。
7. Abstract: 精简到 ≤250 词,拥有 resource+protocol 定位。
