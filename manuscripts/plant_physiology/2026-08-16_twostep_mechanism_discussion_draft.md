# R4 已知控制章节 —— 两步持硫化机制重释（Discussion 草稿）

> **用途**：供用户在下一轮稿件修订（R4 已知控制/机制章节）自行取舍插入。**不自动写入手稿
> tex**（稿件有 wording governance，作者签名处须最终定夺）。全部措辞为可提交级英语，但属于
> 机制重释（narrative）而非已验证结论——插入时请保留文中的 hedge（"in principle"、
> "consistent with"、"we do not claim"），不要升级为因果结论。
>
> **引用**：Corpas et al. 综述处于审稿中（COPLBI-D-26-00068），正式引文尚未定稿；占位
> `\citep{corpas2026inreview}`，投稿前需替换为其正式卷期或改标 preprint/通讯。
>
> 关联记录：`docs/superpowers/plans/2026-08-13-model-improvement-register.md` §2.3；
> 核验审计 `results/diagnostics/oxiptm_sites_v1.json`；先验数据
> `results/diagnostics/{copeptide_structure_separation_v1,sg_chemical_occupancy_v1}.json`。

---

## Draft paragraph 1（主段：为什么共肽分离在原理上受任务约束）

The matched-detection contrast explored in this paper — the co-peptide
diagnostic, in which a modified cysteine is compared with unmodified
cysteines that share its peptide, spectrum, digestion, and enrichment —
returned a null result across both the sequence model (3 species, 17
protein–peptide instances; observed separation below the permutation
expectation) and every structure feature we computed (7 features across
14 groups; all $p > 0.05$, directions inconsistent). A mechanistic
constraint in the persulfidation literature suggests why this null is
expected rather than merely uninformative. Persulfidation of a cysteine
thiol is contingent on prior activation: the thiol must first be
oxidized to sulfenic acid ($-$SOH), a disulfide, or an S-nitrosated
form before it can react with H$_2$S (Corpas et al., in review). Under
this two-step view, a co-peptide unmodified cysteine is a latent mixture
of two distinct failures — a site that was never oxidized (gate 1) and a
site that was oxidized but not persulfidated (gate 2) — and a site of
the second kind is, from the perspective of any oxidation-gate feature,
structurally indistinguishable from a modified site. A single-stage
binary label therefore compresses two physically distinct processes, and
no single feature family is required to separate positives from such a
mixed negative class. We offer this as a reinterpretation of the
co-peptide null rather than a demonstrated mechanism: our data cannot
resolve oxidation state, so the two-stage account is not directly
testable on this axis. Its value is that it bounds the claim we do make —
that within the current feature resolution the co-peptide contrast is a
region where neither sequence nor structure should be asserted to
discriminate — and it motivates the next-generation formulation below.

## Draft paragraph 2（机制重释 + 下一代方向）

The two-step view also reconciles two otherwise paradoxical observations.
First, our most robust structure feature — packing density, on which
known functional sites are significantly more buried than same-protein
background (permutation $p = 0.018$, leave-one-out) — is coherent with
persulfidation reactivity only under the oxidation-gate reading: a
coordinated thiolate is a low-$pK_\mathrm{a}$, strongly nucleophilic
cysteine that is simultaneously the preferred oxidation and H$_2$S
target, whereas the naive assumption that a reactive cysteine must be
surface-exposed would predict the opposite sign. Second, the single
co-peptide case that resolved at the structure level — the BRG3 RING
cluster, where the two modified cysteines are respectively the most
exposed and the most Zn-coordinated of the cluster while the unmodified
cysteine is buried and uncoordinated — is exactly the disjunction the
single-feature contrast could not see: either exposure (gate 1 for
Cys206) or coordination-dependent activation (gate 1 for Cys212) suffices,
and neither alone is necessary. These are case-level consistencies, not
aggregate significance, and we register them as such. They nonetheless
point to the concrete next-generation directions: a two-stage model
$P(\mathrm{persulfidation}) = P(\mathrm{oxidizable}\mid\mathrm{structure})
\times P(\mathrm{persulfidated}\mid\mathrm{oxidized})$, in which the
first stage borrows the abundant, chemically conserved data on cysteine
oxidation across species and the second stage consumes the scarce plant
persulfidation labels; and a subcellular-compartment conditioning term,
since the thiolate fraction of a cysteine of fixed $pK_\mathrm{a}$
changes by roughly an order of magnitude over the pH range spanned by the
cytosol (${\sim}7.2$) and the illuminated chloroplast stroma (${\sim}8.0$).
We emphasize that both are preregistration-bound changes to be evaluated
only after the present blind cohort, and that the oxiPTM-competition
benchmark this framing suggests — discriminating persulfidation from
S-nitrosylation at reactive cysteines, for which "the other modification"
is a more reliable negative than a putatively unmodified site — is now assemblable
from published tables, after per-site primary-source resolution. The
reference review's table spans species, and its six "Arabidopsis-looking"
S-nitrosylation symbols are in fact tomato genes; with each site mapped to
its primary accession, six of ten S-nitrosylation sites and eleven of
thirteen persulfidation sites survive coordinate verification (the
remaining four S-nitrosylation sites fall into reference-proteome coverage
gaps, and the two tomato persulfidation sites whose accessions were never
resolved are documented as unmappable). On this verified subset, the frozen
model showed no persulfidation-specific ranking advantage: known
persulfidation sites did not rank above known S-nitrosylation sites within
their own proteins — if anything the direction ran the other way
(per-site mean within-protein rank 3.8 vs 2.0, with three of the six
S-nitrosylation sites ranked first in their proteins; permutation
$p = 0.976$ for the predicted direction, $n = 8$ vs $6$ proteins). This is
a small, descriptive result, but its direction is the one the two-step view
predicts: a model whose signal is largely generic oxidation-gate
(reactive-Cys) rather than persulfidation-specific would not be expected to
separate the two oxiPTM classes. We register it as such, and as a caution
against treating an unmodified cysteine as a clean negative for any oxiPTM
whose gate-1 chemistry overlaps the modification of interest.

---

## 位置建议

- 草稿段 1 可放 Discussion 已知控制/共肽诊断小节（与现有"Pre-blind calibration"后的
  结构讨论区相邻），或作为共肽诊断结果的机制收尾句。
- 草稿段 2 的前半（埋藏信号自洽 + BRG3 OR 逻辑）可并入"结构承载区分"叙事；后半
  （两步架构 + 区室 + oxiPTM 基准 + SNO 核验失败）放 Discussion 未来方向。
- 如需更保守的短版，仅保留段 1 + 段 2 的 "These are case-level consistencies..." 到
  "...rather than a demonstrated mechanism" 部分。
