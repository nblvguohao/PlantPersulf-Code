# P0-C pilot: non-Cys PTM confound-collapse generalization feasibility

**Date**: 2026-08-17
**Status**: literature/data reconnaissance only — no predictor run, no data downloaded, no code written.
**Gate**: plan requires evidence that the L0→L1 (within-protein) confound collapse is plausible outside cysteine PTMs, since this is what would push the paper from a plant/Cys-specific NC story to an NMI-relevant, field-wide methods claim (`nmi-nc-zippy-kernighan.md` §2, §5 Phase 0).

## Method

Literature search (PubMed via MCP tools, Consensus.app via WebSearch) for existing
documentation of PTM-predictor performance collapse under within-protein or
same-peptide held-out evaluation, plus a survey of public phosphorylation-site
datasets and runnable predictors that could support a fast pilot.

## Headline finding: the precedent already exists and is close to our argument

**Zuallaert et al. 2024 (bioRxiv), "A study on experimental bias in
post-translational modification predictors"** (PhosphoLingo tool; code+data at
https://github.com/jasperzuallaert/PhosphoLingo) shows that current dataset-compilation
and evaluation practice for phosphorylation (and other PTM) predictors creates an
**artefactual protease-specific training bias** that inflates reported accuracy, and
proposes restricting negative candidates to residues actually present in matched
detected peptides — i.e., essentially the same-peptide/co-peptide negative-sampling
fix our L2 level already implements for persulfidation. This is the closest existing
precedent to our claim, frames it via digestion/detectability bias rather than
within-protein/homology bias specifically, and is directly citable as prior art we
extend and generalize (must be distinguished, not just cited, in the paper).

Two further, weaker precedents:
- **Piovesan et al. 2020, PLoS Comput Biol**, PMID 32569263, [doi:10.1371/journal.pcbi.1007967](https://doi.org/10.1371/journal.pcbi.1007967) — benchmarked 7 hydroxylation-site predictors against newly collected independent sites; self-reported performance was "widely overestimated," no predictor beat random on new examples. Temporal/novel-protein holdout, not within-protein/same-peptide, but the standard citation for "self-reported PTM predictor performance is untrustworthy."
- **Esmaili et al. 2022, Genomics Proteomics Bioinformatics** — online phosphorylation predictors perform notably worse on proteins unseen at training time than on their original benchmarks (a generalization-gap result, not a within-protein-negative ablation).

**No existing paper runs the precise within-protein or same-tryptic-peptide negative
construction we propose for phospho/ubiquitination/acetylation.** The genericity claim
is plausible and has adjacent documented support, but is not yet demonstrated in this
exact form — which is exactly the paper's opening contribution, not a redundant replication.

## Data and tooling feasibility

| Resource | Fit | Friction |
|---|---|---|
| **PhosphoLingo's own bundled datasets** (Zuallaert et al., public GitHub) | Best fit — already peptide/spectrum-aware, pre-built for exactly this bias analysis | None significant; pretrained models + code public |
| PhosphoSitePlus | Largest curated resource, excellent multi-site-per-protein coverage | Bulk/programmatic download requires academic registration + terms of use, not fully open API |
| dbPTM (Academia Sinica) | Aggregates PhosphoSitePlus + UniProt + others | Flat-file downloads more open than PSP directly |
| Ochoa et al. 2020 kinase-regulation atlas (PMID 31819260) | Good multi-site coverage, functional scores, public supplementary tables | — |
| CPTAC phosphoproteomics | Spectrum-level evidence exists | Heavier engineering — raw MS reprocessing needed for a true co-peptide split |
| MusiteDeep | Public, widely used, runnable | Packaged benchmarks less peptide-annotated than PhosphoLingo's |

## Verdict

**P0-C gate: PASS.** A pilot re-splitting PhosphoLingo's existing datasets into
within-protein negatives and rerunning its own evaluation code is feasible in an
estimated **2–4 days** of engineering — the peptide-level data plumbing is already
solved by that paper, which is a stronger starting position than the plan assumed.
Building an independent same-peptide/co-peptide split from PhosphoSitePlus or dbPTM
from scratch (license handling, peptide mapping) is a fallback at **1–2 weeks**, only
needed if PhosphoLingo's own data proves too narrow (e.g. single organism/kinase family)
to support a convincing cross-PTM claim.

## What changes in the plan

- Phase 2 (§2 of the main plan, "跨 PTM 类型普适化") should start from PhosphoLingo's
  code and data rather than building a phosphorylation pipeline from scratch.
- The paper must explicitly position itself relative to Zuallaert et al. 2024: same
  underlying phenomenon (detectability/protease bias corrupting PTM benchmarks), but
  this project's contribution is (a) extending it to persulfidation and the
  within-protein/homology axis specifically, not just the co-peptide axis, and (b)
  pairing the diagnostic with a constructive double-axis model rather than only a
  bias correction to negative sampling.
- Recommend citing Piovesan et al. 2020 and Esmaili et al. 2022 as the "field already
  suspects this" citations in the Introduction, parallel to how the manuscript plan
  already treats Corpas et al.'s prose acknowledgment of the gate-1/gate-2 confound.
