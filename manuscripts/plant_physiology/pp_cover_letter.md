# Cover letter — Plant Physiology submission

Manuscript title: Cross-kingdom conservation of persulfidation targeting and a frozen, preregistered framework for ranking candidate cysteine persulfidation sites in tomato (*Solanum lycopersicum*)

Manuscript type: Research Article

---

Dear Editors,

We are pleased to submit our Research Article, "Cross-kingdom conservation of persulfidation targeting and a frozen, preregistered framework for ranking candidate cysteine persulfidation sites in tomato (*Solanum lycopersicum*)", for consideration in Plant Physiology.

Hydrogen sulfide (H2S) signaling in plants acts largely through persulfidation of cysteine residues. Two open questions organize this manuscript. First: is persulfidation targeting itself a conserved biological signal, or an idiosyncratic overlay on each organism's redox chemistry? We answer with the first direct cross-kingdom test of the question. Second: given that site-level mapping is costly and published computational tools rest almost exclusively on within-dataset splits, how should a candidate-ranking resource be built and evaluated so that its claims mean what they say? We answer with a frozen, preregistered design whose primary result is decided by a blind experiment.

The conservation analysis integrates four independent laboratory/chemistry/species datasets spanning two kingdoms (Arabidopsis, rice, tomato, and *Magnaporthe oryzae*). Co-persulfidation of orthologous PANTHER subfamilies is enriched beyond chance in two of the three deeply sampled pairwise comparisons after multiple-testing correction (2.96x and 2.22x), and eight subfamilies — six of them core redox or central-carbon-metabolism enzymes — are persulfidated in Arabidopsis, rice, and *Magnaporthe* alike (8.6x the independence expectation, permutation p = 1.0e-4). Every non-significant cell, including all tomato comparisons, is reported with an explicit detectability floor, so the absence of p < 0.05 cannot be read as an absence of conservation. In a structural-context analysis we found an apparent burial/exposure signature of persulfidated cysteines that did not survive a model-confidence control: persulfidated and other cysteines differ systematically in AlphaFold pLDDT, and restricting the comparison to confident residues removes the signal in all four species. We report this correction in the manuscript because we believe it is a useful caution for structure-based reanalyses of proteomics-derived site sets.

The resource strand reports a provenance-tracked four-species panel (389,609 sites, 2,334 positives), a literature-comparable benchmark in which a gated-fusion positive-unlabeled ranker leads five alternatives (macro AP 0.386 vs 0.047), a one-shot homology-cluster frozen test scored under audit (+0.140 pooled AP difference over the strongest baseline, 95% CI 0.118-0.170), a SHA256-registered frozen release of 179,736 scored tomato candidate sites with matched controls, and a pre-blind calibration against published sites. The blind-cohort design (200 candidates, 393 matched controls, one-sided Fisher exact test, full 144-cell power tabulation) is co-signed with the wet-laboratory partner and registered before any blind data exist; its outcome will be reported regardless of direction.

We believe the manuscript fits Plant Physiology on three grounds: the conservation result speaks directly to the plant redox and H2S-signaling community; the tomato resource addresses a crop whose persulfidome has only a single published dataset; and the two-step discipline the paper practices — state the test before seeing its answer, and report what survives scrutiny rather than what would be most publishable — offers the field a template for evaluation hygiene in post-translational-modification candidate ranking. All artifacts, exclusion ledgers, and negative results are openly archived, and the preregistration is publicly verifiable.

One-sentence summary: The first cross-kingdom test of persulfidation targeting, a self-corrected structural-context analysis, and a frozen, hash-registered, preregistered tomato candidate resource whose blind outcome will be reported regardless of direction.

This manuscript has not been published and is not under consideration elsewhere. All authors have approved the submission and declare no competing interests. We are glad to suggest qualified reviewers on request.

Thank you for your consideration.

Sincerely,

Guohao Lu, on behalf of all authors

Corresponding authors: Ailian Zhou (zhouailian@caas.cn), Lichuan Gu (glc@ahau.edu.cn)

Affiliations: Anhui Agricultural University; Chinese Academy of Agricultural Sciences
