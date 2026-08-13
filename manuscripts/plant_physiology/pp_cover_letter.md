# Cover letter — Plant Physiology submission

Manuscript title: A frozen, preregistered framework for proteome-wide ranking of candidate cysteine persulfidation sites in tomato (*Solanum lycopersicum*)

Manuscript type: Research Article

---

Dear Editors,

We are pleased to submit our Research Article, "A frozen, preregistered framework for proteome-wide ranking of candidate cysteine persulfidation sites in tomato (*Solanum lycopersicum*)", for consideration in Plant Physiology.

Hydrogen sulfide (H2S) signaling in plants acts largely through persulfidation of cysteine residues, yet site-level experimental mapping remains costly and low-throughput, and computational tools for candidate ranking are scarce, with published performance claims resting almost exclusively on within-dataset splits. Our manuscript responds to this situation with a transparency-first framework: a resource and protocol rather than a new performance claim.

The study reports five connected results. First, a provenance-tracked, four-species dataset of 389,609 cysteine sites with 2,334 experimentally supported positives, integrating tomato, Arabidopsis, rice, and *Magnaporthe oryzae* persulfidomes. Second, a within-dataset benchmark of six models under literature-comparable splits, in which a gated-fusion positive-unlabeled ranker leads (macro average precision 0.386; next best 0.047). Third, a one-shot, homology-cluster frozen test scored under audit (pooled delta +0.140 over the strongest baseline, 95% CI 0.118-0.170); the tomato-specific outcome is reported without inflation: the frozen model does not yet concentrate previously unseen tomato sites at the top of the list, and this is stated explicitly. Fourth, a SHA256-registered frozen release of 179,736 scored tomato candidate sites with matched controls, accompanied by a pre-blind calibration against published sites. Fifth, a co-signed, preregistered blind-cohort design (200 candidates, 393 matched controls; one-sided Fisher exact test; full power tabulation) whose outcome will be reported regardless of direction.

We believe the manuscript fits Plant Physiology because it speaks directly to the plant redox community (H2S signaling, persulfidation, tomato ripening), and because the preregistered resource-and-protocol structure offers the field a template for evaluation hygiene in post-translational-modification candidate ranking. All artifacts, exclusion ledgers, and negative results are openly archived; the preregistration will be public before the blind cohort is assayed.

One-sentence summary: A frozen, hash-registered, preregistered resource for ranking persulfidation candidates in tomato, with a co-signed blind cohort whose outcome will be reported regardless of direction.

This manuscript has not been published and is not under consideration elsewhere. All authors have approved the submission and declare no competing interests. We are glad to suggest qualified reviewers on request.

Thank you for your consideration.

Sincerely,

Guohao Lu, on behalf of all authors

Corresponding authors: Ailian Zhou (zhouailian@caas.cn), Lichuan Gu (glc@ahau.edu.cn)

Affiliations: Anhui Agricultural University; Chinese Academy of Agricultural Sciences
