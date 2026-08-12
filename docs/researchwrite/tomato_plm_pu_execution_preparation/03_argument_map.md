# Argument map

## Scientific tension

The current project has a leakage-controlled PU ranking framework and a real tomato site panel, but the small single-study label set limits model capacity. Protein language models may add useful residue context, yet a naive deep model can memorize homologs or study-specific detection patterns and produce an impressive but non-transferable score.

## Central research question

Does a provenance-locked frozen PLM residual improve tomato panel-arena Top-K ranking over the admitted v2 model under identical repeated homology-cluster evaluation, without using hard negatives or weakening stability gates?

## Central hypothesis

A low-capacity PLM residual will add complementary sequence context to the interpretable PU score, but it should replace the v2 model only if all paired admission metrics and leave-one-repetition stability checks pass.

## Supporting arguments

1. The panel arena suppresses protein-level detection shortcuts; limitation: it still represents one study and chemistry.
2. Frozen PLM embeddings increase representation capacity without end-to-end fine-tuning; limitation: pretrained homology information can still inflate weak splits.
3. The residual form preserves a zero-complexity fallback; limitation: residual scale and projection must be selected inside training folds.
4. A prospective blind release measures the actual utility of Top-K prioritization; limitation: its sample size cannot be chosen from a favorable score distribution.

## Counterarguments and responses

- A BERT-GRU stack may be more expressive. Response: include a bounded recurrent ablation only after the linear and residual PLM arms; do not make architectural depth the primary claim.
- The proteome arena is closer to deployment. Response: retain it as an audit, but select models on the panel arena to avoid protein-observation shortcuts.
- A single strong hit may be enough for a mechanism paper. Response: it may support a focused mechanism story, but it cannot establish prospective ranking enrichment or generality.

## Final move

The project should finish and audit v2 first, pre-register one bounded PLM-residual comparison, select a final model once, and only then determine and freeze the blind wet-lab package.

