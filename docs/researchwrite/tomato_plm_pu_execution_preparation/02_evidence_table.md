# Evidence table

| Claim | Evidence/source | Strength | Usable section | Risk | Status |
|---|---|---|---|---|---|
| v2 is a development-only tomato ranking run | Frozen config claim class and Gate 2 firewall test | evidence-backed | Phase A | Mislabeling it as external validation | evidence-backed |
| Panel arena is the correct primary arena | Frozen config; tomato-local analysis showing protein-detection confounding in proteome arena | evidence-backed | Evaluation | Single-study supervision remains | evidence-backed |
| A complex model must beat a simple admitted baseline on all required paired metrics | Existing `admit_candidate_model` implementation | evidence-backed | Admission | Low power may prefer the simpler model | evidence-backed |
| Frozen PLM features may improve residue representation | Existing ESM baseline plus current PTM literature | plausible-inference | Phase B | Homology memorization and small-sample overfit | plausible-inference |
| A gated residual is safer than replacing the biological baseline outright | Small-sample architecture argument; residual can be set to zero | hypothesis | Phase B | Still requires nested model selection | hypothesis |
| BiGRU is not automatically useful after a full-sequence PLM | Architectural redundancy argument | hypothesis | Optional ablation | Could miss directional local patterns | hypothesis |
| Structure should remain disabled unless it adds matched-coverage panel value | Existing structure-confound diagnosis and structure admission policy | evidence-backed | Optional branch | Coverage availability can become a shortcut | evidence-backed |
| Blind wet-lab enrichment is needed for prospective value | Project scientific-integrity and candidate-release policies | evidence-backed | Phase C | Collaborator capacity and protocol not yet frozen | evidence-backed |
| One or two post-hit mechanisms can deepen the story but cannot rescue a failed blind primary endpoint | Candidate-release policy and project claim hierarchy | evidence-backed | Phase C | Cherry-picking individual hits | evidence-backed |

