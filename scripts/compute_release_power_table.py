#!/usr/bin/env python
"""Gate 0 — pre-registered power table for the blind primary endpoint.

Simulates one-sided Fisher's exact tests (alpha=0.05, 10,000 replicates,
seed 20260813) across candidate count K in {50, 100, 200, 300}, control
ratio r in {1, 2, 5} (total assayed sites = (1+r)*K), background
confirmation rate p0 in {0.01, 0.02, 0.05} among matched controls, and
assumed odds ratios OR in {1.5, 2.0, 3.0, 5.0}.

This table is an INPUT to the SAP co-signature: the wet-lab partner
chooses K and the control ratio (assayable site count) BEFORE
unblinding; the full table is reported, not filtered. The model has
never seen blind tomato data, so the table states an assumption grid,
not a claim.

Key reading: for a fixed assay budget, a larger K with fewer controls
per candidate dominates a smaller K with more controls (the tested
group's size drives enrichment power).

Writes power_table.json into the release package.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.stats import fisher_exact

RELEASE_DIR = Path(__file__).resolve().parents[1] / "results" / "candidates" / "multispecies_v2_candidate_release_v1"
N_REPLICATES = 10_000
SEED = 20260813
ALPHA = 0.05
K_VALUES = (50, 100, 200, 300)
P0_VALUES = (0.01, 0.02, 0.05)
OR_VALUES = (1.5, 2.0, 3.0, 5.0)
CONTROL_RATIOS = (1, 2, 5)


def _odds_to_rate(p0: float, odds_ratio: float) -> float:
    """Convert background rate + odds ratio to the candidate-group rate."""
    odds0 = p0 / (1.0 - p0)
    odds_c = odds0 * odds_ratio
    return odds_c / (1.0 + odds_c)


def _simulate_power(
    p0: float, odds_ratio: float, k: int, ratio: int, rng: np.random.Generator
) -> float:
    p_c = _odds_to_rate(p0, odds_ratio)
    n_controls = ratio * k
    rejections = 0
    for _ in range(N_REPLICATES):
        n_cand_hits = rng.binomial(k, p_c)
        n_ctrl_hits = rng.binomial(n_controls, p0)
        table = [[n_cand_hits, k - n_cand_hits], [n_ctrl_hits, n_controls - n_ctrl_hits]]
        _odds_ratio, p_value = fisher_exact(table, alternative="greater")
        if p_value < ALPHA:
            rejections += 1
    return rejections / N_REPLICATES


def main() -> None:
    rng = np.random.default_rng(SEED)
    cells: list[dict[str, float | int]] = []
    for k in K_VALUES:
        for ratio in CONTROL_RATIOS:
            for p0 in P0_VALUES:
                for odds_ratio in OR_VALUES:
                    power = _simulate_power(p0, odds_ratio, k, ratio, rng)
                    cells.append(
                        {
                            "k": k,
                            "controls_per_candidate": ratio,
                            "total_assayed_sites": (1 + ratio) * k,
                            "p0_control_rate": p0,
                            "assumed_odds_ratio": odds_ratio,
                            "p_candidate_rate": round(_odds_to_rate(p0, odds_ratio), 6),
                            "expected_candidate_hits": round(
                                k * _odds_to_rate(p0, odds_ratio), 2
                            ),
                            "power": round(power, 4),
                        }
                    )
    output = {
        "document": "power_table",
        "release_id": "multispecies-v2-candidate-release-v1",
        "status": "draft_pending_co_signature",
        "test": "one-sided Fisher exact, alpha=0.05, alternative='greater'",
        "n_replicates": N_REPLICATES,
        "simulation_seed": SEED,
        "design": {
            "candidate_group_n": "K",
            "control_group_n": "ratio*K (ratio in {1,2,5})",
            "total_assayed_sites": "(1+ratio)*K",
            "table": "2x2: group x confirmed/not-confirmed",
            "recommended_design": {
                "k": 200,
                "controls_per_candidate": 2,
                "total_assayed_sites": 600,
                "rationale": (
                    "best power per assayed site: at a 600-site budget, K=200 "
                    "with 1:2 controls (OR=3 power 0.690) dominates K=100 with "
                    "1:5 controls (0.526); the candidate group size drives "
                    "enrichment power more than control precision"
                ),
                "status": "draft_pending_lab_confirmation",
            },
        },
        "decision_rule": (
            "K and the control ratio are fixed BEFORE unblinding as the "
            "largest assayable design; this table is reported in full "
            "regardless of the chosen design"
        ),
        "cells": cells,
        "created": datetime.now(timezone.utc).isoformat(),
    }
    path = RELEASE_DIR / "power_table.json"
    path.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {path}")
    for cell in cells:
        if cell["p0_control_rate"] == 0.02 and cell["assumed_odds_ratio"] in (2, 3, 5):
            print(
                f"  K={cell['k']:>3} ratio={cell['controls_per_candidate']} "
                f"total={cell['total_assayed_sites']:>4} OR={cell['assumed_odds_ratio']} "
                f"power={cell['power']:.3f} exp_hits={cell['expected_candidate_hits']}"
            )


if __name__ == "__main__":
    main()
