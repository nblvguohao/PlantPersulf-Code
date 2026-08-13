#!/usr/bin/env python
"""Gate 0 — pre-registered power table for the blind primary endpoint.

Simulates one-sided Fisher's exact tests (alpha=0.05, 10,000 replicates,
seed 20260813) across K in {50, 100, 200}, background confirmation rate
p0 in {0.01, 0.02, 0.05} among matched controls (5 controls per
candidate), and assumed odds ratios OR in {1.5, 2.0, 3.0, 5.0}.

This table is an INPUT to the SAP co-signature: the wet-lab partner
chooses K (assayable candidate count) BEFORE unblinding; the full table
is reported, not filtered. The model has never seen blind tomato data,
so the table states an assumption grid, not a claim.

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
K_VALUES = (50, 100, 200)
P0_VALUES = (0.01, 0.02, 0.05)
OR_VALUES = (1.5, 2.0, 3.0, 5.0)
CONTROLS_PER_CANDIDATE = 5


def _odds_to_rate(p0: float, odds_ratio: float) -> float:
    """Convert background rate + odds ratio to the candidate-group rate."""
    odds0 = p0 / (1.0 - p0)
    odds_c = odds0 * odds_ratio
    return odds_c / (1.0 + odds_c)


def _simulate_power(p0: float, odds_ratio: float, k: int, rng: np.random.Generator) -> float:
    p_c = _odds_to_rate(p0, odds_ratio)
    n_controls = CONTROLS_PER_CANDIDATE * k
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
        for p0 in P0_VALUES:
            for odds_ratio in OR_VALUES:
                power = _simulate_power(p0, odds_ratio, k, rng)
                cells.append(
                    {
                        "k": k,
                        "controls_per_candidate": CONTROLS_PER_CANDIDATE,
                        "p0_control_rate": p0,
                        "assumed_odds_ratio": odds_ratio,
                        "p_candidate_rate": round(_odds_to_rate(p0, odds_ratio), 6),
                        "expected_candidate_hits": round(k * _odds_to_rate(p0, odds_ratio), 2),
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
            "control_group_n": "5*K (CONTROLS_PER_CANDIDATE=5)",
            "table": "2x2: group x confirmed/not-confirmed",
        },
        "decision_rule": (
            "K is fixed BEFORE unblinding as the largest assayable K; this "
            "table is reported in full regardless of the chosen K"
        ),
        "cells": cells,
        "created": datetime.now(timezone.utc).isoformat(),
    }
    path = RELEASE_DIR / "power_table.json"
    path.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {path}")
    for cell in cells:
        if cell["k"] in (50, 100, 200) and cell["p0_control_rate"] == 0.02:
            print(
                f"  K={cell['k']:>3} p0={cell['p0_control_rate']} OR={cell['assumed_odds_ratio']} "
                f"power={cell['power']:.3f} expected_hits={cell['expected_candidate_hits']}"
            )


if __name__ == "__main__":
    main()
