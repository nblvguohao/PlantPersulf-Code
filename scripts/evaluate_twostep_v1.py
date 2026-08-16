#!/usr/bin/env python
"""Two-step model v1: does stage-1 oxidation signal help within-protein ranking?

The two-step decomposition P(persulfidation) = P(oxidizable) x P(persulfidated |
oxidized) predicts that gate-1 (oxidation) dominates the variance — the
co-peptide null and the oxiPTM discrimination both point there. This
diagnostic builds the first stage-1 classifier and asks, on the honest
within-protein evidence (12 gold-standard controls + clean untrained kiae271
sites), whether the oxidation score adds ranking signal to the persulfidation
score:

- ``persulf_v2``: the v2 per-species-scaled bundle's within-protein ranking
  (the current best persulfidation signal).
- ``oxid_only``: the stage-1 oxidation classifier's within-protein ranking —
  if gate-1 dominates, this alone should surface persulfidation sites.
- ``combined``: rank by persulf_v2 score x stage-1 oxidation probability.

Stage-1 is a logistic regression over the project's sequence features
(hydrophobicity, cys_density, positive charge density) + structure contact/
pLDDT where available, trained on the 2941 verified gate-1 S-nitrosylation
sites vs a background sample of non-gate-1 Cys. Cross-species warning: stage-1
is Arabidopsis-dominated, the test includes tomato — the LOSO result says
transfer may be ~0, which this test measures directly.

diagnostic_only; no frozen artifact changes.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import csv  # noqa: E402

from plantpersulf.evaluation.comparable_track import (  # noqa: E402
    build_registered_structure_features,
)
from plantpersulf.evaluation.known_controls import REGISTERED_CONTROLS  # noqa: E402
from plantpersulf.evaluation.release_scoring import make_site_row  # noqa: E402
from plantpersulf.evaluation.species_structure_scaling import (  # noqa: E402
    load_species_struct_scalers,
    transform_species_struct,
)
from plantpersulf.evaluation.within_protein_ranking import (  # noqa: E402
    within_protein_metrics,
)
from plantpersulf.features.sequence import _load_proteome  # noqa: E402
from plantpersulf.models.structure_ranker import (  # noqa: E402
    BranchFeatures,
    StructureRankerBundle,
    score_structure_ranker_bundle,
)
from plantpersulf.proteomics.multispecies_v2_sources import (  # noqa: E402
    sequence_feature_map,
)

V2_BUNDLE = (
    _REPO_ROOT / "results" / "candidates" / "multispecies_v2_candidate_release_v2"
    / "model_weights" / "structure_ranker_bundle.pt"
)
V2_SCALERS = (
    _REPO_ROOT / "results" / "candidates" / "multispecies_v2_candidate_release_v2"
    / "model_weights" / "species_struct_scalers.json"
)
GATE1_REGISTRY = _REPO_ROOT / "data" / "registry" / "gate1_oxidation_sites_v1.tsv"
OUTPUT = _REPO_ROOT / "results" / "diagnostics" / "twostep_v1.json"

PROTEOME_PATHS = {
    "arabidopsis": (
        _REPO_ROOT / "data" / "raw" / "references" / "arabidopsis_ref_proteome_v2.fasta"
    ),
    "tomato": (
        _REPO_ROOT / "data" / "raw" / "references" / "tomato_ref_proteome_v1.fasta"
    ),
    "rice": (
        _REPO_ROOT / "data" / "raw" / "references" / "rice_proteome_v1"
        / "uniprot_rice_v1.fasta"
    ),
}
STRUCTURE_REGISTRY = (
    _REPO_ROOT / "data" / "registry" / "releases"
    / "alphafold_structures_release_v3.tsv"
)
STRUCTURE_REGISTRY_BASE = _REPO_ROOT / "data" / "registry"

BACKGROUND_PER_POSITIVE = 8


def _feature_matrix(
    rows: list[tuple[str, str, int]],
    proteomes: dict[str, dict[str, str]],
    structure_features: dict[tuple[str, int], tuple[tuple[float, float], bool]],
    species_scalers: dict | None,
) -> tuple[np.ndarray, list[bool]]:
    """(n, 5) feature matrix: 3 sequence + 2 structure (scaled)."""
    raw = sequence_feature_map(
        [make_site_row(species, acc, pos) for species, acc, pos in rows], proteomes
    )
    seq = np.array(
        [[float(v) for v in raw[(f"{s}|{a}", p)]] for s, a, p in rows],
        dtype=float,
    )
    keys = [(f"{s}|{a}", p) for s, a, p in rows]
    raw_struct = [[float(v) for v in structure_features[k][0]] for k in keys]
    masks = [structure_features[k][1] for k in keys]
    struct = np.array(raw_struct, dtype=float)
    if species_scalers is not None:
        species_list = [s for s, _a, _p in rows]
        scaled = transform_species_struct(
            raw_struct, masks, species_list, species_scalers
        )
        struct = np.array(scaled, dtype=float)
    X = np.hstack([seq, struct])
    return X, masks


def main() -> None:
    proteomes = {
        species: _load_proteome(path) for species, path in PROTEOME_PATHS.items()
    }
    v2 = StructureRankerBundle.load(V2_BUNDLE)
    v2_scalers = load_species_struct_scalers(V2_SCALERS)

    # --- gate-1 positives (verified rows) -----------------------------------
    gate1 = []
    with GATE1_REGISTRY.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["status"] == "verified":
                gate1.append(
                    (row["species"], row["accession"], int(row["cys_position"]))
                )
    gate1_set = set(gate1)
    print(f"gate-1 verified positives: {len(gate1)}")

    # --- background: non-gate-1 Cys, sampled across the same proteomes ------
    rng = np.random.RandomState(20260817)
    background: list[tuple[str, str, int]] = []
    gate_by_species: dict[str, list[tuple[str, str, int]]] = {}
    for s, a, p in gate1:
        gate_by_species.setdefault(s, []).append((s, a, p))
    for species, positives in gate_by_species.items():
        proteome = proteomes[species]
        pool = [
            (species, acc, pos)
            for acc, seq in proteome.items()
            for pos in range(1, len(seq) + 1)
            if seq[pos - 1] == "C" and (species, acc, pos) not in gate1_set
        ]
        n_bg = min(len(positives) * BACKGROUND_PER_POSITIVE, len(pool))
        idx = rng.choice(len(pool), size=n_bg, replace=False)
        background.extend(pool[int(i)] for i in idx)
    print(f"background sampled: {len(background)}")

    # --- within-protein test proteins (honest evidence) ---------------------
    # Collected FIRST so their structure keys are in the single feature build.
    targets: dict[tuple[str, str], dict] = {}
    for control in REGISTERED_CONTROLS:
        if control.status != "mapped":
            continue
        species = (
            "tomato" if control.control_species == "Solanum lycopersicum"
            else "arabidopsis"
        )
        targets.setdefault((species, control.uniprot_accession), {
            "label": control.gene, "positions": set(), "source": "registered_control"
        })["positions"].add(control.cys_position)
    rows_path = (
        _REPO_ROOT / "results" / "known_controls"
        / "kiae271_release_bundle_recovery_v1.rows.tsv"
    )
    with rows_path.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["in_candidate_table"] != "True":
                continue
            key = ("tomato", row["protein_accession"])
            targets.setdefault(
                key,
                {
                    "label": row["protein_accession"],
                    "positions": set(),
                    "source": "kiae271_clean",
                },
            )["positions"].add(int(row["cys_position"]))

    # --- structure features (one audit for ALL keys: train + test) -----------
    test_keys = set()
    for (species, accession), _bucket in targets.items():
        sequence = proteomes[species][accession]
        test_keys.update(
            (f"{species}|{accession}", p)
            for p in range(1, len(sequence) + 1)
            if sequence[p - 1] == "C"
        )
    all_keys = test_keys | set(
        (f"{s}|{a}", p) for s, a, p in [*gate1, *background]
    )
    structure_features = build_registered_structure_features(
        all_keys,
        registry_path=STRUCTURE_REGISTRY,
        registry_base=STRUCTURE_REGISTRY_BASE,
    )

    X_pos, _mask_pos = _feature_matrix(
        gate1, proteomes, structure_features, v2_scalers
    )
    X_bg, _mask_bg = _feature_matrix(
        background, proteomes, structure_features, v2_scalers
    )
    X = np.vstack([X_pos, X_bg])
    y = np.array([1] * len(gate1) + [0] * len(background), dtype=int)
    print(f"stage-1 training matrix: {X.shape}")

    # --- stage-1 logistic regression ----------------------------------------
    from sklearn.linear_model import LogisticRegression

    model = LogisticRegression(max_iter=2000, class_weight="balanced")
    model.fit(X, y)
    train_auc = _auc(X, y, model)
    print(f"stage-1 train AUC: {train_auc:.4f}")

    results: list[dict] = []
    for (species, accession), bucket in sorted(targets.items()):
        sequence = proteomes[species][accession]
        positions = tuple(
            i + 1 for i, residue in enumerate(sequence) if residue == "C"
        )
        if not positions:
            continue
        rows = [(species, accession, p) for p in positions]
        X_test, _masks = _feature_matrix(
            rows, proteomes, structure_features, v2_scalers
        )
        oxid_prob = model.predict_proba(X_test)[:, 1]

        # persulfidation scores from the v2 bundle
        keys = [(f"{species}|{accession}", p) for p in positions]
        raw = sequence_feature_map(
            [make_site_row(species, accession, p) for p in positions], proteomes
        )
        structure = {k: structure_features[k] for k in keys}
        raw_struct = [[float(v) for v in structure[k][0]] for k in keys]
        masks = [structure[k][1] for k in keys]
        scaled = transform_species_struct(
            raw_struct, masks, [species] * len(keys), v2_scalers
        )
        feats = BranchFeatures(
            sequence=[[float(v) for v in raw[k]] for k in keys],
            esm=[[0.0] for _ in keys],
            structure=scaled,
            structure_mask=masks,
            study_ids=None,
        )
        out = score_structure_ranker_bundle(v2, feats, device_name="cpu")
        persulf_scores = {
            p: float(s) for p, s in zip(positions, out.scores, strict=True)
        }

        combined = {
            p: persulf_scores[p] * float(oxid_prob[i])
            for i, p in enumerate(positions)
        }
        oxid_only = {p: float(oxid_prob[i]) for i, p in enumerate(positions)}
        true = bucket["positions"]
        metrics = {
            name: within_protein_metrics(
                positions=positions, scores=scores, true_positions=true
            )
            for name, scores in (
                ("persulf_v2", persulf_scores),
                ("oxid_only", oxid_only),
                ("combined", combined),
            )
        }
        results.append(
            {
                "label": bucket["label"],
                "species": species,
                "accession": accession,
                "source": bucket["source"],
                "n_cys": len(positions),
                "true_positions": sorted(true),
                "hit2": {
                    name: m["hit_at_2"] for name, m in metrics.items()
                },
                "best_rank": {
                    name: min(m["true_site_ranks"].values())
                    for name, m in metrics.items()
                },
            }
        )
        print(
            f"{bucket['source'][:16]:<18} {bucket['label']:<12} n={len(positions)} "
            f"hit2: persulf={metrics['persulf_v2']['hit_at_2']} "
            f"oxid={metrics['oxid_only']['hit_at_2']} "
            f"combined={metrics['combined']['hit_at_2']}"
        )

    def agg(rows: list[dict], name: str) -> dict:
        return {
            "n": len(rows),
            "hit2_persulf_v2": sum(1 for r in rows if r["hit2"]["persulf_v2"]),
            "hit2_oxid_only": sum(1 for r in rows if r["hit2"]["oxid_only"]),
            "hit2_combined": sum(1 for r in rows if r["hit2"]["combined"]),
        }

    document = {
        "track": "twostep_v1",
        "claim_class": "diagnostic_only",
        "note": (
            "Stage-1 oxidation classifier (logistic regression, 3 seq + 2 struct "
            "features) trained on 2941 verified gate-1 SNO sites vs background. "
            "Within-protein Hit@2 of persulfidation-v2 score, oxidation-only "
            "score, and their product, on honest evidence."
        ),
        "stage1": {
            "n_positives": len(gate1),
            "n_background": len(background),
            "train_auc": round(float(train_auc), 4),
            "features": [
                "hydrophobicity", "cys_density", "pos_charge", "contact", "plddt",
            ],
        },
        "all": agg(results, "all"),
        "registered_controls": agg(
            [r for r in results if r["source"] == "registered_control"], "controls"
        ),
        "kiae271_clean": agg(
            [r for r in results if r["source"] == "kiae271_clean"], "kiae271"
        ),
        "per_protein": results,
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(
        json.dumps(document, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(f"\nwrote {OUTPUT}")
    print("all:", json.dumps(document["all"], ensure_ascii=False))


def _auc(X: np.ndarray, y: np.ndarray, model) -> float:
    from sklearn.metrics import roc_auc_score

    return float(roc_auc_score(y, model.predict_proba(X)[:, 1]))


if __name__ == "__main__":
    main()
