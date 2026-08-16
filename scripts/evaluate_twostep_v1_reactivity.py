#!/usr/bin/env python
"""Two-step v1 with reactivity features: does pKa/structure reactivity transfer?

The generic-feature stage-1 (twostep_v1) trained on Arabidopsis SNO sites did
not transfer to tomato within-protein ranking. Hypothesis: the generic
features capture within-species distribution memory, not universal reaction
chemistry (LOSO). This test rebuilds stage-1 on *reactivity* features that
should be chemically universal:

- 7 structure features per Cys from the AFDB model (pLDDT, RSA, S-gamma
  geometry, local positive-residue count, S-gamma Coulomb potential, packing)
  — the feature set validated in the P2/structure-regime diagnostics;
- local positive charge density (the existing thiol-pKa proxy).

Stage-1 is trained on gate-1 SNO sites that HAVE an AFDB structure vs
structure-bearing background sites, then applied to the within-protein test
proteins. If pKa/electrostatics are universal, the reactivity stage-1 should
transfer cross-species better than the generic stage-1 (7/29 Hit@2).

diagnostic_only.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from plantpersulf.evaluation.copeptide_structure import FEATURE_NAMES  # noqa: E402
from plantpersulf.evaluation.known_controls import REGISTERED_CONTROLS  # noqa: E402
from plantpersulf.evaluation.release_scoring import make_site_row  # noqa: E402
from plantpersulf.evaluation.structure_features import (  # noqa: E402
    cys_structure_features,
    parse_pdb,
)
from plantpersulf.evaluation.within_protein_ranking import (  # noqa: E402
    within_protein_metrics,
)
from plantpersulf.features.sequence import (  # noqa: E402
    _flanking_window,
    _load_proteome,
    _local_positive_charge_density,
)
from plantpersulf.proteomics.multispecies_v2_sources import (  # noqa: E402
    sequence_feature_map,
)

GATE1_REGISTRY = _REPO_ROOT / "data" / "registry" / "gate1_oxidation_sites_v1.tsv"
ALPHAFOLD_DIR = _REPO_ROOT / "data" / "raw" / "alphafold"
OUTPUT = _REPO_ROOT / "results" / "diagnostics" / "twostep_v1_reactivity.json"

PROTEOME_PATHS = {
    "arabidopsis": (
        _REPO_ROOT / "data" / "raw" / "references" / "arabidopsis_ref_proteome_v2.fasta"
    ),
    "tomato": (
        _REPO_ROOT / "data" / "raw" / "references" / "tomato_ref_proteome_v1.fasta"
    ),
}
BACKGROUND_PER_POSITIVE = 6
PKA_FEATURE = "local_pos_charge_density"


def _pdb_path(accession: str) -> Path | None:
    candidate = ALPHAFOLD_DIR / f"AF-{accession}-F1-model.pdb"
    return candidate if candidate.is_file() else None


def _reactivity_features(
    species: str,
    accession: str,
    positions: list[int],
    proteomes: dict[str, dict[str, str]],
) -> tuple[list[dict[str, float]], bool]:
    """Per-site reactivity feature dicts; ``has_structure`` flags the protein."""
    sequence = proteomes[species][accession]
    pdb = _pdb_path(accession)
    struct: dict[int, dict[str, float]] = {}
    if pdb is not None:
        try:
            residues = parse_pdb(pdb.read_text(encoding="utf-8"))
            if len(residues) == len(sequence):
                struct = cys_structure_features(residues, positions)
        except (ValueError, KeyError, IndexError):
            struct = {}
    has_structure = len(struct) == len(positions) and bool(struct)
    out: list[dict[str, float]] = []
    for pos in positions:
        flank = _flanking_window(sequence, pos, 10)
        row: dict[str, float] = {PKA_FEATURE: _local_positive_charge_density(flank)}
        if has_structure:
            for name in FEATURE_NAMES:
                row[name] = float(struct.get(pos, {}).get(name, 0.0))
        else:
            for name in FEATURE_NAMES:
                row[name] = 0.0
        out.append(row)
    return out, has_structure


def main() -> None:
    proteomes = {
        species: _load_proteome(path) for species, path in PROTEOME_PATHS.items()
    }

    gate1 = []
    with GATE1_REGISTRY.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["status"] == "verified":
                gate1.append(
                    (row["species"], row["accession"], int(row["cys_position"]))
                )
    gate1_set = set(gate1)

    # test proteins (honest evidence)
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

    # --- reactivity features for gate-1 (only structure-bearing sites) -------
    pos_rows, pos_ok = [], 0
    for species, acc, pos in gate1:
        if _pdb_path(acc) is None:
            continue
        feats, _hs = _reactivity_features(species, acc, [pos], proteomes)
        pos_rows.append(feats[0])
        pos_ok += 1

    rng = np.random.RandomState(20260817)
    bg_rows: list[dict[str, float]] = []
    bg_ok = 0
    gate_by_species: dict[str, list[tuple[str, str, int]]] = {}
    for s, a, p in gate1:
        gate_by_species.setdefault(s, []).append((s, a, p))
    for species, positives in gate_by_species.items():
        proteome = proteomes[species]
        pool = [
            (acc, pos)
            for acc, seq in proteome.items()
            if _pdb_path(acc) is not None
            for pos in range(1, len(seq) + 1)
            if seq[pos - 1] == "C" and (species, acc, pos) not in gate1_set
        ]
        n = min(len(positives) * BACKGROUND_PER_POSITIVE, len(pool))
        idx = rng.choice(len(pool), size=n, replace=False)
        for i in idx:
            acc, pos = pool[int(i)]
            feats, _hs = _reactivity_features(species, acc, [pos], proteomes)
            bg_rows.append(feats[0])
            bg_ok += 1

    features = list(pos_rows[0])
    X_pos = np.array([[r[f] for f in features] for r in pos_rows], dtype=float)
    X_bg = np.array([[r[f] for f in features] for r in bg_rows], dtype=float)
    X = np.vstack([X_pos, X_bg])
    y = np.array([1] * len(pos_rows) + [0] * len(bg_rows), dtype=int)
    print(
        f"reactivity stage-1: {len(pos_rows)} gate-1 (+structure) "
        f"vs {len(bg_rows)} bg"
    )

    from sklearn.linear_model import LogisticRegression

    model = LogisticRegression(max_iter=2000, class_weight="balanced")
    model.fit(X, y)
    from sklearn.metrics import roc_auc_score

    auc = roc_auc_score(y, model.predict_proba(X)[:, 1])
    print(f"reactivity stage-1 train AUC: {auc:.4f}")

    # --- persulfidation baseline assets (load once) ---------------------------
    from plantpersulf.evaluation.comparable_track import (  # noqa: E402
        build_registered_structure_features,
    )
    from plantpersulf.evaluation.species_structure_scaling import (  # noqa: E402
        load_species_struct_scalers,
        transform_species_struct,
    )
    from plantpersulf.models.structure_ranker import (  # noqa: E402
        BranchFeatures,
        StructureRankerBundle,
        score_structure_ranker_bundle,
    )

    v2 = StructureRankerBundle.load(
        _REPO_ROOT / "results" / "candidates"
        / "multispecies_v2_candidate_release_v2" / "model_weights"
        / "structure_ranker_bundle.pt"
    )
    scalers = load_species_struct_scalers(
        _REPO_ROOT / "results" / "candidates"
        / "multispecies_v2_candidate_release_v2" / "model_weights"
        / "species_struct_scalers.json"
    )
    structure_registry = (
        _REPO_ROOT / "data" / "registry" / "releases"
        / "alphafold_structures_release_v3.tsv"
    )
    structure_registry_base = _REPO_ROOT / "data" / "registry"

    # --- apply to test proteins ----------------------------------------------
    results = []
    for (species, accession), bucket in sorted(targets.items()):
        sequence = proteomes[species][accession]
        positions = tuple(
            i + 1 for i, residue in enumerate(sequence) if residue == "C"
        )
        if not positions:
            continue
        feats, has_struct = _reactivity_features(
            species, accession, list(positions), proteomes
        )
        X_test = np.array([[r[f] for f in features] for r in feats], dtype=float)
        oxid_prob = model.predict_proba(X_test)[:, 1]

        # persulfidation baseline (v2 bundle, structure-masked where no model)
        raw = sequence_feature_map(
            [make_site_row(species, accession, p) for p in positions], proteomes
        )
        keys = [(f"{species}|{accession}", p) for p in positions]
        structure = build_registered_structure_features(
            set(keys),
            registry_path=structure_registry,
            registry_base=structure_registry_base,
        )
        raw_struct = [[float(v) for v in structure[k][0]] for k in keys]
        masks = [structure[k][1] for k in keys]
        scaled = transform_species_struct(
            raw_struct, masks, [species] * len(keys), scalers
        )
        feats_b = BranchFeatures(
            sequence=[[float(v) for v in raw[k]] for k in keys],
            esm=[[0.0] for _ in keys],
            structure=scaled,
            structure_mask=masks,
            study_ids=None,
        )
        out = score_structure_ranker_bundle(v2, feats_b, device_name="cpu")
        persulf = {p: float(s) for p, s in zip(positions, out.scores, strict=True)}

        combined = {
            p: persulf[p] * float(oxid_prob[i]) for i, p in enumerate(positions)
        }
        oxid_only = {p: float(oxid_prob[i]) for i, p in enumerate(positions)}
        true = bucket["positions"]
        metrics = {
            name: within_protein_metrics(
                positions=positions, scores=scores, true_positions=true
            )
            for name, scores in (
                ("persulf_v2", persulf), ("oxid_only", oxid_only),
                ("combined", combined),
            )
        }
        results.append({
            "label": bucket["label"], "species": species, "accession": accession,
            "source": bucket["source"], "n_cys": len(positions),
            "has_structure": has_struct,
            "hit2": {name: m["hit_at_2"] for name, m in metrics.items()},
        })
        print(
            f"{bucket['source'][:16]:<18} {bucket['label']:<12} "
            f"struct={'Y' if has_struct else 'N'} "
            f"hit2: persulf={metrics['persulf_v2']['hit_at_2']} "
            f"oxid={metrics['oxid_only']['hit_at_2']} "
            f"combined={metrics['combined']['hit_at_2']}"
        )

    def agg(rows, name):
        return {
            "n": len(rows),
            "n_with_structure": sum(1 for r in rows if r["has_structure"]),
            "hit2_persulf_v2": sum(1 for r in rows if r["hit2"]["persulf_v2"]),
            "hit2_oxid_only": sum(1 for r in rows if r["hit2"]["oxid_only"]),
            "hit2_combined": sum(1 for r in rows if r["hit2"]["combined"]),
        }

    document = {
        "track": "twostep_v1_reactivity",
        "claim_class": "diagnostic_only",
        "note": (
            "Stage-1 rebuilt on reactivity features (7 structure features from "
            "AFDB + local positive charge density as thiol-pKa proxy), trained "
            "on gate-1 SNO sites with structures, applied to within-protein test. "
            "Tests whether pKa/electrostatics transfer cross-species better than "
            "the generic features (twostep_v1: oxid_only 7/29)."
        ),
        "stage1": {"n_pos_with_struct": len(pos_rows), "n_bg_with_struct": len(bg_rows),
                   "train_auc": round(float(auc), 4), "features": features},
        "all": agg(results, "all"),
        "registered_controls": agg(
            [r for r in results if r["source"] == "registered_control"], "c"
        ),
        "kiae271_clean": agg(
            [r for r in results if r["source"] == "kiae271_clean"], "k"
        ),
        "per_protein": results,
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(
        json.dumps(document, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(f"\nwrote {OUTPUT}")
    print("all:", json.dumps(document["all"], ensure_ascii=False))


if __name__ == "__main__":
    main()
