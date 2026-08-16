#!/usr/bin/env python
"""Within-protein localization: v1 (masked) vs v2 (per-species scaled) bundle.

The stated task is site LOCALIZATION within a known protein — "which Cys of
THIS protein" — not proteome-wide Top-K. This diagnostic compares the two
bundles on the metric that matters for that task: the within-protein rank of
published sites, for (a) the 12 registered known controls and (b) the
kiae271 tomato sites.

Both bundles score the same target proteins' full Cys sets; v1 uses the
masked-structure semantics (its release behaviour), v2 uses the v3 structure
registry + per-species scalers. The comparison is paired per protein:
within-protein rank of the published site under each bundle, plus
aggregate Hit@2 / mean-rank / first-hit-burden.
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

V1_BUNDLE = (
    _REPO_ROOT / "results" / "candidates" / "multispecies_v2_candidate_release_v1"
    / "model_weights" / "structure_ranker_bundle.pt"
)
V2_RELEASE = (
    _REPO_ROOT / "results" / "candidates" / "multispecies_v2_candidate_release_v2"
)
V2_BUNDLE = V2_RELEASE / "model_weights" / "structure_ranker_bundle.pt"
V2_SCALERS = V2_RELEASE / "model_weights" / "species_struct_scalers.json"
OUTPUT = _REPO_ROOT / "results" / "diagnostics" / "within_protein_v1_v2.json"

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


def _load_kiae271_sites() -> dict[str, list[int]]:
    """protein_accession -> published positions, from the recovery rows TSV."""
    rows_path = (
        _REPO_ROOT
        / "results"
        / "known_controls"
        / "kiae271_release_bundle_recovery_v1.rows.tsv"
    )
    by_protein: dict[str, list[int]] = {}
    with rows_path.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            by_protein.setdefault(row["protein_accession"], []).append(
                int(row["cys_position"])
            )
    return by_protein


def _score_protein(
    bundle: StructureRankerBundle,
    proteomes: dict[str, dict[str, str]],
    species: str,
    accession: str,
    positions: tuple[int, ...],
    *,
    species_scalers: dict | None,
    structure_features: dict[tuple[str, int], tuple[tuple[float, float], bool]],
) -> dict[int, float]:
    keys = [(f"{species}|{accession}", p) for p in positions]
    raw_features = sequence_feature_map(
        [make_site_row(species, accession, p) for p in positions], proteomes
    )
    raw_struct = [[float(v) for v in structure_features[k][0]] for k in keys]
    masks = [structure_features[k][1] for k in keys]
    struct_values = raw_struct
    if species_scalers is not None:
        struct_values = transform_species_struct(
            raw_struct, masks, [species] * len(keys), species_scalers
        )
    features = BranchFeatures(
        sequence=[[float(v) for v in raw_features[k]] for k in keys],
        esm=[[0.0] for _ in keys],
        structure=struct_values,
        structure_mask=masks,
        study_ids=None,
    )
    output = score_structure_ranker_bundle(bundle, features, device_name="cpu")
    return {
        p: float(score) for p, score in zip(positions, output.scores, strict=True)
    }


def main() -> None:
    proteomes = {
        species: _load_proteome(path) for species, path in PROTEOME_PATHS.items()
    }
    v1 = StructureRankerBundle.load(V1_BUNDLE)
    v2 = StructureRankerBundle.load(V2_BUNDLE)
    v2_scalers = load_species_struct_scalers(V2_SCALERS)

    targets: list[dict] = []
    for control in REGISTERED_CONTROLS:
        if control.status != "mapped":
            continue
        species = (
            "tomato" if control.control_species == "Solanum lycopersicum"
            else "arabidopsis"
        )
        targets.append(
            {
                "label": control.gene,
                "species": species,
                "accession": control.uniprot_accession,
                "positions": [control.cys_position],
                "source": "registered_control",
            }
        )
    for accession, positions in _load_kiae271_sites().items():
        targets.append(
            {
                "label": accession,
                "species": "tomato",
                "accession": accession,
                "positions": positions,
                "source": "kiae271",
            }
        )

    # dedupe by (species, accession)
    seen: dict[tuple[str, str], dict] = {}
    for target in targets:
        key = (target["species"], target["accession"])
        bucket = seen.setdefault(key, dict(target))
        for p in target["positions"]:
            if p not in bucket["positions"]:
                bucket["positions"].append(p)
    targets = list(seen.values())

    # One structure-registry audit for ALL target keys (builds every feature
    # once, not once per protein per bundle).
    all_keys = {
        (f"{t['species']}|{t['accession']}", p)
        for t in targets
        for p in range(1, len(proteomes[t["species"]][t["accession"]]) + 1)
        if proteomes[t["species"]][t["accession"]][p - 1] == "C"
    }
    structure_features = build_registered_structure_features(
        all_keys,
        registry_path=STRUCTURE_REGISTRY,
        registry_base=STRUCTURE_REGISTRY_BASE,
    )
    print(f"structure features built for {len(all_keys)} Cys sites "
          f"({sum(1 for v in structure_features.values() if v[1])} with structure)")

    results: list[dict] = []
    for target in sorted(targets, key=lambda t: (t["species"], t["accession"])):
        sequence = proteomes[target["species"]][target["accession"]]
        positions = tuple(
            i + 1 for i, residue in enumerate(sequence) if residue == "C"
        )
        if not positions:
            continue
        true = set(target["positions"])
        s1 = _score_protein(
            v1, proteomes, target["species"], target["accession"], positions,
            species_scalers=None, structure_features=structure_features,
        )
        s2 = _score_protein(
            v2, proteomes, target["species"], target["accession"], positions,
            species_scalers=v2_scalers, structure_features=structure_features,
        )
        m1 = within_protein_metrics(
            positions=positions, scores=s1, true_positions=true
        )
        m2 = within_protein_metrics(
            positions=positions, scores=s2, true_positions=true
        )
        record = {
            "label": target["label"],
            "species": target["species"],
            "accession": target["accession"],
            "source": target["source"],
            "n_cys": len(positions),
            "true_positions": sorted(true),
            "v1_ranks": m1["true_site_ranks"],
            "v2_ranks": m2["true_site_ranks"],
            "v1_hit2": m1["hit_at_2"],
            "v2_hit2": m2["hit_at_2"],
            "v1_burden": m1["first_hit_burden"],
            "v2_burden": m2["first_hit_burden"],
            "v1_random_burden": m1["random_baseline_burden"],
        }
        results.append(record)
        print(
            f"{target['source']:<18} {target['label']:<12} "
            f"({target['species'][:4]}|{target['accession']}, n={len(positions)}): "
            f"v1={m1['true_site_ranks']} v2={m2['true_site_ranks']}"
        )

    def agg(rows: list[dict]) -> dict:
        hit2 = sum(1 for r in rows if r["v1_hit2"]), sum(
            1 for r in rows if r["v2_hit2"]
        )
        best_v1 = [
            min(r["v1_ranks"].values()) for r in rows if r["v1_ranks"]
        ]
        best_v2 = [
            min(r["v2_ranks"].values()) for r in rows if r["v2_ranks"]
        ]
        return {
            "n": len(rows),
            "hit2_v1_vs_v2": hit2,
            "mean_best_rank_v1": round(float(np.mean(best_v1)), 3) if best_v1 else None,
            "mean_best_rank_v2": round(float(np.mean(best_v2)), 3) if best_v2 else None,
        }

    document = {
        "track": "within_protein_v1_v2",
        "claim_class": "diagnostic_only",
        "note": (
            "Within-protein site-localization comparison (the stated task): rank"
            " of the published site among its protein's Cys, v1 (masked) "
            "vs v2 (per-species scaling + v3 registry). Paired per protein."
        ),
        "all": agg(results),
        "registered_controls": agg(
            [r for r in results if r["source"] == "registered_control"]
        ),
        "kiae271": agg([r for r in results if r["source"] == "kiae271"]),
        "per_protein": results,
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(
        json.dumps(document, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(f"\nwrote {OUTPUT}")
    print("aggregate:", json.dumps(document["all"], ensure_ascii=False))
    print("controls:", json.dumps(document["registered_controls"], ensure_ascii=False))
    print("kiae271:", json.dumps(document["kiae271"], ensure_ascii=False))


if __name__ == "__main__":
    main()
