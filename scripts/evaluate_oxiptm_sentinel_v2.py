#!/usr/bin/env python
"""OxiPTM sentinel on the v2 per-species-scaled bundle.

The v1 discrimination diagnostic (``evaluate_oxiptm_sites.py``, structure
branch masked — the v1 release semantics) found the frozen model ranks
S-nitrosylation sites within their proteins at least as high as
persulfidation sites (mean rank 2.0 vs 3.8, one-sided p=0.976 for
persulfidation-better). The v2 release changes the structure branch
(v3 registry + per-species scaling); this sentinel re-runs the same
within-protein comparison with the NEW bundle and the SAME verified site
list (6 S-nitrosylation + 8 persulfidation proteins), structure branch
ACTIVE with the v2 species scalers where the registry provides a model.

Interpretation: if the v2 model still cannot rank persulfidation sites
above S-nitrosylation sites within their proteins, the gate-2 gap is
unchanged (the change helped the benchmark but not the oxiPTM contrast);
if it separates, the scaling/structure activation changed what the model
learned. Small n either way — descriptive only, not a registered endpoint.
"""

from __future__ import annotations

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
from plantpersulf.evaluation.release_scoring import make_site_row  # noqa: E402
from plantpersulf.evaluation.species_structure_scaling import (  # noqa: E402
    load_species_struct_scalers,
    transform_species_struct,
)
from plantpersulf.evaluation.within_protein_ranking import (  # noqa: E402
    within_protein_metrics,
)
from plantpersulf.evidence.oxiptm_sites import (  # noqa: E402
    ARABIDOPSIS,
    PERSULFIDATION,
    RICE,
    S_NITROSYLATION,
    SITE_TABLE,
    TOMATO,
    VERIFIED,
    verify_site,
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

RELEASE_DIR = (
    _REPO_ROOT / "results" / "candidates" / "multispecies_v2_candidate_release_v2"
)
BUNDLE_PATH = RELEASE_DIR / "model_weights" / "structure_ranker_bundle.pt"
SCALERS_PATH = RELEASE_DIR / "model_weights" / "species_struct_scalers.json"
OUTPUT = _REPO_ROOT / "results" / "diagnostics" / "oxiptm_sentinel_v2.json"

PROTEOME_PATHS = {
    ARABIDOPSIS: (
        _REPO_ROOT / "data" / "raw" / "references" / "arabidopsis_ref_proteome_v2.fasta"
    ),
    TOMATO: (
        _REPO_ROOT / "data" / "raw" / "references" / "tomato_ref_proteome_v1.fasta"
    ),
    RICE: (
        _REPO_ROOT / "data" / "raw" / "references" / "rice_proteome_v1"
        / "uniprot_rice_v1.fasta"
    ),
}
STRUCTURE_REGISTRY = (
    _REPO_ROOT
    / "data"
    / "registry"
    / "releases"
    / "alphafold_structures_release_v3.tsv"
)
STRUCTURE_REGISTRY_BASE = _REPO_ROOT / "data" / "registry"


def main() -> None:
    proteomes = {
        species: _load_proteome(path) for species, path in PROTEOME_PATHS.items()
    }
    bundle = StructureRankerBundle.load(BUNDLE_PATH)
    species_scalers = load_species_struct_scalers(SCALERS_PATH)

    verified = [
        site for site in SITE_TABLE
        if verify_site(site, proteomes=proteomes).status == VERIFIED
    ]
    by_protein: dict[tuple[str, str], dict] = {}
    for site in verified:
        result = verify_site(site, proteomes=proteomes)
        assert result.accession is not None
        key = (site.species, result.accession)
        bucket = by_protein.setdefault(
            key, {"oxiptm": site.oxiptm, "positions": set(), "genes": set()}
        )
        bucket["positions"].add(site.cys_position)
        bucket["genes"].add(site.gene)

    rows: dict[str, list[dict]] = {PERSULFIDATION: [], S_NITROSYLATION: []}
    for (species, accession), bucket in sorted(by_protein.items()):
        sequence = proteomes[species][accession]
        positions = tuple(
            i + 1 for i, residue in enumerate(sequence) if residue == "C"
        )
        assert positions
        keys = [(f"{species}|{accession}", p) for p in positions]
        raw_features = sequence_feature_map(
            [make_site_row(species, accession, p) for p in positions],
            proteomes,
        )
        sequence_cols = [
            [float(v) for v in raw_features[key]]
            for key in keys
        ]
        structure = build_registered_structure_features(
            set(keys),
            registry_path=STRUCTURE_REGISTRY,
            registry_base=STRUCTURE_REGISTRY_BASE,
        )
        raw_struct = [
            [float(v) for v in structure[key][0]] for key in keys
        ]
        masks = [structure[key][1] for key in keys]
        scaled = transform_species_struct(
            raw_struct, masks, [species] * len(keys), species_scalers
        )
        features = BranchFeatures(
            sequence=sequence_cols,
            esm=[[0.0] for _ in keys],
            structure=scaled,
            structure_mask=masks,
            study_ids=None,
        )
        output = score_structure_ranker_bundle(bundle, features, device_name="cpu")
        scores = {
            p: float(score) for p, score in zip(positions, output.scores, strict=True)
        }
        metrics = within_protein_metrics(
            positions=positions,
            scores=scores,
            true_positions=bucket["positions"],
        )
        record = {
            "genes": sorted(bucket["genes"]),
            "accession": accession,
            "species": species,
            "n_cys": len(positions),
            "n_with_structure": sum(masks),
            "verified_positions": sorted(bucket["positions"]),
            "ranks": metrics["true_site_ranks"],
            "first_hit_burden": metrics["first_hit_burden"],
            "random_baseline_burden": metrics["random_baseline_burden"],
        }
        rows[bucket["oxiptm"]].append(record)
        print(
            f"[{bucket['oxiptm']}] {','.join(sorted(bucket['genes']))} "
            f"({accession}, n={len(positions)}, struct={sum(masks)}): "
            f"ranks={metrics['true_site_ranks']} "
            f"burden={metrics['first_hit_burden']}"
        )

    persulf_ranks = [
        rank
        for record in rows[PERSULFIDATION]
        for rank in record["ranks"].values()
    ]
    sno_ranks = [
        rank
        for record in rows[S_NITROSYLATION]
        for rank in record["ranks"].values()
    ]
    observed_diff = float(np.mean(persulf_ranks) - np.mean(sno_ranks))
    rng = np.random.RandomState(20260816)
    pooled = persulf_ranks + sno_ranks
    n_p = len(persulf_ranks)
    null_diffs = [
        float(np.mean(perm[:n_p]) - np.mean(perm[n_p:]))
        for perm in (rng.permutation(pooled) for _ in range(999))
    ]
    p_perm = (sum(1 for d in null_diffs if d <= observed_diff) + 1) / 1000

    document = {
        "track": "oxiptm_sentinel_v2",
        "claim_class": "diagnostic_only",
        "note": (
            "Same verified oxiPTM site list and within-protein comparison as "
            "oxiptm_sites_v1, but scored with the v2 per-species-scaled bundle "
            "(structure branch ACTIVE where the v3 registry provides a model; "
            "v1 masked everything)."
        ),
        "bundle": {
            "path": str(BUNDLE_PATH),
            "seed": bundle.seed,
            "sha256": None,
        },
        "persulfidation": persulf_ranks,
        "s_nitrosylation": sno_ranks,
        "persulfidation_mean_rank": round(float(np.mean(persulf_ranks)), 3),
        "s_nitrosylation_mean_rank": round(float(np.mean(sno_ranks)), 3),
        "observed_mean_diff_persulf_minus_sno": round(observed_diff, 3),
        "permutation_p_lower": round(p_perm, 4),
        "per_protein": rows,
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(
        json.dumps(document, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(f"\npersulfidation mean rank {document['persulfidation_mean_rank']} vs "
          f"SNO {document['s_nitrosylation_mean_rank']} (diff {observed_diff:+.3f}, "
          f"p={document['permutation_p_lower']})")
    print(f"wrote {OUTPUT}")


if __name__ == "__main__":
    main()
