"""Target-label-free cross-crop scoring with a fail-closed source firewall."""

from __future__ import annotations

import json
import statistics
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from plantpersulf.evaluation.metrics import average_precision
from plantpersulf.evaluation.rank_ensemble import ApplicabilityEnvelope
from plantpersulf.features.site_biology import BIOLOGY_FEATURE_NAMES
from plantpersulf.models.additive_pu_ranker import (
    AdditivePuConfig,
    fit_additive_pu_ranker,
)
from plantpersulf.models.traditional import pu_logistic_regression_scores
from plantpersulf.provenance.hashing import hash_file
from plantpersulf.workflows.tomato_ranker_v2 import (
    _lower_quantile,
    _percentiles,
    sensitivity_priors,
)


@dataclass(frozen=True)
class SourceBatch:
    species: str
    study_accession: str
    feature_names: tuple[str, ...]
    features: tuple[tuple[float, ...], ...]
    labels: tuple[str, ...]
    protein_ids: tuple[str, ...]
    source_sha256: str


@dataclass(frozen=True)
class TargetCandidateBatch:
    site_keys: tuple[str, ...]
    feature_names: tuple[str, ...]
    features: tuple[tuple[float, ...], ...]
    source_sha256: str


@dataclass(frozen=True)
class SourceTransferDecision:
    admitted: bool
    reason: str


def validate_target_label_free_sources(
    sources: tuple[SourceBatch, ...], target_species: str
) -> None:
    if not sources:
        raise RuntimeError("at least one registered source batch is required")
    keys = [(source.species, source.study_accession) for source in sources]
    if len(keys) != len(set(keys)):
        raise RuntimeError("duplicate source species-study batch")
    for source in sources:
        if source.species == target_species:
            raise RuntimeError("target label leakage: target species in source batches")
        if not source.source_sha256 or not source.study_accession:
            raise RuntimeError("source batch lacks registered accession or SHA256")
        if source.feature_names != BIOLOGY_FEATURE_NAMES:
            raise RuntimeError("source feature order differs from biological contract")
        if set(source.labels) - {"positive", "unlabeled"}:
            raise RuntimeError("source batch contains a forbidden label")
        if not source.features or not (
            len(source.features) == len(source.labels) == len(source.protein_ids)
        ):
            raise RuntimeError("source batch row count mismatch")
        if any(len(row) != len(source.feature_names) for row in source.features):
            raise RuntimeError("source feature width mismatch")


def source_transfer_gate(sources: tuple[SourceBatch, ...]) -> SourceTransferDecision:
    if len({source.species for source in sources}) < 2:
        return SourceTransferDecision(False, "fewer_than_two_source_species")
    for species in {source.species for source in sources}:
        if len({s.study_accession for s in sources if s.species == species}) < 2:
            return SourceTransferDecision(
                False, f"single_study_source_species:{species}"
            )
    species_holdouts = [
        (f"species:{species}", tuple(s for s in sources if s.species == species))
        for species in sorted({source.species for source in sources})
    ]
    study_holdouts = [(f"study:{s.study_accession}", (s,)) for s in sources]
    holdouts = [*species_holdouts, *study_holdouts]
    for name, held_out in holdouts:
        training = tuple(source for source in sources if source not in held_out)
        train_x = [list(row) for batch in training for row in batch.features]
        train_y = [label for batch in training for label in batch.labels]
        train_proteins = [
            f"{batch.species}|{batch.study_accession}|{protein}"
            for batch in training
            for protein in batch.protein_ids
        ]
        test_x = [list(row) for batch in held_out for row in batch.features]
        test_y = [label for batch in held_out for label in batch.labels]
        priors = sensitivity_priors(
            train_y.count("positive"), len(train_y), (1.0, 1.5, 2.0), 0.5
        )
        candidate, baseline = [], []
        for seed in range(5):
            baseline_value = average_precision(
                list(
                    zip(
                        pu_logistic_regression_scores(train_x, train_y, test_x, seed),
                        test_y,
                        strict=True,
                    )
                )
            )
            if baseline_value is None:
                raise RuntimeError("held-out source batch has no positive")
            baseline.append(baseline_value)
            for prior in priors:
                model = fit_additive_pu_ranker(
                    train_x,
                    train_y,
                    train_proteins,
                    training[0].feature_names,
                    AdditivePuConfig(class_prior=prior, seed=seed),
                )
                value = average_precision(
                    list(zip(model.score(test_x), test_y, strict=True))
                )
                if value is None:
                    raise RuntimeError("held-out source batch has no positive")
                candidate.append(value)
        if min(candidate) <= max(baseline):
            return SourceTransferDecision(
                False,
                f"unstable_source_transfer:{name}:candidate_median={statistics.median(candidate):.6f}",
            )
    return SourceTransferDecision(True, "all_source_holdouts_stably_better")


def run_cross_crop_target_label_free(
    config_path: Path,
    sources: tuple[SourceBatch, ...],
    target_candidates: TargetCandidateBatch,
    output_dir: Path,
) -> dict[str, object]:
    cfg: dict[str, Any] = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if cfg.get("claim_class") != "target_label_free_transfer_not_gate2":
        raise RuntimeError("cross-crop workflow cannot emit a Gate 2 claim")
    target_species = str(cfg["target_species"])
    evaluation = dict(cfg["evaluation"])
    wetlab = dict(cfg["wetlab"])
    if (
        evaluation["target_labels_visible_during_fit"] is not False
        or evaluation["target_labels_visible_during_tuning"] is not False
        or wetlab["enabled"] is not False
        or float(wetlab["max_noncontrol_fraction"]) > 0.25
    ):
        raise RuntimeError("target-label or wet-lab policy differs from frozen config")
    validate_target_label_free_sources(sources, target_species)
    if {source.species for source in sources} != set(cfg["source_species"]):
        raise RuntimeError("source species differ from frozen config")
    if output_dir.exists():
        raise FileExistsError(output_dir)
    if (
        not target_candidates.source_sha256
        or target_candidates.feature_names != BIOLOGY_FEATURE_NAMES
    ):
        raise RuntimeError(
            "target candidate batch lacks registered biological contract"
        )
    if (
        not target_candidates.features
        or len(target_candidates.site_keys) != len(target_candidates.features)
        or len(set(target_candidates.site_keys)) != len(target_candidates.site_keys)
    ):
        raise RuntimeError("target candidate row count or keys are invalid")
    if any(
        len(row) != len(BIOLOGY_FEATURE_NAMES) for row in target_candidates.features
    ):
        raise RuntimeError("target feature width mismatch")
    decision = source_transfer_gate(sources)
    train_x = [list(row) for batch in sources for row in batch.features]
    train_y = [label for batch in sources for label in batch.labels]
    proteins = [
        f"{batch.species}|{batch.study_accession}|{protein}"
        for batch in sources
        for protein in batch.protein_ids
    ]
    target_x = [list(row) for row in target_candidates.features]
    priors = sensitivity_priors(
        train_y.count("positive"), len(train_y), (1.0, 1.5, 2.0), 0.5
    )
    columns = [
        _percentiles(
            fit_additive_pu_ranker(
                train_x,
                train_y,
                proteins,
                BIOLOGY_FEATURE_NAMES,
                AdditivePuConfig(class_prior=prior, seed=seed),
            ).score(target_x)
        )
        for seed in range(5)
        for prior in priors
    ]
    envelope = ApplicabilityEnvelope.fit(train_x)
    in_domain = [envelope.contains(row) for row in target_x]
    fraction = sum(in_domain) / len(in_domain)
    scores = {
        key: {
            "median_percentile": statistics.median(column[i] for column in columns),
            "conservative_percentile": _lower_quantile(
                [column[i] for column in columns], 0.2
            ),
            "in_domain": in_domain[i],
        }
        for i, key in enumerate(target_candidates.site_keys)
    }
    output_dir.mkdir(parents=True, exist_ok=False)
    (output_dir / "scores.json").write_text(
        json.dumps({"scores": scores}, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    result: dict[str, object] = {
        "complete": True,
        "claim_class": cfg["claim_class"],
        "target_species": target_species,
        "sources": [
            {
                "species": s.species,
                "study_accession": s.study_accession,
                "sha256": s.source_sha256,
            }
            for s in sources
        ],
        "target_candidate_sha256": target_candidates.source_sha256,
        "target_label_hash": None,
        "source_gate_admitted": decision.admitted,
        "source_gate_reason": decision.reason,
        "target_in_domain_fraction": fraction,
        "applicability_passed": fraction
        >= float(dict(cfg["applicability"])["min_target_in_domain_fraction"]),
        "wetlab_eligible": False,
        "scores_sha256": hash_file(output_dir / "scores.json", "sha256"),
    }
    result["wetlab_eligible"] = bool(
        result["source_gate_admitted"] and result["applicability_passed"]
    )
    (output_dir / "manifest.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return result
