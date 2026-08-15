"""Deterministic two-model stacking helpers for the accuracy campaign.

W2 (Sul-BertGRU stacking evaluation) needs two combination schemes:

- ``rank_aggregate_scores`` — training-free mean/median rank aggregation.
  Lower rank is better; the combined score is the negated mean rank so a
  higher combined score means "ranked higher by both models on average".
  With two models the mean and median coincide.
- ``fit_score_combiner`` — a deterministic grid search over the convex
  combination ``alpha * scores_a + (1 - alpha) * scores_b``, selecting the
  ``alpha`` that maximizes average precision (with ``label == "positive"``)
  on the supplied rows. Callers pass validation-partition rows so the choice
  never sees test rows.

Both helpers are pure and deterministic; they never touch scientific inputs
directly — every score and label arrives as a plain aligned sequence.
"""

from __future__ import annotations

from plantpersulf.evaluation.rank_ensemble import summarize_rank_runs


def format_site_key(global_protein_id: str, cys_position: int) -> str:
    """Render a comparison site key as the ``"protein|position"`` string used
    by the frozen-test and candidate-table ``site_key`` columns."""
    return f"{global_protein_id}|{cys_position}"


def rank_aggregate_scores(
    site_keys: tuple[str, ...],
    model_scores: tuple[tuple[float, ...], ...],
) -> tuple[float, ...]:
    """Combine per-model score vectors by mean rank across models.

    ``model_scores`` holds one score tuple per model, each aligned with
    ``site_keys``. Ranks are computed per model with the same stable
    tie-break as ``rank_ensemble.summarize_rank_runs`` (rank 1 = best).
    The returned combined score is ``-mean_rank`` so higher is better.
    """
    if not site_keys:
        raise ValueError("site_keys must not be empty")
    if len(model_scores) < 2:
        raise ValueError("rank aggregation requires at least two models")
    if any(len(scores) != len(site_keys) for scores in model_scores):
        raise ValueError("every model score vector must align with site_keys")
    summaries = summarize_rank_runs(site_keys, model_scores)
    return tuple(-summary.median_rank for summary in summaries)


def fit_score_combiner(
    scores_a: tuple[float, ...],
    scores_b: tuple[float, ...],
    labels: tuple[str, ...],
    *,
    n_grid: int = 21,
) -> float:
    """Grid-search the convex combination weight of model A on the given rows.

    Returns the ``alpha`` in ``[0, 1]`` maximizing average precision of
    ``alpha * scores_a + (1 - alpha) * scores_b`` against
    ``label == "positive"``. Ties keep the smallest alpha. Requires at least
    one positive and one non-positive row.
    """
    if not scores_a or not scores_b or not labels:
        raise ValueError("scores and labels must not be empty")
    if not (len(scores_a) == len(scores_b) == len(labels)):
        raise ValueError("scores and labels must align")
    if n_grid < 2:
        raise ValueError("n_grid must be at least 2")
    if "positive" not in labels:
        raise ValueError("combiner fitting requires at least one positive row")
    if all(label == "positive" for label in labels):
        raise ValueError("combiner fitting requires at least one non-positive row")
    best_alpha = 0.0
    best_ap = -1.0
    for step in range(n_grid):
        alpha = step / (n_grid - 1)
        combined = [
            alpha * a + (1.0 - alpha) * b
            for a, b in zip(scores_a, scores_b, strict=True)
        ]
        ap = _average_precision(combined, labels)
        if ap > best_ap:
            best_ap = ap
            best_alpha = alpha
    return best_alpha


def align_site_scores(
    primary: dict[str, float], secondary: dict[str, float]
) -> tuple[tuple[str, ...], tuple[float, ...], tuple[float, ...]]:
    """Align two per-site score maps on exactly the same site keys.

    Returns ``(site_keys, primary_scores, secondary_scores)`` sorted by site
    key. Raises if either map covers keys the other lacks — silent dropping
    of mismatched panel rows would corrupt any stacking evaluation, so the
    mismatch is an error, not a warning.
    """
    if set(primary) != set(secondary):
        raise ValueError(
            "score tables do not cover exactly the same site keys "
            f"(primary-only={len(set(primary) - set(secondary))}, "
            f"secondary-only={len(set(secondary) - set(primary))})"
        )
    keys = tuple(sorted(primary))
    return (
        keys,
        tuple(primary[key] for key in keys),
        tuple(secondary[key] for key in keys),
    )


def _average_precision(scores: list[float], labels: tuple[str, ...]) -> float:
    """sklearn-compatible average precision on ``label == "positive"``."""
    import sklearn.metrics  # type: ignore[import-untyped]  # deferred: optional models dependency

    return float(
        sklearn.metrics.average_precision_score(
            [label == "positive" for label in labels], scores
        )
    )
