"""RED (Task 8): hyperparameter/model selection must never see test data.

``select_hyperparameters`` accepts only train and validation data — its
signature has no test parameter at all, so it is structurally impossible to
select hyperparameters using test labels. ``TestEvaluationGuard`` wraps the
test set and raises if evaluated more than once, enforcing "test only runs
once" (Codex Task 8 TDD focus).

Expected RED: ``plantpersulf.models.traditional`` does not exist yet.
"""

from __future__ import annotations

import inspect

import pytest

from plantpersulf.models.traditional import (  # RED: module missing
    TestEvaluationGuard,
    select_hyperparameters,
)


def test_select_hyperparameters_signature_has_no_test_parameter() -> None:
    signature = inspect.signature(select_hyperparameters)
    parameter_names = set(signature.parameters)
    assert not any("test" in name.lower() for name in parameter_names)
    assert {"candidates", "train_X", "train_y", "val_X", "val_y", "score_fn"} <= (
        parameter_names
    )


def test_select_hyperparameters_picks_best_validation_score() -> None:
    candidates = [{"c": 0.1}, {"c": 1.0}, {"c": 10.0}]

    def score_fn(params: dict, train_X, train_y, val_X, val_y) -> float:  # type: ignore[no-untyped-def]
        # Deterministic stand-in: score is just the candidate's "c" value,
        # so the best candidate is the one with the largest c.
        return float(params["c"])

    best = select_hyperparameters(
        candidates=candidates,
        train_X=[[0.0]],
        train_y=["positive"],
        val_X=[[0.0]],
        val_y=["positive"],
        score_fn=score_fn,
    )
    assert best == {"c": 10.0}


def test_evaluation_guard_allows_exactly_one_test_evaluation() -> None:
    guard = TestEvaluationGuard(test_X=[[0.0]], test_y=["positive"])

    def score_fn(test_X, test_y) -> float:  # type: ignore[no-untyped-def]
        return 1.0

    result = guard.evaluate(score_fn)
    assert result == 1.0


def test_evaluation_guard_raises_on_second_evaluation() -> None:
    guard = TestEvaluationGuard(test_X=[[0.0]], test_y=["positive"])

    def score_fn(test_X, test_y) -> float:  # type: ignore[no-untyped-def]
        return 1.0

    guard.evaluate(score_fn)
    with pytest.raises(RuntimeError, match="once"):
        guard.evaluate(score_fn)
