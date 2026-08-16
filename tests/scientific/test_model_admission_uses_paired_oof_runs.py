"""Model admission must be conservative over paired out-of-fold runs."""

from plantpersulf.evaluation.model_admission import admit_candidate_model


def test_candidate_must_improve_every_required_metric() -> None:
    candidate = {
        f"r{index}f0": (0.3, 1.5, 0.1 if index == 0 else 0.3)
        for index in range(5)
    }
    baseline = {f"r{index}f0": (0.1, 1.2, 0.2) for index in range(5)}
    blocks = {f"r{index}f0": f"r{index}" for index in range(5)}

    decision = admit_candidate_model(candidate, baseline, blocks, n_boot=500)

    assert decision.admitted is False
    assert "mrr" in decision.reason
