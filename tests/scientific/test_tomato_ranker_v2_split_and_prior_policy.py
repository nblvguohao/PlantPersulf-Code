"""Scientific policy checks for frozen tomato-v2 training inputs."""

from plantpersulf.workflows.tomato_ranker_v2 import sensitivity_priors


def test_prior_grid_starts_at_observed_positive_fraction() -> None:
    assert sensitivity_priors(10, 100, (1.0, 1.5, 2.0), 0.5) == (0.1, 0.15, 0.2)
