"""Regenerative scenarios: small-sample uncertainty and zero event credit."""

from __future__ import annotations

import numpy as np
import pytest

from ihd_solutions.evidence import CreditMode
from ihd_solutions.regeneration import (
    DEFAULT_REGENERATION_ANCHOR,
    RegenerationAnchor,
    anchor_posteriors,
    beta_binomial_posterior,
    event_credit,
    simulate_regeneration,
    surveillance_requirements,
)


def test_uniform_prior_posterior_mean_is_correct():
    posterior = beta_binomial_posterior("lvef", 2, 8)
    assert posterior.alpha == 3.0
    assert posterior.beta == 7.0
    assert posterior.mean == pytest.approx(3.0 / 10.0)
    assert posterior.observed_fraction == pytest.approx(0.25)


def test_posterior_interval_brackets_the_mean_and_is_wide():
    posterior = beta_binomial_posterior("lvef", 2, 8)
    low, high = posterior.interval()
    assert 0.0 <= low < posterior.mean < high <= 1.0
    assert high - low > 0.4


def test_all_responders_still_leaves_substantial_uncertainty():
    posterior = beta_binomial_posterior("nyha", 8, 8)
    low, high = posterior.interval()
    assert low < 0.95
    assert posterior.mean == pytest.approx(9.0 / 10.0)


def test_more_responders_shifts_the_posterior_upwards():
    means = [beta_binomial_posterior("x", k, 8).mean for k in (1, 2, 4, 8)]
    assert all(later > earlier for earlier, later in zip(means, means[1:], strict=False))


def test_predictive_counts_are_bounded_and_reproducible():
    posterior = beta_binomial_posterior("lvef", 2, 8)
    first = posterior.predictive_responders(75, 500, seed=7)
    second = posterior.predictive_responders(75, 500, seed=7)
    assert np.array_equal(first, second)
    assert first.min() >= 0
    assert first.max() <= 75
    assert not np.array_equal(first, posterior.predictive_responders(75, 500, seed=8))


@pytest.mark.parametrize("mode", list(CreditMode))
def test_no_mode_credits_events_for_cell_sheets(mode):
    credit, reason = event_credit(mode)
    assert credit == 0.0
    assert "eight uncontrolled patients" in reason


def test_surveillance_requirements_are_explicit():
    requirements = surveillance_requirements()
    assert requirements["approval_type"] == "conditional_time_limited"
    assert requirements["confirmation_target_cases"] == 75
    assert requirements["remuscularization_proven"] is False
    assert requirements["mechanism_emphasised"] == "paracrine"
    assert len(requirements["required_monitoring"]) >= 4


def test_simulation_table_reports_zero_credit_and_wide_intervals():
    frame = simulate_regeneration(draws=400, seed=11)
    assert len(frame) == 3
    assert (frame["mace_credit"] == 0.0).all()
    assert (frame["mortality_credit"] == 0.0).all()
    assert not frame["remuscularization_proven"].any()
    assert (frame["posterior_p95"] > frame["posterior_p5"]).all()
    assert (frame["predicted_responders_p95"] >= frame["predicted_responders_median"]).all()
    assert (frame["predicted_responders_median"] >= frame["predicted_responders_p5"]).all()


def test_anchor_posteriors_cover_the_three_functional_endpoints():
    posteriors = anchor_posteriors()
    assert set(posteriors) == {"lvef", "peak_vo2", "nyha"}
    assert posteriors["peak_vo2"].evaluable == 7


@pytest.mark.parametrize(
    "changes",
    [
        {"n_treated": 0},
        {"lvef_evaluable": 9},
        {"lvef_responders": 9},
        {"vo2_responders": -1},
        {"lvef_threshold_pct": 0.0},
        {"vo2_threshold_pct": 0.0},
    ],
)
def test_anchor_validation(changes):
    with pytest.raises(ValueError):
        RegenerationAnchor(**changes)


def test_posterior_rejects_impossible_counts_and_priors():
    with pytest.raises(ValueError):
        beta_binomial_posterior("x", 3, 2)
    with pytest.raises(ValueError):
        beta_binomial_posterior("x", 1, 0)
    with pytest.raises(ValueError):
        beta_binomial_posterior("x", 1, 2, prior=(0.0, 1.0))


def test_interval_and_predictive_argument_validation():
    posterior = beta_binomial_posterior("x", 2, 8)
    with pytest.raises(ValueError):
        posterior.interval(0.0)
    with pytest.raises(ValueError):
        posterior.interval(1.5)
    with pytest.raises(ValueError):
        posterior.predictive_responders(0, 100, seed=1)
    with pytest.raises(ValueError):
        posterior.predictive_responders(10, 0, seed=1)


def test_default_anchor_matches_the_published_pivotal_counts():
    anchor = DEFAULT_REGENERATION_ANCHOR
    assert anchor.n_treated == 8
    assert (anchor.lvef_responders, anchor.lvef_evaluable) == (2, 8)
    assert (anchor.vo2_responders, anchor.vo2_evaluable) == (4, 7)
    assert (anchor.nyha_responders, anchor.nyha_evaluable) == (8, 8)


@pytest.mark.parametrize("kwargs", [{"cohort_size": 0}, {"draws": 0}])
def test_simulation_argument_validation(kwargs):
    with pytest.raises(ValueError):
        simulate_regeneration(**kwargs)
