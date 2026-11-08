"""Monte Carlo propagation, sensitivity analysis and rank stability."""

from __future__ import annotations

import math

import numpy as np
import pytest

from ihd_solutions.evidence import CreditMode
from ihd_solutions.uncertainty import (
    AdherencePrior,
    MonteCarloResult,
    lognormal_from_ci,
    monte_carlo_regimen,
    one_at_a_time,
    rank_stability,
    spearman_sensitivity,
    tornado,
    truncated_normal,
)


def test_truncated_normal_respects_the_reported_interval():
    rng = np.random.default_rng(1)
    draws = truncated_normal(rng, 0.5, 0.4, 0.6, 5000)
    assert draws.min() >= 0.4
    assert draws.max() <= 0.6
    assert draws.mean() == pytest.approx(0.5, abs=0.01)


def test_degenerate_interval_returns_a_constant():
    rng = np.random.default_rng(1)
    assert np.allclose(truncated_normal(rng, 0.5, 0.5, 0.5, 10), 0.5)
    assert np.allclose(lognormal_from_ci(rng, 0.8, 0.8, 0.8, 10), 0.8)


def test_lognormal_draws_reproduce_the_confidence_interval():
    rng = np.random.default_rng(2)
    draws = lognormal_from_ci(rng, 0.85, 0.78, 0.93, 200000)
    assert np.quantile(draws, 0.025) == pytest.approx(0.78, rel=0.02)
    assert np.quantile(draws, 0.975) == pytest.approx(0.93, rel=0.02)
    assert math.exp(np.mean(np.log(draws))) == pytest.approx(0.85, rel=0.01)


def test_draw_helpers_validate_arguments():
    rng = np.random.default_rng(0)
    with pytest.raises(ValueError):
        truncated_normal(rng, 0.5, 0.6, 0.7, 10)
    with pytest.raises(ValueError):
        truncated_normal(rng, 0.5, 0.4, 0.6, 0)
    with pytest.raises(ValueError):
        lognormal_from_ci(rng, 0.85, 0.9, 1.0, 10)
    with pytest.raises(ValueError):
        lognormal_from_ci(rng, 0.85, 0.78, 0.93, 0)


def test_adherence_prior_alpha_beta_are_consistent():
    prior = AdherencePrior(mean=0.8, concentration=10.0)
    assert prior.alpha == pytest.approx(8.0)
    assert prior.beta == pytest.approx(2.0)
    assert prior.alpha / (prior.alpha + prior.beta) == pytest.approx(0.8)


@pytest.mark.parametrize("kwargs", [{"mean": 0.0}, {"mean": 1.0}, {"concentration": 0.0}])
def test_adherence_prior_validation(kwargs):
    with pytest.raises(ValueError):
        AdherencePrior(**kwargs)


def test_monte_carlo_is_reproducible_for_a_fixed_seed(
    catalog, population, anchors, with_background, short_days
):
    regimen = with_background("evolocumab")
    kwargs = dict(draws=25, seed=123, mode=CreditMode.STRICT)
    first = monte_carlo_regimen(regimen, catalog, population, anchors, short_days, **kwargs)
    second = monte_carlo_regimen(regimen, catalog, population, anchors, short_days, **kwargs)
    assert np.array_equal(first.events_avoided_per_1000, second.events_avoided_per_1000)


def test_monte_carlo_varies_with_the_seed(
    catalog, population, anchors, with_background, short_days
):
    regimen = with_background("evolocumab")
    first = monte_carlo_regimen(regimen, catalog, population, anchors, short_days, draws=25, seed=1)
    second = monte_carlo_regimen(
        regimen, catalog, population, anchors, short_days, draws=25, seed=2
    )
    assert not np.array_equal(first.events_avoided_per_1000, second.events_avoided_per_1000)


def test_monte_carlo_summary_is_internally_consistent(
    catalog, population, anchors, with_background, short_days
):
    result = monte_carlo_regimen(
        with_background("evolocumab"), catalog, population, anchors, short_days, draws=40, seed=5
    )
    summary = result.summary()
    assert summary["draws"] == 40
    assert summary["events_avoided_p5"] <= summary["events_avoided_median"]
    assert summary["events_avoided_median"] <= summary["events_avoided_p95"]
    assert 0.0 <= summary["probability_of_benefit"] <= 1.0
    assert summary["probability_of_benefit"] + summary["probability_of_harm"] <= 1.0
    assert summary["ldl_median_mg_dl"] > 0


def test_neutral_therapy_has_a_benefit_probability_near_one_half(
    catalog, populations, anchors, with_background, short_days
):
    result = monte_carlo_regimen(
        with_background("ziltivekimab"),
        catalog,
        populations["ascvd_ckd_inflammation"],
        anchors,
        short_days,
        draws=60,
        seed=17,
    )
    extra_benefit = result.events_avoided_per_1000
    reference = monte_carlo_regimen(
        with_background(),
        catalog,
        populations["ascvd_ckd_inflammation"],
        anchors,
        short_days,
        draws=60,
        seed=17,
    )
    difference = extra_benefit - reference.events_avoided_per_1000
    assert abs(float(np.median(difference))) < 5.0


def test_monte_carlo_rejects_zero_draws(catalog, population, anchors, with_background, short_days):
    with pytest.raises(ValueError):
        monte_carlo_regimen(with_background(), catalog, population, anchors, short_days, draws=0)


def test_one_at_a_time_reports_direction_and_span():
    row = one_at_a_time(lambda x: 10.0 * x, base=1.0, low=0.5, high=1.5, name="adherence")
    assert row["base_output"] == pytest.approx(10.0)
    assert row["span"] == pytest.approx(10.0)
    assert row["direction"] == "increasing"
    decreasing = one_at_a_time(lambda x: -x, 0.0, -1.0, 1.0)
    assert decreasing["direction"] == "decreasing"


def test_one_at_a_time_requires_bracketing_interval():
    with pytest.raises(ValueError, match="bracket"):
        one_at_a_time(lambda x: x, base=2.0, low=0.0, high=1.0)


def test_tornado_is_sorted_by_span():
    frame = tornado(
        lambda x: x,
        {"small": (1.0, 0.9, 1.1), "large": (1.0, 0.0, 2.0), "medium": (1.0, 0.5, 1.5)},
    )
    assert frame["parameter"].tolist() == ["large", "medium", "small"]
    assert frame["span"].is_monotonic_decreasing


def test_spearman_sensitivity_identifies_the_dominant_parameter():
    rng = np.random.default_rng(3)
    dominant = rng.normal(size=2000)
    weak = rng.normal(size=2000)
    noise = rng.normal(size=2000)
    output = 5.0 * dominant + 0.2 * weak + 0.01 * noise
    frame = spearman_sensitivity({"dominant": dominant, "weak": weak, "noise": noise}, output)
    assert frame["parameter"].iloc[0] == "dominant"
    assert frame["spearman_rho"].iloc[0] > 0.9


def test_spearman_handles_constant_inputs_and_shape_errors():
    output = np.arange(10.0)
    frame = spearman_sensitivity({"constant": np.ones(10)}, output)
    assert frame["spearman_rho"].iloc[0] == 0.0
    flat = spearman_sensitivity({"x": output}, np.ones(10))
    assert flat["spearman_rho"].iloc[0] == 0.0
    with pytest.raises(ValueError, match="sample shape"):
        spearman_sensitivity({"x": np.ones(3)}, output)


def make_result(label: str, values: list[float]) -> MonteCarloResult:
    return MonteCarloResult(
        regimen_label=label,
        population_id="pop",
        mode=CreditMode.STRICT,
        events_avoided_per_1000=np.array(values, dtype=float),
        time_averaged_ldl_mg_dl=np.full(len(values), 50.0),
    )


def test_rank_stability_probabilities_are_bounded():
    results = [
        make_result("always_best", [30.0, 31.0, 32.0, 33.0]),
        make_result("middle", [20.0, 21.0, 19.0, 22.0]),
        make_result("worst", [1.0, 2.0, 3.0, 0.0]),
    ]
    frame = rank_stability(results, top_n=1)
    assert frame["regimen"].iloc[0] == "always_best"
    assert frame["probability_top_rank"].iloc[0] == pytest.approx(1.0)
    assert frame["probability_top_rank"].between(0.0, 1.0).all()


def test_rank_stability_handles_more_ranks_than_regimens():
    frame = rank_stability([make_result("a", [1.0, 2.0]), make_result("b", [2.0, 1.0])], top_n=5)
    assert frame["probability_top_rank"].eq(1.0).all()


def test_rank_stability_validates_inputs():
    with pytest.raises(ValueError, match="no Monte Carlo results"):
        rank_stability([])
    with pytest.raises(ValueError):
        rank_stability([make_result("a", [1.0])], top_n=0)
    with pytest.raises(ValueError, match="same number of draws"):
        rank_stability([make_result("a", [1.0, 2.0]), make_result("b", [1.0])])
