"""Unit conversions, guards and integration helpers."""

from __future__ import annotations

import numpy as np
import pytest

from ihd_solutions import units


def test_cholesterol_conversion_round_trip():
    values = np.array([10.0, 55.0, 70.0, 130.0, 190.0])
    assert np.allclose(units.mmoll_to_mgdl(units.mgdl_to_mmoll(values)), values)


def test_cholesterol_conversion_scalar_matches_known_value():
    assert units.mgdl_to_mmoll(38.67) == pytest.approx(1.0)
    assert units.mmoll_to_mgdl(1.0) == pytest.approx(38.67)


@pytest.mark.parametrize("value", [0.0, -1.0, float("nan"), float("inf")])
def test_require_positive_rejects_bad_values(value):
    with pytest.raises(ValueError):
        units.require_positive("x", value)


@pytest.mark.parametrize("value", [-0.1, 1.0, 1.5, float("nan")])
def test_require_fraction_exclusive_upper_bound(value):
    with pytest.raises(ValueError):
        units.require_fraction("x", value)


def test_require_fraction_can_allow_one():
    assert units.require_fraction("x", 1.0, allow_one=True) == 1.0


def test_require_nonnegative_accepts_zero_and_rejects_negative():
    assert units.require_nonnegative("x", 0.0) == 0.0
    with pytest.raises(ValueError):
        units.require_nonnegative("x", -1e-9)


def test_rate_probability_round_trip():
    for rate in (0.001, 0.03, 0.08, 0.5):
        probability = units.annual_rate_to_probability(rate, 5.0)
        assert 0.0 <= probability < 1.0
        assert units.probability_to_annual_rate(probability, 5.0) == pytest.approx(rate)


def test_zero_rate_gives_zero_probability():
    assert units.annual_rate_to_probability(0.0, 10.0) == 0.0


def test_probability_of_one_is_rejected_by_inverse():
    with pytest.raises(ValueError):
        units.probability_to_annual_rate(1.0)


def test_trapz_cumulative_matches_analytic_integral():
    x = np.linspace(0.0, 4.0, 401)
    cumulative = units.trapz_cumulative(x, x)
    assert cumulative[0] == 0.0
    assert np.allclose(cumulative, 0.5 * x**2, atol=1e-4)


def test_trapz_cumulative_of_constant_is_exact():
    x = np.linspace(0.0, 5.0, 51)
    cumulative = units.trapz_cumulative(np.full_like(x, 0.04), x)
    assert np.allclose(cumulative, 0.04 * x)


def test_trapz_rejects_shape_mismatch_and_nonmonotonic_axis():
    with pytest.raises(ValueError, match="shape mismatch"):
        units.trapz_cumulative(np.ones(3), np.ones(4))
    with pytest.raises(ValueError, match="strictly increasing"):
        units.trapz_cumulative(np.ones(3), np.array([0.0, 1.0, 1.0]))
    with pytest.raises(ValueError, match="one-dimensional"):
        units.trapz_cumulative(np.ones((2, 2)), np.ones((2, 2)))
    with pytest.raises(ValueError, match="non-finite"):
        units.trapz_cumulative(np.array([1.0, np.nan]), np.array([0.0, 1.0]))


def test_trapz_handles_empty_input():
    assert units.trapz_cumulative(np.empty(0), np.empty(0)).size == 0


def test_time_grid_covers_horizon_inclusively():
    grid = units.time_grid(5.0, 30.0)
    assert grid[0] == 0.0
    assert grid[-1] == pytest.approx(5.0 * units.DAYS_PER_YEAR)
    assert np.all(np.diff(grid) > 0)


@pytest.mark.parametrize("horizon,step", [(0.0, 30.0), (5.0, 0.0), (-1.0, 30.0)])
def test_time_grid_rejects_invalid_arguments(horizon, step):
    with pytest.raises(ValueError):
        units.time_grid(horizon, step)


def test_check_all_finite_accepts_lists_and_rejects_nan():
    assert units.check_all_finite("x", [1.0, 2.0]).shape == (2,)
    with pytest.raises(ValueError, match="non-finite"):
        units.check_all_finite("x", [1.0, float("inf")])
