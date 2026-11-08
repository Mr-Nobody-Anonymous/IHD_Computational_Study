"""Unit conversions and numeric guards shared by every simulation module."""

from __future__ import annotations

import math
from collections.abc import Iterable

import numpy as np

#: Cholesterol conversion factor: 1 mmol/L = 38.67 mg/dL.
MG_DL_PER_MMOL_L = 38.67
DAYS_PER_YEAR = 365.25
Z_95 = 1.959963984540054


def mgdl_to_mmoll(value: float | np.ndarray) -> float | np.ndarray:
    """Convert LDL cholesterol from mg/dL to mmol/L."""
    return (
        np.asarray(value, dtype=float) / MG_DL_PER_MMOL_L
        if isinstance(value, np.ndarray)
        else float(value) / MG_DL_PER_MMOL_L
    )


def mmoll_to_mgdl(value: float | np.ndarray) -> float | np.ndarray:
    """Convert LDL cholesterol from mmol/L to mg/dL."""
    return (
        np.asarray(value, dtype=float) * MG_DL_PER_MMOL_L
        if isinstance(value, np.ndarray)
        else float(value) * MG_DL_PER_MMOL_L
    )


def require_finite(name: str, value: float) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite, got {value!r}")
    return value


def require_positive(name: str, value: float) -> float:
    value = require_finite(name, value)
    if value <= 0:
        raise ValueError(f"{name} must be positive, got {value!r}")
    return value


def require_fraction(name: str, value: float, *, allow_one: bool = False) -> float:
    value = require_finite(name, value)
    upper_ok = value <= 1.0 if allow_one else value < 1.0
    if value < 0.0 or not upper_ok:
        limit = "1 inclusive" if allow_one else "1 exclusive"
        raise ValueError(f"{name} must lie in [0, {limit}], got {value!r}")
    return value


def require_nonnegative(name: str, value: float) -> float:
    value = require_finite(name, value)
    if value < 0:
        raise ValueError(f"{name} must be non-negative, got {value!r}")
    return value


def annual_rate_to_probability(rate: float, years: float = 1.0) -> float:
    """Convert a constant hazard into a cumulative-incidence probability."""
    rate = require_nonnegative("annual rate", rate)
    years = require_nonnegative("years", years)
    return float(-np.expm1(-rate * years))


def probability_to_annual_rate(probability: float, years: float = 1.0) -> float:
    """Inverse of :func:`annual_rate_to_probability`."""
    probability = require_fraction("probability", probability)
    years = require_positive("years", years)
    return float(-np.log1p(-probability) / years)


def trapz_cumulative(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Cumulative trapezoidal integral with a leading zero (shape preserved)."""
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    if y.shape != x.shape:
        raise ValueError(f"shape mismatch: {y.shape} versus {x.shape}")
    if y.ndim != 1:
        raise ValueError("trapz_cumulative expects one-dimensional input")
    if y.size == 0:
        return np.empty(0, dtype=float)
    if np.any(np.diff(x) <= 0):
        raise ValueError("x must be strictly increasing")
    if not (np.isfinite(y).all() and np.isfinite(x).all()):
        raise ValueError("non-finite values in integration input")
    segments = 0.5 * (y[1:] + y[:-1]) * np.diff(x)
    return np.concatenate(([0.0], np.cumsum(segments)))


def time_grid(horizon_years: float, step_days: float = 30.0) -> np.ndarray:
    """Uniform day grid from 0 to ``horizon_years`` inclusive."""
    horizon_years = require_positive("horizon_years", horizon_years)
    step_days = require_positive("step_days", step_days)
    total_days = horizon_years * DAYS_PER_YEAR
    n_steps = max(1, int(math.ceil(total_days / step_days)))
    grid = np.linspace(0.0, total_days, n_steps + 1)
    if grid.size < 2:
        raise ValueError("time grid must contain at least two points")
    return grid


def check_all_finite(name: str, values: Iterable[float] | np.ndarray) -> np.ndarray:
    array = np.asarray(list(values) if not isinstance(values, np.ndarray) else values, dtype=float)
    if not np.isfinite(array).all():
        raise ValueError(f"{name} contains non-finite values")
    return array
