"""Dosing kinetics: how much of a therapy's effect is present at a given time."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ihd_solutions.catalog import DosingProfile, Therapy
from ihd_solutions.units import require_fraction, require_nonnegative


@dataclass(frozen=True)
class Adherence:
    """Fraction of prescribed effect taken, and how persistence decays over time."""

    taking_fraction: float = 1.0
    discontinuation_half_life_days: float = 0.0
    start_day: float = 0.0

    def __post_init__(self) -> None:
        require_fraction("taking_fraction", self.taking_fraction, allow_one=True)
        require_nonnegative("discontinuation_half_life_days", self.discontinuation_half_life_days)
        require_nonnegative("start_day", self.start_day)

    def persistence(self, days: np.ndarray) -> np.ndarray:
        days = np.asarray(days, dtype=float)
        elapsed = np.maximum(days - self.start_day, 0.0)
        if self.discontinuation_half_life_days == 0:
            return np.ones_like(elapsed)
        return np.exp(-np.log(2.0) * elapsed / self.discontinuation_half_life_days)


PERFECT_ADHERENCE = Adherence()


def ramp_fraction(days: np.ndarray, therapy: Therapy, start_day: float = 0.0) -> np.ndarray:
    """Linear onset from ``onset_days`` to ``plateau_days`` after the start day."""
    days = np.asarray(days, dtype=float)
    elapsed = days - start_day
    active = np.clip(elapsed - therapy.onset_days, 0.0, None)
    span = max(therapy.plateau_days - therapy.onset_days, 1e-9)
    ramp = np.clip(active / span, 0.0, 1.0)
    return np.where(elapsed < therapy.onset_days, 0.0, ramp)


def interval_modulation(days: np.ndarray, therapy: Therapy, start_day: float = 0.0) -> np.ndarray:
    """Peak-to-trough sawtooth for interval-dosed injectables."""
    if therapy.profile is not DosingProfile.INTERVAL_SC or therapy.interval_days <= 0:
        return np.ones_like(np.asarray(days, dtype=float))
    elapsed = np.maximum(np.asarray(days, dtype=float) - start_day, 0.0)
    phase = np.mod(elapsed, therapy.interval_days) / therapy.interval_days
    return 1.0 - (1.0 - therapy.trough_fraction) * phase


def waning_fraction(days: np.ndarray, therapy: Therapy, start_day: float = 0.0) -> np.ndarray:
    """Durability decay for one-time interventions (0 half-life means permanent)."""
    if therapy.waning_half_life_days <= 0:
        return np.ones_like(np.asarray(days, dtype=float))
    elapsed = np.maximum(np.asarray(days, dtype=float) - start_day - therapy.plateau_days, 0.0)
    return np.exp(-np.log(2.0) * elapsed / therapy.waning_half_life_days)


def effect_fraction(
    days: np.ndarray,
    therapy: Therapy,
    adherence: Adherence = PERFECT_ADHERENCE,
) -> np.ndarray:
    """Fraction of the therapy's maximal effect present on each day in ``days``."""
    days = np.asarray(days, dtype=float)
    if days.ndim != 1:
        raise ValueError("days must be one-dimensional")
    if not np.isfinite(days).all():
        raise ValueError("days must be finite")
    if np.any(np.diff(days) < 0):
        raise ValueError("days must be non-decreasing")
    values = (
        ramp_fraction(days, therapy, adherence.start_day)
        * interval_modulation(days, therapy, adherence.start_day)
        * waning_fraction(days, therapy, adherence.start_day)
        * adherence.persistence(days)
        * adherence.taking_fraction
    )
    values = np.where(days < adherence.start_day, 0.0, values)
    if not np.isfinite(values).all():
        raise ValueError(f"{therapy.therapy_id}: non-finite effect fraction")
    return np.clip(values, 0.0, 1.0)


def mean_effect(
    days: np.ndarray, therapy: Therapy, adherence: Adherence = PERFECT_ADHERENCE
) -> float:
    """Time-weighted mean effect fraction across the grid."""
    values = effect_fraction(days, therapy, adherence)
    if days.size < 2:
        return float(values[0]) if values.size else 0.0
    weights = np.diff(days)
    total = float(weights.sum())
    if total <= 0:
        raise ValueError("time grid must span a positive duration")
    midpoints = 0.5 * (values[1:] + values[:-1])
    return float(np.dot(midpoints, weights) / total)
