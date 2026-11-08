"""Dosing kinetics, adherence and persistence."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from ihd_solutions.pharmacology import (
    Adherence,
    effect_fraction,
    interval_modulation,
    mean_effect,
    ramp_fraction,
    waning_fraction,
)


def test_effect_is_zero_before_start_and_bounded_afterwards(catalog):
    therapy = catalog["evolocumab"]
    days = np.arange(0.0, 200.0, 5.0)
    adherence = Adherence(start_day=60.0)
    values = effect_fraction(days, therapy, adherence)
    assert np.all(values[days < 60.0] == 0.0)
    assert values.min() >= 0.0
    assert values.max() <= 1.0


def test_oral_ramp_is_monotone_until_plateau(catalog):
    therapy = catalog["statin_high_intensity"]
    days = np.arange(0.0, 60.0)
    values = ramp_fraction(days, therapy)
    assert np.all(np.diff(values) >= 0)
    assert values[0] == 0.0
    assert values[-1] == pytest.approx(1.0)
    assert values[int(therapy.plateau_days)] == pytest.approx(1.0)


def test_interval_dosing_produces_a_sawtooth_between_peak_and_trough(catalog):
    therapy = catalog["evolocumab"]
    days = np.array([0.0, 7.0, 13.9, 14.0, 21.0])
    modulation = interval_modulation(days, therapy)
    assert modulation[0] == pytest.approx(1.0)
    assert modulation[3] == pytest.approx(1.0)
    assert modulation.min() >= therapy.trough_fraction
    assert modulation[1] == pytest.approx(1.0 - 0.5 * (1.0 - therapy.trough_fraction))


def test_oral_therapy_has_no_interval_modulation(catalog):
    days = np.arange(0.0, 100.0, 7.0)
    assert np.allclose(interval_modulation(days, catalog["enlicitide_oral"]), 1.0)


def test_permanent_editing_does_not_wane(catalog):
    therapy = catalog["verve102_base_editing"]
    days = np.linspace(0.0, 3650.0, 100)
    assert np.allclose(waning_fraction(days, therapy), 1.0)
    assert effect_fraction(days, therapy)[-1] == pytest.approx(1.0)


def test_waning_half_life_halves_the_effect(catalog):
    therapy = replace(catalog["verve102_base_editing"], waning_half_life_days=365.0)
    days = np.array([therapy.plateau_days, therapy.plateau_days + 365.0])
    values = waning_fraction(days, therapy)
    assert values[0] == pytest.approx(1.0)
    assert values[1] == pytest.approx(0.5)


def test_persistence_halves_at_the_discontinuation_half_life():
    adherence = Adherence(discontinuation_half_life_days=100.0)
    values = adherence.persistence(np.array([0.0, 100.0, 200.0]))
    assert np.allclose(values, [1.0, 0.5, 0.25])


def test_perfect_adherence_has_flat_persistence():
    assert np.allclose(Adherence().persistence(np.linspace(0, 1000, 11)), 1.0)


def test_taking_fraction_scales_the_whole_curve(catalog):
    therapy = catalog["statin_high_intensity"]
    days = np.linspace(0.0, 365.0, 50)
    full = effect_fraction(days, therapy)
    half = effect_fraction(days, therapy, Adherence(taking_fraction=0.5))
    assert np.allclose(half, 0.5 * full)


def test_mean_effect_matches_time_weighted_average(catalog):
    therapy = catalog["statin_high_intensity"]
    days = np.linspace(0.0, 365.0, 366)
    value = mean_effect(days, therapy)
    assert 0.9 < value <= 1.0


def test_mean_effect_handles_single_point_grid(catalog):
    therapy = catalog["statin_high_intensity"]
    assert mean_effect(np.array([100.0]), therapy) == pytest.approx(1.0)
    assert mean_effect(np.empty(0), therapy) == 0.0


def test_mean_effect_rejects_zero_span_multi_point_grid(catalog):
    with pytest.raises(ValueError, match="positive duration"):
        mean_effect(np.array([5.0, 5.0]), catalog["statin_high_intensity"])


@pytest.mark.parametrize(
    "days",
    [np.array([[0.0, 1.0]]), np.array([0.0, np.nan]), np.array([2.0, 1.0])],
)
def test_effect_fraction_rejects_bad_grids(catalog, days):
    with pytest.raises(ValueError):
        effect_fraction(days, catalog["statin_high_intensity"])


@pytest.mark.parametrize(
    "kwargs",
    [{"taking_fraction": 1.5}, {"discontinuation_half_life_days": -1.0}, {"start_day": -5.0}],
)
def test_adherence_validation(kwargs):
    with pytest.raises(ValueError):
        Adherence(**kwargs)


def test_every_catalogued_therapy_stays_within_bounds(catalog):
    days = np.linspace(0.0, 1825.0, 200)
    for therapy in catalog.values():
        values = effect_fraction(days, therapy)
        assert np.isfinite(values).all(), therapy.therapy_id
        assert values.min() >= 0.0
        assert values.max() <= 1.0
