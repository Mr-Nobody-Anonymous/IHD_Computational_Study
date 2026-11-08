"""Biomarker trajectories for a regimen: LDL-C, ApoB, Lp(a), hsCRP and pressure."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from ihd_solutions.catalog import Regimen, Therapy, validate_combination
from ihd_solutions.pharmacology import PERFECT_ADHERENCE, Adherence, effect_fraction
from ihd_solutions.population import Population
from ihd_solutions.units import (
    DAYS_PER_YEAR,
    MG_DL_PER_MMOL_L,
    mgdl_to_mmoll,
    require_positive,
    trapz_cumulative,
)

#: Physiological floor so multiplicative stacking cannot produce absurd LDL values.
LDL_FLOOR_MG_DL = 10.0
APOB_FLOOR_MG_DL = 15.0
LPA_FLOOR_NMOL_L = 2.0
HSCRP_FLOOR_MG_L = 0.2


@dataclass(frozen=True)
class LipidTrajectory:
    """Time-resolved biomarker state plus cumulative exposure summaries."""

    days: np.ndarray
    ldl_mg_dl: np.ndarray
    apob_mg_dl: np.ndarray
    lpa_nmol_l: np.ndarray
    hscrp_mg_l: np.ndarray
    systolic_bp_mmhg: np.ndarray
    ldl_reduction_fraction: np.ndarray
    cumulative_ldl_mmol_years: np.ndarray
    components: dict[str, np.ndarray] = field(default_factory=dict)

    @property
    def years(self) -> np.ndarray:
        return self.days / DAYS_PER_YEAR

    @property
    def time_averaged_ldl_mg_dl(self) -> float:
        span = self.days[-1] - self.days[0]
        if span <= 0:
            return float(self.ldl_mg_dl[0])
        return float(self.cumulative_ldl_mmol_years[-1] * DAYS_PER_YEAR / span * MG_DL_PER_MMOL_L)

    @property
    def total_ldl_exposure_mmol_years(self) -> float:
        return float(self.cumulative_ldl_mmol_years[-1])

    def fraction_time_below(self, threshold_mg_dl: float) -> float:
        require_positive("threshold_mg_dl", threshold_mg_dl)
        below = self.ldl_mg_dl < threshold_mg_dl
        if self.days.size < 2:
            return float(below[0])
        weights = np.diff(self.days)
        mid = below[1:] & below[:-1]
        return float(np.dot(mid.astype(float), weights) / weights.sum())


def _stack_remaining(
    days: np.ndarray,
    therapies: list[Therapy],
    attribute: str,
    adherence: dict[str, Adherence],
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Multiplicative combination: each agent acts on the remaining fraction."""
    remaining = np.ones_like(days, dtype=float)
    components: dict[str, np.ndarray] = {}
    for therapy in therapies:
        potency = float(getattr(therapy, attribute))
        if potency <= 0:
            continue
        fraction = effect_fraction(
            days, therapy, adherence.get(therapy.therapy_id, PERFECT_ADHERENCE)
        )
        components[therapy.therapy_id] = fraction * potency
        remaining *= 1.0 - fraction * potency
    return remaining, components


def simulate_biomarkers(
    regimen: Regimen,
    catalog: dict[str, Therapy],
    population: Population,
    days: np.ndarray,
    adherence: dict[str, Adherence] | None = None,
) -> LipidTrajectory:
    """Simulate biomarkers for a validated regimen on a monotone day grid."""
    days = np.asarray(days, dtype=float)
    if days.ndim != 1 or days.size < 2:
        raise ValueError("days must be a one-dimensional grid with at least two points")
    if not np.isfinite(days).all() or np.any(np.diff(days) <= 0):
        raise ValueError("days must be finite and strictly increasing")
    validate_combination(catalog, regimen.therapy_ids)
    adherence = adherence or {}
    unknown = [key for key in adherence if key not in regimen.therapy_ids]
    if unknown:
        raise ValueError(f"adherence supplied for therapies outside the regimen: {unknown}")
    therapies = [catalog[therapy_id] for therapy_id in regimen.therapy_ids]

    ldl_remaining, ldl_components = _stack_remaining(days, therapies, "ldl_pct", adherence)
    apob_remaining, _ = _stack_remaining(days, therapies, "apob_pct", adherence)
    lpa_remaining, _ = _stack_remaining(days, therapies, "lpa_pct", adherence)
    hscrp_remaining, _ = _stack_remaining(days, therapies, "hscrp_pct", adherence)

    sbp_drop = np.zeros_like(days)
    for therapy in therapies:
        if therapy.sbp_reduction_mmhg <= 0:
            continue
        fraction = effect_fraction(
            days, therapy, adherence.get(therapy.therapy_id, PERFECT_ADHERENCE)
        )
        sbp_drop += fraction * therapy.sbp_reduction_mmhg

    ldl = np.maximum(population.ldl_c_mg_dl * ldl_remaining, LDL_FLOOR_MG_DL)
    apob = np.maximum(population.apob_mg_dl * apob_remaining, APOB_FLOOR_MG_DL)
    lpa = np.maximum(population.lp_a_nmol_l * lpa_remaining, LPA_FLOOR_NMOL_L)
    hscrp = np.maximum(population.hscrp_mg_l * hscrp_remaining, HSCRP_FLOOR_MG_L)
    sbp = np.maximum(population.systolic_bp_mmhg - sbp_drop, 80.0)

    reduction = 1.0 - ldl / population.ldl_c_mg_dl
    exposure = trapz_cumulative(mgdl_to_mmoll(ldl), days / DAYS_PER_YEAR)
    for name, array in (
        ("LDL-C", ldl),
        ("ApoB", apob),
        ("Lp(a)", lpa),
        ("hsCRP", hscrp),
        ("systolic BP", sbp),
    ):
        if not np.isfinite(array).all():
            raise ValueError(f"non-finite {name} trajectory")

    return LipidTrajectory(
        days=days,
        ldl_mg_dl=ldl,
        apob_mg_dl=apob,
        lpa_nmol_l=lpa,
        hscrp_mg_l=hscrp,
        systolic_bp_mmhg=sbp,
        ldl_reduction_fraction=reduction,
        cumulative_ldl_mmol_years=exposure,
        components=ldl_components,
    )


def trajectory_frame(
    trajectory: LipidTrajectory, regimen_label: str, population_id: str
) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "regimen": regimen_label,
            "population_id": population_id,
            "day": trajectory.days,
            "year": trajectory.years,
            "ldl_c_mg_dl": trajectory.ldl_mg_dl,
            "ldl_c_mmol_l": mgdl_to_mmoll(trajectory.ldl_mg_dl),
            "ldl_reduction_pct": 100.0 * trajectory.ldl_reduction_fraction,
            "apob_mg_dl": trajectory.apob_mg_dl,
            "lpa_nmol_l": trajectory.lpa_nmol_l,
            "hscrp_mg_l": trajectory.hscrp_mg_l,
            "systolic_bp_mmhg": trajectory.systolic_bp_mmhg,
            "cumulative_ldl_mmol_years": trajectory.cumulative_ldl_mmol_years,
        }
    )
