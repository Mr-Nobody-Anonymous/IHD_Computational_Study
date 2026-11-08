"""Therapy-driven plaque remodelling scenarios (thin-wall proxy, not FEA)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from ihd_solutions.catalog import Regimen, Therapy, validate_combination
from ihd_solutions.lipids import simulate_biomarkers
from ihd_solutions.pharmacology import PERFECT_ADHERENCE, Adherence, effect_fraction
from ihd_solutions.population import Population
from ihd_solutions.units import require_nonnegative, require_positive

#: Literature heuristics with documented exceptions - not universal rupture criteria.
CAP_THRESHOLD_UM = 65.0
STRESS_THRESHOLD_KPA = 300.0


@dataclass(frozen=True)
class PlaqueParameters:
    reference_radius_mm: float = 1.5
    initial_core_area_mm2: float = 0.40
    core_to_plaque_area: float = 2.2
    initial_cap_um: float = 120.0
    wall_thickness_mm: float = 0.30
    mean_pressure_kpa: float = 13.3
    lipid_influx_coeff: float = 0.010
    cap_loss_coeff: float = 0.90
    cap_repair_coeff: float = 0.30
    mmp_baseline: float = 1.0
    core_floor_mm2: float = 0.26
    cap_ceiling_um: float = 180.0

    def __post_init__(self) -> None:
        for name in (
            "reference_radius_mm",
            "initial_core_area_mm2",
            "core_to_plaque_area",
            "initial_cap_um",
            "wall_thickness_mm",
            "mean_pressure_kpa",
            "cap_loss_coeff",
            "mmp_baseline",
            "core_floor_mm2",
            "cap_ceiling_um",
        ):
            require_positive(name, getattr(self, name))
        require_nonnegative("lipid_influx_coeff", self.lipid_influx_coeff)
        require_nonnegative("cap_repair_coeff", self.cap_repair_coeff)
        if self.core_floor_mm2 > self.initial_core_area_mm2:
            raise ValueError("core floor cannot exceed the initial core area")
        if self.initial_core_area_mm2 >= np.pi * self.reference_radius_mm**2:
            raise ValueError("initial core area must fit inside the reference lumen")

    @property
    def initial_stenosis_pct(self) -> float:
        plaque_area = self.core_to_plaque_area * self.initial_core_area_mm2
        return float(100.0 * plaque_area / (np.pi * self.reference_radius_mm**2))

    def lumen_radius_mm(self, core_area_mm2: float) -> float:
        plaque_area = self.core_to_plaque_area * core_area_mm2
        remaining = np.pi * self.reference_radius_mm**2 - plaque_area
        if remaining <= 0:
            raise ValueError("plaque area occludes the reference lumen")
        return float(np.sqrt(remaining / np.pi))

    def hoop_stress_kpa(self, core_area_mm2: float, cap_um: float) -> float:
        radius = self.lumen_radius_mm(core_area_mm2)
        thickness_mm = max(cap_um, 1.0) / 1000.0
        return float(self.mean_pressure_kpa * radius / thickness_mm)


DEFAULT_PLAQUE = PlaqueParameters()


@dataclass(frozen=True)
class PlaqueOutcome:
    regimen_label: str
    population_id: str
    first_crossing_month: int | None
    min_cap_um: float
    max_stress_kpa: float
    final_core_area_mm2: float
    months_simulated: int
    intact_months: int
    timeline: pd.DataFrame

    @property
    def crossed(self) -> bool:
        return self.first_crossing_month is not None


def simulate_plaque(
    regimen: Regimen,
    catalog: dict[str, Therapy],
    population: Population,
    months: int = 60,
    params: PlaqueParameters = DEFAULT_PLAQUE,
    adherence: dict[str, Adherence] | None = None,
) -> PlaqueOutcome:
    """Monthly remodelling with an explicit intact horizon and counterfactual flag."""
    if months <= 0:
        raise ValueError("months must be positive")
    validate_combination(catalog, regimen.therapy_ids)
    adherence = adherence or {}
    days = np.arange(months + 1, dtype=float) * 30.4375
    biomarkers = simulate_biomarkers(regimen, catalog, population, days, adherence)
    therapies = [catalog[therapy_id] for therapy_id in regimen.therapy_ids]

    mmp = np.full(months + 1, params.mmp_baseline)
    for therapy in therapies:
        if therapy.mmp_factor >= 1.0:
            continue
        exposure = effect_fraction(
            days, therapy, adherence.get(therapy.therapy_id, PERFECT_ADHERENCE)
        )
        mmp *= 1.0 - exposure * (1.0 - therapy.mmp_factor)

    core = np.empty(months + 1)
    cap = np.empty(months + 1)
    stress = np.empty(months + 1)
    stenosis = np.empty(months + 1)
    core[0] = params.initial_core_area_mm2
    cap[0] = params.initial_cap_um
    stress[0] = params.hoop_stress_kpa(core[0], cap[0])
    stenosis[0] = params.initial_stenosis_pct

    crossing: int | None = None
    for month in range(months):
        ldl_drive = biomarkers.ldl_mg_dl[month] / population.ldl_c_mg_dl
        influx = params.lipid_influx_coeff * ldl_drive
        with np.errstate(over="ignore", invalid="ignore"):
            next_core = max(core[month] + influx, params.core_floor_mm2)
            loss = params.cap_loss_coeff * mmp[month] * (next_core / params.initial_core_area_mm2)
            repair = params.cap_repair_coeff * (1.0 - ldl_drive)
            next_cap = min(cap[month] - loss + repair, params.cap_ceiling_um)
            next_cap = max(next_cap, 5.0)
        if not np.isfinite([next_core, next_cap]).all():
            raise ValueError(f"non-finite plaque state at month {month + 1}")
        core[month + 1] = next_core
        cap[month + 1] = next_cap
        stress[month + 1] = params.hoop_stress_kpa(next_core, next_cap)
        plaque_area = params.core_to_plaque_area * next_core
        stenosis[month + 1] = 100.0 * plaque_area / (np.pi * params.reference_radius_mm**2)
        if crossing is None and (
            next_cap <= CAP_THRESHOLD_UM and stress[month + 1] >= STRESS_THRESHOLD_KPA
        ):
            crossing = month + 1

    intact = months if crossing is None else crossing
    timeline = pd.DataFrame(
        {
            "regimen": regimen.label,
            "population_id": population.population_id,
            "month": np.arange(months + 1),
            "ldl_c_mg_dl": biomarkers.ldl_mg_dl,
            "mmp_index": mmp,
            "core_area_mm2": core,
            "cap_thickness_um": cap,
            "hoop_stress_kpa": stress,
            "stenosis_pct": stenosis,
            "post_threshold_counterfactual": np.arange(months + 1) > intact,
        }
    )
    return PlaqueOutcome(
        regimen_label=regimen.label,
        population_id=population.population_id,
        first_crossing_month=crossing,
        min_cap_um=float(cap[: intact + 1].min()),
        max_stress_kpa=float(stress[: intact + 1].max()),
        final_core_area_mm2=float(core[intact]),
        months_simulated=months,
        intact_months=intact,
        timeline=timeline,
    )
