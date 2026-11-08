"""Enumerate and compare every admissible combination of candidate solutions."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from itertools import combinations

import numpy as np
import pandas as pd

from ihd_solutions.catalog import Regimen, Therapy, TherapyClass, validate_combination
from ihd_solutions.evidence import CreditMode, TrialAnchor
from ihd_solutions.pharmacology import Adherence
from ihd_solutions.population import Population, background_regimen_ids
from ihd_solutions.risk import (
    DEFAULT_RISK_PARAMETERS,
    EventSimulation,
    RiskParameters,
    simulate_events,
)


@dataclass(frozen=True)
class EnumerationSpec:
    """Which combinations to build on top of mandatory background therapy."""

    max_optional: int = 2
    background: tuple[str, ...] = background_regimen_ids()
    exclude_classes: tuple[TherapyClass, ...] = (TherapyClass.REGENERATIVE,)
    candidate_ids: tuple[str, ...] | None = None

    def __post_init__(self) -> None:
        if self.max_optional < 0:
            raise ValueError("max_optional must be non-negative")


def enumerate_regimens(
    catalog: dict[str, Therapy], spec: EnumerationSpec = EnumerationSpec()
) -> list[Regimen]:
    """Background alone plus every admissible combination up to ``max_optional``."""
    for therapy_id in spec.background:
        if therapy_id not in catalog:
            raise ValueError(f"background therapy not in catalog: {therapy_id}")
    validate_combination(catalog, tuple(spec.background))
    background_groups = {catalog[tid].exclusive_group for tid in spec.background}
    if spec.candidate_ids is None:
        candidates = sorted(
            therapy.therapy_id
            for therapy in catalog.values()
            if therapy.therapy_id not in spec.background
            and therapy.therapy_class not in spec.exclude_classes
            and therapy.exclusive_group not in background_groups
        )
    else:
        unknown = [tid for tid in spec.candidate_ids if tid not in catalog]
        if unknown:
            raise ValueError(f"unknown candidate ids: {unknown}")
        candidates = sorted(
            tid
            for tid in spec.candidate_ids
            if tid not in spec.background and catalog[tid].exclusive_group not in background_groups
        )
    regimens: list[Regimen] = [validate_combination(catalog, tuple(spec.background))]
    seen: set[tuple[str, ...]] = {regimens[0].therapy_ids}
    for size in range(1, spec.max_optional + 1):
        for combo in combinations(candidates, size):
            groups = [catalog[tid].exclusive_group for tid in combo]
            if len(set(groups)) != len(groups):
                continue
            regimen = validate_combination(catalog, tuple(spec.background) + combo)
            if regimen.therapy_ids in seen:
                continue
            seen.add(regimen.therapy_ids)
            regimens.append(regimen)
    return regimens


def regimen_burden(regimen: Regimen, catalog: dict[str, Therapy]) -> float:
    """Additive burden proxy: dosing frequency, injections and procedures."""
    return float(sum(catalog[therapy_id].burden_score for therapy_id in regimen.therapy_ids))


def evaluate_regimens(
    regimens: Sequence[Regimen],
    catalog: dict[str, Therapy],
    population: Population,
    anchors: dict[str, TrialAnchor],
    days: np.ndarray,
    modes: Iterable[CreditMode] = (CreditMode.STRICT, CreditMode.SURROGATE, CreditMode.EXPLORATORY),
    params: RiskParameters = DEFAULT_RISK_PARAMETERS,
    adherence: dict[str, Adherence] | None = None,
) -> tuple[pd.DataFrame, list[EventSimulation]]:
    """Score every regimen in every requested evidence mode."""
    modes = tuple(modes)
    if not modes:
        raise ValueError("at least one credit mode is required")
    if not regimens:
        raise ValueError("no regimens to evaluate")
    rows: dict[tuple[str, ...], dict[str, object]] = {}
    simulations: list[EventSimulation] = []
    for regimen in regimens:
        key = regimen.therapy_ids
        rows.setdefault(
            key,
            {
                "regimen": regimen.label,
                "n_therapies": len(regimen),
                "therapy_ids": ";".join(regimen.therapy_ids),
                "population_id": population.population_id,
                "burden_score": regimen_burden(regimen, catalog),
            },
        )
        for mode in modes:
            simulation = simulate_events(
                regimen, catalog, population, anchors, days, mode, params, adherence
            )
            simulations.append(simulation)
            prefix = mode.value
            rows[key][f"{prefix}_events_avoided_per_1000"] = simulation.events_avoided_per_1000
            rows[key][f"{prefix}_rrr_pct"] = 100.0 * simulation.relative_risk_reduction
            rows[key][f"{prefix}_nnt"] = simulation.number_needed_to_treat
            rows[key][f"{prefix}_n_credited"] = sum(1 for d in simulation.decisions if d.credited)
            rows[key][f"{prefix}_warnings"] = ";".join(
                sorted({w for d in simulation.decisions for w in d.warnings})
            )
            rows[key]["time_averaged_ldl_mg_dl"] = simulation.biomarkers.time_averaged_ldl_mg_dl
            rows[key]["ldl_exposure_mmol_years"] = (
                simulation.biomarkers.total_ldl_exposure_mmol_years
            )
    frame = pd.DataFrame(rows.values())
    sort_column = f"{modes[0].value}_events_avoided_per_1000"
    return (
        frame.sort_values(
            [sort_column, "burden_score", "regimen"],
            ascending=[False, True, True],
            kind="stable",
        ).reset_index(drop=True),
        simulations,
    )


def pareto_front(
    frame: pd.DataFrame, benefit_column: str, cost_column: str = "burden_score"
) -> pd.DataFrame:
    """Non-dominated regimens: more benefit at no greater burden."""
    if benefit_column not in frame or cost_column not in frame:
        raise ValueError("frame is missing the benefit or cost column")
    benefit = frame[benefit_column].to_numpy(dtype=float)
    cost = frame[cost_column].to_numpy(dtype=float)
    keep = np.ones(len(frame), dtype=bool)
    for index in range(len(frame)):
        dominates = (benefit >= benefit[index]) & (cost <= cost[index])
        strictly_better = (benefit > benefit[index]) | (cost < cost[index])
        if np.any(dominates & strictly_better):
            keep[index] = False
    front = frame.loc[keep].copy()
    return front.sort_values(cost_column, kind="stable").reset_index(drop=True)


def compare_populations(
    regimens: Sequence[Regimen],
    catalog: dict[str, Therapy],
    populations: dict[str, Population],
    anchors: dict[str, TrialAnchor],
    days: np.ndarray,
    mode: CreditMode = CreditMode.STRICT,
    population_ids: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Same regimens across populations, exposing population-specific evidence."""
    ids = list(population_ids) if population_ids is not None else sorted(populations)
    rows = []
    for population_id in ids:
        if population_id not in populations:
            raise ValueError(f"unknown population: {population_id}")
        population = populations[population_id]
        for regimen in regimens:
            simulation = simulate_events(regimen, catalog, population, anchors, days, mode=mode)
            rows.append(simulation.summary_row())
    return pd.DataFrame(rows)


def expected_regimen_count(catalog: dict[str, Therapy], spec: EnumerationSpec) -> int:
    """Independent combinatorial count used to verify the enumerator in tests."""
    background_groups = {catalog[tid].exclusive_group for tid in spec.background}
    candidates = [
        therapy
        for therapy in catalog.values()
        if therapy.therapy_id not in spec.background
        and therapy.therapy_class not in spec.exclude_classes
        and therapy.exclusive_group not in background_groups
    ]
    if spec.candidate_ids is not None:
        allowed = set(spec.candidate_ids)
        candidates = [therapy for therapy in candidates if therapy.therapy_id in allowed]
    groups: dict[str, int] = {}
    for therapy in candidates:
        groups[therapy.exclusive_group] = groups.get(therapy.exclusive_group, 0) + 1
    sizes = sorted(groups.values())
    total = 1
    for size in range(1, spec.max_optional + 1):
        for combo in combinations(sizes, size):
            product = 1
            for value in combo:
                product *= value
            total += product
    return total
