"""Inflammatory cascade scenarios (biomarkers only - never an event model)."""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
import pandas as pd

from ihd_solutions.catalog import Therapy
from ihd_solutions.units import check_all_finite, require_fraction, require_positive

STATE_NAMES = ("oxLDL", "NLRP3", "IL1B", "IL6", "hsCRP", "MMP")


@dataclass(frozen=True)
class InflammationParameters:
    """Dimensionless cascade rates; illustrative, not fitted pharmacology."""

    oxldl_influx: float = 1.0
    oxldl_clearance: float = 0.5
    nlrp3_activation: float = 0.8
    nlrp3_decay: float = 0.6
    il1_production: float = 0.7
    il1_decay: float = 0.5
    il6_production: float = 0.6
    il6_decay: float = 0.5
    crp_production: float = 0.9
    crp_decay: float = 0.4
    mmp_production: float = 0.5
    mmp_decay: float = 0.3

    def __post_init__(self) -> None:
        for name in self.__dataclass_fields__:
            require_positive(name, getattr(self, name))


@dataclass(frozen=True)
class InflammatoryIntervention:
    """Target-specific inhibition fractions in [0, 1)."""

    therapy_id: str = "untreated"
    nlrp3_inhibition: float = 0.0
    il1_inhibition: float = 0.0
    il6_inhibition: float = 0.0
    oxldl_reduction: float = 0.0
    mmp_suppression: float = 0.0

    def __post_init__(self) -> None:
        for name in (
            "nlrp3_inhibition",
            "il1_inhibition",
            "il6_inhibition",
            "oxldl_reduction",
            "mmp_suppression",
        ):
            require_fraction(name, getattr(self, name))


UNTREATED = InflammatoryIntervention()

#: Explicit mechanistic mapping. Methotrexate deliberately engages nothing,
#: reproducing its role as the CIRT negative control.
_TARGET_MAP: dict[str, dict[str, float]] = {
    "colchicine": {"nlrp3_inhibition": 0.30, "mmp_suppression": 0.20},
    "il1_axis": {"il1_inhibition": 0.70, "mmp_suppression": 0.25},
    "il6_axis": {"il6_inhibition": 0.80},
    "nlrp3": {"nlrp3_inhibition": 0.60, "mmp_suppression": 0.25},
    "mtx": {},
}


def intervention_for(therapy: Therapy) -> InflammatoryIntervention:
    """Map a catalogue therapy onto cascade targets using its mechanism group."""
    kwargs = dict(_TARGET_MAP.get(therapy.exclusive_group, {}))
    if therapy.ldl_pct > 0:
        kwargs["oxldl_reduction"] = therapy.ldl_pct
    return InflammatoryIntervention(therapy_id=therapy.therapy_id, **kwargs)


def initial_state() -> np.ndarray:
    return np.array([1.0, 0.20, 0.15, 0.10, 2.0, 1.0], dtype=float)


def _derivatives(
    state: np.ndarray,
    params: InflammationParameters,
    intervention: InflammatoryIntervention,
) -> np.ndarray:
    oxldl, nlrp3, il1, il6, crp, mmp = state
    influx = params.oxldl_influx * (1.0 - intervention.oxldl_reduction)
    d_oxldl = influx - params.oxldl_clearance * oxldl
    d_nlrp3 = (
        params.nlrp3_activation * oxldl * (1.0 - intervention.nlrp3_inhibition)
        - params.nlrp3_decay * nlrp3
    )
    d_il1 = (
        params.il1_production * nlrp3 * (1.0 - intervention.il1_inhibition) - params.il1_decay * il1
    )
    d_il6 = (
        params.il6_production * il1 * (1.0 - intervention.il6_inhibition) - params.il6_decay * il6
    )
    d_crp = params.crp_production * il6 - params.crp_decay * crp
    d_mmp = (
        params.mmp_production * (il1 + oxldl) * (1.0 - intervention.mmp_suppression)
        - params.mmp_decay * mmp
    )
    return np.array([d_oxldl, d_nlrp3, d_il1, d_il6, d_crp, d_mmp], dtype=float)


def simulate_inflammation(
    intervention: InflammatoryIntervention = UNTREATED,
    params: InflammationParameters = InflammationParameters(),
    days: float = 180.0,
    step_days: float = 0.25,
    state0: np.ndarray | None = None,
) -> pd.DataFrame:
    """Integrate the cascade with simultaneous RK4 steps and strict finite checks."""
    days = require_positive("days", days)
    step_days = require_positive("step_days", step_days)
    if step_days > days:
        raise ValueError("step_days must not exceed the simulated horizon")
    n_steps = max(1, int(round(days / step_days)))
    state = np.array(initial_state() if state0 is None else state0, dtype=float)
    if state.shape != (len(STATE_NAMES),):
        raise ValueError(f"state0 must have shape ({len(STATE_NAMES)},)")
    check_all_finite("state0", state)
    if np.any(state < 0):
        raise ValueError("state0 must be non-negative")
    history = np.empty((n_steps + 1, len(STATE_NAMES)), dtype=float)
    history[0] = state
    h = days / n_steps
    for index in range(n_steps):
        with np.errstate(over="ignore", invalid="ignore"):
            k1 = _derivatives(state, params, intervention)
            k2 = _derivatives(state + 0.5 * h * k1, params, intervention)
            k3 = _derivatives(state + 0.5 * h * k2, params, intervention)
            k4 = _derivatives(state + h * k3, params, intervention)
            candidate = state + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        if not np.isfinite(candidate).all():
            raise ValueError(f"non-finite cascade state at step {index}")
        state = np.maximum(candidate, 0.0)
        history[index + 1] = state
    frame = pd.DataFrame(history, columns=STATE_NAMES)
    frame.insert(0, "day", np.linspace(0.0, days, n_steps + 1))
    frame.insert(0, "therapy_id", intervention.therapy_id)
    return frame


def steady_state(
    params: InflammationParameters = InflammationParameters(),
    intervention: InflammatoryIntervention = UNTREATED,
) -> np.ndarray:
    """Analytic equilibrium used to verify the integrator."""
    oxldl = params.oxldl_influx * (1.0 - intervention.oxldl_reduction) / params.oxldl_clearance
    nlrp3 = (
        params.nlrp3_activation * oxldl * (1.0 - intervention.nlrp3_inhibition) / params.nlrp3_decay
    )
    il1 = params.il1_production * nlrp3 * (1.0 - intervention.il1_inhibition) / params.il1_decay
    il6 = params.il6_production * il1 * (1.0 - intervention.il6_inhibition) / params.il6_decay
    crp = params.crp_production * il6 / params.crp_decay
    mmp = (
        params.mmp_production
        * (il1 + oxldl)
        * (1.0 - intervention.mmp_suppression)
        / params.mmp_decay
    )
    return np.array([oxldl, nlrp3, il1, il6, crp, mmp], dtype=float)


def with_parameters(**changes: float) -> InflammationParameters:
    return replace(InflammationParameters(), **changes)
