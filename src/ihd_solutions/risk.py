"""Event-risk scenarios: only crediting effects the evidence actually supports.

Two crediting pathways exist and they are never double counted:

* **Fixed anchor** - a therapy with its own randomised MACE hazard ratio contributes
  that hazard ratio (which already contains its lipid effect).
* **Per-unit pathway** - therapies without their own outcomes anchor contribute their
  biomarker change to a meta-analytic relationship (CTT per mmol/L LDL-C, BPLTTC per
  5 mmHg systolic). In ``strict`` mode this relationship is only used when a therapy
  in the regimen carries it; in ``surrogate`` mode it is transported to unproven
  LDL-lowering modalities and explicitly flagged.

hsCRP and Lp(a) changes never create event credit in any mode.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from ihd_solutions.catalog import Regimen, Therapy, validate_combination
from ihd_solutions.evidence import (
    CreditDecision,
    CreditMode,
    EffectUnit,
    EvidenceTier,
    TrialAnchor,
    decide_credit,
    is_treatment_evidence,
    may_transport_biomarker,
    select_anchor,
)
from ihd_solutions.lipids import LDL_FLOOR_MG_DL, LipidTrajectory, simulate_biomarkers
from ihd_solutions.pharmacology import PERFECT_ADHERENCE, Adherence, effect_fraction
from ihd_solutions.population import Population
from ihd_solutions.units import (
    DAYS_PER_YEAR,
    MG_DL_PER_MMOL_L,
    require_fraction,
    require_positive,
    trapz_cumulative,
)


@dataclass(frozen=True)
class RiskParameters:
    """Scenario guards and benefit-accrual assumptions."""

    year1_benefit_fraction: float = 0.5
    benefit_ramp_years: float = 1.0
    max_combined_rrr: float = 0.85
    sbp_unit_mmhg: float = 5.0
    #: ``null_effect`` holds trials whose published interval spans 1.0 at HR = 1,
    #: so a neutral trial can never manufacture benefit (or harm) from the noise
    #: in its point estimate. ``point_estimate`` uses the published HR verbatim.
    neutral_ci_policy: str = "null_effect"

    def __post_init__(self) -> None:
        require_fraction("year1_benefit_fraction", self.year1_benefit_fraction, allow_one=True)
        require_positive("benefit_ramp_years", self.benefit_ramp_years)
        require_fraction("max_combined_rrr", self.max_combined_rrr)
        require_positive("sbp_unit_mmhg", self.sbp_unit_mmhg)
        if self.neutral_ci_policy not in NEUTRAL_CI_POLICIES:
            raise ValueError(
                f"neutral_ci_policy must be one of {sorted(NEUTRAL_CI_POLICIES)}, "
                f"got {self.neutral_ci_policy!r}"
            )


NEUTRAL_CI_POLICIES = frozenset({"null_effect", "point_estimate"})
DEFAULT_RISK_PARAMETERS = RiskParameters()


def effective_hazard_ratio(anchor: TrialAnchor, params: RiskParameters) -> float:
    """Hazard ratio actually applied to the event curve for a fixed anchor."""
    if anchor.hr is None:
        raise ValueError(f"{anchor.anchor_id}: no effect estimate")
    if params.neutral_ci_policy == "null_effect" and anchor.neutral_at_report:
        return 1.0
    return float(anchor.hr)


@dataclass(frozen=True)
class EventSimulation:
    """Cumulative-incidence scenario for one regimen in one population."""

    regimen_label: str
    population_id: str
    mode: CreditMode
    days: np.ndarray
    relative_risk: np.ndarray
    cumulative_hazard: np.ndarray
    cumulative_incidence: np.ndarray
    baseline_incidence: np.ndarray
    decisions: tuple[CreditDecision, ...]
    biomarkers: LipidTrajectory
    capped: bool = False
    notes: tuple[str, ...] = field(default=())

    @property
    def years(self) -> np.ndarray:
        return self.days / DAYS_PER_YEAR

    @property
    def horizon_years(self) -> float:
        return float(self.days[-1] / DAYS_PER_YEAR)

    @property
    def event_risk(self) -> float:
        return float(self.cumulative_incidence[-1])

    @property
    def baseline_risk(self) -> float:
        return float(self.baseline_incidence[-1])

    @property
    def absolute_risk_reduction(self) -> float:
        return float(self.baseline_risk - self.event_risk)

    @property
    def relative_risk_reduction(self) -> float:
        if self.baseline_risk <= 0:
            return 0.0
        return float(self.absolute_risk_reduction / self.baseline_risk)

    @property
    def events_avoided_per_1000(self) -> float:
        return float(1000.0 * self.absolute_risk_reduction)

    @property
    def number_needed_to_treat(self) -> float:
        arr = self.absolute_risk_reduction
        return float(1.0 / arr) if arr > 1e-12 else float("inf")

    def summary_row(self) -> dict[str, object]:
        credited = [d.therapy_id for d in self.decisions if d.credited]
        denied = [d.therapy_id for d in self.decisions if not d.credited]
        warnings = sorted({w for d in self.decisions for w in d.warnings})
        return {
            "regimen": self.regimen_label,
            "population_id": self.population_id,
            "mode": self.mode.value,
            "horizon_years": self.horizon_years,
            "baseline_risk_pct": 100.0 * self.baseline_risk,
            "scenario_risk_pct": 100.0 * self.event_risk,
            "absolute_risk_reduction_pct": 100.0 * self.absolute_risk_reduction,
            "relative_risk_reduction_pct": 100.0 * self.relative_risk_reduction,
            "events_avoided_per_1000": self.events_avoided_per_1000,
            "number_needed_to_treat": self.number_needed_to_treat,
            "time_averaged_ldl_mg_dl": self.biomarkers.time_averaged_ldl_mg_dl,
            "ldl_exposure_mmol_years": self.biomarkers.total_ldl_exposure_mmol_years,
            "credited_therapies": ";".join(credited),
            "uncredited_therapies": ";".join(denied),
            "evidence_warnings": ";".join(warnings),
            "combined_rrr_capped": self.capped,
            "scenario_notes": ";".join(self.notes),
        }


def benefit_ramp(years: np.ndarray, params: RiskParameters) -> np.ndarray:
    """CTT benefit is smaller in the first year and accrues thereafter."""
    years = np.asarray(years, dtype=float)
    progress = np.clip(years / params.benefit_ramp_years, 0.0, 1.0)
    return params.year1_benefit_fraction + (1.0 - params.year1_benefit_fraction) * progress


def per_unit_relative_risk(
    magnitude: np.ndarray | float,
    hr_per_unit: float,
    years: np.ndarray | float,
    params: RiskParameters = DEFAULT_RISK_PARAMETERS,
) -> np.ndarray:
    """Apply a per-unit hazard ratio with first-year attenuation."""
    magnitude = np.asarray(magnitude, dtype=float)
    if np.any(magnitude < 0):
        raise ValueError("magnitude must be non-negative")
    require_positive("hr_per_unit", hr_per_unit)
    if hr_per_unit > 1.0:
        raise ValueError("per-unit protective hazard ratio must not exceed 1")
    ramp = benefit_ramp(np.asarray(years, dtype=float), params)
    return np.power(hr_per_unit, magnitude * ramp)


def find_per_unit_anchor(anchors: dict[str, TrialAnchor], unit: EffectUnit) -> TrialAnchor | None:
    matches = [
        a
        for a in anchors.values()
        if a.effect_unit is unit
        and a.has_effect_estimate
        and is_treatment_evidence(a.evidence_tier)
    ]
    return sorted(matches, key=lambda a: a.anchor_id)[0] if matches else None


@dataclass(frozen=True)
class CreditPartition:
    """Which therapies may move the event curve, and through which pathway."""

    fixed: tuple[str, ...]
    per_unit: tuple[str, ...]
    uncredited: tuple[str, ...]
    decisions: tuple[CreditDecision, ...]

    def ldl_participants(self, catalog: dict[str, Therapy], mode: CreditMode) -> tuple[str, ...]:
        """Per-unit therapies, plus unproven LDL agents only in looser modes."""
        participants = list(self.per_unit)
        if mode in {CreditMode.SURROGATE, CreditMode.EXPLORATORY}:
            participants.extend(
                tid
                for tid in self.uncredited
                if catalog[tid].ldl_pct > 0 and may_transport_biomarker(catalog[tid].evidence_tier)
            )
        return tuple(sorted(participants))

    def sbp_participants(self, catalog: dict[str, Therapy], mode: CreditMode) -> tuple[str, ...]:
        participants = list(self.per_unit)
        if mode in {CreditMode.SURROGATE, CreditMode.EXPLORATORY}:
            participants.extend(
                tid
                for tid in self.uncredited
                if catalog[tid].sbp_reduction_mmhg > 0
                and may_transport_biomarker(catalog[tid].evidence_tier)
            )
        return tuple(sorted(participants))


def classify_therapies(
    regimen: Regimen,
    population: Population,
    anchors: dict[str, TrialAnchor],
    mode: CreditMode,
) -> CreditPartition:
    """Split a regimen into fixed-anchor, per-unit and uncredited therapies.

    An uncredited therapy still changes biomarkers, but in ``strict`` mode it is
    excluded from every event pathway. Only ``surrogate``/``exploratory`` modes
    transport a meta-analytic relationship onto it, and that is flagged.
    """
    fixed: list[str] = []
    per_unit: list[str] = []
    uncredited: list[str] = []
    decisions: list[CreditDecision] = []
    for therapy_id in regimen.therapy_ids:
        decision = decide_credit(therapy_id, population.population_id, mode, anchors)
        decisions.append(decision)
        if not decision.credited:
            uncredited.append(therapy_id)
            continue
        anchor = anchors[decision.anchor_id] if decision.anchor_id else None
        if anchor is not None and anchor.effect_unit is not EffectUnit.FIXED:
            per_unit.append(therapy_id)
        else:
            fixed.append(therapy_id)
    return CreditPartition(tuple(fixed), tuple(per_unit), tuple(uncredited), tuple(decisions))


def simulate_events(
    regimen: Regimen,
    catalog: dict[str, Therapy],
    population: Population,
    anchors: dict[str, TrialAnchor],
    days: np.ndarray,
    mode: CreditMode = CreditMode.STRICT,
    params: RiskParameters = DEFAULT_RISK_PARAMETERS,
    adherence: dict[str, Adherence] | None = None,
) -> EventSimulation:
    """Simulate a scenario event curve for ``regimen`` versus no therapy."""
    days = np.asarray(days, dtype=float)
    if days.ndim != 1 or days.size < 2:
        raise ValueError("days must be a one-dimensional grid with at least two points")
    validate_combination(catalog, regimen.therapy_ids)
    adherence = adherence or {}
    biomarkers = simulate_biomarkers(regimen, catalog, population, days, adherence)
    years = days / DAYS_PER_YEAR
    partition = classify_therapies(regimen, population, anchors, mode)
    decisions = partition.decisions
    fixed = partition.fixed
    pathway = partition.per_unit
    decision_by_id = {d.therapy_id: d for d in decisions}
    rr = np.ones_like(days)
    notes: list[str] = []
    held_at_null: list[str] = []
    negative_trials: list[str] = []
    superseded: list[str] = []

    # 1. Fixed randomised anchors: the published hazard ratio, ramped by exposure.
    for therapy_id in fixed:
        decision = decision_by_id[therapy_id]
        anchor = anchors[decision.anchor_id]
        therapy = catalog[therapy_id]
        exposure = effect_fraction(days, therapy, adherence.get(therapy_id, PERFECT_ADHERENCE))
        if anchor.evidence_tier is EvidenceTier.RANDOMIZED_TOPLINE:
            negative_trials.append(therapy_id)
        applied_hr = effective_hazard_ratio(anchor, params)
        if applied_hr == 1.0 and float(anchor.hr) != 1.0:
            held_at_null.append(therapy_id)
        if (
            applied_hr == 1.0
            and mode is not CreditMode.STRICT
            and (therapy.ldl_pct > 0 or therapy.sbp_reduction_mmhg > 0)
        ):
            # Admitting the therapy's own (neutral) randomised evidence displaces
            # biomarker transport, so a looser mode can legitimately score lower.
            superseded.append(therapy_id)
        rr *= np.power(applied_hr, exposure)

    if held_at_null:
        notes.append("neutral_trials_held_at_null:" + ";".join(sorted(held_at_null)))
    if negative_trials:
        notes.append("negative_outcomes_trial_no_benefit:" + ";".join(sorted(negative_trials)))
    if superseded:
        notes.append(
            "neutral_anchor_supersedes_biomarker_transport:" + ";".join(sorted(superseded))
        )

    # 2. LDL pathway: biomarker change from per-unit therapies, plus unproven
    #    LDL agents only when the mode explicitly permits surrogate transport.
    ldl_participants = partition.ldl_participants(catalog, mode)
    ldl_remaining = np.ones_like(days)
    for therapy_id in ldl_participants:
        contribution = biomarkers.components.get(therapy_id)
        if contribution is not None:
            ldl_remaining *= 1.0 - contribution
    pathway_ldl = np.maximum(population.ldl_c_mg_dl * ldl_remaining, LDL_FLOOR_MG_DL)
    delta_mmol = (population.ldl_c_mg_dl - pathway_ldl) / MG_DL_PER_MMOL_L
    ldl_anchor = _pathway_anchor(pathway, anchors, population, mode, EffectUnit.PER_MMOL_L_LDL)
    transported = [tid for tid in ldl_participants if tid in partition.uncredited]
    if ldl_anchor is not None and np.any(delta_mmol > 0):
        rr *= per_unit_relative_risk(delta_mmol, float(ldl_anchor.hr), years, params)
        if transported:
            notes.append("ldl_relationship_transported_to:" + ";".join(transported))
    elif np.any(delta_mmol > 0):
        notes.append("ldl_lowering_not_credited_no_relationship_available")
    uncredited_ldl = [tid for tid in partition.uncredited if catalog[tid].ldl_pct > 0]
    if uncredited_ldl and not transported:
        notes.append("ldl_lowering_without_outcomes_evidence:" + ";".join(uncredited_ldl))

    # 3. Blood-pressure pathway, again excluding therapies with their own anchor.
    sbp_drop = np.zeros_like(days)
    sbp_participants = partition.sbp_participants(catalog, mode)
    for therapy_id in sbp_participants:
        therapy = catalog[therapy_id]
        if therapy.sbp_reduction_mmhg <= 0:
            continue
        exposure = effect_fraction(days, therapy, adherence.get(therapy_id, PERFECT_ADHERENCE))
        sbp_drop += exposure * therapy.sbp_reduction_mmhg
    sbp_anchor = _pathway_anchor(pathway, anchors, population, mode, EffectUnit.PER_5_MMHG_SBP)
    if sbp_anchor is not None and np.any(sbp_drop > 0):
        rr *= per_unit_relative_risk(
            sbp_drop / params.sbp_unit_mmhg, float(sbp_anchor.hr), years, params
        )
    if any(catalog[tid].hscrp_pct > 0 for tid in partition.uncredited):
        notes.append("biomarker_only_inflammation_not_credited")
    if any(catalog[tid].lpa_pct > 0 for tid in partition.uncredited):
        notes.append("lpa_lowering_not_outcomes_validated")

    if np.any(~np.isfinite(rr)) or np.any(rr <= 0):
        raise ValueError("non-finite or non-positive relative risk")
    floor_rr = 1.0 - params.max_combined_rrr
    capped = bool(np.any(rr < floor_rr))
    if capped:
        rr = np.maximum(rr, floor_rr)
        notes.append("combined_relative_risk_reduction_capped")

    hazard = population.baseline_annual_mace_rate * rr
    cumulative_hazard = trapz_cumulative(hazard, years)
    incidence = -np.expm1(-cumulative_hazard)
    baseline_hazard = np.full_like(days, population.baseline_annual_mace_rate)
    baseline_incidence = -np.expm1(-trapz_cumulative(baseline_hazard, years))
    return EventSimulation(
        regimen_label=regimen.label,
        population_id=population.population_id,
        mode=mode,
        days=days,
        relative_risk=rr,
        cumulative_hazard=cumulative_hazard,
        cumulative_incidence=incidence,
        baseline_incidence=baseline_incidence,
        decisions=decisions,
        biomarkers=biomarkers,
        capped=capped,
        notes=tuple(notes),
    )


def _pathway_anchor(
    pathway: list[str],
    anchors: dict[str, TrialAnchor],
    population: Population,
    mode: CreditMode,
    unit: EffectUnit,
) -> TrialAnchor | None:
    """Use a carried per-unit relationship, else transport it in looser modes."""
    for therapy_id in pathway:
        anchor = select_anchor(anchors, therapy_id, population.population_id, mode)
        if anchor is None or anchor.effect_unit is not unit or not anchor.has_effect_estimate:
            continue
        decision = decide_credit(therapy_id, population.population_id, mode, anchors)
        if decision.credited:
            return anchor
    if mode in {CreditMode.SURROGATE, CreditMode.EXPLORATORY}:
        candidate = find_per_unit_anchor(anchors, unit)
        if candidate is not None and candidate.applies_to(population.population_id, mode):
            return candidate
    return None


def decisions_frame(simulations: list[EventSimulation], unique: bool = True) -> pd.DataFrame:
    """Tabulate credit decisions.

    A decision depends only on the therapy, population and mode - never on the
    rest of the regimen - so the default output is deduplicated.
    """
    rows = []
    for simulation in simulations:
        for decision in simulation.decisions:
            rows.append(
                {
                    "regimen": simulation.regimen_label,
                    "population_id": simulation.population_id,
                    "mode": simulation.mode.value,
                    "therapy_id": decision.therapy_id,
                    "credited": decision.credited,
                    "anchor_id": decision.anchor_id or "",
                    "hazard_ratio": decision.hr if decision.hr is not None else np.nan,
                    "reason": decision.reason,
                    "warnings": ";".join(decision.warnings),
                }
            )
    frame = pd.DataFrame(rows)
    if not unique or frame.empty:
        return frame
    keys = ["therapy_id", "population_id", "mode"]
    counts = frame.groupby(keys, sort=True).size().rename("n_regimens_observed").reset_index()
    deduped = frame.drop_duplicates(subset=keys, keep="first").drop(columns=["regimen"])
    return (
        deduped.merge(counts, on=keys, how="left")
        .sort_values(keys, kind="stable")
        .reset_index(drop=True)
    )
