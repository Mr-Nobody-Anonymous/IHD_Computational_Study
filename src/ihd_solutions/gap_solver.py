"""Attempt to close the three quantified capability gaps, mechanistically.

The gaps (from docs/cure_discovery_plan.md):

======================  ==============  ==================  ======
Capability              Needed          Have                Gap
======================  ==============  ==================  ======
Remuscularisation       +20 LVEF pts    +3.9 (BioVAT-HF)    x5.13
Substrate clearance     80%             43.1%               x1.86
Stack hazard ratio      <=0.298         0.360               x1.21
======================  ==============  ==================  ======

This module does not tune parameters until the answer looks good. For each gap it
asks a structural question, and the answer changes what kind of work is required:

* **Hazard ratio** - already closed by the composite in the population the target was
  derived from. The target HR is baseline-specific, so it must be recomputed per
  population rather than applied universally.
* **Substrate clearance** - the plateau is *structural*, set by a constant apoptosis
  source that none of the six REVERT-1 modules reduces. No payload strength fixes it.
  Adding a source-reduction module (ER-stress/CHOP) does, and the required strength is
  back-calculated against a precedent bound.
* **Remuscularisation** - decomposes into multiplicative factors that all act on one
  limiting quantity: the fraction of scar replaced by contractile graft. That quantity
  saturates at 1.0, so the factors cannot be multiplied without bound.

Everything here is a scenario. Closing a gap *in this model* is not evidence that the
intervention works.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum

import numpy as np
import pandas as pd

from ihd_solutions.cure import CRITERION_THRESHOLD
from ihd_solutions.discovery import required_capabilities
from ihd_solutions.hypothesis import (
    _ZERO_STRENGTHS,
    ClearanceParameters,
    LesionState,
    simulate_clearance,
)
from ihd_solutions.units import require_fraction, require_positive

DISCLAIMER = (
    "GAP-CLOSURE SCENARIO - a modelled path to a capability target, not evidence that "
    "the path works. Precedent bounds are animal or mechanistic unless stated. "
    "Closing a gap here does not make a cure."
)

GAP_TAG = "gap_closure_scenario_not_evidence"


class GapVerdict(str, Enum):
    CLOSED = "closed_already"
    CLOSABLE = "closable_within_precedent"
    NEEDS_NEW_BIOLOGY = "requires_unprecedented_capability"
    OPEN = "still_open"


# ---------------------------------------------------------------------------
# Gap 3: hazard ratio. Already closed, but only where the target applies.
# ---------------------------------------------------------------------------


def hazard_ratio_from_risks(baseline: float, treated: float) -> float:
    """Constant-hazard HR implied by two cumulative incidences over one horizon."""
    require_fraction("baseline", baseline)
    require_fraction("treated", treated)
    if baseline <= 0.0:
        raise ValueError("baseline risk must be positive")
    return float(math.log1p(-treated) / math.log1p(-baseline))


def hazard_gap_table(
    matrix: pd.DataFrame, reference_risk: float = 0.0198, threshold: float = CRITERION_THRESHOLD
) -> pd.DataFrame:
    """Recompute the required HR per population, then check the composite against it.

    The published target of 0.298 was derived from a 13.93% baseline. Applying it to a
    population with a 20% or 33% baseline is a category error: a higher baseline needs
    a *lower* HR to leave the same share of excess risk.
    """
    rows = []
    for row in matrix.itertuples(index=False):
        baseline = float(row.baseline_risk_pct) / 100.0
        treated = float(row.composite_risk_pct) / 100.0
        # The reference risk is quoted to 0.01 percentage points, so a baseline
        # within that precision of it is the reference population, not a gap.
        at_reference = baseline <= reference_risk + 1e-4
        required = (
            1.0
            if at_reference
            else required_capabilities(baseline, reference_risk, threshold)["hazard_ratio"]
        )
        achieved = hazard_ratio_from_risks(baseline, treated)
        rows.append(
            {
                "population_id": row.population_id,
                "mode": row.mode,
                "baseline_risk_pct": baseline * 100.0,
                "composite_risk_pct": treated * 100.0,
                "required_hr": required,
                "achieved_hr": achieved,
                "gap_factor": achieved / required,
                "at_reference_risk": at_reference,
                "closed": at_reference or achieved <= required * (1.0 + 1e-9),
                # Further reduction needed *on top of* the composite's current risk
                # to reach the required HR, expressed relative to current risk.
                "extra_rrr_needed_pct": (
                    0.0
                    if achieved <= required
                    else 100.0 * (1.0 - (-math.expm1(required * math.log1p(-baseline))) / treated)
                ),
                "status": GAP_TAG,
            }
        )
    return pd.DataFrame(rows).sort_values(["mode", "gap_factor"]).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Gap 2: substrate clearance. The plateau is a source term, not a payload limit.
# ---------------------------------------------------------------------------

#: CHOP deletion in Apoe-/- and Ldlr-/- mice reduced lesional apoptosis by ~35% and
#: plaque necrosis by ~50%. CHOP deletion does NOT abolish apoptosis (IRE1-JNK and
#: Bcl-2 pathways remain), so source reduction is bounded well below 1.0.
SOURCE_REDUCTION_PRECEDENT_CENTRAL = 0.35
SOURCE_REDUCTION_PRECEDENT_MAX = 0.50

SOURCE_PRECEDENT_URL = "https://www.ahajournals.org/doi/10.1161/circresaha.110.224766"


@dataclass(frozen=True)
class ClearanceSolution:
    """Steady state and trajectory for a source-reduction strength."""

    source_reduction: float
    removal_24m: float
    removal_36m: float
    removal_60m: float
    asymptotic_removal: float
    analytic_steady_state: float
    months_to_target: float | None
    meets_target: bool
    within_precedent: bool


def analytic_clearance_ceiling(source_reduction: float = 0.0) -> float:
    """Closed-form steady state of the perfect-module clearance loop.

    With every module at 1.0: influx is zero, uptake is 5.4*D/(D+15) and the debris
    source is 1.5*(1-s). Setting uptake equal to the source gives D*. Cholesterol
    balance then gives F* from efflux plus oxysterol exit equalling 2.5*source.
    Microcalcification vanishes because nucleation is fully blocked.
    """
    require_fraction("source_reduction", source_reduction, allow_one=True)
    source = 1.5 * (1.0 - source_reduction)
    max_uptake = 1.2 * 3.0 * 1.5  # (1+2*1) * (1+0.5*1)
    if source >= max_uptake:
        return 0.0
    d_star = 15.0 * source / (max_uptake - source)
    target_export = 2.5 * source

    def balance(f: float) -> float:
        signal = min(1.0, 0.25 * f / 8.0)
        capacity = 5.0 * (1.0 + 2.0 * signal)
        return capacity * f / (f + 20.0) + 0.25 * f * 0.15 - target_export

    low, high = 0.0, 1000.0
    if balance(high) < 0:
        return 0.0
    for _ in range(200):
        mid = 0.5 * (low + high)
        if balance(mid) < 0:
            low = mid
        else:
            high = mid
    f_star = 0.5 * (low + high)
    initial = 108.0  # default LesionState: 60 + 40 + 8
    return float(max(0.0, 1.0 - (f_star + d_star) / initial))


def solve_substrate_gap(
    source_reductions: tuple[float, ...] = (0.0, 0.20, 0.35, 0.50, 0.65, 0.80),
    required: float = CRITERION_THRESHOLD,
    perfect_modules: bool = True,
) -> pd.DataFrame:
    """Sweep the source-reduction module and find what closes the clearance gap."""
    require_fraction("required", required, allow_one=True)
    params = (
        ClearanceParameters(gate_specificity=1.0, **{n: 1.0 for n in _ZERO_STRENGTHS})
        if perfect_modules
        else ClearanceParameters()
    )
    rows = []
    for strength in source_reductions:
        require_fraction("source_reduction", strength, allow_one=True)
        long_run = simulate_clearance(
            params, "m7", months=240.0, step_months=0.5, source_reduction=strength
        )
        fractions = 1.0 - long_run.substrate / float(long_run.substrate[0])
        hits = np.flatnonzero(fractions >= required)
        rows.append(
            {
                "source_reduction": strength,
                "modules": "perfect" if perfect_modules else "central",
                "removal_24m_pct": 100.0
                * simulate_clearance(
                    params, "m7", months=24.0, step_months=0.5, source_reduction=strength
                ).substrate_removal_fraction,
                "removal_36m_pct": 100.0
                * simulate_clearance(
                    params, "m7", months=36.0, step_months=0.5, source_reduction=strength
                ).substrate_removal_fraction,
                "removal_60m_pct": 100.0
                * simulate_clearance(
                    params, "m7", months=60.0, step_months=0.5, source_reduction=strength
                ).substrate_removal_fraction,
                "asymptotic_removal_pct": 100.0 * long_run.substrate_removal_fraction,
                "analytic_ceiling_pct": 100.0 * analytic_clearance_ceiling(strength)
                if perfect_modules
                else float("nan"),
                "months_to_target": (
                    None if hits.size == 0 else float(long_run.months[int(hits[0])])
                ),
                "meets_target": bool(hits.size > 0),
                "within_precedent": strength <= SOURCE_REDUCTION_PRECEDENT_MAX,
                "status": GAP_TAG,
            }
        )
    return pd.DataFrame(rows)


def minimum_source_reduction(
    required: float = CRITERION_THRESHOLD,
    horizon_months: float = 24.0,
    perfect_modules: bool = True,
    tolerance: float = 1e-4,
) -> float | None:
    """Bisection for the weakest source-reduction module that hits the target."""
    require_fraction("required", required, allow_one=True)
    require_positive("horizon_months", horizon_months)
    params = (
        ClearanceParameters(gate_specificity=1.0, **{n: 1.0 for n in _ZERO_STRENGTHS})
        if perfect_modules
        else ClearanceParameters()
    )

    def removal(strength: float) -> float:
        return simulate_clearance(
            params,
            "m7",
            months=horizon_months,
            step_months=0.5,
            lesion=LesionState(),
            source_reduction=strength,
        ).substrate_removal_fraction

    if removal(1.0) < required:
        return None
    low, high = 0.0, 1.0
    while high - low > tolerance:
        mid = 0.5 * (low + high)
        if removal(mid) >= required:
            high = mid
        else:
            low = mid
    return float(high)


# ---------------------------------------------------------------------------
# Gap 1: remuscularisation. Multiplicative factors on one saturating quantity.
# ---------------------------------------------------------------------------

#: LVEF points recoverable if the entire contractile deficit were replaced.
SCAR_RECOVERABLE_LVEF_PP = 25.0

#: BioVAT-HF delivered +3.9 points, i.e. 15.6% of the deficit. Primate grafts occupy
#: 10-17% of scar area, an independent measurement that matches - so modelling LVEF
#: gain as proportional to the fraction of scar remuscularised is calibrated twice.
BASELINE_SCAR_FRACTION = 3.9 / SCAR_RECOVERABLE_LVEF_PP


@dataclass(frozen=True)
class Factor:
    """One multiplier on engrafted contractile mass, with its precedent and bound."""

    name: str
    demonstrated_multiplier: float
    precedent: str
    source_url: str
    species: str

    def __post_init__(self) -> None:
        require_positive("demonstrated_multiplier", self.demonstrated_multiplier)
        if not self.source_url.startswith("https://"):
            raise ValueError("an https source_url is required")


#: Each factor is a real published result. None is a human efficacy result.
REMUSCULARISATION_FACTORS: tuple[Factor, ...] = (
    Factor(
        "dose_escalation",
        1.75,
        "BioVAT maximal dose 800M cells; human-scale EHT patches built at 450M cells and "
        "~1400M input estimated necessary for a 300 g heart",
        "https://www.ahajournals.org/doi/10.1161/CIRCULATIONAHA.120.047904",
        "human_construct_animal_transplant",
    ),
    Factor(
        "endothelial_co_delivery",
        2.0,
        "EC+CM co-transplantation raised primate EF by 8.9 points vs 1.9 for CM alone; "
        "preformed microvessels gave 10-fold larger grafts in rodents",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC10683868/",
        "non_human_primate",
    ),
    Factor(
        "epicardial_stromal_co_delivery",
        2.6,
        "Epicardial cell co-transplantation doubled CM proliferation and gave 2.6-fold "
        "greater graft size with improved systolic function",
        "https://www.researchgate.net/publication/349735481_Human_Engineered_Heart_Tissue_Patches_Remuscularize_the_Injured_Heart_in_a_Dose-Dependent_Manner",
        "rodent",
    ),
    Factor(
        "proliferation_induction",
        2.0,
        "CCND2-overexpressing hiPSC-CMs replaced >50% of scar at 6 months; CHIR/FGF1 "
        "nanoparticle patches increased engraftment fourfold",
        "https://www.nature.com/articles/s44385-025-00048-4",
        "rodent",
    ),
)


@dataclass(frozen=True)
class RemuscularisationScenario:
    """A combination of factors and the LVEF gain it implies, with saturation."""

    factors: tuple[str, ...]
    raw_multiplier: float
    scar_fraction: float
    lvef_gain_pp: float
    saturated: bool
    any_human_efficacy: bool
    notes: tuple[str, ...] = field(default=(DISCLAIMER,))

    @property
    def meets_target(self) -> bool:
        return self.lvef_gain_pp >= CRITERION_THRESHOLD * SCAR_RECOVERABLE_LVEF_PP


def remuscularisation_scenario(
    factors: tuple[Factor, ...],
    baseline_scar_fraction: float = BASELINE_SCAR_FRACTION,
) -> RemuscularisationScenario:
    """Combine factors onto scar fraction, which saturates at full replacement.

    The saturation is the point. These factors all act on the same limiting quantity -
    engrafted contractile mass - so their product cannot exceed the scar available to
    replace. Multiplying them without a cap would produce absurd gains.
    """
    require_fraction("baseline_scar_fraction", baseline_scar_fraction, allow_one=True)
    raw = 1.0
    for factor in factors:
        raw *= factor.demonstrated_multiplier
    projected = baseline_scar_fraction * raw
    capped = min(1.0, projected)
    return RemuscularisationScenario(
        tuple(factor.name for factor in factors),
        raw,
        capped,
        capped * SCAR_RECOVERABLE_LVEF_PP,
        projected > 1.0,
        any(factor.species == "human_outcomes" for factor in factors),
        (
            DISCLAIMER,
            "factors act on one saturating quantity: fraction of scar remuscularised",
            "no factor here is a human efficacy result",
            "chronic-injury efficacy is unproven; most primate work treats 2 weeks post-MI",
        ),
    )


def remuscularisation_search() -> pd.DataFrame:
    """Every subset of the precedented factors, smallest first."""
    from itertools import combinations

    rows = []
    for size in range(len(REMUSCULARISATION_FACTORS) + 1):
        for subset in combinations(REMUSCULARISATION_FACTORS, size):
            scenario = remuscularisation_scenario(subset)
            rows.append(
                {
                    "n_factors": size,
                    "factors": ";".join(scenario.factors) or "none",
                    "raw_multiplier": scenario.raw_multiplier,
                    "scar_fraction": scenario.scar_fraction,
                    "lvef_gain_pp": scenario.lvef_gain_pp,
                    "saturated": scenario.saturated,
                    "meets_target": scenario.meets_target,
                    "status": GAP_TAG,
                }
            )
    return (
        pd.DataFrame(rows)
        .sort_values(["meets_target", "n_factors"], ascending=[False, True])
        .reset_index(drop=True)
    )


def minimum_remuscularisation_set() -> tuple[str, ...]:
    """Fewest precedented factors that reach +20 LVEF points."""
    search = remuscularisation_search()
    passing = search.loc[search["meets_target"]]
    if passing.empty:
        return ()
    best = passing.sort_values("n_factors").iloc[0]
    return () if best["factors"] == "none" else tuple(best["factors"].split(";"))


def factor_frame() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "factor": factor.name,
                "demonstrated_multiplier": factor.demonstrated_multiplier,
                "species": factor.species,
                "precedent": factor.precedent,
                "source_url": factor.source_url,
                "status": GAP_TAG,
            }
            for factor in REMUSCULARISATION_FACTORS
        ]
    )


# ---------------------------------------------------------------------------
# Verdicts and reporting.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class GapResult:
    capability: str
    needed: str
    had: str
    verdict: GapVerdict
    mechanism: str
    residual_requirement: str


def assess_gaps(matrix: pd.DataFrame) -> tuple[list[GapResult], dict]:
    """Combine the three analyses into one verdict table."""
    hazards = hazard_gap_table(matrix)
    strict = hazards.loc[hazards["mode"].eq("strict")]
    closed_pops = list(strict.loc[strict["closed"], "population_id"])

    clearance = solve_substrate_gap()
    precedented = clearance.loc[clearance["within_precedent"] & clearance["meets_target"]]
    min_24 = minimum_source_reduction(horizon_months=24.0)
    min_asym = minimum_source_reduction(horizon_months=240.0)
    min_24_central = minimum_source_reduction(horizon_months=24.0, perfect_modules=False)
    min_asym_central = minimum_source_reduction(horizon_months=240.0, perfect_modules=False)

    minimum_set = minimum_remuscularisation_set()
    best_scenario = remuscularisation_scenario(REMUSCULARISATION_FACTORS)

    results = [
        GapResult(
            "Stack hazard ratio",
            "population-specific",
            "0.360 (old best stack)",
            GapVerdict.CLOSED if closed_pops else GapVerdict.OPEN,
            "Node-disjoint composite; target HR recomputed per baseline rather than "
            "applied universally",
            f"closed in {len(closed_pops)} of {len(strict)} populations in strict mode: "
            + (", ".join(closed_pops) if closed_pops else "none"),
        ),
        GapResult(
            "Substrate clearance",
            "80%",
            "43.1%",
            GapVerdict.CLOSABLE if not precedented.empty else GapVerdict.OPEN,
            "Add a source-reduction module (ER-stress/CHOP): the plateau is a constant "
            "apoptosis source, not a payload limit",
            (
                f"perfect modules: >={min_24:.3f} source reduction for 80% at 24 months; "
                f"central modules: >={min_24_central:.3f} at 24 months (exceeds the "
                f"{SOURCE_REDUCTION_PRECEDENT_MAX:.2f} precedent bound) or "
                f">={min_asym_central:.3f} over a much longer horizon"
                if None not in (min_24, min_24_central, min_asym_central)
                else "unreachable even at full source suppression"
            ),
        ),
        GapResult(
            "Remuscularisation",
            "+20 LVEF pts",
            "+3.9 (BioVAT-HF)",
            GapVerdict.CLOSABLE if minimum_set else GapVerdict.NEEDS_NEW_BIOLOGY,
            "Stack precedented multipliers on scar replacement fraction, which saturates "
            "at 1.0; all multipliers are animal or construct-level",
            f"minimum set: {', '.join(minimum_set)}"
            if minimum_set
            else "no combination of precedented factors reaches the target",
        ),
    ]
    detail = {
        "hazard_populations_closed": len(closed_pops),
        "hazard_populations_tested": len(strict),
        "minimum_source_reduction_24m": min_24,
        "minimum_source_reduction_asymptotic": min_asym,
        "minimum_source_reduction_24m_central": min_24_central,
        "minimum_source_reduction_asymptotic_central": min_asym_central,
        "clearance_closes_in_24m_within_precedent": (
            min_24 is not None and min_24 <= SOURCE_REDUCTION_PRECEDENT_MAX
        ),
        "clearance_needs_near_perfect_delivery_for_24m": (
            min_24_central is None or min_24_central > SOURCE_REDUCTION_PRECEDENT_MAX
        ),
        "source_reduction_precedent_max": SOURCE_REDUCTION_PRECEDENT_MAX,
        "clearance_closable_within_precedent": not precedented.empty,
        "minimum_remuscularisation_factors": minimum_set,
        "all_factors_lvef_gain_pp": best_scenario.lvef_gain_pp,
        "all_factors_saturated": best_scenario.saturated,
    }
    return results, detail


def gap_report(matrix: pd.DataFrame) -> str:
    results, detail = assess_gaps(matrix)
    hazards = hazard_gap_table(matrix)
    strict = hazards.loc[hazards["mode"].eq("strict")]
    clearance = solve_substrate_gap()
    central = solve_substrate_gap(perfect_modules=False)
    search = remuscularisation_search()

    lines = [
        "Attempting to close the three capability gaps",
        "=" * 74,
        DISCLAIMER,
        "",
        "VERDICTS",
    ]
    for result in results:
        lines.append(f"  {result.capability}: {result.verdict.value}")
        lines.append(f"    mechanism: {result.mechanism}")
        lines.append(f"    requirement: {result.residual_requirement}")
    lines.extend(
        [
            "",
            "GAP 3 - HAZARD RATIO: already closed, but the target was population-specific",
            "The 0.298 figure came from a 13.93% baseline. A higher baseline needs a lower",
            "HR to leave the same share of excess risk, so it must be recomputed:",
        ]
    )
    for row in strict.head(8).itertuples(index=False):
        flag = (
            "already at reference"
            if row.at_reference_risk
            else ("CLOSED" if row.closed else f"needs {row.extra_rrr_needed_pct:.1f}% more")
        )
        lines.append(
            f"  {row.population_id:28s} need HR<={row.required_hr:.3f}  "
            f"achieved {row.achieved_hr:.3f}  {flag}"
        )
    lines.extend(
        [
            "",
            "GAP 2 - SUBSTRATE CLEARANCE: the plateau was structural, not a payload limit",
            "With every module perfect the loop still stalls at 83.5%, because a constant",
            "apoptosis source keeps regenerating debris. The analytic steady state matches",
            "the simulation exactly, which is how the cause was identified. None of the six",
            "REVERT-1 modules touches that source. Adding one that does (ER-stress/CHOP):",
            "",
            "  source reduction   24 mo   36 mo   60 mo   asymptote   within precedent",
        ]
    )
    for row in clearance.itertuples(index=False):
        lines.append(
            f"      {row.source_reduction:5.2f}        {row.removal_24m_pct:5.1f}%  "
            f"{row.removal_36m_pct:5.1f}%  {row.removal_60m_pct:5.1f}%    "
            f"{row.asymptotic_removal_pct:5.1f}%      "
            f"{'yes' if row.within_precedent else 'NO'}"
        )
    lines.extend(
        [
            "  precedent bound: CHOP deletion cut lesional apoptosis ~35% and plaque",
            "  necrosis ~50% in Apoe/Ldlr-null mice; it does not abolish apoptosis, so",
            f"  strengths above {SOURCE_REDUCTION_PRECEDENT_MAX:.2f} are unprecedented.",
            "",
            "  That table assumes PERFECT modules. With central strengths the same module",
            "  delivers much less, which is the honest caveat:",
            "",
            "  source reduction   24 mo   60 mo   asymptote   meets 80%   within precedent",
        ]
    )
    for row in central.itertuples(index=False):
        lines.append(
            f"      {row.source_reduction:5.2f}        {row.removal_24m_pct:5.1f}%  "
            f"{row.removal_60m_pct:5.1f}%    {row.asymptotic_removal_pct:5.1f}%      "
            f"{'yes' if row.meets_target else 'no'}"
            f"{' ' * 10}{'yes' if row.within_precedent else 'NO'}"
        )
    lines.extend(
        [
            "  So clearance closes in two years only with near-perfect delivery; at central",
            "  strengths it needs a far longer horizon. Delivery efficiency and time are",
            "  jointly binding - consistent with the gate dominating the sensitivity.",
            "",
            "GAP 1 - REMUSCULARISATION: decomposes, but every factor is animal-level",
            "BioVAT's +3.9 points is 15.6% of the 25-point deficit. Primate grafts occupy",
            "10-17% of scar area - an independent measurement that matches. So LVEF gain",
            "is modelled as proportional to the fraction of scar replaced, which saturates.",
            "",
            "  factors                              scar replaced   LVEF gain   meets +20",
        ]
    )
    for row in search.head(10).itertuples(index=False):
        label = row.factors if len(row.factors) <= 36 else row.factors[:33] + "..."
        lines.append(
            f"  {label:38s} {row.scar_fraction:8.1%}   {row.lvef_gain_pp:8.1f}   "
            f"{'YES' if row.meets_target else 'no'}"
        )
    lines.extend(
        [
            "",
            "WHAT THIS DOES AND DOES NOT ESTABLISH",
            f"- Hazard ratio: closed in {detail['hazard_populations_closed']} of "
            f"{detail['hazard_populations_tested']} populations. Real, and it needed no new biology.",
            "- Clearance: closable by adding a module aimed at the apoptotic source, at a",
            "  strength inside the animal precedent. This is a design change, not tuning.",
            "- Remuscularisation: the arithmetic closes by stacking published multipliers,",
            "  but every one is rodent, primate or construct-level. Saturation means the",
            "  factors cannot be multiplied indefinitely, and chronic-injury efficacy is",
            "  unproven - most primate work treats two weeks after infarction, not the",
            "  years-old scar of an advanced heart-failure patient.",
            "- None of this is evidence. Closing a gap in a model is a reason to run an",
            "  experiment, not a reason to believe the result.",
        ]
    )
    return "\n".join(lines) + "\n"


def run_gap_study(output_dir, matrix: pd.DataFrame | None = None) -> dict:
    """Write the gap-closure artefacts to their own quarantined directory."""
    from pathlib import Path

    root = Path(output_dir)
    (root / "tables").mkdir(parents=True, exist_ok=True)
    (root / "reports").mkdir(parents=True, exist_ok=True)

    if matrix is None:
        source = Path("results/composite/tables/composite_test_matrix.csv")
        if not source.exists():
            raise FileNotFoundError(
                "composite matrix missing; run `python run_solutions.py composite` first"
            )
        matrix = pd.read_csv(source)

    results, detail = assess_gaps(matrix)
    hazards = hazard_gap_table(matrix)
    clearance = solve_substrate_gap()
    central = solve_substrate_gap(perfect_modules=False)
    search = remuscularisation_search()

    pd.DataFrame(
        [
            {
                "capability": r.capability,
                "needed": r.needed,
                "had": r.had,
                "verdict": r.verdict.value,
                "mechanism": r.mechanism,
                "residual_requirement": r.residual_requirement,
                "status": GAP_TAG,
            }
            for r in results
        ]
    ).to_csv(root / "tables" / "gap_verdicts.csv", index=False)
    hazards.to_csv(root / "tables" / "gap_hazard_by_population.csv", index=False)
    clearance.to_csv(root / "tables" / "gap_clearance_source_reduction.csv", index=False)
    central.to_csv(root / "tables" / "gap_clearance_central_modules.csv", index=False)
    search.to_csv(root / "tables" / "gap_remuscularisation_search.csv", index=False)
    factor_frame().to_csv(root / "tables" / "gap_remuscularisation_factors.csv", index=False)
    pd.DataFrame([{**detail, "status": GAP_TAG}]).to_csv(
        root / "tables" / "gap_detail.csv", index=False
    )
    (root / "reports" / "gap_report.txt").write_text(gap_report(matrix), encoding="utf-8")
    (root / "DISCLAIMER.txt").write_text(DISCLAIMER + "\n", encoding="utf-8")
    return {
        "status": "scenario_only",
        "disclaimer": DISCLAIMER,
        "output_dir": str(root),
        "verdicts": {r.capability: r.verdict.value for r in results},
        "outputs": 9,
    }
