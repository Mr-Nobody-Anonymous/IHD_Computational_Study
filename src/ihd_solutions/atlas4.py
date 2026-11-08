"""ATLAS-4: the graft model corrected against external precedent, then searched
exhaustively.

A literature search (recorded in :mod:`ihd_solutions.precedent`) found that the
architecture ATLAS-3 called invented is largely already precedented, and that two of
its assumptions were wrong:

* **Attempt limit.** Four doses has human precedent (CATO, HOPE-2), not two. ATLAS-3's
  three attempts were never beyond precedent.
* **Retry operator.** ATLAS-3 used best-of-k. Every preclinical repeat-dose study
  reports *cumulative* increments instead, and shows that splitting a fixed total cell
  number across three doses beats giving it all at once. Best-of-k cannot represent
  that finding.

This module implements both operators so they can be contrasted on the same draws,
gates cumulative gain by remaining headroom (the ceiling effect measured in the one
blinded human repeat trial), and tempers repeat responsiveness by the human responder
fraction rather than by rodent enthusiasm.

It then searches the **whole** design grid and reports every point, and audits what
each of the five cure criteria can actually reach. It does not tune anything toward a
target: a design that reads 100% is reported as a saturation artefact and excluded
from the recommendation, exactly as in ATLAS-3.

NOTHING HERE IS A THERAPY OR A RESULT. These are simulations of unvalidated
architectures, anchored where possible to other people's published measurements.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from enum import Enum
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd

from ihd_solutions.atlas3 import (
    DEFAULT_PER_ATTEMPT_RISK,
    _factor_draws,
)
from ihd_solutions.cure import CRITERION_THRESHOLD
from ihd_solutions.gap_solver import (
    BASELINE_SCAR_FRACTION,
    REMUSCULARISATION_FACTORS,
    SCAR_RECOVERABLE_LVEF_PP,
)
from ihd_solutions.precedent import (
    PRECEDENTS,
    cumulative_increment_ratios,
    human_repeat_efficacy_verdict,
    max_precedented_doses,
    measurement_quality_anchor,
    precedent_frame,
    retry_responder_fraction,
)
from ihd_solutions.units import require_fraction, require_positive

ATLAS4_TAG = "atlas4_precedent_corrected_unvalidated"

DISCLAIMER = (
    "ATLAS-4 IS A SIMULATION OF UNVALIDATED ARCHITECTURES. Repeat cell dosing is "
    "precedented and feasible in humans, but its efficacy is not established: the one "
    "randomised blinded repeat-dose trial with a functional primary endpoint was null. "
    "The cumulative-dosing mechanism is anchored to rodent and porcine work, and the "
    "cells involved do not regenerate myocardium. Nothing here is a therapy, a result, "
    "or a recommendation, and no cure is declared."
)

#: Per-dose effect retained when a fixed total dose is split. Anchored to the rat CPC
#: study: three doses of 12M gave +13.0 pp cumulative, one dose of 36M gave +4.7 pp, so
#: each third-dose delivered (13.0/3)/4.7 = 0.922 of a full dose's effect. Duration of
#: exposure matters more than intensity. Rodent-derived; swept, never asserted.
SPLIT_RETENTION = (13.0 / 3.0) / 4.7

#: Measurement quality that carries no viability information. Iron-oxide MRI localises
#: cells but the signal persists in macrophages after they die, and acute retention
#: does not correlate with engraftment. A feedback loop built on either is blind.
NON_VIABILITY_SIGNAL_QUALITY = 0.0


class DosingMode(str, Enum):
    """How successive doses combine."""

    BEST_OF = "best_of"
    CUMULATIVE = "cumulative"


@dataclass(frozen=True)
class GraftDesign:
    """A point in the design grid."""

    doses: int = 4
    mode: DosingMode = DosingMode.CUMULATIVE
    failure_correlation: float = 0.0
    feedback_learning: float = 0.5
    measurement_quality: float = 0.78
    headroom_gating: bool = True
    repeat_responder_probability: float = 1.0 / 3.0
    per_dose_risk: float = DEFAULT_PER_ATTEMPT_RISK
    translation_floor: float = 1.0

    def __post_init__(self) -> None:
        if self.doses < 1:
            raise ValueError("at least one dose is required")
        if not isinstance(self.mode, DosingMode):
            raise ValueError("mode must be a DosingMode")
        require_fraction("failure_correlation", self.failure_correlation, allow_one=True)
        require_fraction("feedback_learning", self.feedback_learning, allow_one=True)
        require_fraction("measurement_quality", self.measurement_quality, allow_one=True)
        require_fraction(
            "repeat_responder_probability", self.repeat_responder_probability, allow_one=True
        )
        require_fraction("per_dose_risk", self.per_dose_risk, allow_one=True)
        require_fraction("translation_floor", self.translation_floor, allow_one=True)

    @property
    def beyond_precedent(self) -> bool:
        return self.doses > max_precedented_doses()


#: Human-tempered default: repeat responsiveness taken from the blinded human trial
#: (10 of 30 responded to dose 2), measurement quality from the published PET
#: correlation midpoint, four doses because four has human precedent.
ATLAS4 = GraftDesign(
    doses=max_precedented_doses(),
    mode=DosingMode.CUMULATIVE,
    feedback_learning=0.5,
    measurement_quality=float(np.mean(measurement_quality_anchor())),
    headroom_gating=True,
    repeat_responder_probability=retry_responder_fraction(),
    translation_floor=1.0,
)

#: The same architecture with rodent optimism instead of human tempering. Reported
#: beside the default so the size of the translation assumption is visible.
ANIMAL_OPTIMISTIC = replace(ATLAS4, repeat_responder_probability=1.0)


@dataclass(frozen=True)
class DesignResult:
    probability_of_target: float
    median_gain_pp: float
    p05_gain_pp: float
    mean_doses_used: float
    cumulative_adverse_risk: float
    saturated_fraction: float
    doses: int
    mode: str
    beyond_precedent: bool
    draws: int

    @property
    def robust(self) -> bool:
        return self.probability_of_target >= CRITERION_THRESHOLD

    @property
    def monte_carlo_standard_error(self) -> float:
        p = self.probability_of_target
        return float(np.sqrt(max(p * (1.0 - p), 0.0) / max(self.draws, 1)))

    @property
    def robust_with_margin(self) -> bool:
        return self.probability_of_target - 2.0 * self.monte_carlo_standard_error >= (
            CRITERION_THRESHOLD
        )

    @property
    def risk_adjusted_probability(self) -> float:
        return self.probability_of_target * (1.0 - self.cumulative_adverse_risk)

    @property
    def saturation_artefact(self) -> bool:
        """A probability of 1 is a saturation artefact, never a result.

        Mirrors the ATLAS-3 lesson: when the model can make success certain, the
        model is wrong or the operating point is degenerate.
        """
        return self.probability_of_target >= 1.0


def design_draws(
    design: GraftDesign = ATLAS4, draws: int = 4000, seed: int = 20261107
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (gain, doses used, first-dose gain) on shared random numbers.

    The first-dose array makes every contrast against single dosing paired, so a
    lever that should contribute nothing contributes exactly zero.
    """
    require_positive("draws", draws)
    count = int(draws)
    rng = np.random.default_rng(seed)
    shared = rng.random((count, len(REMUSCULARISATION_FACTORS)))
    ceiling = SCAR_RECOVERABLE_LVEF_PP
    required = CRITERION_THRESHOLD * ceiling
    shape = cumulative_increment_ratios()

    gain = np.zeros(count)
    first = np.zeros(count)
    used = np.ones(count)
    done = np.zeros(count, dtype=bool)

    for index in range(design.doses):
        delivered = _factor_draws(rng, count, index, _upgrade_view(design), shared)
        penalty = rng.uniform(design.translation_floor, 1.0, size=count)
        increment = np.minimum(1.0, BASELINE_SCAR_FRACTION * delivered * penalty) * ceiling
        if index == 0:
            first = increment.copy()
            gain = increment.copy()
            done = gain >= required
            continue

        # Only a fraction of patients respond to a repeat dose at all. Anchored to the
        # blinded human trial rather than to the rodent cumulative curves.
        responds = rng.random(count) < design.repeat_responder_probability
        active = (~done) & responds
        if design.mode is DosingMode.BEST_OF:
            candidate = np.maximum(gain, increment)
        else:
            ratio = shape[min(index, len(shape) - 1)] / shape[0]
            scaled = increment * SPLIT_RETENTION * ratio
            if design.headroom_gating:
                # The measured human ceiling effect: a patient already near the
                # recoverable ceiling has little left to gain from another dose.
                scaled = scaled * np.clip((ceiling - gain) / ceiling, 0.0, 1.0)
            candidate = np.minimum(ceiling, gain + scaled)
        gain = np.where(active, candidate, gain)
        used = np.where(active, float(index + 1), used)
        done = done | (gain >= required)

    return gain, used, first


def _upgrade_view(design: GraftDesign):
    """Adapt a design to the ATLAS-3 factor-draw signature.

    Reuses ATLAS-3's per-factor uncertainty and bounded aggregate feedback so the two
    studies differ only in the operator, the headroom gate and the human tempering.
    """
    from ihd_solutions.atlas3 import GraftUpgrades

    return GraftUpgrades(
        attempts=design.doses,
        failure_correlation=design.failure_correlation,
        feedback_learning=design.feedback_learning,
        measurement_quality=design.measurement_quality,
        per_attempt_risk=design.per_dose_risk,
        translation_floor=design.translation_floor,
    )


def evaluate_design(
    design: GraftDesign = ATLAS4, draws: int = 4000, seed: int = 20261107
) -> DesignResult:
    gain, used, _ = design_draws(design, draws, seed)
    required = CRITERION_THRESHOLD * SCAR_RECOVERABLE_LVEF_PP
    mean_used = float(np.mean(used))
    return DesignResult(
        probability_of_target=float(np.mean(gain >= required)),
        median_gain_pp=float(np.median(gain)),
        p05_gain_pp=float(np.quantile(gain, 0.05)),
        mean_doses_used=mean_used,
        cumulative_adverse_risk=float(1.0 - (1.0 - design.per_dose_risk) ** mean_used),
        saturated_fraction=float(np.mean(gain >= SCAR_RECOVERABLE_LVEF_PP - 1e-9)),
        doses=design.doses,
        mode=design.mode.value,
        beyond_precedent=design.beyond_precedent,
        draws=int(draws),
    )


def operator_contrast(draws: int = 4000, seed: int = 20261107) -> pd.DataFrame:
    """Best-of-k against cumulative, on identical draws, across dose counts."""
    rows = []
    for doses, mode in product((1, 2, 3, 4, 6), DosingMode):
        design = replace(ATLAS4, doses=doses, mode=mode)
        result = evaluate_design(design, draws, seed)
        rows.append(
            {
                "doses": doses,
                "mode": mode.value,
                "probability_of_target": result.probability_of_target,
                "median_gain_pp": result.median_gain_pp,
                "mean_doses_used": result.mean_doses_used,
                "cumulative_adverse_risk": result.cumulative_adverse_risk,
                "risk_adjusted_probability": result.risk_adjusted_probability,
                "saturated_fraction": result.saturated_fraction,
                "robust": result.robust,
                "beyond_precedent": result.beyond_precedent,
                "status": ATLAS4_TAG,
            }
        )
    return pd.DataFrame(rows)


def correlation_by_operator(
    correlations: tuple[float, ...] = (0.0, 0.25, 0.5, 0.75, 1.0),
    draws: int = 3000,
    seed: int = 20261107,
) -> pd.DataFrame:
    """The independence requirement is a property of the OPERATOR, not of biology.

    Best-of-k needs independent attempts: correlated attempts are one attempt, so the
    gain collapses to zero. Cumulative dosing does not need independence in the same
    way, because a patient whose biology is consistently favourable accumulates
    consistently favourable increments. ATLAS-3's "must measure independence" conclusion
    was therefore operator-specific, and the operator was the wrong one.
    """
    rows = []
    for mode, rho in product(DosingMode, correlations):
        design = replace(ATLAS4, mode=mode, failure_correlation=rho)
        result = evaluate_design(design, draws, seed)
        single = evaluate_design(replace(design, doses=1), draws, seed)
        rows.append(
            {
                "mode": mode.value,
                "failure_correlation": rho,
                "probability_of_target": result.probability_of_target,
                "single_dose_probability": single.probability_of_target,
                "gain_vs_single_pp": 100.0
                * (result.probability_of_target - single.probability_of_target),
                "robust": result.robust,
                "status": ATLAS4_TAG,
            }
        )
    return pd.DataFrame(rows)


def exhaustive_design_search(
    doses: tuple[int, ...] = (1, 2, 3, 4, 6),
    correlations: tuple[float, ...] = (0.0, 0.5, 1.0),
    learning: tuple[float, ...] = (0.0, 0.5, 1.0),
    headroom: tuple[bool, ...] = (True, False),
    responder: tuple[float, ...] = (1.0 / 3.0, 1.0),
    draws: int = 800,
    seed: int = 20261107,
) -> pd.DataFrame:
    """Every combination in the grid, scored. Nothing is hidden or pre-filtered.

    Saturated points (P = 1) are kept in the table and flagged, because deleting them
    would hide how easily this model can be pushed to a meaningless certainty.
    """
    rows = []
    for mode, k, rho, learn, gate, respond in product(
        DosingMode, doses, correlations, learning, headroom, responder
    ):
        design = GraftDesign(
            doses=k,
            mode=mode,
            failure_correlation=rho,
            feedback_learning=learn,
            measurement_quality=ATLAS4.measurement_quality,
            headroom_gating=gate,
            repeat_responder_probability=respond,
            translation_floor=1.0,
        )
        result = evaluate_design(design, draws, seed)
        rows.append(
            {
                "mode": mode.value,
                "doses": k,
                "failure_correlation": rho,
                "feedback_learning": learn,
                "headroom_gating": gate,
                "repeat_responder_probability": respond,
                "probability_of_target": result.probability_of_target,
                "risk_adjusted_probability": result.risk_adjusted_probability,
                "cumulative_adverse_risk": result.cumulative_adverse_risk,
                "saturated_fraction": result.saturated_fraction,
                "robust": result.robust,
                "robust_with_margin": result.robust_with_margin,
                "beyond_precedent": result.beyond_precedent,
                "saturation_artefact": result.saturation_artefact,
                "status": ATLAS4_TAG,
            }
        )
    return pd.DataFrame(rows)


def search_summary(grid: pd.DataFrame) -> dict[str, object]:
    """What the exhaustive search found, including how much of it is admissible."""
    # A recommendation may not rest on a repeat-response assumption more optimistic
    # than the human data. Without this filter the "best" point is the rodent-optimistic
    # one sitting just under certainty, which is how a grid search launders an
    # assumption into a result.
    human_ceiling = retry_responder_fraction() + 1e-9
    admissible = grid.loc[
        (~grid["beyond_precedent"])
        & (~grid["saturation_artefact"])
        & grid["robust_with_margin"]
        & grid["headroom_gating"]
        & (grid["repeat_responder_probability"] <= human_ceiling)
    ]
    best = None
    if not admissible.empty:
        row = admissible.sort_values("risk_adjusted_probability", ascending=False).iloc[0]
        best = {
            "mode": row["mode"],
            "doses": int(row["doses"]),
            "failure_correlation": float(row["failure_correlation"]),
            "feedback_learning": float(row["feedback_learning"]),
            "repeat_responder_probability": float(row["repeat_responder_probability"]),
            "probability_of_target": float(row["probability_of_target"]),
            "risk_adjusted_probability": float(row["risk_adjusted_probability"]),
        }
    return {
        "configurations_tested": int(len(grid)),
        "robust": int(grid["robust"].sum()),
        "robust_with_margin": int(grid["robust_with_margin"].sum()),
        "within_precedent": int((~grid["beyond_precedent"]).sum()),
        "saturation_artefacts": int(grid["saturation_artefact"].sum()),
        "admissible": int(len(admissible)),
        "reaching_certainty": int((grid["probability_of_target"] >= 1.0).sum()),
        "best_admissible": best,
    }


def measurement_modality_falsifier(draws: int = 3000, seed: int = 20261107) -> pd.DataFrame:
    """Feedback is only worth what the signal measures.

    Reporter-gene PET tracks living cells and predicts later function (rho 0.76-0.80).
    Iron-oxide MRI localises cells but its signal survives their death, and acute
    retention does not correlate with engraftment. A loop built on either of the
    latter is blind, and must gain exactly nothing over no feedback at all.
    """
    rows = []
    modalities = (
        ("reporter_gene_pet", float(np.mean(measurement_quality_anchor())), True),
        ("iron_oxide_mri", NON_VIABILITY_SIGNAL_QUALITY, False),
        ("acute_retention_only", NON_VIABILITY_SIGNAL_QUALITY, False),
    )
    blind = evaluate_design(replace(ATLAS4, feedback_learning=0.0), draws, seed)
    for name, quality, viability_specific in modalities:
        result = evaluate_design(replace(ATLAS4, measurement_quality=quality), draws, seed)
        rows.append(
            {
                "modality": name,
                "measurement_quality": quality,
                "viability_specific": viability_specific,
                "probability_of_target": result.probability_of_target,
                "gain_vs_no_feedback_pp": 100.0
                * (result.probability_of_target - blind.probability_of_target),
                "status": ATLAS4_TAG,
            }
        )
    return pd.DataFrame(rows)


def translation_sensitivity(draws: int = 3000, seed: int = 20261107) -> pd.DataFrame:
    """How much of the result is the human tempering rather than the biology."""
    rows = []
    for respond, floor in product((retry_responder_fraction(), 0.5, 1.0), (1.0, 0.7, 0.4)):
        design = replace(ATLAS4, repeat_responder_probability=respond, translation_floor=floor)
        result = evaluate_design(design, draws, seed)
        rows.append(
            {
                "repeat_responder_probability": respond,
                "translation_floor": floor,
                "probability_of_target": result.probability_of_target,
                "median_gain_pp": result.median_gain_pp,
                "risk_adjusted_probability": result.risk_adjusted_probability,
                "robust": result.robust,
                "status": ATLAS4_TAG,
            }
        )
    return pd.DataFrame(rows)


def criterion_ceiling_audit() -> pd.DataFrame:
    """What each cure criterion can actually reach, and what stops it.

    This is the direct answer to "keep going until all of them are 100%". Some
    criteria already reach 1.0 for catalogued agents. The ones that do not are capped
    by explicit rules in this repository, not by an unexplored corner of the design
    space, so no amount of further searching raises them. The caps are honesty rules:
    removing them would produce the number without producing the evidence.
    """
    from ihd_solutions.cure import (
        MAX_SINGLE_AGENT_RISK_NORMALISATION,
        CureCriterion,
    )

    ceilings = _catalogue_ceilings()
    blockers = {
        CureCriterion.DURABLE_DRIVER_CONTROL.value: (
            "Reachable. One-time editing already scores 1.0 for a catalogued agent; "
            "the open item is a verified guide sequence and human durability follow-up."
        ),
        CureCriterion.BURDEN_REVERSAL.value: (
            "Capped by regression biology: partial plaque regression is limited to 0.6 "
            "for any single agent. ATLAS-2 reaches the burden target in simulation only "
            "by an invented 2.80x potency cap at 60 months."
        ),
        CureCriterion.MYOCARDIAL_RESTORATION.value: (
            "Capped at 0.5 for any myocardial claim without human remuscularisation "
            "evidence. ATLAS-4 raises the simulated probability but cannot raise the "
            "score, because the score is an evidence cap, not a probability."
        ),
        CureCriterion.RISK_NORMALISATION.value: (
            f"Capped at {MAX_SINGLE_AGENT_RISK_NORMALISATION:.1f} per single agent: no "
            "single therapy returns event risk to a low-risk reference. This is a "
            "property of the effect sizes, not a search failure."
        ),
        CureCriterion.RANDOMISED_PROOF.value: (
            "Reachable only for the exact regimen randomised. Proof reaches 1.0 for "
            "outcomes-tier catalogued agents, but does not transfer to an untested "
            "combination, so every invented architecture here scores 0."
        ),
    }
    rows = []
    for criterion in CureCriterion:
        name = criterion.value
        ceiling = ceilings.get(name, 0.0)
        row = {
            "criterion": name,
            "catalogue_ceiling": ceiling,
            "meets_threshold": ceiling >= CRITERION_THRESHOLD,
            "already_at_100pct": ceiling >= 1.0,
            "blocker": blockers[name],
            "status": ATLAS4_TAG,
        }
        rows.append(row)
    return pd.DataFrame(rows)


def _catalogue_ceilings() -> dict[str, float]:
    """Per-criterion maxima over the real catalogue, in the reference population."""
    from ihd_solutions.catalog import load_catalog
    from ihd_solutions.cure import assess_catalog, criterion_ceiling
    from ihd_solutions.evidence import load_anchors
    from ihd_solutions.population import load_populations

    catalog = load_catalog()
    anchors = load_anchors()
    populations = load_populations()
    population = populations["ischemic_cm_heart_failure"]
    days = np.arange(0.0, 365.25 * 2.0, 30.0)
    assessments = assess_catalog(catalog, anchors, population, days)
    return criterion_ceiling(assessments)


def _retained(interaction: pd.DataFrame, mode: DosingMode) -> float:
    """Fraction of the zero-correlation gain that survives full correlation."""
    rows = interaction.loc[interaction["mode"].eq(mode.value)]
    free = float(rows.loc[rows["failure_correlation"].eq(0.0)].iloc[0]["gain_vs_single_pp"])
    shared = float(rows.loc[rows["failure_correlation"].eq(1.0)].iloc[0]["gain_vs_single_pp"])
    if abs(free) < 1e-12:
        return 0.0
    return shared / free


def atlas4_report(
    corrected: DesignResult,
    animal: DesignResult,
    contrast: pd.DataFrame,
    grid: pd.DataFrame,
    summary: dict[str, object],
    modalities: pd.DataFrame,
    translation: pd.DataFrame,
    ceilings: pd.DataFrame,
    interaction: pd.DataFrame,
) -> str:
    single = contrast.loc[
        contrast["doses"].eq(1) & contrast["mode"].eq(DosingMode.CUMULATIVE.value)
    ].iloc[0]
    lines = [
        "ATLAS-4: what the literature already did, and what is actually still missing",
        "=" * 76,
        DISCLAIMER,
        "",
        "WHAT THE SEARCH CHANGED",
        f"  precedents recorded                {len(PRECEDENTS)}",
        f"  precedented dose limit             {max_precedented_doses()} (was 2 in ATLAS-3)",
        "  ATLAS-3's 3 attempts               WITHIN precedent, not beyond it",
        "  retry operator                     cumulative, not best-of-k",
        "  mechanism                      paracrine; these cells do not regenerate myocardium",
        "",
        "  " + human_repeat_efficacy_verdict(),
        "",
        "THE OPERATOR CONTRAST (same draws, human-tempered repeat response)",
        "  doses  best-of-k   cumulative   difference",
    ]
    for doses in sorted(contrast["doses"].unique()):
        best = contrast.loc[
            contrast["doses"].eq(doses) & contrast["mode"].eq(DosingMode.BEST_OF.value)
        ].iloc[0]
        cumul = contrast.loc[
            contrast["doses"].eq(doses) & contrast["mode"].eq(DosingMode.CUMULATIVE.value)
        ].iloc[0]
        lines.append(
            f"    {int(doses)}     {100 * best.probability_of_target:6.2f}%     "
            f"{100 * cumul.probability_of_target:6.2f}%    "
            f"{100 * (cumul.probability_of_target - best.probability_of_target):+6.2f} pp"
        )
    lines += [
        "",
        "RESULT",
        f"  single dose                        P {100 * single.probability_of_target:6.2f}%",
        f"  ATLAS-4 (human-tempered)           P {100 * corrected.probability_of_target:6.2f}%"
        f"   median {corrected.median_gain_pp:5.2f} pp",
        f"  robust with 2-SE margin            {'YES' if corrected.robust_with_margin else 'no'}",
        f"  beyond precedent                   {'YES' if corrected.beyond_precedent else 'no'}",
        f"  cumulative adverse risk            {100 * corrected.cumulative_adverse_risk:5.2f}%",
        f"  risk-adjusted success              {100 * corrected.risk_adjusted_probability:6.2f}%",
        "",
        f"  Same design with rodent optimism   P {100 * animal.probability_of_target:6.2f}%"
        f"   (saturated {100 * animal.saturated_fraction:.1f}%)",
        "  The gap between those two lines is the translation assumption, stated not buried.",
        "",
        "EXHAUSTIVE SEARCH: EVERY COMBINATION TESTED",
        f"  configurations                     {summary['configurations_tested']}",
        f"  robust                             {summary['robust']}",
        f"  robust with margin                 {summary['robust_with_margin']}",
        f"  within precedent                   {summary['within_precedent']}",
        f"  saturation artefacts (P = 100%)     {summary['saturation_artefacts']}",
        f"  admissible after all filters       {summary['admissible']}",
        "",
        "  Points that read 100% are reported and then EXCLUDED. A model that can be",
        "  made certain is a model that has stopped describing the problem.",
    ]
    best = summary["best_admissible"]
    if isinstance(best, dict):
        lines += [
            "",
            "  best admissible point:",
            f"    mode {best['mode']}, {best['doses']} doses, correlation "
            f"{best['failure_correlation']:.2f}, learning {best['feedback_learning']:.2f}, "
            f"repeat response {best['repeat_responder_probability']:.2f}",
            f"    P {100 * best['probability_of_target']:.2f}%,"
            f" risk-adjusted {100 * best['risk_adjusted_probability']:.2f}%",
        ]
    lines += [
        "",
        "THE INDEPENDENCE REQUIREMENT WAS OPERATOR-SPECIFIC",
        "  mode          correlation   P(target)   gain vs single dose",
    ]
    for row in interaction.itertuples(index=False):
        lines.append(
            f"  {row.mode:12s}    {row.failure_correlation:4.2f}      "
            f"{100 * row.probability_of_target:6.2f}%    {row.gain_vs_single_pp:+6.2f} pp"
        )
    lines += [
        "",
        f"  Gain retained at correlation 1.00: best-of-k "
        f"{100 * _retained(interaction, DosingMode.BEST_OF):.1f}%, "
        f"cumulative {100 * _retained(interaction, DosingMode.CUMULATIVE):.1f}%.",
        "",
        "  Best-of-k loses most of its benefit when doses are correlated, because",
        "  correlated attempts are close to one attempt; what survives is the feedback",
        "  uplift on later attempts, not the redundancy. Cumulative dosing keeps almost",
        "  all of it: a patient whose biology is consistently favourable accumulates",
        "  consistently favourable increments. ATLAS-3's headline 'measure independence",
        "  above 0.71' was a property of the wrong operator, not of the biology.",
        "",
        "THE MEASUREMENT FALSIFIER",
        "  modality                 quality   viability-specific   gain vs no feedback",
    ]
    for row in modalities.itertuples(index=False):
        lines.append(
            f"  {row.modality:24s} {row.measurement_quality:5.2f}   "
            f"{'yes' if row.viability_specific else 'no ':3s}                "
            f"{row.gain_vs_no_feedback_pp:+6.2f} pp"
        )
    lines += [
        "",
        "  Iron-oxide MRI and acute retention gain exactly nothing: iron signal survives",
        "  cell death and retention does not predict engraftment. Only a viability-specific",
        "  signal can close a feedback loop. That is a published finding, not a modelling",
        "  convenience.",
        "",
        "CAN EVERY CRITERION REACH 100%? MEASURED, NOT ASSERTED",
        "  criterion                    catalogue ceiling   at 100%   meets 0.8",
    ]
    for row in ceilings.itertuples(index=False):
        lines.append(
            f"  {row.criterion:28s} {row.catalogue_ceiling:6.3f}          "
            f"{'YES' if row.already_at_100pct else 'no ':3s}       "
            f"{'YES' if row.meets_threshold else 'no'}"
        )
    lines += [
        "",
        "  What blocks each one:",
    ]
    for row in ceilings.itertuples(index=False):
        lines.append(f"    {row.criterion}: {row.blocker}")
    lines += [
        "",
        "VERDICT ON 'KEEP GOING UNTIL ALL ARE 100%'",
        "That target is not reachable, and the reason is structural rather than a",
        "failure to search hard enough:",
        "",
        "  1. Two criteria already reach 1.0 for catalogued agents, so they were never",
        "     the binding constraint.",
        "  2. Myocardial restoration and single-agent risk normalisation are capped by",
        "     explicit evidence rules (0.5 and 0.4). Raising them means deleting the",
        "     rules, which would manufacture the number without producing evidence.",
        "  3. Randomised proof attaches to the exact regimen tested. It cannot be",
        "     inherited by an untested combination, so every architecture in this",
        "     repository scores zero on it no matter how well it simulates.",
        "  4. Any operating point that does read 100% in the graft model is a saturation",
        "     artefact and is excluded above.",
        "",
        "WHAT IS ACTUALLY STILL MISSING (and not already done by someone else)",
        "  Repeat dosing: DONE by others - CATO is randomising patients to four doses now.",
        "  Engraftment measurement: DONE in animals - reporter-gene PET predicts function.",
        "  Cumulative benefit: DONE in rats and pigs - and it is not regeneration.",
        "  Response-guided retry timing: PROPOSED in print in 2019, still NOT TESTED.",
        "  A viability-specific signal in humans: NOT DONE - the animal signal needs a",
        "    reporter gene the cells must be engineered to carry.",
        "",
        "So the genuinely open item is narrower than ATLAS-3 claimed: not 'repeat the",
        "dose', which is already in trial, but 'measure engraftment in a human well",
        "enough to decide when to repeat'. That is the experiment worth running.",
        "",
        "No cure is declared. Randomised proof for every architecture here remains zero.",
    ]
    return "\n".join(lines) + "\n"


def run_atlas4_study(
    output_dir: str | Path = "results/atlas4", draws: int = 3000, seed: int = 20261107
) -> dict[str, object]:
    root = Path(output_dir)
    tables = root / "tables"
    reports = root / "reports"
    tables.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)

    corrected = evaluate_design(ATLAS4, draws, seed)
    animal = evaluate_design(ANIMAL_OPTIMISTIC, draws, seed)
    contrast = operator_contrast(draws, seed)
    grid = exhaustive_design_search(draws=max(400, draws // 4), seed=seed)
    summary = search_summary(grid)
    modalities = measurement_modality_falsifier(draws, seed)
    translation = translation_sensitivity(draws, seed)
    ceilings = criterion_ceiling_audit()
    interaction = correlation_by_operator(draws=draws, seed=seed)
    precedents = precedent_frame()
    precedents["status"] = ATLAS4_TAG

    precedents.to_csv(tables / "atlas4_precedent_registry.csv", index=False)
    contrast.to_csv(tables / "atlas4_operator_contrast.csv", index=False)
    grid.to_csv(tables / "atlas4_exhaustive_grid.csv", index=False)
    modalities.to_csv(tables / "atlas4_measurement_falsifier.csv", index=False)
    translation.to_csv(tables / "atlas4_translation_sensitivity.csv", index=False)
    ceilings.to_csv(tables / "atlas4_criterion_ceilings.csv", index=False)
    interaction.to_csv(tables / "atlas4_correlation_by_operator.csv", index=False)
    pd.DataFrame(
        [
            {
                "construct": "single_dose",
                "probability_of_target": float(
                    contrast.loc[contrast["doses"].eq(1)].iloc[0]["probability_of_target"]
                ),
                "status": ATLAS4_TAG,
            },
            {
                "construct": "atlas4_human_tempered",
                "probability_of_target": corrected.probability_of_target,
                "status": ATLAS4_TAG,
            },
            {
                "construct": "atlas4_animal_optimistic",
                "probability_of_target": animal.probability_of_target,
                "status": ATLAS4_TAG,
            },
        ]
    ).to_csv(tables / "atlas4_headline.csv", index=False)
    pd.DataFrame(
        [
            {
                **{k: v for k, v in summary.items() if k != "best_admissible"},
                "best_admissible": str(summary["best_admissible"]),
                "status": ATLAS4_TAG,
            }
        ]
    ).to_csv(tables / "atlas4_search_summary.csv", index=False)
    (reports / "atlas4_report.txt").write_text(
        atlas4_report(
            corrected,
            animal,
            contrast,
            grid,
            summary,
            modalities,
            translation,
            ceilings,
            interaction,
        )
    )
    (root / "DISCLAIMER.txt").write_text(DISCLAIMER + "\n")
    return {
        "output_dir": str(root),
        "outputs": 11,
        "status": ATLAS4_TAG,
        "precedented_dose_limit": max_precedented_doses(),
        "probability_of_target": corrected.probability_of_target,
        "robust_with_margin": bool(corrected.robust_with_margin),
        "animal_optimistic_probability": animal.probability_of_target,
        "configurations_tested": summary["configurations_tested"],
        "saturation_artefacts": summary["saturation_artefacts"],
        "admissible": summary["admissible"],
        "cures_declared": 0,
        "disclaimer": DISCLAIMER,
    }
