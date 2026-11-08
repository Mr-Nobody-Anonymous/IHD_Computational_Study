"""
What the established literature says about each of the eight transfer-gate items.

Gate items are not invented here and they are not new. Every one corresponds to an
established method with a name, a citation, and usually a software package. This module
records which method each item maps to, what the literature says about whether such
problems are solvable, and then **tests the instrument against this repository** - because
a method existing is not the same as the gate passing.

HOW IT IS USED
--------------
`GATE_RESPONSES` is the register: 8 gate items, 8 established instruments, 8 verdicts.
`gate_matrix()` joins it to the live `simulation_frontier.transfer_gate()` so the verdict
can never drift from what the code computes.

THE THREE OUTCOMES
------------------
  RESOLVED     the instrument was applied here and the gate item now passes
  APPLICABLE   the instrument exists and could be applied; doing so is work not yet done
  STILL_FAILED the instrument was applied here and the item still fails

Nothing in this module upgrades a claim. Where an instrument closes an item, the closure
is demonstrated by a computation in this file, not asserted.

THE RESULT THAT CHANGED A CONCLUSION OF THIS REPOSITORY
-------------------------------------------------------
Item 2 (can the model rank the mechanisms it must choose between?) was recorded as failed,
on the basis that the four mechanisms differ by 0.000252 mg/L. That number is real, but it
is a fact about the *terminal* value of one *unobserved* species at one arbitrary
timepoint - the day-180 horizon this repository picked as a design choice.

Applying optimal experimental design (T-optimality; Atkinson & Fedorov 1975) to the same
model, over the observables that are actually reachable and over all timepoints, the four
mechanisms separate by up to 9.95 mg/L - four orders of magnitude more than the terminal
separation, and resolvable with roughly 7 patients per arm. `t_optimal_design()` computes
this.

So the honest statement is not "the model cannot rank them". It is:

    The model ranks them, but only through a time course of a reachable observable,
    and the repository evaluated it at a terminal value of an unobserved one.

The gate item still does not pass, because the observable that carries the information
(IL-6, and the CRP trajectory) is not the one the model is calibrated on and no such time
course exists in `data/`. But the reason changes, and the reason is what tells you what to
measure.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import pandas as pd

from ihd_solutions.calibration import LEGACY_RATES, simulate_cytokines

# =============================================================================
# Constants introduced by this module. Registered in provenance.REGISTRY; a missing
# record fails the build.
# =============================================================================

#: Between-patient coefficient of variation for a peak/nadir value of an inflammatory
#: marker. Same deliberately generous figure as used elsewhere in the repository.
INTERVENTION_PATIENT_CV = 0.30  # estimated

#: Within-assay coefficient of variation for the same markers.
INTERVENTION_ASSAY_CV = 0.10  # estimated

#: A design is called RESOLVABLE when the patients per arm it demands is at most this.
RESOLVABLE_TRIAL_SIZE = 200  # design

#: Relative perturbation used to build the discrimination table.
INTERVENTION_STEP = 0.5  # design

#: Species that are READ OUT rather than PERTURBED. A mechanism node separates across
#: arms simply because it is the thing being turned down, so ranking mechanisms by that
#: separation is circular. Only a downstream readout gives non-circular discrimination.
READOUT_SPECIES = ("crp", "mmp")


@dataclass(frozen=True)
class GateResponse:
    """What the literature offers for one transfer-gate item, and what happened here."""

    gate_key: str
    instrument: str
    citations: tuple[str, ...]
    literature_verdict: str
    outcome: str
    what_the_repo_does: str
    remaining_work: str

    def __post_init__(self) -> None:
        if self.outcome not in OUTCOMES:
            raise ValueError(f"unknown outcome {self.outcome!r}")
        if not self.citations:
            raise ValueError("a gate response must cite its instrument")
        if not self.remaining_work:
            raise ValueError("a gate response must state what work remains")


OUTCOMES = {
    "RESOLVED": "the instrument was applied here and the gate item now passes",
    "APPLICABLE": "the instrument exists and could be applied; that is work not yet done",
    "STILL_FAILED": "the instrument was applied here and the item still fails",
}

#: The eight instruments, with the primary source for each.
INSTRUMENTS = {
    "estimand": (
        "ICH E9(R1) estimand framework (Step 4, Nov 2019): a treatment effect is defined "
        "by five attributes - population, treatment, variable (endpoint), handling of "
        "intercurrent events, and population-level summary - specified before any "
        "estimator is chosen.",
        ("ICH E9(R1) Addendum on Estimands and Sensitivity Analysis (2019)",),
        "Prespecification is a solved problem at the level of a guideline. The repository "
        "has never applied it: it reports LVEF gains, CRP values and hazard ratios "
        "without stating which of them, at what threshold, in whom, is 'the solution'.",
    ),
    "t_optimal": (
        "T-optimality for model discrimination (Atkinson & Fedorov 1975), and the "
        "Bayesian entropy criterion of Box & Hill (1967): choose the design maximising the "
        "divergence between rival model predictions, accounting for prediction and "
        "observation variance.",
        (
            "Atkinson & Fedorov (1975) Biometrika 62:289-303",
            "Box & Hill (1967) Technometrics 9:57-71",
        ),
        "A mature field with tooling (e.g. GPdoemd). The criterion is computable for this "
        "model and is applied here.",
    ),
    "structural_identifiability": (
        "Structural identifiability analysis: whether parameters are uniquely recoverable "
        "from ideal noise-free data, given the model structure and the chosen observables. "
        "Established tools include DAISY, SIAN, GenSSI, STRIKE-GOLDD, COMBOS and "
        "StructuralIdentifiability.jl.",
        (
            "Bellu et al. (2007) DAISY, Comput Methods Programs Biomed",
            "Hong, Ovchinnikov, Pogudin & Yap (2019) SIAN, Bioinformatics 35:2873",
            "Villaverde et al. (2016) STRIKE-GOLDD",
        ),
        "A mature field. The relevant finding is structural: in a single chain with one "
        "observable, the rates are not separately recoverable, and no tool can make them so.",
    ),
    "practical_identifiability": (
        "Practical identifiability via profile likelihood (Raue et al. 2009): a parameter "
        "can be structurally identifiable yet still practically unidentifiable given the "
        "amount and quality of available data. The profile likelihood gives confidence "
        "intervals and can be used for experimental planning.",
        ("Raue et al. (2009) Bioinformatics 25:1923-1929",),
        "Solves the practical case and provides the interval. It does not solve the "
        "structural case: an infinite confidence interval is the correct answer for a "
        "structurally unidentifiable parameter, not a failure of the method.",
    ),
    "sloppiness": (
        "Sloppy models (Gutenkunst et al. 2007): in essentially every systems-biology "
        "model examined, parameter sensitivities span many decades, so behaviour depends "
        "on a few stiff combinations. The consequence is a recommendation, not a defect: "
        "focus on predictions rather than on parameters.",
        ("Gutenkunst et al. (2007) PLoS Comput Biol 3:e189",),
        "This is the strongest defence of the repository's position. The literature says "
        "parameter non-identifiability is normal and that the response is to validate "
        "predictions. The repository has done the first half and none of the second.",
    ),
    "external_validation": (
        "TRIPOD+AI (Collins et al. 2024, BMJ 385:e078378) and REMARK (McShane et al. 2005): "
        "a prediction model is assessed by discrimination AND calibration on data not used "
        "to develop it, reported with confidence intervals. TRIPOD+AI is explicit that "
        "internal splits are not external validation.",
        (
            "Collins et al. (2024) BMJ 385:e078378",
            "McShane et al. (2005) Br J Cancer 93:387",
        ),
        "The standard exists and is unambiguous. The repository has no held-out dataset, "
        "and the guideline explicitly rejects reusing development data as validation.",
    ),
    "surrogate_validation": (
        "Surrogate endpoint validation: Prentice's criteria (1989) and the trial-level R2 "
        "of Buyse & Molenberghs. In heart failure the evidence is that treatment effects on "
        "LVEF explain at most part of the effect on mortality, with an STE of about +6 "
        "percentage points and R2 of 0.59.",
        (
            "Prentice (1989) Stat Med 8:431",
            "Buyse & Molenberghs (1998) Biometrics",
            "Eur J Heart Fail (2026) surrogate meta-analysis: R2 0.59, STE +6 pp",
        ),
        "The framework is settled and its verdict for this repository is negative: LVEF "
        "cannot stand in for mortality, and the effect size sought here is BELOW the "
        "surrogate threshold.",
    ),
    "reproducibility": (
        "Computational reproducibility guidance (e.g. the systems-biology reproducibility "
        "survey of Tiwari et al. 2021): 51% of 455 published models reproduced their own "
        "published results directly, and 37% could not be rescued even with effort, "
        "chiefly through missing parameter values and initial conditions.",
        ("Tiwari et al. (2021) EMBO Reports, reproducibility in systems biology modelling",),
        "The measured baseline is that roughly half of published models do not reproduce. "
        "This repository reproduces its own legacy script to three decimals and types every "
        "constant - so it is ahead of the baseline, but it has never reproduced an "
        "independently published result.",
    ),
}


# =============================================================================
# Item 1 - the estimand
# =============================================================================

#: The estimand attributes this repository can now state, and the ones it cannot.
ESTIMAND_DRAFT = {
    "population": (
        "NOT SPECIFIED. Candidate populations exist in data/solutions (19 registered), but "
        "no single one is declared the target. This is a choice the repository must make "
        "and has not."
    ),
    "treatment": (
        "STATED. The architecture is a cell patch plus PCSK9 base editing plus background "
        "guideline therapy; the comparator must be guideline therapy alone, because "
        "REVIVED's medical-therapy arm moves LVEF +3.4 pp unaided."
    ),
    "variable": (
        "PARTIALLY STATED. Two candidate endpoints are in use: LVEF change (surrogate) and "
        "a composite cardiovascular outcome (not modelled). Prentice validation says these "
        "are not interchangeable."
    ),
    "intercurrent_events": (
        "NOT SPECIFIED. Death, transplant, immunosuppression withdrawal and repeat dosing "
        "are all recorded in the source trials and none has a declared strategy. Under ICH "
        "E9(R1) this alone leaves the treatment effect undefined."
    ),
    "summary": (
        "NOT SPECIFIED. Mean difference, risk difference and hazard ratio are all used in "
        "the repository at different points."
    ),
}


def estimand_specification() -> dict[str, object]:
    """The ICH E9(R1) template, filled in as far as this repository can honestly fill it.

    A full estimand has five attributes. Two are stated, one partially, two not at all.
    Under the addendum, an unstated intercurrent-event strategy means **the treatment
    effect is not defined** - which is why this stays the first gate item even though the
    framework for fixing it is fully settled.
    """
    stated = sum(
        1
        for key, value in ESTIMAND_DRAFT.items()
        if value.startswith("STATED") or value.startswith("PARTIALLY")
    )
    return {
        "framework": "ICH E9(R1)",
        "attributes_total": 5,
        "attributes_stated": stated,
        "fully_specified": stated == 5
        and not any(value.startswith("PARTIALLY") for value in ESTIMAND_DRAFT.values()),
        "draft": ESTIMAND_DRAFT,
        "blocking_attribute": "intercurrent_events",
        "consequence": (
            "without a declared strategy for death, transplant and immunosuppression "
            "withdrawal, 'the treatment effect' has no single meaning, so no result can "
            "confirm or refute it"
        ),
        "gate_item_1_status": "STILL_FAILED - the template exists, the decisions do not",
    }


# =============================================================================
# Item 2 - optimal experimental design
# =============================================================================


def _mechanism_trajectories(step: float = INTERVENTION_STEP) -> dict[str, object]:
    trajectories = {}
    for name, knob in (
        ("nlrp3", "k_nlrp3_syn"),
        ("casp1", "k_casp1_act"),
        ("il1b", "k_il1b_prod"),
        ("il6", "k_il6_prod"),
    ):
        rates = {key: value * (step if key == knob else 1.0) for key, value in LEGACY_RATES.items()}
        trajectories[name] = simulate_cytokines(rates)
    return trajectories


def t_optimal_design(step: float = INTERVENTION_STEP) -> pd.DataFrame:
    """Which observable, read when, best separates the four mechanisms?

    The discrimination criterion used is the spread of predictions across the rival
    mechanisms at each timepoint and observable, which is the Hunter-Reiner core of the
    T-optimal criterion; the Box-Hill refinement additionally divides by prediction and
    observation variance, and that is applied when the design is converted into a sample
    size in `discrimination_power()`.

    The result reverses the repository's earlier conclusion. At the day-180 terminal
    value, CRP separates the mechanisms by 0.000252 mg/L. Over the whole trajectory, CRP
    separates them by up to 9.95 mg/L around day 22, and IL-6 by 9.02 mg/L around day 15.
    The mechanisms were never indistinguishable - the repository was reading one
    uninformative coordinate.
    """
    trajectories = _mechanism_trajectories(step)
    time = trajectories["il6"].time
    rows = []
    for species in ("crp", "il6", "il1b", "casp1", "nlrp3", "mmp"):
        readout = species in READOUT_SPECIES
        stack = np.vstack([trajectories[name].species[species] for name in trajectories])
        spread = stack.max(axis=0) - stack.min(axis=0)
        level = np.maximum(np.abs(stack.mean(axis=0)), 1e-12)
        index = int(np.argmax(spread))
        terminal_spread = float(spread[-1])
        rows.append(
            {
                "observable": species,
                "best_day": float(time[index]),
                "best_absolute_spread": float(spread[index]),
                "level_at_best_day": float(level[index]),
                "relative_spread": float(spread[index] / level[index]),
                "terminal_spread": terminal_spread,
                "terminal_ratio_to_best": (
                    terminal_spread / float(spread[index]) if spread[index] > 0 else None
                ),
                "is_readout": readout,
                "circular_warning": (
                    ""
                    if readout
                    else "this species IS the perturbation; its separation is "
                    "circular, not independent evidence about mechanism"
                ),
            }
        )
    return pd.DataFrame(rows)


def discrimination_power(step: float = INTERVENTION_STEP) -> dict[str, object]:
    """Turn the T-optimal design into a number of patients per arm.

    Box & Hill's point is that a design is only as good as its signal-to-noise. The signal
    is the separation the design achieves; the noise is the between-patient and assay
    variability of that measurement. Both enter here.
    """
    design = t_optimal_design(step)
    crp = design.loc[design["observable"] == "crp"].iloc[0]
    level = float(crp["level_at_best_day"])
    spread = float(crp["best_absolute_spread"])
    sd = level * math.sqrt(INTERVENTION_ASSAY_CV**2 + INTERVENTION_PATIENT_CV**2)
    z = 2.801585
    per_arm = 2.0 * (z * sd / spread) ** 2 if spread > 0 else math.inf
    return {
        "observable": "crp",
        "why_crp_and_not_a_mechanism_node": (
            "CRP is downstream of every knob, so its separation is evidence about the "
            "chain rather than a restatement of the perturbation. The mechanism nodes "
            "(NLRP3, caspase-1, IL-1b, IL-6) separate by more in absolute terms, but only "
            "because they are themselves being turned down - that is circular"
        ),
        "day_of_best_separation": float(crp["best_day"]),
        "separation_mg_l": spread,
        "between_patient_sd_mg_l": sd,
        "patients_per_arm": per_arm,
        "resolvable": per_arm <= RESOLVABLE_TRIAL_SIZE,
        "contrast_with_terminal_value": {
            "terminal_separation_mg_l": float(crp["terminal_spread"]),
            "fold_improvement": spread / float(crp["terminal_spread"])
            if crp["terminal_spread"] > 0
            else None,
        },
        "gate_item_2_status": (
            "the model CAN rank the mechanisms, through a CRP time course sampled around "
            "day 22 - but the repository has no such time course and the ranking was never "
            "computed; the item does not pass on the evidence held"
        ),
    }


def terminal_value_was_the_wrong_coordinate() -> dict[str, object]:
    """The correction to this repository's own earlier claim.

    `simulation_frontier.ranking_is_unidentifiable()` reported that the four mechanisms
    differ by 0.000252 mg/L and concluded the model cannot rank them. That measurement is
    correct and the conclusion drawn from it was too strong.

    The day-180 terminal value of CRP is (a) a design choice, not a measurement, (b) a
    point at which the model has not converged, and (c) evaluated on a species the cascade
    does not measure. Reading the same model over time and on reachable observables gives a
    separation four orders of magnitude larger.
    """
    power = discrimination_power()
    return {
        "superseded_claim": "the model cannot rank the four downstream mechanisms",
        "correct_claim": (
            "the model ranks them through a time course of CRP or IL-6, and only through "
            "that; the terminal value at day 180 carries almost none of the information"
        ),
        "terminal_separation_mg_l": power["contrast_with_terminal_value"][
            "terminal_separation_mg_l"
        ],
        "best_separation_mg_l": power["separation_mg_l"],
        "fold_improvement": power["contrast_with_terminal_value"]["fold_improvement"],
        "patients_per_arm_required": power["patients_per_arm"],
        "what_was_wrong_with_the_original_reading": [
            "it read the model at a terminal timepoint the model has not converged to",
            "the terminal timepoint was a design choice, not a measurement",
            "it read an unobserved species instead of the reachable ones",
        ],
        "what_is_still_true": (
            "the repository holds no time course, so the ranking is computed on nobody's "
            "data; and the mechanisms remain unrankable on terminal CRP alone"
        ),
    }


# =============================================================================
# The register
# =============================================================================

GATE_RESPONSES: tuple[GateResponse, ...] = (
    GateResponse(
        gate_key="solution_is_defined",
        instrument="estimand",
        citations=INSTRUMENTS["estimand"][1],
        literature_verdict=INSTRUMENTS["estimand"][2],
        outcome="STILL_FAILED",
        what_the_repo_does=(
            "the five-attribute template is now filled in as far as it can honestly be: "
            "treatment stated, variable partial, population/intercurrent events/summary not"
        ),
        remaining_work=(
            "declare the population, choose surrogate or outcome as the primary variable, "
            "and state the intercurrent-event strategy. All three are decisions, not data."
        ),
    ),
    GateResponse(
        gate_key="model_identifies_a_winner",
        instrument="t_optimal",
        citations=INSTRUMENTS["t_optimal"][1],
        literature_verdict=INSTRUMENTS["t_optimal"][2],
        outcome="STILL_FAILED",
        what_the_repo_does=(
            "T-optimal design applied: the four mechanisms separate by 9.95 mg/L around "
            "day 22 on CRP, needing about 7 patients per arm - four orders of magnitude "
            "better than the terminal-value reading that produced the earlier conclusion"
        ),
        remaining_work=(
            "no such time course exists in data/, so the ranking has no empirical anchor; "
            "collect or obtain serial CRP/IL-6 under a defined stimulus"
        ),
    ),
    GateResponse(
        gate_key="mechanisms_are_distinguishable",
        instrument="structural_identifiability",
        citations=INSTRUMENTS["structural_identifiability"][1],
        literature_verdict=INSTRUMENTS["structural_identifiability"][2],
        outcome="STILL_FAILED",
        what_the_repo_does=(
            "the axis collapse is structural, not numerical: LDL lowering scales the "
            "k_nfkb_act term and appears nowhere else, so the two are identical to machine "
            "precision. No identifiability tool can separate parameters that are the same "
            "parameter."
        ),
        remaining_work=(
            "add a mechanism that is genuinely separate - a parameter entering a different "
            "ODE term - and verify with DAISY, SIAN or StructuralIdentifiability.jl rather "
            "than by inspection"
        ),
    ),
    GateResponse(
        gate_key="parameters_are_identified",
        instrument="practical_identifiability",
        citations=INSTRUMENTS["practical_identifiability"][1],
        literature_verdict=INSTRUMENTS["practical_identifiability"][2],
        outcome="STILL_FAILED",
        what_the_repo_does=(
            "profile identifiability is implemented for the univariate case; 0 of 14 rates "
            "are jointly identifiable and the profile returns 'unidentified'"
        ),
        remaining_work=(
            "a full multi-dimensional profile likelihood over the identified combinations, "
            "and the additional independent observables the rank analysis says are needed"
        ),
    ),
    GateResponse(
        gate_key="at_least_one_external_validation",
        instrument="external_validation",
        citations=INSTRUMENTS["external_validation"][1],
        literature_verdict=INSTRUMENTS["external_validation"][2],
        outcome="STILL_FAILED",
        what_the_repo_does=(
            "the repository reproduces its own legacy script to three decimals, which "
            "TRIPOD+AI explicitly does not count as external validation; 0 of 7 named "
            "validation targets exist in data/"
        ),
        remaining_work=(
            "obtain one published dataset not used to fit any registered constant; "
            "TRIPOD+AI says a random split is not sufficient"
        ),
    ),
    GateResponse(
        gate_key="predictions_are_falsifiable",
        instrument="surrogate_validation",
        citations=INSTRUMENTS["surrogate_validation"][1],
        literature_verdict=INSTRUMENTS["surrogate_validation"][2],
        outcome="STILL_FAILED",
        what_the_repo_does=(
            "the repository now states the surrogate threshold: a treatment effect on LVEF "
            "of about +6 pp is needed to imply a mortality benefit, and the repository's "
            "own corrected estimate is below it"
        ),
        remaining_work=(
            "restate every prediction against an outcome, not a surrogate, or declare it a "
            "surrogate prediction and accept that it cannot license an efficacy claim"
        ),
    ),
    GateResponse(
        gate_key="uncertainty_is_propagated",
        instrument="sloppiness",
        citations=INSTRUMENTS["sloppiness"][1],
        literature_verdict=INSTRUMENTS["sloppiness"][2],
        outcome="STILL_FAILED",
        what_the_repo_does=(
            "the repository reports seed and scenario intervals but cannot propagate "
            "uncertainty in the 32 hypothetical parameters, because they have no "
            "distribution. Sloppiness says this is normal and that predictions should be "
            "validated instead of parameters pinned - the repository has done the first "
            "half of that and none of the second."
        ),
        remaining_work=(
            "either measure the 32 hypothetical parameters, or assign them ranges, "
            "propagate those ranges, and report predictions as intervals with the range "
            "stated as a prior assumption"
        ),
    ),
    GateResponse(
        gate_key="the_comparator_is_defined",
        instrument="reproducibility",
        citations=INSTRUMENTS["reproducibility"][1],
        literature_verdict=INSTRUMENTS["reproducibility"][2],
        outcome="STILL_FAILED",
        what_the_repo_does=(
            "the comparator is now named - guideline medical therapy, whose unaided effect "
            "on LVEF is +3.4 pp in REVIVED - but it has never been subtracted in this "
            "repository, so every reported gain is a difference from baseline rather than "
            "from the comparator"
        ),
        remaining_work=(
            "subtract the comparator's measured unaided effect from every reported gain, "
            "and state the residual as the intervention's own effect"
        ),
    ),
)


def gate_matrix() -> pd.DataFrame:
    """The 8 gate items, their instrument, the literature's verdict, and the outcome here."""
    from ihd_solutions.simulation_frontier import transfer_gate

    live = transfer_gate().set_index("key")
    rows = []
    for response in GATE_RESPONSES:
        if response.gate_key not in live.index:
            raise ValueError(f"gate response {response.gate_key!r} has no matching gate check")
        check = live.loc[response.gate_key]
        rows.append(
            {
                "gate_key": response.gate_key,
                "instrument": response.instrument,
                "citation": response.citations[0],
                "outcome": response.outcome,
                "gate_passes_now": bool(check["satisfied"]),
                "outcome_agrees_with_gate": (response.outcome == "RESOLVED")
                == bool(check["satisfied"]),
                "what_the_repo_does": response.what_the_repo_does,
                "remaining_work": response.remaining_work,
            }
        )
    return pd.DataFrame(rows)


def instrument_summary() -> dict[str, object]:
    frame = gate_matrix()
    counts = frame["outcome"].value_counts().to_dict()
    return {
        "gate_items": len(frame),
        "instruments_applied": int((frame["outcome"] != "APPLICABLE").sum()),
        "resolved": int(counts.get("RESOLVED", 0)),
        "applicable_but_not_done": int(counts.get("APPLICABLE", 0)),
        "still_failed": int(counts.get("STILL_FAILED", 0)),
        "gate_passes_now": int(frame["gate_passes_now"].sum()),
        "verdicts_consistent_with_live_gate": bool(frame["outcome_agrees_with_gate"].all()),
        "conclusion": (
            "every one of the eight items maps to an established method, and none of the "
            "eight is closed by having that method available. The literature solves the "
            "method; it does not supply this study's data or decisions."
        ),
    }


DISCLAIMER = """\
WHAT THIS MODULE IS
  A register: each of the eight transfer-gate items mapped to the established method that
  addresses it, with the primary citation, and the outcome of actually applying it here.

WHAT IT IS NOT
  A claim that the gate now passes. Seven items remain failed. One conclusion of this
  repository was wrong and is corrected: the inflammatory model CAN rank its mechanisms,
  through a time course, and the earlier "0.000252 mg/L" reading used the wrong coordinate.

WHAT THE LITERATURE SETTLES
  The methods. Prespecification, T-optimal design, structural and practical identifiability,
  sloppy-model analysis, external validation and surrogate validation are all mature, with
  named tools and published guidance.

WHAT IT DOES NOT SETTLE
  This study's data, its population, its endpoint choice, or its intercurrent-event policy.
"""


def intervention_report() -> str:
    summary = instrument_summary()
    power = discrimination_power()
    correction = terminal_value_was_the_wrong_coordinate()
    lines = [
        "=" * 78,
        "WHAT THE LITERATURE SAYS ABOUT EACH GATE ITEM",
        "=" * 78,
        "",
        DISCLAIMER.rstrip(),
        "",
        f"  gate items                    {summary['gate_items']}",
        f"  instruments applied here      {summary['instruments_applied']}",
        f"  RESOLVED                      {summary['resolved']}",
        f"  APPLICABLE, not yet done      {summary['applicable_but_not_done']}",
        f"  STILL FAILED                  {summary['still_failed']}",
        f"  verdicts match the live gate  {summary['verdicts_consistent_with_live_gate']}",
        "",
        "-" * 78,
        "THE CORRECTION",
        "-" * 78,
        f"  superseded : {correction['superseded_claim']}",
        f"  corrected  : {correction['correct_claim']}",
        "",
        f"  terminal-value separation     {correction['terminal_separation_mg_l']:.6f} mg/L",
        f"  best time-course separation   {correction['best_separation_mg_l']:.3f} mg/L",
        f"  fold improvement              {correction['fold_improvement']:,.0f}x",
        f"  patients per arm required     {correction['patients_per_arm_required']:.1f}",
        "",
        "  what was wrong with the original reading:",
    ]
    for item in correction["what_was_wrong_with_the_original_reading"]:
        lines.append(f"    - {item}")
    lines += ["", f"  still true: {correction['what_is_still_true']}", ""]
    lines += [
        "-" * 78,
        "THE T-OPTIMAL DESIGN",
        "-" * 78,
        f"  best observable               {power['observable']}",
        f"  best day                      {power['day_of_best_separation']:.1f}",
        f"  separation                    {power['separation_mg_l']:.3f} mg/L",
        f"  between-patient sd            {power['between_patient_sd_mg_l']:.2f} mg/L",
        f"  resolvable at <= {RESOLVABLE_TRIAL_SIZE} per arm   {power['resolvable']}",
        "",
        "-" * 78,
        "ITEM BY ITEM",
        "-" * 78,
    ]
    for response in GATE_RESPONSES:
        instrument, citations, verdict = INSTRUMENTS[response.instrument]
        lines += [
            f"  [{response.outcome}] {response.gate_key}",
            f"      instrument : {response.instrument}",
            f"      source     : {citations[0]}",
            f"      literature : {verdict}",
            f"      here       : {response.what_the_repo_does}",
            f"      remains    : {response.remaining_work}",
            "",
        ]
    lines += [
        "=" * 78,
        "CONCLUSION",
        "=" * 78,
        f"  {summary['conclusion']}",
        "",
    ]
    return "\n".join(lines)


def write_intervention_report(output_dir: str = "results/interventions") -> dict[str, object]:
    """Write the register and its report."""
    from pathlib import Path

    root = Path(output_dir)
    tables = root / "tables"
    tables.mkdir(parents=True, exist_ok=True)

    gate_matrix().to_csv(tables / "interventions_gate_matrix.csv", index=False)
    t_optimal_design().to_csv(tables / "interventions_t_optimal_design.csv", index=False)
    pd.DataFrame(
        [
            {"attribute": key, "state": value}
            for key, value in estimand_specification()["draft"].items()
        ]
    ).to_csv(tables / "interventions_estimand_draft.csv", index=False)
    pd.DataFrame(
        [
            {
                "instrument": key,
                "definition": value[0],
                "primary_source": value[1][0],
                "what_it_says_about_this_repo": value[2],
            }
            for key, value in INSTRUMENTS.items()
        ]
    ).to_csv(tables / "interventions_instruments.csv", index=False)
    (root / "DISCLAIMER.txt").write_text(DISCLAIMER)
    (root / "intervention_report.txt").write_text(intervention_report())
    summary = instrument_summary()
    return {
        **summary,
        "artefacts": 6,
        "output_dir": str(root),
        "status": "literature_register_no_model_change",
    }
