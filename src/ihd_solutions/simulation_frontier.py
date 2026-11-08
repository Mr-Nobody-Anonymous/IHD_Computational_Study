"""
The simulation-first frontier: what simulation can deliver, and where it stops.

THE QUESTION THIS ANSWERS
-------------------------
"Find the solution in simulation only, then transfer to the lab." What would have to
exist for that plan to work, and which of those things exist here?

The answer is a **transfer frontier**: the set of requirements that must all hold before
a simulation result is worth a laboratory's time. Each requirement is checked in code.
Most fail. The failures are not excuses - they name exactly what to build, in order, and
they are cheaper than the experiments they would save.

THE CENTRAL RESULT
------------------
`ranking_is_unidentifiable()`. The legacy inflammatory model contains six mechanisms that
a laboratory could target. The model can order them into coarse classes, but **it cannot
rank the four downstream mechanisms at all** - on terminal CRP they differ by 0.000252
mg/L, and separating them there would need billions of patients per arm.

This is not a resolution problem that more compute solves: each of these parameters is a
multiplicative factor on its own ODE term, so the terminal value is set by the product of
the chain and is invariant to which factor is adjusted.

**SUPERSEDED IN PART** by `interventions.terminal_value_was_the_wrong_coordinate()`. The
terminal reading is the DAY-180 value of an UNOBSERVED species at a DESIGN-CHOSEN horizon.
Over the trajectory of a reachable observable the same four mechanisms separate by 9.95
mg/L around day 22 - about 39,500 times more, resolvable with roughly 7 patients per arm.
The model CAN rank them through a time course; it cannot rank them on the terminal reading,
and the repository holds no time course. The arithmetic below is retained unchanged so the
correction is visible rather than silently edited away.

WHAT THAT DOES NOT MEAN
-----------------------
It does not mean simulation is useless. It means the deliverable changes shape. A
simulation that cannot name a winner CAN do something a laboratory cannot: state precisely
which measurement would discriminate, and at what resolution. The output is an
**experiment specification**, not a cure. `experiment_specification()` produces it.

No model is changed here, nothing is tuned, and no claim is upgraded.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path as _Path

import numpy as np
import pandas as pd

from ihd_solutions.calibration import LEGACY_RATES, simulate_cytokines

# =============================================================================
# Constants introduced by this module. Every one is registered in
# `provenance.REGISTRY`; a missing record fails the build.
# =============================================================================

#: Within-assay coefficient of variation for high-sensitivity CRP. Clinical assays quote
#: 5-10%; the upper end is used so the resolution requirement is not flattered.
CRP_ASSAY_CV = 0.10  # estimated

#: Between-patient coefficient of variation of CRP. Inflammatory markers are notoriously
#: variable between individuals; 30% is a deliberately generous figure, chosen because a
#: smaller one would make the resolution requirement look easier than it is.
CRP_PATIENT_CV = 0.30  # estimated

#: A discrimination is called STRUCTURAL when mechanisms are separated by less than this
#: fraction of the observable's own value - i.e. below any plausible measurement.
STRUCTURAL_BLINDNESS = 1e-3  # design

#: Two-arm 80% power, two-sided 5%, as a z-sum.
FRONTIER_POWER_Z = 2.801585  # estimated

#: Relative perturbation used when probing mechanism additivity.
FRONTIER_STEP = 0.5  # design


# =============================================================================
# 1. CAN THE MODEL NAME A WINNER?
# =============================================================================


def mechanism_effect_table(step: float = FRONTIER_STEP) -> pd.DataFrame:
    """Terminal and peak CRP under each single mechanism, perturbed by `step`.

    Six mechanisms a laboratory could plausibly target are expressed as the model's own
    parameters. Crucially, **LDL lowering and NF-kB inhibition are the same perturbation
    in this model** (`oxldl_factor` multiplies the NF-kB activation rate and nothing
    else), so the table reports them as one row with a note - reporting two would imply
    the model distinguishes two mechanisms it cannot tell apart.
    """
    if not 0.0 < step < 1.0:
        raise ValueError("step must be a fraction strictly inside (0, 1)")

    def run(**kwargs) -> tuple[float, float]:
        trajectory = simulate_cytokines(
            {
                name: value * (step if name in kwargs else 1.0)
                for name, value in LEGACY_RATES.items()
            },
            oxldl_factor=kwargs.pop("oxldl_factor", 1.0),
        )
        return trajectory.final("crp"), float(trajectory.species["crp"].max())

    names = {
        "untreated": {},
        "ldl_lowering_or_nfkb_inhibition": {"oxldl_factor": step},
        "nlrp3_inflammasome": {"k_nlrp3_syn": step},
        "caspase1": {"k_casp1_act": step},
        "il1b": {"k_il1b_prod": step},
        "il6": {"k_il6_prod": step},
    }
    rows = []
    for label, kwargs in names.items():
        final, peak = run(**kwargs)
        rows.append(
            {
                "mechanism": label,
                "final_crp_mg_l": final,
                "peak_crp_mg_l": peak,
                "aliases_not_distinguished": (
                    "LDL lowering == NF-kB inhibition exactly; one row, not two"
                    if label == "ldl_lowering_or_nfkb_inhibition"
                    else ""
                ),
            }
        )
    frame = pd.DataFrame(rows)
    untreated = float(frame.loc[frame["mechanism"] == "untreated", "final_crp_mg_l"].iloc[0])
    frame["reduction_fraction"] = 1.0 - frame["final_crp_mg_l"] / untreated
    return frame


def ldl_and_nfkb_are_the_same_perturbation() -> dict[str, object]:
    """The model's LDL axis and its NF-kB axis are one axis wearing two names.

    `oxldl_factor` multiplies `k_nfkb_act` and appears nowhere else, so scaling one is
    algebraically identical to scaling the other. This is proven by evaluating both and
    taking the difference across every species and every time point.

    Anyone who reads this model as "we can compare an LDL drug against an anti-inflammatory
    drug" is comparing a parameter against itself.
    """
    factors = (0.25, 0.5, 0.75, 1.5, 2.0)
    worst = 0.0
    for factor in factors:
        via_ldl = simulate_cytokines(LEGACY_RATES, oxldl_factor=factor)
        via_rate = simulate_cytokines(
            {**LEGACY_RATES, "k_nfkb_act": LEGACY_RATES["k_nfkb_act"] * factor}
        )
        for species in via_ldl.species:
            worst = max(
                worst, float(np.abs(via_ldl.species[species] - via_rate.species[species]).max())
            )
    return {
        "status": "TWO_PARAMETERS_ONE_AXIS",
        "max_absolute_difference_over_all_species": worst,
        "identical_to_machine_precision": worst < 1e-12,
        "consequence": (
            "the model cannot compare an LDL-lowering strategy against an "
            "anti-inflammatory one; the comparison is a parameter against itself"
        ),
        "what_survives": ("upstream-versus-downstream ordering is still visible; identity is not"),
    }


def ranking_is_unidentifiable(step: float = FRONTIER_STEP) -> dict[str, object]:
    """The result that decides whether simulation can name a winner.

    Four mechanistically distinct anti-inflammatory targets - NLRP3, caspase-1, IL-1b and
    IL-6 - are compared in terminal CRP. They differ by **0.000252 mg/L at most**. The
    model is not merely imprecise about which one is better: it is structurally silent,
    because each of these parameters is a multiplicative factor on its own ODE term and
    the terminal value is set by the product of the chain.

    Converting that separation into a trial size, at a between-patient CV of 30% plus a
    within-assay CV of 10% on a 44 mg/L value, gives billions of patients per arm. The
    absolute figure should not be read too literally - it is an order-of-magnitude
    statement - but its conclusion is robust to any plausible change in the CVs.
    """
    table = mechanism_effect_table(step)
    downstream = table.loc[
        table["mechanism"].isin(["nlrp3_inflammasome", "caspase1", "il1b", "il6"])
    ]
    values = downstream["final_crp_mg_l"].to_numpy()
    spread = float(values.max() - values.min())
    reference = float(table["final_crp_mg_l"].mean())
    sd = reference * math.sqrt(CRP_ASSAY_CV**2 + CRP_PATIENT_CV**2)
    per_arm = 2.0 * (FRONTIER_POWER_Z * sd / spread) ** 2 if spread > 0.0 else math.inf
    return {
        "mechanisms_compared": list(downstream["mechanism"]),
        "max_terminal_crp_spread_mg_l": spread,
        "spread_as_fraction_of_value": spread / reference,
        "structurally_blind": spread / reference < STRUCTURAL_BLINDNESS,
        "between_patient_sd_mg_l": sd,
        "patients_per_arm_to_separate": per_arm,
        "upstream_vs_downstream_is_visible": True,
        "within_class_ranking_is_not": True,
        "why_more_compute_does_not_help_on_this_coordinate": (
            "on this coordinate the separation is not a numerical artefact that finer "
            "integration would resolve; it is a property of a model in which each "
            "mechanism enters as a multiplicative factor on a linear chain, so the chain "
            "product is invariant to which factor is adjusted. It IS an artefact of the "
            "coordinate: read over time, the same model separates them by four orders of "
            "magnitude more (interventions.t_optimal_design)"
        ),
        "why_this_blocks_the_plan": (
            "on this coordinate, choosing WHICH mechanism matters requires the very "
            "experiment the choice is meant to design. That circularity applies to the "
            "terminal reading only, and it is BROKEN by the T-optimal design in "
            "interventions, which names the measurement that would separate them and the "
            "resolution it needs - a better design, not more compute"
        ),
        "status": "ranking_unidentifiable_on_terminal_value_only",
    }


def what_simulation_can_still_do() -> pd.DataFrame:
    """The deliverables that survive the identifiability result.

    Being unable to name a winner is not being unable to contribute. Each row is something
    a simulation can settle that a laboratory cannot, with the caveat attached.
    """
    rows = [
        {
            "deliverable": "falsifiable prediction with a stated resolution",
            "what_it_settles": "which measurement would refute the mechanism, and how precise",
            "why_a_lab_cannot": "a lab must first choose what to measure, and that is the question",
            "caveat": "the prediction is about the model, not the patient",
        },
        {
            "deliverable": "minimum detectable difference for a proposed trial",
            "what_it_settles": "whether a planned study can answer its own question at all",
            "why_a_lab_cannot": "a lab learns this only after the study is over",
            "caveat": "requires an assumed between-patient SD, which is swept",
        },
        {
            "deliverable": "identification of redundant mechanisms",
            "what_it_settles": "which combinations cannot add benefit because they share a pathway",
            "why_a_lab_cannot": "a lab would have to run the combination to discover this",
            "caveat": "redundancy is a property of THIS model's structure",
        },
        {
            "deliverable": "arithmetic limits on what any trial of this size could show",
            "what_it_settles": "whether a published result is distinguishable from noise",
            "why_a_lab_cannot": "it is arithmetic, and nobody runs it",
            "caveat": "none - this one is pure arithmetic and is the strongest deliverable",
        },
        {
            "deliverable": "an ordered list of experiments by information per cost",
            "what_it_settles": "which measurement would collapse the most uncertainty",
            "why_a_lab_cannot": "a lab cannot try every experiment; simulation can rank them",
            "caveat": "the ranking is only as good as the assumed structure",
        },
    ]
    return pd.DataFrame(rows)


# =============================================================================
# 2. WHAT DATA WOULD MAKE EACH MODEL FALSIFIABLE
# =============================================================================


@dataclass(frozen=True)
class ValidationTarget:
    """One measurement that would turn a model from unfalsifiable into testable."""

    module: str
    quantity: str
    resolution_required: str
    dataset_type: str
    why_it_is_the_target: str
    currently_available: bool

    def __post_init__(self) -> None:
        if not self.why_it_is_the_target:
            raise ValueError("a validation target must say why it is the right one")
        if not self.resolution_required:
            raise ValueError("a validation target must state the resolution required")


def validation_targets() -> pd.DataFrame:
    """The data that would move each module up the ladder, and none of it is present.

    This is the concrete answer to "what is missing, computationally". It is not more
    model. It is **measurement**. Every row is a table someone has published or could
    publish, and no row exists in `data/`.
    """
    targets = [
        ValidationTarget(
            module="inflammation.py",
            quantity="time course of IL-6 and CRP after a defined inflammatory stimulus, "
            "with the dose recorded",
            resolution_required=">= 6 time points per subject; inter-individual SD reported",
            dataset_type="human endotoxin challenge or ACS cohort with serial sampling",
            why_it_is_the_target=(
                "the model predicts a SHAPE, and shape is what a single terminal value "
                "cannot test; without a time course the ODE is unfalsifiable"
            ),
            currently_available=False,
        ),
        ValidationTarget(
            module="inflammation.py",
            quantity="effect of a specific NLRP3 or IL-1b inhibitor on the same time course",
            resolution_required="randomised, >= 3 doses, placebo-controlled",
            dataset_type="a target-engagement study (e.g. canakinumab-class)",
            why_it_is_the_target=(
                "the model has ONE knobs-per-mechanism structure with no dose-response; "
                "a randomised dose-response is the only data that can reject it"
            ),
            currently_available=False,
        ),
        ValidationTarget(
            module="plaque.py",
            quantity="paired cap thickness and LDL on the same patients over time",
            resolution_required=">= 2 imaging time points; >= 200 patients",
            dataset_type="serial OCT or IVUS cohort with paired lipids",
            why_it_is_the_target=(
                "the plaque model predicts cap thinning as a function of LDL; no dataset "
                "here pairs the two longitudinally, so the relation is asserted"
            ),
            currently_available=False,
        ),
        ValidationTarget(
            module="revival.py",
            quantity="randomised LVEF and outcome data for a cell-therapy arm versus "
            "guideline therapy",
            resolution_required=">= 26 per arm for LVEF; a mortality endpoint for outcome",
            dataset_type="randomised controlled trial with an outcome endpoint",
            why_it_is_the_target=(
                "the model's central claim is about an uncontrolled observation; only a "
                "control arm can price it, and none exists"
            ),
            currently_available=False,
        ),
        ValidationTarget(
            module="recombine.py",
            quantity="measured efficacy of two of the combined mechanisms administered together",
            resolution_required="factorial or add-on design, >= 2 mechanisms",
            dataset_type="a combination trial",
            why_it_is_the_target=(
                "the interaction terms are the entire justification for refusing additivity; "
                "they are measured in THIS model and nowhere in humans"
            ),
            currently_available=False,
        ),
        ValidationTarget(
            module="precedent.py",
            quantity="adverse-event rate for the delivery procedure with a defined denominator",
            resolution_required="registry-grade denominator, adjudicated events",
            dataset_type="a prospective registry",
            why_it_is_the_target=(
                "28 patients across two trials give no denominator; a rate cannot be "
                "estimated from single-arm counts, only assumed"
            ),
            currently_available=False,
        ),
        ValidationTarget(
            module="calibration.py (legacy cascade)",
            quantity="independent measurements of >= 3 species in the chain simultaneously",
            resolution_required="3 independent observables, which would give 3 identifiable "
            "directions against the 1 currently available",
            dataset_type="multiplex cytokine panel with absolute quantification",
            why_it_is_the_target=(
                "the cascade has 14 rates and 1 usable observable; each additional "
                "independent species removes one degree of structural blindness"
            ),
            currently_available=False,
        ),
    ]
    return pd.DataFrame([target.__dict__ for target in targets])


def validation_arithmetic() -> dict[str, object]:
    """How many independent observables the cascade needs, from its own rank.

    A linear chain of 14 rates with 7 state variables does not have 14 free directions.
    The count of usable directions is bounded by the number of INDEPENDENT observables,
    and the repository has one (terminal CRP, plus its integral, which is not independent).
    """
    frame = validation_targets()
    return {
        "observables_currently_available": 1,
        "observables_used": 2,
        "rates_total": len(LEGACY_RATES),
        "rates_jointly_identifiable": 0,
        "validation_targets_listed": len(frame),
        "validation_targets_present": int(frame["currently_available"].sum()),
        "minimum_observables_to_identify_all_rates": len(LEGACY_RATES),
        "honest_note": (
            "14 independent observables would be needed for 14 rates; two are available. "
            "This is why the axis is one-dimensional and why the four downstream "
            "mechanisms are indistinguishable."
        ),
    }


# =============================================================================
# 3. THE TRANSFER GATE
# =============================================================================


@dataclass(frozen=True)
class GateCheck:
    """One precondition for transferring a simulation result to a laboratory."""

    key: str
    requirement: str
    satisfied: bool
    evidence: str
    if_failed: str

    def __post_init__(self) -> None:
        if not self.requirement:
            raise ValueError("a gate check must state its requirement")
        if not self.satisfied and not self.if_failed:
            raise ValueError("a failed gate check must say what happens if it fails")


def transfer_gate() -> pd.DataFrame:
    """Every precondition for "simulation done, now transfer to lab", checked.

    A laboratory transfer is expensive and exposes patients. The gate is the list of things
    that must be true before it is worth doing. **Seven of eight fail.** That is the
    answer to the question, expressed as a checklist rather than a discussion.
    """
    ranking = ranking_is_unidentifiable()
    identical = ldl_and_nfkb_are_the_same_perturbation()
    targets = validation_targets()
    checks = [
        GateCheck(
            key="solution_is_defined",
            requirement=(
                "a pre-specified endpoint, effect size and population defining what counts "
                "as 'the solution'"
            ),
            satisfied=False,
            evidence=(
                "the repository reports LVEF gains, CRP values and hazard ratios, but has "
                "never pre-specified which of these, at what threshold, in whom, would "
                "constitute success"
            ),
            if_failed=(
                "without a definition, any result can be narrated as success after the "
                "fact; this is the single most important missing item"
            ),
        ),
        GateCheck(
            key="model_identifies_a_winner",
            requirement="the model can rank the candidate mechanisms it is asked to choose between",
            satisfied=not ranking["within_class_ranking_is_not"],
            evidence=(
                f"four downstream mechanisms differ by {ranking['max_terminal_crp_spread_mg_l']:.6f} "
                f"mg/L; separating them needs {ranking['patients_per_arm_to_separate']:,.0f} "
                "patients per arm"
            ),
            if_failed=(
                "the simulation cannot select what to test; it can only exclude what is "
                "structurally impossible"
            ),
        ),
        GateCheck(
            key="mechanisms_are_distinguishable",
            requirement="the model distinguishes the interventions it claims to compare",
            satisfied=not identical["identical_to_machine_precision"],
            evidence=(
                "LDL lowering and NF-kB inhibition are the same perturbation to machine "
                f"precision (max difference {identical['max_absolute_difference_over_all_species']:.1e})"
            ),
            if_failed=(
                "any 'LDL drug versus anti-inflammatory drug' comparison in this model is "
                "a parameter compared with itself"
            ),
        ),
        GateCheck(
            key="parameters_are_identified",
            requirement="the fitted parameters are determined by the data rather than chosen",
            satisfied=False,
            evidence="0 of 14 cytokine rates are jointly identifiable; 32 of 92 registered "
            "parameters are hypothetical",
            if_failed=(
                "a laboratory would be testing the modeller's assumptions, not the "
                "hypothesis; a negative result would be uninterpretable"
            ),
        ),
        GateCheck(
            key="at_least_one_external_validation",
            requirement="the model predicts data it was not built from",
            satisfied=int(targets["currently_available"].sum()) > 0,
            evidence=f"{int(targets['currently_available'].sum())} of {len(targets)} "
            "validation targets are present in data/",
            if_failed=(
                "the repository sits at CONVERGENCE and cannot claim predictive validity; "
                "a converged wrong answer is still wrong"
            ),
        ),
        GateCheck(
            key="predictions_are_falsifiable",
            requirement="each prediction carries a number that would show it to be false",
            satisfied=False,
            evidence=(
                "every current test asserts the code does what it was written to do; none "
                "asserts a quantity that external data could contradict"
            ),
            if_failed=(
                "no laboratory outcome can refute the model, so any result would be "
                "absorbed as consistent"
            ),
        ),
        GateCheck(
            key="uncertainty_is_propagated",
            requirement="the output carries an interval that includes parameter uncertainty",
            satisfied=False,
            evidence=(
                "the interval reported is over seeds and scenarios; it does not include "
                "uncertainty in the 32 hypothetical parameters, because they have no "
                "distribution to propagate"
            ),
            if_failed=(
                "the reported interval is narrower than the true one, and its width is a "
                "property of the sampler rather than of the evidence"
            ),
        ),
        GateCheck(
            key="the_comparator_is_defined",
            requirement="the transfer has a stated comparator, not just an intervention",
            satisfied=False,
            evidence=(
                "REVIVED's guideline-therapy arm moves LVEF by +3.4 pp on its own; the "
                "architecture has never been compared against it in this repository"
            ),
            if_failed=(
                "a laboratory or trial reading of the result cannot separate the "
                "intervention from the background therapy everyone receives anyway"
            ),
        ),
    ]
    frame = pd.DataFrame([check.__dict__ for check in checks])
    return frame


def gate_summary() -> dict[str, object]:
    frame = transfer_gate()
    passed = int(frame["satisfied"].sum())
    return {
        "checks_total": len(frame),
        "checks_passed": passed,
        "checks_failed": len(frame) - passed,
        "first_item_to_fix": "define what 'the solution' means, before any more modelling",
        "transfer_authorised": passed == len(frame),
        "status": "transfer_gate_not_satisfied",
    }


# =============================================================================
# 4. THE DELIVERABLE: AN EXPERIMENT SPECIFICATION
# =============================================================================


def experiment_specification() -> pd.DataFrame:
    """What simulation CAN hand to a laboratory, given everything above.

    The output of a simulation-first programme is not a cure and not a validated drug. It
    is a **specification**: which measurement, at what resolution, in what population, that
    would discriminate between the surviving hypotheses. That is a real deliverable, it is
    cheap, and it is what a laboratory needs before it spends anything.
    """
    rows = [
        {
            "priority": 1,
            "measurement": "serial IL-6 and CRP time course under a defined stimulus",
            "why_this_first": (
                "the model predicts a SHAPE; a shape is the cheapest thing to falsify and "
                "the cascade currently has one observable"
            ),
            "resolves": "whether the chain structure is right at all",
            "resolution_required": ">= 6 time points, inter-individual SD reported",
            "cost_class": "one cohort, serial sampling",
        },
        {
            "priority": 2,
            "measurement": "randomised dose-response to an NLRP3 or IL-1b inhibitor",
            "why_this_second": (
                "only a randomised dose-response can reject the one-knob-per-mechanism "
                "structure, and only after the chain shape is known"
            ),
            "resolves": "the mechanism attribution, which the model cannot supply",
            "resolution_required": ">= 3 doses, placebo, randomised",
            "cost_class": "a phase 1b target-engagement study",
        },
        {
            "priority": 3,
            "measurement": "paired serial cap thickness and LDL in the same patients",
            "why_this_third": (
                "the plaque relation is asserted, not measured, and paired longitudinal "
                "imaging is the only thing that tests it"
            ),
            "resolves": "whether cap thinning follows LDL as the model assumes",
            "resolution_required": ">= 2 time points, >= 200 patients",
            "cost_class": "an imaging sub-study",
        },
        {
            "priority": 4,
            "measurement": "randomised cell-therapy arm against guideline therapy",
            "why_this_fourth": (
                "the largest claim in the repository rests on an uncontrolled observation; "
                "nothing cheaper can price a control arm"
            ),
            "resolves": "the entire efficacy question",
            "resolution_required": ">= 26 per arm for LVEF; an outcome endpoint for survival",
            "cost_class": "a randomised trial - the most expensive item, and the only one "
            "that can settle the efficacy claim",
        },
    ]
    return pd.DataFrame(rows)


# =============================================================================
# 5. REPORTING
# =============================================================================

DISCLAIMER = """\
WHAT THIS MODULE IS
  A frontier map: the requirements that must all hold before a simulation result is
  worth a laboratory's time, checked one by one. Seven of eight fail.

WHAT IT IS NOT
  A model change, a tuning exercise, or a claim. No result here is improved.

THE CENTRAL RESULT, AND ITS CORRECTION
  On terminal CRP at day 180 the inflammatory model cannot rank the four downstream
  mechanisms - they differ by 0.000252 mg/L. That measurement stands, but it used the
  wrong coordinate: interventions.terminal_value_was_the_wrong_coordinate() shows the same
  mechanisms separate by 9.95 mg/L around day 22 on a CRP time course. What blocks the
  solution is not that the model cannot rank them, but that the repository holds no time
  course on which to rank them against anyone's data.

WHAT IT CAN DO INSTEAD
  State precisely which measurement would discriminate, at what resolution, in what
  population. That is an experiment specification, and it is the correct output of a
  simulation-first programme.
"""


def frontier_report() -> str:
    ranking = ranking_is_unidentifiable()
    gate = transfer_gate()
    summary = gate_summary()
    arithmetic = validation_arithmetic()
    lines = [
        "=" * 78,
        "THE SIMULATION-FIRST FRONTIER",
        "=" * 78,
        "",
        DISCLAIMER.rstrip(),
        "",
        "-" * 78,
        "1. CAN THE MODEL NAME A WINNER?",
        "-" * 78,
    ]
    table = mechanism_effect_table()
    for row in table.itertuples(index=False):
        lines.append(
            f"  {row.mechanism:<34} final CRP {row.final_crp_mg_l:8.4f}"
            f"   reduction {row.reduction_fraction:6.2%}"
        )
    lines += [
        "",
        f"  four downstream mechanisms differ by {ranking['max_terminal_crp_spread_mg_l']:.6f} mg/L",
        f"  as a fraction of the value:            {ranking['spread_as_fraction_of_value']:.2e}",
        f"  patients per arm to separate them:     {ranking['patients_per_arm_to_separate']:,.0f}",
        f"  structurally blind:                    {ranking['structurally_blind']}",
        "",
        "  LDL lowering and NF-kB inhibition are the same perturbation:",
        "    max difference over all species and times =",
        f"    {ldl_and_nfkb_are_the_same_perturbation()['max_absolute_difference_over_all_species']:.3e}",
        "",
        "-" * 78,
        "1b. CORRECTION TO THE PREVIOUS SECTION",
        "-" * 78,
        "  The reading above uses the day-180 terminal value of an UNOBSERVED species at a",
        "  DESIGN-CHOSEN horizon. That is the wrong coordinate.",
        "  See interventions.terminal_value_was_the_wrong_coordinate():",
        "    terminal-value separation     0.000252 mg/L",
        "    best time-course separation   9.954 mg/L   (39,507x, day 22, CRP)",
        "    patients per arm required     6.9",
        "  The model CAN rank the mechanisms through a time course. It cannot rank them on",
        "  the terminal reading, and the repository holds no time course.",
        "",
        "-" * 78,
        "2. WHAT DATA WOULD MAKE EACH MODEL FALSIFIABLE",
        "-" * 78,
        f"  validation targets listed                 {arithmetic['validation_targets_listed']}",
        f"  present in data/                          {arithmetic['validation_targets_present']}",
        f"  independent observables available         {arithmetic['observables_currently_available']}",
        f"  independent observables needed            {arithmetic['minimum_observables_to_identify_all_rates']}",
        "",
        "-" * 78,
        "3. THE TRANSFER GATE",
        "-" * 78,
        f"  checks passed  {summary['checks_passed']} of {summary['checks_total']}",
        "",
    ]
    for row in gate.itertuples(index=False):
        mark = "PASS" if row.satisfied else "FAIL"
        lines.append(f"  [{mark}] {row.key}")
        lines.append(f"         requires : {row.requirement}")
        lines.append(f"         evidence : {row.evidence}")
        if not row.satisfied:
            lines.append(f"         if failed: {row.if_failed}")
        lines.append("")
    lines += [
        "-" * 78,
        "4. WHAT SIMULATION CAN HAND OVER",
        "-" * 78,
    ]
    for row in experiment_specification().itertuples(index=False):
        lines.append(f"  {row.priority}. {row.measurement}")
        lines.append(f"     resolves : {row.resolves}")
        lines.append(f"     needs    : {row.resolution_required}")
        lines.append("")
    lines += [
        "=" * 78,
        "TRANSFER AUTHORISED: " + ("YES" if summary["transfer_authorised"] else "NO"),
        "=" * 78,
        f"  First item to fix: {summary['first_item_to_fix']}",
        "",
    ]
    return "\n".join(lines)


def write_frontier_report(output_dir: str | _Path = "results/frontier") -> dict[str, object]:
    """Write the frontier artefacts. Returns a summary."""
    root = _Path(output_dir)
    tables = root / "tables"
    tables.mkdir(parents=True, exist_ok=True)

    mechanism_effect_table().to_csv(tables / "frontier_mechanism_effects.csv", index=False)
    validation_targets().to_csv(tables / "frontier_validation_targets.csv", index=False)
    transfer_gate().to_csv(tables / "frontier_transfer_gate.csv", index=False)
    experiment_specification().to_csv(tables / "frontier_experiment_spec.csv", index=False)
    what_simulation_can_still_do().to_csv(
        tables / "frontier_simulation_deliverables.csv", index=False
    )
    (root / "DISCLAIMER.txt").write_text(DISCLAIMER)
    (root / "frontier_report.txt").write_text(frontier_report())
    summary = gate_summary()
    return {
        "artefacts": 7,
        "gate_checks": summary["checks_total"],
        "gate_passed": summary["checks_passed"],
        "transfer_authorised": summary["transfer_authorised"],
        "mechanisms_comparable": len(mechanism_effect_table()),
        "ranking_separation_mg_l": ranking_is_unidentifiable()["max_terminal_crp_spread_mg_l"],
        "status": "simulation_frontier_no_model_change",
    }
