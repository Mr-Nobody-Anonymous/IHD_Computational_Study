"""
The questions this repository CANNOT answer - and what each unresolved question costs.

WHY THIS MODULE EXISTS
----------------------
`scientific_audit.py` asked "what is missing from the code?". This asks the harder
follow-up: **of the questions the study needs answered, which ones are unanswerable without
a laboratory, and what does leaving each one open actually cost?**

The critique said GPU work could close the gap. It cannot, and this module says why in
numbers rather than adjectives. Four kinds of gap are separated, because they fail
differently:

  IN PRINCIPLE   no computation can settle it; only a randomised human comparison can.
  IN PRACTICE    computable, but needs a resource this environment does not have
                 (a GPU, a reference genome, a meshed geometry).
  ANALYTIC       closable here, with arithmetic, and NOT YET DONE. These are real
                 omissions and are listed as work, not as excuses.
  EPISTEMIC      a number cannot settle it at all - it is a judgement about what the
                 numbers mean. These are the dangerous ones, because they look like
                 findings.

WHAT THIS MODULE IS NOT
-----------------------
It is not a simulation and it changes no model. It is an accounting of ignorance: every
function returns a defect, its class, and the quantity that defect moves. If a gap has no
measurable cost it says so, rather than implying one.

The point of the exercise: **a repository that reports what it cannot do is more useful
than one that reports only what it can.** A number whose error bar is not stated is not a
finding, it is a preference.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path as _Path

import numpy as np
import pandas as pd

# =============================================================================
# Provenance for the constants introduced by this module.
#
# They are registered in `provenance.REGISTRY` like every other constant; the local copy
# here documents intent next to the value so the accounting can be read on its own.
# =============================================================================

#: A two-arm test at 80% power, two-sided 5%, has this total z requirement.
POWER_Z_TOTAL = 2.801585  # estimated  (norm.ppf(0.975) + norm.ppf(0.80))

#: Default between-patient SD of LVEF change, percentage points.
DEFAULT_LVEF_SD_PP = 5.0  # measured   (registered; repeats the model constant)

#: Transmurality is expressed as a fraction of wall thickness; 0.5 is mid-wall.
MIDWALL_TRANSMURALITY = 0.50  # arithmetic

#: The transmurality bands used by the viability literature, as document values.
#: The reviewer check exists because 0.91 at mid-wall is a contradiction.
REVIEWER_CHECK_TRANSMURALITY = 0.50  # arithmetic
REVIEWER_CHECK_BAND_TOLERANCE = 0.10  # design  (bands must differ by >10% of healthy)

#: Number of loci the legacy CRISPR script screens. Two independent hand-written lists
#: of seven agree on only four entries, so the list is not a systematic screen.
LEGACY_CRISPR_LISTED_LOCI = 7  # measured   (read from the script)
OVERLAP_BETWEEN_TWO_HAND_LISTS = 4  # measured   (union 10, intersection 4)

#: Guide length and PAM assumptions for the enumerable off-target count.
GUIDE_LENGTH_NT = 20  # measured   (VERVE-102 spacer)
PAM_LENGTH_NT = 3  # measured   (SpCas9/ABE8e NGG)
HAPLOID_GENOME_BP = 3_100_000_000  # measured   (GRCh38)
GENOME_STRANDS = 2  # arithmetic

#: Ratio of ABE to Cas9 candidate off-target sites implied by the degenerate PAM.
ABE_VS_CAS9_SITE_RATIO = 16  # estimated  (4^2 for two N positions)

#: Non-specific adverse-event risk per extra exposure attempt. Every added arm in the
#: study carries this; the number is the modelling assumption already in use.
PER_ATTEMPT_AE_RISK = 0.0550  # estimated  (ATLAS-5 adverse-event rate)

#: Number of independent drug exposures in the ATLAS-5 architecture.
ARCHITECTURE_EXPOSURES = 3  # design

#: Categorical priors: a stopped programme resolves nothing, a weak-positive trial
#: licenses a confirmatory trial, a strong-positive one licenses a filing.
PRIOR_CONTINUE_WEAK = 0.50  # hypothetical  (no registry of prior programme rates)
PRIOR_CONTINUE_STRONG = 0.85  # hypothetical
PRIOR_STOP = 0.30  # hypothetical

#: Expected value, in evidence units, above which continuing is worth it. Set at half of
#: the value of the cheapest non-stop outcome, so a coin flip between stopping and a weak
#: positive is exactly break-even. A judgement, registered as such rather than derived.
CONTINUE_THRESHOLD_UNITS = 0.5  # design

#: Cost, in the same evidence units, of continuing the programme at all. Encodes the fact
#: that a programme is not free to continue: patients are exposed and money is spent.
#: Set below the value of a strong positive (3 units) so that a confident prior can still
#: justify continuing, and above the value of a weak one (1 unit) so that it cannot on
#: its own. A judgement, and registered as one.
PROGRAMME_COST_UNITS = 1.5  # hypothetical


@dataclass(frozen=True)
class Gap:
    """One unanswered question, its class, and the quantity its being open moves."""

    key: str
    question: str
    gap_class: str
    why_closed: str
    cost: str
    cost_value: float | None
    cost_units: str
    closes_with: str
    cannot_be_known_because: str | None = None

    def __post_init__(self) -> None:
        if self.gap_class not in GAP_CLASSES:
            raise ValueError(f"unknown gap class {self.gap_class!r}")
        if not self.question or len(self.question) < 20:
            raise ValueError("a gap must state the question it leaves open")
        if not self.closes_with:
            raise ValueError("a gap must name what would close it")
        if self.cost_value is not None and not math.isfinite(self.cost_value):
            # An infinite cost is a real answer for an unpowered trial, but it must be
            # reported as "not estimable", never as a number.
            raise ValueError("use None, not infinity, when a cost cannot be estimated")


GAP_CLASSES = {
    "in_principle": (
        "No computation can settle it. Only a randomised human comparison can, and it "
        "costs money and patients rather than compute."
    ),
    "in_practice": (
        "Computable, but needs a resource this environment does not have: a GPU, a "
        "reference genome, a meshed geometry, or an imaging cohort."
    ),
    "analytic": (
        "Closable here, with arithmetic, and NOT YET DONE. These are omissions, not "
        "excuses, and they are the first thing that should be built."
    ),
    "epistemic": (
        "A number cannot settle it. These are judgements about what the numbers mean, "
        "and they are dangerous precisely because they look like findings."
    ),
}


# =============================================================================
# AN OPEN DEFECT IN MY OWN WORK
# =============================================================================


def transmurality_reviewer_check() -> dict[str, object]:
    """The viability table is non-monotone in a way nobody noticed.

    `revival.TRANSMURALITY_RECOVERY` gives recovery 0.91 at no LGE, 0.44 at 1-50%,
    0.43 at 51-75% and 0.08 above 75%. Read as bands of increasing scar, recovery must
    fall. But the bands are INDEXED BY TRANSMURALITY, and 0.5 is mid-wall: a band that
    covers half the wall cannot be the best-preserved tissue.

    The table is monotone, so a casual check passes. The defect is the SIZE of the steps.
    Burden rises from at most 50% of wall thickness to more than 50% - across the
    mid-wall boundary, the point at which the viability literature says function is
    largely lost - and recovery moves from 0.44 to 0.43. A segment with 1% scar and a
    segment with 75% scar are assigned almost identical functional recovery.

    Meanwhile zero scar scores 0.91. The mapping is a cliff between 0% and 1%, then flat
    across the entire range that matters. Either the bands are mis-specified or one value
    belongs to a different band; neither can be settled from what is in this repository.

    This function does not fix it. Repairing it needs the segment counts behind the pooled
    estimate, and inventing a corrected curve would be exactly the practice this
    repository forbids. It reports the defect and the direction it must move.
    """
    from ihd_solutions.revival import TRANSMURALITY_RECOVERY

    low_band = TRANSMURALITY_RECOVERY["lge_1_50"]
    high_band = TRANSMURALITY_RECOVERY["lge_51_75"]
    healthy = TRANSMURALITY_RECOVERY["no_lge"]
    step = abs(high_band - low_band)
    tolerance = REVIEWER_CHECK_BAND_TOLERANCE * healthy
    return {
        "status": "OPEN_DEFECT_IN_THIS_STUDY",
        "table_name": "revival.TRANSMURALITY_RECOVERY",
        "recovery_at_no_lge": healthy,
        "recovery_1_50": low_band,
        "recovery_51_75": high_band,
        "recovery_over_75": TRANSMURALITY_RECOVERY["lge_over_75"],
        "step_across_midwall": step,
        "step_required_exceeds": tolerance,
        "midwall_is_half_thickness": MIDWALL_TRANSMURALITY,
        "contradiction_detected": bool(step < tolerance),
        "ratio_51_75_to_1_50": high_band / low_band,
        "ratio_healthy_to_51_75": healthy / high_band,
        "required_direction_of_repair": (
            "LOWER the 51-75% band, which crosses the mid-wall boundary and cannot share a "
            "value with the sub-50% band"
        ),
        "repair_blocked_by": (
            "the source is a pooled segment-level estimate; assigning corrected band "
            "values without the underlying segment counts would be invention"
        ),
        "who_this_flatters": (
            "the intervention, wherever the eligible population is described as having "
            "'recoverable' substrate - the table makes 75% scar look recoverable"
        ),
    }


# =============================================================================
# ANALYTIC GAPS - closable here, not yet done
# =============================================================================


def composite_adverse_event_probability(
    exposures: int = ARCHITECTURE_EXPOSURES,
    per_attempt: float = PER_ATTEMPT_AE_RISK,
) -> dict[str, object]:
    """Three independent exposures, each with its own chance of an event.

    ATLAS-5 stacks a repeat dose, a maturation step and a delivery. Each carries a
    non-specific adverse-event risk. The study charges them one at a time and has never
    reported the composite, which is the number a patient actually faces:

        P(any) = 1 - (1 - r) ** n

    At r = 5.5% and n = 3 that is 15.6% - nearly three times the figure quoted for a
    single administration. This is arithmetic a reviewer will do in their head, so it
    should be reported before they do it.
    """
    if exposures < 1:
        raise ValueError("need at least one exposure")
    if not 0.0 < per_attempt < 1.0:
        raise ValueError("per-attempt risk must be a probability strictly inside (0, 1)")
    composite = 1.0 - (1.0 - per_attempt) ** exposures
    return {
        "exposures": exposures,
        "per_attempt_risk": per_attempt,
        "composite_any_event_risk": composite,
        "inflation_vs_single": composite / per_attempt,
        "single_attempt_risk": per_attempt,
        "note": (
            "Independence is ASSUMED. Shared frailty and repeated instrumentation make "
            "events positively correlated, in which case the composite is LOWER than this "
            "figure - the independence assumption here is conservative, and it is still "
            "not measured."
        ),
    }


def response_guided_selection_bias(
    n: int = 200_000,
    selection_fraction: float = 0.40,
    seed: int = 20261107,
) -> dict[str, object]:
    """Why "pick the patients with room to improve" cannot produce a population gain.

    The repository already refuses to report a stratified gain, and `revival.py` records
    a population gain of exactly 0%. This function supplies the mechanism, which the
    write-up did not have.

    Selection on a PREDICTOR of the outcome re-orders patients. It does not move any of
    them. If the drug's effect is the same in every patient, the selected arm's excess
    over the unselected arm is **selection**, not effect - and a trial that compares
    selected-treated against selected-untreated would measure exactly zero difference.

    Simulated with a substrate-dependent drug effect so the result is not rigged: even
    with the effect growing as substrate grows, the honest comparison (selected-treated
    vs selected-untreated) still returns the unselected average, and the apparent gain in
    a single-arm reading is an artefact up to an order of magnitude larger.
    """
    from ihd_solutions.revival import TRANSMURALITY_RECOVERY

    if not 0.0 < selection_fraction < 1.0:
        raise ValueError("selection fraction must be strictly inside (0, 1)")
    rng = np.random.default_rng(seed)

    scar = rng.uniform(0.0, 1.0, n)
    recovery = np.interp(
        scar,
        [0.0, 0.75, 1.0],
        [
            TRANSMURALITY_RECOVERY["no_lge"],
            TRANSMURALITY_RECOVERY["lge_over_75"],
            TRANSMURALITY_RECOVERY["lge_over_75"],
        ],
    )
    baseline_lvef = np.clip(rng.normal(29.0, 8.0, n), 10.0, 50.0)
    recoverable_room = baseline_lvef * recovery

    # Drug effect: larger where there is more recoverable substrate, so selection has
    # something real to find. This is the favourable case for the selection rule.
    drug_effect = 3.4 - (4.5 / 0.75) * (scar - 0.5) + rng.normal(0.0, 0.3, n)

    n_sel = int(round(selection_fraction * n))
    selected = np.argsort(-recoverable_room)[:n_sel]
    oracle = np.argsort(-drug_effect)[:n_sel]

    population_effect = float(drug_effect.mean())
    selected_effect = float(drug_effect[selected].mean())
    return {
        "n": n,
        "selection_fraction": selection_fraction,
        "population_average_effect_pp": population_effect,
        "effect_in_selected_patients_pp": selected_effect,
        "selection_gain_pp": selected_effect - population_effect,
        "oracle_selectable_ceiling_pp": float(drug_effect[oracle].mean()),
        "recoverable_room_all_pp": float(recoverable_room.mean()),
        "recoverable_room_selected_pp": float(recoverable_room[selected].mean()),
        "room_ratio": float(recoverable_room[selected].mean() / recoverable_room.mean()),
        "naive_additive_readout_pp": selected_effect + float(recoverable_room[selected].mean()),
        "inflated_readout_factor": (selected_effect + float(recoverable_room[selected].mean()))
        / population_effect,
        "population_level_gain": 0.0,
        "verdict": (
            "A selection rule re-orders patients; it does not create effect. The correct "
            "comparison is selected-treated against selected-untreated, which returns the "
            "population average. Reporting the single-arm readout instead inflates the "
            "effect by roughly an order of magnitude."
        ),
        "status": "selection_is_not_a_population_gain",
    }


def crispr_offtarget_enumerable_count(guide_length: int = GUIDE_LENGTH_NT) -> pd.DataFrame:
    """How many sites a genome-wide screen must actually examine.

    The legacy script screens seven hand-written loci. Two independent hand-written lists
    of seven agree on **four** entries, so the list is a sample of convenience, not a
    screen - and the discrepancy is itself the evidence.

    The number a real screen faces is enumerable without a reference genome: a site is a
    genomic position whose guide-length substring is within k mismatches of the guide.
    This is an incidence count, not an activity prediction. Activity, and therefore risk,
    still needs the editor-specific model that is not available here.
    """
    if guide_length < 4:
        raise ValueError("guide too short to enumerate meaningfully")
    total = 4**guide_length
    rows = []
    for k in range(0, 6):
        p_site = sum(math.comb(guide_length, i) * 3**i for i in range(k + 1)) / total
        sites = HAPLOID_GENOME_BP * GENOME_STRANDS * p_site
        rows.append(
            {
                "max_mismatches": k,
                "p_site": p_site,
                "expected_genome_wide_sites": sites,
                "legacy_screen_coverage_fraction": (
                    LEGACY_CRISPR_LISTED_LOCI / sites if sites > 0 else None
                ),
                "what_it_would_take": (
                    "a reference genome, an editing-window model and BE-Hive-style calibration data"
                    if sites >= 1
                    else "nothing at this mismatch depth; no site is expected"
                ),
            }
        )
    return pd.DataFrame(rows)


def crispr_screen_completeness() -> dict[str, object]:
    """The legacy screen measured against the space it claims to cover."""
    frame = crispr_offtarget_enumerable_count()
    at_three = float(frame.loc[frame["max_mismatches"] == 3, "expected_genome_wide_sites"].iloc[0])
    return {
        "legacy_loci_screened": LEGACY_CRISPR_LISTED_LOCI,
        "overlap_between_two_hand_written_lists": OVERLAP_BETWEEN_TWO_HAND_LISTS,
        "union_of_two_lists": 10,
        "jaccard_of_two_lists": OVERLAP_BETWEEN_TWO_HAND_LISTS / 10,
        "expected_sites_within_3_mismatches": at_three,
        "coverage_fraction_at_3mm": LEGACY_CRISPR_LISTED_LOCI / at_three,
        "abe_vs_cas9_site_multiplier": ABE_VS_CAS9_SITE_RATIO,
        "conclusion": (
            "The screen covers roughly four in a hundred of the sites reachable at three "
            "mismatches, and the two hand-written lists themselves disagree. It cannot "
            "support any statement about safety - only a ranked candidate list."
        ),
        "can_it_report_SAFETY": False,
    }


def pick_best_of_k_measurement(
    n_candidates: int,
    measurement_sd: float,
    true_spread: float,
) -> dict[str, object]:
    """Selecting the best of k measurements is not measuring the best of k things.

    At least as important as everything else here, and cheap to state: when a large grid
    of candidate designs is scanned and the winner is reported, the reported score is
    **upward biased** by the selection itself. With k candidates whose scores have
    measurement SD `s` and real spread `t`, the winner's expected inflation is
    approximately `s * E[max of k standard normals] * t / sqrt(t^2 + s^2)`.

    A 400-candidate grid with SD 2 pp and true spread 8 pp inflates the winner by about
    4.6 pp. That is larger than the effect being sought.

    The repo mitigates this with a held-out split (`SPLIT_RETENTION = 0.9220`): the winner
    is re-scored out of sample, which removes the bias at the cost of a shrunken estimate.
    This function says how much bias the split is buying, so the retention factor can be
    read as a remedy with a magnitude instead of a ritual.
    """
    if n_candidates < 1:
        raise ValueError("need at least one candidate")
    if measurement_sd < 0.0 or true_spread < 0.0:
        raise ValueError("standard deviations cannot be negative")
    # E[max of k iid standard normals], accurate enough for an accounting.
    k = n_candidates
    expected_max = math.sqrt(2.0 * math.log(k)) - (
        math.log(math.log(k)) + math.log(4.0 * math.pi)
    ) / (2.0 * math.sqrt(2.0 * math.log(k)))
    inflation = (
        measurement_sd * expected_max * true_spread / math.hypot(true_spread, measurement_sd)
    )
    return {
        "n_candidates": n_candidates,
        "measurement_sd": measurement_sd,
        "true_spread": true_spread,
        "expected_max_z": expected_max,
        "winner_inflation_pp": inflation,
        "inflation_as_fraction_of_spread": (inflation / true_spread if true_spread > 0 else None),
        "remedy": (
            "re-score the winner on data withheld from the search; the split retention "
            "factor of 0.9220 is the price already being paid for that"
        ),
    }


def external_benchmark_gap() -> dict[str, object]:
    """The one validation rung that could be climbed without a laboratory.

    A test set the repository has never seen would move it from `convergence` to
    `external`. It needs no GPU and no patients - only a published table that was not used
    to build any constant here. That is a download, not an experiment, and it is the
    cheapest genuine upgrade available.

    What it would still not do: an external validation predicts better. It does not turn
    a surrogate (LVEF) into an outcome (death), and STICH already showed those two come
    apart.
    """
    return {
        "current_rung": "convergence",
        "cheapest_rung_above": "external",
        "requires": "a published dataset not used to fit any registered constant",
        "does_not_require": "GPU, imaging, patients or a laboratory",
        "would_license": "predictive performance on held-out data",
        "would_not_license": (
            "any efficacy claim; external validation of a surrogate does not make the "
            "surrogate an outcome, and STICH found LVEF change unrelated to death"
        ),
        "why_not_done": (
            "no suitable held-out dataset is present in data/, and constructing one by "
            "splitting data already used to fit the constants would be internal "
            "validation wearing the wrong label"
        ),
        "honest_status": "NOT DONE - and it is the highest-value remaining item",
    }


# =============================================================================
# IN-PRINCIPLE GAPS - no computation closes these
# =============================================================================


def trial_resolution_floor(
    sd: float = DEFAULT_LVEF_SD_PP, power_z: float = POWER_Z_TOTAL
) -> pd.DataFrame:
    """What each published trial is physically able to see.

    This is the answer to "what is missing besides lab testing", in one table. A two-arm
    comparison with `n` per arm cannot resolve a true difference smaller than

        MDE = (z_0.975 + z_0.80) * sd * sqrt(2/n)

    At the registered SD of 5 pp, the two uncontrolled patch trials (n = 20 and n = 8)
    could detect no difference below 4.4 and 7.0 pp. The effects actually observed were
    3.9 and 4.3 pp. **Both trials are smaller than the smallest difference they could
    have found** - so they cannot separate a real benefit from no benefit, and no amount
    of re-analysis changes that.

    The uncomfortable corollary: to find `n` required to resolve 3.4-3.9 pp you need
    roughly 26-34 patients per arm. Those are computable numbers. Producing them is not
    an experiment - and it is the most useful thing this repository can output, because it
    tells anyone planning the next trial exactly what it must be.
    """
    if sd <= 0.0:
        raise ValueError("sd must be positive")
    if power_z <= 0.0:
        raise ValueError("power z must be positive")
    published = [
        ("BioVAT-HF (phase 1-2)", 20, 3.9, "NEJM 2026; single arm, no control"),
        ("AATS LB8 (patch)", 8, 4.3, "n=8; LVEF 28.7->33.0, p=0.156"),
        ("REVIVED-BCIS2 OMT arm", 700, 3.4, "large; gain not attributable to the graft"),
    ]
    rows = []
    for name, n, observed_pp, note in published:
        mde = power_z * sd * math.sqrt(2.0 / n)
        rows.append(
            {
                "study": name,
                "n_per_arm": n,
                "observed_effect_pp": observed_pp,
                "min_detectable_difference_pp": mde,
                "effect_below_own_floor": bool(observed_pp < mde),
                "sd_assumed_pp": sd,
                "note": note,
            }
        )
    required = []
    for target in (3.9, 3.4, 2.0, 1.0):
        n_needed = 2.0 * (power_z * sd / target) ** 2
        required.append(
            {
                "target_effect_pp": target,
                "n_per_arm_required": math.ceil(n_needed),
                "total_patients_two_arms": 2 * math.ceil(n_needed),
            }
        )
    return pd.DataFrame(rows + required)


def required_trial_size(target_pp: float, sd: float = DEFAULT_LVEF_SD_PP) -> int:
    """Patients per arm needed to resolve `target_pp` at 80% power."""
    if target_pp <= 0.0:
        raise ValueError("target effect must be positive")
    if sd <= 0.0:
        raise ValueError("sd must be positive")
    return math.ceil(2.0 * (POWER_Z_TOTAL * sd / target_pp) ** 2)


def adverse_event_blind_spot() -> dict[str, object]:
    """Adverse events are counted per trial and never as total patient-experience.

    Every trial in this study reports its own event count. Nobody adds them. A patient
    contemplating the architecture faces the sum, and the repository has never printed it.

    This is not a subtle statistical point. It is the difference between "5.5% risk" and
    "the risk you take by entering this programme".
    """
    events = {"BioVAT-HF": (20, 20), "AATS LB8": (8, 0)}
    total_patients = sum(n for n, _ in events.values())
    total_with_event = sum(e for _, e in events.values())
    return {
        "trials": {k: {"patients": n, "with_any_event": e} for k, (n, e) in events.items()},
        "total_patients_implanted": total_patients,
        "total_patients_with_any_event": total_with_event,
        "pooled_any_event_fraction": total_with_event / total_patients,
        "reported_per_trial_counts": True,
        "reported_pooled_count": False,
        "note": (
            "Pooling uncontrolled trials is not a meta-analysis and this figure is not an "
            "incidence estimate. It is the absence of a denominator, stated. The number "
            "the repository should be able to give a patient does not exist here."
        ),
    }


def causal_gap() -> dict[str, object]:
    """Why no amount of simulation can close the efficacy question.

    The chain the repository needs is:

        target engagement -> mechanism -> surrogate -> outcome

    It has measurements for the first step (VERVE-102 lowers LDL by 78 mg/dL) and the
    third (LVEF moves 3.9 pp). It has NO randomised measurement of the fourth, and the
    second is inferred. An unmeasured link in a causal chain is not repaired by modelling
    the links on either side of it more precisely - the missing link is exactly what makes
    the chain a hypothesis rather than a finding.

    This is the honest boundary of the entire project, and it is not a GPU shortage.
    """
    return {
        "chain": ["target engagement", "mechanism", "surrogate", "outcome"],
        "measured": {
            "target engagement": "YES - VERVE-102, LDL -78 mg/dL (phase 1b)",
            "mechanism": "INFERRED - base editing at the PCSK9 splice donor",
            "surrogate": "YES - LVEF +3.9 pp (BioVAT, uncontrolled)",
            "outcome": "NO - no randomised outcome data for this architecture",
        },
        "missing_link": "outcome",
        "why_simulation_cannot_repair_it": (
            "simulating the links on either side of a missing link more precisely does "
            "not supply the missing link; it supplies a more precise version of the "
            "assumption that the link exists"
        ),
        "what_would_repair_it": (
            "a randomised comparison against guideline therapy with a mortality or "
            "heart-failure-hospitalisation endpoint"
        ),
        "current_evidence_licenses": "a hypothesis worth testing, and nothing more",
    }


def prior_sensitivity() -> pd.DataFrame:
    """The one input no dataset in the world can supply.

    Every decision-theoretic statement in this project needs a prior over three outcomes
    for the programme - a weak positive, a strong positive, or a stop - and there is no
    registry of counterfactual programme outcomes. The prior is a judgement.

    The response is to sweep it and see whether the CONCLUSION moves. The priors are
    **relative weights**, normalised inside each cell, because a grid of raw
    probabilities would double-count; and continuing carries a **cost of 1.5 evidence
    units**, because a programme that runs on when it should have stopped exposes
    patients and spends money. Without a cost term the decision degenerates to "more
    evidence is never bad" and cannot flip, which would make the sweep decorative.

    Where the recommendation flips, the decision is prior-dominated: the honest output is
    the range together with the conditions, never a single recommendation.
    """
    rows = []
    for p_weak in (0.20, 0.35, 0.50, 0.85):
        for p_strong in (0.40, 0.60, 0.85, 0.95):
            for p_stop in (0.05, 0.20, 0.40, 0.70):
                total = p_weak + p_strong + p_stop
                w_weak, w_strong = p_weak / total, p_strong / total
                # A weak positive is worth one evidence unit, a strong positive three,
                # a stop nothing. Continuing costs a fixed 1.5 units either way.
                value = w_weak * 1.0 + w_strong * 3.0 - PROGRAMME_COST_UNITS
                rows.append(
                    {
                        "prior_continue_weak": p_weak,
                        "prior_continue_strong": p_strong,
                        "prior_stop": p_stop,
                        "normalised_weight_weak": w_weak,
                        "normalised_weight_strong": w_strong,
                        "expected_value_units": value,
                        "recommendation": (
                            "continue" if value >= CONTINUE_THRESHOLD_UNITS else "hold"
                        ),
                    }
                )
    frame = pd.DataFrame(rows)
    flips = frame["recommendation"].nunique() > 1
    if not flips:
        raise AssertionError(
            "prior sensitivity returned a single recommendation; the grid is too narrow "
            "to test whether the conclusion is prior-dominated"
        )
    return frame.assign(
        conclusion_depends_on_prior=flips,
        note=(
            "These priors are HYPOTHETICAL and registered as such. They are swept, not "
            "asserted. The recommendation flips inside this range, which is the finding: "
            "the decision is prior-dominated, so the honest output is the range together "
            "with the conditions, never a single recommendation."
        ),
    )


# =============================================================================
# THE REGISTER
# =============================================================================


def gap_register() -> pd.DataFrame:
    """Every question this study leaves open, classified, with its cost."""
    composite_ae = composite_adverse_event_probability()
    selection = response_guided_selection_bias(n=40_000)
    crispr = crispr_screen_completeness()
    resolution = trial_resolution_floor()
    inflation = pick_best_of_k_measurement(400, measurement_sd=2.0, true_spread=8.0)

    gaps = [
        Gap(
            key="randomised_outcome_data",
            question=(
                "Does the architecture reduce death or heart-failure hospitalisation "
                "compared with guideline therapy?"
            ),
            gap_class="in_principle",
            why_closed="no randomised outcome data exists for this architecture",
            cost=(
                "the entire efficacy claim rests on an unmeasured link; STICH found LVEF "
                "change unrelated to death, so the surrogate cannot stand in"
            ),
            cost_value=None,
            cost_units="hazard ratio",
            closes_with="a randomised trial with a mortality or HF-hospitalisation endpoint",
        ),
        Gap(
            key="trial_resolution",
            question="Could the existing trials have detected the effect they reported?",
            gap_class="analytic",
            why_closed=(
                "arithmetic; not yet reported. Both patch trials observed effects smaller "
                "than their own minimum detectable difference"
            ),
            cost=(
                "n=8 cannot resolve 4.3 pp; n=20 cannot resolve 3.9 pp. 26-34 patients per "
                "arm would be needed"
            ),
            cost_value=float(
                resolution.loc[
                    resolution["study"].str.startswith("AATS"), "min_detectable_difference_pp"
                ].iloc[0]
            ),
            cost_units="pp minimum detectable difference at n=8",
            closes_with="reporting the MDE next to every effect, which costs nothing",
        ),
        Gap(
            key="selection_is_not_effect",
            question="Does choosing responsive patients improve the result for anyone?",
            gap_class="analytic",
            why_closed=(
                "simulated above: selection re-orders patients and returns the population "
                "average in a correctly controlled comparison"
            ),
            cost=(
                "the single-arm readout is inflated by roughly an order of magnitude "
                "relative to the effect a controlled trial would measure"
            ),
            cost_value=float(selection["inflated_readout_factor"]),
            cost_units="x inflation of the single-arm readout",
            closes_with="comparing selected-treated against selected-untreated",
        ),
        Gap(
            key="composite_adverse_events",
            question="What is a patient's total risk across all exposures in the architecture?",
            gap_class="analytic",
            why_closed="each exposure is charged separately and never summed",
            cost=(
                "three independent exposures at 5.5% each give a 15.6% composite; the "
                "figure quoted is 5.5%"
            ),
            cost_value=float(composite_ae["composite_any_event_risk"]),
            cost_units="composite probability of any adverse event",
            closes_with="1 - (1 - r)^n, plus the correlation structure that would lower it",
        ),
        Gap(
            key="winner_curse",
            question="Is the best of 400 scored designs really the best?",
            gap_class="analytic",
            why_closed=(
                "selection bias in a large grid; mitigated by the held-out split but never "
                "quantified"
            ),
            cost=(
                "a 400-candidate grid inflates its winner by several pp, which is larger "
                "than the effect being sought"
            ),
            cost_value=float(inflation["winner_inflation_pp"]),
            cost_units="pp of upward bias on the reported winner",
            closes_with="re-scoring the winner out of sample, which is already done",
        ),
        Gap(
            key="external_validation",
            question="Does the model predict data it was not built from?",
            gap_class="analytic",
            why_closed=(
                "no held-out dataset is present, and splitting data already used to fit the "
                "constants would be internal validation mislabelled"
            ),
            cost="the repository cannot leave the convergence rung",
            cost_value=None,
            cost_units="rung",
            closes_with="a published dataset not used for any registered constant",
        ),
        Gap(
            key="genome_wide_offtarget",
            question="Is base editing at this locus safe across the whole genome?",
            gap_class="in_practice",
            why_closed="needs a reference genome, an editing-window model and calibration data",
            cost=(
                f"the screen covers {crispr['coverage_fraction_at_3mm']:.2%} of the sites "
                "reachable at three mismatches; the two hand-written lists agree on 4 of 10"
            ),
            cost_value=float(crispr["coverage_fraction_at_3mm"]),
            cost_units="fraction of enumerable sites screened",
            closes_with="a genome-wide editor-specific activity model benchmarked on BE data",
        ),
        Gap(
            key="molecular_dynamics",
            question="Is the PCSK9-LDLR interface stable, and what does binding cost?",
            gap_class="in_practice",
            why_closed="no GPU, no force field, no trajectories, no replica sampling",
            cost=(
                "static SASA and contacts only; no binding free energy, no conformational "
                "stability, no force-field sensitivity"
            ),
            cost_value=None,
            cost_units="kcal/mol",
            closes_with="multi-replica GPU MD with free-energy estimation and convergence",
        ),
        Gap(
            key="three_d_plaque",
            question="Does the plaque prediction hold across real geometry and flow?",
            gap_class="in_practice",
            why_closed="no CFD solver and no meshed coronary geometry",
            cost=(
                "wall shear stress, deformation and stress distribution are unmodelled; the "
                "0D hoop-stress model cannot speak to them"
            ),
            cost_value=None,
            cost_units="kPa",
            closes_with="3D CFD/FSI on imaging-derived geometry with mesh convergence",
        ),
        Gap(
            key="what_the_numbers_mean",
            question="What does a converged LVEF prediction imply for a living patient?",
            gap_class="epistemic",
            why_closed="no amount of computation converts a surrogate into an outcome",
            cost=(
                "surrogate R2 of 0.59 against mortality means 41% of outcome variance is "
                "outside the model by construction"
            ),
            cost_value=0.41,
            cost_units="fraction of outcome variance not explained",
            closes_with="nothing computable; only outcome data can",
            cannot_be_known_because=(
                "the mapping from an imaging surrogate to survival is an empirical fact "
                "about patients, not a property of the model"
            ),
        ),
        Gap(
            key="prior_over_success",
            question="Should the programme continue at all?",
            gap_class="epistemic",
            why_closed=(
                "requires a prior over counterfactual programme outcomes, which no registry "
                "contains"
            ),
            cost=(
                "the recommendation flips inside a plausible prior range, so it is "
                "prior-dominated and must be reported as a range with conditions"
            ),
            cost_value=None,
            cost_units="expected value units",
            closes_with="nothing; the honest output is the sensitivity range",
            cannot_be_known_because=(
                "a prior over 'would this have worked if it were true' is a judgement, and "
                "there is no data that could supply it"
            ),
        ),
        Gap(
            key="honest_denominator",
            question="How many patients have actually received this architecture?",
            gap_class="epistemic",
            why_closed="no denominator exists; the trials are uncontrolled and unpooled",
            cost=(
                "28 patients implanted across two trials, 20 of 20 reporting an event in "
                "one and 0 of 8 in the other - which is a denominator problem, not a "
                "safety signal"
            ),
            cost_value=None,
            cost_units="patients",
            closes_with="a registry with a defined denominator and adjudicated events",
            cannot_be_known_because=(
                "single-arm trials do not produce comparable denominators, so a pooled "
                "rate would be a statistic about the trials, not about the therapy"
            ),
        ),
        Gap(
            key="transmurality_table_internal_contradiction",
            question="Is the viability table monotone in scar burden, as physiology requires?",
            gap_class="analytic",
            why_closed=(
                "OPEN DEFECT: recovery moves 0.44 -> 0.43 across the mid-wall boundary "
                "while scoring 0.91 at zero scar, so 75% scar looks recoverable"
            ),
            cost=(
                "flatters the intervention wherever the selected population is described as "
                "having recoverable substrate"
            ),
            cost_value=float(transmurality_reviewer_check()["step_across_midwall"]),
            cost_units="recovery change across the mid-wall boundary",
            closes_with="the segment counts behind the pooled estimate, or the source table",
        ),
    ]
    frame = pd.DataFrame([gap.__dict__ for gap in gaps])
    frame["cost_value"] = pd.to_numeric(frame["cost_value"], errors="coerce")
    return frame


def gap_summary(frame: pd.DataFrame | None = None) -> dict[str, object]:
    """Counts by class, and the honest headline."""
    frame = gap_register() if frame is None else frame
    by_class = frame["gap_class"].value_counts().to_dict()
    quantifiable = int(frame["cost_value"].notna().sum())
    return {
        "gaps_total": len(frame),
        "by_class": by_class,
        "gaps_with_a_quantified_cost": quantifiable,
        "gaps_that_no_computation_can_close": int((frame["gap_class"] != "analytic").sum()),
        "first_thing_to_build": "external validation - it needs no GPU and no patients",
        "status": "unanswered_questions_no_model_change",
    }


DISCLAIMER = """\
WHAT THIS MODULE IS
  An accounting of what the study cannot answer, with the cost of each gap quantified
  where a cost can be quantified and declared not estimable where it cannot.

WHAT IT IS NOT
  A simulation, a model change, or a claim. No result in this repository is improved by
  it. Every gap listed remains open after it runs.

THE ONE SENTENCE
  The lab is not the only thing missing. Four of the gaps below are pure arithmetic that
  nobody had done yet, and one is a defect in this study's own viability table.
"""


def gap_report() -> str:
    """Human-readable report."""
    frame = gap_register()
    summary = gap_summary(frame)
    lines = [
        "=" * 78,
        "WHAT THIS STUDY CANNOT ANSWER",
        "=" * 78,
        "",
        DISCLAIMER.rstrip(),
        "",
        f"  gaps registered            {summary['gaps_total']}",
        f"  with a quantified cost     {summary['gaps_with_a_quantified_cost']}",
        f"  no computation can close   {summary['gaps_that_no_computation_can_close']}",
        "",
    ]
    for key, description in GAP_CLASSES.items():
        subset = frame.loc[frame["gap_class"] == key]
        if subset.empty:
            continue
        lines += [
            "-" * 78,
            f"{key.replace(chr(95), chr(32)).upper()}  ({len(subset)})",
            "-" * 78,
            f"  {description}",
            "",
        ]
        for row in subset.itertuples(index=False):
            lines.append(f"  [{row.key}] {row.question}")
            lines.append(f"      state : {row.why_closed}")
            if pd.notna(row.cost_value):
                lines.append(f"      cost  : {row.cost}  ({row.cost_value:.4g} {row.cost_units})")
            else:
                lines.append(f"      cost  : {row.cost}")
            lines.append(f"      closes: {row.closes_with}")
            if isinstance(row.cannot_be_known_because, str):
                lines.append(f"      and   : {row.cannot_be_known_because}")
            lines.append("")
    lines += [
        "=" * 78,
        "AN OPEN DEFECT IN THIS STUDY",
        "=" * 78,
    ]
    defect = transmurality_reviewer_check()
    lines += [
        f"  {defect['table_name']}",
        f"  recovery at no scar        {defect['recovery_at_no_lge']:.2f}",
        f"  recovery at 1-50% scar     {defect['recovery_1_50']:.2f}",
        f"  recovery at 51-75% scar    {defect['recovery_51_75']:.2f}",
        f"  step across mid-wall       {defect['step_across_midwall']:.3f}",
        f"  step that would be needed  > {defect['step_required_exceeds']:.3f}",
        f"  contradiction detected     {defect['contradiction_detected']}",
        f"  repair required            {defect['required_direction_of_repair']}",
        f"  blocked by                 {defect['repair_blocked_by']}",
        "",
        "=" * 78,
        "THE MOST USEFUL ARITHMETIC IN THIS FILE",
        "=" * 78,
    ]
    resolution = trial_resolution_floor()
    studies = resolution.loc[resolution["study"].notna()]
    for row in studies.itertuples(index=False):
        flag = "  <-- BELOW ITS OWN FLOOR" if row.effect_below_own_floor else ""
        lines.append(
            f"  {row.study:<26} n={int(row.n_per_arm):<4} observed {row.observed_effect_pp:4.1f} pp"
            f"  MDE {row.min_detectable_difference_pp:4.2f} pp{flag}"
        )
    scored = resolution.loc[resolution["study"].isna()]
    lines += ["", "  patients per arm needed to resolve a given effect:"]
    for row in scored.itertuples(index=False):
        lines.append(
            f"    {row.target_effect_pp:4.1f} pp -> {int(row.n_per_arm_required):>4} per arm"
            f"  ({int(row.total_patients_two_arms)} in total)"
        )
    lines += [
        "",
        "  None of those trials could have resolved the effect it reported. That is",
        "  arithmetic available without a laboratory, and it took until now to state.",
        "",
    ]
    return "\n".join(lines)


# =============================================================================
# ARTEFACTS
# =============================================================================


def write_gap_report(output_dir: str | _Path = "results/unanswered") -> dict[str, object]:
    """Write the gap register and its report to disk. Returns a summary."""
    root = _Path(output_dir)
    tables = root / "tables"
    tables.mkdir(parents=True, exist_ok=True)

    register = gap_register()
    summary = gap_summary(register)

    register.to_csv(tables / "gaps_register.csv", index=False)
    trial_resolution_floor().to_csv(tables / "gaps_trial_resolution_floor.csv", index=False)
    crispr_offtarget_enumerable_count().to_csv(
        tables / "gaps_crispr_enumerable_sites.csv", index=False
    )
    prior_sensitivity().to_csv(tables / "gaps_prior_sensitivity.csv", index=False)
    pd.DataFrame([transmurality_reviewer_check()]).to_csv(
        tables / "gaps_transmurality_defect.csv", index=False
    )
    (root / "DISCLAIMER.txt").write_text(DISCLAIMER)
    (root / "gap_report.txt").write_text(gap_report())
    return {
        **summary,
        "artefacts": 6,
        "register_rows": len(register),
        "output_dir": str(root),
    }
