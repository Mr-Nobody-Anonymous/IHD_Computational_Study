"""Solutions for the six open science gaps, searched then quantified.

The previous turn listed six things "nobody has run". Searching each one properly
found that **three were wrong or overstated**:

1. "Nothing regrows heart muscle" - **false**. Jebran et al. (Nature 2025) examined an
   explanted human heart three months after engineered-heart-muscle implantation and
   found retained, vascularised, viable graft tissue. BioVAT-HF (NEJM 2026) reports
   +4.5 mm target wall thickness and +3.9 pp LVEF at three months in 16 patients at the
   maximal dose. Human remuscularisation has happened. The gap is **magnitude**, not
   existence - and this repository already carried the 3.9 pp figure.

3. "A viability-specific signal needs a reporter gene never imaged in humans" -
   **overstated**. Human NIS is non-immunogenic and has already been imaged in patients:
   MV-NIS trials confirmed NIS expression by 123I SPECT/CT in myeloma (8 of 32 scans
   positive) and ovarian cancer (3 patients, positive across repeat cycles). The
   reporter works in humans; it has simply never been put into a cardiac graft.

6. "PCSK9 editing durability is assumed" - **now measured**. VERVE-102 Heart-2
   (NEJM 2026) reports -88% PCSK9 and -62% LDL-C (-78 mg/dL) at 1.0 mg/kg, durable to
   18 months across hepatocyte turnover. Durability is no longer an assumption. What is
   still missing is outcomes - and that trial can now be *powered* on a measured effect
   rather than a guess.

For what remains genuinely open this module does the arithmetic a protocol needs:
how much graft mass or maturation time the LVEF target requires, whether the trial
that is supposed to answer the repeat-dosing question is even powered to do so, how
much a response-guided selection rule buys, and how many events an outcomes trial
needs.

NOTHING HERE IS A RESULT. Power calculations are design arithmetic under stated
assumptions, not evidence. No cure is declared.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

from ihd_solutions.units import require_fraction, require_positive

SOLUTIONS_TAG = "gap_solutions_design_arithmetic_not_evidence"

DISCLAIMER = (
    "THESE ARE EXPERIMENTAL DESIGNS AND POWER CALCULATIONS, NOT RESULTS. Sample sizes "
    "and detectable differences are arithmetic under stated assumptions; if an "
    "assumption is wrong the design is wrong. Efficacy is nowhere claimed. External "
    "findings are other people's, tagged by design and species. No cure is declared."
)

# --- anchors from the external search, with their sources -----------------------

#: Jebran et al., Nature 2025;639:503-511 and BioVAT-HF, NEJM 2026 (NEJMoa2513525).
BIOVAT_MAX_EHM_UNITS = 20
BIOVAT_CARDIOMYOCYTES_PER_UNIT = 40_000_000  # 800M across 20 units
BIOVAT_LVEF_GAIN_3MO_PP = 3.9
BIOVAT_LVEF_GAIN_17MO_PP = 6.9
BIOVAT_WALL_GAIN_3MO_MM = 4.5
BIOVAT_WALL_GAIN_12MO_MM = 2.9
BIOVAT_BASELINE_LVEF_PCT = 24.0
BIOVAT_N_AT_MAX_DOSE = 16

#: The repository's myocardial criterion, restated here so the arithmetic is local.
REQUIRED_LVEF_GAIN_PP = 20.0

#: Vrtovec/REMEDIUM: 10 of 30 responded (>=5 pp) to dose 2; responders had lower
#: pre-dose LVEF than non-responders (34.9% vs 40.0%, p=0.02).
RETRY_RESPONDERS = 10
RETRY_TOTAL = 30
RESPONDER_PREDose_LVEF = 34.9
NONRESPONDER_PREDose_LVEF = 40.0

#: CATO NCT06145035: 60 patients, 1:1:1, primary change in LVEF by MRI at 12 months.
CATO_TOTAL_N = 60
CATO_ARMS = 3

#: VERVE-102 Heart-2, NEJM 2026 (NEJMoa2601283), 1.0 mg/kg cohort.
VERVE_LDL_ABS_REDUCTION_MG_DL = 78.0
VERVE_LDL_PCT_REDUCTION = 0.62
VERVE_PCSK9_PCT_REDUCTION = 0.88
VERVE_DURABILITY_MONTHS = 18
MG_DL_PER_MMOL_L = 38.67

#: Cholesterol Treatment Trialists: 1 mmol/L LDL-C reduction over ~5 years gives
#: about 22% fewer major vascular events. Used ONLY to power a trial, never to claim
#: the benefit.
CTT_RR_PER_MMOL_L = 0.78

#: Precedented clearance potency cap, from the ATLAS-2 work.
PRECEDENTED_POTENCY_CAP = 1.284


class SolutionStatus(str, Enum):
    CLOSED_BY_OTHERS = "closed_by_others"
    NARROWED = "narrowed"
    STILL_OPEN = "still_open"
    AWAITING_READOUT = "awaiting_readout"


@dataclass(frozen=True)
class GapSolution:
    """One gap, what the search did to it, and the concrete next action."""

    gap_id: int
    prior_claim: str
    search_finding: str
    status: SolutionStatus
    solution: str
    primary_endpoint: str
    blockers: str
    falsifier: str
    source: str

    @property
    def was_overstated(self) -> bool:
        return self.status in (SolutionStatus.CLOSED_BY_OTHERS, SolutionStatus.NARROWED)


# --- design arithmetic ---------------------------------------------------------


def detectable_difference(
    n_per_arm: int, sd: float, alpha: float = 0.05, power: float = 0.8, comparisons: int = 1
) -> float:
    """Smallest true mean difference a two-arm comparison can detect.

    Bonferroni-adjusted across `comparisons` so a three-arm trial is not credited with
    a two-arm alpha.
    """
    require_positive("n_per_arm", n_per_arm)
    require_positive("sd", sd)
    require_fraction("alpha", alpha)
    require_fraction("power", power)
    if n_per_arm < 2:
        raise ValueError("need at least two per arm")
    z_alpha = float(norm.ppf(1.0 - alpha / (2.0 * comparisons)))
    z_power = float(norm.ppf(power))
    return float((z_alpha + z_power) * sd * np.sqrt(2.0 / n_per_arm))


def power_at_difference(
    n_per_arm: int,
    sd: float,
    difference: float,
    alpha: float = 0.05,
    comparisons: int = 1,
) -> float:
    """Power to detect a given true difference."""
    require_positive("n_per_arm", n_per_arm)
    require_positive("sd", sd)
    z_alpha = float(norm.ppf(1.0 - alpha / (2.0 * comparisons)))
    noncentrality = abs(difference) / (sd * np.sqrt(2.0 / n_per_arm))
    return float(norm.cdf(noncentrality - z_alpha))


def n_per_arm_for_difference(
    sd: float, difference: float, alpha: float = 0.05, power: float = 0.8, comparisons: int = 1
) -> int:
    """Arm size needed to detect a difference, rounded up."""
    require_positive("sd", sd)
    require_positive("difference", difference)
    z_alpha = float(norm.ppf(1.0 - alpha / (2.0 * comparisons)))
    z_power = float(norm.ppf(power))
    return int(np.ceil(2.0 * ((z_alpha + z_power) * sd / difference) ** 2))


def remuscularisation_scaling() -> pd.DataFrame:
    """Gap 1: the gap is magnitude, not existence. What closes the magnitude?

    Two candidate levers, both computed from the same human dataset:

    * **Mass** - scale the graft. At the three-month effect this needs 5.13x the
      demonstrated cell dose.
    * **Maturation time** - the same graft improved from +3.9 pp at three months to
      +6.9 pp at a mean 17 months. Waiting buys 1.77x for free and cuts the mass
      requirement to 2.90x.

    The discordance is the interesting part and is reported rather than smoothed:
    target wall thickness *fell* from +4.5 mm to +2.9 mm over the same interval while
    LVEF *rose*. Function improved as mass declined, which points at graft maturation
    and mechanical entrainment rather than bulk. If that reading is right, mass
    scaling is the wrong lever and the table's mass column is an upper bound.
    """
    rows = []
    for label, gain, months in (
        ("3_month_effect", BIOVAT_LVEF_GAIN_3MO_PP, 3.0),
        ("17_month_effect", BIOVAT_LVEF_GAIN_17MO_PP, 17.0),
    ):
        multiplier = REQUIRED_LVEF_GAIN_PP / gain
        rows.append(
            {
                "anchor": label,
                "follow_up_months": months,
                "demonstrated_lvef_gain_pp": gain,
                "required_lvef_gain_pp": REQUIRED_LVEF_GAIN_PP,
                "required_multiplier": multiplier,
                "ehm_units_if_mass_scaled": multiplier * BIOVAT_MAX_EHM_UNITS,
                "cardiomyocytes_if_mass_scaled": multiplier
                * BIOVAT_MAX_EHM_UNITS
                * BIOVAT_CARDIOMYOCYTES_PER_UNIT,
                "status": SOLUTIONS_TAG,
            }
        )
    return pd.DataFrame(rows)


def maturation_benefit() -> dict[str, float]:
    """How much of the magnitude gap closes by waiting instead of scaling."""
    early = BIOVAT_LVEF_GAIN_3MO_PP
    late = BIOVAT_LVEF_GAIN_17MO_PP
    return {
        "maturation_multiplier": late / early,
        "mass_multiplier_at_3_months": REQUIRED_LVEF_GAIN_PP / early,
        "mass_multiplier_at_17_months": REQUIRED_LVEF_GAIN_PP / late,
        "gap_closed_by_time_alone_fraction": 1.0
        - (REQUIRED_LVEF_GAIN_PP / late) / (REQUIRED_LVEF_GAIN_PP / early),
        "wall_thickness_change_mm": BIOVAT_WALL_GAIN_12MO_MM - BIOVAT_WALL_GAIN_3MO_MM,
        "lvef_change_pp": late - early,
    }


def cato_power_analysis() -> pd.DataFrame:
    """Gap 4: is the trial that answers repeat dosing actually powered to answer it?

    CATO randomises 60 patients 1:1:1, so 20 per arm, with two active-versus-control
    comparisons. Cardiac-MRI change in LVEF has a between-patient SD of roughly 5-8 pp;
    that range is swept rather than asserted. The repeat-versus-single increment
    reported across the human literature is small, so the honest question is whether
    20 per arm can see it.
    """
    n_per_arm = CATO_TOTAL_N // CATO_ARMS
    rows = []
    for sd in (5.0, 6.0, 7.0, 8.0):
        detectable = detectable_difference(n_per_arm, sd, comparisons=CATO_ARMS - 1)
        rows.append(
            {
                "n_per_arm": n_per_arm,
                "assumed_sd_pp": sd,
                "comparisons": CATO_ARMS - 1,
                "minimum_detectable_difference_pp": detectable,
                "power_for_3pp": power_at_difference(n_per_arm, sd, 3.0, comparisons=CATO_ARMS - 1),
                "power_for_5pp": power_at_difference(n_per_arm, sd, 5.0, comparisons=CATO_ARMS - 1),
                "n_per_arm_needed_for_3pp": n_per_arm_for_difference(
                    sd, 3.0, comparisons=CATO_ARMS - 1
                ),
                "status": SOLUTIONS_TAG,
            }
        )
    return pd.DataFrame(rows)


def response_guided_enrichment() -> pd.DataFrame:
    """Gap 2: what a pre-dose LVEF selection rule buys, and the trial it implies.

    Vrtovec's responder analysis is the whole basis: patients who gained >=5 pp from a
    second dose had lower LVEF immediately before it (34.9% vs 40.0%). That is a
    selection signal available in any clinic today. The between-patient SD is not
    reported, so it is swept; the enrichment is therefore a range, not a number.

    The design this implies is not a new molecule. It is the existing repeat-dose
    protocol with an entry rule at the second dose.
    """
    base_rate = RETRY_RESPONDERS / RETRY_TOTAL
    rows = []
    for sd in (3.0, 4.0, 5.0):
        for cutoff in (36.0, 37.5, 39.0):
            # Responder and non-responder pre-dose LVEF distributions, normal with a
            # common swept SD. Enrichment = P(responder | LVEF <= cutoff).
            p_screen = base_rate * float(norm.cdf((cutoff - RESPONDER_PREDose_LVEF) / sd)) + (
                1.0 - base_rate
            ) * float(norm.cdf((cutoff - NONRESPONDER_PREDose_LVEF) / sd))
            if p_screen <= 0.0:
                continue
            enriched = (
                base_rate * float(norm.cdf((cutoff - RESPONDER_PREDose_LVEF) / sd)) / p_screen
            )
            rows.append(
                {
                    "assumed_sd_pp": sd,
                    "cutoff_lvef_pp": cutoff,
                    "baseline_responder_rate": base_rate,
                    "fraction_passing_screen": p_screen,
                    "enriched_responder_rate": enriched,
                    "enrichment_fold": enriched / base_rate,
                    "screen_failures_per_responder_enrolled": (1.0 - p_screen)
                    / max(p_screen, 1e-9),
                    "status": SOLUTIONS_TAG,
                }
            )
    return pd.DataFrame(rows)


def pcsk9_outcomes_power() -> pd.DataFrame:
    """Gap 6: power the outcomes trial on the measured LDL effect, not a guess.

    VERVE-102 at 1.0 mg/kg lowered LDL-C by 78 mg/dL, which is 2.02 mmol/L. Applying the
    CTT meta-analytic slope log-linearly gives an assumed hazard ratio; that assumption
    is the weak link and is swept, with a deliberately pessimistic half-effect arm.

    This is design arithmetic. It does NOT assert that editing delivers this event
    reduction - asserting that is exactly what the outcomes trial exists to test.
    """
    mmol = VERVE_LDL_ABS_REDUCTION_MG_DL / MG_DL_PER_MMOL_L
    rows = []
    for label, fraction_of_ctt in (
        ("ctt_slope_full", 1.0),
        ("ctt_slope_three_quarters", 0.75),
        ("ctt_slope_half", 0.5),
    ):
        hazard = float(CTT_RR_PER_MMOL_L ** (mmol * fraction_of_ctt))
        events = int(
            np.ceil(
                4.0 * (float(norm.ppf(0.975)) + float(norm.ppf(0.8))) ** 2 / np.log(hazard) ** 2
            )
        )
        rows.append(
            {
                "assumption": label,
                "ldl_reduction_mmol_l": mmol * fraction_of_ctt,
                "assumed_hazard_ratio": hazard,
                "events_required": events,
                "n_at_3pct_annual_event_rate_3yr": int(np.ceil(events / (0.03 * 3.0))),
                "n_at_2pct_annual_event_rate_3yr": int(np.ceil(events / (0.02 * 3.0))),
                "status": SOLUTIONS_TAG,
            }
        )
    return pd.DataFrame(rows)


def reporter_signal_options() -> pd.DataFrame:
    """Gap 3: the reporter already works in humans. What is missing is the cardiac graft.

    Human NIS is the option that matters because it is human-derived and therefore
    expected to be non-immunogenic, and because it has already been imaged in patients
    with clinical-grade tracers. HSV1-tk is viral and gave lower cardiac contrast.
    """
    return pd.DataFrame(
        [
            {
                "reporter": "human_NIS",
                "tracers": "123I | 124I | 99mTc-pertechnetate | 18F-tetrafluoroborate",
                "viability_specific": True,
                "human_imaging_precedent": (
                    "MV-NIS: 123I SPECT/CT positive in 8 of 32 myeloma patients; "
                    "NIS expression confirmed in 3 ovarian-cancer patients across "
                    "repeat cycles"
                ),
                "immunogenicity": "human-derived, expected non-immunogenic",
                "cardiac_graft_precedent": "none - this is the gap",
                "status": SOLUTIONS_TAG,
            },
            {
                "reporter": "HSV1-tk",
                "tracers": "18F-FHBG",
                "viability_specific": True,
                "human_imaging_precedent": "used in human gene-transfer studies",
                "immunogenicity": "viral protein, immunogenic",
                "cardiac_graft_precedent": (
                    "rodent only; gave LOWER cardiac signal and contrast than hNIS "
                    "in a direct comparison (n=34 animals)"
                ),
                "status": SOLUTIONS_TAG,
            },
            {
                "reporter": "iron_oxide_MRI",
                "tracers": "none (contrast agent)",
                "viability_specific": False,
                "human_imaging_precedent": "yes",
                "immunogenicity": "not applicable",
                "cardiac_graft_precedent": (
                    "signal persists in macrophages after cell death - blind to viability"
                ),
                "status": SOLUTIONS_TAG,
            },
        ]
    )


def potency_redundancy_tradeoff(draws: int = 800, seed: int = 20261107) -> pd.DataFrame:
    """Gap 5: can redundancy buy down the beyond-precedent potency requirement?

    ATLAS-2 needed a 2.80x cap against a precedented 1.284x. Its ablation showed
    titration was worth more than any other lever, which suggests the potency
    requirement is partly substitutable by architecture. This sweeps the cap against
    node and gate redundancy to find out how far, rather than assuming it is not.
    """
    from ihd_solutions.upgrade import SOLVED_MONTHS, Upgrades, evaluate_upgrades

    rows = []
    for cap in (PRECEDENTED_POTENCY_CAP, 1.6, 1.8, 2.0, 2.4, 2.8):
        for nodes in (1, 2, 3, 4):
            for gates in (1, 2, 3):
                config = Upgrades(
                    gate_paths=gates,
                    node_paths=nodes,
                    titrate=True,
                    verified_source=False,
                    max_potency=cap,
                )
                result = evaluate_upgrades(config, months=SOLVED_MONTHS, draws=draws, seed=seed)
                rows.append(
                    {
                        "max_potency": cap,
                        "node_paths": nodes,
                        "gate_paths": gates,
                        "probability_of_target": result.probability_of_target,
                        "robust": result.probability_of_target >= 0.8,
                        "robust_with_margin": (
                            result.probability_of_target
                            - 2.0
                            * float(
                                np.sqrt(
                                    max(
                                        result.probability_of_target
                                        * (1.0 - result.probability_of_target),
                                        0.0,
                                    )
                                    / max(draws, 1)
                                )
                            )
                            >= 0.8
                        ),
                        "precedent_gap_fold": cap / PRECEDENTED_POTENCY_CAP,
                        "redundancy_beyond_precedent": config.redundancy_beyond_precedent,
                        "status": SOLUTIONS_TAG,
                    }
                )
    return pd.DataFrame(rows)


GAP_SOLUTIONS: tuple[GapSolution, ...] = (
    GapSolution(
        gap_id=1,
        prior_claim="Nothing regrows heart muscle; nobody has replaced it in a human.",
        search_finding=(
            "FALSE. Jebran et al. Nature 2025;639:503-511 examined an explanted human "
            "heart 3 months after EHM implantation and found retained, vascularised, "
            "viable graft. BioVAT-HF (NEJM 2026, NEJMoa2513525): 16 patients at 20 EHM "
            "units / 800M iPSC-CMs, target wall +4.5 mm and LVEF +3.9 pp at 3 months, "
            "+6.9 pp at mean 17 months. No graft-related arrhythmia or tumour."
        ),
        status=SolutionStatus.NARROWED,
        solution=(
            "Stop asking whether remuscularisation is possible and ask how to reach "
            "+20 pp. Two levers, both from the same human data: mass scaling needs "
            "5.13x the demonstrated dose at 3 months, but maturation alone lifts the "
            "effect 1.77x and cuts the mass requirement to 2.90x. Because wall "
            "thickness FELL (4.5 to 2.9 mm) while LVEF ROSE (3.9 to 6.9 pp), maturation "
            "and mechanical entrainment look like the real lever. Protocol: randomise "
            "maximal-dose EHM recipients to standard versus maturation-optimised "
            "grafts (longer in-vitro conditioning, paced load training) and follow to "
            "24 months, with the existing patch as active control."
        ),
        primary_endpoint="change in LVEF by cardiac MRI at 24 months, active-controlled",
        blockers=(
            "Lifelong immunosuppression; epicardial surface area limits coverage; "
            "vascularisation depth; open-label uncontrolled precedent with n=16 and "
            "surrogate endpoints; 196 AEs / 57 serious in 20 patients."
        ),
        falsifier=(
            "If LVEF gain does not keep rising after 12 months, maturation is not the "
            "lever and the 5.13x mass requirement stands."
        ),
        source="https://www.nejm.org/doi/abs/10.1056/NEJMoa2513525",
    ),
    GapSolution(
        gap_id=2,
        prior_claim="Response-guided retry timing was proposed in 2019 and never tested.",
        search_finding=(
            "STILL TRUE, and the selection signal is already published: Vrtovec's "
            "responders to a second dose had lower pre-dose LVEF (34.9% vs 40.0%, "
            "p=0.02). No trial has used it as an entry rule."
        ),
        status=SolutionStatus.STILL_OPEN,
        solution=(
            "No new molecule and no new device. Take the existing repeat-dose protocol "
            "and add an entry rule at the second dose: re-dose only if pre-dose LVEF "
            "is below a threshold, and randomise the threshold. This converts an "
            "arbitrary calendar interval into a measured decision, and the enrichment "
            "arithmetic here sizes it."
        ),
        primary_endpoint=(
            "proportion gaining >=5 pp LVEF after the second dose, guided versus fixed-interval"
        ),
        blockers=(
            "The pre-dose LVEF SD is unreported, so enrichment is a swept range; "
            "screening failures cost patients and time; the rule needs prospective "
            "validation before it can guide care."
        ),
        falsifier=(
            "If enrichment does not materialise prospectively, the Vrtovec subgroup "
            "was chance and fixed-interval dosing is as good."
        ),
        source="https://pmc.ncbi.nlm.nih.gov/articles/PMC6400294/",
    ),
    GapSolution(
        gap_id=3,
        prior_claim=(
            "A viability-specific engraftment signal in humans is not done; the animal "
            "version needs an engineered reporter gene."
        ),
        search_finding=(
            "OVERSTATED. Human NIS has already been imaged in patients: MV-NIS gave "
            "123I SPECT/CT positives in 8 of 32 myeloma patients and confirmed NIS "
            "expression in 3 ovarian-cancer patients across repeat cycles. hNIS is "
            "human-derived and expected non-immunogenic, and beat HSV1-tk for cardiac "
            "signal and contrast in a direct comparison."
        ),
        status=SolutionStatus.NARROWED,
        solution=(
            "The graft cells are already an engineered allogeneic iPSC line made under "
            "GMP, so adding a stable hNIS cassette is a manufacturing change, not a new "
            "platform. Image with 123I SPECT/CT or 18F-tetrafluoroborate PET. Both "
            "components are individually precedented in humans; the combination in a "
            "cardiac graft is not. Protocol: hNIS-tagged versus untagged graft in a "
            "small safety cohort, with signal-versus-function correlation as the "
            "endpoint, replicating the rodent rho=0.76-0.80 in humans."
        ),
        primary_endpoint=(
            "correlation between day-7 NIS signal and 6-month LVEF change; safety of "
            "the tagged product"
        ),
        blockers=(
            "Thyroid, stomach and salivary glands express NIS natively, so background "
            "must be blocked with perchlorate; reporter silencing over time; a second "
            "genetic modification adds regulatory burden."
        ),
        falsifier=(
            "If signal does not track function in humans as it does in rodents, the "
            "feedback loop has nothing to close on."
        ),
        source="https://cancerres.aacrjournals.org/content/75/1/22",
    ),
    GapSolution(
        gap_id=4,
        prior_claim="Whether repeat dosing works in humans awaits CATO; nothing to do but wait.",
        search_finding=(
            "STILL AWAITING, but waiting may not settle it. CATO has 20 per arm and two "
            "active-versus-control comparisons. At a plausible MRI LVEF-change SD of "
            "6 pp it can only detect roughly a 6 pp difference, and has well under 50% "
            "power for the ~3 pp increment the human literature suggests."
        ),
        status=SolutionStatus.AWAITING_READOUT,
        solution=(
            "Read CATO out as planned, but pre-register the interpretation now: a null "
            "result at 20 per arm cannot distinguish 'no effect' from 'an effect too "
            "small for this sample'. Pool CATO with the other repeat-dose trials in a "
            "prospective meta-analysis, and use the enrichment rule from gap 2 to "
            "power the successor trial rather than repeating an underpowered design."
        ),
        primary_endpoint="CATO: change in LVEF by cardiac MRI at 12 months (already specified)",
        blockers=(
            "Follow-up completes March 2027; the SD assumption drives the power "
            "estimate and is swept, not known; pooling needs harmonised endpoints."
        ),
        falsifier=(
            "If the observed between-patient SD is at the low end of the sweep, CATO is "
            "better powered than this analysis suggests and the concern is overstated."
        ),
        source="https://clinicaltrials.gov/study/NCT06145035",
    ),
    GapSolution(
        gap_id=5,
        prior_claim="The clearance solve needs 2.18x beyond-precedent potency across every module.",
        search_finding=(
            "Partly addressable by architecture rather than potency. ATLAS-2's own "
            "ablation found titration worth +72.8 pp, more than any other lever, which "
            "implies the potency requirement is substitutable. The trade-off is swept "
            "here instead of assumed."
        ),
        status=SolutionStatus.STILL_OPEN,
        solution=(
            "Buy the requirement down with redundancy instead of potency: sweep the "
            "potency cap against node and gate redundancy and report the cheapest "
            "combination that stays robust. Whatever residual gap remains is the honest "
            "specification, and it is smaller than 2.18x if redundancy substitutes at "
            "all. Note that VERVE-102 now gives -62% LDL-C from one infusion, which "
            "raises the precedented ceiling for the driver module specifically."
        ),
        primary_endpoint="simulated probability of the clearance target, with the precedent gap stated",
        blockers=(
            "Redundant routes mean more agents, more interactions and more cost; "
            "simulation cannot establish that any of these combinations is safe."
        ),
        falsifier=(
            "If no combination below the 2.80x cap is robust, potency is not "
            "substitutable and the 2.18x gap stands as written."
        ),
        source="https://www.nejm.org/doi/full/10.1056/NEJMoa2601283",
    ),
    GapSolution(
        gap_id=6,
        prior_claim="A PCSK9-editing outcomes trial is needed; durability is assumed.",
        search_finding=(
            "DURABILITY IS NOW MEASURED. VERVE-102 Heart-2 (NEJM 2026): -88% PCSK9, "
            "-62% LDL-C, -78 mg/dL absolute at 1.0 mg/kg, durable to 18 months with 15 "
            "participants beyond a year - through hepatocyte turnover of 200-300 days. "
            "No dose-limiting toxicity. Still no outcomes data."
        ),
        status=SolutionStatus.NARROWED,
        solution=(
            "Power the outcomes trial on the measured 78 mg/dL rather than a guess. "
            "That is 2.02 mmol/L; on the CTT slope it implies a hazard ratio near 0.61, "
            "which needs on the order of 130 events - a few thousand patient-years, not "
            "a mega-trial. Sweep the slope down to half effect so the design survives "
            "being wrong. Randomise editing versus continued standard lipid lowering "
            "in established coronary disease."
        ),
        primary_endpoint="MACE on hard endpoints, editing versus standard care",
        blockers=(
            "The CTT slope is an assumption imported from statins, not an editing "
            "result; irreversible editing demands long safety follow-up; the highest "
            "dose cohort has the shortest follow-up."
        ),
        falsifier=(
            "If the event curve does not separate on the assumed slope, the "
            "surrogate-to-outcome translation fails for editing - which is precisely "
            "what the trial is for."
        ),
        source="https://www.nejm.org/doi/full/10.1056/NEJMoa2601283",
    ),
)


def solutions_frame() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "gap_id": s.gap_id,
                "status_of_claim": s.status.value,
                "prior_claim_was_overstated": s.was_overstated,
                "prior_claim": s.prior_claim,
                "search_finding": s.search_finding,
                "solution": s.solution,
                "primary_endpoint": s.primary_endpoint,
                "blockers": s.blockers,
                "falsifier": s.falsifier,
                "source": s.source,
                "tag": SOLUTIONS_TAG,
            }
            for s in GAP_SOLUTIONS
        ]
    )


def solutions_report(
    scaling: pd.DataFrame,
    maturation: dict[str, float],
    cato: pd.DataFrame,
    enrichment: pd.DataFrame,
    outcomes: pd.DataFrame,
    tradeoff: pd.DataFrame,
) -> str:
    # A potency breach may not be paid for with an unflagged redundancy breach.
    # The admissible set requires BOTH to be within precedent on the redundancy axis.
    robust = tradeoff.loc[
        tradeoff["robust_with_margin"] & (~tradeoff["redundancy_beyond_precedent"])
    ]
    cheapest = None
    if not robust.empty:
        # Cheapest cap first, then the MOST headroom at that cap. Sorting by fewest
        # paths instead selects a point sitting on the 80% bar, which is the flaky
        # "solve exactly on the threshold" failure this repository already retracted.
        cheapest = robust.sort_values(
            ["max_potency", "probability_of_target"], ascending=[True, False]
        ).iloc[0]
    any_robust = tradeoff.loc[tradeoff["robust_with_margin"]]
    cheating = None
    if not any_robust.empty:
        cheating = any_robust.sort_values(
            ["max_potency", "probability_of_target"], ascending=[True, False]
        ).iloc[0]
    mid_cato = cato.loc[cato["assumed_sd_pp"].eq(6.0)].iloc[0]
    best_enrichment = enrichment.sort_values("enrichment_fold", ascending=False).iloc[0]
    full_slope = outcomes.loc[outcomes["assumption"].eq("ctt_slope_full")].iloc[0]
    half_slope = outcomes.loc[outcomes["assumption"].eq("ctt_slope_half")].iloc[0]
    lines = [
        "SOLUTIONS TO THE SIX OPEN GAPS",
        "=" * 74,
        DISCLAIMER,
        "",
        "THE SEARCH OVERTURNED THREE OF MY SIX CLAIMS",
    ]
    for entry in GAP_SOLUTIONS:
        marker = "OVERSTATED" if entry.was_overstated else "stands"
        lines.append(f"  gap {entry.gap_id}: {entry.status.value:18s} ({marker})")
    lines += [
        "",
        "GAP 1 - remuscularisation already happened in a human. The gap is magnitude.",
        f"  demonstrated at 3 months           +{BIOVAT_LVEF_GAIN_3MO_PP:.1f} pp LVEF",
        f"  demonstrated at ~17 months         +{BIOVAT_LVEF_GAIN_17MO_PP:.1f} pp LVEF",
        f"  required                           +{REQUIRED_LVEF_GAIN_PP:.1f} pp",
    ]
    for row in scaling.itertuples(index=False):
        lines.append(
            f"  {row.anchor:20s} needs {row.required_multiplier:5.2f}x  "
            f"= {row.ehm_units_if_mass_scaled:6.1f} EHM units  "
            f"= {row.cardiomyocytes_if_mass_scaled / 1e9:.2f} billion cardiomyocytes"
        )
    lines += [
        "",
        f"  Maturation alone buys {maturation['maturation_multiplier']:.2f}x and closes "
        f"{100 * maturation['gap_closed_by_time_alone_fraction']:.1f}% of the mass requirement.",
        f"  But wall thickness FELL {maturation['wall_thickness_change_mm']:+.1f} mm while LVEF "
        f"ROSE {maturation['lvef_change_pp']:+.1f} pp.",
        "  Function improved as mass declined. That points at maturation and mechanical",
        "  entrainment, not bulk - so mass scaling is probably the wrong lever and the",
        "  units above are an upper bound.",
        "",
        "GAP 2 - response-guided retry: no new molecule, just an entry rule.",
        f"  baseline responder rate            {RETRY_RESPONDERS}/{RETRY_TOTAL} = "
        f"{100 * RETRY_RESPONDERS / RETRY_TOTAL:.1f}%",
        f"  best swept enrichment              {best_enrichment.enrichment_fold:.2f}x "
        f"at LVEF <= {best_enrichment.cutoff_lvef_pp:.1f} pp",
        f"  enriched responder rate            {100 * best_enrichment.enriched_responder_rate:.1f}%",
        f"  fraction passing the screen        {100 * best_enrichment.fraction_passing_screen:.1f}%",
        "  The SD behind this is unreported, so the whole table is swept. It is a range,",
        "  not a promise.",
        "",
        "GAP 3 - the reporter already works in humans; put it in the graft.",
        "  Human NIS imaged in patients: 123I SPECT/CT positive in 8 of 32 myeloma",
        "  patients; NIS expression confirmed in 3 ovarian-cancer patients across repeat",
        "  cycles. hNIS is human-derived and beat HSV1-tk for cardiac contrast.",
        "  The graft is already an engineered GMP iPSC line, so this is a manufacturing",
        "  change. Block native uptake with perchlorate; watch for reporter silencing.",
        "",
        "GAP 4 - the trial meant to answer this may not be able to.",
        f"  CATO has {CATO_TOTAL_N // CATO_ARMS} per arm, {CATO_ARMS - 1} comparisons vs control.",
        f"  At SD 6 pp it can detect only {mid_cato.minimum_detectable_difference_pp:.2f} pp.",
        f"  Power for a 3 pp increment:        {100 * mid_cato.power_for_3pp:.1f}%",
        f"  Power for a 5 pp increment:        {100 * mid_cato.power_for_5pp:.1f}%",
        f"  Arm size needed for 3 pp:          {int(mid_cato.n_per_arm_needed_for_3pp)}",
        "",
        "  So a null CATO cannot distinguish 'no effect' from 'too small to see here'.",
        "  Pre-register that interpretation now, pool prospectively, and power the",
        "  successor with the gap-2 enrichment rule.",
        "",
        "GAP 5 - buy the potency requirement down with redundancy.",
    ]
    if cheapest is not None:
        lines += [
            f"  cheapest robust point: cap {cheapest.max_potency:.3f}x, "
            f"{int(cheapest.node_paths)} node paths, {int(cheapest.gate_paths)} gate paths",
            f"  P(target) {100 * cheapest.probability_of_target:.2f}%, "
            f"precedent gap {cheapest.precedent_gap_fold:.2f}x",
        ]
    else:
        lines.append("  no robust point with precedented redundancy - the potency gap stands")
    if cheapest is not None:
        lines += [
            f"  potency precedent gap reduced      {cheapest.precedent_gap_fold:.2f}x "
            f"(was {2.80 / PRECEDENTED_POTENCY_CAP:.2f}x at the ATLAS-2 cap)",
        ]
    if cheating is not None and bool(cheating.redundancy_beyond_precedent):
        lines += [
            "",
            f"  A cheaper-looking point exists at cap {cheating.max_potency:.3f}x with "
            f"{int(cheating.node_paths)} node paths, but that",
            "  redundancy is ITSELF beyond precedent. Paying for a potency breach with an",
            "  unflagged redundancy breach is not a saving, so it is excluded above.",
        ]
    lines += [
        "",
        "GAP 6 - durability is measured; power the outcomes trial on it.",
        f"  VERVE-102 at 1.0 mg/kg: -{VERVE_LDL_ABS_REDUCTION_MG_DL:.0f} mg/dL = "
        f"{full_slope.ldl_reduction_mmol_l:.2f} mmol/L, durable {VERVE_DURABILITY_MONTHS} months",
        f"  full CTT slope:   assumed HR {full_slope.assumed_hazard_ratio:.3f}, "
        f"{int(full_slope.events_required)} events, "
        f"~{int(full_slope.n_at_3pct_annual_event_rate_3yr)} patients at 3%/yr over 3 yr",
        f"  half CTT slope:   assumed HR {half_slope.assumed_hazard_ratio:.3f}, "
        f"{int(half_slope.events_required)} events, "
        f"~{int(half_slope.n_at_3pct_annual_event_rate_3yr)} patients",
        "",
        "  These are event counts under an imported statin slope, NOT a claim that",
        "  editing delivers them. Testing that translation is the point of the trial.",
        "",
        "WHAT CHANGED IN THE ANSWER",
        "  Before: six open gaps, the deepest being that nothing regrows heart muscle.",
        "  After:  human remuscularisation is proven, the reporter is proven in humans,",
        "          and PCSK9-editing durability is measured. Three of six were wrong or",
        "          overstated.",
        "",
        "  What is genuinely left is narrower and mostly arithmetic:",
        "    - roughly 2.9x more graft effect than demonstrated, and maturation may",
        "      supply part of it for free",
        "    - an entry rule nobody has tried, using a signal already in the clinic",
        "    - a reporter cassette already imaged in patients, never put in a heart",
        "    - a readout that may be underpowered, and a successor design if so",
        "    - an outcomes trial needing on the order of a hundred-odd events",
        "",
        "None of this is a cure. Every number above is design arithmetic under stated",
        "assumptions, and every external finding belongs to someone else.",
    ]
    return "\n".join(lines) + "\n"


def run_gap_solutions_study(
    output_dir: str | Path = "results/gap_solutions", draws: int = 1200, seed: int = 20261107
) -> dict[str, object]:
    root = Path(output_dir)
    tables = root / "tables"
    reports = root / "reports"
    tables.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)

    scaling = remuscularisation_scaling()
    maturation = maturation_benefit()
    cato = cato_power_analysis()
    enrichment = response_guided_enrichment()
    outcomes = pcsk9_outcomes_power()
    options = reporter_signal_options()
    tradeoff = potency_redundancy_tradeoff(draws=draws, seed=seed)
    solutions = solutions_frame()

    scaling.to_csv(tables / "gap1_remuscularisation_scaling.csv", index=False)
    pd.DataFrame([{**maturation, "status": SOLUTIONS_TAG}]).to_csv(
        tables / "gap1_maturation_benefit.csv", index=False
    )
    enrichment.to_csv(tables / "gap2_response_guided_enrichment.csv", index=False)
    options.to_csv(tables / "gap3_reporter_options.csv", index=False)
    cato.to_csv(tables / "gap4_cato_power.csv", index=False)
    tradeoff.to_csv(tables / "gap5_potency_redundancy_tradeoff.csv", index=False)
    outcomes.to_csv(tables / "gap6_outcomes_power.csv", index=False)
    solutions.to_csv(tables / "gap_solutions_register.csv", index=False)
    (reports / "gap_solutions_report.txt").write_text(
        solutions_report(scaling, maturation, cato, enrichment, outcomes, tradeoff)
    )
    (root / "DISCLAIMER.txt").write_text(DISCLAIMER + "\n")
    return {
        "output_dir": str(root),
        "outputs": 10,
        "status": SOLUTIONS_TAG,
        "gaps_addressed": len(GAP_SOLUTIONS),
        "claims_overturned": sum(1 for s in GAP_SOLUTIONS if s.was_overstated),
        "still_open": sum(1 for s in GAP_SOLUTIONS if s.status is SolutionStatus.STILL_OPEN),
        "maturation_multiplier": maturation["maturation_multiplier"],
        "mass_multiplier_after_maturation": maturation["mass_multiplier_at_17_months"],
        "cato_power_for_3pp_at_sd6": float(
            cato.loc[cato["assumed_sd_pp"].eq(6.0)].iloc[0]["power_for_3pp"]
        ),
        "outcomes_events_full_slope": int(
            outcomes.loc[outcomes["assumption"].eq("ctt_slope_full")].iloc[0]["events_required"]
        ),
        "cures_declared": 0,
        "disclaimer": DISCLAIMER,
    }
