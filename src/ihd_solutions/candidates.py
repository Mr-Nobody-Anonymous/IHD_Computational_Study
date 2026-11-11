"""
Which IHD candidate is best supported - the literature's answer, ranked.

THE QUESTION
------------
"What is the best candidate for IHD, from the internet and from our results?"

Answering it requires a criterion, and the transfer gate already says what the criterion
must be: an **outcome**, not a surrogate (gate item 1), measured in a **randomised
comparison** (gate item 6). So candidates are ranked by the strength of their OUTCOME
evidence, and everything else is reported as what it is - a mechanism, a hypothesis, or a
marker.

This module does not pick a favourite. It orders candidates by a stated rule, records what
the rule is, and states what the ranking does NOT mean.

THE RULE
--------
1. Randomised outcome evidence, positive  >  2. randomised outcome evidence, neutral/refuted
>  3. randomised surrogate evidence only  >  4. uncontrolled surrogate evidence  >
5. mechanism or magnitude claim with no clinical data.

Within a tier, candidates are ordered by effect magnitude on the endpoint that tier uses.

THE RESULT, AND IT IS NOT FLATTERING TO THIS REPOSITORY
------------------------------------------------------
The best-supported candidate for IHD is **not a new architecture**. It is deep LDL-C
lowering delivered by an agent that already exists, because that is the only mechanism with
repeated, randomised, mortality-level evidence.

And the repository's own central model has been overtaken by a trial. ZEUS randomised
6,385 patients to an IL-6 ligand inhibitor and found **HR 0.99 (0.88-1.11)** for
cardiovascular death, myocardial infarction or stroke, despite hsCRP falling 64% and IL-6
falling 76%. The repository's inflammatory cascade exists to predict exactly that
intervention's effect on exactly that marker. `crp_is_no_longer_a_surrogate()` records the
consequence: the calibration target of `calibration.py` is a marker whose pharmacological
manipulation has now been shown not to move outcomes.

No model is changed here. Nothing is tuned. No candidate is declared a cure.

SOURCES
-------
Every entry carries its trial identifier and the reported effect. Where a result is a press
release rather than a peer-reviewed paper, the entry says so.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

# =============================================================================
# Constants. Registered in provenance.REGISTRY; a missing record fails the build.
# =============================================================================

#: Relative risk reduction in major vascular events per 1 mmol/L (38.7 mg/dL) LDL-C
#: reduction, from the Cholesterol Treatment Trialists' meta-analysis. This is the single
#: best-established quantitative link between a modifiable quantity and IHD outcomes.
CTT_RRR_PER_MMOL = 0.22  # measured

#: Surrogate threshold for LVEF: the treatment effect on LVEF needed to imply a mortality
#: benefit, from the heart-failure surrogate meta-analysis.
LVEF_SURROGATE_THRESHOLD_PP = 6.0  # measured

#: ZEUS hazard ratio point estimate for 3-point MACE. A value at unity means the effect is
#: not distinguishable from none.
ZEUS_HR = 0.99  # measured
ZEUS_CI_LOW = 0.88  # measured
ZEUS_CI_HIGH = 1.11  # measured
ZEUS_N = 6385  # measured
ZEUS_CRP_REDUCTION = 0.64  # measured
ZEUS_IL6_REDUCTION = 0.76  # measured
#: Serious-infection rates in ZEUS, drug arm then comparator. Two scalars rather than one
#: tuple, because the provenance audit scans numeric constants and a tuple would slip past it.
ZEUS_SERIOUS_INFECTION = 0.20  # measured
ZEUS_SERIOUS_INFECTION_CONTROL = 0.18  # measured

#: Enlicitide placebo-adjusted LDL-C reduction at week 24 (CORALreef Lipids).
ENLICITIDE_LDL_REDUCTION = 0.56  # measured

#: VERVE-102 mean LDL-C reduction at 1.0 mg/kg (Heart-2 interim, n=35, open-label).
VERVE102_LDL_REDUCTION = 0.62  # measured
VERVE102_N = 35  # measured

#: CORALreef Outcomes enrols this many; used to state how quickly an answer could arrive.
CORALREEF_OUTCOMES_N = 14500  # measured

#: A candidate is called "supported by outcome evidence" only at or above this n, because
#: below it a null and a benefit are not separable at any plausible effect size.
MIN_OUTCOME_TRIAL_N = 1000  # estimated


@dataclass(frozen=True)
class Candidate:
    """One intervention, its best available evidence, and what that evidence licenses."""

    key: str
    name: str
    mechanism: str
    best_trial: str
    n: int
    endpoint: str
    endpoint_class: str
    result: str
    effect: float | None
    #: Unit of ``effect``. Two candidates with different units are not comparable by the
    #: number; ranking therefore stops at the tier and never averages across units.
    effect_unit: str
    evidence_tier: str
    #: Which way ``effect`` points. Ranking must never invert a refutation into a win:
    #: ZEUS is a *measured* 0.99 hazard ratio, not a 1% benefit.
    direction: str
    what_it_licenses: str
    what_it_does_not: str
    source: str

    def __post_init__(self) -> None:
        if self.evidence_tier not in TIERS:
            raise ValueError(f"unknown evidence tier {self.evidence_tier!r}")
        if self.endpoint_class not in ENDPOINT_CLASSES:
            raise ValueError(f"unknown endpoint class {self.endpoint_class!r}")
        if not self.source:
            raise ValueError("a candidate must cite a source")
        if not self.what_it_does_not:
            raise ValueError("a candidate must state what its evidence does not license")
        if self.direction not in DIRECTIONS:
            raise ValueError(f"unknown direction {self.direction!r}")
        if self.effect is not None and not math.isfinite(self.effect):
            raise ValueError("use None, not infinity, when no effect is estimable")


#: Ordered from strongest to weakest. The rank IS the claim.
#:
#: The key doubles as the rank key, so this dict is the single source of truth: a
#: candidate cannot claim a tier that is not defined here, and reordering the tiers
#: reorders the table without touching any entry.
TIERS = {
    "tier_1_positive_outcome": "randomised outcome evidence, positive",
    "tier_2_refuted": "randomised outcome evidence, neutral or refuted",
    "tier_3_randomised_none_surrogate": "randomised evidence on a surrogate only",
    "tier_4_uncontrolled_surrogate": "uncontrolled evidence on a surrogate",
    "tier_5_no_clinical_data": "mechanism or magnitude claim with no clinical outcome data",
}

#: Rank order, derived from TIERS so the two can never disagree.
TIER_ORDER = {tier: rank for rank, tier in enumerate(TIERS, start=1)}

#: Which way a reported effect points. A refuted mechanism reports a real number that
#: argues against itself, so the number alone cannot be used to rank.
DIRECTIONS = {
    "benefit": "the reported effect favours the intervention",
    "harm": "the reported effect favours the comparator",
    "neutral": "the reported effect is an estimate at or straddling unity",
    "none": "no effect was estimable because no controlled comparison exists",
}

ENDPOINT_CLASSES = {
    "outcome": "death, myocardial infarction, stroke or heart-failure hospitalisation",
    "surrogate": "a measured physiological or biochemical quantity, not events",
}


CANDIDATES: tuple[Candidate, ...] = (
    Candidate(
        key="ldl_lowering_approved",
        name="Deep LDL-C lowering to guideline target with approved agents",
        mechanism="hepatic LDL-receptor availability; the causal quantity is cumulative "
        "LDL-C exposure",
        best_trial="CTT meta-analysis + FOURIER/ODYSSEY; oral PCSK9 inhibitor now approved",
        n=170_000,
        endpoint="major vascular events",
        endpoint_class="outcome",
        result="22% relative risk reduction per 1 mmol/L LDL-C reduction, reproducible "
        "across statins, ezetimibe and PCSK9 inhibitors",
        effect=CTT_RRR_PER_MMOL,
        effect_unit="relative reduction in major vascular events per 1 mmol/L LDL-C",
        evidence_tier="tier_1_positive_outcome",
        direction="benefit",
        what_it_licenses=(
            "the strongest quantitative claim available in IHD: lowering LDL-C by a given "
            "amount lowers event rate by a predictable proportion, and the effect is "
            "additive with baseline risk"
        ),
        what_it_does_not=(
            "it does not reverse established scar or restore myocardium; it slows the "
            "process. It also does not remove residual risk - roughly two thirds of events "
            "in treated patients are not prevented"
        ),
        source="CTT Collaboration; NEJM 2026 CORALreef Lipids (doi:10.1056/NEJMoa2511002)",
    ),
    Candidate(
        key="colchicine",
        name="Low-dose colchicine 0.5 mg daily",
        mechanism="tubulin binding; broadly anti-inflammatory, upstream of the "
        "inflammasome, NOT an IL-6 pathway drug",
        best_trial="COLCOT, LoDoCo2, CLEAR SYNERGY; Cochrane review 2025",
        n=21_425,
        endpoint="major adverse cardiovascular events",
        endpoint_class="outcome",
        result="about 9 fewer myocardial infarctions and 8 fewer strokes per 1,000 treated; "
        "no effect on cardiovascular mortality; CLEAR SYNERGY neutral in acute coronary "
        "syndrome",
        effect=0.017,
        effect_unit="absolute risk difference at 1 year",
        direction="benefit",
        evidence_tier="tier_1_positive_outcome",
        what_it_licenses=(
            "guideline backing: ESC 2024 chronic coronary syndrome Class IIa Level A; "
            "ACC/AHA Class IIb. The only anti-inflammatory with a positive outcome signal "
            "and a guideline place"
        ),
        what_it_does_not=(
            "it does not reduce cardiovascular mortality, it is a Class IIa or IIb "
            "recommendation rather than Class I, the recent large ACS trial was neutral, "
            "and gastrointestinal adverse events increase"
        ),
        source="ESC 2024 CCS guidelines; Cochrane 2025; CLEAR SYNERGY; Rossello 2025",
    ),
    Candidate(
        key="il6_inhibition",
        name="IL-6 ligand inhibition (ziltivekimab)",
        mechanism="IL-6 to C-reactive protein; the exact axis this repository's "
        "inflammatory model simulates",
        best_trial="ZEUS (NCT05021835), plus HERMES and ATHENA stopped for futility",
        n=ZEUS_N,
        endpoint="cardiovascular death, non-fatal myocardial infarction, non-fatal stroke",
        endpoint_class="outcome",
        result=f"NEUTRAL: HR {ZEUS_HR:.2f} (95% CI {ZEUS_CI_LOW:.2f}-{ZEUS_CI_HIGH:.2f}) "
        f"despite hsCRP falling {ZEUS_CRP_REDUCTION:.0%} and IL-6 falling "
        f"{ZEUS_IL6_REDUCTION:.0%}; serious infections "
        f"{ZEUS_SERIOUS_INFECTION:.0%} vs {ZEUS_SERIOUS_INFECTION_CONTROL:.0%}",
        effect=ZEUS_HR,
        effect_unit="hazard ratio",
        direction="neutral",
        evidence_tier="tier_2_refuted",
        what_it_licenses=(
            "a refutation, in 6,385 patients, of the proposition that suppressing IL-6 and "
            "CRP in high-risk ASCVD reduces events. Target engagement was demonstrated and "
            "benefit was not"
        ),
        what_it_does_not=(
            "it does not prove inflammation is irrelevant - the population had advanced "
            "disease and the intervention was downstream. It does show that lowering the "
            "marker is not sufficient"
        ),
        source="ZEUS, reported 31 July 2026, published NEJM; HERMES/ATHENA termination "
        "September 2026",
    ),
    Candidate(
        key="cell_patch",
        name="iPSC-derived cardiomyocyte patch (BioVAT / AATS)",
        mechanism="mechanical remuscularisation and possible wall-stress reduction",
        best_trial="BioVAT-HF (NCT04396899); AATS late-breaking LB8",
        n=20,
        endpoint="left ventricular ejection fraction",
        endpoint_class="surrogate",
        result=f"LVEF +3.9 pp (90% CI 0.9-6.8) uncontrolled; second trial n=8 +4.3 pp, "
        f"p=0.156. Both below the {LVEF_SURROGATE_THRESHOLD_PP:.0f} pp surrogate threshold",
        effect=3.9,
        effect_unit="percentage points of LVEF",
        evidence_tier="tier_4_uncontrolled_surrogate",
        direction="benefit",
        what_it_licenses=(
            "a signal worth testing in a randomised trial with an outcome endpoint, and "
            "nothing beyond that"
        ),
        what_it_does_not=(
            "it cannot separate a real effect from natural history or regression to the "
            "mean: 20 patients, no control, and an effect below the threshold at which LVEF "
            "change predicts mortality benefit"
        ),
        source="NEJM 2026 (BioVAT-HF); AATS 2026 late-breaking; REVIVED BCIS2",
    ),
    Candidate(
        key="pcsk9_base_editing",
        name="VERVE-102 in vivo PCSK9 base editing",
        mechanism="permanent PCSK9 inactivation in hepatocytes, mimicking loss-of-function "
        "variants",
        best_trial="Heart-2 (NCT06164730), phase 1b, open-label, single ascending dose",
        n=VERVE102_N,
        endpoint="LDL-C",
        endpoint_class="surrogate",
        result=f"LDL-C -{VERVE102_LDL_REDUCTION:.0%} at 1.0 mg/kg (-78 mg/dL), PCSK9 -88%, "
        "durable to 18 months, no treatment-related serious adverse events",
        effect=VERVE102_LDL_REDUCTION,
        effect_unit="relative LDL-C reduction",
        evidence_tier="tier_3_randomised_none_surrogate",
        direction="benefit",
        what_it_licenses=(
            "the largest durable LDL-C reduction from a single intervention, on a mechanism "
            "that already has outcome support through other agents"
        ),
        what_it_does_not=(
            "its own outcome benefit is unproven and the edit is irreversible. An approved "
            "oral PCSK9 inhibitor now produces a similar LDL-C reduction without a "
            "permanent edit, which removes most of the practical case for the editing route"
        ),
        source="NEJM 2026 (Heart-2, n=35, open-label); Lilly press release 25 May 2026",
    ),
    Candidate(
        key="oral_pcsk9",
        name="Oral PCSK9 inhibitor (enlicitide, Lipfendra)",
        mechanism="macrocyclic peptide binding circulating PCSK9, preserving LDL-receptor "
        "recycling",
        best_trial="CORALreef Lipids and HeFH; CORALreef Outcomes running",
        n=2_909,
        endpoint="LDL-C",
        endpoint_class="surrogate",
        result=f"LDL-C -{ENLICITIDE_LDL_REDUCTION:.0%} vs placebo at week 24, sustained to "
        "week 52; FDA approved 16 July 2026",
        effect=ENLICITIDE_LDL_REDUCTION,
        effect_unit="relative LDL-C reduction",
        evidence_tier="tier_3_randomised_none_surrogate",
        direction="benefit",
        what_it_licenses=(
            "attains injectable-antibody LDL-C lowering as a once-daily tablet, which is "
            "the adherence problem solved without an irreversible edit"
        ),
        what_it_does_not=(
            "the outcome trial has completed enrolment and not reported; approval rests on "
            "the surrogate, exactly as the transfer gate warns"
        ),
        source="NEJM 2026 (CORALreef Lipids); FDA approval 16 July 2026; CORALreef Outcomes "
        "NCT06008756",
    ),
    Candidate(
        key="atlas5_combination",
        name="Combined architecture (cell patch + base editing + OMT)",
        mechanism="additive or synergistic combination of the above",
        best_trial="none - never administered together in a controlled setting",
        n=0,
        endpoint="none",
        endpoint_class="surrogate",
        result="no clinical data; corrected 20,000-draw estimates are 57.70% for ATLAS-5 "
        "and 57.81% for the re-optimized simulation-only ATLAS-5F; neither is robust",
        effect=None,
        effect_unit="not estimable - no controlled comparison",
        evidence_tier="tier_5_no_clinical_data",
        direction="none",
        what_it_licenses="a hypothesis, and an ordered experiment list",
        what_it_does_not=(
            "it cannot claim the sum of its parts. The repository measured 15 pairwise "
            "interactions and refused 14.08 pp of double-counted benefit, and the wider "
            "re-search found no robust composition"
        ),
        source="this repository: results/revival/, results/recombine/",
    ),
)


#: Direction, ordered so that a refutation never sorts above a benefit at the same tier.
DIRECTION_ORDER = {"benefit": 0, "neutral": 1, "harm": 2, "none": 3}


def _format_effect(effect: float | None, unit: str) -> str:
    """Render an effect with its unit, so two numbers are never shown as comparable.

    The trap this avoids: a 3.9 percentage-point LVEF gain and a 0.22 relative risk
    reduction are both "effects", and printing them side by side invites a comparison that
    means nothing.
    """
    if effect is None:
        return "not estimable - no controlled comparison"
    if unit.startswith("hazard"):
        return f"HR {effect:.2f}"
    if unit.startswith("percentage points"):
        return f"{effect:+.1f} {unit}"
    return f"{effect:.1%} {unit}"


def candidate_table() -> pd.DataFrame:
    """Every candidate, ordered by the stated evidence rule."""
    rows = []
    for candidate in CANDIDATES:
        rows.append(
            {
                "rank_key": TIER_ORDER[candidate.evidence_tier],
                "candidate": candidate.name,
                "key": candidate.key,
                "tier": TIERS[candidate.evidence_tier],
                "endpoint_class": candidate.endpoint_class,
                "best_trial": candidate.best_trial,
                "n": candidate.n,
                "result": candidate.result,
                "effect": candidate.effect,
                "effect_unit": candidate.effect_unit,
                "direction": candidate.direction,
                "direction_rank": DIRECTION_ORDER[candidate.direction],
                "effect_display": _format_effect(candidate.effect, candidate.effect_unit),
                "licenses": candidate.what_it_licenses,
                "does_not_license": candidate.what_it_does_not,
                "source": candidate.source,
            }
        )
    # Ordering never mixes units: tier, then whether the reported direction favours the
    # intervention, then trial size, then the magnitude - the last only ever compared
    # between rows that share a unit, which the sort's stability does not disturb because
    # same-tier rows in this register are always same-endpoint.
    frame = pd.DataFrame(rows).sort_values(
        ["rank_key", "direction_rank", "n", "effect"],
        ascending=[True, True, False, False],
        kind="stable",
    )
    frame.insert(0, "rank", range(1, len(frame) + 1))
    return frame.reset_index(drop=True)


def best_candidate() -> dict[str, object]:
    """The answer, with the reasoning that produced it and its limits.

    The ranking's top entry rests on the Cholesterol Treatment Trialists' meta-analysis,
    which is the most replicated quantitative relationship in cardiovascular prevention.
    The recommendation is deliberately unglamorous: it is the thing already known to work.
    """
    table = candidate_table()
    top = table.iloc[0]
    # "Backed by outcome evidence" is stricter than "measured on an outcome": a refuted
    # outcome trial is outcome evidence against a candidate, not for it, and a trial too
    # small to separate a null from a modest benefit does not back anything.
    outcome_backed = table.loc[
        (table["endpoint_class"] == "outcome")
        & (table["direction"] == "benefit")
        & (table["n"] >= MIN_OUTCOME_TRIAL_N)
    ]
    refuted = table.loc[(table["endpoint_class"] == "outcome") & (table["direction"] != "benefit")]
    return {
        "best_candidate": str(top["candidate"]),
        "mechanism": "cumulative LDL-C exposure",
        "why": (
            "it is the only mechanism with repeated randomised evidence at the mortality "
            "level, and the effect is quantitative and predictable rather than directional"
        ),
        "quantitative_effect": f"{CTT_RRR_PER_MMOL:.0%} relative reduction in major "
        "vascular events per 1 mmol/L LDL-C reduction",
        "second_choice": (
            "low-dose colchicine - the only anti-inflammatory with a positive outcome "
            "signal and a guideline place, and notably NOT an IL-6 inhibitor"
        ),
        "what_the_repository_contributes": (
            "it rules things out: the combined architecture is not robust in its own model, "
            "the cell patch is below the LVEF surrogate threshold, and the inflammatory "
            "arm's marker has now been shown in 6,385 patients not to move outcomes when "
            "manipulated"
        ),
        "candidates_with_outcome_evidence": int(len(outcome_backed)),
        "candidates_refuted_by_outcome": int(len(refuted)),
        "candidates_total": len(table),
        "this_is_not_a_cure": (
            "no candidate here is a cure. None reverses established coronary disease or "
            "restores infarcted myocardium. The best-supported option slows progression and "
            "reduces events; it does not undo the disease"
        ),
        "status": "ranking_by_stated_evidence_rule",
    }


def inflammation_arm_refutation() -> dict[str, object]:
    """What ZEUS does to this repository's central model.

    The repository's inflammatory cascade runs oxLDL -> NF-kB -> NLRP3 -> IL-1b -> IL-6 ->
    CRP. ZEUS intervened at the last step with a targeted antibody in 6,385 patients: hsCRP
    fell 64% and IL-6 fell 76%, and the hazard ratio for cardiovascular death, myocardial
    infarction or stroke was 0.99.

    That is not a model error. It is a refutation of an inference the repository was
    entitled to make only conditionally: **that moving the cascade moves the outcome.**
    """
    return {
        "trial": "ZEUS, NCT05021835",
        "n": ZEUS_N,
        "hazard_ratio": ZEUS_HR,
        "ci": (ZEUS_CI_LOW, ZEUS_CI_HIGH),
        "target_engagement": {
            "hsCRP_reduction": ZEUS_CRP_REDUCTION,
            "IL6_reduction": ZEUS_IL6_REDUCTION,
        },
        "intervened_step_in_our_model": "k_il6_prod downward",
        "our_model_would_predict": (
            "a large fall in CRP, which matches what ZEUS observed - our model gets the "
            "marker dynamics right and has nothing to say about events"
        ),
        "the_refuted_inference": (
            "that a mechanism which moves a marker therefore moves the outcome"
        ),
        "action_for_this_repository": (
            "demote the inflammatory model from 'mechanism of benefit' to 'marker "
            "dynamics', and stop using CRP as a proxy for benefit anywhere"
        ),
        "what_would_change_it_back": (
            "a positive outcome trial of an upstream inflammasome intervention in a less "
            "advanced population; that trial has not reported"
        ),
        "safety_note": (
            "serious infections were more frequent with IL-6 inhibition (20% vs 18%), so "
            "the intervention was not even free"
        ),
    }


def crp_is_no_longer_a_surrogate() -> dict[str, object]:
    """The repository's calibration target has failed as a surrogate, by trial.

    `calibration.py` fits the legacy cascade to terminal CRP and reports which rates could
    be recovered from it. That work remains valid as a demonstration of non-identifiability.

    But CRP can no longer be used as a *target for benefit*, and this is now an empirical
    finding rather than a caveat: manipulating CRP pharmacologically at scale did not move
    outcomes. The relevant lesson is not "do not measure CRP" - it remains a strong risk
    marker - but "a risk marker is not a treatment target", which is exactly the
    distinction the transfer gate item 6 exists to enforce.
    """
    return {
        "previous_use": "calibration target for the legacy 14-rate cascade",
        "still_valid_as": "a demonstration of structural non-identifiability",
        "no_longer_valid_as": "a proxy for treatment benefit",
        "evidence": f"ZEUS, n={ZEUS_N}, HR {ZEUS_HR:.2f} "
        f"({ZEUS_CI_LOW:.2f}-{ZEUS_CI_HIGH:.2f}), hsCRP -{ZEUS_CRP_REDUCTION:.0%}",
        "the_distinction_that_matters": (
            "a marker can be strongly prognostic and still be causally irrelevant to "
            "treatment benefit; ZEUS separates those two roles for CRP"
        ),
        "consequence_for_the_audit": (
            "the calibration work stands, but any result phrased as 'this reduces CRP "
            "therefore it should help' is now unsupported and must be removed"
        ),
    }


DISCLAIMER = """\
WHAT THIS MODULE IS
  A ranking of IHD candidates by a stated evidence rule: randomised outcome evidence
  first, surrogates after, mechanisms last. Every entry cites its trial and its effect.

WHAT IT IS NOT
  A treatment recommendation for any individual, and not a cure for anything. The top
  entry slows a process; it does not reverse disease.

THE HEADLINE
  The best-supported candidate is deep LDL-C lowering with approved agents - the
  unglamorous answer. The repository's own novel architecture ranks last, because it has
  no clinical data, and its inflammatory arm has been refuted at scale by ZEUS.
"""


def candidate_report() -> str:
    table = candidate_table()
    best = best_candidate()
    refutation = inflammation_arm_refutation()
    surrogate = crp_is_no_longer_a_surrogate()
    lines = [
        "=" * 78,
        "BEST IHD CANDIDATE - RANKED BY OUTCOME EVIDENCE",
        "=" * 78,
        "",
        DISCLAIMER.rstrip(),
        "",
        f"  ranked candidates            {len(table)}",
        f"  backed by outcome evidence   {best['candidates_with_outcome_evidence']}",
        f"  refuted by outcome evidence  {best['candidates_refuted_by_outcome']}",
        f"  best candidate               {best['best_candidate']}",
        f"  second                       {best['second_choice']}",
        f"  quantitative effect          {best['quantitative_effect']}",
        "",
        "-" * 78,
        "THE RULE - STATED BEFORE THE ANSWER",
        "-" * 78,
        "  Candidates are ranked by the strength of OUTCOME evidence, because that is the",
        "  criterion the transfer gate already sets: item 1 demands an outcome rather than",
        "  a surrogate, item 6 a randomised comparison. Within a tier, a refutation is",
        "  never counted as a win, and magnitudes are compared only inside a shared unit.",
        "",
    ]
    for rank, key in enumerate(TIERS, start=1):
        count = int((table["rank_key"] == TIER_ORDER[key]).sum())
        lines.append(f"  {rank}. {TIERS[key]}  [{count}]")
    lines += [
        "",
        "-" * 78,
        "THE RANKING",
        "-" * 78,
    ]
    for row in table.itertuples(index=False):
        lines += [
            f"  {row.rank}. {row.candidate}",
            f"     tier     : {row.tier}",
            f"     trial    : {row.best_trial}  (n={row.n:,})",
            f"     endpoint : {row.endpoint_class}",
            f"     effect   : {row.effect_display}  ({row.direction})",
            f"     result   : {row.result}",
            f"     licenses : {row.licenses}",
            f"     does not : {row.does_not_license}",
            "",
        ]
    lines += [
        "-" * 78,
        "THE REPOSITORY'S INFLAMMATORY ARM HAS BEEN REFUTED AT SCALE",
        "-" * 78,
        f"  trial                      {refutation['trial']}",
        f"  n                          {refutation['n']:,}",
        f"  hazard ratio               {refutation['hazard_ratio']:.2f} "
        f"({refutation['ci'][0]:.2f}-{refutation['ci'][1]:.2f})",
        f"  hsCRP reduction            {refutation['target_engagement']['hsCRP_reduction']:.0%}",
        f"  IL-6 reduction             {refutation['target_engagement']['IL6_reduction']:.0%}",
        "",
        "  Our model would have predicted the marker fall correctly. It has nothing to",
        "  say about events, and the trial now says the marker fall was not enough.",
        f"  Refuted inference: {refutation['the_refuted_inference']}",
        f"  Action: {refutation['action_for_this_repository']}",
        f"  Safety: {refutation['safety_note']}",
        "",
        "-" * 78,
        "CONSEQUENCE FOR THE CALIBRATION WORK",
        "-" * 78,
        f"  still valid as        {surrogate['still_valid_as']}",
        f"  no longer valid as    {surrogate['no_longer_valid_as']}",
        f"  evidence              {surrogate['evidence']}",
        f"  {surrogate['the_distinction_that_matters']}",
        "",
        "-" * 78,
        "WHAT THIS IS NOT",
        "-" * 78,
        f"  {best['this_is_not_a_cure']}",
        "",
    ]
    return "\n".join(lines)


def write_candidate_report(output_dir: str = "results/candidates") -> dict[str, object]:
    """Write the ranking and its supporting tables."""
    root = Path(output_dir)
    tables = root / "tables"
    tables.mkdir(parents=True, exist_ok=True)

    candidate_table().to_csv(tables / "candidates_ranking.csv", index=False)
    pd.DataFrame([inflammation_arm_refutation()]).to_csv(
        tables / "candidates_inflammation_refutation.csv", index=False
    )
    pd.DataFrame([crp_is_no_longer_a_surrogate()]).to_csv(
        tables / "candidates_crp_surrogate_status.csv", index=False
    )
    pd.DataFrame([best_candidate()]).to_csv(tables / "candidates_best.csv", index=False)
    (root / "DISCLAIMER.txt").write_text(DISCLAIMER)
    (root / "candidate_report.txt").write_text(candidate_report())
    table = candidate_table()
    best = best_candidate()
    return {
        "artefacts": 6,
        "candidates_ranked": len(table),
        "with_outcome_evidence": best["candidates_with_outcome_evidence"],
        "best_candidate": best["best_candidate"],
        "repository_architecture_rank": int(
            table.loc[table["key"] == "atlas5_combination", "rank"].iloc[0]
        ),
        "inflammation_arm_status": "refuted by ZEUS",
        "output_dir": str(root),
        "status": "candidate_ranking_no_model_change",
    }
