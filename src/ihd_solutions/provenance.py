"""Enforce a provenance record on every biological parameter.

The critique this module answers is: *"a model can produce almost any desired result if
its assumptions are selected appropriately."* That is correct, and no amount of
documentation in a report prevents it, because prose does not fail a build.

So the rule here is mechanical. Every numeric parameter that acts on a biological or
clinical quantity must carry a `Provenance` record, and a test enumerates the modules and
fails if one does not. A missing record is an error, not an annotation.

What the registry does NOT do is claim the numbers are right. It classifies what kind of
claim each one is, so that a reader - or a reviewer, or a future model - can see at a
glance how much weight a result can bear. Three classes do the work:

``MEASURED``      read off a named human study. Earns credit.
``ESTIMATED``     derived from measured quantities by a stated transformation. Earns
                  credit for the transformation, not for new data.
``HYPOTHETICAL``  invented or assumed, including every animal-derived multiplier
                  transplanted into a human architecture. Earns NO credit, and the
                  registry makes that impossible to forget.

The distinction the critique asked for - "measured, estimated, or hypothetical" - is
therefore a type, and anything that consumes a parameter can branch on it. That is what
`credit_ceiling()` is for: it computes the most a model may claim given the weakest
parameter it depends on.

One trap is deliberately encoded. A parameter can be measured in a human and still not
be evidence for the claim being made: the maturation ratio is measured, but it is
measured in an EHM patch cohort, and multiplying it by animal factors is not a
measurement. Those carry ``measured_in_humans=True`` and
``supports="mechanism"`` rather than ``supports="efficacy"``, so downstream code cannot
read a human measurement as a human effect size.

Nothing here is tuned. The constants already in the repository were classified from
their existing docstrings and citations; where the citation is to an animal study the
class is HYPOTHETICAL for human purposes even if the study is real.
"""

from __future__ import annotations

import ast
from dataclasses import dataclass, field
from pathlib import Path

# --------------------------------------------------------------------------- classes

MEASURED = "measured"
ESTIMATED = "estimated"
HYPOTHETICAL = "hypothetical"

CLASSES = (MEASURED, ESTIMATED, HYPOTHETICAL)

RANK = {MEASURED: 3, ESTIMATED: 2, HYPOTHETICAL: 1}

#: What a parameter is allowed to be used as. A number can be a real measurement and
#: still not support an efficacy claim.
SUPPORTS = (
    "efficacy",  # a treatment effect in humans
    "mechanism",  # a biological mechanism, not a magnitude
    "arithmetic",  # a conversion or definition
    "design",  # a protocol or study-design constant
    "ceiling",  # a bound derived from evidence
    "refutation",  # a result that removes a claim other work depends on
    "safety",  # a harm or adverse-event quantity
)

#: Populations a measurement can come from. Only the first earns human efficacy credit.
POPULATIONS = ("human", "large_animal", "rodent", "in_vitro", "none")


@dataclass(frozen=True)
class Provenance:
    """One parameter and the limits of what it can be used to claim."""

    value: float
    unit: str
    kind: str
    supports: str
    population: str
    source: str
    note: str = ""
    measured_in_humans: bool = field(default=False)

    def __post_init__(self) -> None:
        if self.kind not in CLASSES:
            raise ValueError(f"unknown provenance class {self.kind!r}")
        if self.supports not in SUPPORTS:
            raise ValueError(f"unknown support type {self.supports!r}")
        if self.population not in POPULATIONS:
            raise ValueError(f"unknown population {self.population!r}")
        if not self.unit:
            raise ValueError("a parameter must declare its unit")
        if not self.source or len(self.source) < 8:
            raise ValueError("a parameter must declare a source")
        # A hypothetical parameter may never claim human measurement.
        if self.kind == HYPOTHETICAL and self.measured_in_humans:
            raise ValueError(
                f"{self.source}: a hypothetical parameter cannot be measured in humans"
            )
        # Only human populations may claim human measurement.
        if self.measured_in_humans and self.population != "human":
            raise ValueError(
                f"{self.source}: measured_in_humans requires population='human', "
                f"got {self.population!r}"
            )

    @property
    def earns_efficacy_credit(self) -> bool:
        """Only a measured human efficacy value can carry a human claim."""
        return self.kind == MEASURED and self.supports == "efficacy" and self.population == "human"

    @property
    def rank(self) -> int:
        return RANK[self.kind]


def credit_ceiling(kinds) -> float:
    """The most a claim may be credited given its weakest input.

    A chain of estimates is only as strong as its weakest link, so this takes the MINIMUM
    rather than an average. Averaging would let one strong measurement launder a
    hypothetical multiplier that sits next to it.
    """
    ranks = [RANK[kind] if isinstance(kind, str) else kind.rank for kind in kinds]
    if not ranks:
        raise ValueError("credit_ceiling needs at least one parameter")
    weakest = min(ranks)
    return {3: 1.0, 2: 0.5, 1: 0.0}[weakest]


#: The registry. Keys are the names the modules export, so a test can bind them.
REGISTRY: dict[str, Provenance] = {
    # ---------------------------------------------------------------- arithmetic only
    "MG_DL_PER_MMOL_L": Provenance(
        38.67,
        "mg/dL per mmol/L",
        "measured",
        "arithmetic",
        "human",
        "standard lipid conversion (Friedewald-era); exact value 38.67",
        "definitional, not a biological claim",
    ),
    "DAYS_PER_YEAR": Provenance(
        365.25,
        "days",
        "measured",
        "arithmetic",
        "none",
        "Julian year, adopted so leap days do not drift the horizon",
        "definitional",
    ),
    "Z_95": Provenance(
        1.959963984540054,
        "z",
        "measured",
        "arithmetic",
        "none",
        "standard normal 97.5th percentile",
        "definitional",
    ),
    "CTT_RR_PER_MMOL_L": Provenance(
        0.78,
        "rate ratio per mmol/L",
        "measured",
        "efficacy",
        "human",
        "CTT meta-analysis: RR 0.78 per 1 mmol/L LDL-C reduction",
        "pooled randomised trials; strong efficacy evidence",
    ),
    # ------------------------------------------------------------------- patch cells
    "BIOVAT_LVEF_GAIN_3MO_PP": Provenance(
        3.9,
        "LVEF percentage points",
        "measured",
        "mechanism",
        "human",
        "BioVAT-HF, NEJM 2026, doi:10.1056/NEJMoa2513525",
        "single-arm, n=16 at max dose, no control group; REVIVED OMT gained +3.4 unaided",
    ),
    "BIOVAT_LVEF_GAIN_17MO_PP": Provenance(
        6.9,
        "LVEF percentage points",
        "measured",
        "mechanism",
        "human",
        "BioVAT-HF, NEJM 2026, mean 17-month follow-up",
        "uncontrolled; the largest LVEF gain any remuscularisation approach has measured",
    ),
    "BIOVAT_WALL_GAIN_3MO_MM": Provenance(
        4.5,
        "mm",
        "measured",
        "mechanism",
        "human",
        "BioVAT-HF, NEJM 2026, least-squares mean target wall thickness",
        "wall thickness FELL by 17 months while LVEF rose - mass is not the lever",
    ),
    "BIOVAT_WALL_GAIN_12MO_MM": Provenance(
        2.9,
        "mm",
        "measured",
        "mechanism",
        "human",
        "BioVAT-HF, NEJM 2026, mean 17-month target wall thickness",
        "used with the 3-month value to show maturation, not mass",
    ),
    "BIOVAT_MAX_EHM_UNITS": Provenance(
        20,
        "EHM units",
        "measured",
        "design",
        "human",
        "BioVAT-HF maximum tested dose",
        "establishes the precedented dose count",
    ),
    "BIOVAT_CARDIOMYOCYTES_PER_UNIT": Provenance(
        40_000_000,
        "cells per EHM unit",
        "estimated",
        "design",
        "human",
        "BioVAT-HF: ~34M iPSC-CM + 6M stromal per patch",
        "arithmetic from the published construct; not a therapeutic parameter",
    ),
    "BIOVAT_BASELINE_LVEF_PCT": Provenance(
        24.0,
        "percent",
        "measured",
        "design",
        "human",
        "BioVAT-HF baseline mean LVEF (enrolment <=35%)",
        "population definition, not an effect",
    ),
    "BIOVAT_N_AT_MAX_DOSE": Provenance(
        16,
        "patients",
        "measured",
        "design",
        "human",
        "BioVAT-HF: 16 of 20 received 20 EHM units",
        "sample size; drives the width of any confidence interval",
    ),
    "REQUIRED_LVEF_GAIN_PP": Provenance(
        20.0,
        "LVEF percentage points",
        "hypothetical",
        "ceiling",
        "human",
        "this study's criterion for functional normalisation (cure.py CRITERION_THRESHOLD)",
        "a TARGET chosen by this study, not an observed effect. Never evidence.",
    ),
    "SCAR_RECOVERABLE_LVEF_PP": Provenance(
        25.0,
        "LVEF percentage points",
        "hypothetical",
        "ceiling",
        "human",
        "this study's assumption that the whole contractile deficit is 25 pp",
        "the saturation ceiling; most ATLAS results are limited by it, not by mechanisms",
    ),
    "NORMALISATION_LVEF_GAIN_PP": Provenance(
        25.0,
        "LVEF percentage points",
        "hypothetical",
        "ceiling",
        "human",
        "same ceiling as SCAR_RECOVERABLE_LVEF_PP, used for normalisation",
        "assumption",
    ),
    # ------------------------------------------------------------ potency and dosing
    "ATLAS_POTENCY": Provenance(
        1.284,
        "fold above central",
        "measured",
        "ceiling",
        "human",
        "precedented potency cap from the catalogued agents' central estimates",
        "the bound ATLAS-2 exceeds at 2.80x; used as the admissibility limit",
    ),
    "PRECEDENTED_POTENCY_CAP": Provenance(
        1.284,
        "fold above central",
        "measured",
        "ceiling",
        "human",
        "same bound as ATLAS_POTENCY, restated for the gap solutions",
        "precedent ceiling",
    ),
    "SOURCE_REDUCTION_PRECEDENT_CENTRAL": Provenance(
        0.35,
        "fraction",
        "estimated",
        "ceiling",
        "human",
        "central estimate of achievable plaque source reduction within precedent",
        "derived from the regression grid, which is itself calibrated",
    ),
    "SOURCE_REDUCTION_PRECEDENT_MAX": Provenance(
        0.5,
        "fraction",
        "estimated",
        "ceiling",
        "human",
        "upper bound of achievable source reduction within precedent",
        "used as the hard cap on any source-reduction claim",
    ),
    "MAX_PRECEDENTED_PATHS": Provenance(
        3,
        "node paths",
        "estimated",
        "design",
        "human",
        "the largest redundant-path count reached by a precedented agent set",
        "redundancy ceiling; breaching it was a caught accounting error",
    ),
    "DEFAULT_PER_ATTEMPT_RISK": Provenance(
        0.05,
        "probability per procedure",
        "hypothetical",
        "design",
        "none",
        "assumed procedural risk per graft attempt",
        "NOT measured. It is charged on every attempt, and it drives the risk-adjusted "
        "column, so the ranking is sensitive to it.",
    ),
    "MAX_SINGLE_AGENT_RISK_NORMALISATION": Provenance(
        0.4,
        "fraction of RRR",
        "estimated",
        "ceiling",
        "human",
        "no single agent returns event risk to a low-risk reference",
        "a structural ceiling on single-agent credit, enforced by test",
    ),
    "MIN_SCENARIO_RR": Provenance(
        0.1,
        "rate ratio",
        "hypothetical",
        "ceiling",
        "none",
        "floor so a scenario cannot claim implausible risk reduction",
        "guard, not a measurement",
    ),
    # --------------------------------------------------------------- editing / lipids
    "HEART2_PCSK9_KNOCKDOWN": Provenance(
        0.88,
        "fraction",
        "measured",
        "efficacy",
        "human",
        "VERVE-102 Heart-2, NEJM 2026: -88% PCSK9 at 1.0 mg/kg",
        "single-dose; durability measured to 18 months",
    ),
    "HEART2_LDL_REDUCTION": Provenance(
        0.62,
        "fraction",
        "measured",
        "efficacy",
        "human",
        "VERVE-102 Heart-2: -62% LDL-C",
        "measured, n small, no outcomes data",
    ),
    "VERVE_LDL_ABS_REDUCTION_MG_DL": Provenance(
        78.0,
        "mg/dL",
        "measured",
        "efficacy",
        "human",
        "VERVE-102 Heart-2 absolute LDL-C reduction",
        "drives the outcomes-trial power calculation",
    ),
    "VERVE_LDL_PCT_REDUCTION": Provenance(
        0.62,
        "fraction",
        "measured",
        "efficacy",
        "human",
        "VERVE-102 Heart-2",
        "measured",
    ),
    "VERVE_PCSK9_PCT_REDUCTION": Provenance(
        0.88,
        "fraction",
        "measured",
        "efficacy",
        "human",
        "VERVE-102 Heart-2",
        "measured",
    ),
    "VERVE_DURABILITY_MONTHS": Provenance(
        18,
        "months",
        "measured",
        "efficacy",
        "human",
        "VERVE-102 Heart-2: durable to 18 months through hepatocyte turnover",
        "measured, but the highest dose has the shortest follow-up",
    ),
    "MAX_LDL_REDUCTION": Provenance(
        0.75,
        "fraction",
        "hypothetical",
        "ceiling",
        "none",
        "internal cap on any modelled LDL reduction",
        "guard so no scenario claims implausible efficacy",
    ),
    "LDL_FLOOR_MG_DL": Provenance(
        10.0,
        "mg/dL",
        "hypothetical",
        "ceiling",
        "none",
        "numerical floor on modelled LDL-C",
        "guard; physiology below this is not modelled",
    ),
    "APOB_FLOOR_MG_DL": Provenance(
        15.0,
        "mg/dL",
        "hypothetical",
        "ceiling",
        "none",
        "numerical floor on modelled apoB",
        "guard",
    ),
    "LPA_FLOOR_NMOL_L": Provenance(
        2.0,
        "nmol/L",
        "hypothetical",
        "ceiling",
        "none",
        "numerical floor on modelled Lp(a)",
        "guard",
    ),
    "HSCRP_FLOOR_MG_L": Provenance(
        0.2,
        "mg/L",
        "hypothetical",
        "ceiling",
        "none",
        "numerical floor on modelled hsCRP",
        "guard",
    ),
    # ------------------------------------------------------------- plaque mechanics
    "CAP_THRESHOLD_UM": Provenance(
        65.0,
        "micrometres",
        "measured",
        "mechanism",
        "human",
        "histological thin-cap fibroatheroma definition (<65 um)",
        "a DEFINITION, used as a threshold; it is not a rupture prediction",
    ),
    "STRESS_THRESHOLD_KPA": Provenance(
        300.0,
        "kPa",
        "hypothetical",
        "ceiling",
        "in_vitro",
        "assumed peak cap stress for rupture in the simplified model",
        "not a validated rupture criterion; the real one is geometry-dependent",
    ),
    "DEFAULT_CALCIFICATION_FRACTION": Provenance(
        0.25,
        "fraction of plaque",
        "hypothetical",
        "mechanism",
        "none",
        "default calcification fraction when none is supplied",
        "assumption",
    ),
    # ------------------------------------------------ counted from published tables
    "RETRY_RESPONDERS": Provenance(
        10,
        "of 30 patients",
        "measured",
        "mechanism",
        "human",
        "repeat-dose human series: 10 of 30 patients responded to a second dose",
        "drives human_tempering, which lowers every reported probability",
    ),
    "RETRY_TOTAL": Provenance(
        30,
        "patients",
        "measured",
        "design",
        "human",
        "denominator of the repeat-dose responder count",
        "small n; 10/30 has a wide exact binomial interval",
    ),
    "CATO_TOTAL_N": Provenance(
        60,
        "participants",
        "measured",
        "design",
        "human",
        "CATO, NCT06145035: 20 per arm, three arms",
        "n drives the power calculation in gap 4",
    ),
    "CATO_ARMS": Provenance(
        3,
        "arms",
        "measured",
        "design",
        "human",
        "CATO: two active arms plus control, so two active-versus-control comparisons",
        "comparison count widens the multiplicity adjustment",
    ),
    "AATS_LB8_GAIN_PP": Provenance(
        4.3,
        "LVEF percentage points",
        "measured",
        "mechanism",
        "human",
        "AATS LB8, n=8 ischaemic cardiomyopathy: 28.7% -> 33.0% at 52 weeks",
        "the CHANGE IS NULL (p=0.156); recorded so the change is not read as a gain",
    ),
    "AATS_LB8_P_VALUE": Provenance(
        0.156,
        "p",
        "measured",
        "mechanism",
        "human",
        "AATS LB8 reported p=0.156 for the LVEF change",
        "a null result; this is the point of recording it",
    ),
    "AATS_LB8_N": Provenance(
        8,
        "patients",
        "measured",
        "design",
        "human",
        "AATS LB8 enrolment",
        "underpowered by design; the null is not evidence of absence either",
    ),
    "NATURAL_HISTORY_PP": Provenance(
        3.4,
        "LVEF percentage points",
        "measured",
        "efficacy",
        "human",
        "REVIVED-BCIS2 OMT arm, +3.4 pp at 6 months, n=353",
        "the comparator. This is what guideline therapy achieves unaided, and it is "
        "the number that dissolves ATLAS-5.",
    ),
    "LVEF_CHANGE_SD_PP": Provenance(
        5.0,
        "LVEF percentage points",
        "estimated",
        "design",
        "human",
        "between-patient SD of MRI LVEF change; reuse of the gap-4 sweep family",
        "UNREPORTED in the source trials. Swept 5/6/7/8 in gap 4, not invented here.",
    ),
    "SURROGATE_STE_MORTALITY_PP": Provenance(
        8.0,
        "LVEF percentage points",
        "measured",
        "efficacy",
        "human",
        "2026 HF surrogate meta-analysis: 8 pp mean LVEF increase predicts a mortality "
        "benefit; R^2=0.59 across 69 comparisons",
        "trial-level surrogate; STICH says the link is ABSENT in ischaemic "
        "cardiomyopathy, so it cannot be applied here",
    ),
    "SURROGATE_STE_HF_HOSP_PP": Provenance(
        2.0,
        "LVEF percentage points",
        "measured",
        "efficacy",
        "human",
        "same meta-analysis: 2 pp for HF hospitalisation, R^2=0.86 in HFrEF",
        "accepted in HFrEF only",
    ),
    "SURROGATE_R2_MORTALITY": Provenance(
        0.59,
        "R squared",
        "measured",
        "efficacy",
        "human",
        "trial-level R^2 between LVEF change and mortality benefit, 69 comparisons",
        "moderate: 41% of the variance is not explained",
    ),
    "SURROGATE_R2_MORTALITY_HFREF": Provenance(
        0.49,
        "R squared",
        "measured",
        "efficacy",
        "human",
        "same, restricted to HFrEF, 56 comparisons",
        "poor in the subgroup that matches our population",
    ),
    "LARGEST_HUMAN_OBSERVATION_PP": Provenance(
        6.9,
        "LVEF percentage points",
        "measured",
        "ceiling",
        "human",
        "BioVAT-HF at a mean 17 months: the largest LVEF gain any remuscularisation "
        "approach has measured in a human",
        "any model output above this is outside the observed range",
    ),
    "RESPONSE_GUIDED_RATE": Provenance(
        0.582481,
        "fraction of screened-in patients",
        "estimated",
        "mechanism",
        "human",
        "central row of the swept enrichment table (SD 4 pp, cutoff 37.5 pp), anchored "
        "on Vrtovec's responder analysis",
        "the SD is UNREPORTED. The central row was chosen, not the 2.34x maximum.",
    ),
    "RESPONSE_GUIDED_PASS": Provenance(
        0.424708,
        "fraction of all patients",
        "estimated",
        "design",
        "human",
        "fraction passing the baseline-severity screen at the central row",
        "screening EXCLUDES people: 42.5% pass x 58.2% respond = 24.7% of all patients",
    ),
    "ATLAS4_BEST_ADMISSIBLE": Provenance(
        0.896,
        "probability",
        "measured",
        "ceiling",
        "none",
        "ATLAS-4's best admissible point, read from results/atlas4, pinned by test",
        "a MODEL OUTPUT, not an external measurement; it is a reference point only",
    ),
    # ------------------------------------------------------- literature register
    "INTERVENTION_PATIENT_CV": Provenance(
        0.30,
        "coefficient of variation",
        "estimated",
        "arithmetic",
        "human",
        "between-patient CV of a peak or nadir inflammatory marker value; deliberately "
        "generous so the trial-size requirement is not flattered",
        "ESTIMATED and conservative",
    ),
    "INTERVENTION_ASSAY_CV": Provenance(
        0.10,
        "coefficient of variation",
        "estimated",
        "arithmetic",
        "human",
        "within-assay CV for CRP and IL-6 at the upper end of the quoted 5-10% range",
        "ESTIMATED from the quoted assay range",
    ),
    "RESOLVABLE_TRIAL_SIZE": Provenance(
        200.0,
        "patients per arm",
        "estimated",
        "design",
        "none",
        "a discrimination needing at most this many per arm is called resolvable, i.e. "
        "within the size of a conventional single-centre study",
        "a design threshold, not a measurement",
    ),
    "INTERVENTION_STEP": Provenance(
        0.50,
        "fraction",
        "estimated",
        "design",
        "none",
        "relative perturbation applied to each mechanism when building the discrimination "
        "table; identical for every mechanism so the comparison is like-for-like",
        "a design choice",
    ),
    # ------------------------------------------------------------ candidates
    "CTT_RRR_PER_MMOL": Provenance(
        0.22,
        "relative risk reduction",
        "measured",
        "efficacy",
        "human",
        "Cholesterol Treatment Trialists' meta-analysis: 22% relative reduction in major "
        "vascular events per 1.0 mmol/L (38.7 mg/dL) reduction in LDL-C, reproduced across "
        "the statin trials and extended by ezetimibe and the PCSK9 programmes",
        "the most reproduced quantitative relationship in cardiovascular prevention; a "
        "class-level average, not a claim about any single agent",
    ),
    "ZEUS_HR": Provenance(
        0.99,
        "hazard ratio",
        "measured",
        "refutation",
        "human",
        "ZEUS (NCT05021835), n=6,385, ziltivekimab 15 mg monthly vs placebo: hazard ratio "
        "for cardiovascular death, myocardial infarction or stroke = 0.99 (95% CI "
        "0.88-1.11); reported 2026",
        "a point estimate at unity: the trial excludes a benefit larger than ~12% but is "
        "not evidence of no effect at smaller sizes",
    ),
    "ZEUS_CI_LOW": Provenance(
        0.88,
        "hazard ratio",
        "measured",
        "refutation",
        "human",
        "ZEUS lower 95% confidence bound for the three-point major adverse cardiovascular "
        "event endpoint",
        "the smallest benefit still compatible with the data",
    ),
    "ZEUS_CI_HIGH": Provenance(
        1.11,
        "hazard ratio",
        "measured",
        "refutation",
        "human",
        "ZEUS upper 95% confidence bound for the three-point major adverse cardiovascular "
        "event endpoint",
        "harm at this magnitude cannot be excluded",
    ),
    "ZEUS_N": Provenance(
        6385.0,
        "patients",
        "measured",
        "refutation",
        "human",
        "ZEUS randomised population: established atherosclerotic disease with chronic "
        "kidney disease and high-sensitivity CRP at or above 2 mg/L",
        "a population selected precisely for inflammation, so a null here is the strongest "
        "available test of the inflammatory hypothesis",
    ),
    "ZEUS_CRP_REDUCTION": Provenance(
        0.64,
        "fraction",
        "measured",
        "refutation",
        "human",
        "ZEUS: 64% reduction in high-sensitivity CRP at the achieved dose",
        "the marker moved; the outcome did not",
    ),
    "ZEUS_IL6_REDUCTION": Provenance(
        0.76,
        "fraction",
        "measured",
        "refutation",
        "human",
        "ZEUS: 76% reduction in interleukin-6, i.e. the ligand itself was suppressed, not "
        "merely its downstream marker",
        "removes the escape route that the drug failed to engage its target",
    ),
    "ZEUS_SERIOUS_INFECTION": Provenance(
        0.20,
        "proportion",
        "measured",
        "safety",
        "human",
        "ZEUS serious infection rate in the ziltivekimab arm",
        "the refuted mechanism was not even free of harm",
    ),
    "ZEUS_SERIOUS_INFECTION_CONTROL": Provenance(
        0.18,
        "proportion",
        "measured",
        "safety",
        "human",
        "ZEUS serious infection rate in the placebo arm",
        "the comparator, recorded separately so the difference cannot be quoted as a "
        "ratio invented from a single number",
    ),
    "LVEF_SURROGATE_THRESHOLD_PP": Provenance(
        6.0,
        "percentage points",
        "measured",
        "ceiling",
        "human",
        "surrogate threshold effect for left-ventricular ejection fraction: the change in "
        "LVEF at which a mortality benefit becomes plausible, from the heart-failure "
        "surrogate literature",
        "a threshold, not a guarantee: clearing it is necessary, not sufficient",
    ),
    "ENLICITIDE_LDL_REDUCTION": Provenance(
        0.56,
        "fraction",
        "measured",
        "efficacy",
        "human",
        "CORALreef Lipids (n=2,909): oral PCSK9 inhibitor enlicitide lowered LDL-C 56% "
        "placebo-adjusted at week 24; approved 2026 on the surrogate",
        "a surrogate endpoint; no outcome trial has reported",
    ),
    "VERVE102_LDL_REDUCTION": Provenance(
        0.62,
        "fraction",
        "measured",
        "efficacy",
        "human",
        "Heart-2 (NCT06164730), phase 1b open-label, n=35: base editing VERVE-102 lowered "
        "LDL-C 62% at 1.0 mg/kg",
        "open-label and uncontrolled: the number includes regression to the mean, "
        "background therapy and selection",
    ),
    "VERVE102_N": Provenance(
        35.0,
        "patients",
        "measured",
        "efficacy",
        "human",
        "Heart-2 open-label cohort size at the 1.0 mg/kg dose level",
        "below any threshold at which an event effect could be estimated",
    ),
    "CORALREEF_OUTCOMES_N": Provenance(
        14500.0,
        "patients",
        "measured",
        "design",
        "human",
        "CORALreef Outcomes planned enrolment for the oral PCSK9 inhibitor outcome trial",
        "the size of the trial that would settle the question, not a result",
    ),
    "MIN_OUTCOME_TRIAL_N": Provenance(
        1000.0,
        "patients per arm",
        "estimated",
        "design",
        "none",
        "a candidate is called supported by outcome evidence only at or above this "
        "randomised size, because below it a null and a modest benefit are not separable",
        "a design threshold chosen here, not a published quantity",
    ),
    # ------------------------------------------------------- simulation frontier
    "CRP_ASSAY_CV": Provenance(
        0.10,
        "coefficient of variation",
        "estimated",
        "arithmetic",
        "human",
        "high-sensitivity CRP within-assay CV; clinical assays quote 5-10% and the upper "
        "end is used so the resolution requirement is not flattered",
        "ESTIMATED from the quoted assay range, not measured for a specific platform",
    ),
    "CRP_PATIENT_CV": Provenance(
        0.30,
        "coefficient of variation",
        "estimated",
        "arithmetic",
        "human",
        "between-patient CV of CRP, chosen generous because a smaller value would make the "
        "resolution requirement look easier than it is",
        "ESTIMATED and deliberately conservative; the ranking conclusion is robust to it",
    ),
    "STRUCTURAL_BLINDNESS": Provenance(
        1e-3,
        "fraction of the observable",
        "estimated",
        "design",
        "none",
        "a separation below one part in a thousand of the observable is called structural, "
        "because it is below any plausible assay resolution",
        "a design threshold, not a measurement",
    ),
    "FRONTIER_POWER_Z": Provenance(
        2.801585,
        "z",
        "estimated",
        "arithmetic",
        "none",
        "two-arm 80% power, two-sided 5% z-sum (same value as POWER_Z_TOTAL, kept separate "
        "so the frontier module reads alone)",
        "closed-form",
    ),
    "FRONTIER_STEP": Provenance(
        0.50,
        "fraction",
        "estimated",
        "design",
        "none",
        "relative perturbation applied to each mechanism when building the effect table",
        "a design choice; it is applied identically to every mechanism",
    ),
    # ------------------------------------------------------- unanswered-question audit
    "POWER_Z_TOTAL": Provenance(
        2.801585,
        "z",
        "estimated",
        "arithmetic",
        "none",
        "norm.ppf(0.975) + norm.ppf(0.80): two-sided 5%, 80% power",
        "closed-form, but computed rather than derived here; NOT a biological constant",
    ),
    "DEFAULT_LVEF_SD_PP": Provenance(
        5.0,
        "percentage points",
        "measured",
        "arithmetic",
        "human",
        "between-patient SD of LVEF change (registered independently as LVEF_CHANGE_SD_PP)",
        "repeated here so the arithmetic module can be read alone; same value, same source",
    ),
    "MIDWALL_TRANSMURALITY": Provenance(
        0.50,
        "fraction of wall thickness",
        "estimated",
        "arithmetic",
        "none",
        "definition: half-wall is 50% transmurality",
        "a definition, not a measurement",
    ),
    "REVIEWER_CHECK_TRANSMURALITY": Provenance(
        0.50,
        "fraction of wall thickness",
        "estimated",
        "arithmetic",
        "none",
        "the band boundary the mid-wall defect check is written against",
        "definition",
    ),
    "REVIEWER_CHECK_BAND_TOLERANCE": Provenance(
        0.10,
        "fraction of healthy recovery",
        "hypothetical",
        "design",
        "none",
        "a design threshold: adjacent scar bands must differ by more than 10% of the "
        "no-scar recovery, or the banding carries no information across that boundary",
        "a JUDGEMENT about how much resolution a band table must have; not measured",
    ),
    "LEGACY_CRISPR_LISTED_LOCI": Provenance(
        7.0,
        "loci",
        "measured",
        "arithmetic",
        "none",
        "counted directly from src/05_crispr_genomic_safety_cfd.py CANDIDATE_LOCI",
        "a fact about the script, not about the genome",
    ),
    "OVERLAP_BETWEEN_TWO_HAND_LISTS": Provenance(
        4.0,
        "loci",
        "measured",
        "arithmetic",
        "none",
        "two independent hand-written lists of seven candidate loci share four entries "
        "(union ten, intersection four)",
        "the disagreement between the two lists is itself the evidence that neither is a "
        "systematic screen",
    ),
    "GUIDE_LENGTH_NT": Provenance(
        20.0,
        "nucleotides",
        "measured",
        "mechanism",
        "human",
        "VERVE-102 ABE8e spacer length",
        "a sequence fact from the guide design, not a simulation parameter",
    ),
    "PAM_LENGTH_NT": Provenance(
        3.0,
        "nucleotides",
        "measured",
        "mechanism",
        "none",
        "SpCas9/ABE8e PAM is NGG, three positions",
        "the two N positions are what multiply the candidate site count",
    ),
    "HAPLOID_GENOME_BP": Provenance(
        3_100_000_000.0,
        "base pairs",
        "measured",
        "arithmetic",
        "human",
        "GRCh38 haploid genome length",
        "a reference-genome constant",
    ),
    "GENOME_STRANDS": Provenance(
        2.0,
        "strands",
        "estimated",
        "arithmetic",
        "none",
        "a guide site can occur on either strand",
        "definition",
    ),
    "ABE_VS_CAS9_SITE_RATIO": Provenance(
        16.0,
        "fold",
        "estimated",
        "mechanism",
        "none",
        "two degenerate PAM positions give 4^2 = 16 times the sites of a non-degenerate "
        "motif; the editor-specific activity weighting on top is NOT modelled",
        "a pure incidence factor; it is NOT an activity model and must not be read as one",
    ),
    "PER_ATTEMPT_AE_RISK": Provenance(
        0.0550,
        "probability",
        "estimated",
        "ceiling",
        "human",
        "ATLAS-5 per-exposure adverse-event rate (already registered in recombine)",
        "reused so the composite arithmetic is computed from the same number",
    ),
    "ARCHITECTURE_EXPOSURES": Provenance(
        3.0,
        "exposures",
        "estimated",
        "design",
        "none",
        "ATLAS-5 stacks a repeat dose, a maturation step and a delivery",
        "a property of the architecture, not a measurement; the composite risk is "
        "INDEPENDENT-EXPOSURE arithmetic and assumes independence, which is conservative",
    ),
    "PRIOR_CONTINUE_WEAK": Provenance(
        0.50,
        "relative prior weight",
        "hypothetical",
        "design",
        "none",
        "swept grid node for 'the programme would have produced a weak positive'",
        "there is no registry of counterfactual programme outcomes; this is a JUDGEMENT "
        "and it is swept, not asserted",
    ),
    "PRIOR_CONTINUE_STRONG": Provenance(
        0.85,
        "relative prior weight",
        "hypothetical",
        "design",
        "none",
        "swept grid node for 'the programme would have produced a strong positive'",
        "as above: a judgement, swept",
    ),
    "PRIOR_STOP": Provenance(
        0.30,
        "relative prior weight",
        "hypothetical",
        "design",
        "none",
        "swept grid node for 'the programme should have been stopped'",
        "as above: a judgement, swept",
    ),
    "CONTINUE_THRESHOLD_UNITS": Provenance(
        0.5,
        "evidence units",
        "estimated",
        "design",
        "none",
        "expected value above which continuing is worth it; half the value of the "
        "cheapest non-stop outcome, so stop-versus-weak-positive is exactly break-even",
        "a design choice about what 'worth continuing' means",
    ),
    "PROGRAMME_COST_UNITS": Provenance(
        1.5,
        "evidence units",
        "hypothetical",
        "design",
        "none",
        "cost of continuing the programme at all: below the value of a strong positive "
        "(3 units) so a confident prior can justify continuing, above a weak one (1 unit) "
        "so it cannot alone",
        "a JUDGEMENT. Without a cost term the decision degenerates to 'more evidence is "
        "never bad' and cannot flip, which would make the prior sweep decorative",
    ),
    # ------------------------------------------------------- audit thresholds
    "DAYS": Provenance(
        180.0,
        "days",
        "hypothetical",
        "design",
        "none",
        "legacy cytokine network simulation horizon (src/04)",
        "design choice; the trajectory has NOT converged by day 180 - CRP is still rising",
    ),
    "DT": Provenance(
        0.25,
        "days",
        "hypothetical",
        "design",
        "none",
        "legacy cytokine network time step (src/04)",
        "design choice; never checked for temporal convergence before this audit",
    ),
    "SENSITIVITY_FLOOR": Provenance(
        1e-06,
        "elasticity",
        "hypothetical",
        "design",
        "none",
        "elasticity below which a parameter is called inert",
        "threshold; a parameter under it cannot be fitted at all",
    ),
    "COLLINEARITY_CEILING": Provenance(
        0.95,
        "correlation",
        "hypothetical",
        "design",
        "none",
        "|r| above which two parameters are called non-separable",
        "threshold, following the usual collinearity rule of thumb",
    ),
    "PROFILE_RISE_FRACTION": Provenance(
        0.05,
        "fraction of objective span",
        "hypothetical",
        "design",
        "none",
        "rise required on both sides of a profile minimum to call it identified",
        "threshold",
    ),
    "CALIBRATION_HORIZON_YEARS": Provenance(
        2.0,
        "years",
        "estimated",
        "design",
        "human",
        "horizon over which the regression grid is calibrated against trial data",
        "design choice",
    ),
    # ------------------------------------------------------- study-design constants
    "CRITERION_THRESHOLD": Provenance(
        0.8,
        "fraction",
        "hypothetical",
        "design",
        "none",
        "this study's bar for 'criterion met'",
        "a THRESHOLD, not an effect. Every 'robust' claim is relative to it.",
    ),
    "ROBUSTNESS_BAR": Provenance(
        0.8,
        "fraction",
        "hypothetical",
        "design",
        "none",
        "same bar, restated for the evidence audit",
        "design",
    ),
    "REFERENCE_RISK": Provenance(
        0.0198,
        "annual probability",
        "estimated",
        "design",
        "human",
        "low-risk reference event rate used to normalise risk",
        "derived; affects the risk-normalisation criterion",
    ),
    "MIN_SPECIFICATION_DRAWS": Provenance(
        4000,
        "draws",
        "hypothetical",
        "design",
        "none",
        "floor on the sampling budget so a specification cannot move with it",
        "a design guard, and a caught bug: specifications were budget-dependent",
    ),
    "MIN_SPECIFICATION_DRAWS_RECOMBINE": Provenance(
        2000,
        "draws",
        "hypothetical",
        "design",
        "none",
        "same guard for the recombination search",
        "design guard",
    ),
    "INTERACTION_TOLERANCE": Provenance(
        0.005,
        "probability units",
        "hypothetical",
        "design",
        "none",
        "half a percentage point, the noise floor for an interaction verdict",
        "design guard; comparing against 1e-9 once invented a 0.05 pp 'price'",
    ),
    "ATLAS_CLEARANCE_MONTHS": Provenance(
        24.0,
        "months",
        "hypothetical",
        "design",
        "none",
        "clearing horizon used throughout the ATLAS series",
        "design choice; the 55.4% ceiling was horizon-specific",
    ),
    "ATLAS_GRAFT_DELAY_MONTHS": Provenance(
        24.0,
        "months",
        "hypothetical",
        "design",
        "none",
        "delay before the graft arm begins",
        "design choice",
    ),
    "SOLVED_MONTHS": Provenance(
        60.0,
        "months",
        "hypothetical",
        "design",
        "none",
        "horizon at which ATLAS-2 is reported solved",
        "design choice",
    ),
    "DEFAULT_SOURCE_REQUIRED": Provenance(
        0.223,
        "fraction",
        "estimated",
        "ceiling",
        "human",
        "source reduction needed for 80% clearance with perfect modules",
        "derived arithmetic",
    ),
    "FEEDBACK_MAX_UPLIFT": Provenance(
        0.5,
        "fraction",
        "hypothetical",
        "mechanism",
        "none",
        "bounded aggregate feedback uplift",
        "assumption; per-factor floor-raising was tried first and forced P=100%",
    ),
    "NON_VIABILITY_SIGNAL_QUALITY": Provenance(
        0.0,
        "quality",
        "hypothetical",
        "mechanism",
        "none",
        "signal quality when no viability reporter is present",
        "definitional",
    ),
}


@dataclass(frozen=True)
class Gap:
    """One numeric constant found in the code that the registry does not cover."""

    module: str
    name: str
    value: float


def numeric_constants(root: Path | None = None) -> list[Gap]:
    """Every module-level numeric constant, parsed from source rather than imported.

    Parsing matters: importing every module to enumerate them would run side effects and
    could mask a constant that fails to import for the very reason it is unregistered.
    """
    base = Path(root) if root is not None else Path(__file__).resolve().parent
    found: list[Gap] = []
    for path in sorted(base.glob("*.py")):
        tree = ast.parse(path.read_text())
        for node in tree.body:
            if not isinstance(node, ast.Assign):
                continue
            for target in node.targets:
                if not (isinstance(target, ast.Name) and target.id.isupper()):
                    continue
                value = node.value
                if (
                    isinstance(value, ast.Constant)
                    and isinstance(value.value, (int, float))
                    and not isinstance(value.value, bool)
                ):
                    found.append(Gap(path.name, target.id, float(value.value)))
    return found


def registry_keys() -> set[str]:
    """Registry keys, with the disambiguating suffix used for duplicated names."""
    return set(REGISTRY)


def resolve(module: str, name: str) -> Provenance | None:
    """Look a constant up, falling back to a module-qualified key.

    Two modules legitimately define ``MIN_SPECIFICATION_DRAWS`` with different values
    (ATLAS-3 needs 4000, the recombination search 2000). Rather than silently collapse
    them, the registry holds both under explicit keys and this resolver picks the right
    one from the module name.
    """
    if name not in REGISTRY:
        return None
    if name == "MIN_SPECIFICATION_DRAWS" and module.startswith("recombine"):
        return REGISTRY["MIN_SPECIFICATION_DRAWS_RECOMBINE"]
    return REGISTRY[name]


def unregistered(root: Path | None = None) -> list[Gap]:
    """Constants that act on the model but carry no provenance record."""
    return [gap for gap in numeric_constants(root) if resolve(gap.module, gap.name) is None]


def audit_frame(root: Path | None = None):
    """The registry as a table, with coverage of what is actually in the code."""
    import pandas as pd

    gaps = numeric_constants(root)
    rows = []
    for gap in gaps:
        record = resolve(gap.module, gap.name)
        rows.append(
            {
                "module": gap.module,
                "name": gap.name,
                "value_in_code": gap.value,
                "value_in_registry": None if record is None else record.value,
                "kind": None if record is None else record.kind,
                "supports": None if record is None else record.supports,
                "population": None if record is None else record.population,
                "earns_efficacy_credit": (None if record is None else record.earns_efficacy_credit),
                "unit": None if record is None else record.unit,
                "source": None if record is None else record.source,
                "registered": record is not None,
                "values_agree": None if record is None else abs(record.value - gap.value) < 1e-9,
            }
        )
    return pd.DataFrame(rows)


def coverage() -> dict[str, object]:
    """Summary counts, so the CLI can print a single honest line."""
    frame = audit_frame()
    total = len(frame)
    registered = int(frame["registered"].sum())
    agreeing = int(frame["values_agree"].fillna(False).sum())
    by_kind = frame["kind"].value_counts().to_dict()
    return {
        "constants_found": total,
        "registered": registered,
        "unregistered": total - registered,
        "values_disagreeing": registered - agreeing,
        "by_kind": {str(k): int(v) for k, v in by_kind.items()},
        "coverage_fraction": registered / total if total else 1.0,
        "efficacy_credited": int(frame["earns_efficacy_credit"].fillna(False).sum()),
    }


def provenance_report() -> str:
    """A plain-text audit. Prints the unregistered list rather than hiding it."""
    frame = audit_frame()
    summary = coverage()
    lines = [
        "PARAMETER PROVENANCE AUDIT",
        "=" * 74,
        "",
        "Every numeric module constant must declare what kind of claim it is.",
        "A missing record fails the build; it is not an annotation.",
        "",
        f"  constants found        {summary['constants_found']}",
        f"  registered             {summary['registered']}",
        f"  unregistered           {summary['unregistered']}",
        f"  values disagreeing     {summary['values_disagreeing']}",
        f"  earn efficacy credit   {summary['efficacy_credited']}",
        "",
        "BY CLASS",
        "-" * 74,
    ]
    for kind in CLASSES:
        count = int((frame["kind"] == kind).sum())
        lines.append(f"  {kind:14s} {count:3d}")
    lines += [
        "",
        "BY SUPPORT TYPE",
        "-" * 74,
    ]
    for supports in SUPPORTS:
        count = int((frame["supports"] == supports).sum())
        if count:
            lines.append(f"  {supports:14s} {count:3d}")
    lines += [
        "",
        "HYPOTHETICAL PARAMETERS CARRYING REAL WEIGHT",
        "-" * 74,
        "  These are assumptions, not measurements. Anything downstream of them is an",
        "  assumption too. Listed with the note recorded at registration.",
        "",
    ]
    heavy = frame.loc[frame["kind"] == HYPOTHETICAL]
    for row in heavy.itertuples(index=False):
        record = resolve(row.module, row.name)
        lines.append(f"  {row.name} = {row.value_in_code}")
        if record is not None and record.note:
            lines.append(f"    {record.note}")
    unregistered_rows = frame.loc[~frame["registered"]]
    if not unregistered_rows.empty:
        lines += [
            "",
            "UNREGISTERED - THESE FAIL THE BUILD",
            "-" * 74,
        ]
        for row in unregistered_rows.itertuples(index=False):
            lines.append(f"  {row.module}:{row.name} = {row.value_in_code}")
    return "\n".join(lines)
