"""Is any of this a *cure*? An explicit, auditable answer.

"Cure" is not a biomarker. This module states five criteria, scores every
catalogued intervention against them from the encoded evidence, then searches for
the smallest combination of pathways that could satisfy all five. It also
quantifies the one genuinely cure-adjacent idea in the literature - removing
cumulative LDL exposure early and permanently - using the Mendelian-randomisation
contrast rather than pretending a five-year trial says something about a lifetime.

The five criteria
-----------------
1. ``durable_driver_control``   - causal driver suppressed without ongoing dosing
2. ``burden_reversal``          - existing atherosclerotic burden removed
3. ``myocardial_restoration``   - lost contractile myocardium replaced
4. ``risk_normalisation``       - event risk returned to a low-risk reference
5. ``randomised_proof``         - benefit proven on hard endpoints in humans

A criterion counts as satisfied at >= 0.8. Scores are capped by what the evidence
actually shows: partial plaque regression cannot score 1.0, and no single agent
normalises risk, so the ceiling matters more than any individual score.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum
from itertools import combinations

import numpy as np
import pandas as pd

from ihd_solutions.catalog import DosingProfile, Regimen, Therapy, TherapyClass
from ihd_solutions.evidence import (
    CreditMode,
    EndpointClass,
    EvidenceTier,
    TrialAnchor,
)
from ihd_solutions.population import Population
from ihd_solutions.risk import simulate_events
from ihd_solutions.units import MG_DL_PER_MMOL_L, require_fraction, require_positive

#: Score at which a criterion counts as met.
CRITERION_THRESHOLD = 0.8

#: Equal weights, stated rather than hidden. Change them and re-run; the binding
#: constraints below are ceilings, so no reasonable weighting removes them.
CRITERION_WEIGHTS: dict[str, float] = {
    "durable_driver_control": 0.2,
    "burden_reversal": 0.2,
    "myocardial_restoration": 0.2,
    "risk_normalisation": 0.2,
    "randomised_proof": 0.2,
}


class CureCriterion(str, Enum):
    DURABLE_DRIVER_CONTROL = "durable_driver_control"
    BURDEN_REVERSAL = "burden_reversal"
    MYOCARDIAL_RESTORATION = "myocardial_restoration"
    RISK_NORMALISATION = "risk_normalisation"
    RANDOMISED_PROOF = "randomised_proof"


CRITERION_DESCRIPTIONS: dict[CureCriterion, str] = {
    CureCriterion.DURABLE_DRIVER_CONTROL: (
        "The causal driver stays suppressed without repeated dosing (one-time editing is the "
        "only modality that plausibly reaches this)."
    ),
    CureCriterion.BURDEN_REVERSAL: (
        "Existing atherosclerotic burden is removed, not merely stabilised. Imaging shows "
        "partial lipid-core regression with calcified volume rising, so this is capped."
    ),
    CureCriterion.MYOCARDIAL_RESTORATION: (
        "Lost contractile myocardium is replaced. Engineered heart muscle shows a few "
        "percentage points of LVEF, not restored normal myocardium."
    ),
    CureCriterion.RISK_NORMALISATION: (
        "Event risk returns to a low-risk reference. No single therapy does this; the best "
        "stacks remove roughly half of the excess risk."
    ),
    CureCriterion.RANDOMISED_PROOF: (
        "Benefit proven on hard human endpoints. Biomarker, animal and pending results "
        "score low or zero."
    ),
}

#: Structural reversal, scored from the imaging anchors only. Ceilings are set by
#: the fact that regression is partial and calcified burden persists.
BURDEN_REVERSAL_EVIDENCE: dict[str, tuple[float, str]] = {
    "statin_high_intensity": (0.50, "ASTEROID: PAV -0.98 pp, 63.6% of patients regressed"),
    "evolocumab": (0.55, "GLAGOV: PAV -0.95 pp, 64.3% regressed (81.2% if baseline LDL <70)"),
    "ezetimibe": (
        0.60,
        "FLOWPROMOTE: with rosuvastatin, noncalcified -28.2%, 98% of patients regressed; "
        "calcified volume unchanged",
    ),
}

#: Myocardial restoration, scored from human functional/structural data only.
MYOCARDIAL_RESTORATION_EVIDENCE: dict[str, tuple[float, str]] = {
    "biovat_hf_patch": (
        0.50,
        "BioVAT-HF: wall thickness +4.5 mm, LVEF +3.9 pp at 3 months; remuscularisation "
        "confirmed in one explanted heart; n=16, immunosuppression required",
    ),
    "ipsc_cm_sheet": (0.20, "RiHEART: paracrine mechanism, LVEF response in 2/8 patients"),
    "yap101_sav1_gene_therapy": (0.15, "Pig MI model EF +14%; first-in-human trial only"),
    "cd34_cell_therapy": (0.10, "Microvascular/angiogenic effect, not new myocardium"),
}

_PROOF_SCORES: dict[EvidenceTier, float] = {
    EvidenceTier.META_ANALYSIS: 1.0,
    EvidenceTier.RANDOMIZED_OUTCOMES: 1.0,
    EvidenceTier.RANDOMIZED_TOPLINE: 0.5,
    EvidenceTier.RANDOMIZED_ACTIVE_COMPARATOR: 0.4,
    EvidenceTier.RANDOMIZED_SURROGATE: 0.3,
    EvidenceTier.OBSERVATIONAL: 0.15,
    EvidenceTier.CONDITIONAL_APPROVAL: 0.15,
    EvidenceTier.MENDELIAN_RANDOMIZATION: 0.2,
    EvidenceTier.MECHANISTIC: 0.0,
}

#: No single agent normalises risk, so this criterion is capped for one therapy.
MAX_SINGLE_AGENT_RISK_NORMALISATION = 0.4

#: Hard floor on any scenario relative risk, matching the event-model guard.
MIN_SCENARIO_RR = 0.1


def _durable_driver_score(therapy: Therapy) -> tuple[float, str]:
    """How much of the causal driver stays controlled without repeat dosing."""
    drives_lipids = therapy.ldl_pct > 0 or therapy.apob_pct > 0 or therapy.lpa_pct > 0
    drives_inflammation = therapy.hscrp_pct > 0
    if therapy.profile is DosingProfile.ONE_TIME_EDIT and drives_lipids:
        if therapy.is_permanent_scenario:
            score, reason = 1.0, "one-time edit modelled as permanent in the scenario"
        else:
            score, reason = 0.6, "one-time edit with an explicit waning half-life"
    elif not (drives_lipids or drives_inflammation):
        return 0.0, "does not control the causal driver"
    elif therapy.profile is DosingProfile.INTERVAL_SC and therapy.interval_days >= 90:
        score, reason = 0.35, "infrequent dosing (>=90 days) but still repeated"
    elif therapy.profile is DosingProfile.INTERVAL_SC:
        score, reason = 0.20, "repeated parenteral dosing"
    elif therapy.therapy_class is TherapyClass.LIFESTYLE:
        score, reason = 0.15, "requires sustained behaviour, not a durable biological change"
    elif therapy.profile is DosingProfile.ONE_TIME_PROCEDURE:
        score, reason = 0.30, "one-time procedure with a sustained biological effect"
    else:
        score, reason = 0.10, "daily dosing; effect ends when the drug stops"
    if therapy.evidence_tier is EvidenceTier.MECHANISTIC:
        score *= 0.3
        reason += "; animal/in-vitro evidence only"
    return score, reason


@dataclass(frozen=True)
class CureAssessment:
    """One therapy scored against the five cure criteria."""

    therapy_id: str
    scores: dict[str, float]
    rationale: dict[str, str]

    def __post_init__(self) -> None:
        expected = {criterion.value for criterion in CureCriterion}
        if set(self.scores) != expected:
            raise ValueError(f"scores must cover {sorted(expected)}")
        for name, value in self.scores.items():
            require_fraction(name, value, allow_one=True)
        if not math.isclose(sum(CRITERION_WEIGHTS.values()), 1.0, abs_tol=1e-9):
            raise ValueError("criterion weights must sum to 1")

    @property
    def readiness_index(self) -> float:
        return float(sum(CRITERION_WEIGHTS[name] * value for name, value in self.scores.items()))

    @property
    def satisfied(self) -> tuple[str, ...]:
        return tuple(name for name, value in self.scores.items() if value >= CRITERION_THRESHOLD)

    @property
    def blocking(self) -> tuple[str, ...]:
        return tuple(name for name, value in self.scores.items() if value < CRITERION_THRESHOLD)


def assess_therapy(
    therapy: Therapy,
    anchors: dict[str, TrialAnchor],
    strict_rrr: float,
) -> CureAssessment:
    """Score one therapy. ``strict_rrr`` is its simulated strict-mode RRR."""
    scores: dict[str, float] = {}
    rationale: dict[str, str] = {}

    scores["durable_driver_control"], rationale["durable_driver_control"] = _durable_driver_score(
        therapy
    )

    reversal, reversal_reason = BURDEN_REVERSAL_EVIDENCE.get(
        therapy.therapy_id, (0.0, "no human plaque-regression evidence for this agent")
    )
    scores["burden_reversal"] = reversal
    rationale["burden_reversal"] = reversal_reason

    restoration, restoration_reason = MYOCARDIAL_RESTORATION_EVIDENCE.get(
        therapy.therapy_id, (0.0, "does not replace lost myocardium")
    )
    scores["myocardial_restoration"] = restoration
    rationale["myocardial_restoration"] = restoration_reason

    normalisation = MAX_SINGLE_AGENT_RISK_NORMALISATION * max(0.0, strict_rrr)
    scores["risk_normalisation"] = float(normalisation)
    rationale["risk_normalisation"] = (
        f"strict-mode relative risk reduction {strict_rrr:.1%}; a single agent cannot "
        "return risk to a low-risk reference"
    )

    proof = _PROOF_SCORES[therapy.evidence_tier]
    hard = [
        anchor
        for anchor in anchors.values()
        if anchor.therapy_id == therapy.therapy_id
        and anchor.endpoint_class in {EndpointClass.MACE, EndpointClass.CV_DEATH}
    ]
    if proof >= 1.0 and strict_rrr <= 0.02:
        proof = 0.5
        rationale["randomised_proof"] = (
            "randomised hard-endpoint evidence exists but shows no benefit (neutral or negative)"
        )
    else:
        rationale["randomised_proof"] = f"tier {therapy.evidence_tier.value}" + (
            f"; anchors {', '.join(sorted(a.anchor_id for a in hard))}" if hard else ""
        )
    scores["randomised_proof"] = float(proof)
    return CureAssessment(therapy.therapy_id, scores, rationale)


def assess_catalog(
    catalog: dict[str, Therapy],
    anchors: dict[str, TrialAnchor],
    population: Population,
    days: np.ndarray,
) -> pd.DataFrame:
    """Score every catalogued intervention, simulating each one on its own."""
    rows = []
    for therapy_id, therapy in sorted(catalog.items()):
        try:
            simulation = simulate_events(
                Regimen((therapy_id,)),
                catalog,
                population,
                anchors,
                days,
                mode=CreditMode.STRICT,
            )
            rrr = max(0.0, simulation.relative_risk_reduction)
        except ValueError:
            rrr = 0.0
        assessment = assess_therapy(therapy, anchors, rrr)
        rows.append(
            {
                "therapy_id": therapy_id,
                "therapy_class": therapy.therapy_class.value,
                "regulatory_status": therapy.regulatory_status,
                "strict_rrr_pct": 100.0 * rrr,
                **assessment.scores,
                "cure_readiness_index": assessment.readiness_index,
                "criteria_satisfied": ";".join(assessment.satisfied),
                "criteria_blocking": ";".join(assessment.blocking),
                "rationale": " | ".join(
                    f"{name}: {assessment.rationale[name]}" for name in assessment.scores
                ),
                "source_url": therapy.source_url,
            }
        )
    return (
        pd.DataFrame(rows)
        .sort_values(["cure_readiness_index", "therapy_id"], ascending=[False, True])
        .reset_index(drop=True)
    )


@dataclass(frozen=True)
class CoverResult:
    """Outcome of searching for a combination that satisfies every criterion."""

    threshold: float
    satisfied: bool
    selected: tuple[str, ...]
    coverage: dict[str, float]
    missing: tuple[str, ...]
    ceiling: dict[str, float]
    notes: tuple[str, ...] = field(default=())

    @property
    def binding_constraints(self) -> tuple[str, ...]:
        """Criteria that no catalogued therapy can meet, however combined."""
        return tuple(name for name, value in self.ceiling.items() if value < self.threshold)


def criterion_ceiling(assessments: pd.DataFrame) -> dict[str, float]:
    """Best score any single catalogued therapy achieves per criterion."""
    return {
        criterion.value: float(assessments[criterion.value].max()) for criterion in CureCriterion
    }


def search_cure_combination(
    assessments: pd.DataFrame,
    threshold: float = CRITERION_THRESHOLD,
    max_size: int = 5,
) -> CoverResult:
    """Smallest set of therapies whose per-criterion maxima clear ``threshold``.

    Coverage combines by maximum, not by sum: two half-effective regeneration
    procedures do not add up to a restored heart.
    """
    require_positive("threshold", threshold)
    require_fraction("threshold", threshold, allow_one=True)
    if max_size < 1:
        raise ValueError("max_size must be at least 1")
    names = [criterion.value for criterion in CureCriterion]
    ceiling = criterion_ceiling(assessments)
    records = assessments[["therapy_id", *names]].to_dict(orient="records")
    useful = [
        record for record in records if any(float(record[name]) >= threshold for name in names)
    ]

    def coverage(records_subset: list[dict]) -> dict[str, float]:
        return {
            name: max((float(record[name]) for record in records_subset), default=0.0)
            for name in names
        }

    for size in range(1, min(max_size, len(useful)) + 1):
        for subset in combinations(useful, size):
            values = coverage(list(subset))
            if all(values[name] >= threshold for name in names):
                return CoverResult(
                    threshold=threshold,
                    satisfied=True,
                    selected=tuple(sorted(str(record["therapy_id"]) for record in subset)),
                    coverage=values,
                    missing=(),
                    ceiling=ceiling,
                    notes=("criteria combined by maximum, never by summation",),
                )
    best = coverage(useful or records)
    missing = tuple(name for name in names if best[name] < threshold)
    notes = [
        "no combination of catalogued interventions satisfies every criterion",
        "criteria combined by maximum, never by summation",
    ]
    binding = tuple(name for name, value in ceiling.items() if value < threshold)
    if binding:
        notes.append("binding ceilings: " + ", ".join(binding))
    return CoverResult(
        threshold=threshold,
        satisfied=False,
        selected=(),
        coverage=best,
        missing=missing,
        ceiling=ceiling,
        notes=tuple(notes),
    )


@dataclass(frozen=True)
class LifetimeScenario:
    """Lifetime risk under early, durable LDL-C reduction.

    This is a *Mendelian-randomisation scenario*, not a trial result: it asks what
    removing cumulative exposure for a given share of adult life could do, using
    the genetic contrast (RR 0.46 per 1 mmol/L lower LDL-C from early life).
    """

    label: str
    start_age: int
    end_age: int
    exposure_start_age: int
    ldl_reduction_fraction: float
    delta_mmol_l: float
    exposure_fraction_covered: float
    relative_risk: float
    baseline_lifetime_risk: float
    scenario_lifetime_risk: float
    notes: tuple[str, ...] = field(default=())

    @property
    def relative_risk_reduction(self) -> float:
        return 1.0 - self.scenario_lifetime_risk / self.baseline_lifetime_risk

    @property
    def events_avoided_per_1000(self) -> float:
        return 1000.0 * (self.baseline_lifetime_risk - self.scenario_lifetime_risk)


def lifetime_exposure_scenario(
    start_age: int,
    ldl_reduction_fraction: float,
    *,
    label: str = "scenario",
    baseline_ldl_mg_dl: float = 160.0,
    baseline_lifetime_risk: float = 0.20,
    end_age: int = 80,
    exposure_start_age: int = 20,
    genetic_hr_per_mmol_l: float = 0.46,
) -> LifetimeScenario:
    """Interpolate the genetic lifetime contrast by the share of exposure removed.

    ``exposure_fraction_covered`` is the fraction of the adult exposure window
    (``exposure_start_age`` to ``end_age``) spent under treatment. The genetic
    hazard ratio is applied per mmol/L, scaled by that fraction. This is an
    explicit modelling assumption: genetics does not license a treatment claim.
    """
    if not exposure_start_age <= start_age < end_age:
        raise ValueError("require exposure_start_age <= start_age < end_age")
    require_fraction("ldl_reduction_fraction", ldl_reduction_fraction)
    require_fraction("baseline_lifetime_risk", baseline_lifetime_risk)
    require_positive("baseline_ldl_mg_dl", baseline_ldl_mg_dl)
    require_positive("genetic_hr_per_mmol_l", genetic_hr_per_mmol_l)
    if genetic_hr_per_mmol_l >= 1.0:
        raise ValueError("the genetic contrast must be protective (HR < 1)")

    delta_mmol = baseline_ldl_mg_dl * ldl_reduction_fraction / MG_DL_PER_MMOL_L
    covered = (end_age - start_age) / (end_age - exposure_start_age)
    exponent = delta_mmol * covered
    relative_risk = max(MIN_SCENARIO_RR, float(genetic_hr_per_mmol_l**exponent))
    baseline_hazard = -math.log1p(-baseline_lifetime_risk)
    scenario_risk = -math.expm1(-baseline_hazard * relative_risk)
    notes = [
        "Mendelian-randomisation scenario, not a randomised treatment effect",
        "assumes the achieved reduction is sustained for the whole covered window",
        "does not model procedural risk, off-target editing risk or competing mortality",
    ]
    if relative_risk <= MIN_SCENARIO_RR:
        notes.append("relative risk floor applied")
    return LifetimeScenario(
        label=label,
        start_age=int(start_age),
        end_age=int(end_age),
        exposure_start_age=int(exposure_start_age),
        ldl_reduction_fraction=float(ldl_reduction_fraction),
        delta_mmol_l=float(delta_mmol),
        exposure_fraction_covered=float(covered),
        relative_risk=relative_risk,
        baseline_lifetime_risk=float(baseline_lifetime_risk),
        scenario_lifetime_risk=float(scenario_risk),
        notes=tuple(notes),
    )


#: Modalities worth asking the lifetime question about, with their human LDL-C
#: effect sizes. Efficacy comes from trials; the lifetime framing does not.
LIFETIME_MODALITIES: tuple[tuple[str, float], ...] = (
    ("high_intensity_statin", 0.50),
    ("ctx310_angptl3_editing", 0.53),
    ("verve102_base_editing", 0.62),
    ("editing_plus_statin_stack", 0.81),
)


def lifetime_grid(
    start_ages: tuple[int, ...] = (30, 40, 50, 60, 70),
    modalities: tuple[tuple[str, float], ...] = LIFETIME_MODALITIES,
    **kwargs,
) -> pd.DataFrame:
    rows = []
    for label, fraction in modalities:
        for start_age in start_ages:
            scenario = lifetime_exposure_scenario(start_age, fraction, label=label, **kwargs)
            rows.append(
                {
                    "modality": label,
                    "start_age": start_age,
                    "ldl_reduction_pct": 100.0 * fraction,
                    "delta_ldl_mmol_l": scenario.delta_mmol_l,
                    "exposure_window_covered_pct": 100.0 * scenario.exposure_fraction_covered,
                    "scenario_relative_risk": scenario.relative_risk,
                    "baseline_lifetime_risk_pct": 100.0 * scenario.baseline_lifetime_risk,
                    "scenario_lifetime_risk_pct": 100.0 * scenario.scenario_lifetime_risk,
                    "lifetime_rrr_pct": 100.0 * scenario.relative_risk_reduction,
                    "events_avoided_per_1000": scenario.events_avoided_per_1000,
                    "evidence_class": "mendelian_randomisation_scenario",
                }
            )
    return pd.DataFrame(rows)


@dataclass(frozen=True)
class ResidualRisk:
    """How much excess risk survives the best evidence-backed stack."""

    population_id: str
    regimen: str
    baseline_risk_pct: float
    treated_risk_pct: float
    reference_risk_pct: float
    excess_risk_remaining_fraction: float

    @property
    def normalised(self) -> bool:
        return self.excess_risk_remaining_fraction <= 0.0


def residual_risk(
    simulation,
    population: Population,
    reference_annual_rate: float = 0.004,
) -> ResidualRisk:
    """Compare the treated curve with a low-risk reference over the same horizon."""
    require_positive("reference_annual_rate", reference_annual_rate)
    years = float(simulation.days[-1]) / 365.25
    reference_risk = -math.expm1(-reference_annual_rate * years)
    baseline = float(simulation.baseline_incidence[-1])
    treated = float(simulation.cumulative_incidence[-1])
    excess = baseline - reference_risk
    remaining = (treated - reference_risk) / excess if excess > 0 else 0.0
    return ResidualRisk(
        population_id=population.population_id,
        regimen=simulation.regimen_label,
        baseline_risk_pct=100.0 * baseline,
        treated_risk_pct=100.0 * treated,
        reference_risk_pct=100.0 * reference_risk,
        excess_risk_remaining_fraction=float(max(0.0, remaining)),
    )


def cure_search_report(
    assessments: pd.DataFrame,
    cover: CoverResult,
    lifetime: pd.DataFrame,
    residuals: list[ResidualRisk],
) -> str:
    """Plain-text bottom line, written to be quoted without losing the caveats."""
    lines = [
        "IHD cure search: what the evidence can and cannot support",
        "=" * 68,
        "VERDICT: no catalogued intervention, and no combination of them, meets all five",
        "cure criteria. This is a quantitative statement about the encoded evidence,",
        "not a claim that a cure is impossible in principle.",
        "",
        "CRITERIA AND CEILINGS (best score any single therapy achieves)",
    ]
    for criterion in CureCriterion:
        name = criterion.value
        lines.append(
            f"  {name}: ceiling {cover.ceiling[name]:.2f} "
            f"(threshold {cover.threshold:.2f}) - {CRITERION_DESCRIPTIONS[criterion]}"
        )
    if residuals:
        best_removed = max(1.0 - item.excess_risk_remaining_fraction for item in residuals)
        lines.append(
            f"  best observed stack removes {best_removed:.0%} of excess risk "
            f"(risk-normalisation score {best_removed:.2f}, still below threshold)"
        )
    lines.extend(["", "BINDING CONSTRAINTS"])
    if cover.binding_constraints:
        for name in cover.binding_constraints:
            lines.append(f"  {name}: no therapy reaches the threshold, so no stack can either")
    else:
        lines.append("  none individually; the search still found no covering combination")
    lines.extend(
        [
            "",
            "HIGHEST-SCORING INTERVENTIONS (readiness index, equal weights)",
        ]
    )
    for row in assessments.head(10).itertuples(index=False):
        lines.append(
            f"  {row.therapy_id}: {row.cure_readiness_index:.2f} "
            f"[satisfied: {row.criteria_satisfied or 'none'}]"
        )
    lines.extend(
        [
            "",
            "THE STRONGEST CURE-ADJACENT SCENARIO: EARLY, DURABLE EXPOSURE REMOVAL",
            "Mendelian-randomisation scenario, not a treatment effect. Genetic carriers of",
            "lifelong low LDL-C have far lower risk than five-year trials predict, so the",
            "question is when a durable intervention starts, not only how much it lowers.",
        ]
    )
    for row in lifetime.itertuples(index=False):
        if row.modality in {"verve102_base_editing", "editing_plus_statin_stack"}:
            lines.append(
                f"  {row.modality} from age {row.start_age}: lifetime risk "
                f"{row.scenario_lifetime_risk_pct:.1f}% vs {row.baseline_lifetime_risk_pct:.1f}% "
                f"({row.lifetime_rrr_pct:.0f}% lower; {row.events_avoided_per_1000:.0f} "
                "events avoided per 1000)"
            )
    lines.extend(["", "RESIDUAL RISK UNDER THE BEST EVIDENCE-BACKED STACKS"])
    for item in residuals:
        lines.append(
            f"  {item.population_id}: {item.baseline_risk_pct:.1f}% -> {item.treated_risk_pct:.1f}% "
            f"vs reference {item.reference_risk_pct:.1f}%; "
            f"{100.0 * item.excess_risk_remaining_fraction:.0f}% of excess risk remains"
        )
    lines.extend(
        [
            "",
            "WHAT WOULD HAVE TO CHANGE",
            "- Burden reversal: trials showing removal of calcified/total burden, not just",
            "  the lipid-rich fraction. FLOWPROMOTE/GLAGOV show partial regression only.",
            "- Myocardial restoration: controlled trials showing durable, arrhythmia-free",
            "  remuscularisation with normalised function, not a few LVEF points.",
            "- Risk normalisation: evidence that treated risk equals a low-risk reference.",
            "- Randomised proof for one-time editing: outcomes trials, plus off-target",
            "  assays that DSB-based methods cannot provide.",
            "",
            "NON-CLAIMS",
            "- No cure is demonstrated, simulated or implied by this repository.",
            "- Animal and in-vitro results never produce human event credit.",
            "- Active-comparator and Mendelian-randomisation contrasts are never applied",
            "  as placebo-relative treatment effects.",
        ]
    )
    return "\n".join(lines) + "\n"
