"""Evidence tiers, trial anchors, and the rules that decide whether benefit exists.

This module is the safety interlock of the whole simulator. A therapy can lower a
biomarker by any amount and still receive **zero** event-reduction credit: credit
requires an anchor with a hard endpoint, an acceptable evidence tier for the
requested :class:`CreditMode`, and a population match (or an explicitly declared
transportable population).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path

import pandas as pd

from ihd_solutions.units import Z_95, require_positive

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA_DIR = REPO_ROOT / "data" / "solutions"


class EvidenceTier(str, Enum):
    """Strength of the evidence behind a *clinical endpoint*, not a biomarker."""

    META_ANALYSIS = "META_ANALYSIS"
    RANDOMIZED_OUTCOMES = "RANDOMIZED_OUTCOMES"
    RANDOMIZED_SURROGATE = "RANDOMIZED_SURROGATE"
    #: Randomised hard-endpoint result announced as topline only (no HR/CI yet).
    RANDOMIZED_TOPLINE = "RANDOMIZED_TOPLINE"
    #: Randomised against another active drug, so the HR is NOT placebo-relative.
    RANDOMIZED_ACTIVE_COMPARATOR = "RANDOMIZED_ACTIVE_COMPARATOR"
    #: Lifelong genetic exposure contrast, not the effect of starting a treatment.
    MENDELIAN_RANDOMIZATION = "MENDELIAN_RANDOMIZATION"
    OBSERVATIONAL = "OBSERVATIONAL"
    CONDITIONAL_APPROVAL = "CONDITIONAL_APPROVAL"
    MECHANISTIC = "MECHANISTIC"


class EndpointClass(str, Enum):
    """Endpoint families that must never be silently interchanged."""

    MACE = "MACE"
    CV_DEATH = "CV_DEATH"
    ALL_CAUSE_DEATH = "ALL_CAUSE_DEATH"
    HF_COMPOSITE = "HF_COMPOSITE"
    SURROGATE = "SURROGATE"
    FUNCTIONAL = "FUNCTIONAL"
    PENDING = "PENDING"
    PLAQUE_IMAGING = "PLAQUE_IMAGING"


class EffectUnit(str, Enum):
    """How an anchor's hazard ratio must be applied."""

    FIXED = "fixed"
    PER_MMOL_L_LDL = "per_mmol_l_ldl"
    PER_5_MMHG_SBP = "per_5_mmhg_sbp"


class CreditMode(str, Enum):
    """Which evidence classes are allowed to move the event curve."""

    STRICT = "strict"
    SURROGATE = "surrogate"
    EXPLORATORY = "exploratory"


#: Tiers whose hard-endpoint anchors may be used directly.
_RANDOMISED_HARD: frozenset[EvidenceTier] = frozenset(
    {
        EvidenceTier.META_ANALYSIS,
        EvidenceTier.RANDOMIZED_OUTCOMES,
        EvidenceTier.RANDOMIZED_TOPLINE,
    }
)
_DIRECT_TIERS: dict[CreditMode, frozenset[EvidenceTier]] = {
    CreditMode.STRICT: _RANDOMISED_HARD,
    CreditMode.SURROGATE: _RANDOMISED_HARD,
    CreditMode.EXPLORATORY: _RANDOMISED_HARD | {EvidenceTier.OBSERVATIONAL},
}

#: Evidence that describes something other than "start this drug, avoid events".
#: A Mendelian-randomisation contrast is a lifetime genetic exposure difference;
#: an active-comparator HR already contains the comparator's own benefit. Neither
#: may ever be applied as a treatment effect, in any mode.
_NON_TREATMENT_TIERS = frozenset(
    {EvidenceTier.MENDELIAN_RANDOMIZATION, EvidenceTier.RANDOMIZED_ACTIVE_COMPARATOR}
)

#: Tiers that may carry a biomarker effect into a transported per-unit pathway.
#: Animal/in-vitro mechanisms are excluded: mouse data never becomes human benefit.
_HUMAN_BIOMARKER_TIERS = frozenset(set(EvidenceTier) - {EvidenceTier.MECHANISTIC})


def is_treatment_evidence(tier: EvidenceTier) -> bool:
    """True when a tier can describe the effect of starting a therapy."""
    return tier not in _NON_TREATMENT_TIERS


def may_transport_biomarker(tier: EvidenceTier) -> bool:
    """True when human data exist for this therapy's biomarker effect."""
    return tier in _HUMAN_BIOMARKER_TIERS


#: Endpoints that may never be relabelled as MACE, in any mode.
_NEVER_MACE = frozenset(
    {
        EndpointClass.SURROGATE,
        EndpointClass.FUNCTIONAL,
        EndpointClass.PENDING,
        EndpointClass.HF_COMPOSITE,
        EndpointClass.PLAQUE_IMAGING,
    }
)

ANCHOR_COLUMNS = (
    "anchor_id",
    "therapy_id",
    "population_id",
    "endpoint",
    "effect_unit",
    "endpoint_class",
    "hr",
    "ci_low",
    "ci_high",
    "followup_years",
    "n_participants",
    "evidence_tier",
    "transportable_to",
    "source_url",
    "notes",
)


@dataclass(frozen=True)
class TrialAnchor:
    """One published effect estimate, scoped to a population and an endpoint."""

    anchor_id: str
    therapy_id: str
    population_id: str
    endpoint: str
    effect_unit: EffectUnit
    endpoint_class: EndpointClass
    hr: float | None
    ci_low: float | None
    ci_high: float | None
    followup_years: float
    n_participants: int
    evidence_tier: EvidenceTier
    transportable_to: tuple[str, ...] = ()
    source_url: str = ""
    notes: str = ""
    #: Neutrality of the *published* interval. Resampling rewrites the CI bounds,
    #: so the reported verdict is stored separately and survives Monte Carlo.
    neutral_by_report: bool | None = None

    def __post_init__(self) -> None:
        if not self.anchor_id or not self.therapy_id or not self.population_id:
            raise ValueError("anchor_id, therapy_id and population_id are required")
        require_positive("followup_years", self.followup_years)
        if self.n_participants <= 0:
            raise ValueError(f"{self.anchor_id}: n_participants must be positive")
        if self.hr is None:
            if self.ci_low is not None or self.ci_high is not None:
                raise ValueError(f"{self.anchor_id}: CI without a point estimate")
            if self.endpoint_class in {EndpointClass.MACE, EndpointClass.CV_DEATH}:
                raise ValueError(f"{self.anchor_id}: hard endpoint requires a hazard ratio")
            return
        require_positive("hr", self.hr)
        if self.ci_low is None or self.ci_high is None:
            raise ValueError(f"{self.anchor_id}: hazard ratio requires both CI bounds")
        if not (0 < self.ci_low <= self.hr <= self.ci_high):
            raise ValueError(f"{self.anchor_id}: require 0 < ci_low <= hr <= ci_high")

    @property
    def has_effect_estimate(self) -> bool:
        return self.hr is not None

    @property
    def is_neutral(self) -> bool:
        """True when the confidence interval spans the null value of 1.0."""
        if self.hr is None:
            return False
        return bool(self.ci_low <= 1.0 <= self.ci_high)

    @property
    def neutral_at_report(self) -> bool:
        """True when the published interval spans 1.0, regardless of resampling."""
        if self.neutral_by_report is not None:
            return bool(self.neutral_by_report)
        return self.is_neutral

    @property
    def log_sigma(self) -> float:
        """Standard deviation of log(HR) implied by the reported 95% interval."""
        if self.hr is None:
            raise ValueError(f"{self.anchor_id}: no effect estimate to derive uncertainty from")
        return float((math.log(self.ci_high) - math.log(self.ci_low)) / (2.0 * Z_95))

    def applies_to(self, population_id: str, mode: CreditMode) -> bool:
        if population_id == self.population_id:
            return True
        if population_id in self.transportable_to:
            return True
        return mode is CreditMode.EXPLORATORY

    def scaled_hr(self, magnitude: float = 1.0) -> float:
        """Apply a per-unit anchor to ``magnitude`` units, or return a fixed HR."""
        if self.hr is None:
            raise ValueError(f"{self.anchor_id}: anchor has no effect estimate")
        if magnitude < 0:
            raise ValueError("magnitude must be non-negative")
        if self.effect_unit is EffectUnit.FIXED:
            return float(self.hr)
        return float(self.hr**magnitude)


@dataclass(frozen=True)
class CreditDecision:
    """Auditable explanation of why an effect was or was not credited."""

    therapy_id: str
    population_id: str
    mode: CreditMode
    credited: bool
    reason: str
    anchor_id: str | None = None
    hr: float | None = None
    pathway: str = "none"
    warnings: tuple[str, ...] = field(default=())


def _reported_neutrality(row: dict) -> bool:
    hr = _optional_float(row["hr"])
    if hr is None:
        return False
    low = _optional_float(row["ci_low"])
    high = _optional_float(row["ci_high"])
    return bool(low is not None and high is not None and low <= 1.0 <= high)


def load_anchors(path: str | Path | None = None) -> dict[str, TrialAnchor]:
    """Load and validate the anchor registry, keyed by ``anchor_id``."""
    path = Path(path) if path is not None else DEFAULT_DATA_DIR / "trial_anchors.csv"
    if not path.is_file():
        raise FileNotFoundError(f"anchor registry not found: {path}")
    frame = pd.read_csv(path)
    missing = [column for column in ANCHOR_COLUMNS if column not in frame.columns]
    if missing:
        raise ValueError(f"anchor registry missing columns: {missing}")
    anchors: dict[str, TrialAnchor] = {}
    seen_pairs: set[tuple[str, str, str]] = set()
    for row in frame.to_dict(orient="records"):
        transportable = tuple(
            item.strip()
            for item in str(row.get("transportable_to") or "").split(";")
            if item.strip() and item.strip().lower() != "nan"
        )
        anchor = TrialAnchor(
            anchor_id=str(row["anchor_id"]).strip(),
            therapy_id=str(row["therapy_id"]).strip(),
            population_id=str(row["population_id"]).strip(),
            endpoint=str(row["endpoint"]).strip(),
            effect_unit=EffectUnit(str(row["effect_unit"]).strip()),
            endpoint_class=EndpointClass(str(row["endpoint_class"]).strip()),
            hr=_optional_float(row["hr"]),
            ci_low=_optional_float(row["ci_low"]),
            ci_high=_optional_float(row["ci_high"]),
            followup_years=float(row["followup_years"]),
            n_participants=int(row["n_participants"]),
            evidence_tier=EvidenceTier(str(row["evidence_tier"]).strip()),
            transportable_to=transportable,
            neutral_by_report=_reported_neutrality(row),
            source_url=str(row.get("source_url") or ""),
            notes=_clean_text(row.get("notes")),
        )
        if anchor.anchor_id in anchors:
            raise ValueError(f"duplicate anchor_id: {anchor.anchor_id}")
        # One anchor per therapy/population/endpoint-family: a therapy may have
        # both an outcomes trial and an imaging trial in the same population, but
        # never two competing estimates of the same endpoint family.
        pair = (anchor.therapy_id, anchor.population_id, anchor.endpoint_class.value)
        if pair in seen_pairs:
            raise ValueError(f"duplicate therapy/population/endpoint anchor: {pair}")
        seen_pairs.add(pair)
        anchors[anchor.anchor_id] = anchor
    if not anchors:
        raise ValueError("anchor registry is empty")
    return anchors


def _optional_float(value: object) -> float | None:
    if value is None:
        return None
    if isinstance(value, str) and not value.strip():
        return None
    number = float(value)  # type: ignore[arg-type]
    return None if math.isnan(number) else number


def _clean_text(value: object) -> str:
    text = "" if value is None else str(value)
    return "" if text.strip().lower() == "nan" else text


def anchors_for_therapy(
    anchors: dict[str, TrialAnchor], therapy_id: str
) -> tuple[TrialAnchor, ...]:
    return tuple(
        sorted(
            (a for a in anchors.values() if a.therapy_id == therapy_id),
            key=lambda a: (a.population_id, a.anchor_id),
        )
    )


def select_anchor(
    anchors: dict[str, TrialAnchor],
    therapy_id: str,
    population_id: str,
    mode: CreditMode,
) -> TrialAnchor | None:
    """Prefer an exact population match, then a declared transportable match."""
    # Plaque-imaging anchors describe structure, not events. They feed the
    # regression model and must never enter event-credit selection, where they
    # would shadow a therapy's real outcomes anchor.
    candidates = [
        a
        for a in anchors_for_therapy(anchors, therapy_id)
        if a.endpoint_class is not EndpointClass.PLAQUE_IMAGING
    ]

    def preference(anchor: TrialAnchor) -> tuple[int, int, int]:
        return (
            0 if is_treatment_evidence(anchor.evidence_tier) else 1,
            0 if anchor.evidence_tier in _DIRECT_TIERS[mode] else 1,
            0 if anchor.endpoint_class is EndpointClass.MACE else 1,
        )

    def best(subset: list[TrialAnchor]) -> TrialAnchor | None:
        return min(subset, key=preference) if subset else None

    exact = [a for a in candidates if a.population_id == population_id]
    if exact:
        return best(exact)
    transportable = [a for a in candidates if population_id in a.transportable_to]
    if transportable:
        return best(transportable)
    if mode is CreditMode.EXPLORATORY:
        hard = [
            a
            for a in candidates
            if a.endpoint_class is EndpointClass.MACE and is_treatment_evidence(a.evidence_tier)
        ]
        if hard:
            return best(hard)
        return best(candidates)
    return None


def decide_credit(
    therapy_id: str,
    population_id: str,
    mode: CreditMode,
    anchors: dict[str, TrialAnchor],
    *,
    biomarker_only_change: float = 0.0,
) -> CreditDecision:
    """Decide whether a therapy may reduce events, and record the reason.

    ``biomarker_only_change`` is accepted purely so the refusal can be documented:
    a biomarker change never creates credit on its own (the ZEUS/CIRT lesson).
    """
    anchor = select_anchor(anchors, therapy_id, population_id, mode)
    if anchor is None:
        reason = "no trial anchor available for this therapy"
        if biomarker_only_change > 0:
            reason += "; biomarker change alone never credits events"
        return CreditDecision(
            therapy_id, population_id, mode, False, reason, warnings=("no_anchor",)
        )
    if anchor.evidence_tier is EvidenceTier.MENDELIAN_RANDOMIZATION:
        return CreditDecision(
            therapy_id,
            population_id,
            mode,
            False,
            f"{anchor.anchor_id} is a lifelong genetic exposure contrast, not a treatment effect",
            anchor.anchor_id,
            float(anchor.hr) if anchor.hr is not None else None,
            warnings=("mendelian_not_a_treatment_effect",),
        )
    if anchor.evidence_tier is EvidenceTier.RANDOMIZED_ACTIVE_COMPARATOR:
        return CreditDecision(
            therapy_id,
            population_id,
            mode,
            False,
            f"{anchor.anchor_id} uses an active comparator; its HR is not placebo-relative",
            anchor.anchor_id,
            float(anchor.hr) if anchor.hr is not None else None,
            warnings=("active_comparator_not_placebo_relative",),
        )
    if not anchor.has_effect_estimate:
        return CreditDecision(
            therapy_id,
            population_id,
            mode,
            False,
            f"{anchor.anchor_id} reports no effect estimate ({anchor.endpoint_class.value})",
            anchor.anchor_id,
            warnings=("no_effect_estimate",),
        )
    if anchor.endpoint_class in _NEVER_MACE:
        return CreditDecision(
            therapy_id,
            population_id,
            mode,
            False,
            f"{anchor.anchor_id} endpoint class {anchor.endpoint_class.value} "
            "cannot be credited as MACE",
            anchor.anchor_id,
            warnings=("endpoint_class_not_mace",),
        )
    if anchor.endpoint_class is not EndpointClass.MACE and mode is not CreditMode.EXPLORATORY:
        return CreditDecision(
            therapy_id,
            population_id,
            mode,
            False,
            f"{anchor.anchor_id} endpoint {anchor.endpoint_class.value} is not a MACE composite; "
            "credit only in exploratory mode",
            anchor.anchor_id,
            warnings=("non_mace_endpoint",),
        )
    if anchor.evidence_tier not in _DIRECT_TIERS[mode]:
        return CreditDecision(
            therapy_id,
            population_id,
            mode,
            False,
            f"{anchor.anchor_id} tier {anchor.evidence_tier.value} not creditable in "
            f"{mode.value} mode",
            anchor.anchor_id,
            warnings=("tier_not_creditable",),
        )
    if not anchor.applies_to(population_id, mode):
        return CreditDecision(
            therapy_id,
            population_id,
            mode,
            False,
            f"{anchor.anchor_id} population {anchor.population_id} not transportable to "
            f"{population_id}",
            anchor.anchor_id,
            warnings=("population_mismatch",),
        )
    warnings: list[str] = []
    if anchor.is_neutral:
        warnings.append("neutral_confidence_interval")
    if anchor.evidence_tier is EvidenceTier.RANDOMIZED_TOPLINE:
        warnings.append("topline_only_no_published_interval")
    if anchor.population_id != population_id:
        warnings.append("transported_population")
    if anchor.evidence_tier is EvidenceTier.OBSERVATIONAL:
        warnings.append("observational_confounding")
    return CreditDecision(
        therapy_id,
        population_id,
        mode,
        True,
        f"{anchor.anchor_id} credited via {anchor.effect_unit.value} on {anchor.endpoint_class.value}",
        anchor.anchor_id,
        float(anchor.hr),
        pathway="anchor",
        warnings=tuple(warnings),
    )


def assert_no_biomarker_surrogate_credit(decisions: list[CreditDecision]) -> None:
    """Guard used by tests and the pipeline: refusals must stay refusals."""
    for decision in decisions:
        if decision.credited and decision.pathway == "biomarker":
            raise AssertionError(
                f"{decision.therapy_id}: biomarker-only credit is forbidden by design"
            )


def summarise_tiers(anchors: dict[str, TrialAnchor]) -> pd.DataFrame:
    rows = [
        {
            "evidence_tier": tier.value,
            "n_anchors": sum(1 for a in anchors.values() if a.evidence_tier is tier),
            "n_neutral": sum(
                1 for a in anchors.values() if a.evidence_tier is tier and a.is_neutral
            ),
        }
        for tier in EvidenceTier
    ]
    return pd.DataFrame(rows)
