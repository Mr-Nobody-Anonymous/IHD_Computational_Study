"""Base-editing scenarios: pharmacodynamics yes, safety scoring deliberately no."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import pandas as pd

from ihd_solutions.units import DAYS_PER_YEAR, require_fraction, require_positive

#: Heart-2 anchor: about 88% time-averaged PCSK9 reduction produced about 62%
#: time-averaged LDL-C reduction at 1.0 mg/kg (phase 1, n=35).
HEART2_PCSK9_KNOCKDOWN = 0.88
HEART2_LDL_REDUCTION = 0.62
MAX_LDL_REDUCTION = 0.75
DEFAULT_WINDOW = (4, 8)
_DNA = frozenset("ACGT")


def calibrate_sharpness(
    anchor_knockdown: float = HEART2_PCSK9_KNOCKDOWN,
    anchor_reduction: float = HEART2_LDL_REDUCTION,
    max_reduction: float = MAX_LDL_REDUCTION,
) -> float:
    """Solve the saturating curve so it reproduces the published anchor pair."""
    require_fraction("anchor_knockdown", anchor_knockdown)
    require_positive("anchor_knockdown", anchor_knockdown)
    require_fraction("anchor_reduction", anchor_reduction)
    require_fraction("max_reduction", max_reduction)
    if anchor_reduction >= max_reduction:
        raise ValueError("anchor reduction must stay below the saturation ceiling")
    return float(-math.log(1.0 - anchor_reduction / max_reduction) / anchor_knockdown)


def ldl_reduction_from_knockdown(
    knockdown: np.ndarray | float,
    sharpness: float | None = None,
    max_reduction: float = MAX_LDL_REDUCTION,
) -> np.ndarray:
    """Saturating map from PCSK9 protein knockdown to LDL-C reduction."""
    knockdown = np.asarray(knockdown, dtype=float)
    if not np.isfinite(knockdown).all():
        raise ValueError("knockdown must be finite")
    if np.any(knockdown < 0) or np.any(knockdown > 1):
        raise ValueError("knockdown must lie in [0, 1]")
    require_fraction("max_reduction", max_reduction)
    sharpness = calibrate_sharpness(max_reduction=max_reduction) if sharpness is None else sharpness
    require_positive("sharpness", sharpness)
    return max_reduction * (-np.expm1(-sharpness * knockdown))


@dataclass(frozen=True)
class EditingScenario:
    """One-time editing scenario with explicit durability assumptions."""

    scenario_id: str
    editor_version: str = "ABE8e"
    edited_hepatocyte_fraction: float = 0.90
    knockdown_ceiling: float = 0.95
    onset_days: float = 14.0
    plateau_days: float = 28.0
    durability_half_life_years: float | None = None
    redosing_allowed: bool = False

    def __post_init__(self) -> None:
        if not self.scenario_id:
            raise ValueError("scenario_id is required")
        require_fraction(
            "edited_hepatocyte_fraction", self.edited_hepatocyte_fraction, allow_one=True
        )
        require_fraction("knockdown_ceiling", self.knockdown_ceiling, allow_one=True)
        require_positive("plateau_days", self.plateau_days)
        if self.onset_days < 0 or self.onset_days > self.plateau_days:
            raise ValueError("onset_days must lie in [0, plateau_days]")
        if self.durability_half_life_years is not None:
            require_positive("durability_half_life_years", self.durability_half_life_years)

    @property
    def is_permanent(self) -> bool:
        return self.durability_half_life_years is None

    @property
    def plateau_knockdown(self) -> float:
        return float(self.knockdown_ceiling * self.edited_hepatocyte_fraction)


def simulate_editing(
    scenario: EditingScenario,
    horizon_years: float = 10.0,
    step_days: float = 30.0,
    baseline_ldl_mg_dl: float = 130.0,
) -> pd.DataFrame:
    """Project knockdown, LDL-C and durability for a single-infusion scenario."""
    require_positive("horizon_years", horizon_years)
    require_positive("step_days", step_days)
    require_positive("baseline_ldl_mg_dl", baseline_ldl_mg_dl)
    total_days = horizon_years * DAYS_PER_YEAR
    n_steps = max(1, int(math.ceil(total_days / step_days)))
    days = np.linspace(0.0, total_days, n_steps + 1)
    elapsed = np.clip(days - scenario.onset_days, 0.0, None)
    span = max(scenario.plateau_days - scenario.onset_days, 1e-9)
    ramp = np.where(days < scenario.onset_days, 0.0, np.clip(elapsed / span, 0.0, 1.0))
    if scenario.is_permanent:
        durability = np.ones_like(days)
    else:
        decay_days = np.maximum(days - scenario.plateau_days, 0.0)
        durability = np.exp(
            -np.log(2.0) * decay_days / (scenario.durability_half_life_years * DAYS_PER_YEAR)
        )
    knockdown = scenario.plateau_knockdown * ramp * durability
    reduction = ldl_reduction_from_knockdown(knockdown)
    ldl = np.maximum(baseline_ldl_mg_dl * (1.0 - reduction), 10.0)
    return pd.DataFrame(
        {
            "scenario_id": scenario.scenario_id,
            "editor_version": scenario.editor_version,
            "day": days,
            "year": days / DAYS_PER_YEAR,
            "pcsk9_knockdown_pct": 100.0 * knockdown,
            "ldl_reduction_pct": 100.0 * reduction,
            "ldl_c_mg_dl": ldl,
            "durability_assumption": "permanent_scenario" if scenario.is_permanent else "waning",
            "clinical_safety_assessed": False,
        }
    )


def validate_guide(spacer: str, pam: str) -> tuple[str, str]:
    """Strict sequence validation; ambiguous bases and wrong lengths are rejected."""
    spacer = spacer.strip().upper()
    pam = pam.strip().upper()
    if len(spacer) != 20:
        raise ValueError(f"spacer must be 20 nucleotides, got {len(spacer)}")
    if len(pam) != 3:
        raise ValueError(f"PAM must be 3 nucleotides, got {len(pam)}")
    if not set(spacer) <= _DNA:
        raise ValueError("spacer contains ambiguous or non-DNA characters")
    if not set(pam) <= _DNA:
        raise ValueError("PAM contains ambiguous or non-DNA characters")
    if pam[1:] != "GG":
        raise ValueError(f"SpCas9-family PAM must be NGG, got {pam}")
    return spacer, pam


def bystander_report(
    spacer: str,
    pam: str,
    window: tuple[int, int] = DEFAULT_WINDOW,
    intended_position: int | None = None,
) -> dict[str, object]:
    """Count adenines inside the editing window (deterministic sequence analysis)."""
    spacer, pam = validate_guide(spacer, pam)
    start, end = window
    if not (1 <= start <= end <= 20):
        raise ValueError("editing window must lie within protospacer positions 1-20")
    adenines = [index + 1 for index, base in enumerate(spacer) if base == "A"]
    in_window = [position for position in adenines if start <= position <= end]
    if intended_position is not None:
        if intended_position not in adenines:
            raise ValueError(
                f"intended position {intended_position} is not an adenine in this spacer"
            )
        bystanders = [position for position in in_window if position != intended_position]
    else:
        bystanders = in_window
    return {
        "spacer": spacer,
        "pam": pam,
        "window_start": start,
        "window_end": end,
        "adenine_positions": adenines,
        "adenines_in_window": in_window,
        "intended_position": intended_position,
        "bystander_positions": bystanders,
        "n_bystanders": len(bystanders),
        "interpretation": "sequence-level bystander risk only; not an editing-outcome prediction",
    }


def off_target_assessment() -> dict[str, object]:
    """Refuse to emit a numeric safety score: base editors need dedicated assays."""
    return {
        "status": "requires_experimental_assay",
        "numeric_off_target_score": None,
        "reason": (
            "SpCas9 cleavage-based proxies and double-strand-break assays do not measure "
            "adenine-base-editor off-target editing; no numeric safety score is produced"
        ),
        "required_assays": required_safety_assays(),
    }


def required_safety_assays() -> tuple[str, ...]:
    return (
        "base-editor-specific genome-wide DNA off-target assay",
        "RNA off-target transcriptome-wide assay",
        "bystander editing amplicon sequencing at the target locus",
        "indel and structural-variant spectrum at the target locus",
        "germline transmission assessment in an appropriate model",
        "long-term clinical follow-up with malignancy and hepatic monitoring",
    )
