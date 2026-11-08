"""The therapy catalogue: every candidate solution the simulator can combine."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path

import pandas as pd

from ihd_solutions.evidence import DEFAULT_DATA_DIR, EvidenceTier
from ihd_solutions.units import require_fraction, require_nonnegative, require_positive


class TherapyClass(str, Enum):
    LIPID_LOWERING = "LIPID_LOWERING"
    LPA_LOWERING = "LPA_LOWERING"
    ANTI_INFLAMMATORY = "ANTI_INFLAMMATORY"
    METABOLIC = "METABOLIC"
    BLOOD_PRESSURE = "BLOOD_PRESSURE"
    ANTITHROMBOTIC = "ANTITHROMBOTIC"
    LIFESTYLE = "LIFESTYLE"
    REVASCULARIZATION = "REVASCULARIZATION"
    REGENERATIVE = "REGENERATIVE"
    GENE_EDITING = "GENE_EDITING"


class DosingProfile(str, Enum):
    DAILY_ORAL = "DAILY_ORAL"
    INTERVAL_SC = "INTERVAL_SC"
    ONE_TIME_EDIT = "ONE_TIME_EDIT"
    ONE_TIME_PROCEDURE = "ONE_TIME_PROCEDURE"
    CONTINUOUS = "CONTINUOUS"


CATALOG_COLUMNS = (
    "therapy_id",
    "name",
    "therapy_class",
    "mechanism",
    "exclusive_group",
    "profile",
    "onset_days",
    "plateau_days",
    "interval_days",
    "trough_fraction",
    "waning_half_life_days",
    "ldl_pct",
    "ldl_pct_ci_low",
    "ldl_pct_ci_high",
    "apob_pct",
    "lpa_pct",
    "hscrp_pct",
    "sbp_reduction_mmhg",
    "mmp_factor",
    "burden_score",
    "evidence_tier",
    "regulatory_status",
    "requires_fasting",
    "symptom_benefit_only",
    "source_url",
    "notes",
)


@dataclass(frozen=True)
class Therapy:
    """A candidate intervention with mechanistic effects and dosing kinetics."""

    therapy_id: str
    name: str
    therapy_class: TherapyClass
    mechanism: str
    exclusive_group: str
    profile: DosingProfile
    onset_days: float
    plateau_days: float
    interval_days: float
    trough_fraction: float
    waning_half_life_days: float
    ldl_pct: float
    ldl_pct_ci_low: float
    ldl_pct_ci_high: float
    apob_pct: float
    lpa_pct: float
    hscrp_pct: float
    sbp_reduction_mmhg: float
    mmp_factor: float
    burden_score: float
    evidence_tier: EvidenceTier
    regulatory_status: str
    requires_fasting: bool
    symptom_benefit_only: bool
    source_url: str = ""
    notes: str = ""

    def __post_init__(self) -> None:
        if not self.therapy_id or not self.name:
            raise ValueError("therapy_id and name are required")
        require_nonnegative("onset_days", self.onset_days)
        require_positive("plateau_days", self.plateau_days)
        require_nonnegative("interval_days", self.interval_days)
        require_nonnegative("waning_half_life_days", self.waning_half_life_days)
        require_fraction("trough_fraction", self.trough_fraction, allow_one=True)
        if self.trough_fraction <= 0:
            raise ValueError(f"{self.therapy_id}: trough_fraction must be positive")
        for name in ("ldl_pct", "apob_pct", "lpa_pct", "hscrp_pct"):
            require_fraction(name, getattr(self, name))
        require_nonnegative("sbp_reduction_mmhg", self.sbp_reduction_mmhg)
        require_positive("mmp_factor", self.mmp_factor)
        if self.mmp_factor > 1.0:
            raise ValueError(f"{self.therapy_id}: mmp_factor above 1 would increase remodeling")
        require_positive("burden_score", self.burden_score)
        if self.plateau_days < self.onset_days:
            raise ValueError(f"{self.therapy_id}: plateau_days must be at least onset_days")
        if self.profile is DosingProfile.INTERVAL_SC and self.interval_days <= 0:
            raise ValueError(f"{self.therapy_id}: interval dosing requires interval_days > 0")
        if not (0 <= self.ldl_pct_ci_low <= self.ldl_pct <= self.ldl_pct_ci_high < 1):
            raise ValueError(f"{self.therapy_id}: LDL effect interval must bracket the point value")

    @property
    def is_permanent_scenario(self) -> bool:
        return self.profile is DosingProfile.ONE_TIME_EDIT and self.waning_half_life_days == 0

    @property
    def lowers_ldl(self) -> bool:
        return self.ldl_pct > 0

    @property
    def is_procedural(self) -> bool:
        return self.profile is DosingProfile.ONE_TIME_PROCEDURE


@dataclass(frozen=True)
class Regimen:
    """An ordered, duplicate-free combination of therapies."""

    therapy_ids: tuple[str, ...]
    label: str = ""

    def __post_init__(self) -> None:
        if not self.therapy_ids:
            raise ValueError("a regimen needs at least one therapy (use background for control)")
        if len(set(self.therapy_ids)) != len(self.therapy_ids):
            raise ValueError(f"duplicate therapies in regimen: {self.therapy_ids}")
        object.__setattr__(self, "therapy_ids", tuple(sorted(self.therapy_ids)))
        if not self.label:
            object.__setattr__(self, "label", " + ".join(self.therapy_ids))

    def __len__(self) -> int:
        return len(self.therapy_ids)


def _as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"true", "1", "yes"}:
        return True
    if text in {"false", "0", "no"}:
        return False
    raise ValueError(f"cannot interpret boolean value: {value!r}")


def load_catalog(path: str | Path | None = None) -> dict[str, Therapy]:
    """Load and validate the therapy catalogue."""
    path = Path(path) if path is not None else DEFAULT_DATA_DIR / "therapy_catalog.csv"
    if not path.is_file():
        raise FileNotFoundError(f"therapy catalogue not found: {path}")
    frame = pd.read_csv(path)
    missing = [column for column in CATALOG_COLUMNS if column not in frame.columns]
    if missing:
        raise ValueError(f"therapy catalogue missing columns: {missing}")
    catalog: dict[str, Therapy] = {}
    for row in frame.to_dict(orient="records"):
        therapy = Therapy(
            therapy_id=str(row["therapy_id"]).strip(),
            name=str(row["name"]).strip(),
            therapy_class=TherapyClass(str(row["therapy_class"]).strip()),
            mechanism=str(row["mechanism"]).strip(),
            exclusive_group=str(row["exclusive_group"]).strip(),
            profile=DosingProfile(str(row["profile"]).strip()),
            onset_days=float(row["onset_days"]),
            plateau_days=float(row["plateau_days"]),
            interval_days=float(row["interval_days"]),
            trough_fraction=float(row["trough_fraction"]),
            waning_half_life_days=float(row["waning_half_life_days"]),
            ldl_pct=float(row["ldl_pct"]),
            ldl_pct_ci_low=float(row["ldl_pct_ci_low"]),
            ldl_pct_ci_high=float(row["ldl_pct_ci_high"]),
            apob_pct=float(row["apob_pct"]),
            lpa_pct=float(row["lpa_pct"]),
            hscrp_pct=float(row["hscrp_pct"]),
            sbp_reduction_mmhg=float(row["sbp_reduction_mmhg"]),
            mmp_factor=float(row["mmp_factor"]),
            burden_score=float(row["burden_score"]),
            evidence_tier=EvidenceTier(str(row["evidence_tier"]).strip()),
            regulatory_status=str(row["regulatory_status"]).strip(),
            requires_fasting=_as_bool(row["requires_fasting"]),
            symptom_benefit_only=_as_bool(row["symptom_benefit_only"]),
            source_url=str(row.get("source_url") or ""),
            notes=str(row.get("notes") or ""),
        )
        if therapy.therapy_id in catalog:
            raise ValueError(f"duplicate therapy_id: {therapy.therapy_id}")
        catalog[therapy.therapy_id] = therapy
    if not catalog:
        raise ValueError("therapy catalogue is empty")
    return catalog


def therapies_in_class(catalog: dict[str, Therapy], therapy_class: TherapyClass) -> list[Therapy]:
    return sorted(
        (t for t in catalog.values() if t.therapy_class is therapy_class),
        key=lambda t: t.therapy_id,
    )


def validate_combination(catalog: dict[str, Therapy], therapy_ids: tuple[str, ...]) -> Regimen:
    """Reject unknown ids and same-target stacking (e.g. two PCSK9 agents)."""
    unknown = [tid for tid in therapy_ids if tid not in catalog]
    if unknown:
        raise ValueError(f"unknown therapy ids: {unknown}")
    groups: dict[str, str] = {}
    for therapy_id in sorted(therapy_ids):
        group = catalog[therapy_id].exclusive_group
        if group in groups:
            raise ValueError(
                f"therapies {groups[group]} and {therapy_id} share exclusive group {group!r}"
            )
        groups[group] = therapy_id
    return Regimen(tuple(therapy_ids))


def catalog_frame(catalog: dict[str, Therapy]) -> pd.DataFrame:
    rows = []
    for therapy in sorted(catalog.values(), key=lambda t: (t.therapy_class.value, t.therapy_id)):
        rows.append(
            {
                "therapy_id": therapy.therapy_id,
                "name": therapy.name,
                "therapy_class": therapy.therapy_class.value,
                "profile": therapy.profile.value,
                "exclusive_group": therapy.exclusive_group,
                "ldl_pct": therapy.ldl_pct,
                "lpa_pct": therapy.lpa_pct,
                "hscrp_pct": therapy.hscrp_pct,
                "sbp_reduction_mmhg": therapy.sbp_reduction_mmhg,
                "burden_score": therapy.burden_score,
                "evidence_tier": therapy.evidence_tier.value,
                "regulatory_status": therapy.regulatory_status,
                "requires_fasting": therapy.requires_fasting,
                "symptom_benefit_only": therapy.symptom_benefit_only,
                "source_url": therapy.source_url,
            }
        )
    return pd.DataFrame(rows)
