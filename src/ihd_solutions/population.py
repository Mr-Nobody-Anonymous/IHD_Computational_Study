"""Scenario populations: baseline risk and baseline biomarkers, not real patients."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from ihd_solutions.evidence import DEFAULT_DATA_DIR
from ihd_solutions.units import Z_95, require_positive

POPULATION_COLUMNS = (
    "population_id",
    "label",
    "baseline_annual_mace_rate",
    "baseline_rate_ci_low",
    "baseline_rate_ci_high",
    "ldl_c_mg_dl",
    "apob_mg_dl",
    "hscrp_mg_l",
    "systolic_bp_mmhg",
    "lp_a_nmol_l",
    "has_ascvd",
    "has_diabetes",
    "has_ckd",
    "has_heart_failure",
    "is_secondary_prevention",
    "notes",
)


@dataclass(frozen=True)
class Population:
    population_id: str
    label: str
    baseline_annual_mace_rate: float
    baseline_rate_ci_low: float
    baseline_rate_ci_high: float
    ldl_c_mg_dl: float
    apob_mg_dl: float
    hscrp_mg_l: float
    systolic_bp_mmhg: float
    lp_a_nmol_l: float
    has_ascvd: bool
    has_diabetes: bool
    has_ckd: bool
    has_heart_failure: bool
    is_secondary_prevention: bool
    notes: str = ""

    def __post_init__(self) -> None:
        if not self.population_id:
            raise ValueError("population_id is required")
        require_positive("baseline_annual_mace_rate", self.baseline_annual_mace_rate)
        if self.baseline_annual_mace_rate >= 1.0:
            raise ValueError("baseline annual hazard must stay below 1 for these scenarios")
        if not (
            0
            < self.baseline_rate_ci_low
            <= self.baseline_annual_mace_rate
            <= self.baseline_rate_ci_high
        ):
            raise ValueError(f"{self.population_id}: baseline rate interval must bracket the rate")
        for name in (
            "ldl_c_mg_dl",
            "apob_mg_dl",
            "hscrp_mg_l",
            "systolic_bp_mmhg",
            "lp_a_nmol_l",
        ):
            require_positive(name, getattr(self, name))

    @property
    def log_rate_sigma(self) -> float:
        return float(
            (math.log(self.baseline_rate_ci_high) - math.log(self.baseline_rate_ci_low))
            / (2.0 * Z_95)
        )


def _as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"true", "1", "yes"}:
        return True
    if text in {"false", "0", "no"}:
        return False
    raise ValueError(f"cannot interpret boolean value: {value!r}")


def load_populations(path: str | Path | None = None) -> dict[str, Population]:
    path = Path(path) if path is not None else DEFAULT_DATA_DIR / "populations.csv"
    if not path.is_file():
        raise FileNotFoundError(f"population file not found: {path}")
    frame = pd.read_csv(path)
    missing = [column for column in POPULATION_COLUMNS if column not in frame.columns]
    if missing:
        raise ValueError(f"population file missing columns: {missing}")
    populations: dict[str, Population] = {}
    for row in frame.to_dict(orient="records"):
        population = Population(
            population_id=str(row["population_id"]).strip(),
            label=str(row["label"]).strip(),
            baseline_annual_mace_rate=float(row["baseline_annual_mace_rate"]),
            baseline_rate_ci_low=float(row["baseline_rate_ci_low"]),
            baseline_rate_ci_high=float(row["baseline_rate_ci_high"]),
            ldl_c_mg_dl=float(row["ldl_c_mg_dl"]),
            apob_mg_dl=float(row["apob_mg_dl"]),
            hscrp_mg_l=float(row["hscrp_mg_l"]),
            systolic_bp_mmhg=float(row["systolic_bp_mmhg"]),
            lp_a_nmol_l=float(row["lp_a_nmol_l"]),
            has_ascvd=_as_bool(row["has_ascvd"]),
            has_diabetes=_as_bool(row["has_diabetes"]),
            has_ckd=_as_bool(row["has_ckd"]),
            has_heart_failure=_as_bool(row["has_heart_failure"]),
            is_secondary_prevention=_as_bool(row["is_secondary_prevention"]),
            notes=str(row.get("notes") or ""),
        )
        if population.population_id in populations:
            raise ValueError(f"duplicate population_id: {population.population_id}")
        populations[population.population_id] = population
    return populations


def background_regimen_ids() -> tuple[str, ...]:
    """Therapies treated as already-present standard care in every scenario."""
    return ("statin_high_intensity", "aspirin_secondary_prevention")
