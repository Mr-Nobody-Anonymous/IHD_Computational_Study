"""Atherosclerosis regression: how far backwards human plaque actually goes.

This is the closest available quantitative handle on "reversing" coronary disease,
and it is deliberately bounded. Serial-imaging trials show that percent atheroma
volume falls when achieved LDL-C is low enough, and that the *noncalcified*
(lipid-rich, rupture-prone) component shrinks fastest. They also show that
calcified volume does not fall - it rises, as soft plaque stabilises. So the model
has two compartments and can never report eradication.

Calibration (all human, all published):

* REVERSAL      pravastatin 40 mg, LDL-C 110 mg/dL: PAV 40.0 -> 41.8 over 18 months
* statin alone  LDL-C  93 mg/dL: PAV +0.05 pp over 24 months (essentially held)
* ASTEROID      rosuvastatin 40 mg, LDL-C 60.8 mg/dL: PAV -0.98 pp over 24 months
* GLAGOV        + evolocumab, LDL-C 36.5 mg/dL: PAV -0.95 pp over 18 months
* FLOWPROMOTE   statin-naive, 18 months, whole-heart CT:
                LDL-C 42.5 mg/dL -> noncalcified -28.2%, total -21.1%
                LDL-C 69.5 mg/dL -> noncalcified -20.0%, total -13.6%

Cross-trial pooling is a real limitation: these are different populations,
imaging modalities and durations. The fit quality is reported rather than hidden,
and any projection beyond the longest calibration follow-up is flagged.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from ihd_solutions.catalog import Regimen, Therapy, validate_combination
from ihd_solutions.evidence import CreditMode, EndpointClass, TrialAnchor
from ihd_solutions.lipids import simulate_biomarkers
from ihd_solutions.population import Population
from ihd_solutions.units import DAYS_PER_YEAR, require_positive

#: Longest follow-up among the calibration trials, in years.
CALIBRATION_HORIZON_YEARS = 2.0

#: Fraction of lost noncalcified volume assumed to reappear as dense calcium.
#: FLOWPROMOTE reports calcified volume rising in both arms but gives no transfer
#: coefficient, so this is an explicit modelling assumption, not a measurement.
DEFAULT_CALCIFICATION_FRACTION = 0.25


@dataclass(frozen=True)
class ImagingPoint:
    """One published (achieved LDL-C, percent atheroma volume change) pair."""

    trial: str
    ldl_mg_dl: float
    pav_change_pp: float
    years: float

    def __post_init__(self) -> None:
        require_positive("ldl_mg_dl", self.ldl_mg_dl)
        require_positive("years", self.years)

    @property
    def annual_rate_pp(self) -> float:
        return self.pav_change_pp / self.years


PAV_POINTS: tuple[ImagingPoint, ...] = (
    ImagingPoint("REVERSAL-pravastatin", 110.0, 1.8, 1.5),
    ImagingPoint("statin-moderate-achieved-93", 93.0, 0.05, 2.0),
    ImagingPoint("ASTEROID", 60.8, -0.98, 2.0),
    ImagingPoint("GLAGOV", 36.5, -0.95, 1.5),
)

#: (achieved LDL-C, fractional noncalcified volume loss, years) from FLOWPROMOTE.
NONCALCIFIED_POINTS: tuple[tuple[float, float, float], ...] = (
    (42.5, 0.282, 1.5),
    (69.5, 0.200, 1.5),
)


@dataclass(frozen=True)
class RegressionCalibration:
    """Linear fits plus the achieved-LDL threshold where plaque stops growing."""

    pav_slope: float
    pav_intercept: float
    pav_r_squared: float
    nc_decay_slope: float
    nc_decay_intercept: float
    threshold_ldl_mg_dl: float

    def annual_pav_change(self, ldl_mg_dl: float) -> float:
        """Percent-atheroma-volume change per year at a sustained achieved LDL-C."""
        require_positive("ldl_mg_dl", ldl_mg_dl)
        return self.pav_intercept + self.pav_slope * ldl_mg_dl

    def noncalcified_decay_rate(self, ldl_mg_dl: float) -> float:
        """Exponential loss rate of the lipid-rich compartment, per year."""
        require_positive("ldl_mg_dl", ldl_mg_dl)
        return max(0.0, self.nc_decay_intercept + self.nc_decay_slope * ldl_mg_dl)


def fit_calibration() -> RegressionCalibration:
    """Least-squares fit of the published imaging points (documented, not tuned)."""
    x = np.array([point.ldl_mg_dl for point in PAV_POINTS], dtype=float)
    y = np.array([point.annual_rate_pp for point in PAV_POINTS], dtype=float)
    slope, intercept = np.polyfit(x, y, 1)
    predicted = intercept + slope * x
    residual = float(np.sum((y - predicted) ** 2))
    total = float(np.sum((y - float(np.mean(y))) ** 2))
    r_squared = 1.0 - residual / total if total > 0 else 1.0

    nc_x = np.array([ldl for ldl, _, _ in NONCALCIFIED_POINTS], dtype=float)
    nc_y = np.array(
        [-math.log(1.0 - loss) / years for _, loss, years in NONCALCIFIED_POINTS], dtype=float
    )
    nc_slope, nc_intercept = np.polyfit(nc_x, nc_y, 1)
    if slope <= 0:
        raise ValueError("calibration produced a non-positive LDL slope")
    return RegressionCalibration(
        pav_slope=float(slope),
        pav_intercept=float(intercept),
        pav_r_squared=float(r_squared),
        nc_decay_slope=float(nc_slope),
        nc_decay_intercept=float(nc_intercept),
        threshold_ldl_mg_dl=float(-intercept / slope),
    )


DEFAULT_CALIBRATION = fit_calibration()


def calibration_frame(calibration: RegressionCalibration = DEFAULT_CALIBRATION) -> pd.DataFrame:
    """Expose the calibration inputs and fit so a reader can audit them."""
    rows = [
        {
            "trial": point.trial,
            "achieved_ldl_mg_dl": point.ldl_mg_dl,
            "pav_change_pp": point.pav_change_pp,
            "followup_years": point.years,
            "observed_pp_per_year": point.annual_rate_pp,
            "fitted_pp_per_year": calibration.annual_pav_change(point.ldl_mg_dl),
        }
        for point in PAV_POINTS
    ]
    frame = pd.DataFrame(rows)
    frame["residual_pp_per_year"] = frame["observed_pp_per_year"] - frame["fitted_pp_per_year"]
    return frame


def calibration_summary(calibration: RegressionCalibration = DEFAULT_CALIBRATION) -> dict:
    """The fitted parameters, including the LDL-C threshold for net regression."""
    return {
        "pav_slope_pp_per_year_per_mg_dl": calibration.pav_slope,
        "pav_intercept_pp_per_year": calibration.pav_intercept,
        "pav_r_squared": calibration.pav_r_squared,
        "regression_threshold_ldl_mg_dl": calibration.threshold_ldl_mg_dl,
        "noncalcified_decay_slope": calibration.nc_decay_slope,
        "noncalcified_decay_intercept": calibration.nc_decay_intercept,
        "calcification_fraction_assumed": DEFAULT_CALCIFICATION_FRACTION,
        "longest_calibration_followup_years": CALIBRATION_HORIZON_YEARS,
        "caveats": [
            "cross-trial fit: different populations, modalities and durations",
            "calcified volume never regresses in the model, matching the imaging trials",
            "the calcification transfer fraction is an assumption, not a measurement",
            "projections beyond 2 years are extrapolations",
        ],
    }


@dataclass(frozen=True)
class RegressionResult:
    """Two-compartment plaque trajectory under a sustained achieved LDL-C."""

    label: str
    achieved_ldl_mg_dl: float
    years: np.ndarray
    noncalcified_mm3: np.ndarray
    calcified_mm3: np.ndarray
    percent_atheroma_volume: np.ndarray
    notes: tuple[str, ...] = field(default=())

    @property
    def total_mm3(self) -> np.ndarray:
        return self.noncalcified_mm3 + self.calcified_mm3

    @property
    def noncalcified_depletion_fraction(self) -> float:
        start = float(self.noncalcified_mm3[0])
        return 0.0 if start <= 0 else 1.0 - float(self.noncalcified_mm3[-1]) / start

    @property
    def total_reduction_fraction(self) -> float:
        start = float(self.total_mm3[0])
        return 0.0 if start <= 0 else 1.0 - float(self.total_mm3[-1]) / start

    @property
    def calcified_increased(self) -> bool:
        return bool(self.calcified_mm3[-1] > self.calcified_mm3[0])

    @property
    def eradicated(self) -> bool:
        """Always False by construction: calcified burden is never removed."""
        return bool(np.any(self.total_mm3 <= 0.0))

    @property
    def extrapolated_beyond_evidence(self) -> bool:
        return bool(self.years[-1] > CALIBRATION_HORIZON_YEARS)


def simulate_regression(
    achieved_ldl_mg_dl: float,
    years: float = 2.0,
    *,
    label: str = "scenario",
    initial_noncalcified_mm3: float = 120.0,
    initial_calcified_mm3: float = 60.0,
    initial_pav: float = 38.0,
    calcification_fraction: float = DEFAULT_CALCIFICATION_FRACTION,
    step_years: float = 0.25,
    calibration: RegressionCalibration = DEFAULT_CALIBRATION,
) -> RegressionResult:
    """Project plaque compartments at a sustained achieved LDL-C.

    ``initial_*`` values are a representative obstructive-burden scenario, in the
    range reported by the IVUS/CT trials, not a measurement of any patient.
    """
    require_positive("achieved_ldl_mg_dl", achieved_ldl_mg_dl)
    require_positive("years", years)
    require_positive("step_years", step_years)
    require_positive("initial_noncalcified_mm3", initial_noncalcified_mm3)
    require_positive("initial_calcified_mm3", initial_calcified_mm3)
    require_positive("initial_pav", initial_pav)
    if not 0.0 <= calcification_fraction < 1.0:
        raise ValueError("calcification_fraction must be in [0, 1)")

    grid = np.arange(0.0, years + step_years * 0.5, step_years, dtype=float)
    if grid.size < 2:
        raise ValueError("horizon must contain at least two grid points")
    decay = calibration.noncalcified_decay_rate(achieved_ldl_mg_dl)
    noncalcified = initial_noncalcified_mm3 * np.exp(-decay * grid)
    lost = initial_noncalcified_mm3 - noncalcified
    calcified = initial_calcified_mm3 + calcification_fraction * lost
    pav = np.maximum(1.0, initial_pav + calibration.annual_pav_change(achieved_ldl_mg_dl) * grid)

    notes = [
        "two-compartment proxy: lipid-rich volume can regress, calcified volume does not",
        "cross-trial calibration; populations and imaging modalities differ",
    ]
    if years > CALIBRATION_HORIZON_YEARS:
        notes.append(
            f"projection beyond the {CALIBRATION_HORIZON_YEARS:g}-year calibration horizon "
            "is an extrapolation, not an observed result"
        )
    if achieved_ldl_mg_dl > calibration.threshold_ldl_mg_dl:
        notes.append("achieved LDL-C above the fitted regression threshold: plaque still grows")
    if calcification_fraction != DEFAULT_CALCIFICATION_FRACTION:
        notes.append(f"calcification_fraction set to {calcification_fraction:g} (assumption)")
    return RegressionResult(
        label=label,
        achieved_ldl_mg_dl=float(achieved_ldl_mg_dl),
        years=grid,
        noncalcified_mm3=noncalcified,
        calcified_mm3=calcified,
        percent_atheroma_volume=pav,
        notes=tuple(notes),
    )


def years_to_noncalcified_depletion(
    achieved_ldl_mg_dl: float,
    target_fraction: float = 0.95,
    calibration: RegressionCalibration = DEFAULT_CALIBRATION,
) -> float | None:
    """Extrapolated time to remove ``target_fraction`` of the lipid-rich pool.

    Returns ``None`` when the achieved LDL-C produces no net lipid-rich loss.
    This is an extrapolation of an exponential fit far beyond trial follow-up and
    must never be read as "time to cure".
    """
    if not 0.0 < target_fraction < 1.0:
        raise ValueError("target_fraction must be in (0, 1)")
    decay = calibration.noncalcified_decay_rate(achieved_ldl_mg_dl)
    if decay <= 0.0:
        return None
    return float(-math.log(1.0 - target_fraction) / decay)


def regimen_regression(
    regimen: Regimen,
    catalog: dict[str, Therapy],
    population: Population,
    years: float = 2.0,
    *,
    step_years: float = 0.25,
    calibration: RegressionCalibration = DEFAULT_CALIBRATION,
) -> RegressionResult:
    """Use a regimen's time-averaged achieved LDL-C to drive the plaque model."""
    validate_combination(catalog, regimen.therapy_ids)
    days = np.arange(0.0, years * DAYS_PER_YEAR + 1.0, 30.0)
    biomarkers = simulate_biomarkers(regimen, catalog, population, days)
    return simulate_regression(
        max(1.0, biomarkers.time_averaged_ldl_mg_dl),
        years=years,
        label=regimen.label,
        step_years=step_years,
        calibration=calibration,
    )


def imaging_anchor_frame(anchors: dict[str, TrialAnchor]) -> pd.DataFrame:
    """List the imaging evidence separately, so it is never mixed into MACE."""
    rows = [
        {
            "anchor_id": anchor.anchor_id,
            "therapy_id": anchor.therapy_id,
            "population_id": anchor.population_id,
            "endpoint": anchor.endpoint,
            "followup_years": anchor.followup_years,
            "n_participants": anchor.n_participants,
            "mace_credit": 0,
            "notes": anchor.notes,
            "source_url": anchor.source_url,
        }
        for anchor in anchors.values()
        if anchor.endpoint_class is EndpointClass.PLAQUE_IMAGING
    ]
    return pd.DataFrame(rows)


def regression_grid(
    catalog: dict[str, Therapy],
    population: Population,
    background: tuple[str, ...],
    therapy_ids: tuple[str, ...],
    years: float = 2.0,
    mode: CreditMode = CreditMode.STRICT,
) -> pd.DataFrame:
    """Compare single add-ons on structural regression, with no event claim."""
    rows = []
    candidates = [
        ("background_only", background),
        *((tid, (*background, tid)) for tid in therapy_ids),
    ]
    for label, ids in candidates:
        try:
            regimen = validate_combination(catalog, tuple(ids))
        except ValueError:
            continue
        result = regimen_regression(regimen, catalog, population, years=years)
        rows.append(
            {
                "add_on": label,
                "regimen": regimen.label,
                "population_id": population.population_id,
                "mode": mode.value,
                "achieved_ldl_mg_dl": result.achieved_ldl_mg_dl,
                "noncalcified_depletion_pct": 100.0 * result.noncalcified_depletion_fraction,
                "total_burden_reduction_pct": 100.0 * result.total_reduction_fraction,
                "calcified_volume_increased": result.calcified_increased,
                "pav_change_pp": float(
                    result.percent_atheroma_volume[-1] - result.percent_atheroma_volume[0]
                ),
                "years_to_95pct_lipid_depletion": years_to_noncalcified_depletion(
                    result.achieved_ldl_mg_dl
                ),
                "eradicated": result.eradicated,
                "extrapolated_beyond_evidence": result.extrapolated_beyond_evidence,
                "mace_credit_from_imaging": 0,
            }
        )
    frame = pd.DataFrame(rows)
    if not frame.empty:
        frame = frame.sort_values("noncalcified_depletion_pct", ascending=False).reset_index(
            drop=True
        )
    return frame
