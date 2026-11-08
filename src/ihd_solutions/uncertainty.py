"""Uncertainty propagation and sensitivity analysis for the scenario models."""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, replace

import numpy as np
import pandas as pd
from scipy import stats

from ihd_solutions.catalog import Regimen, Therapy
from ihd_solutions.evidence import CreditMode, TrialAnchor
from ihd_solutions.pharmacology import Adherence
from ihd_solutions.population import Population
from ihd_solutions.risk import DEFAULT_RISK_PARAMETERS, RiskParameters, simulate_events
from ihd_solutions.units import DAYS_PER_YEAR, Z_95, require_positive


@dataclass(frozen=True)
class AdherencePrior:
    """Beta prior for the fraction of prescribed therapy actually taken."""

    mean: float = 0.85
    concentration: float = 20.0

    def __post_init__(self) -> None:
        if not 0 < self.mean < 1:
            raise ValueError("adherence mean must lie strictly between 0 and 1")
        require_positive("concentration", self.concentration)

    @property
    def alpha(self) -> float:
        return self.mean * self.concentration

    @property
    def beta(self) -> float:
        return (1.0 - self.mean) * self.concentration


DEFAULT_ADHERENCE_PRIOR = AdherencePrior()


def truncated_normal(
    rng: np.random.Generator, mean: float, low: float, high: float, size: int
) -> np.ndarray:
    """Normal draws clipped to a reported interval, with sigma from that interval."""
    require_positive("size", size)
    if not low <= mean <= high:
        raise ValueError("interval must bracket the mean")
    sigma = (high - low) / (2.0 * Z_95)
    if sigma <= 0:
        return np.full(int(size), float(mean))
    draws = rng.normal(mean, sigma, size=int(size))
    return np.clip(draws, low, high)


def lognormal_from_ci(
    rng: np.random.Generator, point: float, low: float, high: float, size: int
) -> np.ndarray:
    """Log-scale draws reproducing a reported 95% confidence interval."""
    require_positive("point", point)
    require_positive("size", size)
    if not 0 < low <= point <= high:
        raise ValueError("confidence interval must bracket the point estimate")
    sigma = (math.log(high) - math.log(low)) / (2.0 * Z_95)
    if sigma <= 0:
        return np.full(int(size), float(point))
    return rng.lognormal(math.log(point), sigma, size=int(size))


@dataclass(frozen=True)
class MonteCarloResult:
    regimen_label: str
    population_id: str
    mode: CreditMode
    events_avoided_per_1000: np.ndarray
    time_averaged_ldl_mg_dl: np.ndarray
    #: Per-draw untreated and treated cumulative incidence. Optional so existing
    #: callers that only need events avoided keep working. These are what the
    #: risk-normalisation criterion is actually computed from, so scoring that
    #: criterion under uncertainty requires them rather than a point estimate.
    baseline_risk: np.ndarray | None = None
    treated_risk: np.ndarray | None = None

    def risk_normalisation_draws(self, reference_risk: float = 0.0198) -> np.ndarray | None:
        """Per-draw risk-normalisation score, or None if risks were not recorded.

        Mirrors the point-estimate formula in ``composite.score_criteria`` exactly,
        so the nominal score and this distribution are the same quantity.
        """
        if self.baseline_risk is None or self.treated_risk is None:
            return None
        baseline = np.asarray(self.baseline_risk, dtype=float)
        treated = np.asarray(self.treated_risk, dtype=float)
        gap = baseline - reference_risk
        score = np.zeros_like(baseline)
        usable = gap > 0.0
        score[usable] = np.clip(1.0 - (treated[usable] - reference_risk) / gap[usable], 0.0, 1.0)
        return score

    def summary(self) -> dict[str, object]:
        values = self.events_avoided_per_1000
        return {
            "regimen": self.regimen_label,
            "population_id": self.population_id,
            "mode": self.mode.value,
            "draws": int(values.size),
            "events_avoided_median": float(np.median(values)),
            "events_avoided_p5": float(np.quantile(values, 0.05)),
            "events_avoided_p95": float(np.quantile(values, 0.95)),
            "probability_of_benefit": float(np.mean(values > 0)),
            "probability_of_harm": float(np.mean(values < 0)),
            "ldl_median_mg_dl": float(np.median(self.time_averaged_ldl_mg_dl)),
        }


def monte_carlo_regimen(
    regimen: Regimen,
    catalog: dict[str, Therapy],
    population: Population,
    anchors: dict[str, TrialAnchor],
    days: np.ndarray,
    mode: CreditMode = CreditMode.STRICT,
    draws: int = 200,
    seed: int = 20261106,
    params: RiskParameters = DEFAULT_RISK_PARAMETERS,
    adherence_prior: AdherencePrior = DEFAULT_ADHERENCE_PRIOR,
) -> MonteCarloResult:
    """Propagate efficacy, hazard-ratio, adherence and baseline-rate uncertainty."""
    require_positive("draws", draws)
    draws = int(draws)
    rng = np.random.default_rng(seed)
    ldl_draws = {
        therapy_id: truncated_normal(
            rng,
            catalog[therapy_id].ldl_pct,
            catalog[therapy_id].ldl_pct_ci_low,
            catalog[therapy_id].ldl_pct_ci_high,
            draws,
        )
        for therapy_id in regimen.therapy_ids
    }
    adherence_draws = {
        therapy_id: rng.beta(adherence_prior.alpha, adherence_prior.beta, size=draws)
        for therapy_id in regimen.therapy_ids
    }
    rate_draws = lognormal_from_ci(
        rng,
        population.baseline_annual_mace_rate,
        population.baseline_rate_ci_low,
        population.baseline_rate_ci_high,
        draws,
    )
    hr_draws: dict[str, np.ndarray] = {}
    for anchor in anchors.values():
        if anchor.therapy_id in regimen.therapy_ids and anchor.has_effect_estimate:
            # A neutral trial is sampled around the null, not around its noisy
            # point estimate: uncertainty is retained, direction is not invented.
            centre = (
                1.0
                if params.neutral_ci_policy == "null_effect" and anchor.neutral_at_report
                else float(anchor.hr)
            )
            hr_draws[anchor.anchor_id] = lognormal_from_ci(
                rng, centre, float(anchor.ci_low), float(anchor.ci_high), draws
            )

    # The null policy has already been applied to the sampling centre above, so
    # each draw is used verbatim instead of being collapsed back to HR = 1.
    draw_params = replace(params, neutral_ci_policy="point_estimate")

    avoided = np.empty(draws)
    ldl_means = np.empty(draws)
    treated_risk = np.empty(draws)
    years = float(days[-1] - days[0]) / DAYS_PER_YEAR
    for index in range(draws):
        sampled_catalog = dict(catalog)
        for therapy_id, values in ldl_draws.items():
            therapy = catalog[therapy_id]
            ldl_value = float(values[index])
            sampled_catalog[therapy_id] = replace(
                therapy,
                ldl_pct=ldl_value,
                ldl_pct_ci_low=min(therapy.ldl_pct_ci_low, ldl_value),
                ldl_pct_ci_high=max(therapy.ldl_pct_ci_high, ldl_value),
            )
        sampled_anchors = dict(anchors)
        for anchor_id, values in hr_draws.items():
            anchor = anchors[anchor_id]
            hr_value = float(values[index])
            sampled_anchors[anchor_id] = replace(
                anchor,
                hr=hr_value,
                ci_low=min(float(anchor.ci_low), hr_value),
                ci_high=max(float(anchor.ci_high), hr_value),
            )
        sampled_population = replace(
            population,
            baseline_annual_mace_rate=float(rate_draws[index]),
            baseline_rate_ci_low=min(population.baseline_rate_ci_low, float(rate_draws[index])),
            baseline_rate_ci_high=max(population.baseline_rate_ci_high, float(rate_draws[index])),
        )
        sampled_adherence = {
            therapy_id: Adherence(taking_fraction=float(values[index]))
            for therapy_id, values in adherence_draws.items()
        }
        simulation = simulate_events(
            regimen,
            sampled_catalog,
            sampled_population,
            sampled_anchors,
            days,
            mode=mode,
            params=draw_params,
            adherence=sampled_adherence,
        )
        avoided[index] = simulation.events_avoided_per_1000
        ldl_means[index] = simulation.biomarkers.time_averaged_ldl_mg_dl
        treated_risk[index] = float(simulation.cumulative_incidence[-1])
    # Matches the empty-regimen analytic baseline used by composite._incidence,
    # evaluated at the same sampled baseline rate as the treated arm.
    baseline_risk = -np.expm1(-rate_draws * years)
    return MonteCarloResult(
        regimen_label=regimen.label,
        population_id=population.population_id,
        mode=mode,
        events_avoided_per_1000=avoided,
        time_averaged_ldl_mg_dl=ldl_means,
        baseline_risk=baseline_risk,
        treated_risk=treated_risk,
    )


def one_at_a_time(
    model: Callable[[float], float],
    base: float,
    low: float,
    high: float,
    name: str = "parameter",
) -> dict[str, object]:
    """Tornado row for a scalar parameter sweep around a base value."""
    if not low <= base <= high:
        raise ValueError("sweep interval must bracket the base value")
    base_value = float(model(base))
    low_value = float(model(low))
    high_value = float(model(high))
    return {
        "parameter": name,
        "base_input": base,
        "low_input": low,
        "high_input": high,
        "base_output": base_value,
        "low_output": low_value,
        "high_output": high_value,
        "span": abs(high_value - low_value),
        "direction": "increasing" if high_value >= low_value else "decreasing",
    }


def tornado(
    model: Callable[[float], float], sweeps: Mapping[str, tuple[float, float, float]]
) -> pd.DataFrame:
    rows = [
        one_at_a_time(model, base, low, high, name) for name, (base, low, high) in sweeps.items()
    ]
    return (
        pd.DataFrame(rows)
        .sort_values("span", ascending=False, kind="stable")
        .reset_index(drop=True)
    )


def spearman_sensitivity(samples: Mapping[str, np.ndarray], output: np.ndarray) -> pd.DataFrame:
    """Rank-correlation sensitivity of an output to sampled inputs."""
    output = np.asarray(output, dtype=float)
    rows = []
    for name, values in samples.items():
        values = np.asarray(values, dtype=float)
        if values.shape != output.shape:
            raise ValueError(f"{name}: sample shape {values.shape} != output {output.shape}")
        if np.ptp(values) == 0 or np.ptp(output) == 0:
            rho, pvalue = 0.0, 1.0
        else:
            result = stats.spearmanr(values, output)
            rho, pvalue = float(result.statistic), float(result.pvalue)
        rows.append({"parameter": name, "spearman_rho": rho, "p_value": pvalue})
    frame = pd.DataFrame(rows)
    frame["abs_rho"] = frame["spearman_rho"].abs()
    return frame.sort_values("abs_rho", ascending=False, kind="stable").reset_index(drop=True)


def rank_stability(results: Sequence[MonteCarloResult], top_n: int = 3) -> pd.DataFrame:
    """Probability that each regimen lands in the top ``top_n`` across draws."""
    if not results:
        raise ValueError("no Monte Carlo results supplied")
    require_positive("top_n", top_n)
    n_draws = results[0].events_avoided_per_1000.size
    if any(result.events_avoided_per_1000.size != n_draws for result in results):
        raise ValueError("all results must use the same number of draws")
    matrix = np.vstack([result.events_avoided_per_1000 for result in results])
    k = min(int(top_n), matrix.shape[0])
    order = np.argsort(-matrix, axis=0, kind="stable")[:k]
    counts = np.zeros(matrix.shape[0], dtype=int)
    for row in range(k):
        np.add.at(counts, order[row], 1)
    return (
        pd.DataFrame(
            {
                "regimen": [result.regimen_label for result in results],
                "probability_top_rank": counts / n_draws,
                "median_events_avoided": np.median(matrix, axis=1),
            }
        )
        .sort_values("probability_top_rank", ascending=False, kind="stable")
        .reset_index(drop=True)
    )
