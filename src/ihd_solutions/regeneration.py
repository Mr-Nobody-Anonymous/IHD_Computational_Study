"""Regenerative cell-sheet scenarios with honest small-sample uncertainty."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats

from ihd_solutions.evidence import CreditMode
from ihd_solutions.units import require_fraction, require_positive


@dataclass(frozen=True)
class RegenerationAnchor:
    """Responder counts from the eight-patient uncontrolled pivotal experience."""

    n_treated: int = 8
    lvef_evaluable: int = 8
    lvef_responders: int = 2
    vo2_evaluable: int = 7
    vo2_responders: int = 4
    nyha_evaluable: int = 8
    nyha_responders: int = 8
    lvef_threshold_pct: float = 5.0
    vo2_threshold_pct: float = 10.0

    def __post_init__(self) -> None:
        require_positive("n_treated", self.n_treated)
        for evaluable, responders in (
            (self.lvef_evaluable, self.lvef_responders),
            (self.vo2_evaluable, self.vo2_responders),
            (self.nyha_evaluable, self.nyha_responders),
        ):
            if evaluable <= 0 or evaluable > self.n_treated:
                raise ValueError("evaluable counts must lie in [1, n_treated]")
            if not 0 <= responders <= evaluable:
                raise ValueError("responders must lie in [0, evaluable]")
        require_positive("lvef_threshold_pct", self.lvef_threshold_pct)
        require_positive("vo2_threshold_pct", self.vo2_threshold_pct)


DEFAULT_REGENERATION_ANCHOR = RegenerationAnchor()


@dataclass(frozen=True)
class ResponderPosterior:
    endpoint: str
    responders: int
    evaluable: int
    alpha: float
    beta: float

    @property
    def mean(self) -> float:
        return float(self.alpha / (self.alpha + self.beta))

    @property
    def observed_fraction(self) -> float:
        return float(self.responders / self.evaluable)

    def interval(self, mass: float = 0.90) -> tuple[float, float]:
        require_fraction("mass", mass)
        if mass <= 0:
            raise ValueError("mass must be positive")
        tail = (1.0 - mass) / 2.0
        low, high = stats.beta.ppf([tail, 1.0 - tail], self.alpha, self.beta)
        return float(low), float(high)

    def predictive_responders(self, cohort_size: int, draws: int, seed: int) -> np.ndarray:
        require_positive("cohort_size", cohort_size)
        require_positive("draws", draws)
        rng = np.random.default_rng(seed)
        probabilities = rng.beta(self.alpha, self.beta, size=int(draws))
        return rng.binomial(int(cohort_size), probabilities)


def beta_binomial_posterior(
    endpoint: str, responders: int, evaluable: int, prior: tuple[float, float] = (1.0, 1.0)
) -> ResponderPosterior:
    """Uniform-prior Beta posterior; deliberately wide because n is tiny."""
    if evaluable <= 0:
        raise ValueError("evaluable must be positive")
    if not 0 <= responders <= evaluable:
        raise ValueError("responders must lie in [0, evaluable]")
    a_prior, b_prior = prior
    require_positive("prior alpha", a_prior)
    require_positive("prior beta", b_prior)
    return ResponderPosterior(
        endpoint=endpoint,
        responders=int(responders),
        evaluable=int(evaluable),
        alpha=float(a_prior + responders),
        beta=float(b_prior + evaluable - responders),
    )


def anchor_posteriors(
    anchor: RegenerationAnchor = DEFAULT_REGENERATION_ANCHOR,
) -> dict[str, ResponderPosterior]:
    return {
        "lvef": beta_binomial_posterior("lvef", anchor.lvef_responders, anchor.lvef_evaluable),
        "peak_vo2": beta_binomial_posterior(
            "peak_vo2", anchor.vo2_responders, anchor.vo2_evaluable
        ),
        "nyha": beta_binomial_posterior("nyha", anchor.nyha_responders, anchor.nyha_evaluable),
    }


def event_credit(mode: CreditMode) -> tuple[float, str]:
    """Regenerative sheets never receive MACE or mortality credit in this model."""
    return 0.0, (
        f"{mode.value} mode: conditional approval rests on eight uncontrolled patients with "
        "functional endpoints; no event or mortality benefit is credited"
    )


def surveillance_requirements() -> dict[str, object]:
    return {
        "approval_type": "conditional_time_limited",
        "confirmation_target_cases": 75,
        "confirmation_window_years": 7,
        "remuscularization_proven": False,
        "mechanism_emphasised": "paracrine",
        "required_monitoring": (
            "ventricular arrhythmia surveillance",
            "tumorigenicity surveillance",
            "immunological monitoring and graft persistence",
            "all-case post-marketing registry with functional endpoints",
        ),
    }


def simulate_regeneration(
    anchor: RegenerationAnchor = DEFAULT_REGENERATION_ANCHOR,
    cohort_size: int = 75,
    draws: int = 2000,
    seed: int = 20261106,
) -> pd.DataFrame:
    """Posterior-predictive responder counts for the confirmation cohort."""
    require_positive("cohort_size", cohort_size)
    require_positive("draws", draws)
    rows = []
    credit, reason = event_credit(CreditMode.STRICT)
    for index, (endpoint, posterior) in enumerate(anchor_posteriors(anchor).items()):
        counts = posterior.predictive_responders(cohort_size, draws, seed + index)
        low, high = posterior.interval()
        rows.append(
            {
                "endpoint": endpoint,
                "observed_responders": posterior.responders,
                "observed_evaluable": posterior.evaluable,
                "observed_fraction": posterior.observed_fraction,
                "posterior_mean_probability": posterior.mean,
                "posterior_p5": low,
                "posterior_p95": high,
                "cohort_size": cohort_size,
                "predicted_responders_median": float(np.median(counts)),
                "predicted_responders_p5": float(np.quantile(counts, 0.05)),
                "predicted_responders_p95": float(np.quantile(counts, 0.95)),
                "mace_credit": credit,
                "mortality_credit": 0.0,
                "credit_reason": reason,
                "remuscularization_proven": False,
            }
        )
    return pd.DataFrame(rows)
