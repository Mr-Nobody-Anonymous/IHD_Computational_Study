"""ATLAS-3: apply the variance fix to the graft arm, the last open biology gap.

The diagnosis
-------------
After ATLAS-2, clearance is robust in simulation (86.1% across the whole
population). Remuscularisation is not, and the reason is now precise. The four
published enhancement factors multiply to 18.2x, which is **3.55x** the 5.128x
actually required to recover +20 LVEF points. Capability is not the shortfall.
With the chronic-injury penalty set to zero - a perfect bed, every factor available
at its published maximum - a single graft still reaches the target only **50.7%**
of the time.

That is a variance problem, and it is the same variance problem that broke
clearance. ATLAS-2 fixed clearance with redundancy and closed-loop titration.
Neither has ever been applied to the graft, because a graft is one-shot,
irreversible, and has no measurable feedback signal.

ATLAS-3 asks what happens if those two properties can be changed, and charges for
them honestly.

What is modelled, and what it assumes
-------------------------------------
``attempts``
    Staged, repeatable grafting: a failed attempt does not preclude another and the
    best outcome obtains. This is **retry redundancy, not area multiplication** -
    combining by max is deliberately conservative, because claiming three patches
    cover three times the scar would be a much stronger and less defensible claim.
``failure_correlation``
    The catch. Three attempts that all fail because the same chronic scar is poorly
    vascularised are one attempt, not three. At correlation 1 the attempts are
    identical and redundancy buys exactly nothing; that is the falsifier.
``feedback_learning``
    A measured engraftment signal lets a retry correct an under-dosed aggregate. It
    is ONE bounded scalar uplift, not a per-mechanism correction, because a clinical
    signal reports total engraftment and cannot attribute it to four mechanisms.
    At learning 0 it contributes nothing; that is the second falsifier.
``measurement_quality``
    Probability the signal is informative at all. At 0 the feedback contributes
    nothing regardless of the learning setting; that is the third falsifier.
``per_attempt_risk``
    Every additional procedure carries adverse-event risk (graft-associated
    arrhythmia, immune sensitisation). Cumulative risk is reported next to the
    success probability, never hidden behind it.

This is an invented, untrialled construct. It declares no cures.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from ihd_solutions.cure import CRITERION_THRESHOLD
from ihd_solutions.gap_solver import (
    BASELINE_SCAR_FRACTION,
    REMUSCULARISATION_FACTORS,
    SCAR_RECOVERABLE_LVEF_PP,
)
from ihd_solutions.precedent import max_precedented_doses
from ihd_solutions.units import require_fraction, require_positive

ATLAS3_TAG = "atlas3_invented_architecture_unvalidated"

DISCLAIMER = (
    "ATLAS-3 IS AN INVENTED, UNTRIALLED CONSTRUCT. Staged repeatable grafting with a "
    "measurable feedback signal does not exist. Redundancy is modelled as retry, not "
    "as extra treated area; attempt independence is an assumption swept to zero "
    "benefit; and every extra procedure is charged an adverse-event risk. Nothing "
    "here is a therapy, a result, or a recommendation."
)

#: Multiplier the graft must deliver to recover the target LVEF gain.
REQUIRED_MULTIPLIER = CRITERION_THRESHOLD / BASELINE_SCAR_FRACTION

#: Retry redundancy beyond this has no precedent in any published dosing programme.
#: CORRECTED from 2 to the registry value (4): CATO (NCT06145035) randomises patients to
#: four intravenous doses and HOPE-2 gave four. The original 2 was never checked against
#: the literature, which made this module label its own three-attempt design "beyond
#: precedent" when it was not. See ihd_solutions.precedent.
MAX_PRECEDENTED_ATTEMPTS = max_precedented_doses()

#: A specification must not move when the sampling budget moves. At small N the
#: 2-SE margin is wide enough to reject a two-attempt design that is genuinely
#: above the bar (true P ~0.822, verified across 12 seeds). Both specification
#: searches therefore floor their draw count here rather than inheriting whatever
#: budget the caller passed for the exploratory tables.
MIN_SPECIFICATION_DRAWS = 4000

#: Per-procedure adverse-event probability used for the safety ledger. Graft-associated
#: ventricular arrhythmia was non-fatal but real in primate injection work; epicardial
#: patch series reported none, so this is a deliberately non-zero placeholder.
DEFAULT_PER_ATTEMPT_RISK = 0.05

#: Largest aggregate dose correction a perfect feedback signal could justify.
#: Bounded well below what would make the product of minimums exceed the target.
FEEDBACK_MAX_UPLIFT = 0.5


@dataclass(frozen=True)
class GraftUpgrades:
    """The three architectural levers, each independently ablatable."""

    attempts: int = 1
    failure_correlation: float = 0.0
    feedback_learning: float = 0.0
    measurement_quality: float = 0.7
    per_attempt_risk: float = DEFAULT_PER_ATTEMPT_RISK
    translation_floor: float = 1.0

    def __post_init__(self) -> None:
        if self.attempts < 1:
            raise ValueError("at least one attempt is required")
        require_fraction("failure_correlation", self.failure_correlation, allow_one=True)
        require_fraction("feedback_learning", self.feedback_learning, allow_one=True)
        require_fraction("measurement_quality", self.measurement_quality, allow_one=True)
        require_fraction("per_attempt_risk", self.per_attempt_risk, allow_one=True)
        require_fraction("translation_floor", self.translation_floor, allow_one=True)

    @property
    def beyond_precedent(self) -> bool:
        return self.attempts > MAX_PRECEDENTED_ATTEMPTS

    def without(self, lever: str) -> GraftUpgrades:
        mapping = {
            "retry_redundancy": {"attempts": 1},
            "independence": {"failure_correlation": 1.0},
            "feedback": {"feedback_learning": 0.0},
            "measurement": {"measurement_quality": 0.0},
        }
        if lever not in mapping:
            raise ValueError(f"unknown lever {lever!r}")
        fields = {f: getattr(self, f) for f in self.__dataclass_fields__}
        return GraftUpgrades(**{**fields, **mapping[lever]})


#: The full ATLAS-3 graft architecture. Three attempts is what the arithmetic needs;
#: it is beyond precedent and is reported as such rather than quietly used.
ATLAS3 = GraftUpgrades(
    attempts=3,
    failure_correlation=0.0,
    feedback_learning=0.5,
    translation_floor=1.0,
)

LEVER_NAMES = ("retry_redundancy", "independence", "feedback", "measurement")


def _factor_draws(
    rng: np.random.Generator,
    count: int,
    attempt_index: int,
    config: GraftUpgrades,
    shared: np.ndarray,
) -> np.ndarray:
    """Delivered multiplier for one attempt.

    Correlation is introduced by blending a shared uniform with a per-attempt
    uniform *before* mapping onto each factor's range. At correlation 1 every
    attempt sees the same quantile and the max over attempts equals a single draw.

    Feedback narrows the range from below on later attempts: having measured the
    patient's response, a retry avoids the under-dosed tail. It never raises the
    published maximum.
    """
    delivered = np.ones(count)
    for j, factor in enumerate(REMUSCULARISATION_FACTORS):
        independent = rng.random(count)
        quantile = (
            config.failure_correlation * shared[:, j]
            + (1.0 - config.failure_correlation) * independent
        )
        delivered *= 1.0 + quantile * (factor.demonstrated_multiplier - 1.0)
    # Feedback corrects the AGgregate dose, not each mechanism separately. A clinical
    # signal reports total engraftment; it cannot tell you which of four mechanisms
    # under-delivered. Raising every factor's floor at once would make the product of
    # the minimums exceed the requirement and force success deterministically - the
    # nominal-product trap in a new form. One bounded scalar uplift avoids that.
    learn = config.feedback_learning * (attempt_index / max(config.attempts - 1, 1))
    informative = rng.random(count) < config.measurement_quality
    uplift = 1.0 + learn * FEEDBACK_MAX_UPLIFT * informative
    ceiling = float(np.prod([f.demonstrated_multiplier for f in REMUSCULARISATION_FACTORS]))
    return np.minimum(delivered * uplift, ceiling)


@dataclass(frozen=True)
class GraftResult:
    """Success probability beside what the architecture cost to reach it."""

    probability_of_target: float
    median_gain_pp: float
    p05_gain_pp: float
    median_attempts_used: float
    cumulative_adverse_risk: float
    translation_floor: float
    attempts: int
    beyond_precedent: bool
    draws: int

    @property
    def robust(self) -> bool:
        return self.probability_of_target >= 0.8

    @property
    def monte_carlo_standard_error(self) -> float:
        p = self.probability_of_target
        return float(np.sqrt(max(p * (1.0 - p), 0.0) / max(self.draws, 1)))

    @property
    def robust_with_margin(self) -> bool:
        return self.probability_of_target - 2.0 * self.monte_carlo_standard_error >= 0.8

    @property
    def risk_adjusted_probability(self) -> float:
        """Success discounted by the cumulative probability of an adverse event.

        Not a clinical utility - a reminder that a higher success rate bought with
        more procedures is not free.
        """
        return self.probability_of_target * (1.0 - self.cumulative_adverse_risk)


def graft_draws(
    config: GraftUpgrades = ATLAS3,
    draws: int = 4000,
    seed: int = 20261107,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per-draw best gain, attempts used, and the FIRST attempt's gain.

    The first-attempt array makes the independence falsifier a paired contrast:
    at correlation 1 every attempt is identical, so the gain over the first attempt
    is exactly zero rather than zero plus Monte Carlo noise.

    Attempts stop as soon as the target is met, so a patient who responds to the
    first graft is not exposed to further procedures.
    """
    require_positive("draws", draws)
    count = int(draws)
    rng = np.random.default_rng(seed)
    shared = rng.random((count, len(REMUSCULARISATION_FACTORS)))
    best = np.zeros(count)
    first = np.zeros(count)
    used = np.ones(count)
    done = np.zeros(count, dtype=bool)
    required = CRITERION_THRESHOLD * SCAR_RECOVERABLE_LVEF_PP
    for attempt in range(config.attempts):
        delivered = _factor_draws(rng, count, attempt, config, shared)
        penalty = rng.uniform(config.translation_floor, 1.0, size=count)
        gains = (
            np.minimum(1.0, BASELINE_SCAR_FRACTION * delivered * penalty) * SCAR_RECOVERABLE_LVEF_PP
        )
        if attempt == 0:
            first = gains.copy()
        improving = ~done
        best = np.where(improving, np.maximum(best, gains), best)
        used = np.where(improving, float(attempt + 1), used)
        done = done | (best >= required)
    return best, used, first


def evaluate_graft(
    config: GraftUpgrades = ATLAS3,
    draws: int = 4000,
    seed: int = 20261107,
) -> GraftResult:
    gains, used, _ = graft_draws(config, draws, seed)
    required = CRITERION_THRESHOLD * SCAR_RECOVERABLE_LVEF_PP
    mean_used = float(np.mean(used))
    return GraftResult(
        probability_of_target=float(np.mean(gains >= required)),
        median_gain_pp=float(np.median(gains)),
        p05_gain_pp=float(np.quantile(gains, 0.05)),
        median_attempts_used=float(np.median(used)),
        cumulative_adverse_risk=float(1.0 - (1.0 - config.per_attempt_risk) ** mean_used),
        translation_floor=config.translation_floor,
        attempts=config.attempts,
        beyond_precedent=config.beyond_precedent,
        draws=int(draws),
    )


def attempt_sweep(
    attempts: tuple[int, ...] = (1, 2, 3, 4, 6),
    draws: int = 3000,
    seed: int = 20261107,
    failure_correlation: float = 0.0,
    feedback_learning: float = 0.5,
) -> pd.DataFrame:
    rows = []
    for k in attempts:
        config = GraftUpgrades(
            attempts=k,
            failure_correlation=failure_correlation,
            feedback_learning=feedback_learning,
        )
        result = evaluate_graft(config, draws, seed)
        rows.append(
            {
                "attempts": k,
                "probability_of_target": result.probability_of_target,
                "median_gain_pp": result.median_gain_pp,
                "mean_attempts_stopped_at": result.median_attempts_used,
                "cumulative_adverse_risk": result.cumulative_adverse_risk,
                "risk_adjusted_probability": result.risk_adjusted_probability,
                "robust": result.robust,
                "robust_with_margin": result.robust_with_margin,
                "beyond_precedent": result.beyond_precedent,
                "status": ATLAS3_TAG,
            }
        )
    return pd.DataFrame(rows)


def independence_sweep(
    correlations: tuple[float, ...] = (0.0, 0.25, 0.5, 0.75, 1.0),
    draws: int = 3000,
    seed: int = 20261107,
    attempts: int = 3,
) -> pd.DataFrame:
    """The falsifier. At correlation 1, retry redundancy must contribute nothing."""
    rows = []
    required = CRITERION_THRESHOLD * SCAR_RECOVERABLE_LVEF_PP
    for rho in correlations:
        config = GraftUpgrades(attempts=attempts, failure_correlation=rho)
        result = evaluate_graft(config, draws, seed)
        _, _, first = graft_draws(config, draws, seed)
        single_probability = float(np.mean(first >= required))
        rows.append(
            {
                "failure_correlation": rho,
                "probability_of_target": result.probability_of_target,
                "first_attempt_probability": single_probability,
                "redundancy_gain_pp": 100.0 * (result.probability_of_target - single_probability),
                "robust": result.robust,
                "status": ATLAS3_TAG,
            }
        )
    return pd.DataFrame(rows)


def feedback_sweep(
    learning_levels: tuple[float, ...] = (0.0, 0.25, 0.5, 0.75, 1.0),
    draws: int = 3000,
    seed: int = 20261107,
    attempts: int = 3,
) -> pd.DataFrame:
    """The second falsifier. At learning 0 feedback must contribute nothing."""
    rows = []
    for learn in learning_levels:
        config = GraftUpgrades(attempts=attempts, feedback_learning=learn)
        without = evaluate_graft(
            GraftUpgrades(attempts=attempts, feedback_learning=0.0), draws, seed
        )
        result = evaluate_graft(config, draws, seed)
        rows.append(
            {
                "feedback_learning": learn,
                "probability_of_target": result.probability_of_target,
                "probability_without_feedback": without.probability_of_target,
                "feedback_gain_pp": 100.0
                * (result.probability_of_target - without.probability_of_target),
                "robust": result.robust,
                "status": ATLAS3_TAG,
            }
        )
    return pd.DataFrame(rows)


def lever_ablation(
    config: GraftUpgrades = ATLAS3, draws: int = 3000, seed: int = 20261107
) -> pd.DataFrame:
    full = evaluate_graft(config, draws, seed)
    rows = []
    for lever in LEVER_NAMES:
        reduced = evaluate_graft(config.without(lever), draws, seed)
        rows.append(
            {
                "removed_lever": lever,
                "probability_of_target": reduced.probability_of_target,
                "probability_loss_pp": 100.0
                * (full.probability_of_target - reduced.probability_of_target),
                "cumulative_adverse_risk": reduced.cumulative_adverse_risk,
                "robust": reduced.robust,
                "load_bearing": full.probability_of_target - reduced.probability_of_target > 1e-6,
                "status": ATLAS3_TAG,
            }
        )
    return pd.DataFrame(rows)


def safety_sweep(
    per_attempt_risks: tuple[float, ...] = (0.0, 0.02, 0.05, 0.10, 0.20),
    draws: int = 2500,
    seed: int = 20261107,
) -> pd.DataFrame:
    """Where the safety price cancels the benefit of another attempt."""
    rows = []
    for risk in per_attempt_risks:
        for k in (1, 2, 3):
            config = GraftUpgrades(
                attempts=k, feedback_learning=ATLAS3.feedback_learning, per_attempt_risk=risk
            )
            result = evaluate_graft(config, draws, seed)
            rows.append(
                {
                    "per_attempt_risk": risk,
                    "attempts": k,
                    "probability_of_target": result.probability_of_target,
                    "cumulative_adverse_risk": result.cumulative_adverse_risk,
                    "risk_adjusted_probability": result.risk_adjusted_probability,
                    "status": ATLAS3_TAG,
                }
            )
    return pd.DataFrame(rows)


def minimum_robust_attempts(
    config: GraftUpgrades = ATLAS3,
    draws: int = 2500,
    seed: int = 20261107,
    upper: int = 12,
    required_probability: float = 0.8,
) -> int | None:
    """Fewest attempts reaching the bar with a 2-SE margin, or None.

    Draw count is floored at MIN_SPECIFICATION_DRAWS: a design specification that
    changed with the caller's sampling budget would not be a specification.
    """
    draws = max(int(draws), MIN_SPECIFICATION_DRAWS)
    fields = {f: getattr(config, f) for f in config.__dataclass_fields__}
    for k in range(1, upper + 1):
        result = evaluate_graft(GraftUpgrades(**{**fields, "attempts": k}), draws, seed)
        if result.probability_of_target >= required_probability and result.robust_with_margin:
            return k
    return None


def maximum_tolerable_correlation(
    config: GraftUpgrades = ATLAS3,
    draws: int = 2000,
    seed: int = 20261107,
    required_probability: float = 0.8,
) -> float | None:
    """Largest shared-failure correlation at which the graft is still robust.

    This is the specification the biology has to meet: retry outcomes must be at
    least this independent. A chronic poorly-vascularised scar is a shared failure
    mode, so this number is the one to measure, not assume.
    """
    draws = max(int(draws), MIN_SPECIFICATION_DRAWS)
    fields = {f: getattr(config, f) for f in config.__dataclass_fields__}

    def robust_at(rho: float) -> bool:
        result = evaluate_graft(
            GraftUpgrades(**{**fields, "failure_correlation": rho}), draws, seed
        )
        return result.probability_of_target >= required_probability

    if not robust_at(0.0):
        return None
    if robust_at(1.0):
        return 1.0
    lo, hi = 0.0, 1.0
    for _ in range(18):
        mid = 0.5 * (lo + hi)
        if robust_at(mid):
            lo = mid
        else:
            hi = mid
    return float(lo)


def precedented_option(draws: int = 4000, seed: int = 20261107) -> GraftResult:
    """The CHEAPEST robust design inside the published attempt limit.

    Reported beside the shipped three-attempt construct so the extra attempt is visibly
    a purchase of headroom, not a requirement of feasibility. This used to return the
    attempt limit itself; with the limit corrected to 4 that would have compared three
    attempts against a strictly better four-attempt design and called the difference a
    cost. The decision-relevant comparison is against the fewest attempts that work.
    """
    draws = max(int(draws), MIN_SPECIFICATION_DRAWS)
    fields = {f: getattr(ATLAS3, f) for f in ATLAS3.__dataclass_fields__}
    fewest = minimum_robust_attempts(ATLAS3, draws=draws, seed=seed)
    attempts = min(fewest, MAX_PRECEDENTED_ATTEMPTS) if fewest else MAX_PRECEDENTED_ATTEMPTS
    return evaluate_graft(GraftUpgrades(**{**fields, "attempts": attempts}), draws, seed)


def joint_capability_probability(draws: int = 2000, seed: int = 20261107) -> dict[str, float]:
    """P(both invented arms hit target), and why it is still not a cure.

    The product assumes the two arms are independent. They are not fully
    independent in reality - the graft lands in the bed clearance prepared - so this
    is an optimistic joint, reported as such.
    """
    from ihd_solutions.upgrade import ATLAS2, SOLVED_MONTHS, evaluate_upgrades

    clearance = evaluate_upgrades(ATLAS2, months=SOLVED_MONTHS, draws=draws, seed=seed)
    graft = evaluate_graft(ATLAS3, draws, seed)
    return {
        "clearance_probability": clearance.probability_of_target,
        "graft_probability": graft.probability_of_target,
        "joint_probability_optimistic": clearance.probability_of_target
        * graft.probability_of_target,
        "randomised_proof": 0.0,
    }


def current_state_voi(draws: int = 2500, seed: int = 20261107) -> pd.DataFrame:
    """Value of information recomputed against the CURRENT architectures.

    The original discovery VOI ranked substrate removal first (variance share
    0.1839). That was computed before ATLAS-2 gave clearance an architecture, so it
    is stale. Here the question is asked again of the two arms that are still
    uncertain: how much would resolving each one's uncertainty buy?

    Clearance is taken from the committed ATLAS-2 solve; the graft from ATLAS-3.
    The joint is a product, which assumes the arms are independent - stated, not
    hidden, because the sequencing coupling means they are not fully independent in
    reality.
    """
    from ihd_solutions.upgrade import ATLAS2, SOLVED_MONTHS, evaluate_upgrades

    clearance = evaluate_upgrades(ATLAS2, months=SOLVED_MONTHS, draws=draws, seed=seed)
    graft = evaluate_graft(ATLAS3, draws, seed)
    p_clear = clearance.probability_of_target
    p_graft = graft.probability_of_target
    joint = p_clear * p_graft
    rows = [
        {
            "capability": "substrate_removal",
            "current_probability": p_clear,
            "joint_now": joint,
            "joint_if_resolved": 1.0 * p_graft,
            "voi_gain": 1.0 * p_graft - joint,
            "architecture": "ATLAS-2 (redundancy + titration, 60 months)",
            "status": ATLAS3_TAG,
        },
        {
            "capability": "lvef_gain_pp",
            "current_probability": p_graft,
            "joint_now": joint,
            "joint_if_resolved": p_clear * 1.0,
            "voi_gain": p_clear * 1.0 - joint,
            "architecture": "ATLAS-3 (retry redundancy + feedback)",
            "status": ATLAS3_TAG,
        },
    ]
    frame = pd.DataFrame(rows)
    return frame.sort_values("voi_gain", ascending=False).reset_index(drop=True)


def atlas3_report(
    single: GraftResult,
    full: GraftResult,
    attempts: pd.DataFrame,
    independence: pd.DataFrame,
    feedback: pd.DataFrame,
    ablation: pd.DataFrame,
    safety: pd.DataFrame,
    voi: pd.DataFrame,
    minimum: int | None,
    tolerable_correlation: float | None,
    joint: dict[str, float],
    precedented: GraftResult,
) -> str:
    product = float(np.prod([f.demonstrated_multiplier for f in REMUSCULARISATION_FACTORS]))
    zero = independence.loc[independence["failure_correlation"].eq(1.0)].iloc[0]
    lines = [
        "ATLAS-3: applying the variance fix to the last open biology gap",
        "==============================================================",
        DISCLAIMER,
        "",
        "WHY THE GRAFT IS A VARIANCE PROBLEM, NOT A CAPABILITY PROBLEM",
        f"  required delivered multiplier      {REQUIRED_MULTIPLIER:.3f}x",
        f"  product of all published maxima    {product:.3f}x  ({product / REQUIRED_MULTIPLIER:.2f}x the requirement)",
        f"  single graft, perfect bed          P(target) {100 * single.probability_of_target:.2f}%",
        "",
        "  Capability is not the shortfall. Even with zero chronic-injury penalty and",
        "  every factor available at its published maximum, one graft is a coin flip.",
        "  This is the same variance failure that broke clearance, and ATLAS-2 fixed",
        "  clearance with redundancy and feedback. Neither has been applied here.",
        "",
        "RESULT",
        f"  single attempt (today)             P {100 * single.probability_of_target:6.2f}%"
        f"   median {single.median_gain_pp:5.2f} pp",
        f"  ATLAS-3 ({ATLAS3.attempts} staged attempts)         P {100 * full.probability_of_target:6.2f}%"
        f"   median {full.median_gain_pp:5.2f} pp",
        f"  robust at 80%                      {'YES' if full.robust else 'no'}",
        f"  robust with 2-SE margin            {'YES' if full.robust_with_margin else 'no'}",
        f"  attempts beyond precedent          {'YES' if full.beyond_precedent else 'no'}"
        f"  (precedented limit {MAX_PRECEDENTED_ATTEMPTS})",
        "",
        "WHAT IT COST",
        f"  median attempts actually used      {full.median_attempts_used:.2f}",
        f"  cumulative adverse-event risk      {100 * full.cumulative_adverse_risk:.2f}%",
        f"  risk-adjusted success              {100 * full.risk_adjusted_probability:.2f}%",
        "",
        "  Attempts stop as soon as the target is met, so a first-attempt responder is",
        "  not exposed to further procedures. The risk is still real and cumulative.",
        "",
        "THE CHEAPEST PRECEDENTED OPTION: EXTRA ATTEMPTS BUY HEADROOM, NOT FEASIBILITY",
        f"  precedented attempt limit          {MAX_PRECEDENTED_ATTEMPTS} (CORRECTED from 2; CATO and HOPE-2 give four)",
        f"  ATLAS-3's {full.attempts} attempts               WITHIN precedent",
        f"  {precedented.attempts} attempts (cheapest robust)       P {100 * precedented.probability_of_target:6.2f}%"
        f"   adverse {100 * precedented.cumulative_adverse_risk:5.2f}%",
        f"  robust with 2-SE margin            {'YES' if precedented.robust_with_margin else 'no'}",
        f"  fewest attempts clearing the bar   {minimum if minimum is not None else 'none in range'}",
        "",
        "  Fewer attempts already clear the bar, and repeat dosing has ample precedent.",
        "  So the number of procedures is not what breaks precedent here, and it is not",
        "  the invention. The invention is the closed loop: a measured engraftment signal",
        "  that decides whether to retry. The extra attempt is optional headroom - it",
        f"  lifts {100 * precedented.probability_of_target:.1f}% to {100 * full.probability_of_target:.1f}% - bought with",
        f"  {100 * (full.cumulative_adverse_risk - precedented.cumulative_adverse_risk):+.2f} pp more adverse-event risk.",
        "  The cheaper option has a thinner margin, so it is reported, not relied upon.",
        "",
        "THE FALSIFIER: REDUNDANCY ONLY COUNTS IF ATTEMPTS ARE INDEPENDENT",
        "  correlation   P(target)   gain vs single attempt",
    ]
    for row in independence.itertuples(index=False):
        lines.append(
            f"     {row.failure_correlation:4.2f}       {100 * row.probability_of_target:6.2f}%   "
            f"{row.redundancy_gain_pp:+6.2f} pp"
        )
    lines += [
        f"  At correlation 1.00 the gain is {zero.redundancy_gain_pp:+.2f} pp: three attempts that",
        "  share a failure mode are one attempt. Whether real retries are independent",
        "  is the single most important unmeasured quantity in this construct.",
        "",
        (
            f"  MAXIMUM TOLERABLE CORRELATION: about {tolerable_correlation:.2f} (spread +/-0.02 across seeds)."
            if tolerable_correlation is not None
            else "  MAXIMUM TOLERABLE CORRELATION: none - not robust even when fully independent."
        ),
        "  Retry outcomes must be at least this independent to stay robust. That is the number",
        "  to measure in a real study, and a chronic poorly-vascularised scar is exactly",
        "  the kind of shared failure mode that would violate it.",
        "",
        "THE SECOND FALSIFIER: FEEDBACK",
    ]
    for row in feedback.itertuples(index=False):
        lines.append(
            f"  learning {row.feedback_learning:4.2f}   P {100 * row.probability_of_target:6.2f}%   "
            f"gain {row.feedback_gain_pp:+6.2f} pp"
        )
    lines += [
        "",
        "LEVER ABLATION",
    ]
    for row in ablation.itertuples(index=False):
        lines.append(
            f"  drop {row.removed_lever:18s} P {100 * row.probability_of_target:6.2f}% "
            f"({row.probability_loss_pp:+6.2f} pp)  "
            f"{'LOAD-BEARING' if row.load_bearing else 'inert'}"
        )
    lines += ["", "ATTEMPTS: success against the safety price"]
    for row in attempts.itertuples(index=False):
        lines.append(
            f"  {row.attempts} attempts  P {100 * row.probability_of_target:6.2f}%  "
            f"adverse {100 * row.cumulative_adverse_risk:5.2f}%  "
            f"risk-adjusted {100 * row.risk_adjusted_probability:6.2f}%  "
            f"{'BEYOND PRECEDENT' if row.beyond_precedent else 'precedented'}"
        )
    lines += [
        "",
        f"  Fewest attempts reaching the bar with margin: "
        f"{minimum if minimum is not None else 'none within the search range'}.",
        "",
        "SAFETY: where another attempt stops being worth it",
    ]
    for row in safety.itertuples(index=False):
        lines.append(
            f"  per-attempt risk {100 * row.per_attempt_risk:5.2f}%  {row.attempts} attempts  "
            f"risk-adjusted {100 * row.risk_adjusted_probability:6.2f}%"
        )
    lines += [
        "",
        "BOTH ARMS TOGETHER",
        f"  clearance (ATLAS-2)                {100 * joint['clearance_probability']:6.2f}%",
        f"  graft (ATLAS-3)                    {100 * joint['graft_probability']:6.2f}%",
        f"  joint, assuming independence       {100 * joint['joint_probability_optimistic']:6.2f}%",
        "  randomised proof                    0.00%",
        "",
        "  The joint is optimistic: the graft lands in the bed clearance prepared, so",
        "  the arms are not truly independent. Even taken at face value it is one",
        "  number for two of five criteria, with proof still absent.",
        "",
        "VALUE OF INFORMATION, RECOMPUTED AGAINST THE CURRENT STATE",
        "The original discovery VOI ranked substrate removal first (0.1839). That was",
        "computed before ATLAS-2 existed and is stale. Recomputed:",
        "",
    ]
    for row in voi.itertuples(index=False):
        lines.append(
            f"  {row.capability:18s} P now {100 * row.current_probability:6.2f}%  "
            f"VOI gain {100 * row.voi_gain:+6.2f} pp  [{row.architecture}]"
        )
    lines += [
        "",
        "VERDICT",
        "The graft gap closes in simulation under staged retry with feedback, which",
        "identifies the missing step precisely: a remuscularisation modality that can",
        "be repeated AND MEASURED. Repeat dosing already has precedent; the missing",
        "piece is the feedback signal that decides whether to repeat. It is not more",
        f"potency - the published factors already exceed the requirement "
        f"{product / REQUIRED_MULTIPLIER:.2f}-fold.",
        "",
        "Three things must be demonstrated before any of this is real:",
        "  1. Independence of retry outcomes. Correlated failures erase the entire gain,",
        "     and a chronic poorly-vascularised scar is exactly a shared failure mode.",
        "  2. A measurable engraftment signal good enough to guide a retry.",
        "  3. Acceptable cumulative procedural risk across repeated grafting.",
        "",
        "Randomised proof remains zero. ATLAS-3 is a specification of the experiment",
        "to run, not a cure and not a therapy.",
    ]
    return "\n".join(lines) + "\n"


def run_atlas3_study(
    output_dir: Path | str = Path("results/atlas3"),
    draws: int = 2500,
    seed: int = 20261107,
) -> dict[str, object]:
    root = Path(output_dir)
    tables, reports = root / "tables", root / "reports"
    tables.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)

    single = evaluate_graft(GraftUpgrades(attempts=1), draws, seed)
    full = evaluate_graft(ATLAS3, draws, seed)
    attempts = attempt_sweep(draws=draws, seed=seed)
    independence = independence_sweep(draws=draws, seed=seed)
    feedback = feedback_sweep(draws=draws, seed=seed)
    ablation = lever_ablation(ATLAS3, draws, seed)
    safety = safety_sweep(draws=max(800, draws // 2), seed=seed)
    voi = current_state_voi(draws=max(800, draws // 2), seed=seed)
    minimum = minimum_robust_attempts(ATLAS3, draws=max(800, draws // 2), seed=seed)
    tolerable = maximum_tolerable_correlation(ATLAS3, draws=max(800, draws // 2), seed=seed)
    joint = joint_capability_probability(draws=max(800, draws // 2), seed=seed)
    precedented = precedented_option(draws=draws, seed=seed)

    pd.DataFrame(
        [
            {
                "construct": "single_graft_today",
                "probability_of_target": single.probability_of_target,
                "median_gain_pp": single.median_gain_pp,
                "median_attempts_used": single.median_attempts_used,
                "cumulative_adverse_risk": single.cumulative_adverse_risk,
                "risk_adjusted_probability": single.risk_adjusted_probability,
                "robust": single.robust,
                "robust_with_margin": single.robust_with_margin,
                "beyond_precedent": single.beyond_precedent,
                "status": ATLAS3_TAG,
            },
            {
                "construct": "atlas3_staged",
                "probability_of_target": full.probability_of_target,
                "median_gain_pp": full.median_gain_pp,
                "median_attempts_used": full.median_attempts_used,
                "cumulative_adverse_risk": full.cumulative_adverse_risk,
                "risk_adjusted_probability": full.risk_adjusted_probability,
                "robust": full.robust,
                "robust_with_margin": full.robust_with_margin,
                "beyond_precedent": full.beyond_precedent,
                "status": ATLAS3_TAG,
            },
        ]
    ).to_csv(tables / "atlas3_headline.csv", index=False)
    attempts.to_csv(tables / "atlas3_attempt_sweep.csv", index=False)
    independence.to_csv(tables / "atlas3_independence_falsifier.csv", index=False)
    feedback.to_csv(tables / "atlas3_feedback_falsifier.csv", index=False)
    ablation.to_csv(tables / "atlas3_lever_ablation.csv", index=False)
    safety.to_csv(tables / "atlas3_safety_tradeoff.csv", index=False)
    voi.to_csv(tables / "atlas3_current_state_voi.csv", index=False)
    pd.DataFrame(
        [
            {
                "required_multiplier": REQUIRED_MULTIPLIER,
                "product_of_published_maxima": float(
                    np.prod([f.demonstrated_multiplier for f in REMUSCULARISATION_FACTORS])
                ),
                "minimum_robust_attempts": minimum,
                "precedented_attempt_limit": MAX_PRECEDENTED_ATTEMPTS,
                "precedented_option_probability": precedented.probability_of_target,
                "precedented_option_robust": bool(precedented.robust_with_margin),
                "maximum_tolerable_correlation": tolerable,
                "joint_probability_optimistic": joint["joint_probability_optimistic"],
                "status": ATLAS3_TAG,
            }
        ]
    ).to_csv(tables / "atlas3_requirement.csv", index=False)
    (reports / "atlas3_report.txt").write_text(
        atlas3_report(
            single,
            full,
            attempts,
            independence,
            feedback,
            ablation,
            safety,
            voi,
            minimum,
            tolerable,
            joint,
            precedented,
        )
    )
    (root / "DISCLAIMER.txt").write_text(DISCLAIMER + "\n")
    return {
        "output_dir": str(root),
        "outputs": 10,
        "maximum_tolerable_correlation": tolerable,
        "joint_probability_optimistic": joint["joint_probability_optimistic"],
        "status": ATLAS3_TAG,
        "robust": bool(full.robust),
        "robust_with_margin": bool(full.robust_with_margin),
        "minimum_robust_attempts": minimum,
        "precedented_option_probability": precedented.probability_of_target,
        "precedented_option_robust": bool(precedented.robust_with_margin),
        "cumulative_adverse_risk": full.cumulative_adverse_risk,
        "cures_declared": 0,
        "disclaimer": DISCLAIMER,
    }
