"""Self-audit of ATLAS-1: does its own operating point survive delivery uncertainty?

Why this module exists
----------------------
The resolution step solved for a clearance operating point (1.284x potency, source
reduction 0.50, 24 months) that reaches 80.01% substrate removal. ATLAS-1 inherited
that number and scored burden reversal as 0.8001 -- a pass.

But that number is a **point estimate at central module strengths**. Earlier in
this study a remuscularisation "pass" was retracted for exactly this reason:
E[product of uncertain factors] is far below the product of their maxima. The
clearance operating point was never subjected to the same test. This module
subjects it to that test.

It also fixes a second, subtler error. ``atlas.bed_quality`` is derived from the
same nominal clearance, and ``atlas.translation_floor`` is derived from bed
quality. So the sequencing gain -- ATLAS-1's one new idea -- was scored as though
phase 1 reliably succeeds. The two arms are not independent: if the bed is not
actually cleared, there is no prepared bed to graft into. This module scores them
jointly.

Nothing here is evidence. Every module strength is an assumption, and the
distribution placed over it is an assumption too.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from ihd_solutions.atlas import (
    ALL_ARMS,
    ATLAS_CLEARANCE_MONTHS,
    ATLAS_INVENTED_TAG,
    ATLAS_POTENCY,
    ATLAS_SOURCE_REDUCTION,
    DISCLAIMER,
    translation_floor,
)
from ihd_solutions.cure import CRITERION_THRESHOLD
from ihd_solutions.gap_solver import (
    BASELINE_SCAR_FRACTION,
    REMUSCULARISATION_FACTORS,
    SCAR_RECOVERABLE_LVEF_PP,
)
from ihd_solutions.hypothesis import ClearanceParameters, simulate_clearance
from ihd_solutions.resolution import _ALL_STRENGTHS, CENTRAL
from ihd_solutions.units import require_fraction, require_positive

AUDIT_TAG = ATLAS_INVENTED_TAG

AUDIT_DISCLAIMER = (
    "THIS IS A SELF-AUDIT OF AN INVENTED CONSTRUCT, NOT EVIDENCE. It retracts a "
    "nominal pass by the same standard used to retract an earlier one. Module "
    "strengths and the distributions over them are assumptions."
)

#: Each module delivers somewhere between nothing and its assumed central
#: strength. This mirrors the Uniform(1, published) convention already used for
#: the graft factors, so the two arms are judged by the same rule.
DELIVERY_RULE = "uniform_zero_to_central"


def _central_vector() -> np.ndarray:
    return np.array([getattr(CENTRAL, name) for name in _ALL_STRENGTHS], dtype=float)


def clearance_draws(
    potency: float = ATLAS_POTENCY,
    source_reduction: float = ATLAS_SOURCE_REDUCTION,
    months: float = ATLAS_CLEARANCE_MONTHS,
    draws: int = 3000,
    seed: int = 20261107,
    source_uncertain: bool = True,
) -> np.ndarray:
    """Substrate removal under uncertain module delivery.

    Module strengths are drawn once and reused across potencies (common random
    numbers) so that a potency sweep is smooth rather than Monte Carlo noise.
    When ``source_uncertain`` the source reduction is drawn on [0, ceiling],
    because it is no better validated than the modules are.
    """
    require_positive("potency", potency)
    require_positive("months", months)
    require_positive("draws", draws)
    require_fraction("source_reduction", source_reduction, allow_one=True)
    count = int(draws)
    rng = np.random.default_rng(seed)
    unit = rng.random((count, len(_ALL_STRENGTHS)))
    delivered = np.minimum(1.0, unit * _central_vector() * potency)
    sources = (
        rng.uniform(0.0, source_reduction, count)
        if source_uncertain
        else np.full(count, source_reduction)
    )
    out = np.empty(count)
    for i in range(count):
        params = ClearanceParameters(**dict(zip(_ALL_STRENGTHS, delivered[i], strict=True)))
        out[i] = simulate_clearance(
            params,
            "audit",
            months=months,
            step_months=0.5,
            source_reduction=float(sources[i]),
        ).substrate_removal_fraction
    return out


@dataclass(frozen=True)
class ClearanceRobustness:
    """Distributional clearance, alongside the nominal number that misled us."""

    nominal: float
    median: float
    p05: float
    p95: float
    probability_of_target: float
    robust: bool

    @property
    def nominal_is_misleading(self) -> bool:
        return self.nominal >= CRITERION_THRESHOLD and not self.robust


def clearance_robustness(
    potency: float = ATLAS_POTENCY,
    source_reduction: float = ATLAS_SOURCE_REDUCTION,
    months: float = ATLAS_CLEARANCE_MONTHS,
    draws: int = 3000,
    seed: int = 20261107,
    source_uncertain: bool = True,
    target: float = CRITERION_THRESHOLD,
) -> ClearanceRobustness:
    from ihd_solutions.resolution import clearance_at

    require_fraction("target", target, allow_one=True)
    values = clearance_draws(potency, source_reduction, months, draws, seed, source_uncertain)
    probability = float(np.mean(values >= target))
    return ClearanceRobustness(
        nominal=float(clearance_at(potency, source_reduction, months)),
        median=float(np.median(values)),
        p05=float(np.quantile(values, 0.05)),
        p95=float(np.quantile(values, 0.95)),
        probability_of_target=probability,
        robust=probability >= 0.8,
    )


def potency_sweep(
    potencies: tuple[float, ...] = (1.284, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0),
    source_reduction: float = ATLAS_SOURCE_REDUCTION,
    months: float = ATLAS_CLEARANCE_MONTHS,
    draws: int = 2000,
    seed: int = 20261107,
    source_uncertain: bool = True,
) -> pd.DataFrame:
    rows = []
    for potency in potencies:
        values = clearance_draws(potency, source_reduction, months, draws, seed, source_uncertain)
        probability = float(np.mean(values >= CRITERION_THRESHOLD))
        rows.append(
            {
                "potency_multiplier": potency,
                "source_uncertain": source_uncertain,
                "horizon_months": months,
                "median_clearance_pct": 100.0 * float(np.median(values)),
                "p05_clearance_pct": 100.0 * float(np.quantile(values, 0.05)),
                "probability_of_target": probability,
                "robust_at_80pct": probability >= 0.8,
                "status": AUDIT_TAG,
            }
        )
    return pd.DataFrame(rows)


def minimum_robust_potency(
    source_reduction: float = ATLAS_SOURCE_REDUCTION,
    months: float = ATLAS_CLEARANCE_MONTHS,
    draws: int = 1200,
    seed: int = 20261107,
    source_uncertain: bool = True,
    required_probability: float = 0.8,
    upper: float = 256.0,
) -> float | None:
    """Smallest potency reaching the target probability, or None if unreachable.

    Returns None rather than a misleading number when no potency works -- the
    lesson learned when an empty "robust set" was once reported as "no factors
    needed".
    """
    require_fraction("required_probability", required_probability, allow_one=True)
    require_positive("upper", upper)

    def probability_at(potency: float) -> float:
        values = clearance_draws(potency, source_reduction, months, draws, seed, source_uncertain)
        return float(np.mean(values >= CRITERION_THRESHOLD))

    if probability_at(upper) < required_probability:
        return None
    lo, hi = 1.0, upper
    for _ in range(22):
        mid = 0.5 * (lo + hi)
        if probability_at(mid) >= required_probability:
            hi = mid
        else:
            lo = mid
    return float(hi)


def source_threshold_for_perfect_modules(
    months: float = ATLAS_CLEARANCE_MONTHS, target: float = CRITERION_THRESHOLD
) -> float | None:
    """Least source reduction that lets *perfect* modules hit the target.

    Above this source reduction the target is reachable in principle; below it,
    no amount of module potency helps. This sets the asymptotic ceiling on the
    success probability when source reduction is itself uncertain.
    """
    require_positive("months", months)
    require_fraction("target", target, allow_one=True)
    perfect = ClearanceParameters(**dict.fromkeys(_ALL_STRENGTHS, 1.0))

    def reached(source: float) -> bool:
        return (
            simulate_clearance(
                perfect, "audit", months=months, step_months=0.5, source_reduction=source
            ).substrate_removal_fraction
            >= target
        )

    if not reached(1.0):
        return None
    lo, hi = 0.0, 1.0
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if reached(mid):
            hi = mid
        else:
            lo = mid
    return float(hi)


def asymptotic_clearance_probability(
    source_reduction: float = ATLAS_SOURCE_REDUCTION,
    months: float = ATLAS_CLEARANCE_MONTHS,
) -> float:
    """Success probability as potency tends to infinity, source still uncertain.

    With unlimited potency every module saturates at 1.0, so the only remaining
    uncertainty is the source term. The ceiling is therefore the probability that
    the drawn source reduction exceeds the perfect-module threshold.
    """
    require_fraction("source_reduction", source_reduction, allow_one=True)
    threshold = source_threshold_for_perfect_modules(months)
    if threshold is None:
        return 0.0
    if threshold <= 0.0:
        return 1.0
    if source_reduction <= 0.0:
        # No source reduction at all, and perfect modules still need some.
        return 0.0
    return float(max(0.0, (source_reduction - threshold) / source_reduction))


@dataclass(frozen=True)
class CoupledGraft:
    """Graft outcome when the bed it lands in is itself uncertain."""

    probability_of_target: float
    median_gain_pp: float
    decoupled_probability: float
    median_bed_quality: float
    median_translation_floor: float

    @property
    def inflation_pp(self) -> float:
        """How much scoring the arms independently overstated the graft."""
        return 100.0 * (self.decoupled_probability - self.probability_of_target)


def coupled_graft_outcome(
    bed_coupling: float = 0.5,
    base_floor: float = 0.5,
    potency: float = ATLAS_POTENCY,
    source_reduction: float = ATLAS_SOURCE_REDUCTION,
    months: float = ATLAS_CLEARANCE_MONTHS,
    draws: int = 3000,
    seed: int = 20261107,
    required_pp: float = CRITERION_THRESHOLD * SCAR_RECOVERABLE_LVEF_PP,
) -> CoupledGraft:
    """Draw phase 1 and phase 2 together instead of assuming phase 1 worked.

    For each draw: uncertain modules give a clearance, that clearance gives a bed
    quality, the bed quality lifts the translation floor by ``bed_coupling``, and
    the graft factors are then drawn against that floor. At zero coupling this
    reduces exactly to the unsequenced graft.
    """
    require_fraction("bed_coupling", bed_coupling, allow_one=True)
    require_fraction("base_floor", base_floor, allow_one=True)
    require_positive("draws", draws)
    count = int(draws)
    removal = clearance_draws(potency, source_reduction, months, count, seed, True)
    quality = np.clip(removal / CRITERION_THRESHOLD, 0.0, 1.0)
    floors = np.clip(base_floor + bed_coupling * (1.0 - base_floor) * quality, base_floor, 1.0)
    rng = np.random.default_rng(seed + 1)
    delivered = np.ones(count)
    for factor in REMUSCULARISATION_FACTORS:
        delivered *= rng.uniform(1.0, factor.demonstrated_multiplier, size=count)
    # Paired comparison: the same factor draws and the same penalty quantiles are
    # pushed through the uncertain floor and through the nominal floor. The
    # difference is then purely the effect of bed uncertainty, with no Monte Carlo
    # noise in the contrast. At zero coupling the two are identical by construction.
    shared = rng.random(count)
    penalty = floors + (1.0 - floors) * shared
    gains = np.minimum(1.0, BASELINE_SCAR_FRACTION * delivered * penalty) * SCAR_RECOVERABLE_LVEF_PP
    nominal_floor = translation_floor(ALL_ARMS, bed_coupling, base_floor)
    penalty_nominal = nominal_floor + (1.0 - nominal_floor) * shared
    gains_nominal = (
        np.minimum(1.0, BASELINE_SCAR_FRACTION * delivered * penalty_nominal)
        * SCAR_RECOVERABLE_LVEF_PP
    )
    return CoupledGraft(
        probability_of_target=float(np.mean(gains >= required_pp)),
        median_gain_pp=float(np.median(gains)),
        decoupled_probability=float(np.mean(gains_nominal >= required_pp)),
        median_bed_quality=float(np.median(quality)),
        median_translation_floor=float(np.median(floors)),
    )


def coupled_sweep(
    couplings: tuple[float, ...] = (0.0, 0.25, 0.5, 0.75, 1.0),
    draws: int = 2500,
    seed: int = 20261107,
) -> pd.DataFrame:
    """The sequencing idea re-tested with an honestly uncertain bed."""
    rows = []
    for coupling in couplings:
        result = coupled_graft_outcome(coupling, draws=draws, seed=seed)
        rows.append(
            {
                "bed_coupling": coupling,
                "median_bed_quality": result.median_bed_quality,
                "median_translation_floor": result.median_translation_floor,
                "coupled_probability": result.probability_of_target,
                "decoupled_probability": result.decoupled_probability,
                "overstatement_pp": result.inflation_pp,
                "median_gain_pp": result.median_gain_pp,
                "robust_at_80pct": result.probability_of_target >= 0.8,
                "status": AUDIT_TAG,
            }
        )
    return pd.DataFrame(rows)


def audited_criteria(
    bed_coupling: float = 0.5, draws: int = 2500, seed: int = 20261107
) -> pd.DataFrame:
    """ATLAS-1's two contested criteria, scored nominally and distributionally."""
    robustness = clearance_robustness(draws=draws, seed=seed)
    coupled = coupled_graft_outcome(bed_coupling, draws=draws, seed=seed)
    myocardial_median = min(1.0, coupled.median_gain_pp / SCAR_RECOVERABLE_LVEF_PP)
    rows = [
        {
            "criterion": "burden_reversal",
            "nominal_score": robustness.nominal,
            "nominal_passes": robustness.nominal >= CRITERION_THRESHOLD,
            "probability_of_target": robustness.probability_of_target,
            "robust_passes": robustness.robust,
            "status": AUDIT_TAG,
        },
        {
            "criterion": "myocardial_restoration",
            "nominal_score": myocardial_median,
            "nominal_passes": myocardial_median >= CRITERION_THRESHOLD,
            "probability_of_target": coupled.probability_of_target,
            "robust_passes": coupled.probability_of_target >= 0.8,
            "status": AUDIT_TAG,
        },
    ]
    return pd.DataFrame(rows)


def audit_report(
    robustness: ClearanceRobustness,
    sweep: pd.DataFrame,
    coupled: pd.DataFrame,
    criteria: pd.DataFrame,
    robust_potency_fixed: float | None,
    ceiling: float,
    threshold: float | None,
) -> str:
    lines = [
        "ATLAS-1 SELF-AUDIT: does its own operating point survive uncertainty?",
        "===================================================================",
        AUDIT_DISCLAIMER,
        "",
        "THE ERROR BEING CORRECTED",
        "The resolution step solved for 1.284x potency at source reduction 0.50 and",
        "24 months, reaching 80.01% clearance. ATLAS-1 scored that as a burden-reversal",
        "pass. But 80.01% is a point estimate at central module strengths. An earlier",
        "remuscularisation pass was retracted for exactly this reason; clearance had",
        "never been given the same test. It has now.",
        "",
        "RESULT AT THE CLAIMED OPERATING POINT",
        f"  nominal clearance      {100 * robustness.nominal:6.2f}%",
        f"  median clearance       {100 * robustness.median:6.2f}%",
        f"  5th-95th percentile    {100 * robustness.p05:6.2f}% to {100 * robustness.p95:6.2f}%",
        f"  P(clearance >= 80%)    {100 * robustness.probability_of_target:6.2f}%",
        f"  robust                 {'yes' if robustness.robust else 'NO'}",
        "",
        "  The nominal pass does not survive. Burden reversal is retracted.",
        "",
        "WHAT POTENCY WOULD ACTUALLY BE NEEDED",
        "  potency   median    P(>=80%)   robust",
    ]
    for row in sweep.itertuples(index=False):
        lines.append(
            f"  {row.potency_multiplier:7.3f}  {row.median_clearance_pct:6.2f}%  "
            f"{100 * row.probability_of_target:7.2f}%   {'yes' if row.robust_at_80pct else 'no'}"
        )
    lines += [
        "",
        f"  With source reduction uncertain on [0, {ATLAS_SOURCE_REDUCTION:.2f}], perfect modules",
        f"  still need source >= {threshold:.4f} to reach 80% at 24 months."
        if threshold is not None
        else "  Perfect modules cannot reach the target at any source reduction.",
        f"  So as potency tends to infinity the probability ceiling is only {100 * ceiling:.1f}%.",
        "  No potency makes clearance robust while the source term is uncertain.",
        "",
        (
            f"  If source reduction were *guaranteed* at {ATLAS_SOURCE_REDUCTION:.2f}, robustness would need"
            f"  {robust_potency_fixed:.1f}x potency -- not 1.284x."
            if robust_potency_fixed is not None
            else "  Even with source reduction guaranteed, no potency tested is robust."
        ),
        "",
        "THE SECOND ERROR: THE TWO ARMS ARE NOT INDEPENDENT",
        "Bed quality is derived from clearance, and the graft's translation floor is",
        "derived from bed quality. Scoring the graft against a nominal 80% cleared bed",
        "assumed phase 1 worked. Scored jointly:",
        "",
        "  coupling   median bed   coupled P   decoupled P   overstated by",
    ]
    for row in coupled.itertuples(index=False):
        lines.append(
            f"  {row.bed_coupling:8.2f}   {row.median_bed_quality:10.3f}   "
            f"{100 * row.coupled_probability:9.2f}%   {100 * row.decoupled_probability:11.2f}%   "
            f"{row.overstatement_pp:+6.2f} pp"
        )
    lines += [
        "",
        "  Sequencing cannot prepare a bed that phase 1 fails to clear. The idea is",
        "  sound in structure and worthless in practice until clearance itself is robust.",
        "",
        "CORRECTED SCORECARD",
    ]
    for row in criteria.itertuples(index=False):
        lines.append(
            f"  {row.criterion:24s} nominal {row.nominal_score:.3f} "
            f"({'pass' if row.nominal_passes else 'fail'})   "
            f"P(target) {100 * row.probability_of_target:5.2f}%   "
            f"robust {'yes' if row.robust_passes else 'NO'}"
        )
    lines += [
        "",
        "VERDICT",
        "ATLAS-1 does not reach 3 of 5 criteria. Under the same distributional standard",
        "applied to every other construct it reaches 2 of 5, the same as its",
        "predecessors. It remains the best candidate on construction quality -- every",
        "arm is load-bearing and its one new idea is falsifiable -- but the claim that",
        "it advances the criterion count is withdrawn.",
        "",
        "The corrected ranking of open problems is: clearance under delivery",
        "uncertainty first, remuscularisation second, randomised proof third.",
        "The highest-leverage real-world action is unchanged and needs no new biology:",
        "the PCSK9-editing outcomes trial.",
    ]
    return "\n".join(lines) + "\n"


def run_audit_study(
    output_dir: Path | str = Path("results/audit"),
    draws: int = 2500,
    seed: int = 20261107,
) -> dict[str, object]:
    """Write the self-audit. Quarantined outside the evidence tree, like ATLAS-1."""
    root = Path(output_dir)
    tables = root / "tables"
    reports = root / "reports"
    tables.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)

    robustness = clearance_robustness(draws=draws, seed=seed)
    sweep_uncertain = potency_sweep(draws=draws, seed=seed, source_uncertain=True)
    sweep_fixed = potency_sweep(draws=draws, seed=seed, source_uncertain=False)
    sweep = pd.concat([sweep_uncertain, sweep_fixed], ignore_index=True)
    coupled = coupled_sweep(draws=draws, seed=seed)
    criteria = audited_criteria(draws=draws, seed=seed)
    threshold = source_threshold_for_perfect_modules()
    ceiling = asymptotic_clearance_probability()
    robust_fixed = minimum_robust_potency(
        draws=max(600, draws // 3), seed=seed, source_uncertain=False
    )
    robust_uncertain = minimum_robust_potency(
        draws=max(600, draws // 3), seed=seed, source_uncertain=True
    )

    pd.DataFrame(
        [
            {
                "nominal_clearance_pct": 100 * robustness.nominal,
                "median_clearance_pct": 100 * robustness.median,
                "p05_clearance_pct": 100 * robustness.p05,
                "p95_clearance_pct": 100 * robustness.p95,
                "probability_of_target": robustness.probability_of_target,
                "robust": robustness.robust,
                "nominal_is_misleading": robustness.nominal_is_misleading,
                "status": AUDIT_TAG,
            }
        ]
    ).to_csv(tables / "audit_operating_point.csv", index=False)
    sweep.to_csv(tables / "audit_potency_sweep.csv", index=False)
    coupled.to_csv(tables / "audit_coupled_sequencing.csv", index=False)
    criteria.to_csv(tables / "audit_corrected_criteria.csv", index=False)
    pd.DataFrame(
        [
            {
                "source_uncertain": True,
                "perfect_module_source_threshold": threshold,
                "asymptotic_probability_ceiling": ceiling,
                "minimum_robust_potency": robust_uncertain,
                "robust_potency_exists": robust_uncertain is not None,
                "status": AUDIT_TAG,
            },
            {
                "source_uncertain": False,
                "perfect_module_source_threshold": threshold,
                "asymptotic_probability_ceiling": 1.0,
                "minimum_robust_potency": robust_fixed,
                "robust_potency_exists": robust_fixed is not None,
                "status": AUDIT_TAG,
            },
        ]
    ).to_csv(tables / "audit_reachability.csv", index=False)

    text = audit_report(
        robustness, sweep_uncertain, coupled, criteria, robust_fixed, ceiling, threshold
    )
    (reports / "audit_report.txt").write_text(text)
    (root / "DISCLAIMER.txt").write_text(DISCLAIMER + "\n\n" + AUDIT_DISCLAIMER + "\n")

    return {
        "output_dir": str(root),
        "outputs": 7,
        "status": "self_audit_retracts_nominal_pass",
        "burden_reversal_retracted": bool(robustness.nominal_is_misleading),
        "cures_declared": 0,
        "disclaimer": AUDIT_DISCLAIMER,
    }
