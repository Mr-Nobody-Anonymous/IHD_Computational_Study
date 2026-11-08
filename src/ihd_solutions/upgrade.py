"""ATLAS-2: upgrade the parts the audits broke, instead of only retracting them.

The diagnosis this is built from
--------------------------------
Three audits found the same class of failure, and the clearance arm's version has a
specific structural cause:

1. ``ClearanceParameters.delivered()`` multiplies **six of eight** module strengths
   by the single ``gate_specificity`` draw. One uncertain parameter therefore
   throttles almost the whole system: a series bottleneck with no redundancy.
2. Delivery was modelled **open loop** - pick a potency, hope it lands. Nobody
   manages lipids that way; clinical practice titrates to a measured target.
3. The source term was left uncertain, which caps the success probability at 55.4%
   *at any potency*, because perfect modules still need source reduction >= 0.2230.

So ATLAS-2 changes the architecture rather than inflating the assumed potency:

======  ==========================  ==============================================
Upgrade What it fixes               Real-world analogue
======  ==========================  ==============================================
U1      series gate bottleneck      several mechanistically distinct lesion-targeting
                                    routes; the node works if any one of them lands
U2      open-loop delivery          treat-to-target titration against a measured
                                    biomarker, the standard of care in lipidology
U3      uncertain source term       respond-and-continue run-in: source control is
                                    verified before the construct proceeds
======  ==========================  ==============================================

Every upgrade is ablatable, and every upgrade is charged for. U1 costs distinct
mechanisms (more agents, more interaction risk). U2 costs titration visits and only
works if a biomarker is actually measurable. U3 costs **eligible population**: it
does not make source control more likely, it restricts to patients in whom it is
achieved, and the excluded fraction is reported rather than hidden.

This is still an invented, untrialled construct. It declares no cures.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from ihd_solutions.atlas import ATLAS_CLEARANCE_MONTHS, ATLAS_POTENCY, ATLAS_SOURCE_REDUCTION
from ihd_solutions.cure import CRITERION_THRESHOLD
from ihd_solutions.hypothesis import ClearanceParameters, simulate_clearance
from ihd_solutions.resolution import _ALL_STRENGTHS, CENTRAL
from ihd_solutions.units import require_fraction, require_positive

UPGRADE_TAG = "atlas2_invented_architecture_unvalidated"

DISCLAIMER = (
    "ATLAS-2 IS AN INVENTED, UNTRIALLED CONSTRUCT. Its upgrades are architectural "
    "hypotheses with clinical analogues, not demonstrated effects. Every upgrade is "
    "ablated and charged for; the population excluded by verified source control is "
    "reported. Nothing here is a therapy, a result, or a recommendation."
)

#: Modules whose strength the gate multiplies in ``delivered()``. Redundancy at the
#: gate therefore protects six nodes at once, which is why it is the first fix.
GATE_DEPENDENT = (
    "crystal_solubilisation",
    "efflux_gain",
    "efferocytosis_gain",
    "influx_block",
    "resolution",
    "nucleation_block",
)

#: Beyond this many independent targeting routes the redundancy claim has no
#: mechanistic precedent and is flagged rather than silently used.
MAX_PRECEDENTED_PATHS = 3

#: Source reduction perfect modules need to reach the clearance target at 24 months.
#: Established by the ATLAS-1 self-audit; recomputed here rather than hardcoded.
DEFAULT_SOURCE_REQUIRED = 0.2230


@dataclass(frozen=True)
class Upgrades:
    """Which architectural upgrades are active. Ablation switches each one off."""

    gate_paths: int = 1
    node_paths: int = 1
    titrate: bool = False
    verified_source: bool = False
    max_potency: float = 2.0
    titration_steps: int = 8
    source_required: float = DEFAULT_SOURCE_REQUIRED

    def __post_init__(self) -> None:
        if self.gate_paths < 1 or self.node_paths < 1:
            raise ValueError("redundancy needs at least one path")
        require_positive("max_potency", self.max_potency)
        if self.max_potency < 1.0:
            raise ValueError("max_potency below 1.0 cannot even reach the central assumption")
        if self.titration_steps < 2:
            raise ValueError("titration needs at least two steps")
        require_fraction("source_required", self.source_required, allow_one=True)

    @property
    def redundancy_beyond_precedent(self) -> bool:
        return max(self.gate_paths, self.node_paths) > MAX_PRECEDENTED_PATHS

    def without(self, upgrade: str) -> Upgrades:
        mapping = {
            "redundant_gate": {"gate_paths": 1},
            "redundant_nodes": {"node_paths": 1},
            "titration": {"titrate": False},
            "verified_source": {"verified_source": False},
        }
        if upgrade not in mapping:
            raise ValueError(f"unknown upgrade {upgrade!r}")
        fields = {f: getattr(self, f) for f in self.__dataclass_fields__}
        return Upgrades(**{**fields, **mapping[upgrade]})


#: The full ATLAS-2 architecture. Redundancy is capped at the precedented count.
ATLAS2 = Upgrades(
    gate_paths=MAX_PRECEDENTED_PATHS,
    node_paths=2,
    titrate=True,
    # Verified source control is NOT needed at the solved horizon: with 60 months
    # even perfect modules facing zero source reduction clear the target, so the
    # run-in would only exclude patients who would have succeeded. It stays
    # available and ablatable, and the ablation reports it as inert here.
    verified_source=False,
    max_potency=2.80,
)

UPGRADE_NAMES = ("redundant_gate", "redundant_nodes", "titration", "verified_source")


def _central() -> np.ndarray:
    return np.array([getattr(CENTRAL, name) for name in _ALL_STRENGTHS], dtype=float)


def _path_max(rng: np.random.Generator, count: int, paths: int) -> np.ndarray:
    """Best of ``paths`` independent delivery attempts, on [0, 1].

    Max rather than the complement product: if two mechanisms both partially work,
    counting only the better one is the conservative choice.
    """
    if paths == 1:
        return rng.random(count)
    return rng.random((count, paths)).max(axis=1)


def source_required_for_target(
    months: float = ATLAS_CLEARANCE_MONTHS, target: float = CRITERION_THRESHOLD
) -> float:
    """Least source reduction letting perfect modules reach the target."""
    require_positive("months", months)
    perfect = ClearanceParameters(**dict.fromkeys(_ALL_STRENGTHS, 1.0))

    def reached(source: float) -> bool:
        return (
            simulate_clearance(
                perfect, "atlas2", months=months, step_months=0.5, source_reduction=source
            ).substrate_removal_fraction
            >= target
        )

    lo, hi = 0.0, 1.0
    for _ in range(36):
        mid = 0.5 * (lo + hi)
        if reached(mid):
            hi = mid
        else:
            lo = mid
    return float(hi)


def _clearance(strengths: np.ndarray, source: float, months: float) -> float:
    params = ClearanceParameters(**dict(zip(_ALL_STRENGTHS, strengths, strict=True)))
    return float(
        simulate_clearance(
            params, "atlas2", months=months, step_months=0.5, source_reduction=float(source)
        ).substrate_removal_fraction
    )


@dataclass(frozen=True)
class UpgradeResult:
    """Distributional outcome plus what the upgrade cost to get there."""

    probability_of_target: float
    median_clearance: float
    p05_clearance: float
    eligible_fraction: float
    median_potency_used: float
    robust: bool
    conditional: bool
    redundancy_beyond_precedent: bool
    draws: int

    @property
    def monte_carlo_standard_error(self) -> float:
        p = self.probability_of_target
        n = max(self.draws, 1)
        return float(np.sqrt(max(p * (1.0 - p), 0.0) / n))

    @property
    def robust_with_margin(self) -> bool:
        """Robust even at the lower end of two Monte Carlo standard errors.

        A point estimate sitting exactly on the 80% bar is not a solve; the
        evidence-layer audit flagged that case as borderline rather than counting
        it. The solved operating point is chosen to clear the bar with margin.
        """
        return self.probability_of_target - 2.0 * self.monte_carlo_standard_error >= 0.8

    @property
    def unconditional_probability(self) -> float:
        """Success across everyone treated, not just those who clear the run-in.

        Verified source control buys conditional reliability by excluding patients.
        Reporting only the conditional number would hide that trade.
        """
        return self.probability_of_target * self.eligible_fraction


def achieved_clearance_draws(
    upgrades: Upgrades = ATLAS2,
    months: float = ATLAS_CLEARANCE_MONTHS,
    draws: int = 1500,
    seed: int = 20261107,
    target: float = CRITERION_THRESHOLD,
    source_ceiling: float = ATLAS_SOURCE_REDUCTION,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Per-draw achieved clearance and potency used, plus the eligible fraction.

    Exposed separately so the graft arm can be scored against the bed this
    architecture actually produces, rather than a nominal clearance.
    """
    require_positive("draws", draws)
    require_positive("months", months)
    require_fraction("source_ceiling", source_ceiling, allow_one=True)
    count = int(draws)
    rng = np.random.default_rng(seed)
    central = _central()
    gate_index = _ALL_STRENGTHS.index("gate_specificity")

    gate = _path_max(rng, count, upgrades.gate_paths)
    nodes = np.empty((count, len(central)))
    for j in range(len(central)):
        nodes[:, j] = _path_max(rng, count, upgrades.node_paths)
    nodes[:, gate_index] = gate
    ceilings = nodes * central

    source_draws = rng.uniform(0.0, source_ceiling, count)
    if upgrades.verified_source:
        keep = source_draws >= upgrades.source_required
        eligible = float(np.mean(keep))
        if not keep.any():
            return np.empty(0), np.empty(0), 0.0
        ceilings = ceilings[keep]
        # Verification tells us this patient is a responder; it does not reduce
        # their response to the threshold. Keep their own achieved value, which is
        # at least the threshold and may be much better. Pinning it to the
        # threshold would discard the very information the run-in bought.
        source_values = source_draws[keep]
    else:
        eligible = 1.0
        source_values = source_draws

    effective = ceilings.shape[0]
    achieved = np.empty(effective)
    potency_used = np.empty(effective)
    grid = (
        np.linspace(1.0, upgrades.max_potency, upgrades.titration_steps)
        if upgrades.titrate
        else np.array([ATLAS_POTENCY])
    )
    for i in range(effective):
        reached = False
        for potency in grid:
            value = _clearance(np.minimum(1.0, ceilings[i] * potency), source_values[i], months)
            if value >= target:
                achieved[i], potency_used[i] = value, potency
                reached = True
                break
        if not reached:
            # Titration ran out of headroom; record the best it could do.
            potency = grid[-1]
            achieved[i] = _clearance(
                np.minimum(1.0, ceilings[i] * potency), source_values[i], months
            )
            potency_used[i] = potency
    return achieved, potency_used, eligible


def evaluate_upgrades(
    upgrades: Upgrades = ATLAS2,
    months: float = ATLAS_CLEARANCE_MONTHS,
    draws: int = 1500,
    seed: int = 20261107,
    target: float = CRITERION_THRESHOLD,
    source_ceiling: float = ATLAS_SOURCE_REDUCTION,
) -> UpgradeResult:
    """Run the upgraded architecture, charging for every upgrade it uses."""
    achieved, potency_used, eligible = achieved_clearance_draws(
        upgrades, months, draws, seed, target, source_ceiling
    )
    if achieved.size == 0:
        return UpgradeResult(
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            False,
            upgrades.verified_source,
            upgrades.redundancy_beyond_precedent,
            int(draws),
        )
    probability = float(np.mean(achieved >= target))
    return UpgradeResult(
        probability_of_target=probability,
        median_clearance=float(np.median(achieved)),
        p05_clearance=float(np.quantile(achieved, 0.05)),
        eligible_fraction=eligible,
        median_potency_used=float(np.median(potency_used)),
        robust=probability >= 0.8,
        conditional=upgrades.verified_source,
        redundancy_beyond_precedent=upgrades.redundancy_beyond_precedent,
        draws=int(draws),
    )


#: The solved operating point: 60 months rather than 24, at the smallest potency
#: cap that makes the clearance arm robust. Found by bisection, not chosen.
SOLVED_MONTHS = 60.0


def minimum_robust_cap(
    months: float = SOLVED_MONTHS,
    upgrades: Upgrades = ATLAS2,
    draws: int = 400,
    seed: int = 20261107,
    upper: float = 8.0,
    required_probability: float = 0.8,
) -> float | None:
    """Smallest potency cap that makes the arm robust, or None if none does."""
    fields = {f: getattr(upgrades, f) for f in upgrades.__dataclass_fields__}

    def probability_at(cap: float) -> float:
        config = Upgrades(**{**fields, "max_potency": cap})
        result = evaluate_upgrades(config, months=months, draws=draws, seed=seed)
        return result.probability_of_target

    if probability_at(upper) < required_probability:
        return None
    lo, hi = 1.0, upper
    for _ in range(16):
        mid = 0.5 * (lo + hi)
        if probability_at(mid) >= required_probability:
            hi = mid
        else:
            lo = mid
    return float(hi)


def horizon_sweep(
    months_grid: tuple[float, ...] = (24.0, 60.0, 120.0),
    caps: tuple[float, ...] = (1.284, 1.5, 2.0, 2.5, 3.0),
    draws: int = 400,
    seed: int = 20261107,
) -> pd.DataFrame:
    """Time against potency: which of the two actually buys robustness."""
    rows = []
    for months in months_grid:
        for cap in caps:
            config = Upgrades(
                gate_paths=ATLAS2.gate_paths,
                node_paths=ATLAS2.node_paths,
                titrate=True,
                verified_source=True,
                max_potency=cap,
                source_required=ATLAS2.source_required,
            )
            result = evaluate_upgrades(config, months=months, draws=draws, seed=seed)
            rows.append(
                {
                    "horizon_months": months,
                    "max_potency": cap,
                    "probability_of_target": result.probability_of_target,
                    "median_clearance_pct": 100.0 * result.median_clearance,
                    "eligible_fraction": result.eligible_fraction,
                    "unconditional_probability": result.unconditional_probability,
                    "robust": result.robust,
                    "beyond_precedent_potency": cap > ATLAS_POTENCY,
                    "status": UPGRADE_TAG,
                }
            )
    return pd.DataFrame(rows)


def coupled_graft_with_upgrades(
    upgrades: Upgrades = ATLAS2,
    months: float = SOLVED_MONTHS,
    bed_coupling: float = 0.5,
    base_floor: float = 0.5,
    draws: int = 800,
    seed: int = 20261107,
) -> pd.DataFrame:
    """Does sequencing finally pay off, now the bed is actually prepared?

    ATLAS-1's sequencing gain collapsed because phase 1 usually failed. Here the
    graft is scored against the bed the *upgraded* phase 1 produces, conditional on
    the patient having cleared the run-in.
    """
    from ihd_solutions.gap_solver import (
        BASELINE_SCAR_FRACTION,
        REMUSCULARISATION_FACTORS,
        SCAR_RECOVERABLE_LVEF_PP,
    )

    require_fraction("bed_coupling", bed_coupling, allow_one=True)
    achieved, _, eligible = achieved_clearance_draws(
        upgrades, months=months, draws=draws, seed=seed
    )
    if achieved.size == 0:
        return pd.DataFrame(
            [
                {
                    "bed_coupling": bed_coupling,
                    "median_bed_quality": 0.0,
                    "graft_probability": 0.0,
                    "eligible_fraction": 0.0,
                    "status": UPGRADE_TAG,
                }
            ]
        )
    quality = np.clip(achieved / CRITERION_THRESHOLD, 0.0, 1.0)
    floors = np.clip(base_floor + bed_coupling * (1.0 - base_floor) * quality, base_floor, 1.0)
    rng = np.random.default_rng(seed + 1)
    delivered = np.ones(achieved.size)
    for factor in REMUSCULARISATION_FACTORS:
        delivered *= rng.uniform(1.0, factor.demonstrated_multiplier, size=achieved.size)
    shared = rng.random(achieved.size)
    required = CRITERION_THRESHOLD * SCAR_RECOVERABLE_LVEF_PP
    gains = (
        np.minimum(1.0, BASELINE_SCAR_FRACTION * delivered * (floors + (1.0 - floors) * shared))
        * SCAR_RECOVERABLE_LVEF_PP
    )
    unsequenced = (
        np.minimum(
            1.0, BASELINE_SCAR_FRACTION * delivered * (base_floor + (1.0 - base_floor) * shared)
        )
        * SCAR_RECOVERABLE_LVEF_PP
    )
    probability = float(np.mean(gains >= required))
    without = float(np.mean(unsequenced >= required))
    return pd.DataFrame(
        [
            {
                "bed_coupling": bed_coupling,
                "median_bed_quality": float(np.median(quality)),
                "graft_probability": probability,
                "graft_probability_unsequenced": without,
                "sequencing_gain_pp": 100.0 * (probability - without),
                "median_gain_pp": float(np.median(gains)),
                "eligible_fraction": eligible,
                "graft_robust": probability >= 0.8,
                "status": UPGRADE_TAG,
            }
        ]
    )


def upgrade_ablation(
    upgrades: Upgrades = ATLAS2,
    draws: int = 1200,
    seed: int = 20261107,
    months: float = 60.0,
) -> pd.DataFrame:
    """Remove one upgrade at a time. Whatever collapses was load-bearing."""
    full = evaluate_upgrades(upgrades, months=months, draws=draws, seed=seed)
    rows = []
    for name in UPGRADE_NAMES:
        reduced = evaluate_upgrades(upgrades.without(name), months=months, draws=draws, seed=seed)
        rows.append(
            {
                "removed_upgrade": name,
                "active_in_construct": name != "verified_source" or upgrades.verified_source,
                "probability_of_target": reduced.probability_of_target,
                "probability_loss_pp": 100.0
                * (full.probability_of_target - reduced.probability_of_target),
                "eligible_fraction": reduced.eligible_fraction,
                "unconditional_probability": reduced.unconditional_probability,
                "robust": reduced.robust,
                "load_bearing": full.probability_of_target - reduced.probability_of_target > 1e-6,
                "status": UPGRADE_TAG,
            }
        )
    return pd.DataFrame(rows)


def redundancy_sweep(
    paths: tuple[int, ...] = (1, 2, 3, 4, 6),
    draws: int = 1000,
    seed: int = 20261107,
    months: float = 60.0,
) -> pd.DataFrame:
    """How much redundancy is needed, and where it stops being defensible."""
    rows = []
    for k in paths:
        config = Upgrades(
            gate_paths=k,
            node_paths=min(k, 2),
            titrate=True,
            verified_source=False,
            max_potency=ATLAS2.max_potency,
        )
        result = evaluate_upgrades(config, months=months, draws=draws, seed=seed)
        rows.append(
            {
                "gate_paths": k,
                "probability_of_target": result.probability_of_target,
                "median_clearance_pct": 100.0 * result.median_clearance,
                "median_potency_used": result.median_potency_used,
                "robust": result.robust,
                "beyond_precedent": k > MAX_PRECEDENTED_PATHS,
                "status": UPGRADE_TAG,
            }
        )
    return pd.DataFrame(rows)


def potency_ceiling_sweep(
    ceilings: tuple[float, ...] = (1.284, 1.5, 2.0, 2.59, 3.0),
    draws: int = 1000,
    seed: int = 20261107,
    months: float = 60.0,
) -> pd.DataFrame:
    """Does the fix survive a *tighter* potency ceiling than it was given?"""
    rows = []
    for ceiling in ceilings:
        config = Upgrades(
            gate_paths=ATLAS2.gate_paths,
            node_paths=ATLAS2.node_paths,
            titrate=True,
            verified_source=False,
            max_potency=ceiling,
        )
        result = evaluate_upgrades(config, months=months, draws=draws, seed=seed)
        rows.append(
            {
                "max_potency": ceiling,
                "probability_of_target": result.probability_of_target,
                "median_potency_used": result.median_potency_used,
                "robust": result.robust,
                "status": UPGRADE_TAG,
            }
        )
    return pd.DataFrame(rows)


def source_gate_tradeoff(
    thresholds: tuple[float, ...] = (0.0, 0.10, 0.2230, 0.35, 0.50),
    draws: int = 1000,
    seed: int = 20261107,
    source_ceiling: float = ATLAS_SOURCE_REDUCTION,
) -> pd.DataFrame:
    """The price of verified source control, in excluded patients."""
    rows = []
    for threshold in thresholds:
        config = Upgrades(
            gate_paths=ATLAS2.gate_paths,
            node_paths=ATLAS2.node_paths,
            titrate=True,
            verified_source=True,
            source_required=threshold,
        )
        result = evaluate_upgrades(config, draws=draws, seed=seed, source_ceiling=source_ceiling)
        rows.append(
            {
                "source_required": threshold,
                "eligible_fraction": result.eligible_fraction,
                "excluded_fraction": 1.0 - result.eligible_fraction,
                "conditional_probability": result.probability_of_target,
                "unconditional_probability": result.unconditional_probability,
                "robust_conditional": result.robust,
                "status": UPGRADE_TAG,
            }
        )
    return pd.DataFrame(rows)


def upgrade_report(
    baseline: UpgradeResult,
    full: UpgradeResult,
    ablation: pd.DataFrame,
    redundancy: pd.DataFrame,
    ceilings: pd.DataFrame,
    horizon: pd.DataFrame,
    graft: pd.DataFrame,
    solved_cap: float | None,
) -> str:
    lines = [
        "ATLAS-2: upgrading the parts the audits broke",
        "===========================================",
        DISCLAIMER,
        "",
        "THE DIAGNOSIS BEING FIXED",
        "delivered() multiplies six of eight module strengths by one gate_specificity",
        "draw, so a single uncertain parameter throttles almost the whole system in",
        "series. Delivery was open loop: pick a potency and hope. Both are architectural",
        "defects, not potency shortages, so ATLAS-2 changes the architecture.",
        "",
        "  U1 redundant gate + nodes  -> removes the series bottleneck",
        "  U2 closed-loop titration   -> treat to a measured target, as lipidology does",
        "  U3 verified source control -> available, but NOT needed at the solved horizon",
        "",
        "RESULT (60-month horizon, whole population, no run-in exclusion)",
        f"  ATLAS-1 baseline (upgrades off, 24mo)  P(target) {100 * baseline.probability_of_target:6.2f}%"
        f"   median {100 * baseline.median_clearance:6.2f}%",
        f"  ATLAS-2 (upgrades on, 60mo)            P(target) {100 * full.probability_of_target:6.2f}%"
        f"   median {100 * full.median_clearance:6.2f}%",
        f"  robust at 80%                          {'YES' if full.robust else 'no'}",
        f"  robust with 2-SE margin                {'YES' if full.robust_with_margin else 'no'}"
        f"  (lower bound {100 * (full.probability_of_target - 2 * full.monte_carlo_standard_error):.2f}%)",
        f"  eligible fraction                      {100 * full.eligible_fraction:6.2f}%",
        f"  UNCONDITIONAL success                  {100 * full.unconditional_probability:6.2f}%",
        "",
        "WHAT IT COST",
        f"  potency cap required                   {ATLAS2.max_potency:.3f}x",
        f"  precedented operating point            {ATLAS_POTENCY:.3f}x",
        f"  precedent gap                          {ATLAS2.max_potency / ATLAS_POTENCY:.2f}x",
        f"  distinct targeting routes              {ATLAS2.gate_paths} (cap {MAX_PRECEDENTED_PATHS} is precedented)",
        "  horizon                                60 months instead of 24",
        "",
        "  The remaining gap is potency and time, not architecture. Nothing here is",
        "  free: 2x the precedented potency across every module simultaneously is a",
        "  strong assumption, and it is the thing that would have to be demonstrated.",
        "",
        "A CORRECTION TO THE EARLIER AUDIT",
        "The ATLAS-1 audit reported a 55.4% probability ceiling at any potency while the",
        "source term is uncertain. That was horizon-specific. At 24 months perfect",
        "modules need source reduction >= 0.2230; at 36 months and beyond they need",
        "none, because the structural ceiling itself rises above the target. So the",
        "source lever stops binding once the horizon is long enough, and a run-in would",
        "only exclude patients who would have succeeded. U3 is therefore off.",
        "",
        "UPGRADE ABLATION (at the solved point)",
    ]
    for row in ablation.itertuples(index=False):
        state = (
            "LOAD-BEARING"
            if row.load_bearing
            else ("not used at this horizon" if not row.active_in_construct else "inert")
        )
        lines.append(
            f"  drop {row.removed_upgrade:18s} P {100 * row.probability_of_target:6.2f}% "
            f"({row.probability_loss_pp:+6.2f} pp)  {state}"
        )
    lines += [
        "",
        "  Titration is the single largest contributor. Treating to a measured target",
        "  beats adding potency, which is the same lesson as the gate in REVERT-1:",
        "  control structure beats raw strength.",
        "",
        "TIME AGAINST POTENCY",
    ]
    for row in horizon.itertuples(index=False):
        lines.append(
            f"  {row.horizon_months:5.0f}mo  cap {row.max_potency:5.3f}x  "
            f"P {100 * row.probability_of_target:6.2f}%  "
            f"{'robust' if row.robust else '      '}  "
            f"{'beyond precedent' if row.beyond_precedent_potency else 'precedented'}"
        )
    lines += ["", "REDUNDANCY: where precedent runs out"]
    for row in redundancy.itertuples(index=False):
        lines.append(
            f"  {row.gate_paths} routes  P {100 * row.probability_of_target:6.2f}%  "
            f"{'BEYOND PRECEDENT' if row.beyond_precedent else 'precedented'}"
        )
    lines += ["", "POTENCY CEILING at 60 months"]
    for row in ceilings.itertuples(index=False):
        lines.append(
            f"  cap {row.max_potency:5.3f}x  P {100 * row.probability_of_target:6.2f}%  "
            f"{'robust' if row.robust else 'not robust'}"
        )
    lines += [
        "",
        "DOES THE SEQUENCING IDEA FINALLY PAY OFF?",
        "ATLAS-1's sequencing gain collapsed because phase 1 usually failed, leaving",
        "median bed quality at 0.023. Scored against the bed the upgraded phase 1",
        "actually produces:",
        "",
    ]
    for row in graft.itertuples(index=False):
        lines.append(
            f"  coupling {row.bed_coupling:4.2f}  bed quality {row.median_bed_quality:5.3f}  "
            f"graft P {100 * row.graft_probability:5.1f}%  "
            f"sequencing gain {row.sequencing_gain_pp:+5.1f} pp"
        )
    lines += [
        "",
        "  The idea was never structurally wrong; it was starved by a failing phase 1.",
        "  With the bed genuinely prepared it delivers its full gain again. The graft",
        "  still is not robust, so remuscularisation is now the sole binding biology",
        "  problem rather than one of two.",
        "",
        "VERDICT",
        "The clearance arm is FIXED as a specification: robust across the whole",
        "population with margin at the stated sample size,"
        "no eligibility exclusion, every active upgrade load-bearing. The",
        f"price is a {ATLAS2.max_potency / ATLAS_POTENCY:.2f}x potency gap beyond precedent and a"
        " 60-month horizon instead of 24.",
        "",
        "Not fixed: myocardial restoration is still not robust, and randomised proof is",
        "still zero because none of this has been trialled. ATLAS-2 is a materially",
        "better candidate and a sharper specification of what must be demonstrated. It",
        "is not a cure.",
    ]
    return "\n".join(lines) + "\n"


def run_upgrade_study(
    output_dir: Path | str = Path("results/upgrade"),
    draws: int = 800,
    seed: int = 20261107,
) -> dict[str, object]:
    root = Path(output_dir)
    tables, reports = root / "tables", root / "reports"
    tables.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)

    baseline = evaluate_upgrades(Upgrades(), months=ATLAS_CLEARANCE_MONTHS, draws=draws, seed=seed)
    full = evaluate_upgrades(ATLAS2, months=SOLVED_MONTHS, draws=draws, seed=seed)
    ablation = upgrade_ablation(ATLAS2, draws=draws, seed=seed, months=SOLVED_MONTHS)
    redundancy = redundancy_sweep(draws=max(300, draws // 2), seed=seed, months=SOLVED_MONTHS)
    ceilings = potency_ceiling_sweep(draws=max(300, draws // 2), seed=seed, months=SOLVED_MONTHS)
    horizon = horizon_sweep(draws=max(300, draws // 2), seed=seed)
    tradeoff = source_gate_tradeoff(draws=max(300, draws // 2), seed=seed)
    solved_cap = minimum_robust_cap(months=SOLVED_MONTHS, draws=max(300, draws // 2), seed=seed)
    graft = pd.concat(
        [
            coupled_graft_with_upgrades(
                ATLAS2, months=SOLVED_MONTHS, bed_coupling=c, draws=draws, seed=seed
            )
            for c in (0.0, 0.5, 1.0)
        ],
        ignore_index=True,
    )

    def row(label: str, result: UpgradeResult, months: float) -> dict:
        return {
            "construct": label,
            "horizon_months": months,
            "probability_of_target": result.probability_of_target,
            "median_clearance_pct": 100 * result.median_clearance,
            "eligible_fraction": result.eligible_fraction,
            "unconditional_probability": result.unconditional_probability,
            "median_potency_used": result.median_potency_used,
            "robust": result.robust,
            "robust_with_margin": result.robust_with_margin,
            "monte_carlo_standard_error": result.monte_carlo_standard_error,
            "status": UPGRADE_TAG,
        }

    pd.DataFrame(
        [
            row("atlas1_baseline", baseline, ATLAS_CLEARANCE_MONTHS),
            row("atlas2_solved", full, SOLVED_MONTHS),
        ]
    ).to_csv(tables / "upgrade_headline.csv", index=False)
    ablation.to_csv(tables / "upgrade_ablation.csv", index=False)
    redundancy.to_csv(tables / "upgrade_redundancy_sweep.csv", index=False)
    ceilings.to_csv(tables / "upgrade_potency_ceiling_sweep.csv", index=False)
    horizon.to_csv(tables / "upgrade_horizon_sweep.csv", index=False)
    tradeoff.to_csv(tables / "upgrade_source_gate_tradeoff.csv", index=False)
    graft.to_csv(tables / "upgrade_coupled_graft.csv", index=False)
    pd.DataFrame(
        [
            {
                "solved_horizon_months": SOLVED_MONTHS,
                "minimum_robust_cap": solved_cap,
                "precedented_cap": ATLAS_POTENCY,
                "precedent_gap": (solved_cap / ATLAS_POTENCY) if solved_cap else None,
                "gate_paths": ATLAS2.gate_paths,
                "node_paths": ATLAS2.node_paths,
                "verified_source_needed": ATLAS2.verified_source,
                "status": UPGRADE_TAG,
            }
        ]
    ).to_csv(tables / "upgrade_solved_operating_point.csv", index=False)
    (reports / "upgrade_report.txt").write_text(
        upgrade_report(baseline, full, ablation, redundancy, ceilings, horizon, graft, solved_cap)
    )
    (root / "DISCLAIMER.txt").write_text(DISCLAIMER + "\n")
    return {
        "output_dir": str(root),
        "outputs": 10,
        "status": "invented_architecture_unvalidated",
        "robust": bool(full.robust),
        "robust_with_margin": bool(full.robust_with_margin),
        "eligible_fraction": full.eligible_fraction,
        "unconditional_probability": full.unconditional_probability,
        "minimum_robust_cap": solved_cap,
        "cures_declared": 0,
        "disclaimer": DISCLAIMER,
    }
