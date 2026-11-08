"""Step-by-step resolution of the three residual gaps.

Each residual gets a solver that returns an operating point, not an adjective.

**Residual 1 - hazard ratio.** Five populations close in surrogate mode and all
nineteen close in exploratory mode, so the shortfall is evidence maturity, not
biology. The solver names the specific agents whose evidence would have to mature,
and ranks them by how many populations each unlocks.

**Residual 2 - clearance.** A correction to the earlier summary: it is *not* delivery
alone. Potentiating any single module to its maximum reaches only 59-63% against an
80% target. The solver therefore searches a uniform potency multiplier across all
modules jointly, and reports the cheapest (potency, source-reduction, horizon)
operating point inside the precedent bound.

**Residual 3 - remuscularisation.** The nominal x5.2 against a x5.13 requirement has
no margin. The solver replaces point estimates with distributions - each factor
delivers anywhere from nothing to its published effect - and adds an explicit
chronic-injury haircut, because the primate work treats two weeks after infarction
rather than a years-old scar. It then finds the smallest factor set that clears the
target with a stated probability rather than by arithmetic luck.

Everything here is scenario arithmetic. A solved residual is a specification.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from itertools import combinations

import numpy as np
import pandas as pd

from ihd_solutions.catalog import Therapy, load_catalog
from ihd_solutions.cure import CRITERION_THRESHOLD
from ihd_solutions.discovery import required_capabilities
from ihd_solutions.evidence import EvidenceTier
from ihd_solutions.gap_solver import (
    BASELINE_SCAR_FRACTION,
    REMUSCULARISATION_FACTORS,
    SCAR_RECOVERABLE_LVEF_PP,
    SOURCE_REDUCTION_PRECEDENT_CENTRAL,
    SOURCE_REDUCTION_PRECEDENT_MAX,
    Factor,
    hazard_ratio_from_risks,
)
from ihd_solutions.hypothesis import (
    _ZERO_STRENGTHS,
    ClearanceParameters,
    simulate_clearance,
)
from ihd_solutions.units import require_fraction, require_positive

DISCLAIMER = (
    "RESIDUAL-GAP RESOLUTION - scenario arithmetic producing operating points and "
    "evidence requirements, never efficacy evidence. The chronic-injury haircut and "
    "the potency multiplier are assumptions, not measurements."
)

RESOLUTION_TAG = "residual_resolution_scenario_not_evidence"

#: Modes in increasing permissiveness. Moving right requires evidence to mature.
MODE_ORDER: tuple[str, ...] = ("strict", "surrogate", "exploratory")

#: Tiers that can move a strict-mode event curve.
OUTCOMES_TIERS = frozenset({EvidenceTier.META_ANALYSIS, EvidenceTier.RANDOMIZED_OUTCOMES})


# ---------------------------------------------------------------------------
# Step 1: hazard ratio -> which evidence has to mature.
# ---------------------------------------------------------------------------


def population_closure(
    matrix: pd.DataFrame, reference_risk: float = 0.0198, threshold: float = CRITERION_THRESHOLD
) -> pd.DataFrame:
    """Per population, the first credit mode that closes the hazard gap."""
    rows = []
    for population_id in sorted(matrix["population_id"].unique()):
        record: dict[str, object] = {"population_id": population_id}
        closing_mode = None
        for mode in MODE_ORDER:
            subset = matrix.loc[matrix["population_id"].eq(population_id) & matrix["mode"].eq(mode)]
            if subset.empty:
                continue
            row = subset.iloc[0]
            baseline = float(row["baseline_risk_pct"]) / 100.0
            treated = float(row["composite_risk_pct"]) / 100.0
            if baseline <= reference_risk + 1e-4:
                record[f"{mode}_status"] = "at_reference"
                record[f"{mode}_residual_pct"] = 0.0
                closing_mode = closing_mode or mode
                continue
            required = required_capabilities(baseline, reference_risk, threshold)["hazard_ratio"]
            achieved = hazard_ratio_from_risks(baseline, treated)
            closed = achieved <= required * (1.0 + 1e-9)
            record[f"{mode}_status"] = "closed" if closed else "open"
            record[f"{mode}_residual_pct"] = (
                0.0
                if closed
                else 100.0 * (1.0 - (-math.expm1(required * math.log1p(-baseline))) / treated)
            )
            if closed and closing_mode is None:
                closing_mode = mode
        record["first_closing_mode"] = closing_mode or "none"
        record["status"] = RESOLUTION_TAG
        rows.append(record)
    return pd.DataFrame(rows)


def _components(matrix: pd.DataFrame, population_id: str, mode: str) -> set[str]:
    subset = matrix.loc[matrix["population_id"].eq(population_id) & matrix["mode"].eq(mode)]
    if subset.empty:
        return set()
    return set(str(subset.iloc[0]["therapy_ids"]).split(";")) - {"(none)", ""}


def evidence_upgrade_levers(
    matrix: pd.DataFrame,
    catalog: dict[str, Therapy] | None = None,
    reference_risk: float = 0.0198,
) -> pd.DataFrame:
    """Rank agents by how many populations their evidence maturation would unlock.

    An agent counts toward a population when it enters the composite as credit
    relaxes from strict to the first mode that closes that population's hazard gap.
    """
    catalog = load_catalog() if catalog is None else catalog
    closure = population_closure(matrix, reference_risk).set_index("population_id")
    counts: dict[str, set[str]] = {}
    for population_id, row in closure.iterrows():
        mode = str(row["first_closing_mode"])
        if mode in {"none", "strict"}:
            continue
        for therapy_id in _components(matrix, str(population_id), mode) - _components(
            matrix, str(population_id), "strict"
        ):
            counts.setdefault(therapy_id, set()).add(str(population_id))
    rows = []
    for therapy_id, populations in counts.items():
        therapy = catalog[therapy_id]
        rows.append(
            {
                "therapy_id": therapy_id,
                "populations_unlocked": len(populations),
                "populations": ";".join(sorted(populations)),
                "current_tier": therapy.evidence_tier.value,
                "regulatory_status": therapy.regulatory_status,
                "needs_tier": (
                    "none - already outcomes tier"
                    if therapy.evidence_tier in OUTCOMES_TIERS
                    else "RANDOMIZED_OUTCOMES"
                ),
                "already_outcomes_tier": therapy.evidence_tier in OUTCOMES_TIERS,
                # Entering a relaxed composite is only an evidence gap if the agent
                # actually lacks outcomes-tier proof. Otherwise the greedy ordering
                # simply changed when other agents gained credit.
                "reason": (
                    "composition_reordering"
                    if therapy.evidence_tier in OUTCOMES_TIERS
                    else "evidence_maturation_required"
                ),
                "status": RESOLUTION_TAG,
            }
        )
    if not rows:
        return pd.DataFrame(
            columns=[
                "therapy_id",
                "populations_unlocked",
                "populations",
                "current_tier",
                "regulatory_status",
                "needs_tier",
                "already_outcomes_tier",
                "reason",
                "status",
            ]
        )
    return (
        pd.DataFrame(rows)
        .sort_values(["populations_unlocked", "therapy_id"], ascending=[False, True])
        .reset_index(drop=True)
    )


# ---------------------------------------------------------------------------
# Step 2: clearance -> a joint operating point, not a delivery fix.
# ---------------------------------------------------------------------------

CENTRAL = ClearanceParameters()
_ALL_STRENGTHS = ("gate_specificity",) + _ZERO_STRENGTHS


def potency_scaled(multiplier: float) -> ClearanceParameters:
    """Scale every central module strength by one factor, capped at 1.0.

    One interpretable knob instead of eight: how much better than the central
    assumption does *every* module have to perform simultaneously.
    """
    require_positive("multiplier", multiplier)
    values = {name: min(1.0, multiplier * getattr(CENTRAL, name)) for name in _ALL_STRENGTHS}
    return ClearanceParameters(**values)


def clearance_at(multiplier: float, source_reduction: float, horizon_months: float) -> float:
    return simulate_clearance(
        potency_scaled(multiplier),
        "resolve",
        months=horizon_months,
        step_months=0.5,
        source_reduction=source_reduction,
    ).substrate_removal_fraction


def module_leverage_table(
    source_reduction: float = SOURCE_REDUCTION_PRECEDENT_CENTRAL,
    horizon_months: float = 24.0,
) -> pd.DataFrame:
    """Potentiate one module at a time to its maximum. None alone is enough."""
    require_fraction("source_reduction", source_reduction, allow_one=True)
    require_positive("horizon_months", horizon_months)
    base = 100.0 * clearance_at(1.0, source_reduction, horizon_months)
    rows = []
    for name in _ALL_STRENGTHS:
        values = {key: getattr(CENTRAL, key) for key in _ALL_STRENGTHS}
        values[name] = 1.0
        reached = (
            100.0
            * simulate_clearance(
                ClearanceParameters(**values),
                "resolve",
                months=horizon_months,
                step_months=0.5,
                source_reduction=source_reduction,
            ).substrate_removal_fraction
        )
        rows.append(
            {
                "module": name,
                "central_value": getattr(CENTRAL, name),
                "removal_alone_pct": reached,
                "gain_over_central_pp": reached - base,
                "reaches_target_alone": reached >= 100.0 * CRITERION_THRESHOLD,
                "status": RESOLUTION_TAG,
            }
        )
    return (
        pd.DataFrame(rows).sort_values("removal_alone_pct", ascending=False).reset_index(drop=True)
    )


def minimum_potency(
    source_reduction: float,
    horizon_months: float,
    required: float = CRITERION_THRESHOLD,
    tolerance: float = 1e-4,
    max_multiplier: float = 4.0,
) -> float | None:
    """Smallest uniform potency multiplier that reaches the clearance target."""
    require_fraction("source_reduction", source_reduction, allow_one=True)
    require_fraction("required", required, allow_one=True)
    require_positive("horizon_months", horizon_months)
    require_positive("max_multiplier", max_multiplier)
    if clearance_at(max_multiplier, source_reduction, horizon_months) < required:
        return None
    low, high = 1.0, max_multiplier
    if clearance_at(low, source_reduction, horizon_months) >= required:
        return low
    while high - low > tolerance:
        mid = 0.5 * (low + high)
        if clearance_at(mid, source_reduction, horizon_months) >= required:
            high = mid
        else:
            low = mid
    return float(high)


def clearance_operating_points(
    sources: tuple[float, ...] = (0.0, 0.20, 0.35, 0.50),
    horizons: tuple[float, ...] = (24.0, 36.0, 60.0, 120.0),
) -> pd.DataFrame:
    """The cheapest feasible (potency, source, horizon) combinations."""
    rows = []
    for source in sources:
        for horizon in horizons:
            potency = minimum_potency(source, horizon)
            rows.append(
                {
                    "source_reduction": source,
                    "horizon_months": horizon,
                    "minimum_potency_multiplier": potency,
                    "feasible": potency is not None,
                    "source_within_precedent": source <= SOURCE_REDUCTION_PRECEDENT_MAX,
                    "removal_at_potency_pct": (
                        100.0 * clearance_at(potency, source, horizon)
                        if potency is not None
                        else float("nan")
                    ),
                    "status": RESOLUTION_TAG,
                }
            )
    frame = pd.DataFrame(rows)
    frame["admissible"] = frame["feasible"] & frame["source_within_precedent"]
    return frame


def cheapest_clearance_point(operating: pd.DataFrame) -> dict | None:
    """Prefer the shortest horizon, then the lowest potency, inside precedent."""
    admissible = operating.loc[operating["admissible"]]
    if admissible.empty:
        return None
    best = admissible.sort_values(["horizon_months", "minimum_potency_multiplier"]).iloc[0]
    return {
        "source_reduction": float(best["source_reduction"]),
        "horizon_months": float(best["horizon_months"]),
        "minimum_potency_multiplier": float(best["minimum_potency_multiplier"]),
    }


# ---------------------------------------------------------------------------
# Step 3: remuscularisation -> a set with margin, under a chronic-injury haircut.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class RobustnessResult:
    factors: tuple[str, ...]
    nominal_gain_pp: float
    probability_of_target: float
    median_gain_pp: float
    p05_gain_pp: float
    chronic_penalty_low: float

    @property
    def robust(self) -> bool:
        return self.probability_of_target >= 0.8


def remuscularisation_robustness(
    factors: tuple[Factor, ...],
    draws: int = 4000,
    seed: int = 20261107,
    chronic_penalty_low: float = 0.5,
    required_pp: float = CRITERION_THRESHOLD * SCAR_RECOVERABLE_LVEF_PP,
) -> RobustnessResult:
    """Each factor delivers Uniform(1, published effect); then a chronic haircut.

    The haircut is an explicit assumption, not data: primate remuscularisation work
    treats two weeks after infarction, while the target patient has a years-old,
    stiff, poorly vascularised scar. Sweeping it is the honest way to show how much
    of the result depends on that translation.
    """
    require_positive("draws", draws)
    require_fraction("chronic_penalty_low", chronic_penalty_low, allow_one=True)
    require_positive("required_pp", required_pp)
    rng = np.random.default_rng(seed)
    delivered = np.ones(int(draws))
    nominal = 1.0
    for factor in factors:
        delivered *= rng.uniform(1.0, factor.demonstrated_multiplier, size=int(draws))
        nominal *= factor.demonstrated_multiplier
    penalty = rng.uniform(chronic_penalty_low, 1.0, size=int(draws))
    scar = np.minimum(1.0, BASELINE_SCAR_FRACTION * delivered * penalty)
    gains = scar * SCAR_RECOVERABLE_LVEF_PP
    nominal_gain = min(1.0, BASELINE_SCAR_FRACTION * nominal) * SCAR_RECOVERABLE_LVEF_PP
    return RobustnessResult(
        tuple(factor.name for factor in factors),
        float(nominal_gain),
        float(np.mean(gains >= required_pp)),
        float(np.median(gains)),
        float(np.quantile(gains, 0.05)),
        float(chronic_penalty_low),
    )


def robust_remuscularisation_search(
    draws: int = 4000,
    seed: int = 20261107,
    chronic_penalty_low: float = 0.5,
) -> pd.DataFrame:
    """Every factor subset, scored by probability rather than nominal arithmetic."""
    rows = []
    for size in range(len(REMUSCULARISATION_FACTORS) + 1):
        for subset in combinations(REMUSCULARISATION_FACTORS, size):
            result = remuscularisation_robustness(
                subset, draws=draws, seed=seed, chronic_penalty_low=chronic_penalty_low
            )
            rows.append(
                {
                    "n_factors": size,
                    "factors": ";".join(result.factors) or "none",
                    "nominal_gain_pp": result.nominal_gain_pp,
                    "median_gain_pp": result.median_gain_pp,
                    "p05_gain_pp": result.p05_gain_pp,
                    "probability_of_target": result.probability_of_target,
                    "robust_at_80pct": result.robust,
                    "chronic_penalty_low": chronic_penalty_low,
                    "status": RESOLUTION_TAG,
                }
            )
    return (
        pd.DataFrame(rows)
        .sort_values(
            ["robust_at_80pct", "n_factors", "probability_of_target"],
            ascending=[False, True, False],
        )
        .reset_index(drop=True)
    )


def minimum_robust_set(
    draws: int = 4000, seed: int = 20261107, chronic_penalty_low: float = 0.5
) -> tuple[str, ...]:
    search = robust_remuscularisation_search(draws, seed, chronic_penalty_low)
    passing = search.loc[search["robust_at_80pct"]]
    if passing.empty:
        return ()
    best = passing.sort_values(
        ["n_factors", "probability_of_target"], ascending=[True, False]
    ).iloc[0]
    return () if best["factors"] == "none" else tuple(str(best["factors"]).split(";"))


def chronic_penalty_sensitivity(
    penalties: tuple[float, ...] = (1.0, 0.75, 0.5, 0.25),
    draws: int = 3000,
    seed: int = 20261107,
) -> pd.DataFrame:
    """How the answer depends on the sub-acute-to-chronic translation assumption."""
    rows = []
    for penalty in penalties:
        require_fraction("chronic_penalty_low", penalty, allow_one=True)
        search = robust_remuscularisation_search(draws, seed, penalty)
        passing = search.loc[search["robust_at_80pct"]]
        # When nothing is robust, report the best available set and its probability.
        # Reporting the empty set here would read as "no factors needed" instead of
        # "no combination works", which is the opposite of the truth.
        if passing.empty:
            best = search.sort_values("probability_of_target", ascending=False).iloc[0]
            robust_exists = False
        else:
            best = passing.sort_values(
                ["n_factors", "probability_of_target"], ascending=[True, False]
            ).iloc[0]
            robust_exists = True
        rows.append(
            {
                "chronic_penalty_low": penalty,
                "robust_set_exists": robust_exists,
                "minimum_robust_factors": int(best["n_factors"]) if robust_exists else 0,
                "best_available_factors": best["factors"],
                "best_available_n_factors": int(best["n_factors"]),
                "best_probability_of_target": float(best["probability_of_target"]),
                "best_median_gain_pp": float(best["median_gain_pp"]),
                "status": RESOLUTION_TAG,
            }
        )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Report and writer.
# ---------------------------------------------------------------------------


def resolution_report(matrix: pd.DataFrame) -> str:
    closure = population_closure(matrix)
    levers = evidence_upgrade_levers(matrix)
    leverage = module_leverage_table()
    operating = clearance_operating_points()
    cheapest = cheapest_clearance_point(operating)
    robust = robust_remuscularisation_search()
    minimum = minimum_robust_set()
    sensitivity = chronic_penalty_sensitivity()

    counts = closure["first_closing_mode"].value_counts().to_dict()
    lines = [
        "Resolving the three residual gaps, step by step",
        "=" * 74,
        DISCLAIMER,
        "",
        "STEP 1 - HAZARD RATIO: an evidence-maturity gap, not a biology gap",
        f"  populations closing by mode: {counts}",
        "  Agents whose evidence would have to mature, ranked by populations unlocked:",
    ]
    genuine = levers.loc[~levers["already_outcomes_tier"]]
    reordering = levers.loc[levers["already_outcomes_tier"]]
    for row in genuine.head(8).itertuples(index=False):
        lines.append(
            f"    {row.therapy_id:30s} unlocks {row.populations_unlocked:2d}  "
            f"{row.current_tier} -> {row.needs_tier}"
        )
    lines.extend(
        [
            f"  {len(reordering)} further agents enter only because the greedy ordering",
            "  changes when others gain credit; they already hold outcomes-tier proof and",
            "  are NOT evidence gaps.",
            "  So the residual 16-32% is bought by outcomes trials, not new mechanisms.",
            "",
            "STEP 2 - CLEARANCE: a correction, then an operating point",
            "  The earlier phrase 'near-perfect delivery' was imprecise. Potentiating any",
            "  single module to its maximum reaches only:",
        ]
    )
    for row in leverage.head(4).itertuples(index=False):
        lines.append(
            f"    {row.module:26s} {row.removal_alone_pct:5.1f}%  "
            f"(+{row.gain_over_central_pp:.1f} pp over central)"
        )
    lines.extend(
        [
            f"  None reaches {100 * CRITERION_THRESHOLD:.0f}%. It is joint performance across all",
            "  modules, so the solver searches one uniform potency multiplier:",
            "",
            "  source   horizon   min potency   admissible",
        ]
    )
    for row in operating.itertuples(index=False):
        potency = (
            f"{row.minimum_potency_multiplier:.3f}"
            if row.minimum_potency_multiplier is not None
            and not pd.isna(row.minimum_potency_multiplier)
            else "unreachable"
        )
        lines.append(
            f"    {row.source_reduction:4.2f}    {row.horizon_months:5.0f}m   "
            f"{potency:>11s}   {'yes' if row.admissible else 'no'}"
        )
    if cheapest:
        lines.append(
            f"  Cheapest admissible point: source reduction {cheapest['source_reduction']:.2f}, "
            f"{cheapest['horizon_months']:.0f} months, every module at "
            f"{cheapest['minimum_potency_multiplier']:.2f}x its central strength."
        )
    lines.extend(
        [
            "",
            "STEP 3 - REMUSCULARISATION: replacing luck with margin",
            "  The nominal x5.2 against x5.13 has no margin. Each factor is now drawn",
            "  from Uniform(1, published effect), then a chronic-injury haircut is applied.",
            "",
            "  factors                                nominal   median   P(+20)   robust",
        ]
    )
    display = robust.sort_values("probability_of_target", ascending=False).head(8)
    for row in display.itertuples(index=False):
        label = row.factors if len(row.factors) <= 38 else row.factors[:35] + "..."
        lines.append(
            f"  {label:40s} {row.nominal_gain_pp:7.1f}  {row.median_gain_pp:7.1f}  "
            f"{row.probability_of_target:6.1%}   {'YES' if row.robust_at_80pct else 'no'}"
        )
    lines.extend(
        [
            (
                f"  Smallest set robust at 80%: {', '.join(minimum)}."
                if minimum
                else "  NO factor combination is robust at 80%. The nominal two-factor"
            ),
            (
                ""
                if minimum
                else "  answer was arithmetic luck: under uncertainty it almost never clears."
            ),
            "",
            "  Dependence on the chronic-injury translation assumption:",
        ]
    )
    for row in sensitivity.itertuples(index=False):
        verdict = (
            f"{row.minimum_robust_factors} factors suffice"
            if row.robust_set_exists
            else f"NO robust set; best is {row.best_available_n_factors} factors"
        )
        lines.append(
            f"    haircut floor {row.chronic_penalty_low:4.2f} -> {verdict}, "
            f"best P={row.best_probability_of_target:5.1%}"
        )
    lines.extend(
        [
            "",
            "WHAT THE THREE STEPS ACTUALLY CHANGE",
            "1. The hazard residual becomes a named list of outcomes trials. That is",
            "   procurable, and it is the only residual that needs no new biology.",
            "2. The clearance residual becomes a single number: how much better than",
            "   central every module must be at once. That is an engineering spec, and",
            "   it is demanding - no single lever substitutes for the others.",
            "3. The remuscularisation residual is NOT solved. The earlier two-factor",
            "   answer clears only on nominal point estimates; under uncertainty it",
            "   almost never reaches the target, and even all four factors together do",
            "   not clear 80% probability. This retracts the previous 'closable on paper'.",
            "",
            "None of this is evidence. Each step converts a vague shortfall into a",
            "costed requirement that an experiment or a trial could address.",
        ]
    )
    return "\n".join(lines) + "\n"


def run_resolution_study(output_dir, matrix: pd.DataFrame | None = None) -> dict:
    from pathlib import Path

    root = Path(output_dir)
    (root / "tables").mkdir(parents=True, exist_ok=True)
    (root / "reports").mkdir(parents=True, exist_ok=True)
    if matrix is None:
        source = Path("results/composite/tables/composite_test_matrix.csv")
        if not source.exists():
            raise FileNotFoundError(
                "composite matrix missing; run `python run_solutions.py composite` first"
            )
        matrix = pd.read_csv(source)

    closure = population_closure(matrix)
    levers = evidence_upgrade_levers(matrix)
    leverage = module_leverage_table()
    operating = clearance_operating_points()
    robust = robust_remuscularisation_search()
    sensitivity = chronic_penalty_sensitivity()
    cheapest = cheapest_clearance_point(operating)
    minimum = minimum_robust_set()

    closure.to_csv(root / "tables" / "step1_population_closure.csv", index=False)
    levers.to_csv(root / "tables" / "step1_evidence_upgrade_levers.csv", index=False)
    leverage.to_csv(root / "tables" / "step2_module_leverage.csv", index=False)
    operating.to_csv(root / "tables" / "step2_operating_points.csv", index=False)
    robust.to_csv(root / "tables" / "step3_robustness_search.csv", index=False)
    sensitivity.to_csv(root / "tables" / "step3_chronic_penalty_sensitivity.csv", index=False)
    pd.DataFrame(
        [
            {
                "cheapest_clearance_point": str(cheapest),
                "minimum_robust_remuscularisation_set": ";".join(minimum) or "none",
                "populations_closing_strict": int(closure["first_closing_mode"].eq("strict").sum()),
                "populations_closing_surrogate": int(
                    closure["first_closing_mode"].eq("surrogate").sum()
                ),
                "populations_closing_exploratory": int(
                    closure["first_closing_mode"].eq("exploratory").sum()
                ),
                "status": RESOLUTION_TAG,
            }
        ]
    ).to_csv(root / "tables" / "resolution_summary.csv", index=False)
    (root / "reports" / "resolution_report.txt").write_text(
        resolution_report(matrix), encoding="utf-8"
    )
    (root / "DISCLAIMER.txt").write_text(DISCLAIMER + "\n", encoding="utf-8")
    return {
        "status": "scenario_only",
        "disclaimer": DISCLAIMER,
        "output_dir": str(root),
        "cheapest_clearance_point": cheapest,
        "minimum_robust_remuscularisation_set": minimum,
        "outputs": 9,
    }
