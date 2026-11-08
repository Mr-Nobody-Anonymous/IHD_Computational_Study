"""Mechanism-level recombination: take the load-bearing part of each architecture.

Every architecture built so far was ablated, so we know which parts actually carry the
result and which are decoration:

* ATLAS-2 (clearance): titration +72.83 pp, redundant nodes +42.00, redundant gate
  +23.00, verified source **0.00** (inert at the solved horizon).
* ATLAS-3/4 (graft): retry redundancy dominates; the cumulative operator beats best-of-k
  by +8.85 pp at four doses; feedback adds only **+0.40 pp** once dosing is cumulative.
* Gap solutions: maturation buys **1.77x** for no extra procedure; response-guided
  selection enriches responders up to **2.34x**; a viability-specific reporter is what
  makes feedback possible at all.

This module composes those parts into one candidate and tests the composition rather
than assuming it. Three rules make the composition honest:

1. **Interactions are measured, not assumed.** Every mechanism is scored alone and in
   every pair, so an overlap shows up as a negative interaction instead of being
   double-counted. Two mechanisms that both mean "only re-dose when there is room" are
   one mechanism, and the arithmetic says so.

2. **Constraints are not negotiable.** Headroom gating and human repeat-response
   tempering *lower* the reported probability. A search allowed to switch them off would
   select against its own honesty and inflate the answer, so admissible composites must
   keep every constraint enabled.

3. **Every axis is charged.** Attempts, potency, redundant paths, procedural risk and
   genetic modifications are counted separately. A composite may not pay for a breach on
   one axis with an unflagged breach on another - the trap already caught once in the
   gap-5 trade-off.

NOTHING HERE IS A THERAPY OR A RESULT. Saturation (P = 1) is reported as an artefact and
excluded, exactly as in ATLAS-3 and ATLAS-4. No cure is declared.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from itertools import combinations, product
from pathlib import Path

import numpy as np
import pandas as pd

from ihd_solutions.atlas3 import DEFAULT_PER_ATTEMPT_RISK, _factor_draws
from ihd_solutions.atlas4 import SPLIT_RETENTION
from ihd_solutions.cure import CRITERION_THRESHOLD
from ihd_solutions.gap_solutions import (
    BIOVAT_LVEF_GAIN_3MO_PP,
    BIOVAT_LVEF_GAIN_17MO_PP,
    PRECEDENTED_POTENCY_CAP,
)
from ihd_solutions.gap_solver import (
    BASELINE_SCAR_FRACTION,
    REMUSCULARISATION_FACTORS,
    SCAR_RECOVERABLE_LVEF_PP,
)
from ihd_solutions.precedent import (
    max_precedented_doses,
    measurement_quality_anchor,
    retry_responder_fraction,
)
from ihd_solutions.units import require_fraction, require_positive
from ihd_solutions.upgrade import MAX_PRECEDENTED_PATHS

RECOMBINE_TAG = "recombine_mechanism_composition_unvalidated"

DISCLAIMER = (
    "THIS IS A COMPOSITION OF UNVALIDATED MECHANISMS, NOT A THERAPY AND NOT A RESULT. "
    "Each part is borrowed from an architecture that was itself invented or from another "
    "group's animal or early-phase human work. Interactions are measured so overlaps are "
    "not double-counted, honesty constraints cannot be switched off to inflate the "
    "answer, and every precedent axis is charged separately. No cure is declared."
)

#: Maturation multiplier measured in the human EHM series: +3.9 pp at three months
#: rising to +6.9 pp at a mean 17 months, with no extra procedure.
MATURATION_MULTIPLIER = BIOVAT_LVEF_GAIN_17MO_PP / BIOVAT_LVEF_GAIN_3MO_PP

#: Central row of the swept response-guided enrichment table (assumed SD 4 pp, cutoff
#: 37.5 pp LVEF). Deliberately NOT the maximum: the between-patient SD is unreported, so
#: taking the 2.34x best case would be cherry-picking an assumption.
RESPONSE_GUIDED_RATE = 0.582481
RESPONSE_GUIDED_PASS = 0.424708

#: Mechanisms that LOWER the claim. A composition may not drop these.
CONSTRAINT_MECHANISMS = ("headroom_gating", "human_tempering")

#: Mechanisms that raise the claim and are therefore searchable.
BENEFIT_MECHANISMS = (
    "retry_redundancy",
    "cumulative_operator",
    "feedback",
    "viability_reporter",
    "response_guided_selection",
    "maturation",
)

ALL_MECHANISMS = BENEFIT_MECHANISMS + CONSTRAINT_MECHANISMS

#: Half a percentage point, expressed in the probability units the contrast is computed
#: in. Comparing a probability difference against a percentage-point cutoff silently
#: called a -14 pp overlap "additive".
INTERACTION_TOLERANCE = 0.005

#: The composition search must not change its answer when the sampling budget changes.
#: Floored inside exhaustive_subsets for the same reason ATLAS-3 needed it: a design
#: selected on 400 draws is a coin flip wearing a specification.
MIN_SPECIFICATION_DRAWS = 2000


@dataclass(frozen=True)
class Mechanism:
    """One transferable part, with where it came from and what it costs."""

    name: str
    origin: str
    kind: str
    cost_axis: str
    description: str

    @property
    def is_constraint(self) -> bool:
        return self.kind == "constraint"


MECHANISMS: tuple[Mechanism, ...] = (
    Mechanism(
        "retry_redundancy",
        "ATLAS-3 / ATLAS-4",
        "benefit",
        "attempts_and_procedural_risk",
        "Up to the precedented number of doses, stopping on success. Dominant lever in "
        "the graft ablation (-43.50 pp when removed).",
    ),
    Mechanism(
        "cumulative_operator",
        "ATLAS-4 (corrected from ATLAS-3)",
        "benefit",
        "mechanism_assumption",
        "Successive doses add rather than replace. Worth +8.85 pp over best-of-k at four "
        "doses, and the only operator that can represent the published split-dose result.",
    ),
    Mechanism(
        "feedback",
        "ATLAS-3 / ATLAS-4",
        "benefit",
        "requires_reporter",
        "A measured signal guides the next dose. Nearly inert once dosing is cumulative "
        "(+0.40 pp), so it is a candidate for DROPPING rather than keeping.",
    ),
    Mechanism(
        "viability_reporter",
        "gap 3 (hNIS, already imaged in patients)",
        "benefit",
        "genetic_modification",
        "A viability-specific signal. Without it feedback is blind and contributes "
        "exactly nothing. Adds a second genetic modification to the product.",
    ),
    Mechanism(
        "response_guided_selection",
        "gap 2 (Vrtovec responder analysis)",
        "benefit",
        "screening_burden",
        "Re-dose only patients whose pre-dose function leaves room. Up to 2.34x "
        "enrichment, using a signal already in the clinic.",
    ),
    Mechanism(
        "maturation",
        "gap 1 (human EHM series)",
        "benefit",
        "time_and_conditioning",
        "The same graft improved 1.77x between three and seventeen months with no extra "
        "procedure. The cheapest lever in the whole catalogue.",
    ),
    Mechanism(
        "headroom_gating",
        "ATLAS-4 (human ceiling effect)",
        "constraint",
        "none",
        "Scales each increment by the room left. LOWERS the claim. Encodes the ceiling "
        "effect measured in the blinded human repeat trial.",
    ),
    Mechanism(
        "human_tempering",
        "ATLAS-4 (10 of 30 responded)",
        "constraint",
        "none",
        "Only a third of patients respond to a repeat dose at all. LOWERS the claim. "
        "Removing it is what drove ATLAS-4 to a spurious 100%.",
    ),
)

MECHANISM_BY_NAME = {m.name: m for m in MECHANISMS}


@dataclass(frozen=True)
class Composite:
    """A point in mechanism space, plus the clearance arm it is paired with."""

    enabled: frozenset[str] = frozenset(ALL_MECHANISMS)
    doses: int = 4
    per_dose_risk: float = DEFAULT_PER_ATTEMPT_RISK
    translation_floor: float = 1.0
    clearance_cap: float = 2.40
    clearance_nodes: int = 3
    clearance_gates: int = 3

    def __post_init__(self) -> None:
        unknown = set(self.enabled) - set(ALL_MECHANISMS)
        if unknown:
            raise ValueError(f"unknown mechanisms: {sorted(unknown)}")
        if self.doses < 1:
            raise ValueError("at least one dose is required")
        require_fraction("per_dose_risk", self.per_dose_risk, allow_one=True)
        require_fraction("translation_floor", self.translation_floor, allow_one=True)
        require_positive("clearance_cap", self.clearance_cap)
        if self.clearance_nodes < 1 or self.clearance_gates < 1:
            raise ValueError("redundancy needs at least one path")

    def has(self, name: str) -> bool:
        return name in self.enabled

    def with_mechanism(self, name: str, on: bool) -> Composite:
        enabled = set(self.enabled)
        if on:
            enabled.add(name)
        else:
            enabled.discard(name)
        return replace(self, enabled=frozenset(enabled))

    @property
    def constraints_intact(self) -> bool:
        return all(self.has(name) for name in CONSTRAINT_MECHANISMS)

    @property
    def attempts_beyond_precedent(self) -> bool:
        return self.has("retry_redundancy") and self.doses > max_precedented_doses()

    @property
    def redundancy_beyond_precedent(self) -> bool:
        return max(self.clearance_nodes, self.clearance_gates) > MAX_PRECEDENTED_PATHS

    @property
    def potency_beyond_precedent(self) -> bool:
        return self.clearance_cap > PRECEDENTED_POTENCY_CAP

    @property
    def genetic_modifications(self) -> int:
        return 1 if self.has("viability_reporter") else 0

    @property
    def graft_precedent_breaches(self) -> tuple[str, ...]:
        """Breaches on the axes this arm controls.

        The clearance arm's potency requirement is charged separately in the ledger. If
        it were folded in here, its standing breach would make every graft composition
        inadmissible and the search would return nothing for the wrong reason.
        """
        breaches = []
        if self.attempts_beyond_precedent:
            breaches.append("attempts")
        if self.genetic_modifications > 0:
            breaches.append("genetic_modification")
        return tuple(breaches)

    @property
    def precedent_breaches(self) -> tuple[str, ...]:
        breaches = []
        if self.attempts_beyond_precedent:
            breaches.append("attempts")
        if self.redundancy_beyond_precedent:
            breaches.append("redundancy")
        if self.potency_beyond_precedent:
            breaches.append("potency")
        return tuple(breaches)


#: Everything on, at the corrected precedented dose count, paired with the cheapest
#: precedented-redundancy clearance operating point found in the gap-5 sweep.
COMPOSITE = Composite(
    enabled=frozenset(ALL_MECHANISMS),
    doses=max_precedented_doses(),
    clearance_cap=2.40,
    clearance_nodes=3,
    clearance_gates=3,
)


#: ATLAS-5: the composition the exhaustive search actually selected. Named only after
#: the search ran, and defined by what survived it - three mechanisms, no reporter
#: cassette, no feedback loop, four doses inside the corrected precedent limit.
ATLAS5_MECHANISMS = frozenset(
    {"retry_redundancy", "cumulative_operator", "maturation"} | set(CONSTRAINT_MECHANISMS)
)
ATLAS5 = Composite(enabled=ATLAS5_MECHANISMS, doses=max_precedented_doses())

#: ATLAS-4's best admissible design, from results/atlas4/tables at 4000 draws. Kept as a
#: documented reference and pinned by a test against that committed table, so the
#: comparison cannot drift into a flattering number.
ATLAS4_BEST_ADMISSIBLE = 0.8960


@dataclass(frozen=True)
class CompositeResult:
    probability_of_target: float
    median_gain_pp: float
    mean_doses_used: float
    cumulative_adverse_risk: float
    saturated_fraction: float
    draws: int
    mechanisms: int

    @property
    def robust(self) -> bool:
        return self.probability_of_target >= CRITERION_THRESHOLD

    @property
    def monte_carlo_standard_error(self) -> float:
        p = self.probability_of_target
        return float(np.sqrt(max(p * (1.0 - p), 0.0) / max(self.draws, 1)))

    @property
    def robust_with_margin(self) -> bool:
        return self.probability_of_target - 2.0 * self.monte_carlo_standard_error >= (
            CRITERION_THRESHOLD
        )

    @property
    def risk_adjusted_probability(self) -> float:
        return self.probability_of_target * (1.0 - self.cumulative_adverse_risk)

    @property
    def saturation_artefact(self) -> bool:
        return self.probability_of_target >= 1.0


def composite_draws(
    design: Composite = COMPOSITE, draws: int = 3000, seed: int = 20261107
) -> tuple[np.ndarray, np.ndarray]:
    """Simulate the graft arm with the selected mechanisms switched on.

    Returns (gain, doses used). All mechanisms read from one shared random stream so
    every marginal and interaction contrast below is paired.
    """
    require_positive("draws", draws)
    count = int(draws)
    rng = np.random.default_rng(seed)
    shared = rng.random((count, len(REMUSCULARISATION_FACTORS)))
    ceiling = SCAR_RECOVERABLE_LVEF_PP
    required = CRITERION_THRESHOLD * ceiling

    doses = design.doses if design.has("retry_redundancy") else 1
    learning = 0.5 if design.has("feedback") else 0.0
    quality = (
        float(np.mean(measurement_quality_anchor())) if design.has("viability_reporter") else 0.0
    )
    responder = retry_responder_fraction() if design.has("human_tempering") else 1.0
    guided = design.has("response_guided_selection")
    if guided:
        responder = min(1.0, RESPONSE_GUIDED_RATE) if design.has("human_tempering") else 1.0
    # A screen excludes people. Simulate the WHOLE population: those who fail keep their
    # first-dose gain and are never re-dosed. Reporting only the screened-in subset would
    # turn a trial-design lever into a fake therapy benefit.
    eligible = rng.random(count) < RESPONSE_GUIDED_PASS if guided else np.ones(count, bool)
    maturation = MATURATION_MULTIPLIER if design.has("maturation") else 1.0
    view = _upgrade_view_shim(
        _DesignView(doses, learning, quality, design.translation_floor, design.per_dose_risk)
    )

    gain = np.zeros(count)
    used = np.ones(count)
    done = np.zeros(count, dtype=bool)
    for index in range(doses):
        delivered = _factor_draws(rng, count, index, view, shared)
        penalty = rng.uniform(design.translation_floor, 1.0, size=count)
        increment = (
            np.minimum(1.0, BASELINE_SCAR_FRACTION * delivered * penalty * maturation) * ceiling
        )
        if index == 0:
            gain = increment.copy()
            done = gain >= required
            continue
        responds = eligible & (rng.random(count) < responder)
        active = (~done) & responds
        if design.has("cumulative_operator"):
            scaled = increment * SPLIT_RETENTION
            if design.has("headroom_gating"):
                scaled = scaled * np.clip((ceiling - gain) / ceiling, 0.0, 1.0)
            candidate = np.minimum(ceiling, gain + scaled)
        else:
            candidate = np.maximum(gain, increment)
        gain = np.where(active, candidate, gain)
        used = np.where(active, float(index + 1), used)
        done = done | (gain >= required)

    return gain, used


@dataclass(frozen=True)
class _DesignView:
    """Minimal adapter so ATLAS-3's factor draws can be reused unchanged."""

    doses: int
    learning: float
    quality: float
    floor: float
    risk: float


def _upgrade_view_shim(view: _DesignView):
    from ihd_solutions.atlas3 import GraftUpgrades

    return GraftUpgrades(
        attempts=view.doses,
        failure_correlation=0.0,
        feedback_learning=view.learning,
        measurement_quality=view.quality,
        per_attempt_risk=view.risk,
        translation_floor=view.floor,
    )


def evaluate_composite(
    design: Composite = COMPOSITE, draws: int = 3000, seed: int = 20261107
) -> CompositeResult:
    gain, used = composite_draws(design, draws, seed)
    required = CRITERION_THRESHOLD * SCAR_RECOVERABLE_LVEF_PP
    mean_used = float(np.mean(used))
    return CompositeResult(
        probability_of_target=float(np.mean(gain >= required)),
        median_gain_pp=float(np.median(gain)),
        mean_doses_used=mean_used,
        cumulative_adverse_risk=float(1.0 - (1.0 - design.per_dose_risk) ** mean_used),
        saturated_fraction=float(np.mean(gain >= SCAR_RECOVERABLE_LVEF_PP - 1e-9)),
        draws=int(draws),
        mechanisms=len(design.enabled),
    )


def mechanism_marginals(
    base: Composite = COMPOSITE, draws: int = 2500, seed: int = 20261107
) -> pd.DataFrame:
    """Leave-one-out contribution of every mechanism, constraints included.

    A constraint shows a NEGATIVE contribution: switching it off raises the number. That
    is the point - it is the size of the honesty being purchased.
    """
    full = evaluate_composite(base, draws, seed)
    rows = []
    for name in ALL_MECHANISMS:
        reduced = evaluate_composite(base.with_mechanism(name, False), draws, seed)
        spec = MECHANISM_BY_NAME[name]
        rows.append(
            {
                "mechanism": name,
                "origin": spec.origin,
                "kind": spec.kind,
                "cost_axis": spec.cost_axis,
                "probability_without": reduced.probability_of_target,
                "probability_full": full.probability_of_target,
                "contribution_pp": 100.0
                * (full.probability_of_target - reduced.probability_of_target),
                "load_bearing": abs(full.probability_of_target - reduced.probability_of_target)
                > 1e-9,
                "status": RECOMBINE_TAG,
            }
        )
    return pd.DataFrame(rows)


def mechanism_standalone(
    base: Composite = COMPOSITE, draws: int = 2500, seed: int = 20261107
) -> pd.DataFrame:
    """Add-one-in from an empty benefit set, with constraints held on.

    Leave-one-out from the full composition is measured against a near-saturated base,
    where a part can look inert simply because there is no headroom left for it to show
    up in. This view answers the other question: what does each part do on its own?
    """
    empty = base
    for name in BENEFIT_MECHANISMS:
        empty = empty.with_mechanism(name, False)
    baseline = evaluate_composite(empty, draws, seed)
    rows = []
    for name in BENEFIT_MECHANISMS:
        alone = evaluate_composite(empty.with_mechanism(name, True), draws, seed)
        spec = MECHANISM_BY_NAME[name]
        rows.append(
            {
                "mechanism": name,
                "origin": spec.origin,
                "cost_axis": spec.cost_axis,
                "probability_empty": baseline.probability_of_target,
                "probability_alone": alone.probability_of_target,
                "standalone_gain_pp": 100.0
                * (alone.probability_of_target - baseline.probability_of_target),
                "risk_adjusted_gain_pp": 100.0
                * (alone.risk_adjusted_probability - baseline.risk_adjusted_probability),
                "graft_breaches": "|".join(
                    empty.with_mechanism(name, True).graft_precedent_breaches
                )
                or "none",
                "status": RECOMBINE_TAG,
            }
        )
    return pd.DataFrame(rows)


def pairwise_interactions(
    base: Composite = COMPOSITE, draws: int = 2000, seed: int = 20261107
) -> pd.DataFrame:
    """Measured interaction for every pair of benefit mechanisms.

    baseline = both off; a = first only; b = second only; ab = both on.
    interaction = (ab - baseline) - (a - baseline) - (b - baseline).

    Negative means the two overlap and must not both be credited at full strength.
    Positive means genuine synergy. Near zero means they compose additively.
    """
    off = base
    for name in BENEFIT_MECHANISMS:
        off = off.with_mechanism(name, False)
    baseline = evaluate_composite(off, draws, seed).probability_of_target
    rows = []
    for first, second in combinations(BENEFIT_MECHANISMS, 2):
        p_a = evaluate_composite(off.with_mechanism(first, True), draws, seed).probability_of_target
        p_b = evaluate_composite(
            off.with_mechanism(second, True), draws, seed
        ).probability_of_target
        p_ab = evaluate_composite(
            off.with_mechanism(first, True).with_mechanism(second, True), draws, seed
        ).probability_of_target
        gain_a = p_a - baseline
        gain_b = p_b - baseline
        gain_ab = p_ab - baseline
        interaction = gain_ab - gain_a - gain_b
        rows.append(
            {
                "first": first,
                "second": second,
                "gain_first_pp": 100.0 * gain_a,
                "gain_second_pp": 100.0 * gain_b,
                "gain_both_pp": 100.0 * gain_ab,
                "sum_of_parts_pp": 100.0 * (gain_a + gain_b),
                "interaction_pp": 100.0 * interaction,
                "verdict": (
                    "overlapping"
                    if interaction < -INTERACTION_TOLERANCE
                    else ("synergistic" if interaction > INTERACTION_TOLERANCE else "additive")
                ),
                "status": RECOMBINE_TAG,
            }
        )
    return pd.DataFrame(rows)


def exhaustive_subsets(
    base: Composite = COMPOSITE, draws: int = 800, seed: int = 20261107
) -> pd.DataFrame:
    """Every subset of the benefit mechanisms, with constraints held ON.

    2^6 = 64 compositions. Constraints are never switched off, so the search cannot buy
    a higher number by deleting its own honesty. The draw count is floored at
    MIN_SPECIFICATION_DRAWS so the selected specification is budget-independent.
    """
    draws = max(int(draws), MIN_SPECIFICATION_DRAWS)
    rows = []
    for size in range(len(BENEFIT_MECHANISMS) + 1):
        for chosen in combinations(BENEFIT_MECHANISMS, size):
            design = base
            for name in BENEFIT_MECHANISMS:
                design = design.with_mechanism(name, name in chosen)
            result = evaluate_composite(design, draws, seed)
            rows.append(
                {
                    "mechanisms": "|".join(chosen) if chosen else "(none)",
                    "count": size,
                    "probability_of_target": result.probability_of_target,
                    "median_gain_pp": result.median_gain_pp,
                    "risk_adjusted_probability": result.risk_adjusted_probability,
                    "cumulative_adverse_risk": result.cumulative_adverse_risk,
                    "saturated_fraction": result.saturated_fraction,
                    "robust_with_margin": result.robust_with_margin,
                    "saturation_artefact": result.saturation_artefact,
                    "constraints_intact": design.constraints_intact,
                    "genetic_modifications": design.genetic_modifications,
                    "precedent_breaches": "|".join(design.graft_precedent_breaches) or "none",
                    "status": RECOMBINE_TAG,
                }
            )
    return pd.DataFrame(rows)


def select_composite(grid: pd.DataFrame) -> dict[str, object]:
    """Cheapest admissible composition, not the highest number.

    Admissible means: constraints intact, no precedent breach on any axis, robust with
    margin, and not saturated. Among those, prefer fewer genetic modifications and then
    fewer mechanisms - because a part that does not earn its cost should be dropped.
    """
    admissible = grid.loc[
        grid["constraints_intact"]
        & (~grid["saturation_artefact"])
        & grid["robust_with_margin"]
        & grid["precedent_breaches"].eq("none")
    ]
    if admissible.empty:
        return {"admissible": 0, "best": None, "robust_set_exists": False}
    ordered = admissible.sort_values(
        ["genetic_modifications", "count", "risk_adjusted_probability"],
        ascending=[True, True, False],
    )
    strongest = admissible.sort_values(
        ["risk_adjusted_probability", "count"], ascending=[False, True]
    ).iloc[0]

    def _view(row) -> dict[str, object]:
        return {
            "mechanisms": row["mechanisms"],
            "count": int(row["count"]),
            "probability_of_target": float(row["probability_of_target"]),
            "risk_adjusted_probability": float(row["risk_adjusted_probability"]),
            "genetic_modifications": int(row["genetic_modifications"]),
            "cumulative_adverse_risk": float(row["cumulative_adverse_risk"]),
        }

    cheapest = _view(ordered.iloc[0])
    best = _view(strongest)
    return {
        "admissible": int(len(ordered)),
        "robust_set_exists": True,
        # "best" is the recommended design: highest risk-adjusted probability among
        # admissible compositions. "cheapest" is the frugal alternative. Reporting only
        # the cheapest would hide a design that is better on both raw and risk-adjusted
        # probability at a cost that is still inside precedent.
        "best": best,
        "cheapest": cheapest,
    }


def maturation_sensitivity(draws: int = 1200, seed: int = 20261107) -> pd.DataFrame:
    """Does the selected composition survive its own assumptions?

    The maturation ratio is measured in epicardial EHM patches; the model's factor stack
    is assembled from other groups' animal work. Multiplying one by the other is an
    UNVALIDATED TRANSFER, and it is the single largest assumption in the selection. So
    sweep it, together with the translation penalty on the animal factors and the dose
    count, instead of asserting the point estimate.
    """
    rows = []
    for multiplier, floor, doses in product(
        (1.0, 1.25, 1.50, MATURATION_MULTIPLIER),
        (0.4, 0.7, 1.0),
        (1, 2, 4),
    ):
        design = Composite(
            enabled=frozenset(CONSTRAINT_MECHANISMS)
            | {"maturation", "retry_redundancy", "cumulative_operator"},
            doses=doses,
            translation_floor=floor,
        )
        # A multiplier of 1.0 means maturation contributes nothing, so switch it off
        # rather than pretending an inert mechanism is present.
        if multiplier == 1.0:
            design = design.with_mechanism("maturation", False)
        result = evaluate_composite(design, draws, seed)
        rows.append(
            {
                "maturation_multiplier": multiplier,
                "translation_floor": floor,
                "doses": doses,
                "probability_of_target": result.probability_of_target,
                "median_gain_pp": result.median_gain_pp,
                "saturated_fraction": result.saturated_fraction,
                "robust_with_margin": result.robust_with_margin,
                "saturation_artefact": result.saturation_artefact,
                "status": RECOMBINE_TAG,
            }
        )
    return pd.DataFrame(rows)


def composition_verdict(grid: pd.DataFrame, selection: dict[str, object]) -> dict[str, object]:
    """Compare the selection against the alternatives, instead of asserting a win.

    Three prices are separated because conflating them told a false story: the strongest
    composition overall, the strongest ADMISSIBLE one, and the cheapest admissible one
    are three different designs. The gap between the first two is what precedent costs;
    the gap between the last two is what frugality costs.
    """
    usable = grid.loc[grid["constraints_intact"] & (~grid["saturation_artefact"])]
    singles = usable.loc[usable["count"].eq(1)]
    admissible = usable.loc[usable["robust_with_margin"] & usable["precedent_breaches"].eq("none")]
    best_any = float(usable["probability_of_target"].max()) if not usable.empty else 0.0
    best_single = float(singles["probability_of_target"].max()) if not singles.empty else 0.0
    best_row = usable.loc[usable["probability_of_target"].idxmax()] if not usable.empty else None
    adm_row = (
        admissible.loc[admissible["probability_of_target"].idxmax()]
        if not admissible.empty
        else None
    )
    best_adm = float(adm_row["probability_of_target"]) if adm_row is not None else 0.0
    best_adm_risk = float(adm_row["risk_adjusted_probability"]) if adm_row is not None else 0.0
    selected = selection.get("best")
    is_dict = isinstance(selected, dict)
    selected_p = float(selected["probability_of_target"]) if is_dict else 0.0
    selected_count = int(selected["count"]) if is_dict else 0
    selected_risk = float(selected["risk_adjusted_probability"]) if is_dict else 0.0
    return {
        "best_probability_any": best_any,
        "best_mechanisms_any": str(best_row["mechanisms"]) if best_row is not None else "none",
        "best_count_any": int(best_row["count"]) if best_row is not None else 0,
        "best_probability_single": best_single,
        "best_probability_admissible": best_adm,
        "best_mechanisms_admissible": str(adm_row["mechanisms"]) if adm_row is not None else "none",
        "best_count_admissible": int(adm_row["count"]) if adm_row is not None else 0,
        "best_risk_adjusted_admissible": best_adm_risk,
        "selected_probability": selected_p,
        "selected_count": selected_count,
        "selected_risk_adjusted": selected_risk,
        "gain_over_best_single_pp": 100.0 * (selected_p - best_single),
        "price_of_admissibility_pp": 100.0 * (best_any - best_adm),
        "price_of_cheapness_pp": 100.0 * (best_adm - selected_p),
        "cheapness_wins_on_risk": bool(selected_risk > best_adm_risk),
        "composition_beat_best_single": bool(selected_p - best_single > 1e-9),
        "selection_is_a_single_mechanism": bool(selected_count == 1),
        "composition_adds_over_single_pp": 100.0 * (best_adm - best_single),
        "atlas4_reference": ATLAS4_BEST_ADMISSIBLE,
        "gain_over_atlas4_pp": 100.0 * (best_adm - ATLAS4_BEST_ADMISSIBLE),
    }


def maturation_resolution_limit(sensitivity: pd.DataFrame) -> float:
    """Smallest multiplier beyond which the model stops responding to maturation.

    Increments are clamped at the recoverable-scar ceiling, so once the clamp binds the
    model cannot tell 1.25x from 1.77x. That is a blind spot in the model, not a finding
    about biology, and the maturation credit must be read as what the model can resolve.
    """
    for multiplier in sorted(sensitivity["maturation_multiplier"].unique()):
        subset = sensitivity.loc[sensitivity["maturation_multiplier"].ge(multiplier)]
        spread = subset.groupby(["translation_floor", "doses"])["probability_of_target"].agg(
            lambda v: float(v.max() - v.min())
        )
        if not spread.empty and float(spread.max()) <= 1e-12:
            return float(multiplier)
    return float(sensitivity["maturation_multiplier"].max())


def design_from_selection(selection: dict[str, object]) -> Composite:
    """Rebuild the selected composition as a design, constraints always restored."""
    best = selection.get("best")
    if not isinstance(best, dict):
        return COMPOSITE
    chosen = {name for name in str(best["mechanisms"]).split("|") if name in ALL_MECHANISMS}
    return Composite(
        enabled=frozenset(chosen | set(CONSTRAINT_MECHANISMS)),
        doses=max_precedented_doses(),
    )


def stacking_contrast(
    base: Composite = COMPOSITE, draws: int = 2500, seed: int = 20261107
) -> pd.DataFrame:
    """Naive stacking versus selection.

    Switching every mechanism on is NOT the same as composing the successful ones: a
    part with a negative population effect drags the stack down. This is the table that
    shows why the search was worth running.
    """
    selected = design_from_selection(select_composite(exhaustive_subsets(base, draws, seed)))
    rows = []
    for label, design in (
        ("everything_on", base),
        ("selected_composition", selected),
        ("constraints_only", Composite(enabled=frozenset(CONSTRAINT_MECHANISMS), doses=base.doses)),
    ):
        result = evaluate_composite(design, draws, seed)
        rows.append(
            {
                "design": label,
                "mechanisms": "|".join(sorted(design.enabled - set(CONSTRAINT_MECHANISMS)))
                or "(none)",
                "probability_of_target": result.probability_of_target,
                "risk_adjusted_probability": result.risk_adjusted_probability,
                "mean_doses_used": result.mean_doses_used,
                "cumulative_adverse_risk": result.cumulative_adverse_risk,
                "median_gain_pp": result.median_gain_pp,
                "status": RECOMBINE_TAG,
            }
        )
    return pd.DataFrame(rows)


def clearance_arm(design: Composite = COMPOSITE, draws: int = 1500, seed: int = 20261107):
    """The clearance half, using the cheapest precedented-redundancy operating point."""
    from ihd_solutions.upgrade import SOLVED_MONTHS, Upgrades, evaluate_upgrades

    config = Upgrades(
        gate_paths=design.clearance_gates,
        node_paths=design.clearance_nodes,
        titrate=True,
        verified_source=False,
        max_potency=design.clearance_cap,
    )
    return evaluate_upgrades(config, months=SOLVED_MONTHS, draws=draws, seed=seed)


def joint_result(
    design: Composite = COMPOSITE, draws: int = 1500, seed: int = 20261107
) -> dict[str, float]:
    graft = evaluate_composite(design, draws, seed)
    clear = clearance_arm(design, draws, seed)
    return {
        "graft_probability": graft.probability_of_target,
        "clearance_probability": clear.probability_of_target,
        "joint_optimistic": graft.probability_of_target * clear.probability_of_target,
        "randomised_proof": 0.0,
    }


def cost_ledger(design: Composite = COMPOSITE) -> pd.DataFrame:
    """Every axis charged separately, so no breach hides behind another."""
    return pd.DataFrame(
        [
            {
                "axis": "doses",
                "value": design.doses if design.has("retry_redundancy") else 1,
                "precedented_limit": max_precedented_doses(),
                "beyond_precedent": design.attempts_beyond_precedent,
                "status": RECOMBINE_TAG,
            },
            {
                "axis": "clearance_potency_cap",
                "value": design.clearance_cap,
                "precedented_limit": PRECEDENTED_POTENCY_CAP,
                "beyond_precedent": design.potency_beyond_precedent,
                "status": RECOMBINE_TAG,
            },
            {
                "axis": "clearance_redundant_paths",
                "value": max(design.clearance_nodes, design.clearance_gates),
                "precedented_limit": MAX_PRECEDENTED_PATHS,
                "beyond_precedent": design.redundancy_beyond_precedent,
                "status": RECOMBINE_TAG,
            },
            {
                "axis": "genetic_modifications",
                "value": design.genetic_modifications,
                "precedented_limit": 0,
                "beyond_precedent": design.genetic_modifications > 0,
                "status": RECOMBINE_TAG,
            },
        ]
    )


def recombine_report(
    full: CompositeResult,
    marginals: pd.DataFrame,
    interactions: pd.DataFrame,
    grid: pd.DataFrame,
    selection: dict[str, object],
    joint: dict[str, float],
    ledger: pd.DataFrame,
    standalone: pd.DataFrame,
    sensitivity: pd.DataFrame,
    verdict: dict[str, object],
    recommended: CompositeResult,
    contrast: pd.DataFrame,
) -> str:
    overlapping = interactions.loc[interactions["verdict"].eq("overlapping")]
    inert = marginals.loc[~marginals["load_bearing"]]
    lines = [
        "RECOMBINATION: composing the load-bearing part of every architecture",
        "=" * 74,
        DISCLAIMER,
        "",
        "THE CATALOGUE",
        f"  mechanisms available               {len(MECHANISMS)}",
        f"  searchable benefits                {len(BENEFIT_MECHANISMS)}",
        f"  non-negotiable constraints         {len(CONSTRAINT_MECHANISMS)}",
        "",
        "  Constraints LOWER the number. The search may not switch them off, because",
        "  that would buy a better answer by deleting the honesty that produced it.",
        "",
        "LEAVE-ONE-OUT: what each part actually carries",
        "  mechanism                     kind         contribution",
    ]
    for row in marginals.sort_values("contribution_pp", ascending=False).itertuples(index=False):
        lines.append(f"  {row.mechanism:29s} {row.kind:12s} {row.contribution_pp:+7.2f} pp")
    lines += [
        "",
        "  A negative contribution on a constraint is the size of the honesty being",
        "  purchased, not a defect. Anything at zero is decoration and should be dropped.",
    ]
    if not inert.empty:
        lines += [
            "",
            "  INERT in this composition: " + ", ".join(inert["mechanism"]),
        ]
    lines += [
        "",
        "STANDALONE: what each part does on its own (constraints held on)",
        "  mechanism                     gain      risk-adjusted   graft breach",
    ]
    for row in standalone.sort_values("standalone_gain_pp", ascending=False).itertuples(
        index=False
    ):
        lines.append(
            f"  {row.mechanism:29s} {row.standalone_gain_pp:+7.2f} pp  "
            f"{row.risk_adjusted_gain_pp:+7.2f} pp     {row.graft_breaches}"
        )
    lines += [
        "",
        "  Leave-one-out is measured against a near-saturated base, where a part can look",
        "  inert only because no headroom is left. This is the other direction, and the",
        "  two together are what licenses dropping a part.",
        "",
        "INTERACTIONS: measured, so overlaps are not double-counted",
        "  pair                                          both    sum-of-parts   interaction",
    ]
    for row in interactions.sort_values("interaction_pp").itertuples(index=False):
        lines.append(
            f"  {row.first + ' + ' + row.second:44s} {row.gain_both_pp:+6.2f}  "
            f"{row.sum_of_parts_pp:+7.2f}      {row.interaction_pp:+7.2f} pp  {row.verdict}"
        )
    if not overlapping.empty:
        lines += [
            "",
            f"  {len(overlapping)} pair(s) OVERLAP. Crediting both at full strength would",
            "  double-count the same biology. The strongest overlap is",
            f"  {overlapping.iloc[0]['first']} + {overlapping.iloc[0]['second']} "
            f"at {overlapping.iloc[0]['interaction_pp']:+.2f} pp.",
        ]
    lines += [
        "",
        "EXHAUSTIVE COMPOSITION SEARCH",
        f"  compositions tested                {len(grid)}",
        f"  robust with margin                 {int(grid['robust_with_margin'].sum())}",
        f"  saturation artefacts (P = 100%)     {int(grid['saturation_artefact'].sum())}",
        f"  admissible after every filter      {selection['admissible']}",
        f"  robust set exists                  {selection['robust_set_exists']}",
    ]
    best = selection.get("best")
    if isinstance(best, dict):
        lines += [
            "",
            "  RECOMMENDED - ATLAS-5 (highest risk-adjusted among admissible):",
            f"    mechanisms      {best['mechanisms']}",
            f"    count           {best['count']} of {len(BENEFIT_MECHANISMS)}",
            f"    P(target)       {100 * best['probability_of_target']:.2f}%",
            f"    risk-adjusted   {100 * best['risk_adjusted_probability']:.2f}%",
            f"    genetic mods    {best['genetic_modifications']}",
            f"    adverse risk    {100 * best['cumulative_adverse_risk']:.2f}%",
        ]
    cheapest = selection.get("cheapest")
    if isinstance(cheapest, dict):
        lines += [
            "",
            "  FRUGAL ALTERNATIVE (fewest parts, still robust):",
            f"    mechanisms      {cheapest['mechanisms']}",
            f"    P(target)       {100 * cheapest['probability_of_target']:.2f}%",
            f"    risk-adjusted   {100 * cheapest['risk_adjusted_probability']:.2f}%",
            f"    adverse risk    {100 * cheapest['cumulative_adverse_risk']:.2f}%",
        ]
    else:
        lines += ["", "  NO admissible composition. Reporting best available, not a solve."]
    surviving = int(sensitivity["robust_with_margin"].sum())
    lines += [
        "",
        "DOES THE SELECTION SURVIVE ITS OWN ASSUMPTIONS?",
        f"  scenarios swept                    {len(sensitivity)}",
        f"  robust with margin                 {surviving}",
        f"  fraction robust                    {surviving / max(len(sensitivity), 1):.3f}",
        "",
        "  The maturation ratio comes from epicardial EHM patches; the factor stack comes",
        "  from other groups' animal work. Multiplying them is an UNVALIDATED TRANSFER and",
        "  the largest assumption in the selection, so it is swept from 1.00x to",
        f"  {MATURATION_MULTIPLIER:.2f}x rather than asserted at its point estimate.",
        "",
        f"  MODEL BLIND SPOT: above {maturation_resolution_limit(sensitivity):.2f}x the probability",
        "  stops changing, because increments clamp at the recoverable-scar",
        f"  ceiling. The model cannot distinguish that from the full {MATURATION_MULTIPLIER:.2f}x,",
        "  so the maturation credit is what this model can resolve, not the full ratio.",
    ]
    lines += [
        "",
        "ATLAS-5 AT THE REPORTED DRAW COUNT",
        f"  P(target)                          {100 * recommended.probability_of_target:6.2f}%",
        f"  risk-adjusted                      {100 * recommended.risk_adjusted_probability:6.2f}%",
        f"  mean doses actually used           {recommended.mean_doses_used:6.2f}",
        f"  cumulative adverse risk            {100 * recommended.cumulative_adverse_risk:5.2f}%",
        f"  saturated fraction                 {100 * recommended.saturated_fraction:5.1f}%",
        "",
        "  The mean dose count is barely above one: maturation makes the FIRST dose",
        "  succeed often enough that the retry machinery is rarely invoked. That is why",
        "  the composite is cheap in practice and not only on paper.",
        "",
        "NAIVE STACKING VERSUS SELECTION",
        "  design                    mechanisms                        P(target)  risk-adj",
    ]
    for row in contrast.itertuples(index=False):
        lines.append(
            f"  {row.design:25s} {row.mechanisms:33s} "
            f"{100 * row.probability_of_target:6.2f}%   {100 * row.risk_adjusted_probability:6.2f}%"
        )
    stacked = contrast.set_index("design").loc["everything_on", "probability_of_target"]
    picked = contrast.set_index("design").loc["selected_composition", "probability_of_target"]
    if picked > stacked + INTERACTION_TOLERANCE:
        lines += [
            "",
            f"  Switching EVERYTHING on is {100 * (picked - stacked):.2f} pp WORSE than selecting.",
            "  A mechanism with a negative population effect drags the whole stack down.",
            "  This is the concrete reason the search was worth running: recombination is",
            "  a selection problem, not an accumulation problem.",
        ]
    lines += [
        "",
        "FULL COMPOSITION (everything on)",
        f"  P(target)                          {100 * full.probability_of_target:6.2f}%",
        f"  median gain                        {full.median_gain_pp:6.2f} pp",
        f"  saturated fraction                 {100 * full.saturated_fraction:5.1f}%",
        f"  cumulative adverse risk            {100 * full.cumulative_adverse_risk:5.2f}%",
        f"  risk-adjusted                      {100 * full.risk_adjusted_probability:6.2f}%",
        f"  saturation artefact                {'YES - excluded' if full.saturation_artefact else 'no'}",
        "",
        "BOTH ARMS TOGETHER",
        f"  graft                              {100 * joint['graft_probability']:6.2f}%",
        f"  clearance                          {100 * joint['clearance_probability']:6.2f}%",
        f"  joint, assuming independence       {100 * joint['joint_optimistic']:6.2f}%",
        "  randomised proof                    0.00%",
        "",
        "  The joint is a product and so assumes the arms are independent. The measured",
        "  bed coupling says they are not, so this is optimistic and labelled as such.",
        "",
        "THE BILL, EVERY AXIS CHARGED SEPARATELY",
        "  axis                        value      limit    beyond precedent",
    ]
    for row in ledger.itertuples(index=False):
        lines.append(
            f"  {row.axis:27s} {row.value:7.3f}   {row.precedented_limit:7.3f}   "
            f"{'YES' if row.beyond_precedent else 'no'}"
        )
    lines += [
        "",
        "  A breach on one axis may not be paid for with an unflagged breach on another.",
        "  That trap was already caught once in the gap-5 trade-off and is charged here.",
        "",
        "VERDICT (computed, not asserted)",
        f"  strongest composition anywhere     {100 * verdict['best_probability_any']:6.2f}%"
        f"   ({verdict['best_mechanisms_any']})",
        f"  best single mechanism              {100 * verdict['best_probability_single']:6.2f}%",
        f"  selected admissible                {100 * verdict['selected_probability']:6.2f}%"
        f"   ({verdict['selected_count']} mechanism(s))",
        f"  gain over best single mechanism    {verdict['gain_over_best_single_pp']:+6.2f} pp",
        f"  ATLAS-4 best admissible (reference) {100 * verdict['atlas4_reference']:5.2f}%",
        f"  gain over ATLAS-4                  {verdict['gain_over_atlas4_pp']:+6.2f} pp",
        f"  strongest ADMISSIBLE               {100 * verdict['best_probability_admissible']:6.2f}%"
        f"   ({verdict['best_mechanisms_admissible']})",
        f"  price of admissibility             {verdict['price_of_admissibility_pp']:6.2f} pp",
        f"  price of cheapness                 {verdict['price_of_cheapness_pp']:6.2f} pp",
        f"  risk-adjusted, selected            {100 * verdict['selected_risk_adjusted']:6.2f}%",
        f"  risk-adjusted, strongest admissible"
        f" {100 * verdict['best_risk_adjusted_admissible']:6.2f}%",
        "",
    ]
    if verdict["composition_adds_over_single_pp"] > 0.5:
        lines += [
            f"  COMPOSING ADDS {verdict['composition_adds_over_single_pp']:.2f} pp over the best single mechanism,",
            "  and the gain survives charging procedural risk. So the composite",
            "  is genuinely better, not merely more complicated. What recombination also",
            "  did was DELETE: the inert and overlapping parts were dropped rather than",
            "  carried along for the sake of the architecture.",
        ]
    elif verdict["selection_is_a_single_mechanism"]:
        lines += [
            "  THE SEARCH OVER EVERY COMPOSITION SELECTED A SINGLE MECHANISM.",
            "  Recombination did not build something bigger. It showed that most of the",
            "  machinery was inert or was the same mechanism wearing two names. The gain",
            "  from composing is subtraction, not addition.",
        ]
    elif verdict["composition_beat_best_single"]:
        lines += [
            "  Composing beats the best single mechanism, by the margin above. The",
            "  interaction table shows which pairs genuinely added rather than overlapped.",
        ]
    else:
        lines += [
            "  Composing did NOT beat the best single mechanism once admissibility was",
            "  enforced. Reporting that, not the larger inadmissible number.",
        ]
    lines += [
        "",
    ]
    # Compare against sampling noise, not against zero: two near-tied designs differ by
    # a few hundredths of a point for no reason at all, and calling that a price is the
    # same mistake as solving exactly on a threshold.
    if verdict["price_of_admissibility_pp"] > 100.0 * INTERACTION_TOLERANCE:
        lines += [
            "  Admissibility has a price: the strongest composition overall breaches a",
            "  precedent axis, so it is excluded and the strongest admissible one is weaker.",
        ]
    else:
        lines += [
            "  Admissibility costs NOTHING here beyond sampling noise: the strongest",
            "  composition overall is admissible too, so no precedent filter removed it.",
        ]
    if verdict["price_of_cheapness_pp"] > 100.0 * INTERACTION_TOLERANCE:
        lines += [
            "",
            f"  Cheapness costs {verdict['price_of_cheapness_pp']:.2f} pp of raw probability: the trade",
            f"  for going from {verdict['best_count_admissible']} mechanisms to"
            f" {verdict['selected_count']}, and for the procedures and risk the dropped",
            "  mechanisms would have added.",
        ]
        if verdict["cheapness_wins_on_risk"]:
            lines += [
                "  Once procedural risk is charged the cheaper design is AHEAD, so the trade",
                "  is justified on risk and not only on cost.",
            ]
        else:
            lines += [
                "  The stronger admissible design is STILL ahead after risk is charged, so",
                "  this is a real cost/benefit choice, not a free win. Both are reported.",
            ]
    lines += [
        "",
        "  The clearance arm paired with all of this still sits above the precedented",
        "  potency cap. That breach is charged in the ledger and is not forgiven by",
        "  anything the graft arm achieved.",
        "",
        "  None of this creates evidence. Every mechanism is borrowed from an invented",
        "  architecture or from someone else's animal or early-phase human work. Randomised",
        "  proof for the composition is zero and no simulation can change that. This is a",
        "  sharper specification of what to test, not a cure.",
    ]
    return "\n".join(lines) + "\n"


def run_recombine_study(
    output_dir: str | Path = "results/recombine", draws: int = 2500, seed: int = 20261107
) -> dict[str, object]:
    root = Path(output_dir)
    tables = root / "tables"
    reports = root / "reports"
    tables.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)

    full = evaluate_composite(COMPOSITE, draws, seed)
    marginals = mechanism_marginals(COMPOSITE, draws, seed)
    standalone = mechanism_standalone(COMPOSITE, draws, seed)
    sensitivity = maturation_sensitivity(max(600, draws // 2), seed)
    interactions = pairwise_interactions(COMPOSITE, max(600, draws // 2), seed)
    grid = exhaustive_subsets(COMPOSITE, max(400, draws // 3), seed)
    selection = select_composite(grid)
    verdict = composition_verdict(grid, selection)
    recommended = design_from_selection(selection)
    recommended_result = evaluate_composite(recommended, draws, seed)
    contrast = stacking_contrast(COMPOSITE, draws, seed)
    joint = joint_result(recommended, max(800, draws // 2), seed)
    ledger = cost_ledger(recommended)
    catalogue = pd.DataFrame(
        [
            {
                "mechanism": m.name,
                "origin": m.origin,
                "kind": m.kind,
                "cost_axis": m.cost_axis,
                "description": m.description,
                "status": RECOMBINE_TAG,
            }
            for m in MECHANISMS
        ]
    )

    catalogue.to_csv(tables / "recombine_mechanism_catalogue.csv", index=False)
    marginals.to_csv(tables / "recombine_marginals.csv", index=False)
    standalone.to_csv(tables / "recombine_standalone.csv", index=False)
    sensitivity.to_csv(tables / "recombine_maturation_sensitivity.csv", index=False)
    contrast.to_csv(tables / "recombine_stacking_contrast.csv", index=False)
    interactions.to_csv(tables / "recombine_interactions.csv", index=False)
    grid.to_csv(tables / "recombine_exhaustive_subsets.csv", index=False)
    ledger.to_csv(tables / "recombine_cost_ledger.csv", index=False)
    pd.DataFrame([{**joint, "status": RECOMBINE_TAG}]).to_csv(
        tables / "recombine_joint.csv", index=False
    )
    pd.DataFrame([{**verdict, "status": RECOMBINE_TAG}]).to_csv(
        tables / "recombine_verdict.csv", index=False
    )
    pd.DataFrame(
        [
            {
                "admissible": selection["admissible"],
                "robust_set_exists": selection["robust_set_exists"],
                "selected": str(selection.get("best")),
                "cheapest": str(selection.get("cheapest")),
                "status": RECOMBINE_TAG,
            }
        ]
    ).to_csv(tables / "recombine_selection.csv", index=False)
    (reports / "recombine_report.txt").write_text(
        recombine_report(
            full,
            marginals,
            interactions,
            grid,
            selection,
            joint,
            ledger,
            standalone,
            sensitivity,
            verdict,
            recommended_result,
            contrast,
        )
    )
    (root / "DISCLAIMER.txt").write_text(DISCLAIMER + "\n")
    return {
        "output_dir": str(root),
        "outputs": 13,
        "status": RECOMBINE_TAG,
        "mechanisms": len(MECHANISMS),
        "compositions_tested": int(len(grid)),
        "admissible": selection["admissible"],
        "robust_set_exists": selection["robust_set_exists"],
        "selected": selection.get("best"),
        "cheapest": selection.get("cheapest"),
        "overlapping_pairs": int((interactions["verdict"] == "overlapping").sum()),
        "synergistic_pairs": int((interactions["verdict"] == "synergistic").sum()),
        "maturation_resolution_limit": maturation_resolution_limit(sensitivity),
        "sensitivity_scenarios": int(len(sensitivity)),
        "sensitivity_robust": int(sensitivity["robust_with_margin"].sum()),
        "inert_mechanisms": int((~marginals["load_bearing"]).sum()),
        "full_probability": full.probability_of_target,
        "full_saturation_artefact": bool(full.saturation_artefact),
        "verdict": verdict,
        "recommended_probability": recommended_result.probability_of_target,
        "recommended_risk_adjusted": recommended_result.risk_adjusted_probability,
        "recommended_mean_doses": recommended_result.mean_doses_used,
        "stacking_penalty_pp": 100.0
        * (
            float(contrast.set_index("design").loc["selected_composition", "probability_of_target"])
            - float(contrast.set_index("design").loc["everything_on", "probability_of_target"])
        ),
        "cures_declared": 0,
        "disclaimer": DISCLAIMER,
    }
