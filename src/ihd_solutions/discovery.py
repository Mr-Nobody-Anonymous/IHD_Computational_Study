"""In-silico cure discovery: search for a *specification*, never a cure.

Why this module cannot output a cure
------------------------------------
The cure definition in :mod:`ihd_solutions.cure` has five criteria. One of them,
``randomised_proof``, requires randomised human outcome evidence. No simulation can
produce that. It is not a gap in this code; it is a property of the word "cure".

So this module searches the four **simulation-admissible** criteria and returns a
*candidate specification*: the smallest capability vector that would satisfy them,
plus a ranked list of the measurements that would most change the answer. The
output is labelled a specification everywhere and is quarantined from the evidence
tables exactly like :mod:`ihd_solutions.hypothesis`.

The central trick
-----------------
Search the same capability space twice, under two different bound sets:

* ``DEMONSTRATED`` - no capability exceeds what has actually been shown.
* ``CONCEIVABLE``  - generous bounds, allowing mechanisms nobody has built.

If the feasible region is empty under the first and non-empty under the second, the
*difference between the two boxes is precisely the discovery that is missing*. That
turns "what is missing?" into a number instead of an opinion.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum

import numpy as np
import pandas as pd

from ihd_solutions.cure import CRITERION_THRESHOLD
from ihd_solutions.units import require_fraction, require_positive

DISCLAIMER = (
    "SIMULATION-ONLY DISCOVERY - outputs a candidate specification, never a cure. "
    "The randomised_proof criterion is excluded because no simulation can satisfy "
    "it. Nothing here is evidence, a therapy, or a recommendation."
)

DISCOVERY_TAG = "in_silico_specification_not_evidence"

#: Criteria a simulation can legitimately reason about.
SIMULATION_ADMISSIBLE: tuple[str, ...] = (
    "durable_driver_control",
    "burden_reversal",
    "myocardial_restoration",
    "risk_normalisation",
)

#: The criterion no simulation can ever satisfy.
EXCLUDED_CRITERION = "randomised_proof"

#: LVEF points needed to normalise an ischaemic cardiomyopathy (25% -> 50%).
NORMALISATION_LVEF_GAIN_PP = 25.0

#: The four capabilities a specification is built from.
CAPABILITY_NAMES: tuple[str, ...] = (
    "driver_durability",
    "substrate_removal",
    "lvef_gain_pp",
    "hazard_ratio",
)


class BoundSet(str, Enum):
    DEMONSTRATED = "demonstrated"
    CONCEIVABLE = "conceivable"


@dataclass(frozen=True)
class CapabilityBounds:
    """Lower and upper bound for each capability, with provenance."""

    driver_durability: tuple[float, float]
    substrate_removal: tuple[float, float]
    lvef_gain_pp: tuple[float, float]
    hazard_ratio: tuple[float, float]
    provenance: dict[str, str] = field(default_factory=dict)

    def names(self) -> tuple[str, ...]:
        return CAPABILITY_NAMES

    def as_dict(self) -> dict[str, tuple[float, float]]:
        return {name: getattr(self, name) for name in self.names()}


#: What has actually been demonstrated. Substrate removal uses the best human
#: noncalcified-plaque result; LVEF uses BioVAT-HF; hazard ratio uses the best
#: strict stack in the most favourable population.
DEMONSTRATED_BOUNDS = CapabilityBounds(
    driver_durability=(0.0, 1.0),
    substrate_removal=(0.0, 0.431),
    lvef_gain_pp=(0.0, 3.9),
    hazard_ratio=(0.36, 1.0),
    provenance={
        "driver_durability": "one-time hepatic editing modelled as permanent (score 1.0)",
        "substrate_removal": "best human noncalcified plaque reduction 43.1% (regression grid)",
        "lvef_gain_pp": "BioVAT-HF, n=20, +3.9 LVEF points",
        "hazard_ratio": "best strict stack, stable_ccs: 13.93% -> 5.09% five-year risk",
    },
)

#: Generous bounds. Substrate removal is capped by the REVERT-1 architecture's own
#: asymptotic ceiling, so even the optimistic box respects a modelled limit.
CONCEIVABLE_BOUNDS = CapabilityBounds(
    driver_durability=(0.0, 1.0),
    substrate_removal=(0.0, 0.835),
    lvef_gain_pp=(0.0, 25.0),
    hazard_ratio=(0.10, 1.0),
    provenance={
        "driver_durability": "unchanged; already solved in principle",
        "substrate_removal": "REVERT-1 asymptotic ceiling 83.5% with perfect modules",
        "lvef_gain_pp": "full normalisation assumed possible; NO precedent supplies this",
        "hazard_ratio": "event-model floor MIN_SCENARIO_RR = 0.1",
    },
)


def bounds_for(which: BoundSet) -> CapabilityBounds:
    return DEMONSTRATED_BOUNDS if which is BoundSet.DEMONSTRATED else CONCEIVABLE_BOUNDS


def excess_risk_remaining(
    hazard_ratio: float, baseline_risk: float, reference_risk: float
) -> float:
    """Share of the baseline excess risk that survives treatment.

    Constant-hazard arithmetic over a common horizon, matching the event model and
    :func:`ihd_solutions.hypothesis.required_hazard_ratio`.
    """
    require_positive("hazard_ratio", hazard_ratio)
    require_fraction("baseline_risk", baseline_risk)
    require_fraction("reference_risk", reference_risk)
    if reference_risk >= baseline_risk:
        raise ValueError("reference risk must be below baseline risk")
    treated = 1.0 - (1.0 - baseline_risk) ** hazard_ratio
    return float(max(0.0, (treated - reference_risk) / (baseline_risk - reference_risk)))


@dataclass(frozen=True)
class Specification:
    """One point in capability space and the criterion scores it earns."""

    driver_durability: float
    substrate_removal: float
    lvef_gain_pp: float
    hazard_ratio: float
    baseline_risk: float = 0.1393
    reference_risk: float = 0.0198

    def scores(self) -> dict[str, float]:
        return {
            "durable_driver_control": float(self.driver_durability),
            "burden_reversal": float(min(1.0, max(0.0, self.substrate_removal))),
            "myocardial_restoration": float(
                min(1.0, max(0.0, self.lvef_gain_pp / NORMALISATION_LVEF_GAIN_PP))
            ),
            "risk_normalisation": float(
                min(
                    1.0,
                    max(
                        0.0,
                        1.0
                        - excess_risk_remaining(
                            self.hazard_ratio, self.baseline_risk, self.reference_risk
                        ),
                    ),
                )
            ),
        }

    def satisfied(self, threshold: float = CRITERION_THRESHOLD) -> tuple[str, ...]:
        return tuple(name for name, value in self.scores().items() if value >= threshold)

    def blocking(self, threshold: float = CRITERION_THRESHOLD) -> tuple[str, ...]:
        return tuple(name for name, value in self.scores().items() if value < threshold)

    def worst_score(self) -> float:
        return min(self.scores().values())

    @property
    def is_candidate(self) -> bool:
        """True when all four simulation-admissible criteria are met.

        Deliberately NOT called ``is_cure``: the proof criterion is unaddressed.
        """
        return not self.blocking()


def required_capabilities(
    baseline_risk: float = 0.1393,
    reference_risk: float = 0.0198,
    threshold: float = CRITERION_THRESHOLD,
) -> dict[str, float]:
    """Back-calculated minimum capability vector, solved analytically."""
    require_fraction("threshold", threshold, allow_one=True)
    require_fraction("baseline_risk", baseline_risk)
    require_fraction("reference_risk", reference_risk)
    if reference_risk >= baseline_risk:
        raise ValueError("reference risk must be below baseline risk")
    allowed_excess = 1.0 - threshold
    target_risk = reference_risk + allowed_excess * (baseline_risk - reference_risk)
    if not 0.0 < target_risk < 1.0:
        raise ValueError("target risk must lie strictly between 0 and 1")
    return {
        "driver_durability": float(threshold),
        "substrate_removal": float(threshold),
        "lvef_gain_pp": float(threshold * NORMALISATION_LVEF_GAIN_PP),
        "hazard_ratio": float(math.log1p(-target_risk) / math.log1p(-baseline_risk)),
    }


def _feasible_vector(bounds: CapabilityBounds) -> dict[str, float]:
    """Corner of a box that maximises every criterion simultaneously."""
    return {
        "driver_durability": bounds.driver_durability[1],
        "substrate_removal": bounds.substrate_removal[1],
        "lvef_gain_pp": bounds.lvef_gain_pp[1],
        "hazard_ratio": bounds.hazard_ratio[0],
    }


def best_achievable(
    which: BoundSet,
    baseline_risk: float = 0.1393,
    reference_risk: float = 0.0198,
) -> Specification:
    """Best corner of a bound set: the nearest this box can get to the target."""
    bounds = bounds_for(which)
    return Specification(
        **_feasible_vector(bounds),
        baseline_risk=baseline_risk,
        reference_risk=reference_risk,
    )


def capability_gap_table(
    baseline_risk: float = 0.1393,
    reference_risk: float = 0.0198,
    threshold: float = CRITERION_THRESHOLD,
) -> pd.DataFrame:
    """What is missing, per capability, as a ratio rather than an opinion."""
    required = required_capabilities(baseline_risk, reference_risk, threshold)
    bounds = bounds_for(BoundSet.DEMONSTRATED)
    best = _feasible_vector(bounds)
    rows = []
    for name in bounds.names():
        need, have = required[name], best[name]
        lower_is_better = name == "hazard_ratio"
        if lower_is_better:
            gap = have / need if need > 0 else math.inf
            reachable = bounds.as_dict()[name][0] <= need
        else:
            gap = need / have if have > 0 else math.inf
            reachable = bounds.as_dict()[name][1] >= need
        rows.append(
            {
                "capability": name,
                "required": need,
                "best_demonstrated": have,
                "gap_factor": gap,
                "lower_is_better": lower_is_better,
                "reachable_within_demonstrated_bounds": reachable,
                "conceivable_bound": (
                    bounds_for(BoundSet.CONCEIVABLE).as_dict()[name][0]
                    if lower_is_better
                    else bounds_for(BoundSet.CONCEIVABLE).as_dict()[name][1]
                ),
                "provenance": bounds.provenance.get(name, ""),
                "status": DISCOVERY_TAG,
            }
        )
    return pd.DataFrame(rows).sort_values("gap_factor", ascending=False).reset_index(drop=True)


def sample_capabilities(draws: int, bounds: CapabilityBounds, seed: int) -> dict[str, np.ndarray]:
    """Uniform sampling inside a capability box. Priors are deliberately flat."""
    require_positive("draws", draws)
    rng = np.random.default_rng(seed)
    samples = {}
    for name, (low, high) in bounds.as_dict().items():
        if high < low:
            raise ValueError(f"{name}: upper bound {high} below lower bound {low}")
        samples[name] = rng.uniform(low, high, size=int(draws))
    return samples


def feasibility_search(
    draws: int = 4000,
    seed: int = 20261107,
    baseline_risk: float = 0.1393,
    reference_risk: float = 0.0198,
    threshold: float = CRITERION_THRESHOLD,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Search both boxes; return per-draw scores and a per-box summary."""
    per_box = []
    all_rows = []
    for which in BoundSet:
        bounds = bounds_for(which)
        samples = sample_capabilities(draws, bounds, seed)
        frame = pd.DataFrame(samples)
        scores = pd.DataFrame(
            [
                Specification(
                    **row, baseline_risk=baseline_risk, reference_risk=reference_risk
                ).scores()
                for row in frame.to_dict("records")
            ]
        )
        passes = scores.ge(threshold)
        frame = pd.concat([frame, scores], axis=1)
        frame["criteria_met"] = passes.sum(axis=1)
        frame["is_candidate"] = passes.all(axis=1)
        frame["worst_score"] = scores.min(axis=1)
        frame["bound_set"] = which.value
        frame["status"] = DISCOVERY_TAG
        all_rows.append(frame)
        corner = best_achievable(which, baseline_risk, reference_risk)
        per_box.append(
            {
                "bound_set": which.value,
                "draws": int(draws),
                "p_candidate": float(frame["is_candidate"].mean()),
                "any_candidate": bool(frame["is_candidate"].any()),
                "corner_is_candidate": corner.is_candidate,
                "corner_worst_score": corner.worst_score(),
                "corner_blocking": ";".join(corner.blocking(threshold)) or "none",
                "median_criteria_met": float(frame["criteria_met"].median()),
                "status": DISCOVERY_TAG,
            }
        )
    return pd.concat(all_rows, ignore_index=True), pd.DataFrame(per_box)


def feasible_volume_fraction(
    which: BoundSet,
    baseline_risk: float = 0.1393,
    reference_risk: float = 0.0198,
    threshold: float = CRITERION_THRESHOLD,
) -> tuple[float, pd.DataFrame]:
    """Exact feasible volume of a capability box under independent flat priors.

    For a box, no Monte Carlo is needed: the passing fraction is the product of the
    per-capability window fractions. Sampling a region this thin is pure noise, so
    the analytic value is the one to quote.
    """
    bounds = bounds_for(which)
    required = required_capabilities(baseline_risk, reference_risk, threshold)
    rows = []
    product = 1.0
    for name, (low, high) in bounds.as_dict().items():
        span = high - low
        if span <= 0:
            raise ValueError(f"{name}: degenerate bounds")
        if name == "hazard_ratio":
            window = max(0.0, min(high, required[name]) - low)
        elif name == "lvef_gain_pp":
            window = max(0.0, high - max(low, required[name]))
        else:
            # Score equals the capability for durability and substrate removal.
            window = max(0.0, high - max(low, required[name]))
        fraction = window / span
        product *= fraction
        rows.append(
            {
                "capability": name,
                "bound_low": low,
                "bound_high": high,
                "required": required[name],
                "passing_window": window,
                "passing_fraction": fraction,
                "bound_set": which.value,
                "status": DISCOVERY_TAG,
            }
        )
    return float(product), pd.DataFrame(rows)


def decision_value_of_information(
    samples: pd.DataFrame, outcome: str = "worst_score", bins: int = 10
) -> pd.DataFrame:
    """Share of the decision margin's variance each capability resolves.

    First-order Sobol-style index: Var_i(E[Y | X_i]) / Var(Y). The default outcome is
    ``worst_score`` - the binding criterion's score, i.e. the margin that decides
    whether a specification exists. A capability scoring high is one where learning
    the truth would most move that margin, which is what ranks experiments.

    ``is_candidate`` may be passed instead, but a feasible region this thin usually
    yields no positive draws and therefore zero variance; prefer the margin.
    """
    require_positive("bins", bins)
    if outcome not in samples.columns:
        raise ValueError(f"samples must carry a {outcome!r} column")
    values_out = samples[outcome].to_numpy(dtype=float)
    total_variance = float(np.var(values_out))
    names = [name for name in CAPABILITY_NAMES if name in samples.columns]
    rows = []
    for name in names:
        values = samples[name].to_numpy(dtype=float)
        if np.ptp(values) == 0:
            rows.append(
                {
                    "capability": name,
                    "variance_share": 0.0,
                    "outcome": outcome,
                    "status": DISCOVERY_TAG,
                }
            )
            continue
        edges = np.quantile(values, np.linspace(0.0, 1.0, int(bins) + 1))
        edges[0], edges[-1] = -np.inf, np.inf
        index = np.digitize(values, edges[1:-1])
        explained = 0.0
        for group in np.unique(index):
            mask = index == group
            weight = float(mask.mean())
            explained += weight * (float(values_out[mask].mean()) - float(values_out.mean())) ** 2
        share = explained / total_variance if total_variance > 0 else 0.0
        rows.append(
            {
                "capability": name,
                "variance_share": float(min(1.0, max(0.0, share))),
                "outcome": outcome,
                "status": DISCOVERY_TAG,
            }
        )
    return pd.DataFrame(rows).sort_values("variance_share", ascending=False).reset_index(drop=True)


def capability_ablation(
    baseline_risk: float = 0.1393,
    reference_risk: float = 0.0198,
    threshold: float = CRITERION_THRESHOLD,
) -> pd.DataFrame:
    """Start from a passing specification, revert one capability to demonstrated.

    Whatever breaks when reverted is indispensable. Whatever survives was already
    solved and is not part of the discovery problem.
    """
    bounds = bounds_for(BoundSet.CONCEIVABLE)
    ideal = _feasible_vector(bounds)
    demonstrated = _feasible_vector(bounds_for(BoundSet.DEMONSTRATED))
    rows = []
    for name in bounds.names():
        vector = {**ideal, name: demonstrated[name]}
        spec = Specification(**vector, baseline_risk=baseline_risk, reference_risk=reference_risk)
        rows.append(
            {
                "reverted_capability": name,
                "reverted_to": demonstrated[name],
                "still_candidate": spec.is_candidate,
                "indispensable": not spec.is_candidate,
                "blocking_after_revert": ";".join(spec.blocking(threshold)) or "none",
                "worst_score": spec.worst_score(),
                "status": DISCOVERY_TAG,
            }
        )
    return pd.DataFrame(rows)


@dataclass(frozen=True)
class DiscoveryStage:
    """One stage of the simulation-only programme."""

    number: int
    name: str
    question: str
    method: str
    exit_gate: str
    failure_means: str
    executable_now: bool


DISCOVERY_STAGES: tuple[DiscoveryStage, ...] = (
    DiscoveryStage(
        0,
        "Formalise the target",
        "What exactly would count, and which criteria can a simulation address?",
        "Freeze the five criteria; mark randomised_proof as simulation-inadmissible; "
        "set population-specific risk targets",
        "A written target with four admissible criteria and one explicitly excluded",
        "Nothing downstream is interpretable; any 'cure' claim would be definitional fraud",
        True,
    ),
    DiscoveryStage(
        1,
        "Feasibility frontier",
        "Does any capability vector satisfy the admissible criteria?",
        "Search demonstrated and conceivable boxes; compare the feasible regions",
        "Non-empty feasible region under conceivable bounds, empty under demonstrated",
        "If empty under conceivable bounds too, the architecture is wrong - redesign, do not tune",
        True,
    ),
    DiscoveryStage(
        2,
        "Gap quantification",
        "How far is each capability from the requirement?",
        "Back-calculate the minimum vector analytically; divide by best demonstrated",
        "Every gap expressed as a factor with a cited provenance",
        "A gap without provenance is a guess and cannot prioritise work",
        True,
    ),
    DiscoveryStage(
        3,
        "Value of information",
        "Which unknown, if measured, would most change the answer?",
        "Sobol-style variance decomposition of the binary feasibility decision",
        "Capabilities ranked by decision variance resolved",
        "Experiments get chosen by fashion instead of decision impact",
        True,
    ),
    DiscoveryStage(
        4,
        "Minimal capability set",
        "Which capabilities are indispensable versus already solved?",
        "Revert each capability to its demonstrated value inside a passing vector",
        "An irreducible set of required discoveries",
        "Effort gets spent re-solving what is already solved",
        True,
    ),
    DiscoveryStage(
        5,
        "Adversarial falsification",
        "Does the candidate survive hostile but plausible assumptions?",
        "Re-run under pessimistic couplings, added PK, withdrawal and recurrence, "
        "spatial cap mechanics; sweep the REVERT-1 kill criteria",
        "Candidate retains a feasible region across structural perturbations",
        "The specification is an artefact of one model structure and must be discarded",
        False,
    ),
    DiscoveryStage(
        6,
        "Identifiability audit",
        "Are the governing parameters identifiable from obtainable data?",
        "Profile likelihood / sloppy-model analysis on the coupled clearance model; "
        "flag parameters no feasible measurement can pin",
        "Every decision-relevant parameter has an identifiable measurement route",
        "The programme is unfalsifiable; stop and reformulate",
        False,
    ),
    DiscoveryStage(
        7,
        "Specification handoff",
        "What is the cheapest real-world sequence that could convert this to evidence?",
        "Emit a target product profile, kill criteria and an ordered experiment list",
        "A specification a laboratory can act on, with the proof criterion still open",
        "The work stays inert modelling",
        True,
    ),
)


def stages_frame() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "stage": stage.number,
                "name": stage.name,
                "question": stage.question,
                "method": stage.method,
                "exit_gate": stage.exit_gate,
                "failure_means": stage.failure_means,
                "executable_now": stage.executable_now,
                "status": DISCOVERY_TAG,
            }
            for stage in DISCOVERY_STAGES
        ]
    )


def discovery_report(
    summary: pd.DataFrame,
    gaps: pd.DataFrame,
    voi: pd.DataFrame,
    ablation: pd.DataFrame,
    required: dict[str, float],
) -> str:
    lines = [
        "In-silico cure discovery: searching for a specification, not a cure",
        "=" * 68,
        DISCLAIMER,
        "",
        "WHY THIS CANNOT RETURN A CURE",
        f"The '{EXCLUDED_CRITERION}' criterion requires randomised human outcome evidence.",
        "No simulation can generate it. This programme therefore searches the four",
        "simulation-admissible criteria and returns a candidate specification. Calling",
        "that a cure would be a definitional trick, not a discovery.",
        "",
        "REQUIRED CAPABILITY VECTOR (back-calculated)",
    ]
    for name, value in required.items():
        lines.append(f"  {name:22s} {value:.3f}")
    lines.extend(["", "FEASIBILITY, TWO BOUND SETS"])
    for row in summary.itertuples(index=False):
        lines.append(
            f"  {row.bound_set:14s} P(candidate)={row.p_candidate:6.2%}  "
            f"best corner worst-score {row.corner_worst_score:.3f}  "
            f"blocking: {row.corner_blocking}"
        )
    lines.extend(
        [
            "",
            "  Read this pair together: an empty feasible region under demonstrated",
            "  bounds and a non-empty one under conceivable bounds means the difference",
            "  between the boxes IS the missing discovery. A corner that passes while",
            "  P(candidate) is near zero means the feasible region is a vanishingly thin",
            "  corner: essentially every capability must be near its maximum at once.",
            "",
            "EXACT FEASIBLE VOLUME (analytic, no sampling noise)",
        ]
    )
    for which in BoundSet:
        volume, windows = feasible_volume_fraction(which)
        lines.append(f"  {which.value:14s} feasible volume = {volume:.3e}")
        for row in windows.itertuples(index=False):
            lines.append(
                f"      {row.capability:22s} passing window {row.passing_fraction:7.2%} of its range"
            )
    lines.extend(["", "WHAT IS MISSING, PER CAPABILITY"])
    for row in gaps.itertuples(index=False):
        reach = "reachable" if row.reachable_within_demonstrated_bounds else "OUT OF REACH"
        lines.append(
            f"  {row.capability:22s} need {row.required:7.3f}  have {row.best_demonstrated:7.3f}"
            f"  gap x{row.gap_factor:5.2f}  {reach}"
        )
        lines.append(f"      basis: {row.provenance}")
    lines.extend(
        [
            "",
            "WHICH MEASUREMENT WOULD MOST MOVE THE BINDING MARGIN",
            "(variance decomposition of the worst-criterion score)",
        ]
    )
    for row in voi.itertuples(index=False):
        lines.append(
            f"  {row.capability:22s} resolves {row.variance_share:6.1%} of decision variance"
        )
    lines.extend(["", "INDISPENSABLE CAPABILITIES (revert-one ablation)"])
    for row in ablation.itertuples(index=False):
        verdict = "INDISPENSABLE" if row.indispensable else "already solved"
        lines.append(
            f"  revert {row.reverted_capability:22s} -> {verdict}"
            + (f" (blocks {row.blocking_after_revert})" if row.indispensable else "")
        )
    lines.extend(["", "THE PROGRAMME"])
    for stage in DISCOVERY_STAGES:
        flag = "executable now" if stage.executable_now else "requires new model structure"
        lines.append(f"  {stage.number}. {stage.name} [{flag}]")
        lines.append(f"     question: {stage.question}")
        lines.append(f"     exit gate: {stage.exit_gate}")
    lines.extend(
        [
            "",
            "THE HONEST BOTTOM LINE",
            "- Driver durability is already solved in principle by one-time editing.",
            "- Substrate clearance is short by roughly a factor of two, and the modelled",
            "  architecture's own ceiling only just clears the bar.",
            "- Myocardial restoration is the binding discovery: nothing demonstrated is",
            "  close, and no parameter tuning creates a capability that does not exist.",
            "- Randomised proof stays open no matter how good the simulation gets.",
            "- Therefore: a cure cannot be discovered in simulation. A specification can,",
            "  and this is it.",
        ]
    )
    return "\n".join(lines) + "\n"


def run_discovery_study(
    output_dir,
    draws: int = 4000,
    seed: int = 20261107,
    baseline_risk: float = 0.1393,
    reference_risk: float = 0.0198,
) -> dict:
    """Write the discovery artefacts to their own quarantined directory."""
    from pathlib import Path

    root = Path(output_dir)
    (root / "tables").mkdir(parents=True, exist_ok=True)
    (root / "reports").mkdir(parents=True, exist_ok=True)

    samples, summary = feasibility_search(
        draws=draws, seed=seed, baseline_risk=baseline_risk, reference_risk=reference_risk
    )
    gaps = capability_gap_table(baseline_risk, reference_risk)
    volumes = []
    for which in BoundSet:
        volume, windows = feasible_volume_fraction(which, baseline_risk, reference_risk)
        windows["feasible_volume"] = volume
        volumes.append(windows)
    volume_table = pd.concat(volumes, ignore_index=True)
    conceivable = samples.loc[samples["bound_set"].eq(BoundSet.CONCEIVABLE.value)]
    voi = decision_value_of_information(conceivable)
    ablation = capability_ablation(baseline_risk, reference_risk)
    required = required_capabilities(baseline_risk, reference_risk)

    summary.to_csv(root / "tables" / "discovery_feasibility.csv", index=False)
    gaps.to_csv(root / "tables" / "discovery_capability_gaps.csv", index=False)
    voi.to_csv(root / "tables" / "discovery_value_of_information.csv", index=False)
    ablation.to_csv(root / "tables" / "discovery_capability_ablation.csv", index=False)
    stages_frame().to_csv(root / "tables" / "discovery_stages.csv", index=False)
    volume_table.to_csv(root / "tables" / "discovery_feasible_volume.csv", index=False)
    samples.to_csv(root / "tables" / "discovery_samples.csv", index=False)
    pd.DataFrame([{**required, "status": DISCOVERY_TAG}]).to_csv(
        root / "tables" / "discovery_required_vector.csv", index=False
    )
    (root / "reports" / "discovery_report.txt").write_text(
        discovery_report(summary, gaps, voi, ablation, required), encoding="utf-8"
    )
    (root / "DISCLAIMER.txt").write_text(DISCLAIMER + "\n", encoding="utf-8")

    demonstrated = summary.loc[summary["bound_set"].eq(BoundSet.DEMONSTRATED.value)].iloc[0]
    return {
        "status": "specification_only",
        "disclaimer": DISCLAIMER,
        "output_dir": str(root),
        "excluded_criterion": EXCLUDED_CRITERION,
        "p_candidate_demonstrated": float(demonstrated["p_candidate"]),
        "p_candidate_conceivable": float(
            summary.loc[summary["bound_set"].eq(BoundSet.CONCEIVABLE.value)].iloc[0]["p_candidate"]
        ),
        "indispensable": list(ablation.loc[ablation["indispensable"], "reverted_capability"]),
        "feasible_volume_conceivable": feasible_volume_fraction(
            BoundSet.CONCEIVABLE, baseline_risk, reference_risk
        )[0],
        "feasible_volume_demonstrated": feasible_volume_fraction(
            BoundSet.DEMONSTRATED, baseline_risk, reference_risk
        )[0],
        "outputs": 10,
    }
