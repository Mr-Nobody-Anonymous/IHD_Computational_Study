"""Best-of-breed composite: assemble the strongest part per disease node, then test it.

The trap this module is built to avoid
--------------------------------------
"Take the best part of each and combine them" invites two errors:

1. **Summing partial capabilities.** Two half-effective regeneration therapies do not
   make one restored heart. Coverage of a criterion combines by *maximum*, matching
   :mod:`ihd_solutions.cure`.
2. **Double-counting one target.** PCSK9 base editing and a PCSK9 antibody hit the same
   molecule. The catalogue puts them in different ``exclusive_group`` values, so the
   catalogue alone will happily stack them. This module adds an explicit target-overlap
   rule on top.

And a third rule that the single-agent scoring did not need:

3. **An untested combination earns no randomised-proof credit.** Every component can
   have phase-3 outcome evidence and the *combination* still has never been trialled.
   Component proof does not transfer. The composite therefore scores
   ``randomised_proof = 0.0`` with that reason stated, which is why no composite built
   here can be declared a cure.

Selection is greedy by *marginal evidence-graded event benefit* in strict mode, under
node and target constraints - so "best" means measured incremental benefit in the
repository's own event model, not a hand-picked favourite.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

import numpy as np
import pandas as pd

from ihd_solutions.catalog import Regimen, Therapy, load_catalog, validate_combination
from ihd_solutions.cure import (
    BURDEN_REVERSAL_EVIDENCE,
    CRITERION_THRESHOLD,
    MYOCARDIAL_RESTORATION_EVIDENCE,
)
from ihd_solutions.evidence import CreditMode, TrialAnchor, load_anchors
from ihd_solutions.population import Population, load_populations
from ihd_solutions.risk import simulate_events
from ihd_solutions.units import (
    DAYS_PER_YEAR,
    require_nonnegative,
    require_positive,
    time_grid,
)

DISCLAIMER = (
    "COMPOSITE SCENARIO - a best-per-node assembly, not a tested regimen. Every "
    "combination here is untrialled, so randomised proof is zero by construction. "
    "Invented components are tagged separately and never earn evidence credit."
)

COMPOSITE_TAG = "untrialled_composite_scenario"
INVENTED_TAG = "invented_component_unvalidated"

#: Disease nodes, each mapped to the catalogue exclusive groups that can serve it.
#: One component per node: four anti-inflammatories are not four cures.
DISEASE_NODES: dict[str, tuple[str, ...]] = {
    "durable_driver": ("pcsk9_gene", "angptl3_gene", "pcsk9_vaccine"),
    "ldl_synthesis": ("statin", "bempedoic"),
    "ldl_absorption": ("ezetimibe",),
    "residual_pcsk9": ("pcsk9_protein",),
    "ldl_ceil": ("cetp",),
    "inflammation": ("colchicine", "il1_axis", "il6_axis", "nlrp3", "mtx"),
    "thrombosis": ("antiplatelet",),
    "hemodynamic": ("bp",),
    "metabolic": ("incretin", "sGLT2"),
    "lipoprotein_a": ("lpa_agent",),
    "triglyceride": ("triglyceride_agent", "apoc3"),
    "lifestyle": ("diet", "tobacco", "rehab"),
    "myocardium": ("regenerative", "cardiac_gene_therapy", "cd34"),
    "perfusion": ("revascularization", "cs_reducer"),
    "recurrence_tolerance": ("treg_tolerance",),
}

#: Exclusive groups that hit the same molecular target and must never be stacked,
#: even though the catalogue permits it.
TARGET_OVERLAPS: tuple[frozenset[str], ...] = (
    frozenset({"pcsk9_gene", "pcsk9_protein", "pcsk9_vaccine"}),
    frozenset({"statin", "bempedoic"}),
)

#: Nodes whose only available components are invented rather than catalogued.
INVENTED_NODES: dict[str, str] = {
    "plaque_substrate": "REVERT-1 coupled clearance loop (M1-M6, AND-gated)",
    "recurrence_tolerance": "ApoB/oxLDL CAR-Treg durability lock (S2)",
}


#: Criterion order used throughout, matching ihd_solutions.cure.
CRITERION_ORDER: tuple[str, ...] = (
    "durable_driver_control",
    "burden_reversal",
    "myocardial_restoration",
    "risk_normalisation",
    "randomised_proof",
)


class Variant(str, Enum):
    EVIDENCE_ONLY = "evidence_only"
    WITH_INVENTED = "with_invented"


def node_of(exclusive_group: str) -> str | None:
    for node, groups in DISEASE_NODES.items():
        if exclusive_group in groups:
            return node
    return None


def conflicts(
    chosen: tuple[str, ...], candidate: Therapy, catalog: dict[str, Therapy]
) -> str | None:
    """Return the reason a candidate cannot join, or None if it can."""
    group = candidate.exclusive_group
    for therapy_id in chosen:
        other = catalog[therapy_id]
        if other.exclusive_group == group:
            return f"shares exclusive group {group!r} with {therapy_id}"
    node = node_of(group)
    if node is not None:
        for therapy_id in chosen:
            if node_of(catalog[therapy_id].exclusive_group) == node:
                return f"node {node!r} already served by {therapy_id}"
    for overlap in TARGET_OVERLAPS:
        if group in overlap:
            for therapy_id in chosen:
                if catalog[therapy_id].exclusive_group in overlap:
                    return f"same molecular target as {therapy_id} (overlap set {sorted(overlap)})"
    return None


@dataclass(frozen=True)
class SelectionStep:
    """One greedy acceptance, with the marginal benefit that justified it."""

    step: int
    therapy_id: str
    node: str
    marginal_events_avoided_per_1000: float
    incidence_pct: float
    rejected_candidates: int


@dataclass(frozen=True)
class Composite:
    """A node-disjoint assembly and the trail of decisions that built it."""

    variant: Variant
    population_id: str
    therapy_ids: tuple[str, ...]
    steps: tuple[SelectionStep, ...]
    invented_components: tuple[str, ...] = ()
    coverage_additions: tuple[str, ...] = ()
    notes: tuple[str, ...] = field(default=(DISCLAIMER,))

    @property
    def regimen(self) -> Regimen:
        return Regimen(self.therapy_ids)

    @property
    def label(self) -> str:
        return " + ".join(self.therapy_ids) if self.therapy_ids else "(none)"


def _incidence(
    therapy_ids: tuple[str, ...],
    catalog: dict[str, Therapy],
    anchors: dict[str, TrialAnchor],
    population: Population,
    days: np.ndarray,
    mode: CreditMode,
) -> float:
    if not therapy_ids:
        # An empty regimen is not a valid Regimen; the untreated baseline comes
        # straight from the population's constant hazard, matching risk.py.
        years = float(days[-1] - days[0]) / DAYS_PER_YEAR
        return float(-np.expm1(-population.baseline_annual_mace_rate * years))
    regimen = validate_combination(catalog, therapy_ids)
    result = simulate_events(regimen, catalog, population, anchors, days, mode=mode)
    return float(result.cumulative_incidence[-1])


#: Nodes that matter for the cure criteria but are event-neutral in this evidence
#: framework, so greedy event selection can never reach them. Regeneration carries
#: zero event credit by design; that must not silently drop the myocardium node.
CRITERION_NODES: dict[str, dict[str, tuple[float, str]]] = {
    "myocardium": MYOCARDIAL_RESTORATION_EVIDENCE,
    "ldl_absorption": BURDEN_REVERSAL_EVIDENCE,
}


def complete_cure_nodes(
    chosen: tuple[str, ...], catalog: dict[str, Therapy]
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Attach the best component for cure-relevant, event-neutral nodes.

    Returns the completed tuple and the ids added this way. These additions carry
    zero event credit - they are chosen for criterion coverage, and recorded as such
    so they can never be mistaken for event benefit.
    """
    completed = chosen
    added: list[str] = []
    for node, evidence in sorted(CRITERION_NODES.items()):
        served = any(node_of(catalog[tid].exclusive_group) == node for tid in completed)
        if served:
            continue
        candidates = [
            (evidence[tid][0], tid)
            for tid, therapy in catalog.items()
            if tid in evidence
            and node_of(therapy.exclusive_group) == node
            and conflicts(completed, therapy, catalog) is None
        ]
        if not candidates:
            continue
        _, best_id = max(candidates)
        completed = (*completed, best_id)
        added.append(best_id)
    return completed, tuple(added)


def assemble_composite(
    population: Population,
    catalog: dict[str, Therapy] | None = None,
    anchors: dict[str, TrialAnchor] | None = None,
    days: np.ndarray | None = None,
    max_components: int = 8,
    min_marginal_per_1000: float = 0.5,
    mode: CreditMode = CreditMode.STRICT,
    complete_nodes: bool = True,
) -> Composite:
    """Greedy node-disjoint assembly ranked by marginal evidence-graded benefit.

    Selection runs in the requested credit mode (strict by default, so the assembly
    is built only from evidence that can move an event curve). Candidates that add
    nothing are never accepted, which is how neutral and mechanistic agents drop out
    on their own rather than by hand.
    """
    require_positive("max_components", max_components)
    require_nonnegative("min_marginal_per_1000", min_marginal_per_1000)
    catalog = load_catalog() if catalog is None else catalog
    anchors = load_anchors() if anchors is None else anchors
    days = time_grid(5.0, 30.0) if days is None else days

    chosen: tuple[str, ...] = ()
    current = _incidence(chosen, catalog, anchors, population, days, mode)
    steps: list[SelectionStep] = []
    for step in range(1, int(max_components) + 1):
        best_id: str | None = None
        best_incidence = current
        rejected = 0
        for therapy_id, therapy in sorted(catalog.items()):
            if conflicts(chosen, therapy, catalog) is not None:
                rejected += 1
                continue
            trial = _incidence((*chosen, therapy_id), catalog, anchors, population, days, mode)
            if trial < best_incidence - 1e-12:
                best_id, best_incidence = therapy_id, trial
        if best_id is None:
            break
        marginal = (current - best_incidence) * 1000.0
        if marginal < min_marginal_per_1000:
            break
        chosen = (*chosen, best_id)
        steps.append(
            SelectionStep(
                step,
                best_id,
                node_of(catalog[best_id].exclusive_group) or "unmapped",
                marginal,
                best_incidence * 100.0,
                rejected,
            )
        )
        current = best_incidence

    added: tuple[str, ...] = ()
    if complete_nodes:
        chosen, added = complete_cure_nodes(chosen, catalog)
    if chosen:
        validate_combination(catalog, chosen)
    notes = [DISCLAIMER]
    if added:
        notes.append("criterion-coverage additions (zero event credit): " + ", ".join(added))
    if not chosen:
        notes.append(
            "no catalogued agent earns event credit in this population and mode; "
            "an empty composite is the honest result"
        )
    # An empty composite is a real result: in some populations no catalogued agent
    # earns strict-mode event credit, and inventing one would be the error.
    return Composite(
        Variant.EVIDENCE_ONLY,
        population.population_id,
        chosen,
        tuple(steps),
        coverage_additions=added,
        notes=tuple(notes),
    )


def with_invented(composite: Composite) -> Composite:
    """Attach the invented nodes the catalogue cannot serve at all."""
    invented = tuple(
        f"{node}: {description}" for node, description in sorted(INVENTED_NODES.items())
    )
    return Composite(
        Variant.WITH_INVENTED,
        composite.population_id,
        composite.therapy_ids,
        composite.steps,
        invented,
        composite.coverage_additions,
        (
            DISCLAIMER,
            "invented components carry no evidence credit and do not move the event curve",
        ),
    )


def score_criteria(
    composite: Composite,
    catalog: dict[str, Therapy],
    anchors: dict[str, TrialAnchor],
    population: Population,
    days: np.ndarray,
    reference_risk: float = 0.0198,
    mode: CreditMode = CreditMode.STRICT,
) -> dict[str, float | str]:
    """Score the composite. Coverage by maximum; proof is zero because untested."""
    baseline = _incidence((), catalog, anchors, population, days, CreditMode.STRICT)
    treated = _incidence(composite.therapy_ids, catalog, anchors, population, days, mode)

    driver = 0.0
    for therapy_id in composite.therapy_ids:
        if catalog[therapy_id].exclusive_group in {"pcsk9_gene", "angptl3_gene"}:
            driver = 1.0
    reversal = max(
        (BURDEN_REVERSAL_EVIDENCE.get(tid, (0.0, ""))[0] for tid in composite.therapy_ids),
        default=0.0,
    )
    restoration = max(
        (MYOCARDIAL_RESTORATION_EVIDENCE.get(tid, (0.0, ""))[0] for tid in composite.therapy_ids),
        default=0.0,
    )
    if baseline > reference_risk:
        normalisation = max(0.0, 1.0 - (treated - reference_risk) / (baseline - reference_risk))
    else:
        normalisation = 0.0

    invented_bonus = composite.variant is Variant.WITH_INVENTED
    return {
        "durable_driver_control": driver,
        "burden_reversal": float(reversal),
        "myocardial_restoration": float(restoration),
        "risk_normalisation": float(min(1.0, normalisation)),
        # The combination has never been trialled. Component proof does not transfer.
        "randomised_proof": 0.0,
        "baseline_risk_pct": baseline * 100.0,
        "treated_risk_pct": treated * 100.0,
        "invented_nodes_present": invented_bonus,
        "proof_reason": "no trial of this combination; component evidence does not transfer",
    }


def leave_one_out(
    composite: Composite,
    catalog: dict[str, Therapy],
    anchors: dict[str, TrialAnchor],
    population: Population,
    days: np.ndarray,
    mode: CreditMode = CreditMode.STRICT,
) -> pd.DataFrame:
    """Remove one component at a time: what each part actually contributes."""
    full = _incidence(composite.therapy_ids, catalog, anchors, population, days, mode)
    rows = []
    for therapy_id in composite.therapy_ids:
        reduced = tuple(tid for tid in composite.therapy_ids if tid != therapy_id)
        without = _incidence(reduced, catalog, anchors, population, days, mode)
        rows.append(
            {
                "population_id": population.population_id,
                "removed": therapy_id,
                "node": node_of(catalog[therapy_id].exclusive_group) or "unmapped",
                "mode": mode.value,
                "incidence_without_pct": without * 100.0,
                "incidence_full_pct": full * 100.0,
                "contribution_per_1000": (without - full) * 1000.0,
                "indispensable": (without - full) * 1000.0 >= 0.5,
                "status": COMPOSITE_TAG,
            }
        )
    return (
        pd.DataFrame(rows)
        .sort_values("contribution_per_1000", ascending=False)
        .reset_index(drop=True)
    )


def evaluate_matrix(
    catalog: dict[str, Therapy] | None = None,
    anchors: dict[str, TrialAnchor] | None = None,
    populations: dict[str, Population] | None = None,
    population_ids: tuple[str, ...] | None = None,
    horizon_years: float = 5.0,
    max_components: int = 8,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Assemble and test across every population and every credit mode.

    Assembly happens *in the mode being evaluated*. A strict composite cannot contain
    gene editing because editing has no outcomes evidence; a surrogate composite can.
    Reporting one composite for all three modes would hide that.
    """
    catalog = load_catalog() if catalog is None else catalog
    anchors = load_anchors() if anchors is None else anchors
    populations = load_populations() if populations is None else populations
    ids = tuple(sorted(populations)) if population_ids is None else population_ids
    days = time_grid(horizon_years, 30.0)

    assembly_rows: list[dict] = []
    matrix_rows: list[dict] = []
    criteria_rows: list[dict] = []
    for population_id in ids:
        if population_id not in populations:
            raise ValueError(f"unknown population: {population_id}")
        population = populations[population_id]
        baseline = _incidence((), catalog, anchors, population, days, CreditMode.STRICT)
        for mode in CreditMode:
            composite = assemble_composite(
                population, catalog, anchors, days, max_components=max_components, mode=mode
            )
            for step in composite.steps:
                assembly_rows.append(
                    {
                        "population_id": population_id,
                        "mode": mode.value,
                        "step": step.step,
                        "therapy_id": step.therapy_id,
                        "node": step.node,
                        "marginal_events_avoided_per_1000": step.marginal_events_avoided_per_1000,
                        "incidence_pct": step.incidence_pct,
                        "rejected_candidates": step.rejected_candidates,
                        "status": COMPOSITE_TAG,
                    }
                )
            treated = _incidence(composite.therapy_ids, catalog, anchors, population, days, mode)
            matrix_rows.append(
                {
                    "population_id": population_id,
                    "mode": mode.value,
                    "components": len(composite.therapy_ids),
                    "therapy_ids": ";".join(composite.therapy_ids) or "(none)",
                    "coverage_additions": ";".join(composite.coverage_additions) or "(none)",
                    "baseline_risk_pct": baseline * 100.0,
                    "composite_risk_pct": treated * 100.0,
                    "events_avoided_per_1000": (baseline - treated) * 1000.0,
                    "relative_risk_reduction_pct": (
                        (1.0 - treated / baseline) * 100.0 if baseline > 0 else 0.0
                    ),
                    "excess_risk_remaining_fraction": (
                        max(0.0, (treated - 0.0198) / (baseline - 0.0198))
                        if baseline > 0.0198
                        else 0.0
                    ),
                    "status": COMPOSITE_TAG,
                }
            )
            for variant, candidate in (
                (Variant.EVIDENCE_ONLY, composite),
                (Variant.WITH_INVENTED, with_invented(composite)),
            ):
                scores = score_criteria(candidate, catalog, anchors, population, days, mode=mode)
                met = [
                    name for name in CRITERION_ORDER if float(scores[name]) >= CRITERION_THRESHOLD
                ]
                criteria_rows.append(
                    {
                        "population_id": population_id,
                        "mode": mode.value,
                        "variant": variant.value,
                        **{name: float(scores[name]) for name in CRITERION_ORDER},
                        "criteria_met": len(met),
                        "coverage_additions": ";".join(candidate.coverage_additions) or "(none)",
                        "satisfied": ";".join(met) or "none",
                        "is_cure": False,
                        "proof_reason": scores["proof_reason"],
                        "status": (
                            INVENTED_TAG if variant is Variant.WITH_INVENTED else COMPOSITE_TAG
                        ),
                    }
                )
    return (
        pd.DataFrame(assembly_rows),
        pd.DataFrame(matrix_rows),
        pd.DataFrame(criteria_rows),
    )


def component_frequency(assembly: pd.DataFrame) -> pd.DataFrame:
    """Which parts get chosen across populations - the recurring best-of-breed core."""
    if assembly.empty:
        return pd.DataFrame(columns=["therapy_id", "node", "populations_selected", "status"])
    grouped = (
        assembly.groupby(["therapy_id", "node"], as_index=False)
        .agg(
            populations_selected=("population_id", "nunique"),
            mean_marginal_per_1000=("marginal_events_avoided_per_1000", "mean"),
            mean_step=("step", "mean"),
        )
        .sort_values(["populations_selected", "mean_marginal_per_1000"], ascending=False)
        .reset_index(drop=True)
    )
    grouped["status"] = COMPOSITE_TAG
    return grouped


def composite_report(
    assembly: pd.DataFrame, matrix: pd.DataFrame, criteria: pd.DataFrame, frequency: pd.DataFrame
) -> str:
    lines = [
        "Best-of-breed composite: strongest part per disease node, tested everywhere",
        "=" * 74,
        DISCLAIMER,
        "",
        "HOW IT WAS BUILT",
        "Greedy selection by marginal event benefit, one component per disease node, with",
        "an explicit same-target rule so PCSK9 editing and a PCSK9 antibody can never be",
        "stacked. Components that add nothing are refused automatically. Then a coverage",
        "pass attaches the best part for cure-relevant nodes that are event-neutral",
        "(the myocardium), at zero event credit, so optimising events cannot silently",
        "drop the heart.",
        "",
        "THE RECURRING CORE (across populations and modes)",
    ]
    for row in frequency.head(12).itertuples(index=False):
        lines.append(
            f"  {row.therapy_id:32s} {row.node:20s} "
            f"{row.populations_selected:2d} populations  "
            f"mean marginal {row.mean_marginal_per_1000:6.2f}/1000"
        )
    for mode in CreditMode:
        subset = matrix.loc[matrix["mode"].eq(mode.value)]
        evidence = criteria.loc[
            criteria["variant"].eq(Variant.EVIDENCE_ONLY.value) & criteria["mode"].eq(mode.value)
        ]
        lines.extend(["", f"MODE: {mode.value.upper()}", "-" * 74])
        top = subset.sort_values("events_avoided_per_1000", ascending=False).head(6)
        for row in top.itertuples(index=False):
            lines.append(
                f"  {row.population_id:28s} {row.baseline_risk_pct:5.2f}% -> "
                f"{row.composite_risk_pct:5.2f}%  avoided {row.events_avoided_per_1000:6.2f}/1000"
                f"  excess left {row.excess_risk_remaining_fraction:5.1%}"
            )
        empty = int(subset["components"].eq(0).sum())
        lines.append(
            f"  populations with an empty composite: {empty} of {len(subset)} "
            "(no agent earns credit there)"
        )
        lines.append(f"  max criteria met: {int(evidence['criteria_met'].max())} of 5")
        lines.append("  best score per criterion across populations:")
        for name in CRITERION_ORDER:
            value = float(evidence[name].max())
            flag = "MET" if value >= CRITERION_THRESHOLD else "blocking"
            lines.append(f"    {name:24s} {value:.3f}  {flag}")
    lines.extend(
        [
            "",
            "WHAT CHANGES WITH THE CREDIT MODE",
            "Strict mode cannot include gene editing, because editing has no outcomes",
            "evidence - so durable driver control scores 0. Surrogate and exploratory modes",
            "admit it and driver control reaches 1.0. That difference is the honest measure",
            "of how much of the composite depends on evidence we do not yet have.",
            "",
            "WHAT THE INVENTED NODES ADD",
            "Nothing to any score. The REVERT-1 clearance loop and the CAR-Treg lock are",
            "listed as the nodes the catalogue cannot serve, and carry no evidence credit.",
            "Their hypothetical effect lives in results/hypotheses and results/discovery,",
            "never here.",
            "",
            "THE THREE RULES THAT KEEP THIS HONEST",
            "1. Coverage combines by maximum. Two half-restorations are not one heart.",
            "2. The same molecular target is never stacked, even across exclusive groups.",
            "3. An untrialled combination scores zero on randomised proof. Every component",
            "   can have phase-3 evidence and the combination still has none.",
            "",
            "BOTTOM LINE",
            f"Composites declared a cure: {int(criteria['is_cure'].sum())}.",
            "Risk normalisation is the one criterion the composite can clear, and it does -",
            "but burden reversal (0.60), myocardial restoration (0.50) and randomised proof",
            "(0.00) stay below threshold in every population and every mode. Assembling the",
            "best proven parts produces the strongest supportable regimen, not a cure.",
        ]
    )
    return "\n".join(lines) + "\n"


def run_composite_study(
    output_dir,
    population_ids: tuple[str, ...] | None = None,
    horizon_years: float = 5.0,
    max_components: int = 8,
) -> dict:
    """Write the composite artefacts to their own directory."""
    from pathlib import Path

    root = Path(output_dir)
    (root / "tables").mkdir(parents=True, exist_ok=True)
    (root / "reports").mkdir(parents=True, exist_ok=True)

    catalog, anchors, populations = load_catalog(), load_anchors(), load_populations()
    assembly, matrix, criteria = evaluate_matrix(
        catalog, anchors, populations, population_ids, horizon_years, max_components
    )
    frequency = component_frequency(assembly)

    ablations = []
    days = time_grid(horizon_years, 30.0)
    ids = tuple(sorted(populations)) if population_ids is None else population_ids
    for population_id in ids:
        composite = assemble_composite(
            populations[population_id], catalog, anchors, days, max_components=max_components
        )
        if composite.therapy_ids:
            ablations.append(
                leave_one_out(composite, catalog, anchors, populations[population_id], days)
            )
    ablation = pd.concat(ablations, ignore_index=True) if ablations else pd.DataFrame()

    assembly.to_csv(root / "tables" / "composite_assembly.csv", index=False)
    matrix.to_csv(root / "tables" / "composite_test_matrix.csv", index=False)
    criteria.to_csv(root / "tables" / "composite_criteria.csv", index=False)
    frequency.to_csv(root / "tables" / "composite_component_frequency.csv", index=False)
    ablation.to_csv(root / "tables" / "composite_leave_one_out.csv", index=False)
    pd.DataFrame(
        [
            {
                "node": node,
                "exclusive_groups": ";".join(groups),
                "invented_only": node in INVENTED_NODES,
                "status": COMPOSITE_TAG,
            }
            for node, groups in sorted(DISEASE_NODES.items())
        ]
        + [
            {
                "node": node,
                "exclusive_groups": "",
                "invented_only": True,
                "status": INVENTED_TAG,
            }
            for node in sorted(INVENTED_NODES)
            if node not in DISEASE_NODES
        ]
    ).to_csv(root / "tables" / "composite_node_map.csv", index=False)
    (root / "reports" / "composite_report.txt").write_text(
        composite_report(assembly, matrix, criteria, frequency), encoding="utf-8"
    )
    (root / "DISCLAIMER.txt").write_text(DISCLAIMER + "\n", encoding="utf-8")

    evidence = criteria.loc[criteria["variant"].eq(Variant.EVIDENCE_ONLY.value)]
    return {
        "status": "untrialled_composite",
        "disclaimer": DISCLAIMER,
        "output_dir": str(root),
        "populations_tested": len(ids),
        "credit_modes": len(CreditMode),
        "matrix_rows": len(matrix),
        "max_criteria_met": int(evidence["criteria_met"].max()),
        "cures_declared": int(criteria["is_cure"].sum()),
        "outputs": 8,
    }
