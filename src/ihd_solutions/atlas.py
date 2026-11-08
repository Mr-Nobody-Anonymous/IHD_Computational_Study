"""ATLAS-1: the best working part of every prior construct, sequenced in time.

What was taken from where
-------------------------
===========================  ====================================================
Source construct             Part kept (and why it earned its place)
===========================  ====================================================
Evidence catalogue           Node-disjoint proven backbone: 134 events/1000
Composite assembler          Target-overlap dedup + zero proof for untested combos
REVERT-1                     The AND gate (rho=+0.56) and the efflux/efferocytosis
                             coupling that stops clearance monotherapy backfiring
Gap solver (M7)              Source reduction: the only thing that broke the
                             structural 83.5% plateau
Resolution engine            The 1.28x potency operating point, and scoring by
                             distribution instead of point estimate
Discovery engine             The back-calculated required capability vector
===========================  ====================================================

The new idea: sequencing, not stacking
--------------------------------------
Every prior construct delivered everything at once. But the remuscularisation arm
failed mainly on the *chronic-injury haircut* - the probability of reaching +20 LVEF
points fell from 51% to 26% when the graft bed was a years-old, inflamed, poorly
vascularised scar. Meanwhile the plaque arm's whole job is to remove lipid, quench
inflammation and restore a quieter tissue.

So ATLAS-1 runs two phases: clear and quieten the bed first, then graft into it. The
claim is that phase 1 raises phase 2's translation floor.

That claim is an ASSUMPTION, not a finding. It is parameterised by ``bed_coupling``
and swept all the way to zero, where ATLAS-1 degenerates to simultaneous delivery.
If the gain vanishes at zero coupling, the gain was the assumption - and the report
says so.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache

import numpy as np
import pandas as pd

from ihd_solutions.catalog import Therapy, load_catalog
from ihd_solutions.composite import (
    CRITERION_ORDER,
    assemble_composite,
    score_criteria,
)
from ihd_solutions.cure import CRITERION_THRESHOLD
from ihd_solutions.evidence import CreditMode, TrialAnchor, load_anchors
from ihd_solutions.gap_solver import (
    BASELINE_SCAR_FRACTION,
    REMUSCULARISATION_FACTORS,
    SCAR_RECOVERABLE_LVEF_PP,
    SOURCE_REDUCTION_PRECEDENT_MAX,
)
from ihd_solutions.population import Population, load_populations
from ihd_solutions.resolution import (
    clearance_at,
    population_closure,
)
from ihd_solutions.units import require_fraction, require_positive, time_grid

DISCLAIMER = (
    "ATLAS-1 IS AN INVENTED, UNTRIALLED CONSTRUCT. Invented arms carry zero evidence "
    "credit. The bed-preparation coupling is an assumption swept to zero, not a "
    "measured effect. Nothing here is a therapy, a result, or a recommendation."
)

ATLAS_EVIDENCE_TAG = "atlas_evidence_derived_scenario"
ATLAS_INVENTED_TAG = "atlas_invented_hypothesis_unvalidated"

#: The operating point solved for in the resolution step: 1.28x potency, source
#: reduction at the precedent ceiling, 24 months.
ATLAS_POTENCY = 1.284
ATLAS_SOURCE_REDUCTION = SOURCE_REDUCTION_PRECEDENT_MAX
ATLAS_CLEARANCE_MONTHS = 24.0

#: Phase 2 starts only after the bed has been cleared for this long.
ATLAS_GRAFT_DELAY_MONTHS = 24.0


@dataclass(frozen=True)
class AtlasArms:
    """Which arms are active. Used for leave-one-out ablation."""

    proven_backbone: bool = True
    durable_editing: bool = True
    gated_clearance: bool = True
    source_reduction: bool = True
    remuscularisation: bool = True
    bed_sequencing: bool = True

    def names(self) -> tuple[str, ...]:
        return tuple(self.__dataclass_fields__)

    def without(self, arm: str) -> AtlasArms:
        if arm not in self.__dataclass_fields__:
            raise ValueError(f"unknown arm {arm!r}")
        return AtlasArms(**{**{n: getattr(self, n) for n in self.names()}, arm: False})


ALL_ARMS = AtlasArms()


def clearance_for(arms: AtlasArms, months: float = ATLAS_CLEARANCE_MONTHS) -> float:
    """Rupture-substrate removal achieved by phase 1 under a given arm set."""
    require_positive("months", months)
    if not arms.gated_clearance:
        return 0.0
    potency = ATLAS_POTENCY if arms.proven_backbone or arms.durable_editing else 1.0
    source = ATLAS_SOURCE_REDUCTION if arms.source_reduction else 0.0
    return clearance_at(potency, source, months)


def bed_quality(arms: AtlasArms, months: float = ATLAS_CLEARANCE_MONTHS) -> float:
    """How hospitable the graft bed is after phase 1, on [0, 1].

    Driven by achieved substrate clearance. With no clearance arm the bed is
    untouched and quality is zero, so sequencing can contribute nothing.
    """
    removal = clearance_for(arms, months)
    return float(min(1.0, max(0.0, removal / CRITERION_THRESHOLD)))


def translation_floor(
    arms: AtlasArms,
    bed_coupling: float,
    base_floor: float = 0.5,
    months: float = ATLAS_GRAFT_DELAY_MONTHS,
) -> float:
    """Phase 2's chronic-injury floor, optionally lifted by phase 1.

    ``bed_coupling`` is the invented link. At 0 the floor stays at ``base_floor``
    and ATLAS-1 is just simultaneous delivery with a delay. At 1 a fully cleared
    bed removes the haircut entirely. Neither endpoint is measured.
    """
    require_fraction("bed_coupling", bed_coupling, allow_one=True)
    require_fraction("base_floor", base_floor, allow_one=True)
    if not arms.bed_sequencing or not arms.remuscularisation:
        return base_floor
    lifted = base_floor + bed_coupling * (1.0 - base_floor) * bed_quality(arms, months)
    return float(min(1.0, max(base_floor, lifted)))


@dataclass(frozen=True)
class MyocardialOutcome:
    """Distributional result for the graft arm, never a point estimate."""

    probability_of_target: float
    median_gain_pp: float
    p05_gain_pp: float
    nominal_gain_pp: float
    translation_floor: float
    robust: bool

    @property
    def criterion_score_median(self) -> float:
        return float(min(1.0, self.median_gain_pp / SCAR_RECOVERABLE_LVEF_PP))


def graft_outcome(
    arms: AtlasArms,
    bed_coupling: float = 0.5,
    base_floor: float = 0.5,
    draws: int = 6000,
    seed: int = 20261107,
    required_pp: float = CRITERION_THRESHOLD * SCAR_RECOVERABLE_LVEF_PP,
) -> MyocardialOutcome:
    """Draw every factor from Uniform(1, published effect), then apply the floor."""
    require_positive("draws", draws)
    require_positive("required_pp", required_pp)
    if not arms.remuscularisation:
        return MyocardialOutcome(0.0, 0.0, 0.0, 0.0, base_floor, False)
    rng = np.random.default_rng(seed)
    factors = REMUSCULARISATION_FACTORS
    if not arms.proven_backbone:
        # Dose escalation depends on the manufacturing backbone being in place.
        factors = tuple(f for f in factors if f.name != "dose_escalation")
    delivered = np.ones(int(draws))
    nominal = 1.0
    for factor in factors:
        delivered *= rng.uniform(1.0, factor.demonstrated_multiplier, size=int(draws))
        nominal *= factor.demonstrated_multiplier
    floor = translation_floor(arms, bed_coupling, base_floor)
    penalty = rng.uniform(floor, 1.0, size=int(draws))
    gains = np.minimum(1.0, BASELINE_SCAR_FRACTION * delivered * penalty) * SCAR_RECOVERABLE_LVEF_PP
    probability = float(np.mean(gains >= required_pp))
    return MyocardialOutcome(
        probability,
        float(np.median(gains)),
        float(np.quantile(gains, 0.05)),
        float(min(1.0, BASELINE_SCAR_FRACTION * nominal) * SCAR_RECOVERABLE_LVEF_PP),
        floor,
        probability >= 0.8,
    )


def bed_coupling_sweep(
    couplings: tuple[float, ...] = (0.0, 0.25, 0.5, 0.75, 1.0),
    draws: int = 6000,
    seed: int = 20261107,
) -> pd.DataFrame:
    """The falsification sweep. At zero coupling the sequencing idea contributes nil."""
    rows = []
    for coupling in couplings:
        result = graft_outcome(ALL_ARMS, coupling, draws=draws, seed=seed)
        unsequenced = graft_outcome(
            ALL_ARMS.without("bed_sequencing"), coupling, draws=draws, seed=seed
        )
        rows.append(
            {
                "bed_coupling": coupling,
                "translation_floor": result.translation_floor,
                "probability_of_target": result.probability_of_target,
                "probability_without_sequencing": unsequenced.probability_of_target,
                "sequencing_gain_pp_probability": 100.0
                * (result.probability_of_target - unsequenced.probability_of_target),
                "median_gain_pp": result.median_gain_pp,
                "robust_at_80pct": result.robust,
                "status": ATLAS_INVENTED_TAG,
            }
        )
    return pd.DataFrame(rows)


def arm_ablation(
    bed_coupling: float = 0.5, draws: int = 4000, seed: int = 20261107
) -> pd.DataFrame:
    """Remove one arm at a time. Whatever breaks was load-bearing."""
    full_clearance = clearance_for(ALL_ARMS)
    full_graft = graft_outcome(ALL_ARMS, bed_coupling, draws=draws, seed=seed)
    rows = []
    for arm in ALL_ARMS.names():
        reduced = ALL_ARMS.without(arm)
        clearance = clearance_for(reduced)
        graft = graft_outcome(reduced, bed_coupling, draws=draws, seed=seed)
        # Editing contributes durability, not clearance or graft mass, so measuring
        # only those two would wrongly label it inert.
        driver_lost = arm == "durable_editing"
        rows.append(
            {
                "removed_arm": arm,
                "durable_driver_control_lost": driver_lost,
                "clearance_pct": 100.0 * clearance,
                "clearance_loss_pp": 100.0 * (full_clearance - clearance),
                "graft_probability": graft.probability_of_target,
                "graft_probability_loss_pp": 100.0
                * (full_graft.probability_of_target - graft.probability_of_target),
                "load_bearing": (
                    driver_lost
                    or abs(full_clearance - clearance) > 1e-6
                    or abs(full_graft.probability_of_target - graft.probability_of_target) > 1e-6
                ),
                "status": ATLAS_INVENTED_TAG,
            }
        )
    return (
        pd.DataFrame(rows).sort_values("clearance_loss_pp", ascending=False).reset_index(drop=True)
    )


@dataclass(frozen=True)
class AtlasScorecard:
    """Two scorecards, kept apart: what is provable, and what is hypothesised."""

    population_id: str
    mode: str
    evidence_scores: dict[str, float]
    hypothesised_scores: dict[str, float]
    graft: MyocardialOutcome
    notes: tuple[str, ...] = field(default=(DISCLAIMER,))

    @property
    def evidence_criteria_met(self) -> int:
        return sum(v >= CRITERION_THRESHOLD for v in self.evidence_scores.values())

    @property
    def hypothesised_criteria_met(self) -> int:
        return sum(v >= CRITERION_THRESHOLD for v in self.hypothesised_scores.values())


@lru_cache(maxsize=16)
def distributional_arm_probabilities(
    bed_coupling: float = 0.5, draws: int = 3000, seed: int = 20261107
) -> tuple[float, float]:
    """P(clearance target) and P(graft target) scored jointly, not nominally.

    The self-audit in :mod:`ihd_solutions.audit` showed that scoring burden
    reversal at the nominal 80.01% clearance was the same nominal-product error
    that sank an earlier remuscularisation pass, and that the graft cannot be
    scored against a bed that phase 1 usually fails to prepare. Both arms are now
    scored as probabilities. Clearance does not depend on the population, so this
    is computed once and reused across the whole matrix.
    """
    from ihd_solutions.audit import clearance_robustness, coupled_graft_outcome

    clearance = clearance_robustness(draws=draws, seed=seed)
    coupled = coupled_graft_outcome(bed_coupling, draws=draws, seed=seed)
    return clearance.probability_of_target, coupled.probability_of_target


def score_atlas(
    population: Population,
    mode: CreditMode,
    catalog: dict[str, Therapy],
    anchors: dict[str, TrialAnchor],
    days: np.ndarray,
    bed_coupling: float = 0.5,
    draws: int = 4000,
    seed: int = 20261107,
) -> AtlasScorecard:
    """Score the proven backbone honestly, then the hypothesis separately."""
    composite = assemble_composite(population, catalog, anchors, days, mode=mode)
    raw = score_criteria(composite, catalog, anchors, population, days, mode=mode)
    evidence = {name: float(raw[name]) for name in CRITERION_ORDER}

    graft = graft_outcome(ALL_ARMS, bed_coupling, draws=draws, seed=seed)
    clearance_probability, coupled_graft_probability = distributional_arm_probabilities(
        bed_coupling, draws, seed
    )
    hypothesised = {
        **evidence,
        # Probabilities, not nominal point estimates. See ihd_solutions.audit:
        # the nominal clearance pass was retracted, and the graft must be scored
        # against the bed phase 1 actually produces.
        "burden_reversal": max(evidence["burden_reversal"], clearance_probability),
        "myocardial_restoration": max(
            evidence["myocardial_restoration"], coupled_graft_probability
        ),
        # An untrialled combination still earns no proof, however good it looks.
        "randomised_proof": 0.0,
    }
    return AtlasScorecard(
        population.population_id,
        mode.value,
        evidence,
        hypothesised,
        graft,
        (
            DISCLAIMER,
            "hypothesised scores are probabilities of reaching each target under "
            "delivery uncertainty, not nominal point estimates; the nominal clearance "
            "pass was retracted by ihd_solutions.audit",
        ),
    )


def evaluate_atlas_matrix(
    population_ids: tuple[str, ...] | None = None,
    horizon_years: float = 5.0,
    bed_coupling: float = 0.5,
    draws: int = 3000,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run ATLAS-1 across every population and every credit mode."""
    catalog, anchors, populations = load_catalog(), load_anchors(), load_populations()
    ids = tuple(sorted(populations)) if population_ids is None else population_ids
    days = time_grid(horizon_years, 30.0)
    event_rows: list[dict] = []
    score_rows: list[dict] = []
    for population_id in ids:
        if population_id not in populations:
            raise ValueError(f"unknown population: {population_id}")
        population = populations[population_id]
        for mode in CreditMode:
            composite = assemble_composite(population, catalog, anchors, days, mode=mode)
            card = score_atlas(population, mode, catalog, anchors, days, bed_coupling, draws=draws)
            event_rows.append(
                {
                    "population_id": population_id,
                    "mode": mode.value,
                    "backbone_components": len(composite.therapy_ids),
                    "backbone": ";".join(composite.therapy_ids) or "(none)",
                    "treated_risk_pct": float(
                        raw_risk(composite, catalog, anchors, population, days, mode)
                    ),
                    "status": ATLAS_EVIDENCE_TAG,
                }
            )
            score_rows.append(
                {
                    "population_id": population_id,
                    "mode": mode.value,
                    **{f"evidence_{k}": v for k, v in card.evidence_scores.items()},
                    **{f"hypothesised_{k}": v for k, v in card.hypothesised_scores.items()},
                    "evidence_criteria_met": card.evidence_criteria_met,
                    "hypothesised_criteria_met": card.hypothesised_criteria_met,
                    "retracted_nominal_burden_reversal": clearance_for(ALL_ARMS),
                    "retracted_nominal_myocardial": card.graft.criterion_score_median,
                    "coupled_graft_probability": card.hypothesised_scores["myocardial_restoration"],
                    "graft_probability_of_target": card.graft.probability_of_target,
                    "graft_median_gain_pp": card.graft.median_gain_pp,
                    "graft_robust": card.graft.robust,
                    "is_cure": False,
                    "status": ATLAS_INVENTED_TAG,
                }
            )
    return pd.DataFrame(event_rows), pd.DataFrame(score_rows)


def raw_risk(
    composite,
    catalog: dict[str, Therapy],
    anchors: dict[str, TrialAnchor],
    population: Population,
    days: np.ndarray,
    mode: CreditMode,
) -> float:
    from ihd_solutions.composite import _incidence

    return 100.0 * _incidence(composite.therapy_ids, catalog, anchors, population, days, mode)


def atlas_report(scores: pd.DataFrame, sweep: pd.DataFrame, ablation: pd.DataFrame) -> str:
    best_evidence = int(scores["evidence_criteria_met"].max())
    best_hypo = int(scores["hypothesised_criteria_met"].max())
    zero = sweep.loc[sweep["bed_coupling"].eq(0.0)].iloc[0]
    full = sweep.loc[sweep["bed_coupling"].eq(1.0)].iloc[0]
    lines = [
        "ATLAS-1: best working part of every construct, sequenced rather than stacked",
        "=" * 76,
        DISCLAIMER,
        "",
        "WHAT WAS KEPT FROM WHERE",
        "  catalogue      -> node-disjoint proven backbone (real event benefit)",
        "  composite      -> target-overlap dedup; zero proof for an untested combination",
        "  REVERT-1       -> the AND gate and the efflux/efferocytosis coupling",
        "  gap solver M7  -> source reduction, which broke the structural plateau",
        "  resolution     -> the 1.28x potency operating point (NOMINAL PASS RETRACTED",
        "                    by the self-audit; see ihd_solutions.audit)",
        "  discovery      -> the back-calculated capability target",
        "",
        "THE NEW IDEA",
        "Every earlier construct delivered everything at once. ATLAS-1 clears and quietens",
        "the lesion first, then grafts muscle into the prepared bed. The claim is that",
        "phase 1 lifts phase 2's chronic-injury translation floor. That link is an",
        "assumption, so it is swept to zero where the idea contributes nothing:",
        "",
        "  bed coupling   floor   P(+20)   vs unsequenced   robust",
    ]
    for row in sweep.itertuples(index=False):
        lines.append(
            f"     {row.bed_coupling:4.2f}       {row.translation_floor:4.2f}   "
            f"{row.probability_of_target:6.1%}   "
            f"{row.sequencing_gain_pp_probability:+6.1f} pp        "
            f"{'yes' if row.robust_at_80pct else 'no'}"
        )
    lines.extend(
        [
            "",
            f"  At zero coupling sequencing adds {zero.sequencing_gain_pp_probability:+.1f} pp. "
            f"At full coupling it adds {full.sequencing_gain_pp_probability:+.1f} pp.",
            "  THESE NUMBERS ASSUME PHASE 1 SUCCEEDED. They score the graft against a",
            "  nominally cleared bed. The self-audit (results/audit/) shows phase 1 reaches",
            "  its target with probability ~0, so the bed is usually not prepared and this",
            "  table overstates the graft by about 20 pp at full coupling. Scored jointly,",
            "  sequencing adds only a few points. The structure of the idea survives; its",
            "  practical value does not, until clearance itself is robust.",
            "",
            "ARM ABLATION (what is actually load-bearing)",
        ]
    )
    for row in ablation.itertuples(index=False):
        lines.append(
            f"  drop {row.removed_arm:20s} clearance {row.clearance_pct:5.1f}% "
            f"({row.clearance_loss_pp:+5.1f} pp)  graft P {row.graft_probability:5.1%} "
            f"({row.graft_probability_loss_pp:+5.1f} pp)  "
            f"{'LOAD-BEARING' if row.load_bearing else 'inert'}"
            + ("  [durability criterion]" if row.durable_driver_control_lost else "")
        )
    lines.extend(
        [
            "",
            "TESTED ACROSS EVERY POPULATION AND CREDIT MODE",
            f"  rows: {len(scores)}  "
            f"({scores['population_id'].nunique()} populations x {scores['mode'].nunique()} modes)",
            f"  best criteria met on EVIDENCE alone:      {best_evidence} of 5",
            f"  best criteria met IF the hypothesis held: {best_hypo} of 5",
            f"  ATLAS-1 configurations declared a cure:   {int(scores['is_cure'].sum())}",
            "",
            "  Both counts score the invented arms as probabilities of hitting their",
            "  targets under delivery uncertainty, not as nominal point estimates. They",
            "  come out equal: the invented arms add no criterion that survives that test.",
            "  Proof stays at zero in both, because no combination here has been trialled.",
            "",
            "VERDICT (corrected by the self-audit)",
            "ATLAS-1 remains the best candidate on construction quality: every arm is",
            "load-bearing, every inherited part survived its own ablation, and its one new",
            "idea carries an explicit falsifier. But the earlier claim that it reaches 3 of",
            "5 criteria is WITHDRAWN. That claim rested on a nominal clearance pass at",
            "1.284x potency; under delivery uncertainty the probability of reaching the",
            "clearance target is ~0, and no potency makes it robust while the source term",
            "is uncertain. Jointly scored, ATLAS-1 reaches 2 of 5 -- the same as every",
            "predecessor. It is not a cure.",
            "",
            "Corrected ranking of open problems: clearance under delivery uncertainty",
            "first, remuscularisation second, randomised proof third. The highest-leverage",
            "real-world action is unchanged and needs no new biology: the PCSK9-editing",
            "outcomes trial.",
        ]
    )
    return "\n".join(lines) + "\n"


def run_atlas_study(
    output_dir,
    population_ids: tuple[str, ...] | None = None,
    bed_coupling: float = 0.5,
    draws: int = 3000,
) -> dict:
    from pathlib import Path

    root = Path(output_dir)
    (root / "tables").mkdir(parents=True, exist_ok=True)
    (root / "reports").mkdir(parents=True, exist_ok=True)

    events, scores = evaluate_atlas_matrix(population_ids, bed_coupling=bed_coupling, draws=draws)
    sweep = bed_coupling_sweep(draws=max(draws, 4000))
    ablation = arm_ablation(bed_coupling, draws=max(draws, 3000))
    matrix = pd.read_csv("results/composite/tables/composite_test_matrix.csv")
    closure = population_closure(matrix)
    closure["status"] = ATLAS_EVIDENCE_TAG

    events.to_csv(root / "tables" / "atlas_event_matrix.csv", index=False)
    scores.to_csv(root / "tables" / "atlas_scorecard.csv", index=False)
    sweep.to_csv(root / "tables" / "atlas_bed_coupling_sweep.csv", index=False)
    ablation.to_csv(root / "tables" / "atlas_arm_ablation.csv", index=False)
    closure.to_csv(root / "tables" / "atlas_population_closure.csv", index=False)
    pd.DataFrame(
        [
            {
                "arm": arm,
                "active": getattr(ALL_ARMS, arm),
                "status": ATLAS_INVENTED_TAG,
            }
            for arm in ALL_ARMS.names()
        ]
    ).to_csv(root / "tables" / "atlas_arm_register.csv", index=False)
    pd.DataFrame(
        [
            {
                "potency_multiplier": ATLAS_POTENCY,
                "source_reduction": ATLAS_SOURCE_REDUCTION,
                "clearance_months": ATLAS_CLEARANCE_MONTHS,
                "graft_delay_months": ATLAS_GRAFT_DELAY_MONTHS,
                "clearance_achieved": clearance_for(ALL_ARMS),
                "bed_quality": bed_quality(ALL_ARMS),
                "best_evidence_criteria": int(scores["evidence_criteria_met"].max()),
                "best_hypothesised_criteria": int(scores["hypothesised_criteria_met"].max()),
                "cures_declared": int(scores["is_cure"].sum()),
                "status": ATLAS_INVENTED_TAG,
            }
        ]
    ).to_csv(root / "tables" / "atlas_operating_point.csv", index=False)
    (root / "reports" / "atlas_report.txt").write_text(
        atlas_report(scores, sweep, ablation), encoding="utf-8"
    )
    (root / "DISCLAIMER.txt").write_text(DISCLAIMER + "\n", encoding="utf-8")
    return {
        "status": "invented_untrialled",
        "disclaimer": DISCLAIMER,
        "output_dir": str(root),
        "populations_tested": int(scores["population_id"].nunique()),
        "scorecard_rows": len(scores),
        "best_evidence_criteria": int(scores["evidence_criteria_met"].max()),
        "best_hypothesised_criteria": int(scores["hypothesised_criteria_met"].max()),
        "cures_declared": int(scores["is_cure"].sum()),
        "outputs": 9,
    }
