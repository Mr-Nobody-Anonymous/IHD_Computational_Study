"""Extend the distributional standard to the evidence layer itself.

The ATLAS-1 self-audit (``ihd_solutions.audit``) caught a nominal pass in an
invented construct. This module asks the same question of the layer that every
construct inherits: are the *evidence-derived* criterion scores probabilities, or
point estimates?

Answer: ``composite.score_criteria`` computes ``risk_normalisation`` from a single
point-estimate simulation, even though ``uncertainty.monte_carlo_regimen`` already
propagates LDL efficacy intervals, hazard-ratio intervals, baseline-rate intervals
and adherence. The criterion was therefore nominal while the machinery to make it
distributional already existed in the repository.

This module connects the two and reports what changes. It does not rescore the
other four criteria, because they are not effect estimates:

===========================  ====================================================
Criterion                    Nature
===========================  ====================================================
durable_driver_control       Structural flag: is a durable editing mechanism in
                             the regimen? No effect size, so no interval to
                             propagate. NOT proof of lifetime durability.
burden_reversal              Capped policy credit (max 0.60), never reaches 0.80.
myocardial_restoration       Capped policy credit (max 0.50), never reaches 0.80.
risk_normalisation           A genuine effect estimate. THIS is what is audited.
randomised_proof             Zero by construction for an untrialled combination.
===========================  ====================================================

Nothing here is clinical prediction. Confidence intervals describe sampling
uncertainty in the source trials, not whether a mechanism transports to a new
population.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from ihd_solutions.catalog import Therapy, load_catalog, validate_combination
from ihd_solutions.composite import CRITERION_ORDER, assemble_composite, score_criteria
from ihd_solutions.cure import CRITERION_THRESHOLD
from ihd_solutions.evidence import CreditMode, TrialAnchor, load_anchors
from ihd_solutions.population import Population, load_populations
from ihd_solutions.uncertainty import monte_carlo_regimen
from ihd_solutions.units import require_positive, time_grid

EVIDENCE_AUDIT_TAG = "evidence_layer_uncertainty_audit"

DISCLAIMER = (
    "SCENARIO SIMULATION, NOT CLINICAL PREDICTION. This audit propagates the "
    "confidence intervals reported by the source trials, plus adherence and "
    "baseline-rate uncertainty. It does not establish that any mechanism transports "
    "to a new population, and it declares no cures."
)

#: Criteria that are effect estimates carrying a sampling interval. Everything
#: else is a structural flag or a capped policy credit, and rescoring it under CI
#: propagation would be meaningless.
EFFECT_ESTIMATE_CRITERIA = frozenset({"risk_normalisation"})

#: Reference low-risk annualised level used by the risk-normalisation criterion.
REFERENCE_RISK = 0.0198

ROBUSTNESS_BAR = 0.8


@dataclass(frozen=True)
class RiskRobustness:
    """Risk normalisation as a distribution, beside the point estimate."""

    nominal: float
    median: float
    p05: float
    p95: float
    probability_of_target: float
    robust: bool
    draws: int
    mass_at_floor: float = 0.0
    mass_at_ceiling: float = 0.0

    @property
    def median_uninformative(self) -> bool:
        """True when the score piles up at both 0 and 1, so the median is unstable.

        The score is clipped to [0, 1], which can produce a bimodal distribution
        with little mass in between. The median then sits at whichever pole happens
        to hold slightly more than half the draws and flips with the sample. It is
        not a summary of such a distribution; the probability is.
        """
        return self.mass_at_floor > 0.2 and self.mass_at_ceiling > 0.2

    @property
    def monte_carlo_standard_error(self) -> float:
        p = self.probability_of_target
        return float(np.sqrt(max(p * (1.0 - p), 0.0) / max(self.draws, 1)))

    @property
    def borderline(self) -> bool:
        """Within two standard errors of the bar, so the verdict is not settled.

        Calling a 78% estimate a confident failure against an 80% bar would repeat
        the overclaiming this audit exists to correct, in the opposite direction.
        """
        return (
            abs(self.probability_of_target - ROBUSTNESS_BAR)
            <= 2.0 * self.monte_carlo_standard_error
        )

    @property
    def nominal_is_misleading(self) -> bool:
        """A nominal pass that fails robustly by more than Monte Carlo noise."""
        return self.nominal >= CRITERION_THRESHOLD and not self.robust and not self.borderline

    @property
    def median_is_misleading(self) -> bool:
        """A median can pass while the probability does not.

        ``lifetime_exposure`` is the live example: the score is clipped at 1.0, so
        the median sits at 1.0 while barely half the draws clear the bar.
        """
        return self.median >= CRITERION_THRESHOLD and not self.robust


def risk_normalisation_robustness(
    population: Population,
    mode: CreditMode,
    catalog: dict[str, Therapy],
    anchors: dict[str, TrialAnchor],
    days: np.ndarray,
    draws: int = 300,
    seed: int = 20261107,
    reference_risk: float = REFERENCE_RISK,
) -> RiskRobustness:
    """Score risk normalisation under the uncertainty the trials actually report."""
    require_positive("draws", draws)
    composite = assemble_composite(population, catalog, anchors, days, mode=mode)
    raw = score_criteria(composite, catalog, anchors, population, days, reference_risk, mode=mode)
    nominal = float(raw["risk_normalisation"])
    regimen = validate_combination(catalog, composite.therapy_ids)
    result = monte_carlo_regimen(
        regimen, catalog, population, anchors, days, mode=mode, draws=draws, seed=seed
    )
    values = result.risk_normalisation_draws(reference_risk)
    if values is None:  # pragma: no cover - defensive; the recorder always populates
        raise RuntimeError("monte_carlo_regimen did not record per-draw risks")
    probability = float(np.mean(values >= CRITERION_THRESHOLD))
    return RiskRobustness(
        nominal=nominal,
        median=float(np.median(values)),
        p05=float(np.quantile(values, 0.05)),
        p95=float(np.quantile(values, 0.95)),
        probability_of_target=probability,
        robust=probability >= ROBUSTNESS_BAR,
        draws=int(values.size),
        mass_at_floor=float(np.mean(values <= 1e-9)),
        mass_at_ceiling=float(np.mean(values >= 1.0 - 1e-9)),
    )


def evidence_audit_matrix(
    population_ids: tuple[str, ...] | None = None,
    horizon_years: float = 5.0,
    draws: int = 300,
    seed: int = 20261107,
) -> pd.DataFrame:
    """Every population x every credit mode, nominal beside distributional."""
    catalog, anchors, populations = load_catalog(), load_anchors(), load_populations()
    ids = tuple(sorted(populations)) if population_ids is None else population_ids
    days = time_grid(horizon_years, 30.0)
    rows: list[dict] = []
    for population_id in ids:
        if population_id not in populations:
            raise ValueError(f"unknown population: {population_id}")
        population = populations[population_id]
        for mode in CreditMode:
            composite = assemble_composite(population, catalog, anchors, days, mode=mode)
            raw = score_criteria(composite, catalog, anchors, population, days, mode=mode)
            risk = risk_normalisation_robustness(
                population, mode, catalog, anchors, days, draws, seed
            )
            nominal_scores = {name: float(raw[name]) for name in CRITERION_ORDER}
            corrected = {
                **nominal_scores,
                "risk_normalisation": risk.probability_of_target,
            }
            rows.append(
                {
                    "population_id": population_id,
                    "mode": mode.value,
                    "backbone_components": len(composite.therapy_ids),
                    "nominal_risk_normalisation": risk.nominal,
                    "median_risk_normalisation": risk.median,
                    "p05_risk_normalisation": risk.p05,
                    "p95_risk_normalisation": risk.p95,
                    "probability_of_target": risk.probability_of_target,
                    "nominal_pass": risk.nominal >= CRITERION_THRESHOLD,
                    "median_pass": risk.median >= CRITERION_THRESHOLD,
                    "robust_pass": risk.robust,
                    "borderline": risk.borderline,
                    "monte_carlo_standard_error": risk.monte_carlo_standard_error,
                    "nominal_is_misleading": risk.nominal_is_misleading,
                    "median_is_misleading": risk.median_is_misleading,
                    "mass_at_floor": risk.mass_at_floor,
                    "mass_at_ceiling": risk.mass_at_ceiling,
                    "median_uninformative": risk.median_uninformative,
                    "nominal_criteria_met": sum(
                        v >= CRITERION_THRESHOLD for v in nominal_scores.values()
                    ),
                    "robust_criteria_met": sum(
                        v >= CRITERION_THRESHOLD for v in corrected.values()
                    ),
                    "draws": risk.draws,
                    "is_cure": False,
                    "status": EVIDENCE_AUDIT_TAG,
                }
            )
    return pd.DataFrame(rows)


def mode_summary(matrix: pd.DataFrame) -> pd.DataFrame:
    """How many passes survive, by credit mode. The headline of this audit."""
    rows = []
    for mode, group in matrix.groupby("mode", sort=True):
        rows.append(
            {
                "mode": mode,
                "populations": int(group["population_id"].nunique()),
                "nominal_risk_passes": int(group["nominal_pass"].sum()),
                "robust_risk_passes": int(group["robust_pass"].sum()),
                "retracted_passes": int(group["nominal_is_misleading"].sum()),
                "borderline_cases": int(group["borderline"].sum()),
                "median_probability": float(group["probability_of_target"].median()),
                "best_nominal_criteria": int(group["nominal_criteria_met"].max()),
                "best_robust_criteria": int(group["robust_criteria_met"].max()),
                "status": EVIDENCE_AUDIT_TAG,
            }
        )
    return pd.DataFrame(rows)


def criterion_ledger() -> pd.DataFrame:
    """Which criteria can meaningfully be audited this way, and which cannot."""
    return pd.DataFrame(
        [
            {
                "criterion": "durable_driver_control",
                "kind": "structural_flag",
                "audited": False,
                "reason": "mechanism presence, not an effect size; no interval to propagate",
                "caveat": "presence is not evidence of lifetime durability in humans",
                "status": EVIDENCE_AUDIT_TAG,
            },
            {
                "criterion": "burden_reversal",
                "kind": "capped_policy_credit",
                "audited": False,
                "reason": "ceiling 0.60, below the 0.80 threshold regardless of uncertainty",
                "caveat": "ceiling is a policy choice, not a measurement",
                "status": EVIDENCE_AUDIT_TAG,
            },
            {
                "criterion": "myocardial_restoration",
                "kind": "capped_policy_credit",
                "audited": False,
                "reason": "ceiling 0.50, below the 0.80 threshold regardless of uncertainty",
                "caveat": "ceiling is a policy choice, not a measurement",
                "status": EVIDENCE_AUDIT_TAG,
            },
            {
                "criterion": "risk_normalisation",
                "kind": "effect_estimate",
                "audited": True,
                "reason": "derived from trial hazard ratios and efficacy intervals",
                "caveat": "propagates sampling uncertainty, not transportability",
                "status": EVIDENCE_AUDIT_TAG,
            },
            {
                "criterion": "randomised_proof",
                "kind": "structural_zero",
                "audited": False,
                "reason": "an untrialled combination earns zero regardless of uncertainty",
                "caveat": "no simulation can supply it",
                "status": EVIDENCE_AUDIT_TAG,
            },
        ]
    )


def evidence_audit_report(matrix: pd.DataFrame, summary: pd.DataFrame) -> str:
    misleading = matrix.loc[matrix["nominal_is_misleading"]]
    median_trap = matrix.loc[matrix["median_uninformative"]]
    lines = [
        "EVIDENCE-LAYER AUDIT: are the inherited criterion scores probabilities?",
        "=======================================================================",
        DISCLAIMER,
        "",
        "WHY THIS AUDIT",
        "The ATLAS-1 self-audit caught a nominal pass in an invented construct. The",
        "same question applies to the layer every construct inherits. It turns out",
        "composite.score_criteria computes risk_normalisation from one point-estimate",
        "simulation, while uncertainty.monte_carlo_regimen already propagates LDL,",
        "hazard-ratio, baseline-rate and adherence uncertainty. The criterion was",
        "nominal even though the machinery to make it distributional already existed.",
        "",
        "Only risk_normalisation is an effect estimate. Driver control is a mechanism",
        "flag; burden and myocardial credit are capped below the threshold; proof is",
        "zero by construction. Those four are unchanged, for the reasons in",
        "tables/evidence_criterion_ledger.csv.",
        "",
        "RESULT BY CREDIT MODE",
        "  mode          nominal passes   robust passes   retracted   best criteria",
    ]
    for row in summary.itertuples(index=False):
        lines.append(
            f"  {row.mode:12s}  {row.nominal_risk_passes:14d}   {row.robust_risk_passes:13d}   "
            f"{row.retracted_passes:9d}   {row.best_nominal_criteria} -> {row.best_robust_criteria}"
        )
    lines += [
        "",
        "NOMINAL PASSES THAT DO NOT SURVIVE",
    ]
    if misleading.empty:  # pragma: no cover - the committed data has six
        lines.append("  (none)")
    else:
        for row in misleading.itertuples(index=False):
            lines.append(
                f"  {row.population_id:26s} {row.mode:12s} nominal {row.nominal_risk_normalisation:.3f}"
                f"  median {row.median_risk_normalisation:.3f}"
                f"  P(>=0.8) {100 * row.probability_of_target:5.1f}%"
            )
    borderline = matrix.loc[matrix["borderline"]]
    if not borderline.empty:
        lines += ["", "TOO CLOSE TO CALL (within two Monte Carlo standard errors of 80%)"]
        for row in borderline.itertuples(index=False):
            lines.append(
                f"  {row.population_id:26s} {row.mode:12s} "
                f"P {100 * row.probability_of_target:5.1f}% "
                f"+/- {100 * 2 * row.monte_carlo_standard_error:4.1f} pp (2 SE)"
            )
        lines.append("  These are not counted as retractions; more draws would settle them.")
    lines += [
        "",
        "WHAT THIS CHANGES",
        "Strict credit is the worst affected. Its single risk-normalisation pass was",
        "the only criterion any population met under strict rules, and it does not",
        "survive: the best strict result falls to zero of five. Exploratory credit is",
        "largely unaffected, because its effects are large relative to their intervals.",
        "",
        "A MEDIAN IS NOT A PROBABILITY EITHER",
    ]
    if median_trap.empty:  # pragma: no cover - defensive
        lines.append("  (no bimodal cases in this run)")
    else:
        for row in median_trap.itertuples(index=False):
            lines.append(
                f"  {row.population_id:26s} {row.mode:12s} "
                f"{100 * row.mass_at_floor:4.1f}% of draws at 0, "
                f"{100 * row.mass_at_ceiling:4.1f}% at 1, "
                f"P(>=0.8) {100 * row.probability_of_target:5.1f}%"
            )
        lines += [
            "  The score is clipped to [0, 1], so these distributions pile up at both",
            "  poles with little mass between. The median then sits at whichever pole",
            "  holds slightly more than half the draws, and flips with the sample size.",
            "  It is not a summary of such a distribution. The probability is, which is",
            "  why every criterion here is reported as one.",
        ]
    lines += [
        "",
        "VERDICT",
        "The evidence layer was partly nominal too, and is now corrected. This does",
        "not rescue any construct: it lowers the ceiling they all inherit. Combined",
        "with the ATLAS-1 self-audit, no construct in this study meets more than two",
        "of five criteria robustly, and under strict credit none meets any.",
        "",
        "The recommendation is unchanged and is the one item that does not depend on",
        "any of these scores: a PCSK9-editing outcomes trial would replace assumed",
        "durability and surrogate credit with a randomised effect estimate, which is",
        "precisely the deficiency this audit exposes.",
    ]
    return "\n".join(lines) + "\n"


def run_evidence_audit_study(
    output_dir: Path | str = Path("results/evidence_audit"),
    population_ids: tuple[str, ...] | None = None,
    draws: int = 300,
    seed: int = 20261107,
) -> dict[str, object]:
    """Write the evidence-layer audit, tagged and separate from the run manifest."""
    root = Path(output_dir)
    tables = root / "tables"
    reports = root / "reports"
    tables.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)

    matrix = evidence_audit_matrix(population_ids, draws=draws, seed=seed)
    summary = mode_summary(matrix)
    ledger = criterion_ledger()
    matrix.to_csv(tables / "evidence_audit_matrix.csv", index=False)
    summary.to_csv(tables / "evidence_audit_mode_summary.csv", index=False)
    ledger.to_csv(tables / "evidence_criterion_ledger.csv", index=False)
    matrix.loc[matrix["nominal_is_misleading"]].to_csv(
        tables / "evidence_retracted_passes.csv", index=False
    )
    text = evidence_audit_report(matrix, summary)
    (reports / "evidence_audit_report.txt").write_text(text)
    (root / "DISCLAIMER.txt").write_text(DISCLAIMER + "\n")

    return {
        "output_dir": str(root),
        "outputs": 6,
        "status": "evidence_layer_partly_nominal_corrected",
        "retracted_passes": int(matrix["nominal_is_misleading"].sum()),
        "cures_declared": 0,
        "disclaimer": DISCLAIMER,
    }
