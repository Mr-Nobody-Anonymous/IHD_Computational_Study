"""Evidence-aware test of complementary IHD medicines (not a covalent hybrid).

A small molecule cannot safely inherit a useful fragment from another drug by simple
copy/paste: changing a scaffold can destroy binding, exposure, selectivity, and safety.
This module therefore tests the scientifically interpretable version of the request:
add an independently acting, outcome-studied medicine to established LDL-lowering, and
separately check whether a new oral PCSK9 agent earns hard-outcome credit.

All combinations below are untrialled scenarios. Their component evidence does not
transfer to the combination.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from ihd_solutions.catalog import load_catalog, validate_combination
from ihd_solutions.evidence import CreditMode, load_anchors
from ihd_solutions.pipeline import PipelineConfig
from ihd_solutions.population import load_populations
from ihd_solutions.risk import simulate_events
from ihd_solutions.units import time_grid

STATUS = "untrialled_complementary_regimen_simulation"
DEFAULT_POPULATION = "secondary_post_mi"

# LDL-lowering is the outcome-supported backbone. Colchicine adds a complementary
# anti-inflammatory mechanism; enlicitide is included separately to distinguish its
# LDL surrogate effect from hard-outcome credit not yet provided by its own outcomes trial.
REGIMENS: dict[str, tuple[str, ...]] = {
    "ldl_backbone": ("aspirin_secondary_prevention", "statin_high_intensity"),
    "ldl_plus_colchicine": (
        "aspirin_secondary_prevention",
        "statin_high_intensity",
        "colchicine",
    ),
    "ldl_plus_enlicitide": (
        "aspirin_secondary_prevention",
        "statin_high_intensity",
        "enlicitide_oral",
    ),
    "ldl_plus_colchicine_plus_enlicitide": (
        "aspirin_secondary_prevention",
        "statin_high_intensity",
        "colchicine",
        "enlicitide_oral",
    ),
}

DISCLAIMER = (
    "SIMULATION ONLY. These are separate medicines in an untrialled regimen, not a new "
    "chemically fused molecule. Combination efficacy and safety are unproven; component "
    "evidence is not randomized proof for the combination. Results are population- and "
    "evidence-credit-mode dependent, and must not guide individual treatment."
)


def complementarity_table(population_id: str = DEFAULT_POPULATION) -> pd.DataFrame:
    """Compare a statin backbone with colchicine and oral PCSK9 additions.

    Strict mode credits only event-outcome evidence supported by the selected population;
    surrogate mode additionally credits surrogate-only evidence. The delta column is a
    difference between model scenarios, not a clinical estimate of combination synergy.
    """
    catalog = load_catalog()
    anchors = load_anchors()
    populations = load_populations()
    if population_id not in populations:
        raise ValueError(f"unknown population: {population_id}")
    population = populations[population_id]
    config = PipelineConfig(make_figures=False)
    days = time_grid(config.horizon_years, config.step_days)

    rows: list[dict[str, object]] = []
    for label, therapy_ids in REGIMENS.items():
        regimen = validate_combination(catalog, therapy_ids)
        for mode in (CreditMode.STRICT, CreditMode.SURROGATE):
            result = simulate_events(regimen, catalog, population, anchors, days, mode=mode)
            summary = result.summary_row()
            rows.append(
                {
                    "scenario": label,
                    "therapy_ids": ";".join(therapy_ids),
                    "population_id": population_id,
                    "credit_mode": mode.value,
                    "events_avoided_per_1000": result.events_avoided_per_1000,
                    "relative_risk_reduction_pct": 100.0 * result.relative_risk_reduction,
                    "number_needed_to_treat": result.number_needed_to_treat,
                    "time_averaged_ldl_mg_dl": result.biomarkers.time_averaged_ldl_mg_dl,
                    "credited_therapies": summary["credited_therapies"],
                    "uncredited_therapies": summary["uncredited_therapies"],
                    "evidence_warnings": summary["evidence_warnings"],
                    "combination_randomized_proof": False,
                    "status": STATUS,
                }
            )

    frame = pd.DataFrame(rows)
    base = frame.loc[frame["scenario"].eq("ldl_backbone")].set_index("credit_mode")
    frame["increment_vs_ldl_backbone_per_1000"] = frame.apply(
        lambda row: float(
            row["events_avoided_per_1000"] - base.loc[row["credit_mode"], "events_avoided_per_1000"]
        ),
        axis=1,
    )
    return frame


def complementarity_report(frame: pd.DataFrame) -> str:
    """Render the tested scenarios without implying that a hybrid molecule exists."""
    lines = [
        "COMPLEMENTARY IHD MEDICINES — MODEL TEST, NOT A NEW MOLECULE",
        "=" * 76,
        DISCLAIMER,
        "",
        "Scenario results by evidence-credit mode (events avoided per 1,000 over the",
        "configured horizon; deltas are versus the LDL-lowering backbone):",
        "",
        f"{'scenario':38s} {'mode':12s} {'avoided/1000':>12s} {'delta':>9s} {'NNT':>9s}  warnings",
        "-" * 118,
    ]
    for row in frame.itertuples(index=False):
        nnt = (
            "not estimable"
            if not pd.notna(row.number_needed_to_treat)
            else f"{row.number_needed_to_treat:.2f}"
        )
        warnings = row.evidence_warnings or "none"
        lines.append(
            f"{row.scenario:38s} {row.credit_mode:12s} "
            f"{row.events_avoided_per_1000:12.2f} {row.increment_vs_ldl_backbone_per_1000:9.2f} "
            f"{nnt:>9s}  {warnings}"
        )
    lines += [
        "",
        "Interpretation:",
        "- Colchicine is the complementary anti-inflammatory add-on tested here; it is",
        "  a separate drug, not a fragment grafted into the statin. Trial results vary by",
        "  context (CLEAR SYNERGY was neutral in ACS), consistent with the transport warning.",
        "- Enlicitide is credited for LDL surrogate change in surrogate mode, but its",
        "  own hard-outcome benefit is not credited in strict mode without outcomes data.",
        "- No row is a tested fixed-dose combination; interaction, tolerability, and",
        "  combination-specific efficacy remain unknown.",
    ]
    return "\n".join(lines)


def write_complementarity_report(output_dir: str | Path) -> pd.DataFrame:
    """Write reproducible CSV and text outputs for the complementary-regimen test."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    frame = complementarity_table()
    frame.to_csv(output / "complementarity_scenarios.csv", index=False)
    (output / "complementarity_report.txt").write_text(complementarity_report(frame) + "\n")
    (output / "DISCLAIMER.txt").write_text(DISCLAIMER + "\n")
    return frame
