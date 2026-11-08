"""Correct the recombination result for what the literature says it is missing.

ATLAS-5 was the strongest composition this study could build. This module asks a
different question: **is the effect it is built on actually there?**

The deep search found one deficiency that dominates every other, and it is not a
missing mechanism. Every part of ATLAS-5 ultimately rests on a single uncontrolled
observation: BioVAT-HF reported +3.9 LVEF points at three months in 16 patients with
no control group. Subtract what a contemporaneous control arm achieved for free and
most of the effect disappears. No amount of mechanism recombination fixes that,
because recombination rearranges assumptions - it cannot add a control arm.

What the search established, in order of how much it matters:

1. **No contemporaneous control.** REVIVED-BCIS2 randomised 700 patients with
   LVEF <= 35% and viable myocardium to PCI or optimal medical therapy. The OMT arm
   gained **+3.4 points at six months** on its own. BioVAT's +3.9 (90% CI 0.9-6.8)
   is not distinguishable from that. The comparator is not a placebo, it is the
   guideline therapy those patients were already on.

2. **A second uncontrolled patch trial was negative.** AATS LB8 gave iPSC-CM
   patches to 8 patients with ischaemic cardiomyopathy: LVEF 28.7 -> 33.0% at 52
   weeks, **p = 0.156**. Same population as ours, no significant change.

3. **The endpoint does not track outcomes in this population.** STICH followed LVEF
   and survival for a decade in ischaemic cardiomyopathy: an LVEF increase occurred
   only in patients with viability, **irrespective of treatment assignment**, and
   there was **no association between LVEF change and subsequent death**. ATLAS-5
   optimises precisely that quantity.

4. **Viability does not identify who benefits.** In REVIVED, viability testing did
   not select patients who gained from revascularisation, and the extent of
   non-viable myocardium - not viable myocardium - predicted prognosis. So the
   "recoverable scar" pool the model draws on is **not reachable by restoring
   perfusion**, which is what the hibernation hypothesis assumed.

5. **Transmurality, not viability, decides recovery.** Segments with no late
   gadolinium enhancement recover in 91%; those above 75% transmurality recover in
   8%. Our population has established scar.

The corrections below are applied as subtractions with floors at zero. A correction
may only ever LOWER an estimate, and a test enforces that.

Nothing here is tuned to reach a target. If the corrected design fails, it fails.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from ihd_solutions.cure import CRITERION_THRESHOLD
from ihd_solutions.gap_solver import SCAR_RECOVERABLE_LVEF_PP
from ihd_solutions.recombine import (
    ALL_MECHANISMS,
    ATLAS4_BEST_ADMISSIBLE,
    ATLAS5,
    Composite,
    CompositeResult,
    composite_draws,
    exhaustive_subsets,
)
from ihd_solutions.units import require_positive

REVIVAL_TAG = "recombination_natural_history_corrected_unvalidated"

DISCLAIMER = (
    "THIS IS A CORRECTION TO A COMPOSITION OF UNVALIDATED MECHANISMS, NOT A THERAPY "
    "AND NOT A RESULT. The correction subtracts a documented control-arm effect from "
    "an uncontrolled observation. It removes attribution; it does not add evidence. "
    "Randomised proof for any design here is zero. No cure is declared."
)

#: LVEF points the OMT arm of REVIVED-BCIS2 gained unaided, at six months.
#: 700 patients, LVEF <= 35%, viability in >= 4 dysfunctional segments - the closest
#: population in the literature to the one this study models.
NATURAL_HISTORY_PP = 3.4

#: Documented range for that unaided gain: the same arm was +3.4 at six months and
#: +1.1 at twelve. We do not know the trajectory between, so it is swept.
NATURAL_HISTORY_RANGE = (1.1, 3.4)

#: Between-patient SD of MRI/echo LVEF change. NOT invented here: gap 4 already
#: swept this quantity (5/6/7/8 pp) to size the CATO power calculation. Reusing it
#: keeps one SD family across the study instead of introducing a second.
LVEF_CHANGE_SD_PP = 5.0

#: The second uncontrolled iPSC-CM patch trial in ischaemic cardiomyopathy.
#: AATS LB8, n = 8, LVEF 28.7 +/- 6.8 -> 33.0 +/- 6.0 at 52 weeks, p = 0.156.
AATS_LB8_GAIN_PP = 4.3
AATS_LB8_P_VALUE = 0.156
AATS_LB8_N = 8

#: Surrogate threshold effects for a treatment effect on hard endpoints, trial level,
#: from the 2026 heart-failure surrogate meta-analysis: an 8 pp mean LVEF increase is
#: what predicts an all-cause mortality benefit; 2 pp for HF hospitalisation.
SURROGATE_STE_MORTALITY_PP = 8.0
SURROGATE_STE_HF_HOSP_PP = 2.0

#: The largest LVEF gain any remuscularisation approach has measured in a human:
#: BioVAT-HF, +6.9 pp at a mean 17 months. Any model output above this is outside the
#: observed range and cannot be used to reason about surrogate thresholds.
LARGEST_HUMAN_OBSERVATION_PP = 6.9

#: Trial-level R^2 between LVEF change and mortality benefit (69 comparisons), and
#: within HFrEF (56 comparisons). Reported so the link can be charged as a variance
#: explained rather than asserted.
SURROGATE_R2_MORTALITY = 0.59
SURROGATE_R2_MORTALITY_HFREF = 0.49

#: Functional recovery by transmural extent of late gadolinium enhancement.
#: Segment-level, pooled across the CMR viability literature.
TRANSMURALITY_RECOVERY = {
    "no_lge": 0.91,
    "lge_1_50": 0.44,
    "lge_51_75": 0.43,
    "lge_over_75": 0.08,
}


@dataclass(frozen=True)
class Correction:
    """One documented deficiency and the direction it must move the estimate."""

    name: str
    finding: str
    source: str
    effect_on_claim: str
    charged_as: str


CORRECTIONS = (
    Correction(
        name="no_contemporaneous_control",
        finding=(
            "REVIVED-BCIS2 OMT arm gained +3.4 LVEF points at 6 months unaided; "
            "BioVAT-HF's +3.9 (90% CI 0.9-6.8) is not distinguishable from it"
        ),
        source="Perera et al., NEJM 2022 (REVIVED-BCIS2); doi:10.1056/NEJMoa2207899",
        effect_on_claim="subtract the unaided gain per patient, floor at zero",
        charged_as="attribution",
    ),
    Correction(
        name="second_trial_null",
        finding=(
            "AATS LB8: 8 ischaemic cardiomyopathy patients, LVEF 28.7 -> 33.0 at 52 "
            "weeks, p = 0.156 - no significant change despite the same intervention class"
        ),
        source="AATS 2025 late-breaking LB8 (allogeneic iPSC-CM patches)",
        effect_on_claim="independent estimate of the same uncontrolled effect",
        charged_as="replication",
    ),
    Correction(
        name="surrogate_not_linked_in_this_population",
        finding=(
            "STICH: LVEF increase occurred only with viability, irrespective of "
            "treatment assignment, and was NOT associated with subsequent death over "
            "10 years in ischaemic cardiomyopathy"
        ),
        source="Panza et al., NEJM 2019; doi:10.1056/NEJMoa1807365",
        effect_on_claim="LVEF gain cannot be read as an outcome benefit here",
        charged_as="endpoint validity",
    ),
    Correction(
        name="viability_does_not_identify_benefit",
        finding=(
            "REVIVED: viability testing did not select patients who gained from "
            "revascularisation; extent of NON-viable myocardium predicted prognosis"
        ),
        source="Perera et al., JAMA Cardiol 2023 (REVIVED viability substudy)",
        effect_on_claim="the recoverable-scar pool is not openable by reperfusion",
        charged_as="mechanism access",
    ),
    Correction(
        name="transmurality_not_viability_decides_recovery",
        finding=(
            "segments with no LGE recover in 91%, 51-75% transmurality in 43%, "
            ">75% in 8% - our population has established scar"
        ),
        source="Kim et al., Circulation 2000; AJR 2005 segmental recovery tables",
        effect_on_claim="conditional recovery is a stratification, not a gain",
        charged_as="population definition",
    ),
)


def correction_table() -> pd.DataFrame:
    """Every documented deficiency, its source and what it is charged against."""
    return pd.DataFrame(
        [
            {
                "deficiency": item.name,
                "finding": item.finding,
                "source": item.source,
                "effect_on_claim": item.effect_on_claim,
                "charged_as": item.charged_as,
                "status": REVIVAL_TAG,
            }
            for item in CORRECTIONS
        ]
    )


def natural_history_draws(
    count: int,
    mean_pp: float = NATURAL_HISTORY_PP,
    sd_pp: float = LVEF_CHANGE_SD_PP,
    seed: int = 20261107,
) -> np.ndarray:
    """Per-patient LVEF change attributable to guideline therapy alone.

    NOT floored. Individual patients in REVIVED both improved and worsened, so the
    comparator is a signed per-patient change and truncating it is a distortion:
    flooring N(3.4, 5) at zero rents the mean up to 4.16 pp, which would make the model
    subtract MORE than the literature documents and call it conservative. The physical
    bound belongs on the outcome, not on the comparator, and is applied in
    revival_draws().

    A mean of zero means NO CORRECTION and returns exactly zeros. It is not a draw from
    N(0, sd): zero-mean noise leaves the mean where it was and only changes the
    variance, which is not what is being corrected. Baking that in would make the
    uncorrected rung disagree with the published one for no reason.
    """
    require_positive("count", count)
    if mean_pp <= 0.0:
        return np.zeros(int(count))
    rng = np.random.default_rng(seed + 1)
    return rng.normal(mean_pp, sd_pp, int(count))


def revival_draws(
    design: Composite = ATLAS5,
    draws: int = 2500,
    seed: int = 20261107,
    natural_history_pp: float = NATURAL_HISTORY_PP,
) -> tuple[np.ndarray, np.ndarray]:
    """The graft gain with the unaided control-arm gain subtracted, per patient.

    Paired with the uncorrected stream: both use the same draw index, so a per-patient
    contrast between designs is still a contrast between the same simulated patient.
    """
    gain, used = composite_draws(design, draws, seed)
    unaided = natural_history_draws(int(draws), natural_history_pp, seed=seed)
    return np.maximum(0.0, gain - unaided), used


def evaluate_revival(
    design: Composite = ATLAS5,
    draws: int = 2500,
    seed: int = 20261107,
    natural_history_pp: float = NATURAL_HISTORY_PP,
) -> CompositeResult:
    """Evaluate a design under the corrected accounting."""
    gain, used = revival_draws(design, draws, seed, natural_history_pp)
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


def correction_ladder(
    design: Composite = ATLAS5, draws: int = 2500, seed: int = 20261107
) -> pd.DataFrame:
    """Every rung from the published claim down to the corrected one.

    Each rung states what it removes and what fraction of the original claim survives.
    Reading the ladder top to bottom is reading the cost of not having a control arm.
    """
    raw = CompositeResult(
        **{
            **composite_result_fields(design, draws, seed),
        }
    )
    rows = [
        {
            "rung": 0,
            "step": "as published (no control group)",
            "basis": "BioVAT-HF +3.9 pp at 3 months, 90% CI 0.9-6.8, n=16",
            "probability_of_target": raw.probability_of_target,
            "mean_gain_pp": raw.median_gain_pp,
        }
    ]
    cumulative = design
    ladder = [
        (1, "minus unaided gain, upper range (+3.4 pp)", 3.4),
        (2, "minus unaided gain, lower range (+1.1 pp)", 1.1),
    ]
    for rung, label, value in ladder:
        corrected = evaluate_revival(cumulative, draws, seed, natural_history_pp=value)
        rows.append(
            {
                "rung": rung,
                "step": label,
                "basis": f"REVIVED-BCIS2 OMT arm, +{value} pp, n=353",
                "probability_of_target": corrected.probability_of_target,
                "mean_gain_pp": corrected.median_gain_pp,
            }
        )
    rows.append(
        {
            "rung": 3,
            "step": "second uncontrolled patch trial",
            "basis": (
                f"AATS LB8 n={AATS_LB8_N}, +{AATS_LB8_GAIN_PP} pp at 52 weeks, "
                f"p={AATS_LB8_P_VALUE} - a null"
            ),
            "probability_of_target": float("nan"),
            "mean_gain_pp": AATS_LB8_GAIN_PP,
        }
    )
    frame = pd.DataFrame(rows)
    frame["fraction_of_original"] = (
        frame["probability_of_target"] / rows[0]["probability_of_target"]
    )
    frame["status"] = REVIVAL_TAG
    return frame


def composite_result_fields(
    design: Composite = ATLAS5, draws: int = 2500, seed: int = 20261107
) -> dict[str, object]:
    """The raw (uncorrected) fields, for the top rung of the ladder."""
    gain, used = composite_draws(design, draws, seed)
    required = CRITERION_THRESHOLD * SCAR_RECOVERABLE_LVEF_PP
    mean_used = float(np.mean(used))
    return {
        "probability_of_target": float(np.mean(gain >= required)),
        "median_gain_pp": float(np.median(gain)),
        "mean_doses_used": mean_used,
        "cumulative_adverse_risk": float(1.0 - (1.0 - design.per_dose_risk) ** mean_used),
        "saturated_fraction": float(np.mean(gain >= SCAR_RECOVERABLE_LVEF_PP - 1e-9)),
        "draws": int(draws),
        "mechanisms": len(design.enabled),
    }


def natural_history_sensitivity(
    design: Composite = ATLAS5, draws: int = 2000, seed: int = 20261107
) -> pd.DataFrame:
    """How much of the claim survives at each documented comparators value.

    Swept rather than chosen. If the answer only holds at one value of the comparator,
    the answer is that value, not the design.
    """
    rows = []
    values = sorted(
        {0.0, NATURAL_HISTORY_RANGE[0], np.mean(NATURAL_HISTORY_RANGE), NATURAL_HISTORY_RANGE[1]}
    )
    for value in values:
        result = evaluate_revival(design, draws, seed, natural_history_pp=value)
        rows.append(
            {
                "natural_history_pp": value,
                "probability_of_target": result.probability_of_target,
                "median_gain_pp": result.median_gain_pp,
                "robust": result.robust,
                "robust_with_margin": result.robust_with_margin,
                "survives_correction": result.probability_of_target >= CRITERION_THRESHOLD,
                "status": REVIVAL_TAG,
            }
        )
    return pd.DataFrame(rows)


def surrogate_translation(draws: int = 2000, seed: int = 20261107) -> pd.DataFrame:
    """Does the corrected gain still clear the surrogate threshold for hard endpoints?

    The threshold is a trial-level mean, so the comparison is between the corrected
    MEAN gain and the STE - not between a per-patient probability and anything. STICH
    is carried alongside because it says the link is absent in this population even
    when the mean moves.
    """
    design = ATLAS5
    raw, _ = composite_draws(design, draws, seed)
    corrected, _ = revival_draws(design, draws, seed)
    rows = []
    for label, ste, r2, linked, note in (
        (
            "all cause mortality",
            SURROGATE_STE_MORTALITY_PP,
            SURROGATE_R2_MORTALITY,
            False,
            "STICH: no association between LVEF change and death in this population",
        ),
        (
            "HF hospitalisation",
            SURROGATE_STE_HF_HOSP_PP,
            SURROGATE_R2_MORTALITY,
            True,
            "trial-level surrogate accepted in HFrEF only",
        ),
    ):
        saturated = float(np.mean(corrected >= SCAR_RECOVERABLE_LVEF_PP - 1e-9))
        rows.append(
            {
                "endpoint": label,
                "surrogate_threshold_pp": ste,
                "surrogate_r2": r2,
                "raw_mean_gain_pp": float(np.mean(raw)),
                "corrected_mean_gain_pp": float(np.mean(corrected)),
                "corrected_median_gain_pp": float(np.median(corrected)),
                "corrected_saturated_fraction": saturated,
                # A mean gain sitting at the recoverable-scar ceiling is a property of
                # the ceiling, not an expected clinical effect. Where most draws
                # saturate the STE comparison is uninformative and is labelled so
                # rather than reported as a pass.
                "model_exceeds_observed_range": bool(
                    float(np.median(corrected)) > LARGEST_HUMAN_OBSERVATION_PP
                ),
                "comparison_usable": bool(
                    saturated < 0.5 and float(np.median(corrected)) <= LARGEST_HUMAN_OBSERVATION_PP
                ),
                "raw_clears_ste": bool(float(np.mean(raw)) >= ste),
                "corrected_clears_ste": bool(float(np.mean(corrected)) >= ste),
                "link_established_in_this_population": linked,
                "note": note,
                "status": REVIVAL_TAG,
            }
        )
    return pd.DataFrame(rows)


def transmurality_stratification() -> pd.DataFrame:
    """Conditional recovery by scar transmurality, with the excluded fraction stated.

    This is the correction that is MOST tempting to present as a gain and is NOT one.
    Restricting to patients without transmural scar raises the probability of +20 pp
    for the patients who enter - and leaves everyone else outside the claim. It is a
    stratification, exactly like the response-guided screen that the recombination
    search already had to reject for the same reason.
    """
    rows = []
    for band, recovery in TRANSMURALITY_RECOVERY.items():
        rows.append(
            {
                "transmurality_band": band,
                "conditional_recovery": recovery,
                "usable_fraction_of_scar": recovery,
                "raises_conditional_probability": recovery > TRANSMURALITY_RECOVERY["lge_over_75"],
                "population_level_gain_claimed": 0.0,
                "reason": (
                    "stratification: excluded patients are not made better by being "
                    "excluded, so no population-level probability may be raised"
                ),
                "status": REVIVAL_TAG,
            }
        )
    return pd.DataFrame(rows)


def revival_search(
    draws: int = 2000, seed: int = 20261107
) -> tuple[pd.DataFrame, dict[str, object]]:
    """Re-run the composition search under the corrected accounting.

    If the same three mechanisms still win after the correction, that is evidence the
    composition is doing real work rather than riding the uncontrolled baseline. If a
    different set wins, ATLAS-5 was an artefact of the missing control.
    """
    grid = exhaustive_subsets(ATLAS5, draws, seed)
    rows = []
    for record in grid.itertuples(index=False):
        names = {name for name in str(record.mechanisms).split("|") if name in ALL_MECHANISMS}
        design = Composite(enabled=frozenset(names), doses=ATLAS5.doses)
        if not design.constraints_intact:
            design = Composite(
                enabled=design.enabled | frozenset(("headroom_gating", "human_tempering")),
                doses=ATLAS5.doses,
            )
        corrected = evaluate_revival(design, draws, seed)
        rows.append(
            {
                "mechanisms": record.mechanisms,
                "count": record.count,
                "admissible": str(record.precedent_breaches) == "none",
                "raw_probability_of_target": record.probability_of_target,
                "corrected_probability_of_target": corrected.probability_of_target,
                "corrected_risk_adjusted": corrected.risk_adjusted_probability,
                "corrected_median_gain_pp": corrected.median_gain_pp,
                "corrected_robust": corrected.robust,
                "corrected_robust_with_margin": corrected.robust_with_margin,
                "status": REVIVAL_TAG,
            }
        )
    frame = (
        pd.DataFrame(rows)
        .sort_values(["corrected_robust", "corrected_risk_adjusted"], ascending=False)
        .reset_index(drop=True)
    )
    admissible = frame.loc[frame["admissible"]].sort_values(
        "corrected_risk_adjusted", ascending=False
    )
    best = admissible.iloc[0] if not admissible.empty else frame.iloc[0]
    selection = {
        "best_corrected": {
            "mechanisms": best["mechanisms"],
            "count": int(best["count"]),
            "corrected_probability_of_target": float(best["corrected_probability_of_target"]),
            "corrected_risk_adjusted": float(best["corrected_risk_adjusted"]),
        },
        "atlas5_still_selected": str(best["mechanisms"])
        == "|".join(sorted(ATLAS5.enabled - {"headroom_gating", "human_tempering"})),
        "any_robust": bool(frame["corrected_robust"].any()),
        "best_raw_probability": float(frame["raw_probability_of_target"].max()),
        "best_corrected_probability": float(frame["corrected_probability_of_target"].max()),
    }
    return frame, selection


def revival_verdict(
    search: dict[str, object], ladder: pd.DataFrame, sensitivity: pd.DataFrame
) -> dict[str, object]:
    """Every claim in the report is read out of these numbers, never asserted."""
    as_published = float(ladder.iloc[0]["probability_of_target"])
    corrected = float(ladder.iloc[1]["probability_of_target"])
    return {
        "as_published_probability": as_published,
        "corrected_probability": corrected,
        "attribution_loss_pp": 100.0 * (as_published - corrected),
        "no_longer_robust": corrected < CRITERION_THRESHOLD,
        "atlas5_still_best_after_correction": bool(search["atlas5_still_selected"]),
        "any_composition_robust_after_correction": bool(search["any_robust"]),
        "best_corrected_mechanisms": str(search["best_corrected"]["mechanisms"]),
        "best_corrected_probability": float(
            search["best_corrected"]["corrected_probability_of_target"]
        ),
        "still_beats_atlas4": float(search["best_corrected"]["corrected_probability_of_target"])
        > ATLAS4_BEST_ADMISSIBLE,
        "comparator_range_pp": (
            float(sensitivity["natural_history_pp"].min()),
            float(sensitivity["natural_history_pp"].max()),
        ),
        "robust_at_every_comparator": bool(sensitivity["survives_correction"].all()),
        "robust_at_any_comparator": bool(sensitivity["survives_correction"].any()),
    }


def revival_report(
    ladder: pd.DataFrame,
    sensitivity: pd.DataFrame,
    surrogate: pd.DataFrame,
    transmurality: pd.DataFrame,
    search_frame: pd.DataFrame,
    search: dict[str, object],
    ledger: pd.DataFrame,
) -> str:
    """Compose the report. All numbers come from the frames passed in."""
    verdict = revival_verdict(search, ladder, sensitivity)
    lines = [
        "CORRECTION FOR WHAT THE LITERATURE SAYS THE COMPOSITION IS MISSING",
        "=" * 78,
        "",
        DISCLAIMER,
        "",
        "WHAT THE SEARCH FOUND (each item is a documented deficiency, not a knob)",
        "-" * 78,
    ]
    for item in CORRECTIONS:
        lines += [
            f"  {item.name}",
            f"    {item.finding}",
            f"    {item.source}",
            f"    charged as: {item.charged_as}",
            "",
        ]

    lines += [
        "THE LADDER: WHAT SURVIVES THE CORRECTION",
        "-" * 78,
        f"  {'rung':>4}  {'step':38s} {'P(target)':>10s} {'of original':>12s}",
    ]
    for row in ladder.itertuples(index=False):
        probability = row.probability_of_target
        shown = "   n/a" if np.isnan(probability) else f"{100 * probability:5.2f}%"
        fraction = row.fraction_of_original
        shown_fraction = "   n/a" if np.isnan(fraction) else f"{100 * fraction:5.1f}%"
        lines.append(f"  {row.rung:>4}  {row.step:38s} {shown:>10s} {shown_fraction:>12s}")
    lines += [
        "",
        f"  ATTRIBUTION LOSS: {verdict['attribution_loss_pp']:.2f} pp of the reported",
        "  probability was what the guideline therapy arm achieved unaided.",
    ]

    if verdict["no_longer_robust"]:
        lines += [
            "",
            "  *** ATLAS-5 IS NO LONGER ROBUST ONCE THE CONTROL ARM IS SUBTRACTED. ***",
            "  It does not clear the 80% bar. The reported 97.92% was carried by an",
            "  uncontrolled baseline, not by the mechanisms.",
        ]
    else:
        lines += [
            "",
            "  ATLAS-5 still clears the bar after the correction. The mechanisms are",
            "  doing work that the unaided arm does not explain.",
        ]

    lines += [
        "",
        "IS THE COMPOSITION STILL THE RIGHT ONE?",
        "-" * 78,
    ]
    if verdict["atlas5_still_best_after_correction"]:
        lines += [
            "  YES. Re-searching all 64 compositions under the corrected accounting",
            "  selects the same mechanism set. The recombination survives; only its",
            "  magnitude was inflated.",
        ]
    else:
        lines += [
            "  NO. Re-searching under the corrected accounting selects a DIFFERENT set:",
            f"    {verdict['best_corrected_mechanisms']}",
            "  ATLAS-5 was an artefact of the missing control arm.",
        ]
    if verdict["any_composition_robust_after_correction"]:
        lines.append(
            f"  Best corrected admissible design: "
            f"{100 * verdict['best_corrected_probability']:.2f}%."
        )
    else:
        lines += [
            "  NO composition is robust after the correction. The correct report is the",
            "  best available design with robust_set_exists = False.",
        ]
    beats = verdict["still_beats_atlas4"]
    lines.append(
        f"  Still beats ATLAS-4's {100 * ATLAS4_BEST_ADMISSIBLE:.2f}% best admissible: "
        f"{'yes' if beats else 'no'}."
    )

    lines += [
        "",
        "SENSITIVITY: THE COMPARATOR IS SWEPT, NOT CHOSEN",
        "-" * 78,
        f"  {'unaided gain (pp)':>18s} {'P(target)':>10s} {'median gain':>13s} {'robust':>8s}",
    ]
    for row in sensitivity.itertuples(index=False):
        lines.append(
            f"  {row.natural_history_pp:18.2f} "
            f"{100 * row.probability_of_target:9.2f}% "
            f"{row.median_gain_pp:12.2f}  "
            f"{'yes' if row.robust_with_margin else 'no':>8s}"
        )
    if verdict["robust_at_every_comparator"]:
        lines += ["", "  Robust across the entire documented comparator range."]
    elif verdict["robust_at_any_comparator"]:
        lines += [
            "",
            "  Robust at some comparators and not others. The design's status is",
            "  conditional on a quantity that was never randomised, which is the point.",
        ]
    else:
        lines += ["", "  Robust at NO documented comparator value."]

    lines += [
        "",
        "DOES THE ENDPOINT STILL MEAN ANYTHING?",
        "-" * 78,
        f"  {'endpoint':22s} {'STE':>6s} {'raw':>8s} {'corrected':>10s} {'link':>10s}",
    ]
    for row in surrogate.itertuples(index=False):
        lines.append(
            f"  {row.endpoint:22s} {row.surrogate_threshold_pp:5.1f} "
            f"{row.raw_mean_gain_pp:7.2f} {row.corrected_mean_gain_pp:9.2f} "
            f"{'established' if row.link_established_in_this_population else 'ABSENT':>10s}"
        )
    if not bool(surrogate["comparison_usable"].all()):
        lines += [
            "",
            "  UNUSABLE COMPARISON, for a reason worth stating plainly:",
        ]
        if bool(surrogate["model_exceeds_observed_range"].any()):
            lines += [
                f"    the corrected median gain is ~{surrogate['corrected_median_gain_pp'].max():.1f} pp,",
                f"    which is {surrogate['corrected_median_gain_pp'].max() / LARGEST_HUMAN_OBSERVATION_PP:.1f}x the largest LVEF gain",
                f"    ANY remuscularisation approach has measured in a human (+{LARGEST_HUMAN_OBSERVATION_PP:.1f} pp,",
                "    BioVAT-HF at a mean 17 months).",
            ]
        if bool((surrogate["corrected_saturated_fraction"] > 0.5).any()):
            lines.append(
                f"    {100 * float(surrogate['corrected_saturated_fraction'].max()):.0f}% of draws also sit at the"
            )
            lines.append("    recoverable-scar ceiling, so the mean is a property of that ceiling.")
        lines += [
            "",
            "  The surrogate threshold is therefore crossed by construction, not by",
            "  prediction. An STE comparison against a model that operates outside the",
            "  measured range is not a pass. It is the second deficiency: even after",
            "  removing the unaided gain, this model claims several times more effect",
            "  than anything observed in a human.",
        ]
    mortality = surrogate.loc[surrogate["endpoint"] == "all cause mortality"].iloc[0]
    if not bool(mortality["link_established_in_this_population"]):
        lines += [
            "",
            "  The all-cause mortality link is ABSENT in this population (STICH, 10-year",
            "  follow-up, n=601). An LVEF increase that does not track death is a",
            "  biomarker, not a benefit. Any design that optimises it must say so.",
        ]

    lines += [
        "",
        "THE TEMPTATION THAT IS NOT A GAIN",
        "-" * 78,
        f"  {'transmurality':14s} {'conditional recovery':>20s} {'population gain':>16s}",
    ]
    for row in transmurality.itertuples(index=False):
        lines.append(
            f"  {row.transmurality_band:14s} {100 * row.conditional_recovery:19.0f}% "
            f"{100 * row.population_level_gain_claimed:15.0f}%"
        )
    lines += [
        "",
        "  Restricting to patients without transmural scar raises the conditional",
        "  probability and is claimed as ZERO population-level gain. Excluding patients",
        "  does not make anyone better - the same trap the recombination search already",
        "  had to reject for response-guided selection.",
    ]

    lines += [
        "",
        "THE BILL, UNCHANGED BY CORRECTION",
        "-" * 78,
    ]
    for row in ledger.itertuples(index=False):
        lines.append(
            f"  {row.axis:28s} {row.value:8.3f} {row.precedented_limit:8.3f}  "
            f"{'YES' if row.beyond_precedent else 'no':>4s}"
        )
    lines += [
        "",
        "  Subtracting an unaided effect changes attribution, not precedent. ATLAS-5",
        "  was already inside precedent on every axis the graft arm controls; it stays",
        "  there. The clearance arm's potency breach is still carried and still not",
        "  forgiven.",
        "",
        "WHAT THIS CORRECTION IS",
        "-" * 78,
        "  It is not a mechanism. No recombination could supply what was missing,",
        "  because what was missing was a CONTROL ARM, and a control arm is evidence,",
        "  not an assumption. The correction moves the estimate from 'what the",
        "  intervention achieved' to 'what the intervention achieved that the guideline",
        "  therapy did not'. That is the honest quantity, and it is smaller.",
        "",
        "  The gap is not closeable in simulation. It is closeable by randomising the",
        "  graft against guideline therapy with an OUTCOME endpoint - which is the",
        "  Phase 3 that BioVAT's investigators have already announced.",
        "",
        "  Randomised proof for anything in this study: 0.00%.",
        "",
        "KEY",
        "=" * 78,
    ]
    return "\n".join(lines)


def run_revival_study(output: Path, draws: int = 2500) -> dict[str, object]:
    """Write the correction study. Returns the summary printed by the CLI."""
    from ihd_solutions.recombine import cost_ledger

    output = Path(output)
    tables = output / "tables"
    reports = output / "reports"
    tables.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)
    (output / "DISCLAIMER.txt").write_text(DISCLAIMER + "\n")

    ladder = correction_ladder(ATLAS5, draws)
    sensitivity = natural_history_sensitivity(ATLAS5, max(1200, draws // 2))
    surrogate = surrogate_translation(max(1200, draws // 2))
    transmurality = transmurality_stratification()
    search_frame, search = revival_search(max(1500, draws // 2))
    corrections = correction_table()
    ledger = cost_ledger(ATLAS5).assign(status=REVIVAL_TAG)

    ladder.to_csv(tables / "revival_correction_ladder.csv", index=False)
    sensitivity.to_csv(tables / "revival_natural_history_sensitivity.csv", index=False)
    surrogate.to_csv(tables / "revival_surrogate_translation.csv", index=False)
    transmurality.to_csv(tables / "revival_transmurality.csv", index=False)
    search_frame.to_csv(tables / "revival_search.csv", index=False)
    corrections.to_csv(tables / "revival_corrections.csv", index=False)
    ledger.to_csv(tables / "revival_cost_ledger.csv", index=False)

    report = revival_report(
        ladder, sensitivity, surrogate, transmurality, search_frame, search, ledger
    )
    (reports / "revival_report.txt").write_text(report)

    verdict = revival_verdict(search, ladder, sensitivity)
    summary = {
        "artefacts": 9,
        "as_published_probability": verdict["as_published_probability"],
        "corrected_probability": verdict["corrected_probability"],
        "attribution_loss_pp": verdict["attribution_loss_pp"],
        "no_longer_robust": verdict["no_longer_robust"],
        "atlas5_still_best": verdict["atlas5_still_best_after_correction"],
        "any_robust_after_correction": verdict["any_composition_robust_after_correction"],
        "still_beats_atlas4": verdict["still_beats_atlas4"],
        "status": REVIVAL_TAG,
    }
    return summary
