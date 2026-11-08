"""Test whether a fitted parameter can be identified at all, instead of reporting one fit.

The critique names this precisely:

    *"two different combinations of plaque-growth rate and inflammation rate might
    generate nearly identical outputs. If that happens, the model cannot reliably
    identify which biological mechanism caused the observed result. The repository must
    expose this problem instead of silently reporting one fitted parameter set."*

That is correct, and it is the failure mode that makes every 'calibrated' model in the
literature untrustworthy when nobody checks. This module implements the three checks that
actually decide whether a fitted parameter means anything:

1. **Local sensitivity** - does perturbing the parameter move the output at all? A
   parameter with zero sensitivity cannot be fitted and is decorative.
2. **Collinearity** - do two parameters move the output in the same direction? If so,
   only their sum (or ratio) is identified, not either one.
3. **Profile identifiability** - sweep one parameter, re-fit the others at each point,
   and see whether the objective has a single minimum. A flat valley means the data
   cannot separate the parameters.

The verdict is deliberately pessimistic. A parameter is only called ``identified`` if it
passes all three. Everything else is reported as what it is: a number that was chosen, not
measured, and that no amount of downstream simulation can promote.

Everything here is deterministic given a seed; nothing is tuned.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
import pandas as pd

#: Sensitivity below this cannot be distinguished from solver/round-off noise.
SENSITIVITY_FLOOR = 1e-6

#: |correlation| above this means two parameters move together and are not separable.
COLLINEARITY_CEILING = 0.95

#: A profile is "identified" only if the objective rises on BOTH sides of the optimum.
PROFILE_RISE_FRACTION = 0.05


@dataclass(frozen=True)
class ParameterSensitivity:
    """One-at-a-time sensitivity of an output to a parameter."""

    name: str
    nominal: float
    elasticity: float
    direction: str

    @property
    def moves_output(self) -> bool:
        return abs(self.elasticity) > SENSITIVITY_FLOOR

    @property
    def verdict(self) -> str:
        if not self.moves_output:
            return "inert - cannot be fitted, and cannot have caused the result"
        if abs(self.elasticity) > 1.0:
            return "dominant - the result is mostly this parameter's assumption"
        return "sensitive"


def one_at_a_time(
    model: Callable[[dict[str, float]], float],
    nominal: dict[str, float],
    *,
    step: float = 0.01,
) -> list[ParameterSensitivity]:
    """Elasticity of the output to each parameter, holding the others at nominal.

    Elasticity is a relative change divided by a relative change, so parameters on
    different scales (a rate, a fraction, a count) are directly comparable. The step is
    PROPORTIONAL, not absolute: adding 0.01 to a parameter of order 1e6 is a 1e-8
    perturbation and would be lost in floating point.
    """
    base = float(model(dict(nominal)))
    if base == 0.0:
        raise ValueError("nominal output is zero; elasticity is undefined")
    rows = []
    for name, value in nominal.items():
        if value == 0.0:
            raise ValueError(f"{name} is zero; a relative step is undefined")
        bumped = dict(nominal)
        bumped[name] = value * (1.0 + step)
        moved = float(model(bumped))
        elasticity = ((moved - base) / base) / step
        rows.append(
            ParameterSensitivity(
                name=name,
                nominal=value,
                elasticity=float(elasticity),
                direction="increases" if elasticity > 0 else "decreases",
            )
        )
    return sorted(rows, key=lambda row: -abs(row.elasticity))


def sensitivity_matrix(
    model: Callable[[dict[str, float]], float],
    nominal: dict[str, float],
    *,
    step: float = 0.01,
) -> pd.DataFrame:
    """Normalised Jacobian: how each parameter moves the output, per unit relative change."""
    base = float(model(dict(nominal)))
    columns = {}
    for name, value in nominal.items():
        bumped = dict(nominal)
        bumped[name] = value * (1.0 + step)
        columns[name] = ((float(model(bumped)) - base) / base) / step
    return pd.DataFrame([columns])


def jacobian(
    model: Callable[[dict[str, float]], float],
    nominal: dict[str, float],
    *,
    step: float = 0.01,
    points: int = 9,
) -> tuple[np.ndarray, list[str]]:
    """Normalised Jacobian sampled at several points around the nominal.

    One row per sample point, one column per parameter: the relative change in the output
    per unit relative change in that parameter, holding the rest fixed. Sampling at more
    than one point is what makes collinearity detectable at all - at a single point there
    is no variation to correlate, and a one-row matrix carries no covariance information.
    """
    names = list(nominal)
    rng = np.random.default_rng(20261107)
    rows = []
    for _ in range(points):
        # Perturb every parameter by up to +-10%, clipped at positive values.
        centre = {
            name: max(nominal[name] * (1.0 + rng.uniform(-0.1, 0.1)), 1e-12) for name in names
        }
        base = float(model(dict(centre)))
        if base == 0.0:
            continue
        row = np.zeros(len(names))
        for index, name in enumerate(names):
            up, down = dict(centre), dict(centre)
            up[name] = centre[name] * (1.0 + step)
            down[name] = centre[name] * (1.0 - step)
            row[index] = (float(model(up)) - float(model(down))) / (2.0 * step * base)
        rows.append(row)
    if not rows:
        raise ValueError("the model returned zero at every sampled point")
    return np.vstack(rows), names


def effect_correlation(effects: np.ndarray, names: list[str]) -> pd.DataFrame:
    """Correlate each parameter's effect ACROSS sample points.

    Two parameters are collinear when their columns move together across the sampled
    space: whatever one does to the output, the other does proportionally. That is the
    definition that catches an additive cascade, where every production rate has the same
    effect everywhere and can be exchanged for another.
    """
    if effects.shape[0] < 2:
        return pd.DataFrame(np.eye(len(names)), index=names, columns=names)
    with np.errstate(invalid="ignore", divide="ignore"):
        correlation = np.corrcoef(effects, rowvar=False)
    if np.ndim(correlation) == 0:
        correlation = np.ones((1, 1))
    return pd.DataFrame(np.nan_to_num(correlation, nan=0.0), index=names, columns=names)


def collinearity(
    model: Callable[[dict[str, float]], float],
    nominal: dict[str, float],
    *,
    step: float = 0.01,
    points: int = 9,
) -> pd.DataFrame:
    """Correlation between parameter EFFECTS across the sampled parameter space."""
    effects, names = jacobian(model, nominal, step=step, points=points)
    return effect_correlation(effects, names)


def collinear_pairs(
    model: Callable[[dict[str, float]], float],
    nominal: dict[str, float],
    *,
    ceiling: float = COLLINEARITY_CEILING,
    step: float = 0.01,
) -> list[tuple[str, str, float]]:
    """Pairs whose effects cannot be told apart."""
    frame = collinearity(model, nominal, step=step)
    found = []
    names = list(frame.index)
    for i, first in enumerate(names):
        for second in names[i + 1 :]:
            value = float(frame.loc[first, second])
            if abs(value) >= ceiling:
                found.append((first, second, value))
    return found


def profile_identifiability(
    objective: Callable[[dict[str, float]], float],
    nominal: dict[str, float],
    name: str,
    *,
    span: float = 0.5,
    points: int = 21,
) -> pd.DataFrame:
    """Sweep one parameter and RE-FIT the others at each point.

    The sweep is CENTRED ON THE NOMINAL VALUE, so the nominal must be a plausible
    optimum. A test that plants a truth far outside the swept range legitimately reports
    an edge optimum, which is the right answer for that test but not a bug here.

    This is the check that matters. If the objective is flat over a range of ``name``
    while the others are re-optimised, the data do not determine ``name`` - only some
    combination do. Practically it is done here by coordinate descent, which is enough
    to reveal a valley and is deliberately conservative.
    """
    others = [key for key in nominal if key != name]
    values = np.linspace(nominal[name] * (1.0 - span), nominal[name] * (1.0 + span), points)
    rows = []
    for value in values:
        candidate = dict(nominal)
        candidate[name] = float(value)
        # One round of coordinate refinement on the others: enough to descend into the
        # valley, and cheap enough to sweep. A full optimiser would look more rigorous
        # without changing the qualitative reading.
        for _ in range(3):
            for other in others:
                best_value = candidate[other]
                best = objective(candidate)
                for factor in (0.9, 0.95, 1.0, 1.05, 1.1):
                    probe = dict(candidate)
                    probe[other] = candidate[other] * factor
                    score = objective(probe)
                    if score < best:
                        best = score
                        best_value = probe[other]
                candidate[other] = best_value
        rows.append(
            {
                "parameter": name,
                "value": float(value),
                "relative": float(value / nominal[name]),
                "objective": float(objective(candidate)),
                "refitted": "|".join(f"{key}={candidate[key]:.4g}" for key in others),
            }
        )
    return pd.DataFrame(rows)


def profile_verdict(profile: pd.DataFrame) -> str:
    """Decide whether a profile has a real minimum or a flat valley.

    Two degenerate cases return immediately, because inside either of them the notion of
    a minimum is undefined and a shape test over an all-equal array would silently pass:
    a completely flat profile, and an ALL-ZERO profile - which is what a valley looks
    like when the other parameters absorb the swept one exactly.
    """
    values = profile["objective"].to_numpy(dtype=float)
    span = float(values.max() - values.min())
    # An all-zero profile is a valley, not a flat line: the objective is achievable
    # everywhere because the other parameters absorb the swept one exactly. It must be
    # tested BEFORE the flat check, which would otherwise catch it and describe it wrongly.
    if float(np.abs(values).max()) <= 1e-12:
        return "unidentified - a perfect valley: the other parameters absorb this one"
    if span <= 0.0:
        return "unidentified - the objective is completely flat"
    best = int(np.argmin(values))
    at_edge = best == 0 or best == len(values) - 1
    rise = (values[best:].max() if best < len(values) - 1 else 0.0) - values[best]
    fall = (values[: best + 1].max() if best > 0 else 0.0) - values[best]
    threshold = PROFILE_RISE_FRACTION * span
    if at_edge:
        return "unidentified - the optimum is at the edge of the swept range"
    if rise < threshold or fall < threshold:
        return "weakly identified - the objective rises on one side only"
    return "identified - a single interior minimum"


def identifiability_verdict(
    sensitivities: list[ParameterSensitivity],
    pairs: list[tuple[str, str, float]],
    profiles: dict[str, str],
) -> dict[str, object]:
    """Roll the three checks into one classification per parameter."""
    collinear_names = {name for pair in pairs for name in pair[:2]}
    identified, chosen = [], []
    for row in sensitivities:
        ok = (
            row.moves_output
            and row.name not in collinear_names
            and profiles.get(row.name, "").startswith("identified")
        )
        (identified if ok else chosen).append(row.name)
    return {
        "identified": tuple(identified),
        "chosen_not_measured": tuple(chosen),
        "inert": tuple(row.name for row in sensitivities if not row.moves_output),
        "collinear_pairs": tuple(pairs),
        "fraction_identified": len(identified) / max(len(sensitivities), 1),
    }


def identifiability_report(
    sensitivities: list[ParameterSensitivity],
    pairs: list[tuple[str, str, float]],
    profiles: dict[str, str],
    *,
    subject: str = "model",
) -> str:
    """Report the classification, naming every parameter that is chosen not measured."""
    verdict = identifiability_verdict(sensitivities, pairs, profiles)
    lines = [
        f"PARAMETER IDENTIFIABILITY - {subject}",
        "=" * 78,
        "",
        "A parameter is 'identified' only if it moves the output, is not collinear with",
        "another parameter, and its profile has a single interior minimum. Anything else",
        "was CHOSEN, and a chosen parameter cannot be reported as a fitted one.",
        "",
        "SENSITIVITY (elasticity: relative output change per relative parameter change)",
        "-" * 78,
        f"  {'parameter':22s} {'nominal':>10s} {'elasticity':>12s}  verdict",
    ]
    for row in sensitivities:
        lines.append(f"  {row.name:22s} {row.nominal:10.4g} {row.elasticity:12.3f}  {row.verdict}")
    lines += ["", "COLLINEARITY", "-" * 78]
    if pairs:
        for first, second, value in pairs:
            lines.append(
                f"  {first} and {second} correlate at r={value:+.3f}: only their "
                f"combination is determined, not either one"
            )
    else:
        lines.append("  no pair exceeds the collinearity ceiling")
    lines += ["", "PROFILE IDENTIFIABILITY", "-" * 78]
    for name, verdict_text in profiles.items():
        lines.append(f"  {name:22s} {verdict_text}")
    lines += [
        "",
        "CLASSIFICATION",
        "-" * 78,
        f"  identified        {list(verdict['identified'])}",
        f"  chosen, NOT measured  {list(verdict['chosen_not_measured'])}",
        f"  inert             {list(verdict['inert'])}",
        "",
        f"  {100 * verdict['fraction_identified']:.0f}% of parameters are identified.",
    ]
    if verdict["chosen_not_measured"]:
        lines += [
            "",
            "  The parameters in 'chosen, NOT measured' cannot be fitted by the available",
            "  data. Reporting a single fitted value for them would be reporting an",
            "  assumption with a decimal point.",
        ]
    return "\n".join(lines)
