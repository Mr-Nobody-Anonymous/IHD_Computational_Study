"""The evidence interlock: what may and may not move an event curve."""

from __future__ import annotations

import math
from dataclasses import replace

import pytest

from ihd_solutions import evidence
from ihd_solutions.evidence import (
    CreditMode,
    EndpointClass,
    EvidenceTier,
    TrialAnchor,
    decide_credit,
    load_anchors,
    select_anchor,
)


def make_anchor(**changes) -> TrialAnchor:
    base = TrialAnchor(
        anchor_id="TEST",
        therapy_id="drug",
        population_id="pop",
        endpoint="3-point MACE",
        effect_unit=evidence.EffectUnit.FIXED,
        endpoint_class=EndpointClass.MACE,
        hr=0.80,
        ci_low=0.72,
        ci_high=0.90,
        followup_years=3.0,
        n_participants=1000,
        evidence_tier=EvidenceTier.RANDOMIZED_OUTCOMES,
    )
    return replace(base, **changes)


def test_log_sigma_reproduces_a_log_symmetric_interval():
    # Reported intervals are rarely log-symmetric; use one that is by construction.
    anchor = make_anchor(hr=0.80, ci_low=0.72, ci_high=0.80**2 / 0.72)
    sigma = anchor.log_sigma
    assert math.exp(math.log(0.80) - 1.959963984540054 * sigma) == pytest.approx(0.72)
    assert math.exp(math.log(0.80) + 1.959963984540054 * sigma) == pytest.approx(0.80**2 / 0.72)


def test_log_sigma_of_an_asymmetric_interval_centres_on_the_geometric_mean():
    anchor = make_anchor()
    sigma = anchor.log_sigma
    lower = math.exp(math.log(0.80) - 1.959963984540054 * sigma)
    upper = math.exp(math.log(0.80) + 1.959963984540054 * sigma)
    assert math.sqrt(lower * upper) == pytest.approx(0.80)


@pytest.mark.parametrize(
    "hr,low,high,neutral",
    [
        (0.69, 0.57, 0.83, False),
        (0.99, 0.88, 1.11, True),
        (1.01, 0.82, 1.25, True),
        (0.85, 0.9, 1.1, False),
    ],
)
def test_neutral_detection(hr, low, high, neutral):
    if not low <= hr <= high:
        with pytest.raises(ValueError):
            make_anchor(hr=hr, ci_low=low, ci_high=high)
    else:
        assert make_anchor(hr=hr, ci_low=low, ci_high=high).is_neutral is neutral


def test_anchor_requires_ci_with_point_estimate():
    with pytest.raises(ValueError, match="both CI bounds"):
        make_anchor(ci_low=None, ci_high=None)


def test_anchor_rejects_ci_without_point_estimate():
    with pytest.raises(ValueError, match="CI without a point estimate"):
        make_anchor(hr=None)


def test_hard_endpoint_requires_effect_estimate():
    with pytest.raises(ValueError, match="hard endpoint requires"):
        make_anchor(hr=None, ci_low=None, ci_high=None)


def test_surrogate_anchor_may_omit_effect_estimate():
    anchor = make_anchor(hr=None, ci_low=None, ci_high=None, endpoint_class=EndpointClass.SURROGATE)
    assert not anchor.has_effect_estimate
    assert not anchor.is_neutral
    with pytest.raises(ValueError, match="no effect estimate"):
        anchor.scaled_hr()
    with pytest.raises(ValueError, match="no effect estimate"):
        _ = anchor.log_sigma


def test_anchor_validates_identifiers_and_counts():
    with pytest.raises(ValueError):
        make_anchor(anchor_id="")
    with pytest.raises(ValueError):
        make_anchor(n_participants=0)
    with pytest.raises(ValueError):
        make_anchor(followup_years=0.0)
    with pytest.raises(ValueError):
        make_anchor(hr=0.0)


def test_fixed_anchor_ignores_magnitude_but_per_unit_scales():
    assert make_anchor().scaled_hr(5.0) == pytest.approx(0.80)
    per_unit = make_anchor(effect_unit=evidence.EffectUnit.PER_MMOL_L_LDL, hr=0.78)
    assert per_unit.scaled_hr(1.0) == pytest.approx(0.78)
    assert per_unit.scaled_hr(2.0) == pytest.approx(0.78**2)
    assert per_unit.scaled_hr(0.0) == pytest.approx(1.0)
    with pytest.raises(ValueError):
        per_unit.scaled_hr(-1.0)


def test_transport_requires_declaration_outside_exploratory_mode():
    anchor = make_anchor(transportable_to=("other",))
    assert anchor.applies_to("pop", CreditMode.STRICT)
    assert anchor.applies_to("other", CreditMode.STRICT)
    assert not anchor.applies_to("third", CreditMode.STRICT)
    assert not anchor.applies_to("third", CreditMode.SURROGATE)
    assert anchor.applies_to("third", CreditMode.EXPLORATORY)


def test_decide_credit_grants_randomised_mace_anchor():
    anchors = {"TEST": make_anchor()}
    decision = decide_credit("drug", "pop", CreditMode.STRICT, anchors)
    assert decision.credited
    assert decision.anchor_id == "TEST"
    assert decision.hr == pytest.approx(0.80)


def test_decide_credit_denies_when_no_anchor_exists():
    decision = decide_credit("drug", "pop", CreditMode.STRICT, {}, biomarker_only_change=0.9)
    assert not decision.credited
    assert "biomarker change alone never credits events" in decision.reason


@pytest.mark.parametrize(
    "endpoint_class",
    [
        EndpointClass.SURROGATE,
        EndpointClass.FUNCTIONAL,
        EndpointClass.PENDING,
        EndpointClass.HF_COMPOSITE,
    ],
)
def test_non_mace_endpoints_never_credited_in_any_mode(endpoint_class):
    anchors = {"TEST": make_anchor(endpoint_class=endpoint_class)}
    for mode in CreditMode:
        decision = decide_credit("drug", "pop", mode, anchors)
        assert not decision.credited
        assert "cannot be credited as MACE" in decision.reason


@pytest.mark.parametrize("endpoint_class", [EndpointClass.CV_DEATH, EndpointClass.ALL_CAUSE_DEATH])
def test_mortality_endpoints_only_credited_in_exploratory_mode(endpoint_class):
    anchors = {"TEST": make_anchor(endpoint_class=endpoint_class)}
    assert not decide_credit("drug", "pop", CreditMode.STRICT, anchors).credited
    assert not decide_credit("drug", "pop", CreditMode.SURROGATE, anchors).credited
    assert decide_credit("drug", "pop", CreditMode.EXPLORATORY, anchors).credited


def test_observational_evidence_requires_exploratory_mode():
    anchors = {"TEST": make_anchor(evidence_tier=EvidenceTier.OBSERVATIONAL)}
    assert not decide_credit("drug", "pop", CreditMode.STRICT, anchors).credited
    assert not decide_credit("drug", "pop", CreditMode.SURROGATE, anchors).credited
    decision = decide_credit("drug", "pop", CreditMode.EXPLORATORY, anchors)
    assert decision.credited
    assert "observational_confounding" in decision.warnings


def test_mechanistic_tier_is_never_creditable():
    anchors = {"TEST": make_anchor(evidence_tier=EvidenceTier.MECHANISTIC)}
    for mode in CreditMode:
        assert not decide_credit("drug", "pop", mode, anchors).credited


def test_neutral_anchor_is_credited_but_flagged():
    anchors = {"TEST": make_anchor(hr=0.99, ci_low=0.88, ci_high=1.11)}
    decision = decide_credit("drug", "pop", CreditMode.STRICT, anchors)
    assert decision.credited
    assert "neutral_confidence_interval" in decision.warnings


def test_transported_population_is_flagged():
    anchors = {"TEST": make_anchor(transportable_to=("elsewhere",))}
    decision = decide_credit("drug", "elsewhere", CreditMode.STRICT, anchors)
    assert decision.credited
    assert "transported_population" in decision.warnings


def test_select_anchor_prefers_exact_population_match():
    exact = make_anchor(anchor_id="EXACT", population_id="pop")
    other = make_anchor(anchor_id="OTHER", population_id="other", transportable_to=("pop",))
    anchors = {"EXACT": exact, "OTHER": other}
    assert select_anchor(anchors, "drug", "pop", CreditMode.STRICT).anchor_id == "EXACT"
    assert select_anchor(anchors, "drug", "pop", CreditMode.SURROGATE).anchor_id == "EXACT"


def test_biomarker_pathway_credit_is_forbidden():
    decision = evidence.CreditDecision(
        "drug", "pop", CreditMode.STRICT, True, "bad", pathway="biomarker"
    )
    with pytest.raises(AssertionError, match="biomarker-only credit"):
        evidence.assert_no_biomarker_surrogate_credit([decision])
    evidence.assert_no_biomarker_surrogate_credit([replace(decision, pathway="anchor")])


def test_registry_loads_and_is_internally_consistent(anchors, catalog, populations):
    assert len(anchors) >= 25
    for anchor in anchors.values():
        assert anchor.therapy_id in catalog
        assert anchor.population_id in populations
        for transported in anchor.transportable_to:
            assert transported in populations
        assert anchor.source_url.startswith("https://")


def test_registry_contains_the_key_neutral_controls(anchors):
    neutral = {anchor.anchor_id for anchor in anchors.values() if anchor.is_neutral}
    assert {"ZEUS", "CIRT", "CLEAR-SYNERGY", "ISCHEMIA"} <= neutral


def test_registry_rejects_missing_file_and_bad_columns(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_anchors(tmp_path / "missing.csv")
    bad = tmp_path / "bad.csv"
    bad.write_text("anchor_id,therapy_id\nA,B\n", encoding="utf-8")
    with pytest.raises(ValueError, match="missing columns"):
        load_anchors(bad)


def test_registry_rejects_duplicate_anchor_ids(tmp_path, anchors):
    import pandas as pd

    frame = pd.DataFrame(
        [
            {
                "anchor_id": anchor.anchor_id,
                "therapy_id": anchor.therapy_id,
                "population_id": anchor.population_id,
                "endpoint": anchor.endpoint,
                "effect_unit": anchor.effect_unit.value,
                "endpoint_class": anchor.endpoint_class.value,
                "hr": anchor.hr,
                "ci_low": anchor.ci_low,
                "ci_high": anchor.ci_high,
                "followup_years": anchor.followup_years,
                "n_participants": anchor.n_participants,
                "evidence_tier": anchor.evidence_tier.value,
                "transportable_to": ";".join(anchor.transportable_to),
                "source_url": anchor.source_url,
                "notes": anchor.notes,
            }
            for anchor in anchors.values()
        ]
    )
    duplicated = pd.concat([frame, frame.iloc[[0]]], ignore_index=True)
    path = tmp_path / "dup.csv"
    duplicated.to_csv(path, index=False)
    with pytest.raises(ValueError, match="duplicate anchor_id"):
        load_anchors(path)


def test_tier_summary_counts_every_tier(anchors):
    frame = evidence.summarise_tiers(anchors)
    assert set(frame["evidence_tier"]) == {tier.value for tier in EvidenceTier}
    assert int(frame["n_anchors"].sum()) == len(anchors)
