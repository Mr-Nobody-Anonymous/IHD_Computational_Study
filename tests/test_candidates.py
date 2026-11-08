"""The ranking must be reproducible and must not flatter the repository.

The failure mode this file exists to prevent: a module that nominates a favourite and calls
it best. These tests pin the *rule*, the ordering it produces, and the two inconvenient
facts - that the best answer is an old, generic intervention, and that this repository's own
architecture ranks last because it has no clinical outcome data at all.
"""

from __future__ import annotations

import pandas as pd
import pytest

from ihd_solutions import candidates as C

# ------------------------------------------------------------------- structure


def test_seven_candidates_are_registered():
    assert len(C.CANDIDATES) == 7
    assert len({c.key for c in C.CANDIDATES}) == 7


def test_every_candidate_cites_a_source():
    for candidate in C.CANDIDATES:
        assert len(candidate.source) > 10, candidate.key


def test_every_candidate_states_what_it_does_not_license():
    for candidate in C.CANDIDATES:
        assert len(candidate.what_it_does_not) > 60, candidate.key


def test_tier_map_is_the_single_source_of_truth():
    """The rank order must be derived from TIERS, never restated alongside it."""
    assert set(C.TIER_ORDER) == set(C.TIERS)
    assert sorted(C.TIER_ORDER.values()) == list(range(1, len(C.TIERS) + 1))
    assert all(c.evidence_tier in C.TIERS for c in C.CANDIDATES)


def test_unknown_tier_is_rejected():
    with pytest.raises(ValueError, match="unknown evidence tier"):
        C.Candidate(
            key="x",
            name="x",
            mechanism="x",
            best_trial="none",
            n=10,
            endpoint="none",
            endpoint_class="surrogate",
            result="none",
            effect=None,
            effect_unit="percentage points",
            evidence_tier="tier_0_asserted",
            direction="benefit",
            what_it_licenses="nothing",
            what_it_does_not="nothing at all, which is the point of the field",
            source="nowhere",
        )


def test_unknown_direction_is_rejected():
    with pytest.raises(ValueError, match="unknown direction"):
        C.Candidate(
            key="x",
            name="x",
            mechanism="x",
            best_trial="none",
            n=10,
            endpoint="none",
            endpoint_class="surrogate",
            result="none",
            effect=None,
            effect_unit="not estimable - no controlled comparison",
            evidence_tier="tier_5_no_clinical_data",
            direction="hopeful",
            what_it_licenses="nothing",
            what_it_does_not="nothing at all, which is the point of the field",
            source="somewhere",
        )


def test_a_candidate_must_state_what_it_does_not_license():
    with pytest.raises(ValueError, match="does not license"):
        C.Candidate(
            key="x",
            name="x",
            mechanism="x",
            best_trial="none",
            n=10,
            endpoint="none",
            endpoint_class="surrogate",
            result="none",
            effect=None,
            effect_unit="not estimable - no controlled comparison",
            evidence_tier="tier_5_no_clinical_data",
            direction="none",
            what_it_licenses="nothing",
            what_it_does_not="",
            source="somewhere",
        )


def test_unknown_endpoint_class_is_rejected():
    with pytest.raises(ValueError, match="unknown endpoint class"):
        C.Candidate(
            key="x",
            name="x",
            mechanism="x",
            best_trial="none",
            n=10,
            endpoint="creatinine",
            endpoint_class="biomarker",
            result="none",
            effect=None,
            effect_unit="not estimable - no controlled comparison",
            evidence_tier="tier_5_no_clinical_data",
            direction="none",
            what_it_licenses="nothing",
            what_it_does_not="nothing at all, which is the point of the field",
            source="nowhere",
        )


def test_candidate_without_a_source_is_rejected():
    with pytest.raises(ValueError, match="must cite a source"):
        C.Candidate(
            key="x",
            name="x",
            mechanism="x",
            best_trial="none",
            n=10,
            endpoint="none",
            endpoint_class="surrogate",
            result="none",
            effect=None,
            effect_unit="not estimable - no controlled comparison",
            evidence_tier="tier_5_no_clinical_data",
            direction="none",
            what_it_licenses="nothing",
            what_it_does_not="nothing at all, which is the point of the field",
            source="",
        )


def test_infinite_effect_is_rejected():
    with pytest.raises(ValueError, match="not infinity"):
        C.Candidate(
            key="x",
            name="x",
            mechanism="x",
            best_trial="none",
            n=10,
            endpoint="none",
            endpoint_class="surrogate",
            result="none",
            effect=float("inf"),
            effect_unit="relative LDL-C reduction",
            evidence_tier="tier_5_no_clinical_data",
            direction="benefit",
            what_it_licenses="nothing",
            what_it_does_not="nothing at all, which is the point of the field",
            source="somewhere",
        )


# ------------------------------------------------------------------- the table


def test_table_is_ranked_one_to_seven():
    frame = C.candidate_table()
    assert list(frame["rank"]) == list(range(1, 8))
    assert frame["rank_key"].is_monotonic_increasing


def test_three_candidates_are_measured_on_an_outcome():
    frame = C.candidate_table()
    outcome = frame.loc[frame["endpoint_class"] == "outcome"]
    assert set(outcome["key"]) == {"ldl_lowering_approved", "colchicine", "il6_inhibition"}
    assert len(outcome) == 3


def test_only_two_are_backed_by_that_outcome_evidence():
    """Being measured on an outcome and being supported by one are different claims."""
    answer = C.best_candidate()
    assert answer["candidates_with_outcome_evidence"] == 2
    assert answer["candidates_refuted_by_outcome"] == 1


def test_a_refutation_is_not_counted_as_a_benefit():
    """The bug this pins: effect = 1 - HR turned ZEUS into a 1% win."""
    frame = C.candidate_table().set_index("key")
    assert frame.loc["il6_inhibition", "effect"] == pytest.approx(0.99)
    assert frame.loc["il6_inhibition", "direction"] == "neutral"
    assert frame.loc["il6_inhibition", "direction_rank"] > 0


def test_effect_units_are_recorded_for_every_candidate():
    for candidate in C.CANDIDATES:
        assert candidate.effect_unit, candidate.key


def test_effects_are_only_rendered_in_their_own_unit():
    frame = C.candidate_table().set_index("key")
    assert frame.loc["il6_inhibition", "effect_display"].startswith("HR 0.99")
    assert "percentage points" not in frame.loc["colchicine", "effect_display"]
    assert frame.loc["cell_patch", "effect_display"].startswith("+3.9")
    assert frame.loc["atlas5_combination", "effect_display"].startswith("not estimable")


def test_il6_inhibition_is_classified_as_refuted():
    frame = C.candidate_table().set_index("key")
    assert frame.loc["il6_inhibition", "rank_key"] == C.TIER_ORDER["tier_2_refuted"]
    assert frame.loc["il6_inhibition", "direction"] == "neutral"


def test_best_candidate_is_ldl_lowering():
    answer = C.best_candidate()
    assert "LDL" in str(answer["best_candidate"])
    assert answer["candidates_total"] == 7
    assert answer["status"] == "ranking_by_stated_evidence_rule"


def test_repositorys_own_architecture_ranks_last():
    """It has no clinical data. The rule puts that at the bottom, without exception."""
    frame = C.candidate_table()
    last = frame.iloc[-1]
    assert last["key"] == "atlas5_combination"
    assert last["rank"] == 7
    assert last["n"] == 0
    assert last["rank_key"] == max(C.TIER_ORDER.values())
    assert last["direction"] == "none"


def test_no_candidate_is_called_a_cure():
    """The only permitted appearance of the word is inside a refusal of it."""
    answer = C.best_candidate()
    assert "cure" not in str(answer["best_candidate"]).lower()
    refusal = str(answer["this_is_not_a_cure"]).lower()
    assert "no candidate here is a cure" in refusal
    report = C.candidate_report().lower()
    for sentence in report.splitlines():
        if "cure" in sentence:
            assert "not a cure" in sentence or "is a cure" in sentence, sentence


def test_table_is_a_frame_with_the_documented_columns():
    frame = C.candidate_table()
    assert isinstance(frame, pd.DataFrame)
    for column in ("rank", "candidate", "tier", "endpoint_class", "does_not_license"):
        assert column in frame.columns


# ------------------------------------------------------- the ZEUS refutation


def test_inflammation_refutation_carries_the_trial_numbers():
    refutation = C.inflammation_arm_refutation()
    assert refutation["n"] == 6_385
    assert refutation["hazard_ratio"] == pytest.approx(0.99)
    assert refutation["ci"] == pytest.approx((0.88, 1.11))
    assert refutation["target_engagement"]["hsCRP_reduction"] == pytest.approx(0.64)
    assert refutation["target_engagement"]["IL6_reduction"] == pytest.approx(0.76)


def test_the_refuted_inference_is_about_inference_not_measurement():
    refutation = C.inflammation_arm_refutation()
    assert "marker" in refutation["the_refuted_inference"]
    assert "outcome" in refutation["the_refuted_inference"]
    assert "marker dynamics" in refutation["action_for_this_repository"]


def test_crp_is_demoted_but_not_thrown_away():
    verdict = C.crp_is_no_longer_a_surrogate()
    assert "non-identifiability" in verdict["still_valid_as"]
    assert "benefit" in verdict["no_longer_valid_as"]
    assert "6,385" in verdict["evidence"] or "6385" in verdict["evidence"]


# ------------------------------------------------------------------- artefacts


def test_report_contains_every_candidate():
    report = C.candidate_report()
    for candidate in C.CANDIDATES:
        assert candidate.name in report


def test_report_states_the_rule_before_the_ranking():
    report = C.candidate_report()
    assert report.index("THE RULE") < report.index("THE RANKING")


def test_report_lists_every_tier_with_its_count():
    report = C.candidate_report()
    for label in C.TIERS.values():
        assert label in report


def test_write_report_produces_six_artefacts_plus_a_summary(tmp_path):
    summary = C.write_candidate_report(tmp_path)
    assert summary["artefacts"] == 6
    assert summary["candidates_ranked"] == 7
    assert summary["repository_architecture_rank"] == 7
    assert summary["inflammation_arm_status"] == "refuted by ZEUS"
    written = sorted(p.name for p in (tmp_path / "tables").iterdir())
    assert len(written) == 4
    assert (tmp_path / "candidate_report.txt").exists()
    assert (tmp_path / "DISCLAIMER.txt").exists()


def test_disclaimer_refuses_the_cure_framing():
    assert "not a cure" in C.DISCLAIMER.lower()
