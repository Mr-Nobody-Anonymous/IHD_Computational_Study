"""Validate the roadmap's estimand and external-data intake registries.

This is Phase 0/1 planning infrastructure. It validates schemas and cross-references;
source truth, data access, model calibration, and biological validity require human/data
review and are never certified by this module.
"""

from __future__ import annotations

from pathlib import Path
from urllib.parse import urlparse

import pandas as pd

from ihd_solutions.evidence import DEFAULT_DATA_DIR

ESTIMAND_COLUMNS = (
    "estimand_id",
    "disease_stage",
    "population_id",
    "comparator",
    "intervention_scope",
    "primary_endpoint",
    "endpoint_class",
    "effect_measure",
    "horizon_years",
    "model_coverage",
    "status",
    "rationale",
)

VALID_ENDPOINT_CLASSES = {
    "MACE",
    "HF_COMPOSITE",
    "FUNCTIONAL",
    "SYMPTOM",
    "SCENARIO",
    "NOT_DEFINED",
}
VALID_EFFECT_MEASURES = {
    "hazard_ratio",
    "risk_difference",
    "mean_difference",
    "absolute_risk_difference",
    "not_set",
}
VALID_MODEL_COVERAGE = {"scenario_supported", "proxy_only", "not_modelled", "scenario_only"}

DATA_INVENTORY_COLUMNS = (
    "dataset_id",
    "trial_id",
    "design",
    "population",
    "reported_sample_size",
    "ipd_plan_declared",
    "access_status",
    "overlap_with_current_model",
    "possible_validation_use",
    "not_licensed_to_validate",
    "source_url",
    "blocking_action",
)
VALID_IPD_STATUSES = {"yes", "no", "unknown", "no documented plan found"}


def _require_columns(frame: pd.DataFrame, required: tuple[str, ...], label: str) -> None:
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"{label} missing columns: {missing}")


def _valid_http_url(value: object) -> bool:
    parsed = urlparse(str(value).strip())
    return parsed.scheme in {"http", "https"} and bool(parsed.netloc)


def load_estimands(
    path: str | Path | None = None,
    population_ids: set[str] | None = None,
) -> pd.DataFrame:
    """Load planning estimands and fail on ambiguous or cross-reference-invalid rows."""
    path = Path(path) if path is not None else DEFAULT_DATA_DIR / "estimands.csv"
    if not path.is_file():
        raise FileNotFoundError(f"estimand registry not found: {path}")
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    _require_columns(frame, ESTIMAND_COLUMNS, "estimand registry")
    if frame.empty:
        raise ValueError("estimand registry is empty")
    if frame["estimand_id"].str.strip().eq("").any():
        raise ValueError("estimand_id cannot be blank")
    if frame["estimand_id"].duplicated().any():
        duplicates = sorted(frame.loc[frame["estimand_id"].duplicated(), "estimand_id"].tolist())
        raise ValueError(f"duplicate estimand_id: {duplicates}")

    populations = population_ids
    if populations is None:
        from ihd_solutions.population import load_populations

        populations = set(load_populations())

    for row in frame.to_dict(orient="records"):
        estimand_id = row["estimand_id"].strip()
        for required in (
            "disease_stage",
            "population_id",
            "comparator",
            "intervention_scope",
            "primary_endpoint",
            "endpoint_class",
            "effect_measure",
            "model_coverage",
            "status",
            "rationale",
        ):
            if not row[required].strip():
                raise ValueError(f"{estimand_id}: {required} cannot be blank")
        endpoint_class = row["endpoint_class"].strip()
        if endpoint_class not in VALID_ENDPOINT_CLASSES:
            raise ValueError(f"{estimand_id}: unsupported endpoint_class {endpoint_class!r}")
        measure = row["effect_measure"].strip()
        if measure not in VALID_EFFECT_MEASURES:
            raise ValueError(f"{estimand_id}: unsupported effect_measure {measure!r}")
        coverage = row["model_coverage"].strip()
        if coverage not in VALID_MODEL_COVERAGE:
            raise ValueError(f"{estimand_id}: unsupported model_coverage {coverage!r}")

        population_id = row["population_id"].strip()
        blocked = coverage == "not_modelled" and endpoint_class == "NOT_DEFINED"
        if population_id not in populations and not (population_id == "not_registered" and blocked):
            raise ValueError(f"{estimand_id}: unknown population_id {population_id!r}")
        horizon_text = row["horizon_years"].strip()
        if blocked:
            if horizon_text or measure != "not_set":
                raise ValueError(
                    f"{estimand_id}: blocked estimand must not invent horizon or effect measure"
                )
        else:
            try:
                horizon = float(horizon_text)
            except ValueError as exc:
                raise ValueError(f"{estimand_id}: horizon_years must be positive") from exc
            if not (horizon > 0.0):
                raise ValueError(f"{estimand_id}: horizon_years must be positive")

        if endpoint_class in {"MACE", "HF_COMPOSITE"} and measure not in {
            "hazard_ratio",
            "risk_difference",
        }:
            raise ValueError(
                f"{estimand_id}: hard outcomes require a risk or time-to-event measure"
            )
        if endpoint_class in {"FUNCTIONAL", "SYMPTOM"} and measure != "mean_difference":
            raise ValueError(f"{estimand_id}: functional/symptom outcome requires mean_difference")
        if endpoint_class == "NOT_DEFINED" and not blocked:
            raise ValueError(
                f"{estimand_id}: NOT_DEFINED is permitted only for explicitly blocked rows"
            )

    return frame


def load_validation_inventory(path: str | Path | None = None) -> pd.DataFrame:
    """Load candidate external datasets without implying that any dataset is accessed."""
    path = Path(path) if path is not None else DEFAULT_DATA_DIR / "validation_data_inventory.csv"
    if not path.is_file():
        raise FileNotFoundError(f"validation data inventory not found: {path}")
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    _require_columns(frame, DATA_INVENTORY_COLUMNS, "validation data inventory")
    if frame.empty:
        raise ValueError("validation data inventory is empty")
    if frame["dataset_id"].str.strip().eq("").any():
        raise ValueError("dataset_id cannot be blank")
    if frame["dataset_id"].duplicated().any():
        duplicates = sorted(frame.loc[frame["dataset_id"].duplicated(), "dataset_id"].tolist())
        raise ValueError(f"duplicate dataset_id: {duplicates}")

    for row in frame.to_dict(orient="records"):
        dataset_id = row["dataset_id"].strip()
        for required in (
            "trial_id",
            "design",
            "population",
            "reported_sample_size",
            "ipd_plan_declared",
            "access_status",
            "overlap_with_current_model",
            "possible_validation_use",
            "not_licensed_to_validate",
            "source_url",
            "blocking_action",
        ):
            if not row[required].strip():
                raise ValueError(f"{dataset_id}: {required} cannot be blank")
        if row["ipd_plan_declared"].strip().lower() not in VALID_IPD_STATUSES:
            raise ValueError(f"{dataset_id}: invalid ipd_plan_declared value")
        try:
            sample_size = int(row["reported_sample_size"])
        except ValueError as exc:
            raise ValueError(f"{dataset_id}: reported_sample_size must be an integer") from exc
        if sample_size <= 0:
            raise ValueError(f"{dataset_id}: reported_sample_size must be positive")
        if not _valid_http_url(row["source_url"]):
            raise ValueError(f"{dataset_id}: source_url must be an HTTP(S) URL")
    return frame


def audit_anchor_links(
    anchors: dict,
    catalog_ids: set[str],
    population_ids: set[str],
) -> dict[str, object]:
    """Check source/link completeness and registry cross-references, not source truth."""
    invalid_sources = sorted(
        anchor.anchor_id for anchor in anchors.values() if not _valid_http_url(anchor.source_url)
    )
    unknown_therapies = sorted({anchor.therapy_id for anchor in anchors.values()} - catalog_ids)
    unknown_populations = sorted(
        {anchor.population_id for anchor in anchors.values()} - population_ids
    )
    unknown_transport_targets = sorted(
        {
            target
            for anchor in anchors.values()
            for target in anchor.transportable_to
            if target not in population_ids
        }
    )
    if invalid_sources or unknown_therapies or unknown_populations or unknown_transport_targets:
        raise ValueError(
            "evidence registry link audit failed: "
            f"invalid_sources={invalid_sources}, unknown_therapies={unknown_therapies}, "
            f"unknown_populations={unknown_populations}, "
            f"unknown_transport_targets={unknown_transport_targets}"
        )
    return {
        "anchors_with_http_source": len(anchors),
        "therapy_links_valid": True,
        "population_links_valid": True,
        "transport_links_valid": True,
        "source_content_human_verified": False,
    }


def roadmap_input_audit(
    anchors: dict,
    catalog_ids: set[str],
    population_ids: set[str],
) -> dict[str, object]:
    """Validate structured inputs and report what remains unvalidated or inaccessible."""
    anchor_audit = audit_anchor_links(anchors, catalog_ids, population_ids)
    estimands = load_estimands(population_ids=population_ids)
    inventory = load_validation_inventory()
    return {
        "estimands_registered": len(estimands),
        "estimands_blocked_for_scope_or_population": int(
            estimands["model_coverage"].eq("not_modelled").sum()
        ),
        "external_dataset_leads": len(inventory),
        "external_datasets_accessed": 0,
        "external_validation_completed": False,
        "evidence_anchors_pending_source_review": len(anchors),
        **anchor_audit,
    }


def build_endgame_report(
    anchors: dict,
    catalog_ids: set[str],
    population_ids: set[str],
) -> str:
    """Render an honest current-state report; never claim inaccessible phases complete."""
    audit = roadmap_input_audit(anchors, catalog_ids, population_ids)
    estimands = load_estimands(population_ids=population_ids)
    inventory = load_validation_inventory()
    blocked = estimands.loc[estimands["model_coverage"] == "not_modelled"]

    lines = [
        "# IHD endgame execution report",
        "",
        "Generated from the repository registries. This is a computational status report, "
        "not a clinical, laboratory, or external-validation result.",
        "",
        "## Machine-checkable inputs",
        "",
        f"- Therapies: {len(catalog_ids)}; trial anchors: {len(anchors)}; "
        f"populations: {len(population_ids)}.",
        f"- Planning estimands: {audit['estimands_registered']}; explicitly not modelled: "
        f"{audit['estimands_blocked_for_scope_or_population']}.",
        f"- Candidate external-data leads: {audit['external_dataset_leads']}; datasets accessed: "
        f"{audit['external_datasets_accessed']}; held-out external validation completed: "
        f"{str(audit['external_validation_completed']).lower()}.",
        f"- Anchors awaiting human source review: "
        f"{audit['evidence_anchors_pending_source_review']}; no automated source verification is implied.",
        "- Source URL syntax and registry links pass; source contents human-verified: "
        f"{str(audit['source_content_human_verified']).lower()}.",
        "",
        "## Explicitly unmodelled planning scopes",
        "",
        "| Estimand | Population state | Endpoint state | Reason/status |",
        "|---|---|---|---|",
    ]
    for row in blocked.to_dict(orient="records"):
        lines.append(
            f"| `{row['estimand_id']}` | `{row['population_id']}` | "
            f"{row['endpoint_class']} | {row['status']} |"
        )

    lines.extend(
        [
            "",
            "## Candidate data leads (not acquired or approved)",
            "",
            "| Dataset | Design | Access status | Permitted near-term use / limitation |",
            "|---|---|---|---|",
        ]
    )
    for row in inventory.to_dict(orient="records"):
        limitation = row["not_licensed_to_validate"].replace("|", "/")
        lines.append(
            f"| {row['dataset_id']} ({row['trial_id']}) | {row['design']} | "
            f"{row['access_status']} | {limitation} |"
        )

    lines.extend(
        [
            "",
            "## Stop/go sequence",
            "",
            "1. **Proceed with registry/schema maintenance:** automated validation passes. "
            "Human source verification is still required before evidence-credit changes.",
            "2. **Do not claim external validation:** obtain data-use approval and an independent "
            "dataset first; preserve strict study-level holdout separation.",
            "3. **Do not start new efficacy optimization:** first reconcile sources, define the "
            "estimand, and validate/calibrate the relevant model.",
            "4. **Lab and clinical milestones are blocked outside this codebase:** product-specific "
            "preclinical safety/target-engagement work and randomized clinical evidence require "
            "facilities, approvals, participants, and independent oversight.",
            "5. **No cure or 80% target is established:** computational outputs remain model- and "
            "scenario-specific; unavailable observations cannot be replaced by simulations.",
            "",
            "The machine-generated source-review queue is written beside this report as "
            "`endgame_report_evidence_review_queue.csv`; all entries remain pending until a "
            "qualified human reviewer documents source checks.",
            "",
            "See `docs/ihd_endgame_roadmap.md` for phase details, dependencies, and acceptance gates.",
            "",
        ]
    )
    return "\n".join(lines)


def write_endgame_report(
    output_path: str | Path,
    anchors: dict,
    catalog_ids: set[str],
    population_ids: set[str],
) -> Path:
    """Write the current roadmap status report to a requested path."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(build_endgame_report(anchors, catalog_ids, population_ids), encoding="utf-8")
    return path


def build_evidence_review_queue(anchors: dict) -> pd.DataFrame:
    """Create a human-review work queue; URLs and fields are not verified by automation."""
    rows = []
    for anchor in sorted(anchors.values(), key=lambda item: item.anchor_id):
        rows.append(
            {
                "anchor_id": anchor.anchor_id,
                "therapy_id": anchor.therapy_id,
                "population_id": anchor.population_id,
                "endpoint_class": anchor.endpoint_class.value,
                "evidence_tier": anchor.evidence_tier.value,
                "source_url": anchor.source_url,
                "source_review_status": "pending_human_source_review",
                "reviewer": "",
                "reviewed_on_utc": "",
                "verify_design_population_endpoint_effect_ci_followup_and_source_location": "",
                "notes_and_discrepancies": "",
            }
        )
    return pd.DataFrame(rows)


def write_evidence_review_queue(output_path: str | Path, anchors: dict) -> Path:
    """Write a machine-generated human source-review queue, with every anchor pending."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    build_evidence_review_queue(anchors).to_csv(path, index=False)
    return path
