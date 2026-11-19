# IHD endgame roadmap: what is missing, what can be simulated, and what must be measured

## Purpose and honest boundary

This is the integrated plan for advancing the repository from a set of useful simulations to a validated decision-support and translational research programme. It consolidates the existing discovery, evidence, cure-search, gap-closure, unanswered-question, and revival analyses. It is a plan, **not a claim that IHD can be cured, that the present models predict individual patients, or that a simulated design is a treatment**.

“Once and for ever” is not a scientifically achievable guarantee for a heterogeneous, lifelong disease. The meaningful goal is durable reduction of preventable events and mortality, prevention of recurrent ischemia and disease progression, restoration of function where biologically possible, acceptable safety, and evidence that these benefits apply to the intended population. Primary prevention, stable coronary disease, acute MI, post-MI remodeling, and ischemic cardiomyopathy are related but distinct problems; they need distinct cohorts, mechanisms, and outcomes.

A simulation can rank questions, calculate requirements, design experiments, and reject implausible plans. It cannot create an effect, demonstrate safe delivery, or replace randomized human outcome evidence. In this repository, a “cure” still requires all five registered criteria, including randomized proof; the computational `discover` result is a capability specification, not a cure.

## What the repository currently says

Results below are model outputs with the status and limitations documented in their source files; none is a clinical prediction.

- The cure search reports **zero qualifying cures**. Its binding capability gaps are myocardial restoration (required +20 LVEF points in that specification versus +3.9 points observed in an uncontrolled BioVAT-HF cohort), substrate removal (0.80 required versus 0.431 demonstrated in the model's evidence ledger), and risk normalization (required HR ≤0.298 versus best evidence-graded value 0.360). Driver durability is the one capability the search considers solved in principle.
- The revival correction searched 256 composition-dose combinations. At the 3.4-point natural-history comparator its best modeled target probability is **57.81%** and risk-adjusted probability **54.62%**. Re-searching all 256 combinations at each documented nonzero comparator produced a best of **73.91% / 69.82%** at 1.1 points; none reaches 80%. The missing 3.37/4.34-point uniform shifts remain counterfactual thresholds, not treatment effects.
- The broad solution simulator separates prevention/event outcomes from LVEF, functional, biomarker, and heart-failure endpoints; its strict, surrogate, and exploratory modes must remain separate. A risk-reducing regimen is not automatically a myocardial-restoration therapy.
- The repo's own audits flag remaining evidence, identifiability, external-validation, model-structure, delivery, and lab/clinical gaps. Several are not closable in this workspace, regardless of compute budget.

Read together: the immediate computational task is not to raise a score. It is to establish which model outputs are calibrated against independent data, identify the smallest decisive unknowns, and translate those unknowns into safe, falsifiable experiments.

## Execution status on the current branch

This is a staged programme. The first computational prerequisites have now been started; downstream validation is explicitly blocked until the required data and experiments exist.

| Plan phase | Status now | Completed here | Still required / blocker |
|---|---|---|---|
| 0 — estimands | **Partly implemented** | `data/solutions/estimands.csv` defines nine planning estimands; `run_solutions.py estimands` lists them and `validate` checks population, endpoint class, measure, and horizon. Acute MI is deliberately marked blocked because no acute-MI population/protocol is registered. | Domain review is needed before these planning definitions become trial protocols; HF-composite and remodeling models remain uncalibrated. |
| 1 — evidence ledger | **Registry links audited; human source audit still open** | Validation checks source URL syntax and therapy/population/transport cross-references for all 48 anchor rows. It checks schema and references for the estimand and data-intake registries. `run_solutions.py endgame` exports a 48-row pending review queue at `results/endgame/endgame_report_evidence_review_queue.csv`. | URL syntax is not source verification. Each of the 48 sources still needs human review/version/date, study-design details, risk-of-bias and data-sharing metadata; no efficacy credit was upgraded. |
| 2 — external validation | **Blocked: no data accessed** | `data/solutions/validation_data_inventory.csv` records four candidate datasets and labels the possible use, overlap, access status, and what each cannot validate. | DREAM-HF lists an IPD-sharing plan, but actual release/criteria have not been verified; REVIVED's data-sharing supplement must be reviewed. No dataset has been requested, downloaded, approved, or used for held-out validation. |
| 3 — mechanistic model repair | **Partial / not externally validated** | Existing modules and audits remain available; evidence-gated uplift refuses to turn a mean treatment effect into uniform individual gain. | Integrated PK/PD, validated response heterogeneity, spatial mechanics, scar biology, and product-specific safety data are absent. |
| 4 — scenario search | **Specific subproblem done; global programme partial** | The revival search covers all 256 composition-dose pairs at three documented nonzero comparators. | This is not a fully calibrated optimizer across every disease stage, safety/cost/access objective, or validated interaction. |
| 5 — value of information | **Partial** | Existing discovery and gap analyses rank assumed capability uncertainties and calculate some sample-size/experiment requirements. | Priors, costs, and empirical variance inputs need independent review and real data; there is no approved experiment plan or budget. |
| 6 — lab/preclinical | **Not executable in this environment** | Requirements and kill criteria are specified in the roadmap and domain plans. | Product, tissue, animal/lab facilities, safety assays, ethics approvals, and independent replication are required. |
| 7 — clinical evidence | **Not simulatable** | Trial-design considerations and current evidence gaps are documented. | Sponsor, regulatory/ethics approval, trial sites, participants, randomization, outcome follow-up, and replication are indispensable. |
| 8 — governance | **Partly documented** | Provenance, correction ledgers, tests, and evidence tiers exist. | Continue independent source review, data-use governance, conflict disclosures, and versioned revalidation as new evidence appears. |

The data inventory is a **lead register**, not data access. ClinicalTrials.gov's DREAM-HF record reports completed phase 3 randomization, posted results, and an IPD-sharing plan marked “Yes”; this is a plan, not evidence that files are released or a request will be approved. There is a denominator discrepancy that must be reconciled before requesting data: the current v2 API response returned 565 actual enrollment, study-record version 39 lists 566, and the JACC paper reports n=565. Confirm screening, randomized, treated, intention-to-treat, and endpoint-available denominators from the source flow diagram and custodian. The current record/API should be rechecked before any application. REVIVED has already supplied aggregate inputs to the revival correction, so it cannot be a held-out validation set for that specific component; its data-sharing statement and any later investigator-request route still require formal review. The single-arm BioVAT-HF and AATS LB8 cohorts can inform feasibility/safety or exploratory response characterization, but cannot provide a randomized treatment contrast. Data-lead desk review was refreshed on 2026-11-06 for DREAM-HF and REVIVED-BCIS2; see the exact source URLs in `data/solutions/validation_data_inventory.csv`. This is not the pending human source review of the 48 evidence anchors.

The integrated status can be regenerated with `python run_solutions.py endgame`; the current snapshot is `results/endgame/endgame_report.md`, accompanied by a 48-anchor human-review queue. That report only certifies schema and cross-reference checks, not source truth or biological validity. Draft, unsent data-access requests and their mandatory governance checklist are in `docs/data_access_request_templates.md`.

## What is still missing

| Missing item | Why it matters | Can this workspace settle it? | Evidence or deliverable needed |
|---|---|---|---|
| Operational disease scope and success criteria | “IHD cure” conflates prevention, plaque stabilization, acute reperfusion, scar repair, and heart failure | Partly | Freeze separate estimands, populations, time horizons, and hard outcomes before each analysis |
| External validation data | Agreement with fitted inputs is not validation; held-out trials test transport and calibration | Only if accessible public/permissioned data are supplied | Prespecified study-level holdout or independently curated participant-level data; report calibration and error |
| Patient-level heterogeneity and response distribution | Trial mean effects do not imply uniform individual gains; subgroup claims can be unstable | Not from aggregate means alone | IPD or a validated predictive response model, with transportability checks and uncertainty |
| Causal treatment effect in each target population | Biomarker, surrogate, open-label, and uncontrolled results cannot establish hard-outcome benefit | No, except by synthesis of suitable trials | Randomized comparator evidence on relevant populations and endpoints |
| Plaque natural history and spatial mechanics | Scalar compartments may miss lesion location, geometry, cap mechanics, thrombosis, and competing lesions | Some equations are computable; parameter truth is not | Human imaging/pathology longitudinal data; validated CFD/FSI/finite-element inputs |
| Myocardial recovery and scar biology | Viability, transmurality, perfusion, remodeling, and cell integration determine who can recover | Not fully | Segment-level CMR/echo data, controlled treatment contrasts, preclinical tissue integration and safety data |
| Delivery, PK/PD, durability, and target engagement | A potent payload that does not reach the target is not an effective intervention | Partly in models; actual delivery is experimental | Biodistribution, dose/exposure, target-engagement, durability, and withdrawal data |
| Combination interactions and cumulative harms | Independent-effect multiplication can overstate efficacy and understate toxicity | Exploratory simulations only | Factorial or otherwise identifiable combination data; adverse events and competing risks by exposure |
| Long-term safety | Off-target editing, immune effects, arrhythmia, tumors, thrombosis, and repeat-dose accumulation cannot be inferred from efficacy alone | No | Product-specific toxicology, monitoring, follow-up, and independent adjudication |
| Clinical benefit and “cure” status | A surrogate change or model probability is not reduced death/MI/stroke or durable restoration | No | Adequately powered randomized trials, follow-up, replication, and post-market surveillance |
| Equity, feasibility, cost, and access | A technically effective approach may be inaccessible, non-scalable, or harmful at population level | Partly | Resource, implementation, ethics, and distributional-impact analyses with real data |

## Work plan

### Phase 0 — Freeze scope and estimands (repository work)

**Question:** Which problem is being solved in each run?

1. Maintain separate task definitions for: primary prevention; established ASCVD; acute MI/reperfusion; post-MI remodeling; ischemic cardiomyopathy/HFrEF; and recurrent events after intervention.
2. For each task, prespecify population, inclusion/exclusion, baseline risk, comparator, intervention, follow-up, primary endpoint, safety endpoints, and estimand (intention-to-treat unless justified otherwise).
3. Preserve the five-criterion cure specification, but report which criteria are computable and which require trials. Do not let an LVEF threshold silently become a death/MI/stroke endpoint.

**Deliverable:** versioned `data/solutions/estimands.csv` (or equivalent schema), a checked registry, and a test that each run names its estimand.

**Exit gate:** all registered simulations point to a population and endpoint; cure proof stays false until randomized human evidence meets its definition.

### Phase 1 — Evidence ledger and provenance audit (repository work)

1. Reconcile `trial_anchors.csv`, candidate catalogues, cure-search tables, revival corrections, and literature register into one source-of-truth evidence graph.
2. Store DOI/registry ID, exact source URL, trial design, randomization/blinding, sample size, population, comparator, endpoint class, effect estimate, CI, follow-up, missingness, funding, prespecification, subgroup status, and evidence tier.
3. Distinguish direct RCT effects, observational associations, surrogate endpoints, mechanistic measurements, and hypothetical assumptions. Record null and harmful findings as carefully as positive findings.
4. Add a source-review date and reproducible provenance. For future-looking or changing trial statuses, require re-checking the registry and publication before updating credit.

**Deliverable:** validated evidence catalogue and machine-readable evidence-to-credit decisions; no unsourced efficacy values.

**Exit gate:** malformed/missing provenance fails validation; neutral, pending, surrogate, and exploratory anchors cannot leak into strict hard-outcome credit.

### Phase 2 — External validation and calibration (requires independent data)

1. Inventory public, licensed, or collaborator-provided data without assuming access to individual participant records.
2. Pre-register a study-level split: never randomly split records from the same trial across train and test. Reserve independent trials/populations for validation.
3. Compare predicted versus observed treatment effects for hard outcomes separately from LVEF, biomarkers, and functional outcomes. Report calibration-in-the-large, calibration slope, discrimination where appropriate, confidence intervals, and prediction error.
4. Test population transport explicitly (age, sex, diabetes, renal disease, prior MI, scar burden, baseline therapy, geography); flag unsupported transport rather than filling it with a generic average.
5. If no valid external dataset is available, label validation `not assessed`; do not substitute repeated simulation or re-fitting to the same anchors.

**Deliverable:** held-out validation report, input data dictionary, data-use note, and reproducible scripts; no patient data committed unless permission and governance allow it.

**Exit gate:** prespecified performance thresholds met on held-out trials, or model is restricted to descriptive scenario analysis and its poor transport is reported.

### Phase 3 — Mechanistic and pharmacology model repair (some compute; data-dependent)

1. Build separate, coupled but auditable modules for lipid/atheroma progression, plaque cap and thrombosis, ischemia/reperfusion injury, scar/remodeling, and pump function. Specify units, timescales, compartments, and conservation/bound constraints.
2. Add PK/PD, delivery probability, adherence, target engagement, waning, withdrawal, repeat dosing, and dose-dependent harms where supporting data exist. Keep unavailable quantities as explicit uncertainty, not fitted constants.
3. Model patient heterogeneity using evidence-supported covariates. Never map a cohort-average LVEF effect to every individual. The evidence gate should require an externally validated response distribution before patient-level uplift is consumed.
4. Treat combinations conservatively: no multiplication of independent effects unless the causal assumptions are stated and supported; require factorial or otherwise identifiable interaction evidence for synergy. Always model cumulative adverse events and competing risks.
5. Separate epistemic uncertainty (unknown parameters/model form) from aleatoric variation (patient-to-patient variation); publish both, and do not conflate them with Monte Carlo numerical error.

**Deliverable:** versioned model equations, unit tests, parameter provenance, identifiability audit, posterior/predictive checks only where data justify them, and adversarial model-form sensitivity.

**Exit gate:** no decision-driving parameter is hidden, unidentifiable products are exposed, and structural alternatives do not reverse conclusions without the reversal being reported.

### Phase 4 — Reproducible scenario search and robust decision analysis (computable now, evidence constrained)

1. Exhaustively enumerate registered therapies, dose ranges supported by precedent, timing, adherence scenarios, target populations, and comparator assumptions. Keep simulation-only constructs in a separate namespace.
2. Use fixed seeds for reproducibility and multiple seeds for Monte Carlo error checks; calculate confidence/credible intervals for probabilities, not only point estimates.
3. Optimize a prespecified multi-objective: hard-outcome benefit, safety, burden of treatment, access/cost, and robustness. Show the Pareto frontier rather than inventing one scalar “best” score.
4. Run global sensitivity, ablation, threshold, comparator, and worst-case analyses. Correct for selection/winner's curse with held-out simulation seeds or nested selection; never select a comparator after seeing the best score.
5. Include negative controls: zero effect, neutral-trial effect, no-adherence, no-delivery, sham, and code paths where all uncredited components yield zero benefit.
6. Retain explicit labels for `observed`, `transported`, `assumed`, `hypothetical`, and `unvalidated`. A counterfactual may show what would be required but must not change an evidence-backed ranking.

**Deliverable:** command such as `python run_solutions.py endgame` that emits estimand-specific scenario tables, robust frontiers, uncertainty decomposition, and a summary showing which constraints fail.

**Exit gate:** deterministic reproduction, no optimistic credit leakage, all best designs reported with uncertainty and assumptions; if the feasible set is empty, report it as empty.

### Phase 5 — Value-of-information and experiment prioritization (computable with explicit priors)

1. Rank uncertainties by expected change in the decision, not merely parameter variance. Compare the value of measuring delivery, plaque clearance, myocardial recovery, and adverse-event incidence.
2. For each top uncertainty, specify the smallest experiment that discriminates among competing models, expected information, costs, and kill criteria.
3. Sweep priors and costs; if the recommended experiment changes over plausible values, report a conditional recommendation instead of a single answer.
4. Do not claim a simulation has selected the exact laboratory or clinical effect; it selects measurements that would reduce decision uncertainty.

**Deliverable:** experiment-priority ledger with assumptions, minimum detectable effects, sample-size/power calculations, and explicit links to decision-changing uncertainties.

**Exit gate:** every proposed experiment has a falsifiable hypothesis, comparator, endpoint, uncertainty target, and stop/continue rule.

### Phase 6 — Laboratory and preclinical validation (cannot be completed by this repository)

For any regenerative or delivery construct, require sequential evidence for product identity/purity, potency, target engagement, biodistribution, retention, dose response, off-target activity, immune compatibility, tumorigenicity, arrhythmogenicity, thrombosis, and durability. Use appropriate human cells/tissue and disease-relevant models; animal results support mechanism and safety but are not human efficacy.

For plaque-clearing designs, measure human plaque explant penetration, lipid/crystal handling, macrophage survival/efferocytosis, efflux, cap integrity, and pro-thrombotic effects. For gene editing, include editor-specific off-target and structural-variant assays, liver toxicity, immunogenicity, durability, and reversibility/mitigation strategy.

**Deliverable:** preregistered experimental protocol, product characterization, blinded readouts, independent replication, and a safety dossier. No passing simulation may bypass these gates.

**Exit gate:** predefined safety and mechanism thresholds met; otherwise stop, redesign, or narrow the indication.

### Phase 7 — Clinical evidence generation (cannot be simulated into existence)

1. Select the indication where the intervention's mechanism and patient biology match; define the control as contemporary standard of care, not an obsolete baseline.
2. Use randomized, blinded where feasible, independently adjudicated trials. Power the primary hard outcome or a clinically justified endpoint; declare LVEF/biomarkers secondary unless their surrogate validity is established for that population.
3. Prespecify subgroup analyses; treat exploratory subgroup findings as hypothesis-generating until independently confirmed. Avoid post-hoc enrichment claims.
4. Measure patient-level safety and benefit, including cardiovascular death, MI, stroke, HF hospitalization, arrhythmia, serious adverse events, quality of life, and function as appropriate. Follow long enough to measure durability and late harms.
5. For combination regimens, establish component and interaction effects; do not infer that each promising monotherapy will add safely or synergistically.
6. Replicate pivotal findings and continue post-market pharmacovigilance. Only then update evidence tiers and any cure criterion.

**Deliverable:** trial-ready target product profile and protocol synopsis generated from validated evidence and expert review; the repository's role is decision support, not trial authorization.

**Exit gate:** clinically meaningful, safe, durable benefit versus current standard of care in the target population, independently replicated. A “cure” claim additionally requires all five registered criteria, not only a surrogate or modeled probability.

### Phase 8 — Governance, update, and maintenance

1. Add versioned evidence review, independent review of high-impact claims, and a visible correction/retraction ledger.
2. Record model, data, code, seed, and dependency versions for each result. Pin output summaries to the exact commit.
3. Protect privacy, consent, data-use terms, equity, access, and conflict-of-interest disclosures.
4. Re-run validation when a new pivotal trial, safety signal, guideline, or major model correction changes a decision.
5. Maintain an explicit `not solved` register and do not let age of a result convert `unvalidated` into `validated`.

## Simulation implementation map for this repository

| Sequence | Repository action | Existing starting point | Acceptance check |
|---|---|---|---|
| 1 | Audit and align registries; create estimand and evidence schemas | `data/solutions/`, `src/ihd_solutions/evidence.py`, `docs/literature_register.md` | validator rejects missing sources, population, comparator, endpoint, or tier |
| 2 | Reconcile cure/gap metrics with current outputs | `discovery.py`, `cure.py`, `gap_solver.py`, CLI `unanswered` command | same requirements/provenance appear in one machine-readable gap ledger |
| 3 | Add study-level holdout validation harness | `pipeline.py`, `validation.py`, `trial_anchors.csv` | no trial leakage; outputs calibration and uncertainty with data access status |
| 4 | Connect uncertainty and sensitivity modules to scenario orchestration | `uncertainty.py`, `simulation_frontier.py`, `revival.py` | every run exports seeds, draw counts, assumptions, intervals, and Pareto/frontier tables |
| 5 | Add PK/PD, delivery, spatial plaque, and scar-response models only as data permits | `pharmacology.py`, `plaque.py`, `regeneration.py`, legacy structural models | identifiability plus held-out validation precedes efficacy credit |
| 6 | Add experiment/sample-size decision outputs | `discovery.py`, `gap_solutions.py`, `docs/gap_closure_attempt.md`, CLI `unanswered` command | each experiment targets a ranked unknown and states a falsification boundary |
| 7 | Generate a consolidated `endgame` report/CLI command | `cli.py`, `roadmap.py`, current registries | **Implemented (initial)**: `python run_solutions.py endgame`; emits explicit data-access, source-review, scope, and lab/clinical blockers without mixing evidence tiers |
| 8 | Create controlled lab/trial handoff, not simulated efficacy | `docs/`, output schemas | **Not complete**: roadmap checklists exist, but a product-specific protocol package requires independent scientific, safety, ethics, and trial-expert review |

Do these in sequence. Do not start Phase 5 optimization before Phase 1 evidence reconciliation and Phase 2 validation; otherwise the optimizer will confidently search an uncalibrated model.

## First 90-day computational sprint (conditional on data availability)

1. **Weeks 1–2:** freeze estimands and reconcile evidence/source conflicts; publish a versioned gap ledger.
2. **Weeks 3–5:** identify external trial datasets and permissions; implement the study-level holdout harness and run only what data permit.
3. **Weeks 6–8:** repair model interfaces and uncertainty outputs; run negative controls, identifiability, and adversarial sensitivity.
4. **Weeks 9–10:** exhaustive robust scenario and value-of-information search; rank the next measurements, not hypothetical cures.
5. **Weeks 11–12:** issue a reproducible endgame report and a lab/clinical handoff with stop/go criteria; independent scientific review before any claim upgrade.

If external data are unavailable, Weeks 3–5 produce a documented data-acquisition request rather than a fake validation score. Wet-lab and clinical phases are separate, resourced programmes, not deliverables promised by this sprint.

## Non-negotiable success rules

- Never change the requested target, comparator, population, or endpoint after seeing results; changes require a versioned rationale and both old/new results.
- Never label a simulation probability as clinical efficacy, a cure, or individual prognosis.
- Never transfer an effect between products, populations, endpoints, or dose regimens without a stated transportability argument.
- Never apply average treatment effects uniformly to simulated individuals without a validated response distribution.
- Never claim synergy from stacking component effects; show evidence and interactions.
- Never hide adverse events, competing risks, null findings, Monte Carlo error, selection bias, or uncertainty.
- If no evidence-backed combination crosses a threshold, the correct finding is **not reached**, plus the measurements/trials required to test it.

## Existing detailed analyses to preserve and extend

- `docs/cure_discovery_plan.md` — feasibility frontier, capability gaps, and discovery stages.
- `docs/cure_search_2026.md` — candidate search and evidence anchors.
- `docs/solutions_simulation.md` — intervention/evidence modes and code map.
- `docs/unanswered_questions.md` — computation, resource, analytic, and epistemic gaps.
- `docs/revival_correction.md` — control-arm correction, 256-pair search, comparator frontier, and fail-closed evidence gate.
- `docs/scientific_methodology_and_limitations.md` and `docs/scientific_audit.md` — epistemic limits, provenance, and known corrections.

This roadmap is the umbrella. It does not supersede more specific domain plans, and later work should update those source documents rather than duplicate their changing results.
