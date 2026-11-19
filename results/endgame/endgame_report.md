# IHD endgame execution report

Generated from the repository registries. This is a computational status report, not a clinical, laboratory, or external-validation result.

## Machine-checkable inputs

- Therapies: 47; trial anchors: 48; populations: 19.
- Planning estimands: 9; explicitly not modelled: 4.
- Candidate external-data leads: 4; datasets accessed: 0; held-out external validation completed: false.
- Anchors awaiting human source review: 48; no automated source verification is implied.
- Source URL syntax and registry links pass; source contents human-verified: false.

## Explicitly unmodelled planning scopes

| Estimand | Population state | Endpoint state | Reason/status |
|---|---|---|---|
| `acute_mi_reperfusion` | `not_registered` | NOT_DEFINED | blocked_missing_population_and_protocol |
| `post_mi_remodeling` | `post_mi_pci` | FUNCTIONAL | planning_only_not_trial_protocol |
| `advanced_ischemic_hf` | `advanced_ischemic_hf` | HF_COMPOSITE | planning_only_requires_expert_review |
| `refractory_angina` | `refractory_angina_no_option` | SYMPTOM | planning_only_requires_expert_review |

## Candidate data leads (not acquired or approved)

| Dataset | Design | Access status | Permitted near-term use / limitation | Lead/source review note |
|---|---|---|---|---|
| DREAM_HF (NCT02032004) | Phase 3 randomized double-blind sham-controlled | access route and data release not verified | not direct efficacy validation for the modeled iPSC patch or its effects | Denominator discrepancy to resolve before requesting data: current ClinicalTrials.gov v2 response shows 565 actual enrollment, while record version 39 (https://clinicaltrials.gov/ct2/history/NCT02032004?V_39=View) shows 566; the JACC report describes n=565 (https://www.jacc.org/doi/10.1016/j.jacc.2022.11.061). Confirm consented, randomized, treated, ITT, and endpoint-available counts from the flow diagram and sponsor. |
| REVIVED_BCIS2 (NCT01920048) | Randomized open-label PCI plus OMT versus OMT | publication supplement data-sharing route must be reviewed | not external validation of the revival comparator correction because its aggregate result already informed that model | Aggregate REVIVED findings already informed one model component. A 2025 IPD analysis states data may be available on reasonable request (https://pmc.ncbi.nlm.nih.gov/articles/PMC12146480/); verify its terms, scope, custodians, and overlap before treating any analysis as independent validation. |
| BioVAT_HF (NCT04396899) | Phase 1/2 single-arm dose-finding and feasibility | IPD availability not verified | no randomized treatment-effect validation or causal estimate from a single arm | No current participant-level data availability or approval verified. Single-arm design cannot estimate a randomized causal treatment effect. |
| AATS_LB8 (LB8) | Prospective single-arm open-label multicenter feasibility cohort | individual data access not verified | no comparative treatment effect or efficacy validation | No current participant-level data availability or approval verified. Small single-arm feasibility cohort supports descriptive feasibility/safety only, not comparative efficacy. |

## Stop/go sequence

1. **Proceed with registry/schema maintenance:** automated validation passes. Human source verification is still required before evidence-credit changes.
2. **Do not claim external validation:** obtain data-use approval and an independent dataset first; preserve strict study-level holdout separation.
3. **Do not start new efficacy optimization:** first reconcile sources, define the estimand, and validate/calibrate the relevant model.
4. **Lab and clinical milestones are blocked outside this codebase:** product-specific preclinical safety/target-engagement work and randomized clinical evidence require facilities, approvals, participants, and independent oversight.
5. **No cure or 80% target is established:** computational outputs remain model- and scenario-specific; unavailable observations cannot be replaced by simulations.

The machine-generated source-review queue is written beside this report as `endgame_report_evidence_review_queue.csv`; all entries remain pending until a qualified human reviewer documents source checks.

See `docs/ihd_endgame_roadmap.md` for phase details, dependencies, and acceptance gates.
