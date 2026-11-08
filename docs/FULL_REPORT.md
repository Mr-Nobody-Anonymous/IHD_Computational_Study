# FULL REPORT — IHD Computational Study
## Every agent, every compound, every experiment, every result

**Repository:** `Mr-Nobody-Anonymous/IHD_Computational_Study`
**Branch:** `arena/2330d67a-ihd-computational-study`
**Report date:** 2026-11-07

---

## How to read this report

This study contains five *different kinds* of work, and mixing them would be
misleading. Every section is labelled with its epistemic status:

| Label | Meaning | Can it support a clinical claim? |
|---|---|---|
| **EMPIRICAL** | Measured data from published wet-lab work, or deterministic geometry computed from deposited crystal structures | Yes, within the limits of the original measurement |
| **EVIDENCE-GRADED SIMULATION** | Scenario model anchored to published trial effect estimates, with credit rules | No — scenario only, not prediction |
| **EXPLORATORY MODEL** | Textbook dynamical systems with assumed parameters | No — illustrative |
| **INVENTED / UNVALIDATED** | Constructs designed in this study that do not exist | No — explicitly quarantined |
| **NON-USABLE** | Legacy output that fails internal quality checks | **No — flagged, not interpreted** |

**Headline verdict:** no cure was found, and none is claimed. The best candidate is
**ATLAS-2**, an invented construct that makes plaque clearance robust *in
simulation* at a stated price. Randomised proof is zero for every combination in
this study, because no combination has been trialled.

---

# PART 0 — INVENTORY

## 0.1 What is in the repository

| results tree | files |
|---|---|
| `results/atlas/` | 9 |
| `results/atlas3/` | 10 |
| `results/atlas4/` | 11 |
| `results/audit/` | 7 |
| `results/composite/` | 8 |
| `results/discovery/` | 10 |
| `results/evidence_audit/` | 6 |
| `results/figures/` | 3 |
| `results/gap_solutions/` | 10 |
| `results/gaps/` | 9 |
| `results/hypotheses/` | 10 |
| `results/recombine/` | 13 |
| `results/reports/` | 5 |
| `results/resolution/` | 9 |
| `results/solutions/` | 34 |
| `results/tables/` | 10 |
| `results/upgrade/` | 10 |

**Input registries:** 47 therapies · 48 trial anchors · 19 populations
**Documentation:** 19 documents (plus this report)
**Code:** 33 modules in `src/ihd_solutions/`
**Tests:** 32 test files

---

# PART A — THE THERAPEUTIC AGENTS (all 47)

**Status: EVIDENCE-GRADED.** Each row is a real agent with a published mechanism.
Effect sizes are percentage reductions from published trials or meta-analyses;
they are inputs to the scenario model, not predictions.

## A.1 Counts by class, evidence tier and regulatory status

**By therapeutic class**

| class | agents |
|---|---|
| LIPID_LOWERING | 13 |
| ANTI_INFLAMMATORY | 8 |
| LPA_LOWERING | 5 |
| METABOLIC | 5 |
| REGENERATIVE | 5 |
| LIFESTYLE | 3 |
| REVASCULARIZATION | 3 |
| GENE_EDITING | 2 |
| ANTITHROMBOTIC | 2 |
| BLOOD_PRESSURE | 1 |

**By evidence tier**

| evidence tier | agents |
|---|---|
| RANDOMIZED_SURROGATE | 19 |
| RANDOMIZED_OUTCOMES | 15 |
| META_ANALYSIS | 5 |
| MECHANISTIC | 4 |
| OBSERVATIONAL | 2 |
| RANDOMIZED_TOPLINE | 1 |
| CONDITIONAL_APPROVAL | 1 |

**By regulatory status**

| status | agents |
|---|---|
| approved | 18 |
| phase_3 | 6 |
| phase_2 | 4 |
| recommended | 3 |
| phase_1 | 2 |
| investigational | 2 |
| preclinical | 2 |
| not_approved_for_ascvd | 1 |
| not_indicated | 1 |
| conditional_japan | 1 |
| phase_3_outcomes_pending | 1 |
| phase_1b | 1 |
| approved_sHTG | 1 |
| filed_sHTG | 1 |
| phase_2a | 1 |
| post_hoc_hypothesis | 1 |
| approved_limited | 1 |

## A.2 Full agent register — mechanism of action

Every agent, its class, and what it actually does molecularly.

| therapy_id | name | therapy_class | mechanism |
|---|---|---|---|
| statin_high_intensity | High-intensity statin | LIPID_LOWERING | HMG-CoA reductase inhibition |
| statin_moderate_intensity | Moderate-intensity statin | LIPID_LOWERING | HMG-CoA reductase inhibition |
| ezetimibe | Ezetimibe | LIPID_LOWERING | NPC1L1 cholesterol-absorption inhibition |
| bempedoic_acid | Bempedoic acid | LIPID_LOWERING | ATP-citrate lyase inhibition |
| evolocumab | Evolocumab | LIPID_LOWERING | Anti-PCSK9 monoclonal antibody |
| alirocumab | Alirocumab | LIPID_LOWERING | Anti-PCSK9 monoclonal antibody |
| inclisiran | Inclisiran | LIPID_LOWERING | siRNA silencing of hepatic PCSK9 |
| enlicitide_oral | Enlicitide oral macrocyclic peptide | LIPID_LOWERING | Oral macrocyclic peptide blocking PCSK9-LDLR binding |
| azd0780_oral | AZD0780 oral small molecule | LIPID_LOWERING | Small-molecule PCSK9 inhibition |
| verve102_base_editing | VERVE-102 PCSK9 base editing | GENE_EDITING | Single-infusion adenine base editing of PCSK9 splice site |
| pelacarsen | Pelacarsen | LPA_LOWERING | Antisense oligonucleotide against apolipoprotein(a) |
| olpasiran | Olpasiran | LPA_LOWERING | siRNA against apolipoprotein(a) |
| muvalaplin_oral | Muvalaplin | LPA_LOWERING | Oral small-molecule disruptor of apo(a)-apoB assembly |
| colchicine | Low-dose colchicine | ANTI_INFLAMMATORY | Tubulin and NLRP3-mediated inflammatory dampening |
| canakinumab | Canakinumab | ANTI_INFLAMMATORY | IL-1beta monoclonal antibody |
| ziltivekimab | Ziltivekimab | ANTI_INFLAMMATORY | IL-6 ligand monoclonal antibody |
| methotrexate_low_dose | Low-dose methotrexate | ANTI_INFLAMMATORY | Broad anti-inflammatory without IL-6 pathway engagement |
| nlrp3_inhibitor | Selective NLRP3 inhibitor | ANTI_INFLAMMATORY | Direct NLRP3 inflammasome inhibition |
| semaglutide_2_4mg | Semaglutide 2.4 mg weekly | METABOLIC | GLP-1 receptor agonism with weight and metabolic effects |
| empagliflozin | Empagliflozin | METABOLIC | SGLT2 inhibition |
| sglt2_inhibitor_hf | SGLT2 inhibitor for heart failure | METABOLIC | SGLT2 inhibition in heart failure |
| bp_intensification | Blood-pressure intensification | BLOOD_PRESSURE | Additional antihypertensive therapy to lower systolic pressure |
| aspirin_secondary_prevention | Aspirin secondary prevention | ANTITHROMBOTIC | Platelet cyclo-oxygenase inhibition |
| icosapent_ethyl | Icosapent ethyl | ANTITHROMBOTIC | High-dose EPA formulation |
| smoking_cessation | Smoking cessation | LIFESTYLE | Removal of tobacco exposure |
| exercise_cardiac_rehab | Exercise-based cardiac rehabilitation | LIFESTYLE | Structured exercise and risk-factor programme |
| diet_ldl_lowering | LDL-lowering dietary pattern | LIFESTYLE | Saturated-fat replacement and plant-sterol intake |
| pci_stable_cad | PCI for stable coronary disease | REVASCULARIZATION | Percutaneous revascularisation |
| cabg_ischemic_cm | CABG for ischaemic cardiomyopathy | REVASCULARIZATION | Surgical revascularisation |
| ipsc_cm_sheet | iPSC-derived cardiomyocyte sheet | REGENERATIVE | Paracrine epicardial cell-sheet therapy |
| obicetrapib | Obicetrapib 10 mg | LIPID_LOWERING | Selective CETP inhibition |
| ctx310_angptl3_editing | CTX310 ANGPTL3 CRISPR editing | GENE_EDITING | One-time in vivo Cas9 loss-of-function editing of hepatic ANGPTL3 |
| lepodisiran | Lepodisiran | LPA_LOWERING | Ultra-long-acting GalNAc-siRNA against LPA |
| zerlasiran | Zerlasiran | LPA_LOWERING | GalNAc-siRNA against LPA |
| plozasiran | Plozasiran (APOC3 siRNA) | LIPID_LOWERING | RNAi suppression of apolipoprotein C-III |
| olezarsen | Olezarsen (APOC3 ASO) | LIPID_LOWERING | Antisense suppression of apolipoprotein C-III |
| tirzepatide | Tirzepatide | METABOLIC | Dual GIP/GLP-1 receptor agonism |
| retatrutide | Retatrutide | METABOLIC | Triple GIP/GLP-1/glucagon receptor agonism |
| dfv890_nlrp3 | DFV890 (selective NLRP3 inhibitor) | ANTI_INFLAMMATORY | Direct selective NLRP3 inflammasome inhibition |
| canakinumab_tet2_guided | Genotype-guided canakinumab (TET2-CHIP) | ANTI_INFLAMMATORY | IL-1beta blockade restricted to TET2 clonal-haematopoiesis carriers |
| yap101_sav1_gene_therapy | YAP101 (AAV9-SAV1-shRNA) | REGENERATIVE | Cardiomyocyte-restricted Hippo pathway release to trigger proliferation |
| biovat_hf_patch | BioVAT engineered heart muscle patch | REGENERATIVE | Epicardial iPSC-derived cardiomyocyte/stromal tissue for remuscularisation |
| cd34_cell_therapy | Autologous CD34+ cell therapy | REGENERATIVE | Pro-angiogenic cell delivery to ischaemic myocardium/microvasculature |
| coronary_sinus_reducer | Coronary sinus reducer | REVASCULARIZATION | Controlled coronary sinus narrowing to redistribute perfusion |
| pcsk9_vaccine | PCSK9 peptide vaccine (AT04A-class) | LIPID_LOWERING | Active immunisation raising anti-PCSK9 antibodies |
| car_treg_oxldl | CAR-Treg against oxidised LDL | ANTI_INFLAMMATORY | Antigen-specific regulatory T cells restoring tolerance to oxLDL |
| ab1002_gene_therapy | AB-1002 / umiposgene parvec | REGENERATIVE | Intracoronary AAV delivery of constitutively active inhibitor-1 (I-1c) |

## A.3 Full agent register — measured effects

Fractional reductions (0.500 = 50% reduction). `sbp` is mmHg, not a fraction.
These are the model's pharmacodynamic inputs.

| therapy_id | ldl_pct | apob_pct | lpa_pct | hscrp_pct | sbp_reduction_mmhg | burden_score |
|---|---|---|---|---|---|---|
| statin_high_intensity | 0.5 | 0.42 | 0 | 0.15 | 0 | 1 |
| statin_moderate_intensity | 0.3 | 0.25 | 0 | 0.1 | 0 | 1 |
| ezetimibe | 0.18 | 0.15 | 0 | 0 | 0 | 1 |
| bempedoic_acid | 0.21 | 0.17 | 0 | 0 | 0 | 1 |
| evolocumab | 0.59 | 0.5 | 0.2 | 0 | 0 | 3 |
| alirocumab | 0.55 | 0.46 | 0.15 | 0 | 0 | 3 |
| inclisiran | 0.5 | 0.42 | 0.2 | 0 | 0 | 2 |
| enlicitide_oral | 0.558 | 0.503 | 0.282 | 0 | 0 | 1 |
| azd0780_oral | 0.507 | 0.45 | 0 | 0 | 0 | 1 |
| verve102_base_editing | 0.62 | 0.52 | 0 | 0 | 0 | 8 |
| pelacarsen | 0 | 0.1 | 0.8 | 0 | 0 | 3 |
| olpasiran | 0 | 0.1 | 0.95 | 0 | 0 | 2 |
| muvalaplin_oral | 0 | 0.05 | 0.65 | 0 | 0 | 1 |
| colchicine | 0 | 0 | 0 | 0.25 | 0 | 1 |
| canakinumab | 0 | 0 | 0 | 0.37 | 0 | 3 |
| ziltivekimab | 0 | 0 | 0 | 0.6 | 0 | 3 |
| methotrexate_low_dose | 0 | 0 | 0 | 0 | 0 | 1 |
| nlrp3_inhibitor | 0 | 0 | 0 | 0.4 | 0 | 1 |
| semaglutide_2_4mg | 0 | 0.05 | 0 | 0.35 | 5 | 2 |
| empagliflozin | 0 | 0 | 0 | 0 | 4 | 1 |
| sglt2_inhibitor_hf | 0 | 0 | 0 | 0 | 3 | 1 |
| bp_intensification | 0 | 0 | 0 | 0 | 10 | 1 |
| aspirin_secondary_prevention | 0 | 0 | 0 | 0 | 0 | 1 |
| icosapent_ethyl | 0 | 0.05 | 0 | 0.2 | 0 | 1 |
| smoking_cessation | 0 | 0 | 0 | 0.3 | 5 | 2 |
| exercise_cardiac_rehab | 0.03 | 0.02 | 0 | 0.2 | 4 | 2 |
| diet_ldl_lowering | 0.15 | 0.12 | 0 | 0.1 | 3 | 2 |
| pci_stable_cad | 0 | 0 | 0 | 0 | 0 | 9 |
| cabg_ischemic_cm | 0 | 0 | 0 | 0 | 0 | 10 |
| ipsc_cm_sheet | 0 | 0 | 0 | 0 | 0 | 10 |
| obicetrapib | 0.3 | 0.25 | 0.37 | 0 | 0 | 1 |
| ctx310_angptl3_editing | 0.53 | 0.4 | 0 | 0 | 0 | 3 |
| lepodisiran | 0 | 0 | 0.94 | 0 | 0 | 2 |
| zerlasiran | 0 | 0 | 0.8 | 0 | 0 | 2 |
| plozasiran | 0 | 0.08 | 0 | 0 | 0 | 2 |
| olezarsen | 0 | 0.08 | 0 | 0 | 0 | 3 |
| tirzepatide | 0 | 0 | 0 | 0.3 | 2.1 | 3 |
| retatrutide | 0 | 0 | 0 | 0.35 | 14 | 3 |
| dfv890_nlrp3 | 0 | 0 | 0 | 0.4 | 0 | 1 |
| canakinumab_tet2_guided | 0 | 0 | 0 | 0.37 | 0 | 6 |
| yap101_sav1_gene_therapy | 0 | 0 | 0 | 0 | 0 | 9 |
| biovat_hf_patch | 0 | 0 | 0 | 0 | 0 | 10 |
| cd34_cell_therapy | 0 | 0 | 0 | 0 | 0 | 8 |
| coronary_sinus_reducer | 0 | 0 | 0 | 0 | 0 | 7 |
| pcsk9_vaccine | 0.2 | 0 | 0 | 0 | 0 | 2 |
| car_treg_oxldl | 0 | 0 | 0 | 0.3 | 0 | 10 |
| ab1002_gene_therapy | 0 | 0 | 0 | 0 | 0 | 8 |

## A.4 Full agent register — evidence and approval state

| therapy_id | evidence_tier | regulatory_status | onset_days | plateau_days | interval_days |
|---|---|---|---|---|---|
| statin_high_intensity | META_ANALYSIS | approved | 7 | 28 | 1 |
| statin_moderate_intensity | META_ANALYSIS | approved | 7 | 28 | 1 |
| ezetimibe | RANDOMIZED_OUTCOMES | approved | 7 | 14 | 1 |
| bempedoic_acid | RANDOMIZED_OUTCOMES | approved | 7 | 28 | 1 |
| evolocumab | RANDOMIZED_OUTCOMES | approved | 3 | 14 | 14 |
| alirocumab | RANDOMIZED_OUTCOMES | approved | 3 | 14 | 14 |
| inclisiran | RANDOMIZED_SURROGATE | approved | 14 | 60 | 182 |
| enlicitide_oral | RANDOMIZED_SURROGATE | approved | 7 | 28 | 1 |
| azd0780_oral | RANDOMIZED_SURROGATE | phase_3 | 7 | 28 | 1 |
| verve102_base_editing | RANDOMIZED_SURROGATE | phase_1 | 14 | 28 | 0 |
| pelacarsen | RANDOMIZED_TOPLINE | phase_3 | 14 | 60 | 28 |
| olpasiran | RANDOMIZED_SURROGATE | phase_3 | 14 | 90 | 84 |
| muvalaplin_oral | RANDOMIZED_SURROGATE | phase_2 | 7 | 28 | 1 |
| colchicine | RANDOMIZED_OUTCOMES | approved | 7 | 30 | 1 |
| canakinumab | RANDOMIZED_OUTCOMES | not_approved_for_ascvd | 7 | 30 | 90 |
| ziltivekimab | RANDOMIZED_OUTCOMES | investigational | 7 | 30 | 30 |
| methotrexate_low_dose | RANDOMIZED_OUTCOMES | not_indicated | 7 | 30 | 7 |
| nlrp3_inhibitor | MECHANISTIC | investigational | 7 | 30 | 1 |
| semaglutide_2_4mg | RANDOMIZED_OUTCOMES | approved | 14 | 168 | 7 |
| empagliflozin | RANDOMIZED_OUTCOMES | approved | 14 | 90 | 1 |
| sglt2_inhibitor_hf | RANDOMIZED_OUTCOMES | approved | 14 | 90 | 1 |
| bp_intensification | META_ANALYSIS | approved | 7 | 42 | 1 |
| aspirin_secondary_prevention | META_ANALYSIS | approved | 1 | 7 | 1 |
| icosapent_ethyl | RANDOMIZED_OUTCOMES | approved | 14 | 90 | 1 |
| smoking_cessation | OBSERVATIONAL | recommended | 30 | 365 | 1 |
| exercise_cardiac_rehab | META_ANALYSIS | recommended | 30 | 180 | 1 |
| diet_ldl_lowering | RANDOMIZED_SURROGATE | recommended | 14 | 90 | 1 |
| pci_stable_cad | RANDOMIZED_OUTCOMES | approved | 1 | 30 | 0 |
| cabg_ischemic_cm | RANDOMIZED_OUTCOMES | approved | 7 | 90 | 0 |
| ipsc_cm_sheet | CONDITIONAL_APPROVAL | conditional_japan | 30 | 180 | 0 |
| obicetrapib | RANDOMIZED_SURROGATE | phase_3_outcomes_pending | 7 | 28 | 1 |
| ctx310_angptl3_editing | RANDOMIZED_SURROGATE | phase_1b | 14 | 60 | 0 |
| lepodisiran | RANDOMIZED_SURROGATE | phase_3 | 14 | 60 | 182 |
| zerlasiran | RANDOMIZED_SURROGATE | phase_2 | 14 | 60 | 182 |
| plozasiran | RANDOMIZED_SURROGATE | approved_sHTG | 14 | 90 | 91 |
| olezarsen | RANDOMIZED_SURROGATE | filed_sHTG | 14 | 60 | 30 |
| tirzepatide | RANDOMIZED_OUTCOMES | approved | 7 | 180 | 7 |
| retatrutide | RANDOMIZED_SURROGATE | phase_3 | 7 | 240 | 7 |
| dfv890_nlrp3 | RANDOMIZED_SURROGATE | phase_2a | 7 | 28 | 1 |
| canakinumab_tet2_guided | OBSERVATIONAL | post_hoc_hypothesis | 14 | 90 | 90 |
| yap101_sav1_gene_therapy | MECHANISTIC | phase_1 | 30 | 180 | 0 |
| biovat_hf_patch | RANDOMIZED_SURROGATE | phase_2 | 30 | 90 | 0 |
| cd34_cell_therapy | RANDOMIZED_SURROGATE | phase_3 | 30 | 180 | 0 |
| coronary_sinus_reducer | RANDOMIZED_SURROGATE | approved_limited | 30 | 90 | 0 |
| pcsk9_vaccine | MECHANISTIC | preclinical | 60 | 180 | 90 |
| car_treg_oxldl | MECHANISTIC | preclinical | 30 | 180 | 0 |
| ab1002_gene_therapy | RANDOMIZED_SURROGATE | phase_2 | 30 | 180 | 0 |

### Agents that deliberately earn zero credit in specific populations

These are not bugs. They are the credit rules refusing to transport evidence where
it does not belong:

- **Statins and PCSK9 antibodies in `ischemic_cm_heart_failure`** — lipid lowering
  has no demonstrated event benefit in established severe heart failure.
- **Evolocumab in `post_mi_pci`** — no anchor in that population.
- **All Lp(a)-lowering agents** — Lp(a)HORIZON missed its primary endpoint, so the
  class earns no event credit despite large biomarker effects.
- **Colchicine after modern PCI** — CLEAR SYNERGY was neutral (HR 0.99).
- **Methotrexate** — CIRT was neutral (HR 1.01); it does not engage IL-6.
- **AB-1002 / GenePHIT** — tested only in *non-ischaemic* cardiomyopathy, so it
  earns nothing in IHD. A dedicated exclusion-control population proves this.
- **All imaging anchors** (ASTEROID, GLAGOV, FLOWPROMOTE) — plaque imaging never
  enters event credit, only the burden-reversal criterion.

---

# PART B — THE TRIAL EVIDENCE BASE (all 48 anchors)

**Status: EMPIRICAL (published) → EVIDENCE-GRADED (as used).**


## B.1 Anchors with a hazard ratio and confidence interval

These carry randomised or meta-analytic event evidence. A **neutral** trial is one
whose published interval spans 1.0; the model samples such trials around the null
rather than around their noisy point estimate, so uncertainty is retained without
inventing a direction.

| anchor_id | therapy_id | population_id | hr | ci_low | ci_high | followup_years | n_participants | neutral |
|---|---|---|---|---|---|---|---|---|
| CTT-STATIN | statin_high_intensity | any_ascvd | 0.78 | 0.76 | 0.8 | 5 | 170000 |  |
| BPLTTC-SBP | bp_intensification | any_ascvd | 0.9 | 0.88 | 0.92 | 4.2 | 344716 |  |
| IMPROVE-IT | ezetimibe | secondary_post_mi | 0.936 | 0.897 | 0.978 | 6 | 18144 |  |
| CLEAR-OUTCOMES | bempedoic_acid | statin_intolerant_ascvd | 0.87 | 0.79 | 0.96 | 3.4 | 13970 |  |
| FOURIER | evolocumab | stable_ccs | 0.85 | 0.78 | 0.93 | 2.2 | 27564 |  |
| ODYSSEY-OUTCOMES | alirocumab | secondary_post_mi | 0.85 | 0.78 | 0.93 | 2.8 | 18924 |  |
| LODOCO2 | colchicine | stable_ccs | 0.69 | 0.57 | 0.83 | 2.4 | 5522 |  |
| COLCOT | colchicine | secondary_post_mi | 0.77 | 0.61 | 0.96 | 1.9 | 4745 |  |
| CLEAR-SYNERGY | colchicine | post_mi_pci | 0.99 | 0.85 | 1.16 | 3 | 7062 | NEUTRAL |
| CANTOS-150 | canakinumab | post_mi_hscrp_high | 0.85 | 0.74 | 0.98 | 3.7 | 10061 |  |
| ZEUS | ziltivekimab | ascvd_ckd_inflammation | 0.99 | 0.88 | 1.11 | 3 | 6385 | NEUTRAL |
| CIRT | methotrexate_low_dose | stable_ccs | 1.01 | 0.82 | 1.25 | 2.3 | 4786 | NEUTRAL |
| SELECT | semaglutide_2_4mg | ascvd_obesity_no_diabetes | 0.8 | 0.72 | 0.9 | 3.3 | 17604 |  |
| EMPA-REG | empagliflozin | t2d_ascvd | 0.86 | 0.74 | 0.99 | 3.1 | 7020 |  |
| DAPA-HF | sglt2_inhibitor_hf | ischemic_cm_heart_failure | 0.74 | 0.65 | 0.85 | 1.5 | 4744 |  |
| ATT-ASPIRIN | aspirin_secondary_prevention | any_ascvd | 0.81 | 0.77 | 0.86 | 2 | 135000 |  |
| REDUCE-IT | icosapent_ethyl | ascvd_high_triglycerides | 0.75 | 0.68 | 0.83 | 4.9 | 8179 |  |
| QUIT-POSTMI | smoking_cessation | secondary_post_mi | 0.64 | 0.58 | 0.71 | 5 | 12603 |  |
| EXERCISE-REHAB | exercise_cardiac_rehab | secondary_post_mi | 0.74 | 0.64 | 0.86 | 1 | 14486 |  |
| ISCHEMIA | pci_stable_cad | stable_ccs | 0.93 | 0.83 | 1.04 | 3.2 | 5179 | NEUTRAL |
| STICH-CABG | cabg_ischemic_cm | ischemic_cm_heart_failure | 0.84 | 0.73 | 0.97 | 9.8 | 1212 |  |
| LPA-HORIZON | pelacarsen | ascvd_high_lpa | 1 | 1 | 1 | 4 | 8323 | NEUTRAL |
| BROADWAY-EVENTS | obicetrapib | any_ascvd | 0.79 | 0.54 | 1.15 | 1 | 2530 | NEUTRAL |
| SURPASS-CVOT | tirzepatide | t2d_ascvd | 0.92 | 0.83 | 1.01 | 4 | 13299 | NEUTRAL |
| CANTOS-TET2 | canakinumab_tet2_guided | ascvd_tet2_chip | 0.38 | 0.2 | 0.72 | 3.7 | 103 |  |
| MR-LIFELONG-LDL | statin_high_intensity | lifetime_exposure | 0.46 | 0.39 | 0.54 | 50 | 200000 |  |

**26 anchors carry an effect estimate.** Total participants represented: **1,062,161**.

## B.2 Anchors without an event estimate (surrogate, safety or pending)

These support biomarker or mechanistic credit only. They can never earn event
credit, which is why so many agents stall at the same criterion.

| anchor_id | therapy_id | population_id | endpoint | n_participants | evidence_tier |
|---|---|---|---|---|---|
| ORION-11 | inclisiran | stable_ccs | LDL-C percent change at day 510 | 1617 | RANDOMIZED_SURROGATE |
| CORALREEF-LIPIDS | enlicitide_oral | stable_ccs | LDL-C percent change at week 24 | 2909 | RANDOMIZED_SURROGATE |
| PURSUIT | azd0780_oral | stable_ccs | LDL-C percent change at week 12 | 428 | RANDOMIZED_SURROGATE |
| HEART-2 | verve102_base_editing | hefh_primary | Time-averaged LDL-C percent change | 35 | RANDOMIZED_SURROGATE |
| RIHEART-PIVOTAL | ipsc_cm_sheet | ischemic_cm_heart_failure | LVEF and exercise capacity in 8 patients | 8 | CONDITIONAL_APPROVAL |
| OCEAN-A | olpasiran | ascvd_high_lpa | Cardiovascular outcomes pending | 7297 | RANDOMIZED_SURROGATE |
| MUV-Phase2 | muvalaplin_oral | ascvd_high_lpa | Lp(a) percent change | 233 | RANDOMIZED_SURROGATE |
| PREVAIL | obicetrapib | any_ascvd | 3-point MACE, primary completion Nov 2026 | 9541 | RANDOMIZED_SURROGATE |
| CTX310-P1A | ctx310_angptl3_editing | hefh_primary | LDL-C and triglyceride change at 1 year | 15 | RANDOMIZED_SURROGATE |
| ALPACA | lepodisiran | ascvd_high_lpa | Time-averaged Lp(a) change | 320 | RANDOMIZED_SURROGATE |
| ACCLAIM-LPA | lepodisiran | ascvd_high_lpa | MACE, estimated completion 2029 | 17300 | RANDOMIZED_SURROGATE |
| ZERLASIRAN-P2 | zerlasiran | ascvd_high_lpa | Time-averaged Lp(a) change over 36 weeks | 367 | RANDOMIZED_SURROGATE |
| SHASTA-3-4 | plozasiran | severe_hypertriglyceridaemia | Triglycerides and acute pancreatitis | 757 | RANDOMIZED_SURROGATE |
| CORE-TIMI-72 | olezarsen | severe_hypertriglyceridaemia | Triglycerides and acute pancreatitis | 1063 | RANDOMIZED_SURROGATE |
| BIOVAT-HF | biovat_hf_patch | advanced_ischemic_hf | Wall thickness, LVEF and KCCQ at 3 months | 16 | RANDOMIZED_SURROGATE |
| SALVADOR-HF | yap101_sav1_gene_therapy | advanced_ischemic_hf | Safety, phase 1 ongoing | 12 | MECHANISTIC |
| GENEPHIT | ab1002_gene_therapy | non_ischemic_cardiomyopathy | LVEF, NYHA and 6-minute walk | 173 | RANDOMIZED_SURROGATE |
| CD34-ANGINA | cd34_cell_therapy | refractory_angina_no_option | Angina frequency and exercise time | 167 | RANDOMIZED_SURROGATE |
| REMEDY-PILOT | coronary_sinus_reducer | refractory_angina_no_option | Myocardial perfusion reserve versus sham | 54 | RANDOMIZED_SURROGATE |
| ASTEROID | statin_high_intensity | stable_ccs | Percent atheroma volume by IVUS | 507 | RANDOMIZED_SURROGATE |
| GLAGOV | evolocumab | stable_ccs | Percent atheroma volume by IVUS | 968 | RANDOMIZED_SURROGATE |
| FLOWPROMOTE | ezetimibe | stable_ccs | Total and noncalcified plaque volume by CT | 104 | RANDOMIZED_SURROGATE |

## B.3 The neutral and failed trials, explicitly

| Trial | Agent | Result | Consequence in this study |
|---|---|---|---|
| **Lp(a)HORIZON** | pelacarsen | HR 1.00, primary endpoint **not met** (n=8,323) | Entire Lp(a)-lowering class earns no event credit |
| **CLEAR SYNERGY** | colchicine | HR 0.99 (0.85–1.16), n=7,062 | No colchicine credit after modern PCI |
| **CIRT** | methotrexate | HR 1.01 (0.82–1.25), n=4,786 | No credit; no IL-6 engagement |
| **ZEUS** | ziltivekimab | HR 0.99 (0.88–1.11), n=6,385 | No credit despite 60% hsCRP reduction |
| **SURPASS-CVOT** | tirzepatide | HR 0.92 (0.83–1.01) vs dulaglutide | Active-comparator; interval spans 1.0 |
| **BROADWAY** | obicetrapib | HR 0.79 (0.54–1.15), not powered | Surrogate credit only; PREVAIL pending |
| **ISCHEMIA** | PCI | HR 0.93 (0.83–1.04), n=5,179 | No event credit for PCI in stable disease |

This table is the single most important honesty feature of the evidence layer: five
of the seventeen outcomes-tier anchors are **neutral**, and the model treats them
as neutral rather than quietly using their point estimates.

---

# PART C — THE POPULATIONS (all 19)

**Status: EVIDENCE-GRADED SCENARIO.** Baseline rates are scenario inputs with
intervals, not patient estimates.

| population_id | label | baseline_annual_mace_rate | baseline_rate_ci_low | baseline_rate_ci_high |
|---|---|---|---|---|
| secondary_post_mi | Post-myocardial-infarction secondary prevention | 0.045 | 0.038 | 0.053 |
| post_mi_pci | Early post-MI treated with PCI | 0.05 | 0.042 | 0.059 |
| post_mi_hscrp_high | Post-MI with hsCRP at least 2 mg/L | 0.05 | 0.042 | 0.059 |
| stable_ccs | Stable chronic coronary syndrome | 0.03 | 0.025 | 0.036 |
| ascvd_ckd_inflammation | ASCVD with chronic kidney disease and inflammation | 0.06 | 0.051 | 0.07 |
| statin_intolerant_ascvd | ASCVD with documented statin intolerance | 0.04 | 0.034 | 0.047 |
| t2d_ascvd | Type 2 diabetes with established ASCVD | 0.04 | 0.034 | 0.047 |
| ascvd_obesity_no_diabetes | ASCVD with obesity and without diabetes | 0.03 | 0.025 | 0.036 |
| ascvd_high_triglycerides | ASCVD with elevated triglycerides on statin | 0.04 | 0.034 | 0.047 |
| ascvd_high_lpa | ASCVD with elevated lipoprotein(a) | 0.04 | 0.034 | 0.047 |
| hefh_primary | Heterozygous familial hypercholesterolaemia | 0.02 | 0.016 | 0.025 |
| ischemic_cm_heart_failure | Ischaemic cardiomyopathy with severe heart failure | 0.08 | 0.068 | 0.094 |
| any_ascvd | Generic atherosclerotic disease reference | 0.035 | 0.029 | 0.042 |
| ascvd_tet2_chip | ASCVD with TET2 clonal haematopoiesis | 0.055 | 0.035 | 0.08 |
| advanced_ischemic_hf | Advanced ischaemic HF, LVEF <=35% | 0.07 | 0.05 | 0.095 |
| refractory_angina_no_option | Refractory angina, no revascularisation option | 0.045 | 0.03 | 0.065 |
| severe_hypertriglyceridaemia | Severe hypertriglyceridaemia (TG >=500 mg/dL) | 0.03 | 0.018 | 0.048 |
| lifetime_exposure | Primordial prevention: lifelong LDL-C exposure scenario | 0.004 | 0.002 | 0.008 |
| non_ischemic_cardiomyopathy | Non-ischaemic cardiomyopathy (exclusion control) | 0.03 | 0.02 | 0.045 |

## C.2 Baseline risk-factor profile

| population_id | ldl_c_mg_dl | apob_mg_dl | hscrp_mg_l | systolic_bp_mmhg | lp_a_nmol_l |
|---|---|---|---|---|---|
| secondary_post_mi | 100 | 90 | 2.5 | 135 | 50 |
| post_mi_pci | 95 | 88 | 3 | 132 | 50 |
| post_mi_hscrp_high | 100 | 90 | 4.2 | 134 | 55 |
| stable_ccs | 95 | 85 | 2 | 132 | 45 |
| ascvd_ckd_inflammation | 90 | 85 | 4.5 | 138 | 60 |
| statin_intolerant_ascvd | 150 | 125 | 2.4 | 134 | 50 |
| t2d_ascvd | 95 | 90 | 2.6 | 136 | 45 |
| ascvd_obesity_no_diabetes | 100 | 92 | 3 | 133 | 45 |
| ascvd_high_triglycerides | 85 | 88 | 2.4 | 132 | 50 |
| ascvd_high_lpa | 90 | 85 | 2.2 | 130 | 180 |
| hefh_primary | 190 | 150 | 1.8 | 128 | 70 |
| ischemic_cm_heart_failure | 90 | 85 | 3.2 | 120 | 50 |
| any_ascvd | 110 | 95 | 2.4 | 134 | 50 |
| ascvd_tet2_chip | 95 | 85 | 4.5 | 138 | 60 |
| advanced_ischemic_hf | 85 | 75 | 3.5 | 120 | 50 |
| refractory_angina_no_option | 95 | 85 | 2.8 | 135 | 55 |
| severe_hypertriglyceridaemia | 110 | 110 | 3 | 132 | 50 |
| lifetime_exposure | 160 | 115 | 1.2 | 125 | 50 |
| non_ischemic_cardiomyopathy | 95 | 80 | 2.5 | 118 | 40 |

Two populations exist purely as controls:

- **`non_ischemic_cardiomyopathy`** — deliberately non-ischaemic. It holds the
  GenePHIT/AB-1002 anchor and proves that a heart-failure gene therapy tested only
  outside IHD earns no credit inside IHD.
- **`lifetime_exposure`** — a Mendelian-randomisation scenario, not a trial
  population. No treatment anchor is transportable to it.

---

# PART D — MOLECULAR AND STRUCTURAL EXPERIMENTS

## D.1 Crystallographic interface analysis — **EMPIRICAL / deterministic**

Six structures were used, four of them PCSK9-related:

| PDB | Resolution | Content |
|---|---|---|
| `3GCX` | 1.98 Å | PCSK9 + LDLR EGF-A (Kwon 2008) |
| `3BPS` | 2.30 Å | PCSK9 catalytic domain + LDLR EGF-A |
| `2P4E` | 1.97 Å | PCSK9 apo / domain map |
| `3P5C` | 3.30 Å | PCSK9 + full LDLR ectodomain (Lo Surdo 2011) |
| `4NMX` | — | PCSK9 comparison structure |
| `4HHB` | — | Haemoglobin, method control |

**Method.** Shrake–Rupley solvent-accessible surface area, probe radius 1.4 Å,
960 test points per atom. Burial = SASA(isolated) − SASA(bound).

**Result (3BPS):** total interface buried surface **591.48 Å²**; PCSK9 side
2676.02 Å², LDLR side 639.03 Å²; **317 atom contacts** below 5.0 Å across 15 residues.

### Contact residues ranked by number of contacts

| residue | contacts_<5A | min_dist_A | buried_A2 | pct_buried |
|---|---|---|---|---|
| PHE379 | 70 | 2.75 | 77.4 | 99.4 |
| THR377 | 46 | 2.54 | 43 | 60.5 |
| CYS378 | 31 | 3.27 | 37 | 100 |
| SER153 | 28 | 3.11 | 66.2 | 39.1 |
| ILE369 | 24 | 3.59 | 64.9 | 94.9 |
| ASP238 | 23 | 2.88 | 43.2 | 56.9 |
| VAL380 | 20 | 3.44 | 37.5 | 79.9 |
| ARG194 | 19 | 2.95 | 20.6 | 45 |
| PRO155 | 18 | 3.65 | 22.9 | 92.4 |
| ASP374 | 9 | 3.44 | 38.8 | 35.9 |
| ILE154 | 8 | 3.89 |  |  |
| SER381 | 6 | 4.58 | 13.5 | 31.7 |
| SER372 | 5 | 3.77 |  |  |
| CYS375 | 5 | 3.47 | 12.5 | 26.2 |
| ALA239 | 3 | 4.05 | 5.16 | 95.9 |

### The "flat PPI" conclusion

The interface is **591 Å² spread over 15+ residues on discontinuous loops**. A
Lipinski-compliant small molecule (<500 Da) covers roughly 250–350 Å² — under 40% of
the contact zone. This quantitatively explains why 15 years of classical
small-molecule medicinal chemistry produced no approved drug against this
interface, and why the modalities that *did* work are either macrocyclic peptides
(~1.5 kDa, e.g. enlicitide/MK-0616) or upstream genetic silencing (inclisiran siRNA,
VERVE-102 base editing).

**This is the strongest validated finding in the molecular work.**

## D.2 Cross-structure variability — **EMPIRICAL**

The same residue measured in three independent crystal forms. Large standard
deviations mean the residue's burial is context-dependent, not a fixed property.

| ResID | ResName | Burial_3GCX_A2 | Burial_3BPS_A2 | Burial_3P5C_A2 | Mean_Burial_A2 | Std_Burial_A2 |
|---|---|---|---|---|---|---|
| 153 | SER | 57.89 | 66.24 | 26.5 | 50.21 | 17.11 |
| 154 | ILE | 2.92 | 4.553 | 37.88 | 15.12 | 16.11 |
| 155 | PRO | 21.28 | 22.95 | 8.908 | 17.71 | 6.26 |
| 194 | ARG | 33.07 | 20.64 | 26.05 | 26.59 | 5.09 |
| 238 | ASP | 36.74 | 43.17 | 31.87 | 37.26 | 4.63 |
| 367 | ASP | 0 | 13.5 | 14.73 | 9.41 | 6.67 |
| 368 | ILE | 0 | 0 | 0 | 0 | 0 |
| 369 | ILE | 65.16 | 64.91 | 50.7 | 60.26 | 6.76 |
| 374 | ASP | 36.31 | 38.81 | 26.9 | 34.01 | 5.13 |
| 375 | CYS | 12.4 | 12.54 | 6.032 | 10.32 | 3.04 |
| 376 | SER | 0 | 0 | 0 | 0 | 0 |
| 377 | THR | 41.08 | 42.99 | 54.93 | 46.33 | 6.13 |
| 378 | CYS | 36.18 | 36.97 | 29.22 | 34.12 | 3.48 |
| 379 | PHE | 81.92 | 77.43 | 92.34 | 83.89 | 6.24 |
| 380 | VAL | 28.29 | 37.46 | 27.05 | 30.93 | 4.64 |
| 381 | SER | 14.44 | 13.5 | 23.11 | 17.02 | 4.33 |

`SER153` is the cautionary case: 57.9 / 66.2 / 26.5 Å² across the three structures
(SD 17.1 Å²). Any single-structure claim about it would be unreliable.

## D.3 Experimental benchmark: SPR alanine scanning — **EMPIRICAL (wet-lab)**

Source: Kwon et al., *PNAS* 2008, 105(6):1820–1825. Surface plasmon resonance
measurement of PCSK9 binding to LDLR EGF-A. **This is real measured data, not a
model output.**

| Mutation | ResID | KD_uM | ddG_exp_kcal_mol | Assay_Method | Notes |
|---|---|---|---|---|---|
| WT | - | 0.63 | 0 | SPR | Baseline wild-type PCSK9 binding to LDLR EGF-A |
| D374Y | 374 | 0.025 | -1.91 | SPR | Severe gain-of-function mutation causing Familial Hypercholesterolemia (~25-fold higher affinity) |
| D374A | 374 | 2.5 | 0.82 | SPR | Alanine substitution removes negative charge / H-bond network (4-fold reduction) |
| F379A | 379 | 100 | 3 | SPR | Essential hydrophobic anchor; binding virtually abolished (>150-fold loss) |
| I369A | 369 | 14.8 | 1.87 | SPR | Hydrophobic core packing residue (23-fold loss) |
| R194A | 194 | 12 | 1.75 | SPR | Key basic residue interacting with acidic EGF-A loop (19-fold loss) |
| T377A | 377 | 5.7 | 1.31 | SPR | Interface polar contact (9-fold loss) |
| D238A | 238 | 4.2 | 1.13 | SPR | Salt-bridge partner near interface (6.7-fold loss) |
| V380A | 380 | 2.1 | 0.71 | SPR | Adjacent hydrophobic residue (3.3-fold loss) |
| S153A | 153 | 0.75 | 0.1 | SPR | N-terminal mature cleavage site; near-neutral effect (1.2-fold loss) |

**Key readings.** Wild-type Kd = 0.63 µM. The gain-of-function **D374Y** mutation
binds 25× tighter (Kd 0.025 µM, ΔΔG −1.91 kcal/mol) and causes familial
hypercholesterolaemia. **F379A** abolishes binding (Kd 100 µM, ΔΔG +3.00 kcal/mol,
>150-fold loss) — consistent with PHE379 being the most-buried, most-contacted
residue in the geometric analysis (77.4 Å², 99.4% buried, 70 contacts).

**The computational geometry and the wet-lab energetics agree on which residue
matters most.** That agreement is the validation.

## D.4 Computed burial vs measured ΔΔG — the correlation table

| Mutation | KD_uM | ddG_exp_kcal_mol | Mean_Burial_A2 | Std_Burial_A2 |
|---|---|---|---|---|
| D374A | 2.5 | 0.82 | 34.01 | 5.13 |
| F379A | 100 | 3 | 83.89 | 6.24 |
| I369A | 14.8 | 1.87 | 60.26 | 6.76 |
| R194A | 12 | 1.75 | 26.59 | 5.09 |
| T377A | 5.7 | 1.31 | 46.33 | 6.13 |
| D238A | 4.2 | 1.13 | 37.26 | 4.63 |
| V380A | 2.1 | 0.71 | 30.93 | 4.64 |
| S153A | 0.75 | 0.1 | 50.21 | 17.11 |

The ordering is broadly consistent (F379A: largest ΔΔG and largest burial; S153A:
near-neutral ΔΔG), but this is an **8-point comparison**, not a fitted or
cross-validated model. It supports a qualitative claim, not a predictive one.

## D.5 CRISPR base-editing specificity screen — **COMPUTATIONAL, first-pass only**

Target: human *PCSK9* intron-1 splice donor, `GGGCAGCGCCACCCACAGGA`, PAM `TGG`
(the VERVE-102 locus). Method: Cutting Frequency Determination mismatch scoring
against homologous human proprotein convertases plus intergenic matches.

| locus | chromosome | mismatches | CFD score | risk | consequence |
|---|---|---|---|---|---|
| PCSK9 (On-Target) | Chr1p32.3 | 0 | 1 | ON-TARGET THERAPEUTIC | Permanent loss-of-function; therapeutic LDL drop |
| PCSK1 (Convertase 1) | Chr5q15 | 2 | 0.08612 | VERY LOW / BORDERLINE | Endocrine hormone processing |
| PCSK2 (Convertase 2) | Chr20p12 | 4 | 0.000646 | NEGLIGIBLE (<0.1%) | Glucagon/insulin processing |
| FURIN (PCSK3) | Chr15q26 | 4 | 0.001117 | NEGLIGIBLE (<0.1%) | Essential general convertase |
| MBTPS1 (SKI-1 / S1P) | Chr16q24 | 3 | 0.01941 | VERY LOW / BORDERLINE | Sterol regulation |
| Intergenic Off-Target 1 | Chr8q21.1 | 1 | 0.265 | POTENTIAL RISK (Deep-Seq Needed) | No known transcript |
| Intergenic Off-Target 2 | Chr12p11.2 | 1 | 0.4 | POTENTIAL RISK (Deep-Seq Needed) | Intron of uncharacterized lncRNA |

**Two intergenic loci score 0.265 and 0.400 with a single mismatch each — flagged
POTENTIAL RISK.** The on-target CFD is 1.0 by definition.

**Limitation, stated plainly:** in-silico string matching and CFD scoring are
filtering tools. A regulatory dossier requires empirical genome-wide assays
(GUIDE-seq, CIRCLE-seq, DISCOVER-seq) in primary human hepatocytes, targeted deep
sequencing at >100,000× coverage, and long-term registry follow-up. **Nothing here
clears VERVE-102 for safety.** The guide sequence itself is unpublished; the locus
is analysed from the public splice-donor context.

## D.6 Monte Carlo conformational search — **EXPLORATORY**

10,001 steps recorded in `mc_energy_trace.csv`, final energy −9.47 kcal/mol. This is
a sampling trace on a model potential, with no experimental observable attached.
It demonstrates the search ran; it does not establish a binding affinity.

## D.7 Docking and de novo design — **NON-USABLE, FLAGGED**

`results/tables/docking_and_denovo_results.csv` contains 8 compounds: 3 de novo
designs, 3 approved drugs, 2 preclinical leads. **This table must not be used to
rank affinity.** Three independent defects:

| Rank | Compound | Type | dG_bind_kcal_mol | Est_Kd_uM | MW | LogP | TPSA |
|---|---|---|---|---|---|---|---|
| 1 | PCSK9-Inh-D3 | De Novo Macrocyclic Inhibitor | 5.28 | 1e+06 | 524.6 | 0.41 | 173.9 |
| 2 | PCSK9-Inh-D1 | De Novo Designed Small Molecule | 17.47 | 1e+06 | 403.3 | 5.46 | 66.4 |
| 3 | Valsartan | Approved Drug (ARB) | 53.8 | 1e+06 | 370.5 | 5.4 | 55.1 |
| 4 | PF-06446846 | Preclinical Lead | 65.82 | 1e+06 | 315.8 | 4 | 55.4 |
| 5 | Rosuvastatin | Approved Drug (Statin) | 68.7 | 1e+06 | 481.5 | 2.4 | 140.9 |
| 6 | Ezetimibe | Approved Drug (Cholesterol Absorption) | 75.05 | 1e+06 | 365.4 | 4.41 | 49.3 |
| 7 | SBC-115076 | Preclinical Lead | 81.31 | 1e+06 | 246.7 | 3.3 | 42 |
| 8 | PCSK9-Inh-D2 | De Novo Designed Small Molecule | 90.81 | 1e+06 | 453.5 | 3.15 | 116.1 |

1. **Every ΔG is positive** (+5.28 to +90.81 kcal/mol). A positive binding free
   energy means *no favourable binding predicted*. The "Rank" column therefore
   orders compounds from least-unfavourable to most-unfavourable; it is not an
   affinity ranking, and rank 1 is not a hit.
2. **Every `Est_Kd_uM` is identical: 1,000,000 µM (1 M).** The conversion saturated
   at its clamp for all eight compounds, so the column carries zero discriminating
   information.
3. **Molecular weights disagree with the screening table** for four shared
   compounds (see D.8), so at least one of the two tables has incorrect structures.

**Conclusion: the docking exercise is reported for completeness and is not
interpreted.** The de novo designs (PCSK9-Inh-D1/D2/D3) are unvalidated sketches.
The sound conclusion about this target comes from the interface geometry (D.1),
not from these scores.

## D.8 Virtual screening library — **NON-USABLE for biologics, heuristic overall**

28 compounds were scored on a combined drug-likeness + PCSK9-affinity heuristic.

| Rank | Name | MW | logP | HBD | HBA | DrugLikeness | PCSK9Affinity | CombinedScore |
|---|---|---|---|---|---|---|---|---|
| 1 | Anakinra_mimic | 428.5 | 1.2 | 5 | 8 | 100 | 100 | 100 |
| 2 | Empagliflozin | 450.9 | 1.6 | 3 | 8 | 100 | 100 | 100 |
| 3 | Ezetimibe | 409.4 | 3.9 | 2 | 4 | 100 | 100 | 100 |
| 4 | Alirocumab_frag | 402.8 | 4.5 | 1 | 4 | 100 | 100 | 100 |
| 5 | Inclisiran | 337.3 | 1.5 | 2 | 6 | 100 | 100 | 100 |
| 6 | Valsartan | 435.5 | 3.6 | 2 | 5 | 100 | 98 | 98.9 |
| 7 | Rosuvastatin | 481.5 | 0.3 | 4 | 10 | 100 | 96 | 97.8 |
| 8 | Ramipril | 416.5 | 1.2 | 2 | 7 | 100 | 96 | 97.8 |
| 9 | PF-06446846 | 357.8 | 3.5 | 2 | 5 | 100 | 95 | 97.2 |
| 10 | Berberine | 336.4 | 1.8 | 0 | 7 | 100 | 95 | 97.2 |
| 11 | Quercetin | 302.2 | 1.5 | 5 | 7 | 100 | 95 | 97.2 |
| 12 | Clopidogrel | 321.8 | 2.5 | 1 | 4 | 100 | 95 | 97.2 |
| 13 | Ziltivekimab | 434.5 | 4.1 | 2 | 6 | 100 | 93 | 96.2 |
| 14 | Atorvastatin | 558.6 | 4.1 | 3 | 7 | 95 | 96 | 95.6 |
| 15 | Ticagrelor | 522.6 | 2.5 | 3 | 7 | 95 | 93 | 93.9 |
| 16 | SBC-115076 | 261.7 | 2.9 | 2 | 3 | 100 | 86 | 92.3 |
| 17 | Colchicine | 399.4 | 1.3 | 1 | 7 | 100 | 85 | 91.8 |
| 18 | Resveratrol | 228.2 | 2.8 | 3 | 3 | 100 | 85 | 91.8 |
| 19 | Simvastatin | 418.6 | 4.7 | 1 | 6 | 100 | 80 | 89 |
| 20 | Lomitapide | 363.3 | 5.2 | 1 | 3 | 100 | 80 | 89 |
| 21 | Evolocumab_frag | 296.4 | 1 | 3 | 5 | 100 | 76 | 86.8 |
| 22 | Aspirin | 180.2 | 1.2 | 1 | 4 | 100 | 73 | 85.2 |
| 23 | Nicotinamide | 122.1 | -0.4 | 2 | 2 | 100 | 66 | 81.3 |
| 24 | Canakinumab | 518.6 | 1.8 | 7 | 9 | 85 | 78 | 81.2 |
| 25 | Doxycycline | 444.4 | -0.2 | 7 | 11 | 80 | 81 | 80.6 |
| 26 | Metformin | 129.2 | -1.4 | 3 | 2 | 100 | 58 | 76.9 |
| 27 | Semaglutide | 283.5 | 5.2 | 1 | 2 | 85 | 55 | 68.5 |
| 28 | Inclisiran_analog | 591.7 | 7.2 | 1 | 7 | 70 | 60 | 64.5 |

**Defect 1 — saturation.** 18 of 28 compounds score ≥90 and five score exactly
100.0. A heuristic that puts two-thirds of a library in the top decile is not
discriminating.

**Defect 2 — biologics given small-molecule molecular weights.** These rows are
physically wrong and their drug-likeness scores are meaningless:

| name | MW in table | actual modality | true MW |
|---|---|---|---|
| Canakinumab | 519 | monoclonal antibody | ~145,000 Da |
| Semaglutide | 284 | GLP-1 peptide | 4,113.6 Da |
| Inclisiran | 337 | double-stranded siRNA | ~17,300 Da |
| Ziltivekimab | 434 | monoclonal antibody | ~145,000 Da |

A Lipinski-style drug-likeness score cannot be applied to an antibody or an siRNA
at all. `Alirocumab_frag`, `Evolocumab_frag`, `Anakinra_mimic` and
`Inclisiran_analog` are explicitly fragments or mimics, which is defensible;
the four rows above are not.

**Defect 3 — MW inconsistency between the two compound tables:**

| compound | screening MW | docking MW | difference |
|---|---|---|---|
| Ezetimibe | 409 | 365 | 44 |
| Valsartan | 436 | 370 | 65 |
| PF-06446846 | 358 | 316 | 42 |
| SBC-115076 | 262 | 247 | 15 |

For reference, the accepted molecular weights of ezetimibe (409.4 Da) and valsartan
(435.5 Da) match the **screening** table, which suggests the docking table's
structures are the incorrect ones.

**Net position on the chemistry:** the crystallographic interface analysis and its
SPR benchmark are sound and are the basis for the modality conclusion. The docking
and screening tables are retained for auditability but are **excluded from every
conclusion in this study**. They are the weakest artefacts in the repository and
should be regenerated with a validated scoring function and correct structures
before being used for anything.

---

# PART E — DYNAMICAL SYSTEMS EXPERIMENTS

## E.1 Plaque biomechanics, 60-month progression — **EXPLORATORY**

Coupled ODEs for intimal lipid accumulation, necrotic core expansion, fibrous cap
remodelling (VSMC collagen synthesis vs MMP degradation) and wall stress via
Laplace's law for a thin-walled cylinder:

$$\sigma_{hoop} = \frac{P_{lumen} \cdot R_{lumen}}{T_{cap}}$$

Vulnerability thresholds used: cap thickness ≤65 µm and peak tensile stress
≥300 kPa. 44 timeline rows across regimens.

| Regimen | rows | final_cap_um | final_core_mm2 | max_stress_kPa |
|---|---|---|---|---|
| Dual Cure (VERVE-102 + Anti-IL-1β) | 11 | 210.4 | 0.26 | 155.5 |
| High-Intensity Statin | 11 | 56.37 | 0.847 | 426.2 |
| Untreated (High Risk) | 11 | 15 | 1.579 | 1511 |
| VERVE-102 (CRISPR Base Edit) | 11 | 156.3 | 0.26 | 171.7 |

**Limitations.** Axisymmetric cylinder, not a patient-specific 3D vessel; linear
elastic material, whereas real vascular tissue is nonlinear, anisotropic and
hyperelastic (Holzapfel–Gasser–Ogden); uniform pressure instead of fluid–structure
interaction with wall-shear-stress gradients. **The rupture months are illustrative
outputs of chosen parameters, not predictions for any patient.**

## E.2 Vascular cytokine network, 180 days — **EXPLORATORY**

Seven-variable ODE chain:

`oxLDL → NF-κB → NLRP3 inflammasome → Caspase-1 → IL-1β → IL-6 → hs-CRP → MMP-1/9`

65 rows across regimens (untreated, statin, canakinumab, colchicine, combination).

| Regimen | rows | final_IL1b | final_IL6 | final_hsCRP | final_MMP_pct |
|---|---|---|---|---|---|
| Canakinumab (Anti-IL-1β, CANTOS) | 13 | 10.87 | 10.76 | 13.29 | 154.6 |
| Colchicine (NLRP3 Inhibitor) | 13 | 16.31 | 16.14 | 19.94 | 231.9 |
| Curative Combo (PCSK9 silencer + NLRP3i) | 13 | 1.45 | 1.44 | 1.78 | 20.8 |
| High-Intensity Statin | 13 | 30.83 | 30.52 | 37.72 | 438.8 |
| Untreated High Risk | 13 | 72.48 | 71.72 | 88.61 | 1031 |

**What it illustrates.** Lipid lowering alone leaves residual inflammatory risk;
combining lipid suppression with inflammasome inhibition is synergistic in the
model. This is consistent with CANTOS and the colchicine trials, but it is a
7-variable caricature of a system with hundreds of species, spatial gradients,
regulatory T cells and negative feedback (IL-10, SOCS).

---

# PART F — THE EVIDENCE-GRADED SIMULATION PROGRAMME

Nine studies were run. Each writes to its own results tree, and invented work is
kept out of the evidence manifest.

## F.1 Study 1 — single-agent cure search

**Question:** can any single approved or investigational agent satisfy a
five-criterion definition of cure?

**Criteria** (threshold 0.80 each): durable driver control · atherosclerotic burden
reversal · myocardial restoration · risk normalisation · randomised proof.

**Result: all 47 agents satisfy exactly 1 of 5** (randomised proof, where they have
it). Highest readiness indices:

| therapy_id | strict_rrr_pct | durable_driver_control | burden_reversal | myocardial_restoration | risk_normalisation | randomised_proof | cure_readiness_index |
|---|---|---|---|---|---|---|---|
| evolocumab | 12.56 | 0.2 | 0.55 | 0 | 0.05023 | 1 | 0.36 |
| ezetimibe | 5.701 | 0.1 | 0.6 | 0 | 0.0228 | 1 | 0.3446 |
| statin_high_intensity | 23.99 | 0.1 | 0.5 | 0 | 0.09597 | 1 | 0.3392 |
| ctx310_angptl3_editing | 0 | 1 | 0 | 0 | 0 | 0.3 | 0.26 |
| verve102_base_editing | 0 | 1 | 0 | 0 | 0 | 0.3 | 0.26 |
| alirocumab | 12.56 | 0.2 | 0 | 0 | 0.05023 | 1 | 0.25 |
| colchicine | 20.87 | 0.1 | 0 | 0 | 0.08349 | 1 | 0.2367 |
| bempedoic_acid | 11.67 | 0.1 | 0 | 0 | 0.04666 | 1 | 0.2293 |

**Structural ceilings for a single agent:** driver control 1.0 · burden reversal
**0.60** · myocardial restoration **0.50** · risk normalisation **0.096** · proof 1.0.
Three of the five ceilings are below the 0.80 threshold *by construction*, so no
single agent can ever pass. That is a property of the criteria, not a discovery
about drugs.

## F.2 Study 2 — best-of-breed composite

**Question:** if we assemble the best agent per disease node, deduplicating by
mechanism, how close does the combination get?

15 disease nodes, target-overlap deduplication (three PCSK9 modalities count as one
mechanism), maximum 8 components, coverage by **maximum never summation**.
57 population × mode configurations, 114 criterion rows including invented variants.

| Credit mode | Criteria met | Risk normalisation | Best events avoided /1000 |
|---|---|---|---|
| strict | 1 of 5 | 0.819 | 134.3 |
| surrogate | 2 of 5 | 0.974 | 161.9 |
| exploratory | 2 of 5 | 1.000 | 266.9 |

Blocking everywhere: burden reversal 0.60, myocardial restoration 0.50, proof 0.

**Recurring core across populations:** genotype-guided canakinumab (19 populations),
VERVE-102 (19), smoking cessation (18), aspirin (18), icosapent ethyl (18).

**Proof is forced to zero** for every composite: the combination has never been
trialled, and component evidence does not transfer. Adding an invented component
changes no score and no curve — this is asserted by a test.

> **Later retracted in part:** the strict-mode "1 of 5" did not survive uncertainty
> propagation. See F.8.

## F.3 Study 3 — REVERT-1, an invented clearance architecture

**Status: INVENTED, QUARANTINED.** Eight assumed modules, no clinical existence.

| Module | Assumed strength | Role |
|---|---|---|

| module | strength | function |
|---|---|---|
| gate_specificity | 0.9 | lesion-targeted delivery gate; multiplies 6 of 8 modules |
| crystal_solubilisation | 0.6 | dissolve cholesterol crystals in situ |
| efflux_gain | 0.8 | cholesterol export capacity (ABCA1/ABCG1/ApoA-I) |
| efferocytosis_gain | 0.7 | corpse clearance by phagocytes |
| influx_block | 0.5 | block further lipid entry |
| resolution | 0.6 | pro-resolving inflammation programme |
| nucleation_block | 0.6 | prevent new crystal nucleation |
| plasma_driver_reduction | 0.7 | lower the circulating lipid driver |

**Ablation, 24 months (substrate removed):**

| Arm | Substrate removed | Necrotic core |
|---|---|---|
| **REVERT-1 full** | **43.7%** | **−46.5%** |
| efflux + efferocytosis | 18.3% | −23.7% |
| two-pronged published | 12.7% | +35.5% |
| efflux only | 5.3% | +35.5% |
| ungated | −4.7% | +12.7% |
| driver only | −14.2% | +35.5% |
| **efferocytosis only** | **−102.3%** | **+9.8%** |
| untreated | −60.9% | — |

**The central finding: efferocytosis alone is worse than doing nothing.** Engulfing
a corpse loads the macrophage with cholesterol; without matching export capacity
the phagocyte dies and becomes necrotic core itself. The AND gate coupling efflux to
efferocytosis is the highest-correlated parameter in the model (Spearman ρ = +0.555).

**Monte Carlo (400 draws): P(≥80% clearance) = 0%, median −4.9%.** The 43.7% figure
is a central-parameter scenario, not an expected outcome. The model's asymptotic
ceiling with perfect modules is 83.5%, set by a constant apoptosis source term.

**Design traps recorded so they are not repeated:** systemic CD47 blockade →
anaemia; systemic LXR agonism → hepatic lipogenesis; cleavage-resistant MerTK →
liver fibrosis; cyclodextrin ototoxicity caps systemic dosing; removing
macrocalcification can *destabilise* a plaque.

## F.4 Study 4 — back-calculating the required capability

**Question:** instead of asking "does this work?", ask "what would have to be true?"

Required vector: driver durability **0.80** · substrate removal **0.80** ·
LVEF gain **+20 pp** · hazard ratio **≤0.298**.

4,000 draws, seed 20261107. **Demonstrated feasible volume: exactly 0.**
Conceivable feasible volume: 3.686×10⁻⁴ (0.03%).

| Gap | Required improvement |
|---|---|
| LVEF gain | **×5.13** |
| Substrate removal | ×1.86 |
| Hazard ratio | ×1.21 |
| Driver durability | ×0.80 (already exceeded) |

Value of information (on worst-score, since VOI on a binary outcome is degenerate):
substrate 0.1839 > HR 0.1524 ≈ LVEF 0.1523 > durability 0.1409.

## F.5 Study 5 — attempting to close the three gaps

Three verdicts:

1. **Hazard ratio — already closed**, but the 0.298 target was population-specific.
   Recomputed per baseline: `stable_ccs` needs ≤0.298 and achieves 0.282 (closed);
   `secondary_post_mi` needs ≤0.257 and achieves 0.309 (16.4% short);
   `hefh_primary` needs ≤0.355, achieves 0.464 (23.1% short). Applying one universal
   HR was a category error.
2. **Substrate clearance — closable within precedent** via a *source-reduction*
   module (ER-stress/CHOP axis), because the plateau is a constant apoptosis source,
   not a payload limit. Precedent: CHOP deletion in Apoe⁻/ and Ldlr⁻/⁻ mice reduced
   plaque necrosis ~50% and lesional apoptosis ~35%. Bounded at **0.35 central /
   0.50 maximum** — never 1.0, because CHOP deletion does not abolish apoptosis.
3. **Remuscularisation — nominally closable**, later **retracted** (F.6).

### The four remuscularisation factors and their provenance

| factor | multiplier | evidence | level |
|---|---|---|---|
| dose_escalation | 1.75 | BioVAT-HF max dose 20 EHM units = 800M iPSC-CMs in 16 of 20 patients (NEJM 2026, doi:10.1056/NEJMoa2513525); human dose precedent, 1.75x escalation above it unvalidated | human dose precedent; escalation unvalidated |
| endothelial_co_delivery | 2 | EC+CM co-transplant raised primate EF 8.9 pp vs 1.9 pp CM alone; 10-fold larger grafts in rodents | primate + rodent |
| epicardial_stromal_co_delivery | 2.6 | epicardial co-delivery doubled CM proliferation, 2.6-fold graft size | rodent / construct |
| proliferation_induction | 2 | CCND2-overexpressing hiPSC-CMs replaced >50% of scar at 6 months; CHIR/FGF1 patches 4-fold engraftment | rodent / construct |

Baseline scar fraction 3.9/25 = 0.156; recoverable LVEF span 25 pp; so +20 pp
requires replacing 80% of scar — a **×5.13** improvement over baseline.

## F.6 Study 6 — step-by-step resolution, and a retraction

**Step 1 — the binding constraint on paper is evidence maturity, not biology.**
Counting populations that close as credit is relaxed: **strict 2, surrogate 6,
exploratory 11**; all 19 close somewhere. Genuine evidence-upgrade levers:

| Agent | Populations unlocked | Current tier → needed |
|---|---|---|
| **VERVE-102 base editing** | **17** | randomised surrogate → outcomes |
| Smoking cessation | 15 | observational → randomised |
| Genotype-guided canakinumab | 11 | observational → randomised |
| Retatrutide | 6 | surrogate → outcomes |
| Obicetrapib | 4 | surrogate → outcomes |
| LDL-lowering diet | 2 | surrogate → outcomes |

**One PCSK9-editing outcomes trial unlocks 17 of 19 populations and needs no new
biology.** This is the highest-leverage item found in the entire study.

Three agents were correctly *excluded* as evidence gaps (aspirin, bempedoic acid,
icosapent ethyl): they are already outcomes-tier and enter relaxed composites by
composition reordering, not by missing trials.

**Step 2 — clearance operating points** (minimum uniform potency multiplier for 80%
removal):

| Source reduction | 24 mo | 36 mo | 60 mo | 120 mo |
|---|---|---|---|---|
| 0.00 | unreachable | 1.627 | 1.403 | 1.332 |
| 0.20 | unreachable | 1.302 | 1.156 | 1.101 |
| 0.35 | 1.426 | 1.152 | 1.041 | 1.000 |
| 0.50 | **1.284** | 1.058 | **1.000** | 1.000 |

No single module suffices at 24 months: best is crystal solubilisation at 63.3%;
all modules perfect reaches only 82.7%.

> **This step was later retracted as a specification** — the 1.284× point is a
> nominal pass. See F.7.

**Step 3 — REMUSCULARISATION RETRACTED.** Treating each published multiplier as
Uniform(1, m) rather than a point value:

| Factor set | Nominal gain | P(+20 pp), floor 1.0 | P, floor 0.5 |
|---|---|---|---|
| endothelial + epicardial | +20.3 pp | **0.06%** | 0.00% |
| all four factors | +25.0 pp | **50.9%** | 26.5% |

**No subset is robust at 80% under any haircut.** The two-factor "pass" was an
artefact of multiplying point estimates: E[∏U(1,m)] ≪ ∏m.

## F.7 Study 7 — ATLAS-1 synthesis and its self-audit

**ATLAS-1** assembled the best *validated* part of each prior construct into six
arms — proven backbone, durable editing, gated clearance, source reduction,
remuscularisation, and one new idea: **temporal sequencing** (clear the lesion bed
first, then graft into it). Tested across 19 populations × 3 credit modes × 5
criteria = 57 configurations.

The new idea's link was an explicit assumption, `bed_coupling`, swept to zero:

| Bed coupling | P(+20 pp) | Gain vs unsequenced |
|---|---|---|
| 0.00 | 26.8% | **+0.0 pp** |
| 0.50 | 38.7% | +11.9 pp |
| 1.00 | 50.9% | +24.1 pp |

All six arms load-bearing under ablation. Initially reported **3 of 5** criteria.

**Then the self-audit retracted it.** The clearance operating point had been scored
at *central* module strengths while the graft was scored distributionally — the
same nominal-product error already retracted once:

| At the claimed operating point | Value |
|---|---|
| Nominal clearance | **80.01%** |
| Median under delivery uncertainty | **1.83%** |
| 5th–95th percentile | −38.3% to 41.7% |
| **P(clearance ≥ 80%)** | **0.00%** |

Second error: `bed_quality` derives from clearance, so the sequencing gain had been
scored against a bed that phase 1 usually fails to prepare. Scored jointly with a
paired contrast, decoupled scoring **overstated the graft by 20.3 pp**; median bed
quality was **0.023**, not 1.0.

**Corrected: ATLAS-1 reaches 2 of 5, not 3** — the invented arms add no criterion
that survives distributional scoring.

## F.8 Study 8 — auditing the evidence layer itself

The same question applied to the layer every construct inherits.
`composite.score_criteria` computed `risk_normalisation` from a single point
estimate, while `uncertainty.monte_carlo_regimen` already propagated LDL, hazard
ratio, baseline-rate and adherence intervals. The machinery existed; it was never
connected.

Only `risk_normalisation` is an effect estimate. Driver control is a structural
mechanism flag; burden and myocardial credit are capped below threshold; proof is
zero by construction.

| Credit mode | Nominal risk passes | Robust passes | Retracted | Best criteria |
|---|---|---|---|---|
| strict | 1 | **0** | 1 | **1 → 0** |
| surrogate | 7 | **2** | 4 | 2 → 2 |
| exploratory | 19 | **18** | 1 | 2 → 2 |

**Six passes retracted**, including `stable_ccs` strict (0.819 nominal → **11.7%**
probability) — which was the composite's *only* strict criterion. One case at 78.3%
is flagged **borderline** (within 2 Monte Carlo SE of the bar), not retracted.

**A median is not a probability either.** `lifetime_exposure` exploratory has
nominal 1.000 but a **bimodal** distribution: ~50% of draws at 0, ~50% at 1,
P(≥0.8) = 50%. Its median flips between 0.000 and 1.000 with sample size.

## F.9 Study 9 — ATLAS-2, fixing the parts the audits broke

**Diagnosis.** `delivered()` multiplies **six of eight** module strengths by one
`gate_specificity` draw — a series bottleneck with no redundancy — and delivery was
open loop. Both are architecture defects, not potency shortages.

**Three upgrades:**

| | Upgrade | Analogue |
|---|---|---|
| U1 | redundant gate and node paths, combined by max | distinct targeting routes |
| U2 | closed-loop titration within a hard cap | treat-to-target, standard of care in lipidology |
| U3 | verified source control by run-in | enrichment — **not needed at the solved horizon** |

| construct | horizon_months | probability_of_target | median_clearance_pct | eligible_fraction | unconditional_probability | robust | robust_with_margin |
|---|---|---|---|---|---|---|---|
| atlas1_baseline | 24 | 0 | 1.576 | 1 | 0 | False | False |
| atlas2_solved | 60 | 0.8612 | 81.16 | 1 | 0.8612 | True | True |

**Ablation at the solved point**

| removed_upgrade | probability_of_target | probability_loss_pp | active_in_construct | load_bearing |
|---|---|---|---|---|
| redundant_gate | 0.6125 | 24.87 | True | True |
| redundant_nodes | 0.4163 | 44.5 | True | True |
| titration | 0.06125 | 80 | True | True |
| verified_source | 0.8612 | 0 | False | False |

**Result: clearance P(target) 0.00% → 86.12%**, median 81.16%, robust with a 2-SE
margin (lower bound 83.68%), **eligible fraction 100%** so conditional and
unconditional are equal. Titration alone is worth **+72.8 pp**.

**The price:** potency cap 2.80× against a precedented 1.284× — a **2.18× precedent
gap** — and a 60-month horizon instead of 24. Bisection finds the minimum robust cap
at 2.590×, but that sits on the 80% bar within noise, so margin was bought instead.

**A correction to the earlier audit.** The "55.4% ceiling at any potency" was
**horizon-specific**. Perfect modules need source reduction 0.2230 at 24 months but
**none** at 36 months or beyond, because the structural ceiling rises above target
on its own. Time and source control are substitutes. That is why U3 is off: a run-in
at 60 months would only exclude patients who would have succeeded.

**The sequencing idea is recovered.** With the bed genuinely prepared (median
quality 1.000, previously 0.023):

| Bed coupling | Graft P(+20 pp) | Sequencing gain |
|---|---|---|
| 0.00 | 27.9% | **+0.0 pp** |
| 0.50 | 39.0% | +11.1 pp |
| 1.00 | 52.0% | **+24.1 pp** |

The idea was never structurally wrong; it was starved by a failing phase 1. The
graft still is not robust, so **remuscularisation is now the sole binding biology
problem**.

## F.9 ATLAS-3 — the variance fix applied to the graft (`atlas3`, INVENTED)

Full account: [`atlas3_graft.md`](atlas3_graft.md). Artefacts: `results/atlas3/` (10).

ATLAS-2 closed clearance without new potency, by attacking variance. The graft had
never been given the same treatment: it was still one irreversible, unmeasured shot.
ATLAS-3 applies retry redundancy plus a measured feedback loop, and charges every
extra procedure an adverse-event risk.

**The diagnosis.** Required multiplier **5.128×**; the product of all four published
maxima is **18.200×**, a **3.55×** surplus. Yet one graft into a *perfect* bed hits
target only **49.60%** of the time. Capability was never the shortfall — this is the
nominal-product trap again. A single hypothetical 8× mechanism scores **40.9%**, worse
than the four modest factors combined.

| Construct | P(target) | Median gain | Adverse risk | Risk-adjusted |
|---|---|---|---|---|
| Single graft (today) | **49.60%** | 19.92 pp | 5.00% | 47.12% |
| 2 attempts — **within precedent** | **82.50%** | 25.00 pp | 7.42% | 76.37% |
| ATLAS-3, 3 attempts | **93.10%** | 25.00 pp | 8.38% | 85.30% |

**The precedent breach is optional.** Two attempts already clear the bar (verified
across 12 seeds, P = 0.8216, range 0.8117–0.8300) and repeat dosing has precedent. So
the invention is not the number of procedures — it is the **closed loop**, a measured
engraftment signal that decides whether to repeat. The third attempt buys headroom
(82.5% → 93.1%) for **+0.95 pp** of adverse-event risk.

**Falsifier, paired.** Redundancy contributes **+38.88 pp** at zero failure
correlation and **exactly +0.00 pp** at correlation 1. Maximum tolerable correlation
**≈0.71** (±0.02 across seeds) — the quantity a real study must measure, and the one a
chronic poorly-vascularised scar would most likely violate. Feedback contributes
**exactly +0.00 pp** when it cannot learn, rising to +7.05 pp at full learning.

**Ablation — all four levers load-bearing:** retry redundancy **−43.50 pp**,
independence −26.60 pp, feedback −4.63 pp, measurement −4.63 pp.

**Both arms together:** clearance 83.50%, graft 92.45%, optimistic joint **77.20%**,
randomised proof **0.00%**. The joint is a product and therefore assumes arm
independence, which the measured bed coupling contradicts — labelled optimistic.

**VOI recomputed** (the discovery VOI predates ATLAS-2 and is stale): substrate
removal **+15.25 pp** still ranks above the graft **+6.30 pp**. The ranking held, but
the reason changed — clearance now has an architecture that simply delivers a lower
probability.

**Verdict.** The missing step is a remuscularisation modality that can be **repeated
and measured**, not one that is more potent. ATLAS-3 specifies the experiment; it is
not a therapy. Cure count remains **zero**.

## F.10 ATLAS-4 — searching the literature, and three corrections (`atlas4`)

Full account: [`atlas4_precedent_search.md`](atlas4_precedent_search.md).
Artefacts: `results/atlas4/` (11). Registry: `src/ihd_solutions/precedent.py`, 12 sources.

The question asked was whether anyone had **already done** the missing step. Most of
them have. **CATO (NCT06145035)** is a Phase IIA randomised double-blind placebo-controlled
trial randomising 60 patients with ischaemic cardiomyopathy to four placebo, one, or
**four** intravenous UC-MSC doses two months apart, primary endpoint ΔLVEF by cardiac MRI
at 12 months; 47 were enrolled as of 2025-11-30. HOPE-2 gave four doses in humans in
2022. Reporter-gene PET already predicts later cardiac function in animals
(**ρ=0.76** echo, **ρ=0.80** MRI).

That forced **three corrections to ATLAS-3**:

| ATLAS-3 said | Actually | Source |
|---|---|---|
| Precedented dose limit 2; three attempts "beyond precedent" | Limit is **4**; three was never beyond precedent | CATO, HOPE-2 |
| Retries combine as best-of-k | They combine **cumulatively**; splitting a fixed total dose beats one big dose (+13.0 vs +4.7 pp) | Circ Res 2016; JAHA 2018; Cardiovasc Res 2024 (pig) |
| LVEF gain is scar replacement | Cells **do not regenerate myocardium**; effect is paracrine/anti-fibrotic | same, explicit |

**Corrected result.** Cumulative beats best-of-k at every dose count (+8.85 pp at four
doses). Tempered by the human responder rate — only **10 of 30** patients responded to a
second dose in the one blinded repeat trial — four doses reach **82.97%**, robust with
margin, **within precedent**, at 7.98% cumulative adverse risk (76.35% risk-adjusted).
A 0.7 translation floor drops it to 76.40%, no longer robust.

**A consequence.** ATLAS-3's "measure independence above ≈0.71" was operator-specific:
best-of-k retains only **35.8%** of its gain at full correlation, cumulative retains
**90.6%**. The most emphasised unmeasured quantity mattered less than claimed.

**Measurement falsifier.** Iron-oxide MRI and acute-retention imaging gain **exactly
0.00 pp** — iron signal survives cell death and retention does not predict engraftment.
Only a viability-specific signal can close a loop. Reporter-gene PET gains just
+0.40 pp under the corrected operator, so feedback matters far less than ATLAS-3 implied.

**All 360 design combinations tested.** 176 robust, 160 robust with margin, 288 within
precedent, **48 reading exactly 100%** — every one a saturation artefact — leaving
**8 admissible**. Best admissible: cumulative, four doses, ρ=0.50, full learning,
human-tempered: **89.60%**.

**On reaching 100% everywhere.** Measured catalogue ceilings: durable driver control
**1.000** and randomised proof **1.000** already; burden reversal 0.600, myocardial
restoration 0.500, risk normalisation 0.061. The last three are capped by explicit
evidence rules, not by an unexplored design space, and proof attaches to the exact
regimen randomised so it never transfers to an untested combination. **100% across all
five is unreachable by searching harder.** The one design that does read 100% is the
rodent-optimistic one, and it stays at 99.98% even under a harsh translation penalty —
the certainty comes from the responder assumption, not the biology.

**What is genuinely still missing:** response-guided retry timing, proposed in print in
2019 and never tested, and a viability-specific engraftment signal in humans.

## F.11 Solutions to the six open gaps (`gap-solutions`)

Full account: [`gap_solutions.md`](gap_solutions.md). Artefacts: `results/gap_solutions/` (10).

Each of the six gaps listed as "nobody has run this" was searched individually. **Three
of the six claims were wrong or overstated.**

| Gap | Verdict | Basis |
|---|---|---|
| 1 Nothing regrows heart muscle | **FALSE** | Jebran, Nature 2025;639:503-511: explanted human heart 3 months post-EHM showed retained, vascularised, viable graft. BioVAT-HF, NEJM 2026: +4.5 mm wall, +3.9 pp LVEF at 3 months in 16 patients at max dose |
| 2 Response-guided retry | **Stands** | Proposed 2019, never tested; selection signal already published (Vrtovec responders' pre-dose LVEF 34.9 vs 40.0, p=0.02) |
| 3 Viability signal in humans | **Overstated** | hNIS already imaged in patients: MV-NIS 123I SPECT/CT positive in 8/32 myeloma and 3 ovarian-cancer patients across repeat cycles; human-derived, non-immunogenic, beat HSV1-tk for cardiac contrast |
| 4 Repeat dosing awaits CATO | **Stands, but underpowered** | 20/arm, 2 comparisons; at SD 6 pp detects only 5.85 pp and has **25.5% power** for a 3 pp increment |
| 5 2.18x potency gap | **Reduced to 1.87x** | Redundancy substitutes: cap 2.40x with 3 node + 3 gate paths gives P 89.25% |
| 6 PCSK9 durability assumed | **Measured** | VERVE-102 Heart-2, NEJM 2026: -88% PCSK9, -62% LDL-C, -78 mg/dL, durable 18 months through hepatocyte turnover |

**Gap 1 quantified.** The gap is magnitude, not existence. Against the +20 pp criterion
the 3-month effect needs **5.13x** the demonstrated dose (103 EHM units, 4.10 billion
cardiomyocytes); the 17-month effect needs **2.90x**. Maturation alone buys **1.77x** and
closes **43.5%** of the mass requirement. Critically, wall thickness **fell 1.6 mm** while
LVEF **rose 3.0 pp** over that interval - function improved as mass declined, pointing at
maturation and mechanical entrainment rather than bulk. Mass scaling may be the wrong
lever.

**Gap 6 powered.** 78 mg/dL = 2.02 mmol/L. On the CTT slope that implies HR 0.606 and
**125 events** (~1,389 patients at 3%/yr over 3 years); at half slope, 500 events. Design
arithmetic under an imported statin slope, explicitly **not** a claim that editing
delivers it.

**Two accounting traps found and fixed.** The first gap-5 pass reported a 1.56x gap at
cap 2.00x, but that point needs 4 node paths and `MAX_PRECEDENTED_PATHS = 3`, so the
redundancy was itself beyond precedent - paying for a potency breach with an unflagged
redundancy breach is not a saving. Selection now requires both axes precedented, giving
the honest 1.87x. Separately, selecting the fewest paths at the cheapest cap picked a
point at 82.88%, within 2 SE of the bar; selection now takes the most headroom at the
cheapest cap.

Three gaps need no new biology: gap 2 is an entry rule, gap 3 a manufacturing change,
gap 6 a trial design. No cure is declared.

---

## F.12 Mechanism recombination — ATLAS-5 (`recombine`)

Full account: [`mechanism_recombination.md`](mechanism_recombination.md).
Artefacts: `results/recombine/` (13).

The previous arms were ablated, so the load-bearing part of each was already known. This
turn decomposed them into **eight transferable mechanisms**, measured each alone and in
every pair, and composed only the parts that earned their cost.

**Three rules make the composition honest.** Interactions are measured so two mechanisms
doing the same job are not both credited. Two of the eight parts are *constraints* that
lower the reported probability, and the search may not switch them off — otherwise it
would buy a better answer by deleting the honesty that produced it. And every precedent
axis is charged separately.

| mechanism | standalone | leave-one-out | verdict |
|---|---|---|---|
| `maturation` (human EHM series) | **+40.04 pp** | +25.12 pp | kept |
| `retry_redundancy` | +20.52 pp | +3.20 pp | kept |
| `cumulative_operator` | 0.00 pp | +0.08 pp | kept (interaction only) |
| `response_guided_selection` | +0.08 pp | **−3.52 pp** | **dropped** |
| `feedback` | 0.00 pp | 0.00 pp | **dropped — inert** |
| `viability_reporter` | 0.00 pp | 0.00 pp | **dropped — inert, costs a genetic mod** |

**The overlap that mattered.** Maturation and retry redundancy both push toward the same
recoverable-scar ceiling. Credited additively they would claim **61.20 pp**; the measured
joint gain is **47.12 pp**. **14.08 pp of that sum was the same biology counted twice.**
Conversely `retry_redundancy` + `cumulative_operator` is genuinely **synergistic
(+10.48 pp)**: adding increments differs from replacing them only when there is more than
one dose.

**Response-guided selection is a trial-design lever, not a therapy.** Vrtovec's screen
concentrates responders up to 2.34x, but it *excludes patients*. Over the whole population,
using the central swept assumption rather than the maximum: 42.5% pass x 58.2% respond =
**24.7% of all patients**, against 33.3% unscreened. Fewer people benefit, so the effect is
**−3.52 pp** and the mechanism is excluded. Enrichment shrinks the trial needed to detect
an effect; it does not increase the number who get one.

**Naive stacking is worse than selecting.** All 64 subsets, constraints held on:

| design | mechanisms | P(target) | risk-adjusted |
|---|---|---|---|
| everything on | 6 | 94.40% | 89.48% |
| **ATLAS-5 (selected)** | **3** | **97.92%** | **92.53%** |
| constraints only | 0 | 51.68% | 49.10% |

Switching everything on is **3.52 pp worse**. Recombination is a selection problem, not an
accumulation problem.

**ATLAS-5** = `retry_redundancy` + `cumulative_operator` + `maturation`, four doses
available, **no feedback loop and no reporter cassette — zero genetic modifications**.
Mean doses *actually used* is **1.10**, because maturation makes the first dose succeed
often enough that retries are rarely invoked.

| comparison | value |
|---|---|
| best single mechanism | 90.60% |
| gain over best single | **+6.55 pp** |
| ATLAS-4 best admissible | 89.60% |
| **gain over ATLAS-4** | **+7.55 pp** |
| price of admissibility | 0.05 pp (sampling noise) |

The reporter buys 0.05 pp — noise — for a genetic modification with no cardiac precedent,
so it is dropped. The graft arm is now **inside precedent on every axis it controls**; the
clearance arm's 1.87x potency breach is carried forward and **not forgiven**.

**Robustness is 50%, not 97.92%.** The maturation ratio is measured in EHM patches while
the factor stack comes from animal work; multiplying them is an unvalidated transfer. Over
36 scenarios (multiplier x translation penalty x doses) only **18 are robust with margin**.
And a **declared blind spot**: increments clamp at the ceiling, so above **1.25x** the model
cannot resolve a larger multiplier at all — it cannot represent the full 1.77x. The
maturation credit is what the model can resolve, not the full human ratio.

Randomised proof for the composition: **0.00%**. Zero saturation artefacts in 64
compositions. No cure is declared.

---

# PART G — THE HONESTY LEDGER

Every retraction and correction made during this study, in order. This section
exists because a report that only shows final numbers hides how they were reached.

| # | Claim | What was wrong | Corrected to |
|---|---|---|---|
| 1 | Two-factor graft gives +20.3 pp | Multiplied point estimates; E[∏U] ≪ ∏m | **P = 0.06%** |
| 2 | Universal HR target 0.298 | Category error; target depends on baseline | Recomputed per population |
| 3 | Clearance solved at 1.284× | Nominal pass at central strengths | **P = 0.00%**, median 1.83% |
| 4 | ATLAS-1 reaches 3 of 5 | Inherited the nominal clearance pass | **2 of 5** |
| 5 | Sequencing worth +24.1 pp | Scored against a bed phase 1 rarely prepares | Overstated **20.3 pp** |
| 6 | Composite strict 1 of 5 | Risk score was a point estimate | **0 of 5** (P = 11.7%) |
| 7 | 55.4% ceiling at any potency | Fixed the horizon at 24 months | Vanishes at **≥36 months** |
| 8 | Docking ranks affinity | All ΔG positive; Kd clamped at 1 M | **Non-usable, flagged** |
| 9 | Screening drug-likeness valid | Biologics given small-molecule MWs | **4 rows invalid** |
| 10 | Graft needs 3 attempts, beyond precedent | 2-SE margin at small N rejected a design genuinely above the bar | **2 attempts suffice**, within precedent |
| 11 | Feedback raises each factor's floor | Product of minima exceeded the requirement → P = 100% | **One bounded aggregate uplift**; no P = 100% reachable |
| 12 | Correlated retries gain +0.1 pp | Independent streams; noise read as signal | **Paired contrast, exactly 0.00 pp** |
| 13 | Two doses is the precedented limit | Never checked the literature | **Four** (CATO, HOPE-2); three was never beyond precedent |
| 14 | Retries combine as best-of-k | Published repeat-dose studies show cumulative increments | **Cumulative operator**; +8.85 pp at four doses |
| 15 | Graft LVEF gain is regenerated myocardium | These cells do not engraft or differentiate | **Paracrine**; ceiling is functional, not new tissue |
| 16 | Independence above 0.71 must be measured | Property of the wrong operator | Cumulative retains **90.6%** of gain at full correlation |
| 17 | Grid search's best admissible point | Filter admitted rodent optimism sitting just under certainty | Filter now caps repeat response at the **human** rate; 32 → 8 admissible |
| 18 | Nothing regrows heart muscle in a human | Never searched the iPSC-CM patch trials | **Proven** histologically in an explanted human heart (Nature 2025); gap is magnitude |
| 19 | A viability reporter has never been imaged in humans | Confused "not in a cardiac graft" with "not in humans" | hNIS imaged in **myeloma and ovarian-cancer patients** |
| 20 | PCSK9-editing durability is assumed | Pre-dated the Heart-2 readout | **Measured** to 18 months through hepatocyte turnover |
| 21 | Redundancy buys the potency gap down to 1.56x | The cheap point used 4 node paths, itself beyond precedent | Both axes must be precedented: **1.87x** |
| 22 | Interaction of −14.15 pp is "additive" | Verdict compared a probability difference against percentage-point cutoffs | Tolerance in probability units; **3 overlaps, 1 synergy** |
| 23 | Response-guided selection is a therapy benefit | Reported the screened-in subset, hiding the patients the screen excludes | Whole population: **−3.52 pp**; a trial-design lever only |
| 24 | Composing every mechanism is the best composite | A negative-effect part drags the stack down | Selection beats stacking by **3.52 pp** |
| 25 | The strongest composite needs an unprecedented modification | Conflated the price of admissibility with the price of cheapness | It needed **none**; admissibility costs 0.05 pp (noise) |
| 26 | Dose escalation is a construct-scale projection | Pre-dated BioVAT-HF's reported maximum dose | **800M cells given to 16 of 20 patients**; the 1.75x above it stays unvalidated |

**Implementation bugs found and fixed:** run-in pinned responders to the threshold
instead of keeping their achieved value; gate threshold computed at 24 months but
applied at 60; divide-by-zero in the asymptotic-ceiling helper at zero source
reduction; ablation labelled durable editing "inert" because it only measured
clearance and graft mass, not durability; pytest collected a library function named
`test_*`.

**Rules enforced by tests, not by prose:**

- Invented arms earn **zero** evidence credit.
- An untested combination earns **zero** randomised proof.
- Coverage combines by **maximum, never summation**.
- Animal or in-vitro data never earns human proof or biomarker transport.
- Partial regression credit ≤0.6; myocardial credit ≤0.5.
- Single-agent risk normalisation ≤0.4 × RRR.
- Imaging anchors never enter event credit.
- Neutral trials are sampled around the null, not their point estimate.
- Adding an invented component must change **no** score and **no** curve.
- Invented artefacts live outside `results/solutions/` and outside its manifest.

---

# PART H — REPRODUCIBILITY AND QUALITY GATE

| Check | Result |
|---|---|
| Test suite | **1060 passed** (917 s) |
| Coverage of `ihd_solutions` | **98.74%** (`recombine.py` **100%**, `atlas4.py` 99%, `gap_solutions.py` 99%) |
| Lint (`ruff`) | clean, 66 files formatted |
| Registry validation | **47 therapies / 48 anchors / 19 populations** |
| Evidence artefacts hash-verified | **33**, **0 mismatches** |
| Cures declared | **0** |
| Manifest scope | `scenario simulation; not clinical prediction` |

Every command that regenerates a results tree:

```bash
python run_solutions.py run             # evidence-graded pipeline  -> results/solutions/
python run_solutions.py hypothesis      # REVERT-1 (invented)       -> results/hypotheses/
python run_solutions.py discover        # capability back-calc      -> results/discovery/
python run_solutions.py composite       # best-of-breed composite   -> results/composite/
python run_solutions.py gaps            # gap-closure attempt       -> results/gaps/
python run_solutions.py resolve         # step-by-step resolution   -> results/resolution/
python run_solutions.py atlas           # ATLAS-1 synthesis         -> results/atlas/
python run_solutions.py audit           # ATLAS-1 self-audit        -> results/audit/
python run_solutions.py audit-evidence  # evidence-layer audit      -> results/evidence_audit/
python run_solutions.py upgrade         # ATLAS-2 (the fix)         -> results/upgrade/
python run_solutions.py atlas3          # ATLAS-3 (graft retry+loop) -> results/atlas3/
python run_solutions.py atlas4          # ATLAS-4 (precedent + grid) -> results/atlas4/
python run_solutions.py gap-solutions   # designs for the 6 gaps    -> results/gap_solutions/
python run_solutions.py recombine       # ATLAS-5 mechanism composition -> results/recombine/
python run_solutions.py validate        # registries only
```

Monte Carlo is seeded (20261107 for the cure programme, 20261106 for the regimen
uncertainty runs). Common random numbers are used across potency sweeps, and the
coupled-vs-decoupled contrast is paired so its zero-coupling value is exactly 0, as
is the ATLAS-3 retry-independence falsifier. ATLAS-3's specification searches floor
their draw count at 4000 so a specification cannot move with the sampling budget.

---

# PART I — WHAT REMAINS OPEN

## I.1 Ranked open problems

1. **Remuscularisation — the open item is narrower than first thought.** ATLAS-4's
   literature search found repeat dosing is already in human trial (CATO, four doses)
   and engraftment imaging already predicts function in animals. Corrected for the
   cumulative operator and tempered by the human responder rate, four doses reach
   **82.97%** within precedent. What is genuinely untested is **response-guided retry
   timing** (proposed 2019) and a **viability-specific engraftment signal in humans**.
   Note the mechanism is paracrine, not regeneration, and the one blinded human
   repeat-dose trial was **null**.
2. **The ATLAS-2 potency assumption** — 2.18× beyond precedent across every
   clearance module simultaneously. This is the specification to demonstrate, not
   a result.
3. **Randomised proof** — zero for every combination in this study. No simulation
   can supply it.
4. **The legacy chemistry tables** — docking and screening need regeneration with
   a validated scoring function and correct structures before any use.

## I.2 The single highest-leverage real-world action

**A PCSK9-editing outcomes trial.** It unlocks **17 of 19 populations** in the
credit model, requires **no new biology**, and replaces assumed durability and
surrogate credit with a randomised effect estimate — precisely the deficiency all
three audits exposed. Nothing in the invented constructs changes that recommendation.

## I.3 Final statement

This study did not find a cure for ischaemic heart disease and does not claim one.
What it produced is: a validated structural explanation for why small molecules
failed against PCSK9; an evidence-graded map of which agents help which populations
and where the evidence runs out; a quantified specification of what a cure would
require; an invented architecture that makes plaque clearance robust in simulation
at a stated price; a second that does the same for the graft and shows the missing
step is *repeatability plus measurement*, not potency; a third composite assembled from
the measured load-bearing parts of the first two, which is stronger and simpler than
either and refuses 14.08 pp of double-counted benefit; and a documented record of
twenty-six claims that were tested and withdrawn. **The retractions are part of the
result.**

---

# APPENDIX — COMPOUND AND AGENT INDEX

A single consolidated index of every chemical entity named anywhere in this study
is written alongside this report as `FULL_REPORT_chemical_index.csv`, with one row
per entity and a `source` column identifying which register it came from:

- `therapy_catalog` — the 47 evidence-graded therapeutic agents
- `virtual_screening` — the 28 screening-library compounds (**heuristic, partly invalid**)
- `docking_de_novo` — the 8 docked/designed compounds (**non-usable scores**)
- `revert1_module` — the 8 assumed clearance modules (**invented**)
- `remuscularisation_factor` — the 4 graft enhancement factors (**animal/construct precedent**)

Molecular weights are reproduced as recorded, **not** corrected; the `qa_flag`
column marks the rows this report identified as unreliable.
