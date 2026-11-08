# IHD solution simulation — methods, code map and test strategy

**Scope.** This package simulates *scenarios* for candidate interventions against
ischaemic heart disease (IHD). It is **not** a clinical trial, **not** a
patient-level prediction model, **not** a treatment recommendation, and **not** a
cure. Every number it produces is a transparent arithmetic consequence of
published effect sizes plus explicitly stated modelling assumptions. Where the
evidence does not support an event benefit, the code returns **zero credit** and
says why.

Everything described here lives in `src/ihd_solutions/`, `data/solutions/` and
`tests/`. The original numbered scripts (`src/01_*.py` … `src/05_*.py`) and the
older `docs/*.md` narratives are untouched by this work and are **not** validated
by this test suite.

---

## 1. Quick start

```bash
python -m venv .venv && . .venv/bin/activate
pip install -r requirements-dev.txt

python run_solutions.py validate                 # schema + consistency checks
python run_solutions.py therapies                # catalogue with mechanisms
python run_solutions.py anchors                  # every effect estimate + URL
python run_solutions.py populations              # 13 scenario populations
python run_solutions.py regimens --max-optional 2

# one regimen, one population, one credit mode
python run_solutions.py evaluate \
    --therapies statin_high_intensity aspirin_secondary_prevention colchicine \
    --population stable_ccs --mode strict --json

# full run: 6 populations x 337 regimens x 3 modes -> results/solutions/
python run_solutions.py simulate

pytest -q                                        # 647 tests
pytest --cov=ihd_solutions --cov-report=term-missing
ruff check src tests run_solutions.py && ruff format --check src tests run_solutions.py
```

`python -m ihd_solutions …` is equivalent to `python run_solutions.py …`.
A default run takes ~10 s and writes 33 hashed artefacts (~3 MB).

---

## 2. What is modelled

Nine intervention pathways, **47 therapies, 19 populations, 48 trial anchors**
after the 2026 cure search (see [`cure_search_2026.md`](cure_search_2026.md)).

| Pathway | Therapies | Event evidence handled by |
| --- | --- | --- |
| LDL / ApoB lowering | high- and moderate-intensity statin, ezetimibe, bempedoic acid, evolocumab, alirocumab, inclisiran, enlicitide, AZD0780 | CTT per-mmol/L relationship + agent-specific outcome trials |
| Permanent lipid lowering | VERVE-102 *PCSK9* base editing | knockdown → LDL translation, no MACE credit |
| Lp(a) lowering | pelacarsen, olpasiran, muvalaplin | surrogate only; outcomes pending |
| Residual inflammation | colchicine, canakinumab, ziltivekimab, methotrexate, NLRP3 inhibitor | population-specific trials, including three neutral ones |
| Metabolic | semaglutide 2.4 mg, empagliflozin, SGLT2 for HF | SELECT / EMPA-REG; HF composite kept separate |
| Blood pressure | guideline-style intensification | BPLTTC per-5 mmHg relationship |
| Antithrombotic | aspirin, icosapent ethyl | ATT, REDUCE-IT (comparator caveat recorded) |
| Lifestyle | smoking cessation, exercise/cardiac rehab, diet | observational tiers, flagged |
| Revascularisation & regeneration | PCI for stable CAD, CABG for ischaemic CM, iPSC-CM sheet | ISCHEMIA (neutral), STICH (mortality only), RiHEART (functional only) |

### Therapy catalogue (values read from `data/solutions/therapy_catalog.csv`)

The table below covers the 30 core therapies. The 17 candidates added by the 2026
cure search (CETP inhibition, ANGPTL3 editing, APOC3 silencing, next-generation
Lp(a) agents, incretins, genotype-guided IL-1β blockade, engineered heart muscle,
regeneration gene therapy, cell therapy, devices, vaccines and CAR-Treg) are
tabulated with their sources in [`cure_search_2026.md`](cure_search_2026.md).

| Therapy | Class | LDL ↓ % | SBP ↓ mmHg | hsCRP ↓ % | Lp(a) ↓ % | Anchors |
| --- | --- | --- | --- | --- | --- | --- |
| `statin_high_intensity` | lipid | 50 | 0 | 15 | 0 | CTT-STATIN |
| `statin_moderate_intensity` | lipid | 30 | 0 | 10 | 0 | — |
| `ezetimibe` | lipid | 18 | 0 | 0 | 0 | IMPROVE-IT |
| `bempedoic_acid` | lipid | 21 | 0 | 0 | 0 | CLEAR-OUTCOMES |
| `evolocumab` | lipid | 59 | 0 | 0 | 20 | FOURIER |
| `alirocumab` | lipid | 55 | 0 | 0 | 15 | ODYSSEY-OUTCOMES |
| `inclisiran` | lipid | 50 | 0 | 0 | 20 | ORION-11 (LDL only) |
| `enlicitide_oral` | lipid | 56 | 0 | 0 | 28 | CORALREEF-LIPIDS (LDL only) |
| `azd0780_oral` | lipid | 51 | 0 | 0 | 0 | PURSUIT (LDL only) |
| `verve102_base_editing` | gene editing | 62 | 0 | 0 | 0 | HEART-2 (phase 1) |
| `pelacarsen` | Lp(a) | 0 | 0 | 0 | 80 | LPA-HORIZON (pending) |
| `olpasiran` | Lp(a) | 0 | 0 | 0 | 95 | OCEAN-A (pending) |
| `muvalaplin_oral` | Lp(a) | 0 | 0 | 0 | 65 | MUV-Phase2 (surrogate) |
| `colchicine` | anti-inflammatory | 0 | 0 | 25 | 0 | LODOCO2 / COLCOT / CLEAR-SYNERGY |
| `canakinumab` | anti-inflammatory | 0 | 0 | 37 | 0 | CANTOS-150 |
| `ziltivekimab` | anti-inflammatory | 0 | 0 | 60 | 0 | ZEUS (neutral) |
| `methotrexate_low_dose` | anti-inflammatory | 0 | 0 | 0 | 0 | CIRT (neutral) |
| `nlrp3_inhibitor` | anti-inflammatory | 0 | 0 | 40 | 0 | — (preclinical) |
| `semaglutide_2_4mg` | metabolic | 0 | 5 | 35 | 0 | SELECT |
| `empagliflozin` | metabolic | 0 | 4 | 0 | 0 | EMPA-REG |
| `sglt2_inhibitor_hf` | metabolic | 0 | 3 | 0 | 0 | DAPA-HF (HF endpoint) |
| `bp_intensification` | blood pressure | 0 | 10 | 0 | 0 | BPLTTC-SBP |
| `aspirin_secondary_prevention` | antithrombotic | 0 | 0 | 0 | 0 | ATT-ASPIRIN |
| `icosapent_ethyl` | antithrombotic | 0 | 0 | 20 | 0 | REDUCE-IT |
| `smoking_cessation` | lifestyle | 0 | 5 | 30 | 0 | QUIT-POSTMI |
| `exercise_cardiac_rehab` | lifestyle | 3 | 4 | 20 | 0 | EXERCISE-REHAB |
| `diet_ldl_lowering` | lifestyle | 15 | 3 | 10 | 0 | — |
| `pci_stable_cad` | revascularisation | 0 | 0 | 0 | 0 | ISCHEMIA (neutral, angina only) |
| `cabg_ischemic_cm` | revascularisation | 0 | 0 | 0 | 0 | STICH-CABG (all-cause death) |
| `ipsc_cm_sheet` | regenerative | 0 | 0 | 0 | 0 | RIHEART-PIVOTAL (functional) |

Biomarker percentages are *relative* changes and stack multiplicatively, never
additively, with hard floors (LDL-C 10 mg/dL, ApoB 15 mg/dL, Lp(a) 2 nmol/L,
hsCRP 0.2 mg/L, SBP 80 mmHg).

---

## 3. Evidence handling — the part that decides the answer

### 3.1 Three credit modes

| Mode | Admits | Purpose |
| --- | --- | --- |
| `strict` (headline) | meta-analyses and randomised hard-outcome trials whose population matches or is explicitly declared transportable | what the evidence supports today |
| `surrogate` | additionally transports the CTT LDL and BPLTTC pressure relationships to agents that only have biomarker data | "if the LDL hypothesis holds for this modality" |
| `exploratory` | additionally admits observational evidence and CV/all-cause mortality endpoints | upper bound of plausible benefit |

The mode is a first-class argument everywhere; the ranking table reports all
three side by side so the reader can see exactly how much of a result depends on
an assumption rather than a trial.

### 3.2 Rules that are enforced by tests, not by prose

1. **No surrogate may create events.** hsCRP and Lp(a) changes never produce MACE
   credit in any mode. `EndpointClass.SURROGATE`, `FUNCTIONAL`, `PENDING` and
   `HF_COMPOSITE` can never be relabelled as MACE.
2. **Neutral trials stay neutral.** ZEUS (0.99), CIRT (1.01), CLEAR SYNERGY (0.99)
   and ISCHEMIA (0.93, CI 0.83–1.04) have intervals spanning 1.0. The default
   `RiskParameters.neutral_ci_policy="null_effect"` applies HR = 1.0 exactly, so a
   neutral drug changes nothing. `neutral_ci_policy="point_estimate"` is available
   as a sensitivity analysis and is labelled as such. Neutrality is stored as
   *reported* (`TrialAnchor.neutral_by_report`) so Monte Carlo resampling cannot
   accidentally turn a neutral trial into a positive one.
3. **No double counting.** A therapy with its own randomised MACE anchor
   contributes that anchor (which already contains its lipid effect) and is
   excluded from the per-unit LDL/pressure pathway. `CreditPartition` keeps
   `fixed`, `per_unit` and `uncredited` disjoint, and the pathway participant lists
   are mode-gated — this is what prevents an unproven agent from leaking into the
   strict curve through the biomarker route.
4. **Population matters.** Colchicine resolves to LoDoCo2 in chronic coronary
   syndrome, COLCOT in early post-MI, and CLEAR SYNERGY after modern PCI. Statins
   are deliberately **not** transportable to ischaemic cardiomyopathy with heart
   failure, so a statin scores 0 MACE there.
5. **Every decision is auditable.** `credit_decisions.csv` records therapy,
   population, mode, credited flag, anchor, HR, reason string and warnings
   (`neutral_confidence_interval`, `transported_population`,
   `observational_confounding`, …).
6. **Guards.** Benefit accrues at half rate in year 1 (CTT), combined RRR is
   capped at 85 %, and the baseline curve is integrated with the same trapezoid
   routine as the treated curve so an all-uncredited regimen returns exactly 0.0
   events avoided rather than floating-point dust.

### 3.3 Quantitative anchors

| Anchor | Effect | Source |
| --- | --- | --- |
| CTT 2010 | RR 0.78 per 1 mmol/L LDL-C | 170 000 patients, 26 trials |
| BPLTTC 2021 | HR 0.90 per 5 mmHg SBP | 344 716 patients, 48 trials |
| IMPROVE-IT | 0.936 (0.897–0.978) | ezetimibe add-on |
| CLEAR Outcomes | 0.87 (0.79–0.96) | bempedoic acid, statin-intolerant |
| FOURIER / ODYSSEY | 0.85 (0.78–0.93) each | PCSK9 antibodies |
| SELECT | 0.80 (0.72–0.90) | semaglutide, ASCVD + obesity, no diabetes |
| CANTOS 150 mg | 0.85 (0.74–0.98) | IL-1β blockade |
| LoDoCo2 / COLCOT | 0.69 (0.57–0.83) / 0.77 (0.61–0.96) | colchicine |
| CLEAR SYNERGY / ZEUS / CIRT | 0.99 / 0.99 / 1.01 | neutral controls |
| REDUCE-IT | 0.75 (0.68–0.83) | mineral-oil comparator caveat recorded |
| EMPA-REG | 0.86 (0.74–0.99) | T2D + ASCVD |
| ATT aspirin | 0.81 | secondary prevention |
| ISCHEMIA | 0.93 (0.83–1.04) | angina relief only |
| STICH 10 y | 0.84 (0.73–0.97) | all-cause death, not MACE |
| HEART-2 | 88 % PCSK9 knockdown → 62 % LDL-C | n = 35, ≤18 months, phase 1 |
| RiHEART | LVEF +≥5 % in 2/8, peak VO₂ +≥10 % in 4/7 | paracrine, conditional approval |

Full URLs, participant counts, follow-up and endpoint classes are in
`data/solutions/trial_anchors.csv`; `python run_solutions.py anchors` prints them.

---

## 4. Model layers

```
catalog.py ─┐
population ─┼→ pharmacology.py → lipids.py → risk.py → portfolio.py → pipeline.py
evidence.py ┘        ↑                          ↑           ↑
                 uncertainty.py ────────────────┘      editing.py
                 inflammation.py   plaque.py           regeneration.py
```

* **`pharmacology.py`** — exposure profiles: onset, plateau, dosing interval with
  trough fraction (so a 14-day antibody produces a sawtooth, not a step), waning
  half-life for one-time edits, and adherence as a taking-fraction.
* **`lipids.py`** — multiplicative biomarker stacking with floors and per-therapy
  `components`, which is what the risk layer consumes for pathway credit.
* **`risk.py`** — piecewise hazard → cumulative hazard (trapezoid) → incidence
  `1 − exp(−H)`; events avoided per 1000 over the horizon; NNT = `1000 / benefit`
  (∞ when nothing is credited); `notes` explain every non-credit.
* **`inflammation.py`** — 6-state ODE (oxLDL → NLRP3 → IL-1β → IL-6 → hsCRP, plus
  MMP) integrated with RK4 and cross-checked against the analytic steady state.
  Methotrexate has no targets, reproducing CIRT as an exact negative control.
* **`plaque.py`** — lipid-core growth, cap thinning and Laplace hoop stress with
  65 µm / 300 kPa proxy thresholds. Timeline rows after a crossing are explicitly
  flagged `post_threshold_counterfactual`; this is a teaching proxy, **not** a
  rupture prediction.
* **`editing.py`** — PCSK9 knockdown → LDL-C translation calibrated to HEART-2,
  durability scenarios (permanent vs 5-year waning), bystander enumeration inside
  the editor window, and `off_target_assessment()` which returns
  `numeric_off_target_score = None`: DSB-based assays miss adenine-base-editor
  off-targets, so safety requires a dedicated experiment.
* **`regeneration.py`** — beta-binomial posteriors for the 8-patient RiHEART
  pivotal data, `event_credit(mode) == 0.0` in every mode, and a
  conditional/time-limited surveillance plan (75 confirmation cases).
* **`portfolio.py`** — regimen enumeration with exclusivity groups (all
  protein-level PCSK9 agents share one slot; the gene edit has its own;
  the two SGLT2 entries share one), contraindication/eligibility filtering, an
  independent combinatorial count check, burden scoring and a Pareto front.
* **`uncertainty.py`** — Monte Carlo over efficacy CIs, hazard-ratio CIs, a Beta
  adherence prior and the baseline-rate CI; plus one-at-a-time sweeps, tornado
  data, Spearman sensitivity and bootstrap rank stability.
* **`pipeline.py` / `cli.py`** — orchestration, figures, report and a SHA-256 run
  manifest.

---

## 5. Outputs

`results/solutions/` after a default run:

* `tables/` — `therapy_catalog_resolved`, `evidence_tier_summary`,
  `endpoint_specific_effects`, `regimen_ranking`, `credit_decisions`,
  `pareto_front`, `lipid_trajectories`, `risk_curves`, `monte_carlo_summary`,
  `rank_stability`, `sensitivity_tornado`, `plaque_solution_grid`,
  `inflammation_solution_grid`, `editing_scenarios`, `regeneration_scenarios`
* cure-search tables — `cure_criteria_scores`, `lifetime_exposure_scenarios`,
  `regression_calibration`, `plaque_imaging_evidence`, `plaque_regression_grid`,
  `residual_risk`
* `json/` — `bystander_report`, `off_target_status`, `regeneration_surveillance`,
  `cure_combination_search` (+ relaxed threshold), `regression_calibration_fit`
* `figures/` — lipid trajectories, regimen ranking, uncertainty intervals,
  sensitivity tornado
* `reports/solutions_simulation_report.txt` — scope, credit rules, best regimen
  per population, top-10 table, uncertainty intervals, explicit non-claims
* `reports/cure_search_report.txt` — cure criteria, binding ceilings, lifetime
  scenarios, residual risk, and what would have to change
* `run_manifest.json` — config, status, scope string, and SHA-256 + byte size of
  every artefact

---

## 6. Results snapshot (5-year horizon, strict mode)

Events avoided per 1000 people versus no therapy; background = high-intensity
statin + aspirin.

| Population | Background | Best regimen found | Best |
| --- | --- | --- | --- |
| secondary post-MI | 75.2 | + BP intensification + colchicine (COLCOT) | 119.2 |
| post-MI after PCI | 80.0 | + alirocumab (ODYSSEY) | 98.2 |
| chronic coronary syndrome | 52.0 | + BP intensification + colchicine (LoDoCo2) | 88.4 |
| ASCVD + CKD + inflammation | 55.0 | + BP intensification | 88.1 |
| HeFH, primary prevention | 40.4 | + BP intensification | 49.8 |
| ischaemic CM with HF | 0.0 | + BP intensification | 50.0 |

How to read this honestly:

* The two largest levers in every population are the oldest ones — lipid lowering
  and blood pressure — because they are the only ones with meta-analytic
  per-unit evidence covering hundreds of thousands of patients.
* Colchicine helps where a positive trial exists and does **nothing** after modern
  PCI (CLEAR SYNERGY) or in the CKD/inflammation scenario, where no positive
  anchor applies. Ziltivekimab lowers hsCRP by 60 % and still moves events by 0.0.
* The ischaemic-CM background is 0.0 **by design**: no statin or aspirin anchor is
  transportable to that population, and DAPA-HF, STICH and RiHEART measure HF,
  mortality and functional endpoints that are never relabelled as MACE. Their
  effects are reported in `endpoint_specific_effects.csv` instead.
* Unproven modalities are the point of the three modes. For enlicitide in chronic
  coronary syndrome: strict 52.0 (identical to background), surrogate 64.3. For
  VERVE-102: strict 52.0, surrogate 65.6. The gap *is* the assumption.
* Monte Carlo (120 draws) puts the top post-MI regimen at 105.4 [91.9, 129.0]
  events avoided per 1000, with only ~76 % probability of being top-ranked — the
  ranking is not robust enough to be read as a league table of therapies.

---

## 7. Testing

647 tests, 98.5 % statement coverage of `ihd_solutions`, ~10 s wall clock.

| File | Focus |
| --- | --- |
| `test_units.py` | conversions, validators, trapezoid integration |
| `test_evidence.py` | tiers, endpoint gating, neutrality, transportability, loader errors |
| `test_catalog.py` | 30 therapies, exclusivity groups, eligibility, burden, CSV round-trip |
| `test_population.py` | 13 populations, CI → log-σ, boolean parsing, duplicate/missing-column loaders |
| `test_pharmacology.py` | onset/plateau, interval sawtooth, waning edits, adherence |
| `test_lipids.py` | multiplicative stacking, floors, component accounting |
| `test_risk.py` | zero-credit invariance, no double counting, ramp, cap, NNT, notes |
| `test_neutral_policy.py` | neutral drugs change nothing in all 3 modes; resampling cannot flip neutrality |
| `test_inflammation.py` | RK4 vs steady state, monotone dose response, CIRT negative control |
| `test_plaque.py` | thresholds, occlusion guard, counterfactual flagging |
| `test_editing.py` | HEART-2 calibration, durability, bystanders, no numeric off-target score |
| `test_regeneration.py` | beta-binomial posteriors, zero event credit, surveillance plan |
| `test_uncertainty.py` | reproducibility, interval ordering, sensitivity, rank stability |
| `test_portfolio.py` | enumeration count vs independent combinatorics, Pareto front |
| `test_pipeline_cli.py` | end-to-end run, all CLI subcommands, error exit codes |
| `test_scientific_guards.py` | all 13 populations × 30 therapies × 3 modes: invariants that must never break |
| `test_regression.py` | imaging calibration recovery, monotonicity, calcified burden never removed |
| `test_cure.py` | criteria ceilings, max-not-sum coverage, animal/active-comparator/MR refusals |
| `test_results_integrity.py` | every committed artefact matches its manifest hash; tables obey the credit rules |

Two deliberate design choices in the suite:

* `filterwarnings = ["error::RuntimeWarning"]` turns silent numpy overflow or
  divide-by-zero into a test failure.
* `test_scientific_guards.py` sweeps the whole cross-product rather than a few
  hand-picked cases, which is how the original surrogate-leak bug (unproven LDL
  agents being credited through the CTT relationship in strict mode) was caught.

Regenerate artefacts after any model change:

```bash
python run_solutions.py simulate && pytest -q
```

`test_results_integrity.py` fails if the committed outputs no longer match the
code, so results cannot silently drift.

---

## 8. Limitations, and what would change the conclusions

1. **Population-average, not individual.** Constant baseline hazard, no competing
   risks, no age/sex/time-varying hazard, no adherence dynamics beyond a static
   taking-fraction.
2. **Multiplicative independence** of relative risks is an assumption; it is
   capped at 85 % combined RRR but not otherwise validated for combinations that
   were never tested together.
3. **Pending outcomes could move the ranking.** ORION-4 (inclisiran, 2027),
   CORALreef Outcomes (enlicitide, NCT06008756), AZD0780 phase 3 (NCT07000357) and
   Lp(a)HORIZON (pelacarsen) will convert surrogate rows into either credited or
   neutral rows. The registry is built so that adding an anchor is a data edit.
4. **No safety or cost model.** Adverse events, serious infection signals (ZEUS),
   infusion reactions (7/35 in HEART-2), bleeding on aspirin, and cost are not
   optimised; `burden_score` is a crude regimen-complexity proxy, not a
   cost-effectiveness measure.
5. **Mechanistic submodels are proxies.** The 65 µm / 300 kPa plaque thresholds
   have documented exceptions; the inflammation ODE is not fitted to patient
   data; the base-editing off-target question is explicitly left to experiment.
6. **Regeneration is not remuscularisation.** The RiHEART mechanism is
   paracrine, the pivotal series is 8 patients, and no mortality or MACE claim is
   made anywhere in the code.

## 9. Primary sources

CTT Collaboration, *Lancet* 2010 · BPLTTC, *Lancet* 2021;397:1625–36 ·
IMPROVE-IT · CLEAR Outcomes · FOURIER · ODYSSEY OUTCOMES · SELECT (*NEJM* 2023) ·
CANTOS · CIRT · LoDoCo2 · COLCOT · CLEAR SYNERGY (2024) · ZEUS (2025) ·
REDUCE-IT · EMPA-REG OUTCOME · DAPA-HF · ATT Collaboration · ISCHEMIA ·
STICH 10-year extension · CORALreef Lipids (2026) · PURSUIT (AZD0780) ·
VERVE-102 HEART-2, *NEJM* 2026, doi:10.1056/NEJMoa2601283 · RiHEART conditional
approval (Japan, March 2026). Exact URLs are stored per row in
`data/solutions/trial_anchors.csv` and `data/solutions/therapy_catalog.csv`.


## 10. Quarantined invention: REVERT-1

`src/ihd_solutions/hypothesis.py` models an **invented, non-existent** construct
(`docs/invented_therapy_REVERT1.md`). It is deliberately sealed off from everything
above:

- Its IDs appear in neither `therapy_catalog.csv` nor `trial_anchors.csv`.
- It writes to `results/hypotheses/`, outside the hash-verified evidence manifest.
- Every row carries `status = invented_hypothesis_unvalidated`; the report and a
  `DISCLAIMER.txt` lead with the disclaimer.
- `tests/test_hypothesis.py` (51 tests) asserts the quarantine, and separately asserts
  that the coupling claim is genuinely falsifiable — including the arm where
  efferocytosis enhancement alone makes the lesion worse than untreated.

Run with `python run_solutions.py hypothesis`. The model is a four-state toy with
assumed rate constants; it can show that a coupling matters, never an effect size.


## 11. Simulation-only cure discovery

`src/ihd_solutions/discovery.py` searches for a cure **specification**, never a cure.
The `randomised_proof` criterion is excluded by construction because no simulation can
satisfy it; `Specification` exposes `is_candidate` and deliberately has no `is_cure`.

Method: search the four-dimensional capability space under `DEMONSTRATED` and
`CONCEIVABLE` bounds, compute the feasible volume analytically (exact for a box under
flat priors), back-calculate the minimum vector, decompose the variance of the binding
margin to rank measurements, and ablate each capability to find the indispensable set.

Outputs go to `results/discovery/`, tagged `in_silico_specification_not_evidence`, and
are covered by the same quarantine tests as the REVERT-1 hypothesis. See
`docs/cure_discovery_plan.md` for the eight-stage programme and its exit gates.


## 12. Best-of-breed composite

`src/ihd_solutions/composite.py` assembles the strongest component per disease node and
tests it across all 19 populations and all three credit modes. Assembly is greedy by
marginal event benefit, constrained by `DISEASE_NODES` (one per node) and `TARGET_OVERLAPS`
(same molecular target never stacked, even across exclusive groups). A coverage pass adds
the best component for cure-relevant but event-neutral nodes at zero event credit, so
optimising events cannot silently drop the myocardium.

The composite scores `randomised_proof = 0.0` by construction: an untested combination
inherits no proof from its components. Outputs go to `results/composite/`; see
`docs/composite_best_of_breed.md`.


## 13. Gap-closure attempt

`src/ihd_solutions/gap_solver.py` attacks the three quantified capability gaps. It adds an
optional `source_reduction` argument to `simulate_clearance` (default 0.0, so all prior
artefacts are bit-identical) modelling an ER-stress/CHOP module that lowers the apoptotic
source term. Remuscularisation is modelled as multiplicative factors on the fraction of
scar replaced, saturating at 1.0. Hazard targets are recomputed per population baseline.

Outputs go to `results/gaps/`, tagged `gap_closure_scenario_not_evidence`. See
`docs/gap_closure_attempt.md`.


## 14. Residual resolution

`src/ihd_solutions/resolution.py` solves each residual as a requirement rather than an
adjective: per-population closure by credit mode plus ranked evidence-upgrade levers
(distinguishing genuine evidence gaps from greedy `composition_reordering`); a uniform
potency multiplier searched jointly because no single clearance module suffices; and
remuscularisation scored by `Uniform(1, published effect)` draws with a chronic-injury
haircut instead of nominal point estimates.

That last change **retracts** the earlier two-factor "closable on paper" result: its
probability of reaching +20 LVEF points is 0.06%, and no factor subset is robust at 80%.
A test pins the retraction. Outputs: `results/resolution/`, tagged
`residual_resolution_scenario_not_evidence`. See `docs/residual_gap_resolution.md`.

## 15. ATLAS-1 synthesis (invented, untrialled)

`ihd_solutions.atlas` takes the best part of each construct that survived its own
tests — node-disjoint proven backbone, the REVERT-1 AND gate, gap-solver source
reduction, the resolution operating point, the discovery target — and adds one new
idea: **sequence the arms instead of stacking them**. Phase 1 (months 0-24) clears
and quietens the lesion bed; phase 2 grafts muscle into the prepared bed.

The bed-preparation link is an explicit assumption, `bed_coupling`, swept from 0 to
1. At zero the construct degenerates to simultaneous delivery and the idea
contributes exactly +0.0 pp; at full coupling it contributes +24.1 pp. Even at
perfect coupling the graft arm reaches only 50.9%, below the 80% robustness bar.

Tested across 19 populations x 3 credit modes (57 configurations). Two scorecards
are reported separately: evidence-only (best 2 of 5 criteria) and
hypothesis-conditional (also 2 of 5, after the self-audit in section 16 retracted the
nominal clearance pass). Randomised proof is 0.0 in both, myocardial restoration and
burden reversal are blocked in both, and zero configurations are declared a cure.
All six arms are load-bearing under ablation.

**ATLAS-1 is the best candidate in this study and is not a cure.** See
[`atlas_synthesis.md`](atlas_synthesis.md) for the provenance table, the
falsification sweep, the ablation, and the consolidated findings.

## 16. Self-audit: retracting the nominal clearance pass

`ihd_solutions.audit` applies the study's own distributional standard to ATLAS-1's
clearance arm, which had been scored at central module strengths while the graft arm
was scored as a distribution. At the claimed operating point (1.284x potency, source
reduction 0.50, 24 months) the nominal 80.01% becomes a median of 1.83% with
P(clearance >= 80%) = 0.00%.

Reachability: perfect modules still need source reduction >= 0.2230 to hit 80% at 24
months, so with source uniform on [0, 0.50] the probability ceiling as potency tends
to infinity is only 55.4%. `minimum_robust_potency` returns None in that case rather
than a misleading number. With source reduction *guaranteed* at 0.50, robustness
needs about 17.7x potency.

Second correction: `bed_quality` derives from clearance and `translation_floor` from
bed quality, so scoring the graft against a nominally cleared bed assumed phase 1
succeeded. Scored jointly with a paired contrast (exact 0.00 pp difference at zero
coupling), the decoupled scoring overstated the graft by 20.3 pp at full coupling;
median bed quality is 0.023.

`atlas.score_atlas` now scores both invented arms as probabilities. Consequence: the
hypothesised criterion count equals the evidence-only count everywhere, so the
invented arms add no robust criterion. Retracted nominal values are preserved as
`retracted_nominal_*` columns. See [`atlas_self_audit.md`](atlas_self_audit.md).

## 17. Evidence-layer audit

`ihd_solutions.evidence_audit` applies the distributional standard to the layer every
construct inherits. `composite.score_criteria` computed `risk_normalisation` from one
point-estimate simulation while `uncertainty.monte_carlo_regimen` already propagated
LDL, hazard-ratio, baseline-rate and adherence uncertainty; the two were never
connected. `MonteCarloResult` now records per-draw baseline and treated risk and
exposes `risk_normalisation_draws()` using the identical formula.

Only `risk_normalisation` is an effect estimate. Driver control is a structural
mechanism flag, burden and myocardial credit are capped below the threshold, and
proof is zero by construction; `criterion_ledger()` records why each is out of scope.

Across 19 populations x 3 modes: strict 1 nominal risk pass -> 0 robust (best
criteria 1 -> 0), surrogate 7 -> 2, exploratory 19 -> 18. Six passes are retracted,
including the composite's only strict criterion. One case at 78.3% is flagged
borderline (within 2 Monte Carlo SE of the bar) rather than retracted.

`lifetime_exposure` under exploratory credit exposes a second trap: its score is
clipped to [0, 1] and bimodal, with ~50% of draws at each pole, so the median flips
between 0 and 1 with sample size while P(>=0.8) is 50%. The audit reports
`mass_at_floor`, `mass_at_ceiling` and `median_uninformative` instead of a median.

See [`evidence_layer_audit.md`](evidence_layer_audit.md).

## 18. ATLAS-2: fixing the parts the audits broke

`ihd_solutions.upgrade` diagnoses rather than just retracts. Two structural causes:
`ClearanceParameters.delivered()` multiplies six of eight module strengths by one
`gate_specificity` draw (a series bottleneck with no redundancy), and delivery was
open loop. Both are architecture defects, not potency shortages.

Three upgrades: U1 redundant gate/node paths (combine by max, the conservative
choice), U2 closed-loop titration to a measured target within a hard cap, U3 verified
source control via run-in. At the solved operating point (60 months, cap 2.80x)
clearance probability goes from 0.00% to 86.12%, robust with a 2-SE margin (lower
bound 83.68%), eligible fraction 100%, so conditional and unconditional are equal.
Ablation: titration +72.83 pp, redundant nodes +42.00 pp, redundant gate +23.00 pp;
U3 is reported as not used at this horizon. Three gate routes suffice (precedented
limit); four would be beyond precedent.

Cost: 2.80x cap against a precedented 1.284x, a 2.18x precedent gap, and 60 months
instead of 24. Bisection finds the minimum robust cap at 2.590x, but that sits on the
80% bar within Monte Carlo noise, so the shipped cap buys margin instead.

Horizon correction to section 16's upstream audit: `source_required_for_target` is
0.2230 at 24 months and 0.0 at 36 months and beyond. The 55.4% "ceiling at any
potency" was horizon-specific. A gate cannot create success — unconditional
probability is invariant to a non-binding gate threshold, and a tighter gate raises
the conditional number while lowering the unconditional one.

Because the bed is now genuinely prepared (median quality 1.000, previously 0.023),
ATLAS-1's sequencing idea recovers its gain: +11.1 pp at central coupling, +24.1 pp at
full coupling, with the zero-coupling falsifier still exactly 0.0 pp. The graft
reaches 52.0% and is still not robust, so remuscularisation is now the sole binding
biology problem. Proof remains zero. See [`atlas2_upgrade.md`](atlas2_upgrade.md).

## 19. ATLAS-3: the variance fix applied to the graft

Section 18 left remuscularisation as the sole binding biology problem. ATLAS-3 attacks
it the way ATLAS-2 attacked clearance — by changing the architecture, not the assumed
potency. Module: `src/ihd_solutions/atlas3.py`; tests: `tests/test_atlas3.py` (58).

Diagnosis first. The required delivered multiplier is 5.128x and the product of the
four published maxima is 18.200x, a 3.55x surplus. Despite that surplus a single graft
into a perfect bed (translation_floor = 1.0, no chronic-injury penalty) reaches target
only 49.60% of the time. This is the nominal-product trap: E[product of uniforms] is
far below the product of their maxima. A single hypothetical U(1,8) mechanism scores
40.9%, worse than the four modest factors combined, which rules out "more potency" as
the answer.

Three levers, each ablatable: retry redundancy (best-of-k, never area multiplication,
stopping as soon as the target is met), attempt independence, and feedback gated by
measurement quality. Results at 4000 draws: single 49.60%, two attempts 82.50%, three
attempts 93.10%. Every attempt is charged a per-attempt adverse-event risk (default
5%, swept 0-20%) and success is reported both raw and risk-discounted (85.30% at three
attempts).

The precedent breach is optional. Two attempts already clear the bar and repeat dosing
has precedent (MAX_PRECEDENTED_ATTEMPTS = 2), verified across 12 seeds (mean 0.8216,
range 0.8117-0.8300). The invention is therefore the closed loop, not the procedure
count; the third attempt buys headroom for +0.95 pp of extra risk.

Falsifiers are paired. `graft_draws` returns the best gain, the attempts used, and the
first attempt's gain from the same run, so at failure correlation 1 the redundancy
gain is exactly 0.0 pp rather than zero plus noise (an earlier independent-stream
version read +0.1 pp). Maximum tolerable correlation is about 0.71, spread +/-0.02
across seeds. Feedback at zero learning contributes exactly 0.0 pp, rising to +7.05 pp
at full learning. Ablation: retry redundancy -43.50 pp, independence -26.60 pp,
feedback -4.63 pp, measurement -4.63 pp; feedback and measurement coincide by
construction because an unmeasurable signal leaves nothing to learn from.

Two model bugs were found and fixed. Raising each factor's lower bound as "feedback"
made the product of the minima exceed the requirement and forced P = 100%; it was
replaced by one bounded aggregate uplift (FEEDBACK_MAX_UPLIFT = 0.5) gated by
measurement quality, because a clinical signal reports total engraftment rather than
per-mechanism attribution. A test asserts no configuration reaches certainty. Separately,
the first specification search reported three attempts because the 2-SE margin at small
N rejected a two-attempt design that is genuinely above the bar; both specification
searches now floor their draw count at MIN_SPECIFICATION_DRAWS = 4000. A specification
must not move with the sampling budget.

Both arms together: clearance 83.50%, graft 92.45%, optimistic joint 77.20%, proof
0.00%. The joint is a product and so assumes arm independence, which the measured bed
coupling contradicts; it is labelled optimistic rather than presented as a combined
result. The discovery VOI predates ATLAS-2 and is stale; recomputed, substrate removal
still ranks first (+15.25 pp) over the graft (+6.30 pp), but now because its
architecture delivers a lower probability, not because it lacks one.

Proof remains zero and no cure is declared. See [`atlas3_graft.md`](atlas3_graft.md).

## 20. ATLAS-4: precedent search, three corrections, exhaustive design grid

Modules: `src/ihd_solutions/precedent.py` (12 external sources) and
`src/ihd_solutions/atlas4.py`; tests: `tests/test_atlas4.py` (59).

A literature search asked whether the step ATLAS-3 called missing had already been done.
Most of it had. CATO (NCT06145035) is a Phase IIA randomised double-blind placebo-controlled
trial giving 60 patients with ischaemic cardiomyopathy four intravenous UC-MSC doses two
months apart against one dose and placebo, primary endpoint change in LVEF by cardiac MRI
at 12 months (47 enrolled as of 2025-11-30). HOPE-2 gave four doses in humans in 2022.
Reporter-gene PET predicts later function in animals with rho=0.76 against echo and
rho=0.80 against MRI.

Three corrections follow. The precedented dose limit is four, not two, so ATLAS-3's three
attempts were never beyond precedent. Retries combine cumulatively, not as best-of-k:
three 12M CPC doses gave +13.0 pp against +4.7 pp for one 36M dose, and the extra
injections could not be replaced by tripling the single dose. And the cells do not
regenerate myocardium at all - the effect is paracrine, anti-inflammatory and
anti-fibrotic, stated explicitly in the porcine study.

The registry enforces evidence discipline: animal work never earns human efficacy credit,
an ongoing trial supports a dose count but not efficacy, and the null blinded human trial
earns nothing. Split retention is anchored to the rat arithmetic, (13.0/3)/4.7 = 0.922,
and the increment shape is normalised so no rodent magnitude is imported as a human claim.

Cumulative beats best-of-k at every dose count on shared draws (+8.85 pp at four doses).
Repeat responsiveness is tempered to the human observation that 10 of 30 patients
responded to a second dose, and gain is gated by remaining headroom, the ceiling effect
measured in that same trial. Result: 82.97% at four doses, robust with margin, within
precedent, 7.98% cumulative adverse risk. A 0.7 translation floor drops it to 76.40%.

A consequence worth recording: ATLAS-3's independence specification was operator-specific.
Best-of-k retains 35.8% of its gain at full dose correlation; cumulative retains 90.6% and
stays robust. The measurement falsifier is exact: iron-oxide MRI and acute-retention
imaging gain 0.00 pp because iron signal survives cell death and retention does not
predict engraftment; only reporter-gene PET is viability-specific, and under the corrected
operator it adds just +0.40 pp.

The full Cartesian grid of 360 configurations was scored with nothing pre-filtered: 176
robust, 160 robust with margin, 288 within precedent, 48 reading exactly 100%, 8 admissible.
Every point at certainty is a saturation artefact. The first admissibility filter was too
weak and returned a rodent-optimistic point at 99.83% as "best"; it now also requires the
repeat-response assumption to be no more optimistic than the human data, which cut
admissible points from 32 to 8.

On the question of driving every criterion to 100%: measured catalogue ceilings are 1.000
for durable driver control and randomised proof, 0.600 for burden reversal, 0.500 for
myocardial restoration and 0.061 for risk normalisation. The last three are capped by
explicit evidence rules rather than by an unexplored design space, and proof attaches to
the exact regimen randomised so it cannot transfer to an untested combination. Raising
them would mean deleting the caps. The one design that does read 100% is the
rodent-optimistic variant, and it stays at 99.98% under a harsh translation penalty,
showing the certainty comes from the responder assumption rather than the biology.

Genuinely still missing: response-guided retry timing, proposed in print in 2019 and never
tested, and a viability-specific engraftment signal in humans. No cure is declared. See
[`atlas4_precedent_search.md`](atlas4_precedent_search.md).
