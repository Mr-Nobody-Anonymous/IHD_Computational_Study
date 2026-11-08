# Computational Biophysics & Systems Biology of Ischemic Heart Disease (IHD)

> **Full report:** [`docs/FULL_REPORT.md`](docs/FULL_REPORT.md) consolidates every
> agent (47), trial anchor (48), population (19), compound library, molecular
> experiment and simulation study in one document, with an epistemic status label
> on each section and a complete retraction ledger. A machine-readable index of all
> 95 chemical entities is in [`docs/FULL_REPORT_chemical_index.csv`](docs/FULL_REPORT_chemical_index.csv).

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![PDB Structures](https://img.shields.io/badge/PDB-3BPS%20%7C%203GCX%20%7C%203P5C%20%7C%202P4E-brightgreen.svg)](https://www.rcsb.org/)

A reproducible computational biology framework investigating the structural biophysics, vascular biomechanics, inflammatory signaling networks, and candidate therapeutic modalities for **Ischemic Heart Disease (IHD)**.

---

## 🧪 Solution simulation & test suite (`ihd_solutions`)

A tested, evidence-graded scenario simulator covering **every major IHD solution
pathway** — 47 therapies, 19 populations, 48 published trial anchors — lives in
`src/ihd_solutions/`, with data in `data/solutions/` and 783 tests in `tests/`.

```bash
pip install -r requirements-dev.txt
python run_solutions.py validate      # registry checks
python run_solutions.py simulate      # 6 populations x 337 regimens x 3 modes -> results/solutions/
pytest -q                             # 783 tests, ~99% coverage of ihd_solutions
```

Pathways covered: LDL/ApoB lowering (statins, ezetimibe, bempedoic acid, PCSK9
antibodies, inclisiran, oral enlicitide and AZD0780), *PCSK9* base editing
(VERVE-102), Lp(a) lowering, residual inflammation (colchicine, IL-1β/IL-6/NLRP3,
methotrexate), metabolic therapy (GLP-1, SGLT2), blood pressure, antithrombotic
therapy, lifestyle, revascularisation and regenerative cell sheets.

**How the evidence is handled.** Every effect size is data, with a source URL and
an evidence tier. Three credit modes are reported side by side: `strict` (only
randomised hard-outcome evidence in a matching or explicitly transportable
population), `surrogate` (additionally transports the CTT LDL / BPLTTC pressure
relationships to biomarker-only agents) and `exploratory`. hsCRP and Lp(a)
improvements **never** create event credit; heart-failure, functional and
mortality endpoints are never relabelled as MACE; neutral trials (ZEUS, CIRT,
CLEAR SYNERGY, ISCHEMIA) are held at HR = 1.0. Every credit decision is written to
`results/solutions/tables/credit_decisions.csv` with its reason.

### Is any of it a cure?

`cure.py` states five explicit cure criteria (durable driver control, burden
reversal, myocardial restoration, risk normalisation, randomised proof), scores
all 47 interventions, and searches for the smallest combination that satisfies
all five. **None does.** The binding ceilings are burden reversal (0.60),
myocardial restoration (0.50) and risk normalisation (0.74 even for the best
stack). The strongest cure-adjacent lever is starting durable LDL-C reduction
decades earlier, quantified with a Mendelian-randomisation scenario that is
labelled as such and never used as a treatment effect.

📖 **Methods and test strategy: [`docs/solutions_simulation.md`](docs/solutions_simulation.md)**
🔬 **The 2026 cure search, with sources: [`docs/cure_search_2026.md`](docs/cure_search_2026.md)**

### 🧪 An invented therapy, quarantined: REVERT-1

[`docs/invented_therapy_REVERT1.md`](docs/invented_therapy_REVERT1.md) designs a
construct that does **not exist** — a lesion-AND-gated particle coupling cholesterol
crystal dissolution, export, corpse clearance and influx blockade, so the plaque's own
cholesterol fuels its removal. It is modelled, ablated and stress-tested, **not**
credited as evidence:

```bash
python run_solutions.py hypothesis     # -> results/hypotheses/ (10 artefacts)
```

The interesting results are negative. Efferocytosis enhancement *alone* ends up worse
than untreated (the phagocyte eats a cholesterol-laden corpse, cannot export it, dies,
and becomes necrotic core) while looking like a success at 12 months. The delivery gate
outranks every payload in sensitivity. And even with every module set to perfect the
architecture plateaus at 83.5% substrate removal — it never eradicates, and it does not
regenerate myocardium. No hypothesis ID can enter the catalogue, the anchors or the
cure search; tests enforce that.

> **REVERT-1 is an unvalidated hypothesis, not a therapy or a recommendation.**

### 🎯 What a cure would require, searched in simulation: the discovery engine

[`docs/cure_discovery_plan.md`](docs/cure_discovery_plan.md) answers *which candidate is
best*, *what is missing*, and gives an eight-stage simulation-only programme.

```bash
python run_solutions.py discover       # -> results/discovery/ (10 artefacts)
```

It searches the capability space twice — once bounded by what has actually been
demonstrated, once by generous bounds. The feasible region is **exactly zero** under
demonstrated capabilities and **0.037%** under generous ones; the difference between
those boxes *is* the missing discovery. Three capabilities are indispensable:
substrate clearance (gap ×1.86), stack risk normalisation (×1.21) and remuscularisation
(**×5.13**, +20 LVEF points needed versus +3.9 demonstrated). Driver durability is
already solved by one-time editing and is provably not part of the problem.

### 🧩 Best-of-breed composite: strongest part per node, tested everywhere

[`docs/composite_best_of_breed.md`](docs/composite_best_of_breed.md) assembles the best
component per disease node and tests it across **19 populations × 3 credit modes**.

```bash
python run_solutions.py composite      # -> results/composite/ (8 artefacts)
```

Greedy selection by marginal evidence-graded event benefit, under three anti-double-counting
rules: coverage combines by **maximum**, one component per node, and the **same molecular
target is never stacked** (PCSK9 editing + a PCSK9 antibody would otherwise slip through).
A fourth rule: an **untrialled combination earns zero randomised proof**.

Result: up to **134 events avoided per 1000**, and risk normalisation finally crosses the
bar (0.819 strict / 0.974 surrogate). Still **0 cures** — burden reversal (0.60), myocardial
restoration (0.50) and proof (0.00) block everywhere. Invented components provably change
no score.

### 🔧 Attempting to close the three capability gaps

[`docs/gap_closure_attempt.md`](docs/gap_closure_attempt.md) attacks the quantified gaps
structurally rather than by tuning.

```bash
python run_solutions.py gaps           # -> results/gaps/ (9 artefacts)
```

- **Hazard ratio ×1.21 — already closed.** The composite achieves HR 0.282 vs 0.298 required.
  But the target was baseline-specific; recomputed per population, 1 of 19 closes and the
  rest need a further 16–32%.
- **Clearance ×1.86 — the plateau was structural.** A closed-form steady state matched the
  simulation exactly, revealing a constant apoptosis source no module touched. Adding an
  ER-stress/CHOP source-reduction module (bounded by the ~35% apoptosis / ~50% necrosis
  reduction from CHOP deletion) lifts the ceiling to 89.5%. Caveat: closing in 2 years needs
  near-perfect delivery; at central strengths it needs a much longer horizon.
- **Remuscularisation ×5.13 — decomposes into two co-delivery changes.** Endothelial +
  epicardial/stromal co-delivery gives ×5.2 on scar replacement → +20.3 LVEF points. Needs
  neither dose escalation nor gene editing. Every multiplier is animal-level, and the
  margin is nil.

### 🧮 Resolving the residuals step by step

[`docs/residual_gap_resolution.md`](docs/residual_gap_resolution.md) converts each residual
into a costed requirement — and **retracts** one earlier conclusion.

```bash
python run_solutions.py resolve        # -> results/resolution/ (9 artefacts)
```

[`docs/atlas_synthesis.md`](docs/atlas_synthesis.md) assembles the best *validated*
part of every preceding construct into one new candidate, **ATLAS-1**, and adds a
single new idea — clear the lesion bed first, then graft into it — whose coupling
is swept to zero so it can be falsified. Tested across all 19 populations x 3
credit modes x 5 criteria. Best candidate on construction quality; **not a cure**,
and after the self-audit below it reaches **2 of 5** criteria, not the 3 originally
claimed.

```bash
python run_solutions.py atlas          # -> results/atlas/ (9 artefacts)
```

[`docs/atlas_self_audit.md`](docs/atlas_self_audit.md) audits that candidate against
the study's own standard and **retracts its headline**. The clearance operating point
was a nominal pass: under delivery uncertainty P(clearance ≥ 80%) is 0.00%, and no
potency makes it robust while the source term is uncertain (ceiling 55.4%). Because
bed quality derives from clearance, the sequencing gain was overstated by ~20 pp —
the arms of a sequenced construct cannot be scored independently. Both invented arms
are now scored as probabilities, and the retracted values are kept as labelled
columns.

```bash
python run_solutions.py audit          # -> results/audit/ (7 artefacts)
```

[`docs/evidence_layer_audit.md`](docs/evidence_layer_audit.md) asks the same question
of the layer every construct inherits, and finds the same defect: `risk_normalisation`
was a point estimate even though the repository already had CI-propagating Monte
Carlo. Connected and re-run across all 19 populations x 3 modes, **6 passes are
retracted** — including the composite's only strict-mode criterion, so the best strict
result falls to **0 of 5**. Exploratory credit largely survives (18/19). One case at
78.3% is flagged borderline rather than retracted.

```bash
python run_solutions.py audit-evidence # -> results/evidence_audit/ (6 artefacts)
```

[`docs/atlas2_upgrade.md`](docs/atlas2_upgrade.md) then **fixes** what the audits
broke, by changing the architecture rather than inflating assumed potency: redundant
targeting paths remove a series bottleneck (one gate draw throttled six of eight
modules), and closed-loop titration replaces open-loop dosing. Clearance goes from
**0.00% to 86.12%** probability of target, robust with margin, across **100%** of the
population. Titration alone is worth +72.8 pp. The stated price is a **2.18x**
beyond-precedent potency cap and a 60-month horizon. Because the bed is now genuinely
prepared, ATLAS-1's sequencing idea recovers its full gain (+24.1 pp). Remuscularisation
(52.0%) and randomised proof (0) remain unfixed.

```bash
python run_solutions.py upgrade        # -> results/upgrade/ (10 artefacts)
```

[`docs/atlas3_graft.md`](docs/atlas3_graft.md) applies that same variance fix to the
last open biology gap. The graft never needed more potency: the four published factors
already exceed the requirement **3.55-fold**, yet one graft into a *perfect* bed hits
target only **49.6%** of the time. Adding staged retries plus a measured engraftment
signal lifts it to **93.1%** (3 attempts) or **82.5%** (2 attempts, **within precedent**),
with every extra procedure charged an adverse-event risk (**8.38%** cumulative at three).
Redundancy is worth **+38.9 pp** when retries fail independently and **exactly 0.00 pp**
when they share a failure mode; the design tolerates correlation up to **~0.71**. So the
missing step is a modality that can be **repeated and measured** — not a stronger one.
The recomputed value of information still ranks clearance first (**+15.25 pp** vs
**+6.30 pp**). Randomised proof remains **0**; no cure is declared.

```bash
python run_solutions.py atlas3         # -> results/atlas3/ (10 artefacts)
```

[`docs/atlas4_precedent_search.md`](docs/atlas4_precedent_search.md) then asks whether
anyone has **already done** the missing step. Most of them have: **CATO
(NCT06145035)** is randomising patients with ischaemic cardiomyopathy to **four**
intravenous cell doses right now, and reporter-gene PET already predicts later cardiac
function in animals (**rho=0.76-0.80**). That search forced **three corrections** to
ATLAS-3: the precedented dose limit is 4, not 2; retries combine **cumulatively**, not
best-of-k; and these cells are **paracrine, not regenerative**. Corrected and tempered
by the human responder rate (10 of 30 responded to a second dose), four doses reach
**82.97%**, within precedent. All **360** design combinations were tested; **48 read
exactly 100%** and every one is a saturation artefact of assuming patients respond like
rats. Two of the five cure criteria already reach 1.000; the rest are capped by evidence
rules, so **100% across all five is unreachable by searching**. What is genuinely still
missing is narrower than claimed: response-guided retry timing, proposed in 2019 and
still untested.

```bash
python run_solutions.py atlas4         # -> results/atlas4/ (11 artefacts)
```

[`docs/gap_solutions.md`](docs/gap_solutions.md) searches each of the six remaining open
gaps and designs a solution for it. **Three of the six prior claims turned out to be wrong
or overstated:** human remuscularisation is already proven histologically (Jebran,
*Nature* 2025; BioVAT-HF, *NEJM* 2026 - so the gap is **magnitude**, needing ~2.9x more
effect than demonstrated, and maturation supplies part of it free); the hNIS reporter has
already been imaged in **patients**, just never in a cardiac graft; and PCSK9-editing
durability is **measured** to 18 months, not assumed. The most actionable new result is
that **CATO has only ~25% power** for the plausible 3 pp increment at 20 per arm, so a null
readout would be uninformative. The clearance potency gap falls from 2.18x to **1.87x** by
substituting redundancy - but only after excluding a cheaper point whose redundancy was
itself beyond precedent. And the PCSK9 outcomes trial needs on the order of **125 events**,
not a mega-trial. Three gaps need no new biology at all.

```bash
python run_solutions.py gap-solutions  # -> results/gap_solutions/ (10 artefacts)
```

[`docs/mechanism_recombination.md`](docs/mechanism_recombination.md) takes the **successful
part of each architecture** and composes them. Eight mechanisms were extracted, each
measured alone and in every pair, and only the parts that earned their cost were kept.
**Maturation and retry redundancy overlap by 14.08 pp** - crediting both additively would
count the same biology twice - while **retry + cumulative dosing is genuinely synergistic
(+10.48 pp)**. Two parts turned out **inert**: the feedback loop and the reporter cassette
it exists to enable, which removes the design's only genetic modification.
**Response-guided selection is a trial-design lever, not a therapy**: the screen
concentrates responders but excludes patients, so over the whole population it is
**-3.52 pp**. Consequently **switching everything on is 3.52 pp worse than selecting**.
The result, **ATLAS-5**, keeps three mechanisms, reaches **97.92%** (risk-adjusted
**92.53%**) against ATLAS-4's 89.60%, and uses **1.10 procedures on average** because
maturation makes the first dose usually sufficient. Honest limits: only **18 of 36**
assumption scenarios are robust, and the model **cannot resolve maturation above 1.25x**
at all. The ATLAS-3 precedent limit was also corrected in code (**2 -> 4**, from the
registry), so ATLAS-3 was never beyond precedent. Randomised proof: **zero**.

```bash
python run_solutions.py recombine      # -> results/recombine/ (13 artefacts)
```

- **Hazard ratio — solved as procurement.** All 19 populations close by exploratory credit.
  One PCSK9-editing outcomes trial unlocks **17 of 19**. No new biology needed.
- **Clearance — solved as an engineering spec.** Correction: no *single* module suffices
  (best alone 63.3%, gate only 62.3%). Cheapest operating point is **1.28× central potency
  across all modules + 0.50 source reduction at 24 months**, or central potency at 60 months.
- **Remuscularisation — RETRACTED.** The two-factor "×5.2 vs ×5.13" pass was arithmetic
  luck: under Uniform(1, published effect) it reaches +20 only **0.06%** of the time. Even
  all four factors give 50.9%, falling to 26.0% with a chronic-injury haircut. **No
  combination is robust at 80%.**

> **A cure cannot be discovered in simulation.** One of the five criteria requires
> randomised human outcome evidence. The engine returns a *specification* and is
> structurally prevented from calling it a cure — `Specification.is_candidate` exists,
> `is_cure` does not, and a test enforces that.

> **Scope note.** This package produces *scenario arithmetic*, not clinical
> predictions, treatment recommendations, or evidence of a cure. The sections
> below describe the **original** structural/biophysical study
> (`src/01_*.py` … `src/05_*.py`). Those scripts and their earlier "curative"
> framing are left unchanged and are **not** validated by the new test suite;
> where they conflict with the graded evidence above, prefer the graded evidence.

---

## 📌 Scientific Context & The "Flat PPI" Paradox

Ischemic Heart Disease (IHD) remains the leading cause of human mortality worldwide (~9 million deaths annually). Although standard-of-care pharmacotherapies (statins) and revascularization procedures (PCI stents, CABG surgery) prolong event-free survival, they manage symptoms and slow progression without permanently halting atherogenesis or regrowing infarcted myocardium.

This repository provides an empirical, multi-scale computational investigation addressing the core biological bottlenecks:
1. **The Biophysical Paradox of PCSK9:** Why 15+ years of classical medicinal chemistry failed to yield an approved oral small molecule targeting the PCSK9–LDLR interface (quantified via multi-crystal solvent-accessible surface area burial).
2. **Experimental Alanine-Scanning Validation:** Benchmark of computational residue burial against empirical wet-lab Surface Plasmon Resonance (SPR) binding affinities from Kwon et al. (*PNAS* 2008).
3. **Plaque Rupture Biomechanics:** A 60-month finite-difference simulation of fibrous cap remodeling, necrotic core dynamics, and Laplace tensile hoop stress.
4. **Vascular Systems Immunology:** Ordinary differential equation (ODE) modeling of the $\text{oxLDL} \rightarrow \text{NLRP3} \rightarrow \text{IL-1}\beta \rightarrow \text{hs-CRP}$ inflammatory cascade (CANTOS/COLCOT axis).
5. **CRISPR Base-Editing Specificity:** Computational genomic off-target scoring (CFD) of the clinical-stage VERVE-102 ABE8e guide RNA.

---

## 📂 Repository Structure

```
IHD_Computational_Study/
├── data/
│   ├── structures/                  # High-resolution crystallographic coordinate files (.cif)
│   │   ├── 3BPS.cif                 # PCSK9 + LDLR EGF-A (2.30 Å resolution)
│   │   ├── 3GCX.cif                 # PCSK9 + LDLR EGF-A (1.98 Å resolution, Kwon et al. 2008)
│   │   ├── 3P5C.cif                 # PCSK9 + full LDLR extracellular domain (3.30 Å)
│   │   ├── 2P4E.cif                 # High-resolution PCSK9 apo / domain map (1.97 Å)
│   │   ├── 4NMX.cif                 # PCSK9 with engineered inhibitor peptide (1.90 Å)
│   │   └── 4HHB.cif                 # Human hemoglobin tetramer reference (1.74 Å)
│   └── experimental/
│       └── kwon_2008_alanine_scanning.csv  # Peer-reviewed SPR binding data (PNAS 2008)
├── src/
│   ├── 01_multistructure_sasa_benchmark.py # Multi-structure SASA & SPR correlation
│   ├── 02_interface_contact_mapping.py     # Atomic contacts, H-bonds, and docking grid
│   ├── 03_plaque_biomechanics_model.py     # 60-month coronary atherogenesis simulation
│   ├── 04_immune_cytokine_network.py       # 180-day vascular inflammatory ODE network
│   └── 05_crispr_genomic_safety_cfd.py     # Genome-wide CFD off-target screening
├── results/
│   ├── figures/
│   │   ├── experimental_validation_plot.png # Cross-crystal burial & SPR correlation
│   │   ├── plaque_simulation_curves.png     # 4-panel plaque biomechanics trajectories
│   │   └── cytokine_dynamics.png            # 4-panel inflammatory biomarker timecourses
│   ├── tables/
│   │   ├── multi_structure_sasa_variability.csv
│   │   ├── experimental_vs_computational_correlation.csv
│   │   ├── plaque_biomechanics_timeline.csv
│   │   ├── cytokine_network_data.csv
│   │   └── verve102_offtarget_scores.json
│   └── reports/
│       ├── interface_contact_summary.txt
│       ├── rigorous_sasa_report.txt
│       └── verve102_genomic_safety_report.txt
├── docs/
│   ├── solutions_simulation.md                   # Methods, results & test strategy for ihd_solutions
│   ├── cure_search_2026.md                       # 2026 cure search: candidates, verdicts, sources
│   ├── invented_therapy_REVERT1.md               # INVENTED hypothesis, quarantined from evidence
│   ├── cure_discovery_plan.md                    # simulation-only cure discovery programme
│   ├── composite_best_of_breed.md                # best part per node, tested across all pops
│   ├── gap_closure_attempt.md                    # structural attack on the three capability gaps
│   ├── residual_gap_resolution.md                # step-by-step solve; retracts the remuscularisation pass
│   ├── atlas_synthesis.md                        # ATLAS-1: best part of each construct, sequenced; best candidate
│   ├── atlas_self_audit.md                       # self-audit; retracts the nominal clearance pass
│   ├── evidence_layer_audit.md                   # audits the inherited evidence scores; strict 1/5 -> 0/5
│   ├── atlas2_upgrade.md                         # ATLAS-2: fixes clearance (0% -> 86%), states the price
│   ├── atlas3_graft.md                           # ATLAS-3: fixes the graft by retry + feedback, not potency
│   ├── atlas4_precedent_search.md                # ATLAS-4: literature search, 3 corrections, all 360 designs tested
│   ├── gap_solutions.md                          # designs + power calculations for the six open gaps
│   ├── mechanism_recombination.md                # ATLAS-5: mechanism ablation, overlaps, all 64 compositions
│   ├── FULL_REPORT.md                            # consolidated report: all agents, compounds, experiments, results
│   ├── FULL_REPORT_chemical_index.csv            # 95 entities, source-tagged, QA-flagged
│   ├── scientific_methodology_and_limitations.md # Critical reality check & assumptions
│   ├── IHD_master_solution_blueprint.md         # Comprehensive multi-scale dossier
│   └── IHD_cure_deep_research.md                # Research review on VERVE-102 & iPSCs
├── data/solutions/                  # Evidence-graded inputs for the solution simulator
│   ├── trial_anchors.csv            # 48 published effect estimates + tiers + source URLs
│   ├── therapy_catalog.csv          # 47 therapies: pharmacology, biomarker effects, burden
│   └── populations.csv              # 19 scenario populations with baseline rates and CIs
├── src/ihd_solutions/               # Tested simulation package (20 modules)
├── tests/                           # 783 tests incl. cross-product scientific guards
├── results/solutions/               # Generated tables, figures, report + SHA-256 manifest
├── run_solutions.py                 # CLI entry point (validate / simulate / evaluate / ...)
├── pyproject.toml                   # Packaging, pytest, coverage and ruff configuration
├── requirements-dev.txt             # Simulator + test dependencies
├── .github/workflows/tests.yml      # Lint, tests and artefact-hash verification on 3.10-3.12
├── .gitignore
├── requirements.txt
└── README.md
```

---

## 🔬 Key Scientific Modules & Results

### 1. Cross-Crystallographic Interface Burial & Alanine-Scanning Benchmark
Using the **Shrake-Rupley algorithm** ($r_{\text{probe}} = 1.4\text{ \AA}$ water probe, 960 points/atom), we evaluated per-residue solvent-accessible surface area burial across three independent crystal structures (`3BPS`, `3GCX`, `3P5C`):

```
ResID   ResName   3GCX (1.98Å)   3BPS (2.30Å)   3P5C (3.30Å)    Mean ± SD (Å²)   Experimental ΔΔG (kcal/mol)
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
379     PHE           81.9           77.4           92.3          83.9 ± 6.2     +3.00 (Binding abolished)
369     ILE           65.2           64.9           50.7          60.3 ± 6.8     +1.87 (23-fold loss)
153     SER           57.9           66.2           26.5          50.2 ± 17.1    +0.10 (Neutral effect)
377     THR           41.1           43.0           54.9          46.3 ± 6.1     +1.31 (9-fold loss)
238     ASP           36.7           43.2           31.9          37.3 ± 4.6     +1.13 (6.7-fold loss)
378     CYS           36.2           37.0           29.2          34.1 ± 3.5     N/A (Disulfide bridge)
374     ASP           36.3           38.8           26.9          34.0 ± 5.1     +0.82 (4-fold loss / D374Y: -1.91)
380     VAL           28.3           37.5           27.0          30.9 ± 4.6     +0.71 (3.3-fold loss)
194     ARG           33.1           20.6           26.0          26.6 ± 5.1     +1.75 (19-fold loss)
```

#### Why Classical Small Molecules Failed:
* Total buried surface area (BSA) is **$591.5\text{ \AA}^2$**, with the altered PCSK9 footprint covering **$2,676\text{ \AA}^2$**.
* A Rule-of-5 small molecule ($<500\text{ Da}$) covers only **$\sim 250 - 350\text{ \AA}^2$** ($<40\%$ of the interface).
* The only viable molecular modalities are:
  1. **Synthetic Macrocycles (e.g., Merck's MK-0616):** Large cyclic peptides ($\sim 1.5\text{ kDa}$) designed to cover $>700\text{ \AA}^2$.
  2. **Gene Silencing (VERVE-102 base editing / Inclisiran siRNA):** Bypasses the flat interface by arresting protein synthesis upstream in hepatocytes.

---

### 2. Plaque Biomechanics & Vulnerability Dynamics
Simulates a 60-month continuum finite-difference model of human Left Anterior Descending (LAD) coronary artery atherogenesis under physiological pulsatile pressure ($140\text{ mmHg} = 18.66\text{ kPa}$):
* **Untreated High-Risk (LDL 190 mg/dL):** Fibrous cap thins below the $65\text{ \mu m}$ critical threshold, and tensile hoop stress surpasses the $300\text{ kPa}$ material failure threshold, predicting rupture at **Month 23.0**.
* **High-Intensity Statin (LDL 95 mg/dL):** Prolongs event-free survival, delaying rupture to **Month 56.0**; persistent residual vascular inflammation ultimately causes late cap failure.
* **VERVE-102 Base Editing (LDL 38 mg/dL):** Drives net cholesterol efflux, reversing necrotic lipid core volume from $0.40\text{ mm}^2$ to **$0.26\text{ mm}^2$** and maintaining thick, stable caps ($>156\text{ \mu m}$).
* **Dual Curative Strategy:** Combines gene silencing with inflammasome blockade, achieving complete plaque stabilization ($180\text{ \mu m}$ cap).

---

### 3. Systems Biology of Vascular Inflammation (CANTOS/COLCOT)
An ODE network tracking the inflammatory cascade over 180 days:
* While statins moderately lower biomarkers, **residual inflammatory risk** remains high ($\text{hs-CRP} \approx 37.7\text{ mg/L}$ in inflamed plaques).
* Targeted inflammasome blockade (Colchicine / Canakinumab) combined with radical LDL reduction drives hs-CRP down to **$1.78\text{ mg/L}$** (safely below the clinical $<2.0\text{ mg/L}$ threshold) and suppresses destructive matrix metalloproteinase (MMP) activity to **$20.8\%$**.

---

### 4. CRISPR Base-Editing Genomic Specificity (VERVE-102)
Off-target cutting frequency determination (CFD) analysis for the *PCSK9* Intron 1 splice donor guide (`GGGCAGCGCCACCCACAGGA` PAM `TGG`):
* Mismatch penalties across the 10 bp seed region ensure negligible off-target cleavage in homologous proprotein convertases (**PCSK1**, **PCSK2**, **FURIN**, **MBTPS1**; all CFD scores $<0.086$).
* The use of **Adenine Base Editing (ABE8e-nCas9)** avoids double-strand breaks (DSBs), thereby eliminating chromosomal translocations, chromothripsis, and p53 damage responses associated with classical Cas9 nucleases.

---

## 🚀 Quickstart & Reproducibility

### Installation

Clone the repository and install dependencies:

```bash
git clone https://github.com/Mr-Nobody-Anonymous/IHD_Computational_Study.git
cd IHD_Computational_Study
pip install -r requirements.txt
```

### Running the Analysis Suite

Execute the scripts in sequence:

```bash
# 1. Multi-structure SASA calculation & experimental SPR correlation
python src/01_multistructure_sasa_benchmark.py

# 2. Interface contact mapping and docking grid derivation
python src/02_interface_contact_mapping.py

# 3. 60-month coronary plaque biomechanics simulation
python src/03_plaque_biomechanics_model.py

# 4. 180-day vascular inflammatory cytokine ODE network
python src/04_immune_cytokine_network.py

# 5. CRISPR base-editing genomic specificity & CFD scoring
python src/05_crispr_genomic_safety_cfd.py
```

All generated plots are written to `results/figures/`, tables to `results/tables/`, and reports to `results/reports/`.

---

## ⚠️ Limitations & Scientific Transparency

Please read [`docs/scientific_methodology_and_limitations.md`](docs/scientific_methodology_and_limitations.md) for full context:
* **Structural calculations** are based directly on empirical crystallographic ground truth from the RCSB Protein Data Bank.
* **Biomechanics and cytokine models** are idealized continuum and dynamical systems intended to explore mechanistic concepts; they are not patient-specific clinical diagnostic tools.
* **Genomic off-target scores** are in-silico screening estimates; regulatory dossiers require wet-lab assays (GUIDE-seq, CIRCLE-seq) in primary hepatocytes.

---

## 📚 References

1. **Kwon HJ, et al.** (2008). Molecular basis for LDL receptor recognition by PCSK9. *Proc Natl Acad Sci USA*, 105(6):1820–1825.
2. **Lo Surdo P, et al.** (2011). Mechanistic implications for LDL receptor degradation from the PCSK9/LDLR structure at neutral pH. *EMBO Rep*, 12(12):1300–1305.
3. **Ridker PM, et al. (CANTOS)** (2017). Antiinflammatory Therapy with Canakinumab for Atherosclerotic Disease. *N Engl J Med*, 377(12):1119–1131.
4. **Tardif JC, et al. (COLCOT)** (2019). Efficacy and Safety of Low-Dose Colchicine after Myocardial Infarction. *N Engl J Med*, 381(26):2497–2505.
5. **Verve Therapeutics / NEJM** (2025). Interim Safety and Efficacy of VERVE-102 in Heterozygous Familial Hypercholesterolemia. *N Engl J Med*.

---

## 📄 License
This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
