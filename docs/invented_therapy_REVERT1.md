# REVERT-1 — an invented therapy architecture for ischaemic heart disease

**Status: invented hypothesis. Not a therapy. Not a trial result. Not efficacy evidence.**

Everything below is either (a) a published preclinical or mechanistic precedent, cited,
or (b) an assumption I made, labelled as such. The construct itself does not exist. It
has never been built, dosed in an animal, or given to a human. Nothing here is a
treatment recommendation.

REVERT-1 = **RE**gression **V**ia **E**fferocytosis-coupled **R**everse **T**ransport.

Reproduce: `python run_solutions.py hypothesis` → `results/hypotheses/`
Code: `src/ihd_solutions/hypothesis.py` · Tests: `tests/test_hypothesis.py` (51 tests)

---

## 1. Bottom line

I designed the strongest construct I could justify, modelled it, and tried to kill it.
Three results matter more than the design itself:

1. **The rate-limiting step is not lipid lowering, and not corpse clearance — it is the
   coupling between them.** Boosting efferocytosis alone made lesions *worse than
   untreated* in the model, because a macrophage that eats a cholesterol-laden corpse
   and cannot export that cholesterol dies and becomes necrotic core itself. This is
   not my speculation; it is published mouse genetics ([Yancey 2010](https://www.ahajournals.org/doi/10.1161/circresaha.110.217281)).
2. **The delivery gate, not the payloads, is the highest-leverage parameter**
   (Spearman ρ = +0.56, ahead of every payload). The same architecture with poor
   lesional targeting clears nothing.
3. **Even this design does not meet the cure bar.** Under central assumptions it
   plateaus at **66.8%** removal of rupture-prone substrate. With *every module set to
   perfect* — physically unachievable — it reaches **83.5%**, crossing the 80% criterion
   only at ~34.5 months, and still leaves ~16.5% behind. It does not regenerate a
   single cardiomyocyte.

So the honest answer to "invent one that attaches to all of it and fixes it" is: here
is a construct that attacks more nodes simultaneously than anything in the catalogue,
here is the quantified reason it still would not be a cure, and here is the cheapest
experiment that would tell us whether the central claim is even true.

---

## 2. What "all of it" actually is

IHD is not one lesion. From the node map this repository already encodes:

| Node | Current best intervention | What it does *not* do |
|---|---|---|
| Atherogenic lipoprotein supply | One-time PCSK9/ANGPTL3 editing, ~62% LDL-C reduction | Does not remove existing plaque |
| Lipid retention in the wall | Aggressive LDL lowering | Threshold for regression ≈74 mg/dL; slow |
| Defective corpse clearance | Nothing approved | Necrotic core keeps growing |
| Macrophage cholesterol overload | Nothing approved | Foam cells persist |
| Chronic inflammation | Colchicine, IL-6/IL-1β blockade | Partial risk reduction only |
| Cap microcalcification | Statins (indirectly) | Not selectively targeted |
| Residual event risk | Best stack leaves 26–84% of excess risk | Never normalises |
| Dead myocardium | BioVAT-HF: +3.9 LVEF points | Nowhere near +25 needed |

Any single-agent "cure" fails because it addresses one row. REVERT-1 targets rows 2–6
with one lesion-gated particle, rows 1 and 7 with a systemic backbone, and explicitly
leaves row 8 to a different therapy.

---

## 3. Why the obvious combination fails

The naive invention is "put an anti-CD47 efferocytosis enhancer and an LXR agonist in
a nanoparticle." Two published failure modes kill it:

- **Systemic CD47 blockade causes anaemia.** CD47 is the "don't eat me" signal on red
  cells; blocking it systemically destroys them.
- **Systemic LXR agonism causes hepatic lipogenesis.** The first-generation LXR agonists
  raised triglycerides and caused steatosis, which is why they were abandoned.

And a third, subtler failure that motivated the whole design: **efferocytosis without
efflux is counterproductive.** ABCA1/ABCG1-deficient macrophages undergo oxidative
stress and apoptosis *after* engulfing apoptotic cells; the double knockout enlarges
the necrotic core, and cyclodextrin-mediated efflux rescues it
([Yancey et al., Circ Res 2010](https://www.ahajournals.org/doi/10.1161/circresaha.110.217281)).
Foam-cell transformation itself reduces efferocytic capacity, and uncleared apoptotic
bodies nucleate calcification
([Front Immunol 2023](https://doi.org/10.3389/fimmu.2023.1215612)).

That gives the design requirement: **the payload must be released only inside lesional
macrophages, and uptake must be delivered together with export capacity or not at all.**

---

## 4. The architecture

**Chassis:** a nanoparticle requiring *two* coincident recognition events to release —
an oxidised-phospholipid epitope plus a macrophage scavenger receptor — with
pH/enzyme-triggered uncoating in the phagolysosome. The AND gate is the invention's
load-bearing element: it is what makes an efficacious dose reachable without
systemic CD47 or LXR exposure.

Precedent for macrophage-restricted LXR delivery avoiding hepatic lipogenic gene
induction is real
([Biomaterials 2020](https://www.sciencedirect.com/science/article/abs/pii/S0142961220305792),
[PMC7530139](https://pmc.ncbi.nlm.nih.gov/articles/PMC7530139/)).

| Module | Target | Role | Precedent level |
|---|---|---|---|
| **M0** gate | oxPL + macrophage receptor | Confines release to plaque | Mouse pharmacologic |
| **M1** solubilise | Crystalline cholesterol | Dissolves crystals → generates the LXR ligand that fuels M2/M3 | Mouse pharmacologic |
| **M2** efflux | LXRβ → ABCA1/ABCG1 | Exports cholesterol so the efferocyte survives its meal | Mouse genetic |
| **M3** efferocytosis | Cleavage-resistant MerTK + local CD47/SIRPα F(ab')₂ | Clears the apoptotic cells forming the core | Mouse genetic |
| **M4** influx block | SR-A/CD36 siRNA | Stops re-loading the cell while unloading it | Mouse pharmacologic |
| **M5** resolution | NLRP3 / ADAM17 | Breaks the inflammation→MerTK-cleavage feed-forward loop | Mouse genetic |
| **M6** mineral | Microcalcification nucleation | Stops new cap micro-deposits; lets existing mineral mature | Human mechanism |
| **S1** driver (systemic) | Hepatic PCSK9 ± ANGPTL3 ± LPA | One-time editing shuts off supply | Human mechanism |
| **S2** tolerance (systemic) | ApoB/oxLDL CAR-Treg | Prevents recurrence without chronic immunosuppression | Mouse pharmacologic |

### The self-fuelling loop — the actual novelty

This is the part I have not found published as an integrated design:

```
   crystalline cholesterol in the plaque
                 │  M1 dissolves it
                 ▼
            oxysterols  ──► LXR activation
                 │                 │
                 │                 ├──► ABCA1/ABCG1  (M2: export)
                 │                 └──► MerTK        (M3: clearance)
                 ▼                            │
        corpse engulfment ◄───────────────────┘
                 │
                 └──► cholesterol load ──► exported by M2 ──► cell survives
                                                              and eats again
```

The plaque's own crystalline cholesterol becomes the signal that drives its removal.
Cyclodextrin dissolving crystals, raising oxysterols and reprogramming macrophages
toward regression in an LXR-dependent manner is published
([Zimmer et al., Sci Transl Med 2016](https://www.semanticscholar.org/paper/Cyclodextrin-promotes-atherosclerosis-regression-Zimmer-Grebe/4c3bf0237601aa3dc9756240b1aa5c1570868e06));
LXR coordinating ABCA1/ABCG1 with MerTK is published; the *closed loop with an AND gate*
is what I am proposing.

---

## 5. What the model says

An explicit-Euler toy model of four coupled states — free/crystalline cholesterol,
necrotic debris, microcalcification, viable phagocyte pool. Rate constants are
**assumed**, calibrated so an untreated lesion slowly progresses with a viable
phagocyte pool. The coupling is emergent, not asserted: engulfment adds cholesterol,
cholesterol above tolerance kills the cell, dead cells deposit both debris and
cholesterol.

### Ablation, 24 months (`revert1_ablation.csv`)

| Arm | Substrate removed | Necrotic core Δ | Phagocyte depletion |
|---|---|---|---|
| **revert1_full** | **+43.7%** | **−46.5%** | no |
| efflux + efferocytosis | +18.3% | −23.7% | no |
| two-pronged (published style) | +12.7% | +35.5% | no |
| efflux only | +5.3% | +35.5% | no |
| revert1 **ungated** | −4.7% | +12.7% | no |
| driver lowering only | −14.2% | +35.5% | no |
| untreated | −60.9% | +42.7% | no |
| **efferocytosis only** | **−102.3%** | +9.8% | **yes** |

Negative removal means the lesion grew. Two things to read carefully:

- **Efferocytosis alone is the worst arm in the study** — worse than doing nothing.
  The lesion roughly doubles.
- **It looks like a success at 12 months.** Necrotic core is −6.3% at month 12, then
  +9.8% by month 24 as cholesterol overload kills the phagocytes. A one-year study
  with a core-area endpoint would have called this a win. The model reproduces the
  trap, and `revert1_ablation.csv` carries both timepoints for exactly this reason.

### The ceiling (`revert1_ceiling.csv`, `revert1_horizon_sweep.csv`)

| Horizon | Central assumptions | Every module perfect |
|---|---|---|
| 12 mo | 30.1% | 55.3% |
| 24 mo | 43.8% | 74.9% |
| 36 mo | 52.9% | **80.4%** |
| 60 mo | 62.4% | 82.3% |
| 240 mo | 66.8% | **83.5%** |

Perfect modules cross 80% at **34.5 months**. Central assumptions **never** cross it.
Full eradication: **false**, at any horizon.

### Monte Carlo, 400 draws, uniform [0,1] on every module

Deliberately generous priors — this asks how often the architecture hits the bar even
if every module worked as well as it conceivably could.

- P(substrate removal ≥ 80%) = **0.0%**
- median removal **−4.9%**, 95th percentile **35.6%**, best draw **57.0%**

### Sensitivity (`revert1_sensitivity.csv`)

| Module | Spearman ρ |
|---|---|
| gate_specificity | **+0.56** |
| crystal_solubilisation | +0.43 |
| plasma_driver_reduction | +0.39 |
| efflux_gain | +0.31 |
| influx_block | +0.21 |
| resolution | −0.11 |
| nucleation_block | +0.10 |
| efferocytosis_gain | −0.04 |

**Efferocytosis strength is uncorrelated with clearance (ρ ≈ −0.04).** That is the
coupling thesis showing up independently in the sensitivity analysis: more uptake is
not the lever. Delivery and export are. Resolution comes out slightly negative
because it both boosts uptake and keeps loaded cells alive — when export is limiting,
that is mildly harmful. I did not tune this; it fell out.

---

## 6. Calcium: the design deliberately does *not* decalcify

The intuition "dissolve the calcium" is wrong and would destabilise plaques.
Microcalcification (0.5–15 µm, in the cap) raises local stress and tracks rupture;
dense macrocalcification is **stabilising**, and CAC *density* is inversely associated
with events ([Cells 2026](https://doi.org/10.3390/cells15161496),
[IJMS 2026](https://doi.org/10.3390/ijms27115113)).
Statins reduce the microcalcified fraction while overall density rises.

So M6 blocks *nucleation* of new micro-deposits and lets existing mineral mature.
The rupture-substrate metric in this model therefore **excludes mature
macrocalcification**. This is a proposed refinement to criterion 2 of the cure
framework, not a silent overwrite: burden reversal should score rupture-prone
substrate, not total plaque volume.

---

## 7. The heart is a separate problem

REVERT-1 does nothing for dead myocardium. Restoring an ischaemic cardiomyopathy from
LVEF 25% to 50% needs **+25 points**; the best human cell-therapy result is BioVAT-HF's
**+3.9 points** (n=20) — a **6.4× gap**. That belongs to a different arm (engineered
heart muscle, YAP101/SALVADOR-HF-style Hippo modulation), and pretending one particle
fixes both would be the central dishonesty this repository is built to avoid.

---

## 8. What a cure would require (back-calculated)

From `cure_target_product_profile.csv`, using this repo's own simulated numbers:

| Requirement | Target | Best current | Gap |
|---|---|---|---|
| Rupture-substrate removal | 80% | 43.1% | **1.9×** |
| Hazard ratio to normalise risk | **0.298** | best stack ≈0.6 | — |
| LVEF gain for the scar arm | +25 points | +3.9 points | **6.4×** |
| Randomised outcomes proof | required | — | cannot be back-calculated |

The HR of 0.298 is what it takes to leave ≤20% of excess risk in chronic coronary
syndrome (13.9% five-year risk against a 2.0% reference). No existing agent is close.

---

## 9. How to kill it — cheapest experiment first

| # | Claim | Falsified if | Cheapest test |
|---|---|---|---|
| 1 | Uptake without export is counterproductive | Anti-CD47 alone shrinks the core as well as co-delivery | 4-arm mouse study, 12 weeks Western diet; core area + TUNEL |
| 2 | Dissolved crystal cholesterol fuels the loop | No oxysterol/ABCA1/MerTK induction in human plaque macrophages | **Ex vivo human carotid endarterectomy explants, cyclodextrin vs vehicle, 48 h** |
| 3 | Lesion-restricted release is achievable | Lesional delivery ≤25%, or anaemia/steatosis persists | Large-animal biodistribution + Hb, reticulocytes, liver TG |
| 4 | Suppressing nucleation doesn't destabilise the cap | Cap thins or rupture features rise | OCT/NaF-PET + histology, micro-deposit size distribution |
| 5 | Tolerance prevents recurrence | Lesions regrow after withdrawal despite CAR-Treg | Withdrawal-and-rechallenge |

**Experiment 2 is the one to run first.** It is cheap, uses human tissue, needs no
particle engineering, and directly tests the self-fuelling claim on which the whole
architecture rests. If cyclodextrin does not induce oxysterol→LXR→ABCA1/MerTK in human
plaque macrophages, REVERT-1 is dead before anyone builds a nanoparticle.

---

## 10. Known failure modes I am not hiding

- **Cyclodextrin regression is contested.** A separate ApoE-null study found no lesion
  regression. Systemic cyclodextrin is ototoxic.
- **Cleavage-resistant MerTK increases liver fibrosis** in MASH models — a second
  reason the gate must hold, not just an anaemia argument.
- **Scavenger-receptor knockdown may impair host defence.** SR-A/CD36 are not
  atherosclerosis-specific.
- **Sustained inflammasome blockade carries infection risk**, visible in the human
  IL-6 data.
- **CAR-Treg Foxp3 instability** could convert tolerance into pathogenic effector cells.
- **A nine-module AND-gated particle is a manufacturing and regulatory nightmare.**
  CMC for this is plausibly harder than the biology.
- **The model is a toy.** Four state variables, assumed rate constants, no spatial
  structure, no shear stress, no immune trafficking, no pharmacokinetics. It can
  demonstrate that a coupling *matters*; it cannot predict an effect size.

---

## 11. Quarantine

This hypothesis is structurally incapable of contaminating the evidence tables:

- No REVERT-1 ID exists in `therapy_catalog.csv` (47 therapies) or
  `trial_anchors.csv` (48 anchors). Tests assert this.
- Outputs go to `results/hypotheses/`, never `results/solutions/`, and are outside the
  hash-verified evidence manifest.
- Every output row carries `status = invented_hypothesis_unvalidated`; a test reads
  all eight CSVs and asserts it.
- `DISCLAIMER.txt` sits in the output root and the report leads with the disclaimer.
- A test runs the real cure search and asserts no hypothesis ID can be selected.

---

## 12. Reproduce

```bash
python run_solutions.py hypothesis                      # 10 artefacts
python run_solutions.py hypothesis --draws 2000 --months 36
python -m pytest tests/test_hypothesis.py -q             # 51 tests
```

| Artefact | Contents |
|---|---|
| `revert1_ablation.csv` | 8 arms, both timepoints, depletion flag |
| `revert1_modules.csv` | 9 modules, precedent + failure mode + URL |
| `revert1_monte_carlo.csv` | 400 draws over all module strengths |
| `revert1_sensitivity.csv` | Spearman ranking |
| `revert1_horizon_sweep.csv` | Central vs perfect, 12–240 months |
| `revert1_ceiling.csv` | Asymptote, time-to-target, eradication flag |
| `revert1_kill_criteria.csv` | 5 falsifiers, cheapest first |
| `cure_target_product_profile.csv` | Back-calculated specification |

---

## 13. The honest summary

I could not invent a cure, and I do not think one can be invented from the current
evidence base by design alone. What I could do is identify where the leverage is:
**the coupling between corpse clearance and cholesterol export, and the delivery gate
that makes acting on it survivable.** Both are testable this year, the second one
cheaply, in human tissue.

If the coupling holds, the payoff is a plaque that clears itself using its own
cholesterol as the signal — a genuinely different mechanism from lowering LDL further.
If it does not hold, experiment 2 kills the idea for the cost of a few endarterectomy
specimens. Either outcome is worth more than another decimal place on LDL-C.

**This is a hypothesis. It is not a treatment. Do not act on it.**
