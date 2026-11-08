# Finding a cure in simulation only: what is achievable, what is missing, and the full programme

**Read this first.** A cure cannot be discovered in simulation. Not because the models
are weak, but because this repository's own definition of a cure has five criteria and
one of them — `randomised_proof` — requires randomised human outcome evidence. No
amount of compute produces that. Any pipeline that reported "cure found" from
simulation alone would be committing a definitional trick.

What simulation *can* discover is a **specification**: the smallest set of capabilities,
with numbers, that would satisfy every criterion a model can legitimately reason about,
plus a ranked list of the measurements that would most change the answer. That is what
this programme produces, and it is now executable:

```bash
python run_solutions.py discover        # -> results/discovery/ (10 artefacts)
```

Code: `src/ihd_solutions/discovery.py` · Tests: `tests/test_discovery.py` (43 tests)

---

## Part 1 — Which is the best candidate still?

There are three different questions hiding in "best", and they have different answers.

### Best evidence-backed single agent: **evolocumab** (readiness 0.360)

| Therapy | Readiness | Driver control | Burden reversal | Myocardial | Risk normalisation | Randomised proof |
|---|---|---|---|---|---|---|
| **evolocumab** | **0.360** | 0.20 | 0.55 | 0.00 | 0.050 | 1.00 |
| ezetimibe | 0.345 | 0.10 | 0.60 | 0.00 | 0.023 | 1.00 |
| statin (high intensity) | 0.339 | 0.10 | 0.50 | 0.00 | 0.096 | 1.00 |
| ctx310 (ANGPTL3 editing) | 0.260 | **1.00** | 0.00 | 0.00 | 0.000 | 0.30 |
| verve102 (PCSK9 editing) | 0.260 | **1.00** | 0.00 | 0.00 | 0.000 | 0.30 |

Every one of them satisfies exactly **one** criterion. Evolocumab wins on proven
event reduction plus GLAGOV plaque regression. Note what wins instead on the column
that matters most for a cure: **one-time editing is the only modality scoring 1.00 on
durable driver control** — the only capability in the entire catalogue that is already
solved in principle.

### Best actual regimen: the stack, not any drug

Best strict stack in chronic coronary syndrome: **13.93% → 5.09%** five-year risk, i.e.
**74% of excess risk removed**. That is the closest anything gets to risk
normalisation. But it is wildly population-dependent:

| Population | Baseline 5-y risk | Best stack | Excess risk remaining |
|---|---|---|---|
| stable_ccs | 13.93% | 5.09% | **26.0%** |
| hefh_primary | 9.52% | 4.53% | 33.9% |
| secondary_post_mi | 20.15% | 8.23% | 34.4% |
| post_mi_pci | 22.12% | 12.30% | 51.2% |
| ascvd_ckd_inflammation | 25.92% | 17.11% | 63.2% |
| **ischemic_cm_heart_failure** | 32.97% | 27.97% | **83.9%** |

**This table is the real answer.** "Cure IHD" is two problems with very different
difficulty. In a patient with living myocardium we are ~74% of the way to normalising
risk. In a patient whose myocardium is already dead we are ~16% of the way, because
removing the cause does not rebuild the muscle.

### Best cure-directed candidate: **REVERT-1**, and it is not close to a cure

The invented architecture clears 43.7% of rupture-prone substrate under central
assumptions, 83.5% with every module perfect, and never eradicates. It addresses the
plaque, not the heart.

**So the honest answer: the best candidate is a stack for prevention plus a
regeneration capability that does not exist.** No single agent is a cure candidate,
and the highest-leverage validated mechanism (one-time editing) scores zero on the two
criteria that actually bind.

---

## Part 2 — What is missing? Three answers, only one of which is a research problem

### (a) A missing *capability* — this is the real blocker

Back-calculated minimum vector, now computed rather than asserted
(`discovery_required_vector.csv`):

| Capability | Required | Best demonstrated | Gap | Reachable today? |
|---|---|---|---|---|
| **lvef_gain_pp** | **20.0 points** | **3.9 points** (BioVAT-HF, n=20) | **×5.13** | **no** |
| substrate_removal | 0.800 | 0.431 | ×1.86 | no |
| hazard_ratio | ≤0.298 | 0.360 | ×1.21 | no |
| driver_durability | 0.800 | 1.000 | ×0.80 | **yes — solved** |

Myocardial restoration is **5.1× beyond anything demonstrated**, and no parameter
tuning creates a capability that does not exist. This is not a modelling gap.

### (b) A missing *measurement* — which numbers we do not know

The coupled-clearance model runs on assumed rate constants. The decision-relevant
unknowns, in order of how much they move the binding margin
(`discovery_value_of_information.csv`, variance decomposition of the worst-criterion
score):

| Unknown | Share of decision variance |
|---|---|
| substrate_removal (clearance coupling + gate efficiency) | **18.4%** |
| hazard_ratio (stack-level risk normalisation) | 15.2% |
| lvef_gain_pp (remuscularisation ceiling) | 15.2% |
| driver_durability | 14.1% |

Note the split personality, which is the most useful finding here: **the heart is the
biggest capability gap (×5.13), but plaque clearance is the tightest engineering
tolerance.** Only **4.19%** of the conceivable substrate-removal range passes, because
REVERT-1's own asymptotic ceiling (83.5%) sits barely above the 80% bar. There is
almost no margin for error in the plaque arm.

### (c) A missing *evidence type* — structurally unobtainable in simulation

`randomised_proof` cannot be satisfied in silico. Full stop.

### The feasibility frontier — the headline number

Searching the same capability space under two bound sets:

| Bound set | Exact feasible volume | Best corner passes? | Blocking at that corner |
|---|---|---|---|
| **Demonstrated** | **0** (exactly) | no | burden reversal, myocardial restoration, risk normalisation |
| **Conceivable** | **3.686 × 10⁻⁴** (0.0369%) | **yes** | none |

The feasible region is **empty** using only demonstrated capabilities and **non-empty**
using generous ones. *The difference between those two boxes is precisely the missing
discovery.* And even in the generous box the passing region is a vanishingly thin
corner: every capability must sit near its theoretical maximum **simultaneously**.

Per-capability passing windows in the conceivable box: substrate removal **4.19%**,
LVEF gain 20.0%, hazard ratio 21.99%, durability 20.0%.

### Which capabilities are indispensable?

Revert-one ablation inside a passing specification (`discovery_capability_ablation.csv`):

- revert `driver_durability` → **still passes** → already solved, not part of the problem
- revert `substrate_removal` → fails → **indispensable**
- revert `lvef_gain_pp` → fails → **indispensable**
- revert `hazard_ratio` → fails → **indispensable**

**Three discoveries are required. Durability is not one of them.**

---

## Part 3 — The full simulation-only programme

Eight stages, each with an explicit exit gate and an explicit statement of what failure
means. Six are executable today; two need model structure that does not exist yet.
Machine-readable in `discovery_stages.csv`.

### Stage 0 — Formalise the target ✅ executable

**Question.** What would count, and which criteria can a model address?
**Method.** Freeze the five criteria. Mark `randomised_proof` simulation-inadmissible.
Set population-specific risk targets rather than one global number.
**Exit gate.** A written target with four admissible criteria and one explicitly excluded.
**If it fails.** Nothing downstream is interpretable; any cure claim is definitional fraud.

> Already done. The `Specification` class deliberately exposes `is_candidate`, never
> `is_cure`, and a test asserts no `is_cure` attribute exists.

### Stage 1 — Feasibility frontier ✅ executable

**Question.** Does any capability vector satisfy the admissible criteria?
**Method.** Search the demonstrated and conceivable boxes; compute the feasible volume
analytically (exact for a box under flat priors — sampling a 0.037% region is noise).
**Exit gate.** Non-empty under conceivable bounds, empty under demonstrated.
**If it fails.** If empty under *both*, the architecture is wrong. Redesign; do not tune.

> Result: 0 versus 3.686 × 10⁻⁴. Gate passed.

### Stage 2 — Gap quantification ✅ executable

**Question.** How far is each capability from its requirement?
**Method.** Solve the minimum vector analytically; divide by the best demonstrated
value; attach a cited provenance to every number.
**Exit gate.** Every gap is a factor with a source, not an adjective.
**If it fails.** A gap without provenance cannot prioritise work.

> Result: ×5.13, ×1.86, ×1.21, ×0.80. Gate passed.

### Stage 3 — Value of information ✅ executable

**Question.** Which unknown, if measured, would most change the answer?
**Method.** First-order Sobol-style variance decomposition of the binding margin:
Var_i(E[Y|X_i]) / Var(Y), on `worst_score`.
**Exit gate.** Capabilities ranked by decision variance resolved.
**If it fails.** Experiments get chosen by fashion instead of decision impact.

> Result: clearance coupling 18.4% > risk 15.2% ≈ myocardium 15.2% > durability 14.1%.

### Stage 4 — Minimal capability set ✅ executable

**Question.** Which capabilities are indispensable versus already solved?
**Method.** Start from a passing vector; revert one capability at a time to its
demonstrated value; whatever breaks is indispensable.
**Exit gate.** An irreducible set of required discoveries.
**If it fails.** Effort gets spent re-solving what is solved.

> Result: {substrate clearance, remuscularisation, stack risk normalisation}. Durability excluded.

### Stage 5 — Adversarial falsification ⚠️ needs new model structure

**Question.** Does the candidate survive hostile but plausible assumptions?
**Method.** Re-run the specification under:
1. pessimistic coupling (efflux capacity halved, load tolerance halved);
2. **added pharmacokinetics** — clearance half-life, lesion residence, repeat-dose accumulation;
3. **withdrawal and recurrence** — remove the chassis at month 24 and follow to month 60;
4. **spatial cap mechanics** — finite-element stress on the micro-deposit distribution
   instead of a scalar burden;
5. immune trafficking and competing mortality;
6. the five REVERT-1 kill criteria as parameter perturbations.

**Exit gate.** A feasible region survives across structural perturbations, not just one
model form.
**If it fails.** The specification is an artefact of one model structure. Discard it.

> **This is the stage most likely to kill the current answer**, because the 4.19%
> passing window on substrate clearance leaves almost no robustness margin.

### Stage 6 — Identifiability audit ⚠️ needs new model structure

**Question.** Are the governing parameters identifiable from obtainable data?
**Method.** Profile-likelihood and sloppy-model analysis on the coupled clearance model.
For each decision-relevant parameter, ask whether any feasible measurement pins it.
Flag structural non-identifiability (parameters that only ever appear as a product).
**Exit gate.** Every decision-relevant parameter has an identifiable measurement route.
**If it fails.** The programme is unfalsifiable. Stop and reformulate.

> Honest expectation: several rate constants in the current four-state model are
> sloppy. This stage will probably force a reduced model, which is a good outcome.

### Stage 7 — Specification handoff ✅ executable

**Question.** What is the cheapest real sequence that could convert this to evidence?
**Method.** Emit the target product profile, the kill criteria, and the experiment list
ordered by Stage 3's information ranking.
**Exit gate.** A specification a laboratory can act on, with `randomised_proof` still open.
**If it fails.** The work stays inert modelling.

---

## Part 4 — What the programme concludes, and what it refuses to conclude

**Concludes:**

1. Driver durability is solved in principle. Stop spending discovery effort there.
2. Three capabilities are indispensable; the feasible region requiring all three is
   0.037% of the generous design space. A cure is not one clever idea away — it needs
   three simultaneous near-maximal capabilities.
3. The biggest *capability* gap is remuscularisation (×5.13). The tightest *tolerance* is
   plaque clearance (4.19% passing window). These demand different kinds of work.
4. Prevention and restoration are different problems: 74% of excess risk removable in
   stable coronary disease versus 16% in ischaemic cardiomyopathy.
5. The first measurement to buy is the clearance coupling — human plaque explants,
   cyclodextrin versus vehicle, oxysterol and ABCA1/MerTK readout. Cheap, human tissue,
   no particle engineering, highest decision impact.

**Refuses to conclude:**

- That a cure has been found. `randomised_proof` is untouched and untouchable here.
- That REVERT-1 works. It is an invented hypothesis with assumed rate constants.
- That the feasible corner is reachable. It requires near-perfect delivery, near-perfect
  payloads, and a remuscularisation capability with **no mechanistic precedent at all**.
- Anything patient-specific. This is scenario arithmetic, not prediction.

**The one-sentence version:** simulation can tell you exactly what a cure would have to
do, prove that nothing demonstrated can do it, and rank the experiments that would
change the answer — but the cure itself has to be found in tissue and in trials.

---

## Reproduce

```bash
python run_solutions.py discover                                  # 10 artefacts
python run_solutions.py discover --baseline-risk 0.3297           # ischaemic cardiomyopathy
python run_solutions.py discover --draws 20000 --seed 42
python -m pytest tests/test_discovery.py -q                        # 43 tests
```

| Artefact | Contents |
|---|---|
| `discovery_required_vector.csv` | Analytic minimum capability vector |
| `discovery_feasibility.csv` | Both bound sets, corner scores, blocking criteria |
| `discovery_feasible_volume.csv` | Exact volumes and per-capability passing windows |
| `discovery_capability_gaps.csv` | Gap factors with provenance |
| `discovery_value_of_information.csv` | Decision variance resolved per capability |
| `discovery_capability_ablation.csv` | Indispensable versus already-solved |
| `discovery_stages.csv` | The eight stages, gates and failure meanings |
| `discovery_samples.csv` | 8 000 sampled capability vectors, tagged |
| `reports/discovery_report.txt` | Plain-text bottom line |
| `DISCLAIMER.txt` | The disclaimer, in the output root |

Everything is quarantined exactly like the REVERT-1 hypothesis: separate directory,
outside the hash-verified evidence manifest, every row tagged
`in_silico_specification_not_evidence`, and tests asserting the output can never be
called a cure.
