# The best-of-breed composite: strongest part per disease node, tested everywhere

> **STRICT-MODE HEADLINE RETRACTED — see
> [`evidence_layer_audit.md`](evidence_layer_audit.md).** The "strict 1/5" result
> below rested on a `risk_normalisation` score of 0.819 for `stable_ccs`, computed as
> a point estimate. Under the confidence intervals the source trials actually
> report, P(score ≥ 0.8) is **11.7%**, so that criterion does not survive and the
> best strict result is **0 of 5**. Four surrogate-mode and one exploratory-mode
> risk pass are also retracted; exploratory credit otherwise holds (18/19).


**Status: untrialled composite scenario.** Every combination here is assembled, not
studied. No composite is a cure, and the code cannot declare one.

```bash
python run_solutions.py composite      # -> results/composite/ (8 artefacts)
```

Code: `src/ihd_solutions/composite.py` · Tests: `tests/test_composite.py` (39 tests)
Tested across **19 populations × 3 credit modes × 2 variants = 114 criterion scorings**.

---

## 1. How it was assembled

For each population and each credit mode, greedily pick the component with the largest
**marginal evidence-graded event benefit**, subject to three constraints. Then run a
coverage pass for nodes that events cannot see.

The three rules that make "combine the best parts" legitimate rather than double counting:

| Rule | Why it exists |
|---|---|
| **Coverage combines by maximum, never summation** | Two half-effective regeneration therapies are not one restored heart |
| **One component per disease node** | Four anti-inflammatories are not four cures; infection risk compounds |
| **Same molecular target never stacked** | The catalogue puts PCSK9 *editing* and a PCSK9 *antibody* in different exclusive groups, so it would happily stack them. An explicit overlap rule blocks it |

And a fourth rule the single-agent scoring never needed:

> **An untrialled combination earns zero randomised-proof credit.** Every component can
> have phase-3 outcome evidence and the combination still has never been tested.
> Component proof does not transfer. This is why no composite here can be a cure.

### The coverage pass — a bug this design would otherwise have had

Selecting by event benefit is **structurally blind to the myocardium**: regeneration
carries zero event credit in this framework by design, so a greedy event optimiser never
picks it — silently dropping the node that binds the cure. The coverage pass attaches the
best available component for cure-relevant, event-neutral nodes at **zero event credit**,
recorded separately. A test asserts these additions do not move the event curve at all.

---

## 2. The recurring core

Components chosen most often across all populations and modes:

| Component | Node | Populations | Mean marginal benefit |
|---|---|---|---|
| canakinumab_tet2_guided | inflammation | 19 | 99.95/1000 |
| verve102_base_editing | durable_driver | 19 | 32.18/1000 |
| smoking_cessation | lifestyle | 18 | 19.42/1000 |
| aspirin | thrombosis | 18 | 10.98/1000 |
| icosapent_ethyl | triglyceride | 18 | 10.48/1000 |
| statin_high_intensity | ldl_synthesis | 12 | 27.43/1000 |
| retatrutide | metabolic | 11 | 27.05/1000 |
| bp_intensification | hemodynamic | 11 | 17.37/1000 |
| ezetimibe | ldl_absorption | 9 | 2.51/1000 |
| biovat_hf_patch | myocardium | coverage pass | **0.00/1000** |

The top two only appear in relaxed modes — that is the point of the next section.

---

## 3. Results across all 19 populations and all 3 modes

| Mode | Max criteria met | Driver control | Burden reversal | Myocardial | Risk normalisation | Proof | Best avoided |
|---|---|---|---|---|---|---|---|
| **strict** | **1 / 5** | 0.000 | 0.600 | 0.500 | **0.819 ✓** | 0.000 | 134.3/1000 |
| **surrogate** | **2 / 5** | **1.000 ✓** | 0.600 | 0.500 | **0.974 ✓** | 0.000 | 161.9/1000 |
| **exploratory** | **2 / 5** | **1.000 ✓** | 0.600 | 0.500 | **1.000 ✓** | 0.000 | 266.9/1000 |

**Composites declared a cure: 0.** Three criteria block in every population and every mode:
burden reversal (0.60), myocardial restoration (0.50), randomised proof (0.00).

### Best strict-mode results

| Population | Components | Baseline | Composite | Avoided /1000 |
|---|---|---|---|---|
| secondary_post_mi | 7 | 20.15% | **6.71%** | **134.34** |
| post_mi_hscrp_high | 7 | 22.12% | 8.80% | 133.21 |
| t2d_ascvd | 7 | 18.13% | 6.77% | 113.57 |
| statin_intolerant_ascvd | 5 | 18.13% | 7.66% | 104.71 |
| ascvd_high_triglycerides | 6 | 18.13% | 7.98% | 101.50 |
| post_mi_pci | 5 | 22.12% | 12.30% | 98.18 |

### The same population, three modes — watch the composition change

**stable_ccs**, baseline 13.93%:

- **strict → 4.14%**: colchicine + statin + aspirin + bp + alirocumab + ezetimibe + biovat patch
- **surrogate → 2.84%**: colchicine + **verve102 editing** + retatrutide + aspirin + bp + bempedoic + diet + ezetimibe + biovat patch
- **exploratory → 2.45%**: **canakinumab (TET2-guided)** + smoking cessation + verve102 + icosapent + aspirin + ezetimibe + biovat patch

Strict mode **cannot** include gene editing — no outcomes evidence — so durable driver
control scores 0. Relaxed modes admit it and driver control reaches 1.0. **That difference
is the honest measure of how much of the composite depends on evidence we do not have.**

Note what the assembler refused on its own: **no Lp(a) agent was ever selected**, because
Lp(a)HORIZON failed and the neutral anchor holds it at null. No hand-curation needed.

---

## 4. What genuinely improved

**Risk normalisation now crosses the bar.** The previous best stack removed 74% of excess
risk in stable coronary disease (score 0.74, below the 0.8 threshold). The composite reaches
**0.819 strict / 0.974 surrogate / 1.000 exploratory** — the first time any configuration in
this repository has cleared that criterion. It got there by adding alirocumab and ezetimibe
on top of the four-drug backbone, each accepted only for measured marginal benefit.

That is a real gain, and it is the *only* criterion that moved. Assembling the best proven
parts does not create burden reversal or new myocardium.

---

## 5. What the invented parts add

**Nothing to any score.** The REVERT-1 clearance loop and the CAR-Treg lock are recorded as
the nodes the catalogue cannot serve, tagged `invented_component_unvalidated`, and carry
zero evidence credit. A test asserts the `with_invented` variant produces *identical*
criterion scores and an identical event curve. Their hypothetical effect lives in
`results/hypotheses/` and `results/discovery/`, never here.

This is deliberate: it would be easy to let an invented module quietly lift a score into
passing territory. The structure makes that impossible.

---

## 6. Honest limitations

- **Greedy is not optimal.** Forward selection can miss a pair that only helps together.
  The repo's exhaustive enumeration (753 regimens) remains the check on that.
- **The composite is untrialled.** Nine components with nine safety profiles; drug-drug
  interactions, cumulative immunosuppression and adherence at this complexity are not
  modelled. `burden_score` is a crude complexity proxy, not tolerability.
- **Same-target rules are hand-specified.** `TARGET_OVERLAPS` encodes PCSK9 and
  statin/bempedoic overlaps. A target pair not listed there would slip through.
- **Node assignment is a modelling choice.** `DISEASE_NODES` is explicit and tested for
  catalogue coverage, but a different granularity gives a different composite.
- **Empty composites are real.** In `non_ischemic_cardiomyopathy` and similar populations no
  agent earns strict credit; the honest output is an empty regimen, not a fabricated one.

---

## 7. Reproduce

```bash
python run_solutions.py composite
python run_solutions.py composite --populations stable_ccs secondary_post_mi
python run_solutions.py composite --horizon-years 10 --max-components 6
python -m pytest tests/test_composite.py -q
```

| Artefact | Contents |
|---|---|
| `composite_assembly.csv` | Every greedy acceptance, its marginal benefit, rejections |
| `composite_test_matrix.csv` | 57 rows: 19 populations × 3 modes, with composition |
| `composite_criteria.csv` | 114 rows: × 2 variants, all five criteria, `is_cure=False` |
| `composite_component_frequency.csv` | The recurring best-of-breed core |
| `composite_leave_one_out.csv` | What each part actually contributes |
| `composite_node_map.csv` | Node → exclusive groups, invented-only nodes flagged |
| `reports/composite_report.txt` | Plain-text bottom line |
| `DISCLAIMER.txt` | The disclaimer, in the output root |

---

## 8. Bottom line

Taking the best working part from each node produces a measurably better regimen —
**134 events avoided per 1000** in secondary prevention, and the first configuration in this
repository to normalise event risk. It is still not a cure, and the reason is precise:
**the parts that would make it one do not exist.** Burden reversal caps at 0.60, myocardial
restoration at 0.50, and an untested combination has no proof at all.

The composite is the strongest thing the current evidence supports. The gap between it and
a cure is exactly the gap quantified in `docs/cure_discovery_plan.md`.
