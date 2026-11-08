# The 2026 cure search: everything I could find, encoded, simulated and tested

**Bottom line up front.** I searched the current literature and trial record for
every plausible route to *curing* ischaemic heart disease, added all of them to the
simulator, and ran them. **No intervention, and no combination of interventions,
satisfies a defensible definition of a cure.** Three requirements are not met by
anything in human evidence today: removing existing plaque burden (not just
softening it), replacing lost myocardium (not just gaining a few LVEF points), and
returning event risk to a low-risk reference. The strongest genuinely cure-adjacent
lever is not a new drug target at all — it is **starting durable LDL-C reduction
decades earlier**.

Everything below is reproducible: `python run_solutions.py simulate` regenerates
all 33 artefacts, and `pytest -q` runs 425 tests over the result.

---

## 1. What changed in the evidence since the last build

The most consequential finding is a **failure**, not a success.

**Lp(a)HORIZON missed its primary endpoint.** On 4 September 2026 Novartis
announced that pelacarsen lowered Lp(a) but did not reduce 4-point MACE in 8,323
patients with established cardiovascular disease
[4](https://www.biopharmawatch.com/news/NVS/novartis-lpa-horizon-phase-iii-results-pelacarsen).
The published commentary is explicit that the level fell and the risk did not, and
that hazard ratios and confidence intervals were still unpublished at the time of
writing [1](https://pubmed.ncbi.nlm.nih.gov/42880684/). Roughly 80% Lp(a) reduction
was not enough
[5](https://www.unibestpharm.com/news/novartis-pelacarsen-lp-a-horizon-phase-3-miss-cardiovascular-2026).

I re-encoded pelacarsen from "pending" to **negative randomised topline**, held at
HR = 1.0. It now behaves as a negative control, exactly like CIRT and ZEUS. The
deeper siRNA agents remain unproven: olpasiran lowers Lp(a) ~94–97% with
OCEAN(a)-Outcomes (n = 7,297) due around December 2026, and lepodisiran's
ACCLAIM-Lp(a) (n = 17,300) reads out ~2029
[3](https://www.fishtownmedicine.com/cardiovascular-risk/lowering-lipoprotein-a-new-drugs).
Zerlasiran is phase-3-ready with no active outcomes trial
[1](https://elevatedcholesterol.com/lpa-lowering-drugs.html). **None of them earns
event credit in any mode**, and after HORIZON that refusal is no longer merely
precautionary.

Other updates now encoded:

| Finding | Status |
| --- | --- |
| CTX310, one-time CRISPR editing of *ANGPTL3* | Phase 1a, n = 15: ANGPTL3 −79%, LDL-C −53%, TG −48%, sustained to one year [2](https://crisprtx.com/about-us/press-releases-and-presentations/crispr-therapeutics-presents-phase-1a-data-for-ctx310-demonstrating-deep-and-durable-angptl3-editing-triglyceride-and-ldl-lowering-at-esc-congress-2026) |
| VERVE-102 *PCSK9* base editing | Phase 1b, n = 35: PCSK9 −88%, LDL-C −62%, durable to 18 months; no DLTs [1](https://www.nejm.org/doi/full/10.1056/NEJMoa2601283) |
| Enlicitide (Lipfendra) | FDA-approved 16 July 2026, first oral PCSK9 inhibitor; LDL-C −56%; outcomes trial still running [3](https://www.pharmacytimes.com/view/fda-approves-enlicitide-the-first-oral-pcsk9-inhibitor-for-high-cholesterol) |
| Inclisiran outcomes | ORION-4 follow-up completes 2026, results **early 2027** [5](https://pubmed.ncbi.nlm.nih.gov/42567431/) |
| Obicetrapib (CETP) | BROADWAY n = 2,530, LDL-C −29.9%; events HR 0.79 (0.54–1.15), unpowered; PREVAIL n = 9,541 primary completion Nov 2026 [1](https://www.radcliffecardiology.com/news/obicetrapib-significantly-reduces-ldl-cholesterol-high-risk-patients) [3](https://cetpinhibition.org/) |
| Tirzepatide | SURPASS-CVOT n = 13,299: MACE-3 HR 0.92 (0.83–1.01) **versus dulaglutide**, non-inferior, not superior [5](https://www.hcplive.com/view/tirzepatide-gets-fda-approval-to-cut-cardiovascular-risk-in-type-2-diabetes) |
| Plozasiran / olezarsen (APOC3) | TG −79/81% and up to −72%; pancreatitis RR 0.22 and −85%. **Pancreatitis is not MACE** [1](https://www.escardio.org/news/press/press-releases/plozasiran-significantly-reduces-triglyceride-levels-in-patients-with-severe-hypertriglyceridaemia/) [3](https://news.gastro.org/issues/2026/february/olezarsen-shows-promise-for-severe-hypertriglyceridemia-management) |
| BioVAT-HF engineered heart muscle | NEJM 2026, n = 20: wall +4.5 mm, LVEF +3.9 pp, KCCQ +6.7 at 3 months; remuscularisation confirmed in an explanted heart; 3 deaths and 3 VT judged unrelated [3](https://www.springermedicine.com/chronic-heart-failure/heart-surgery/biovat-bioengineered-heart-muscle-tissue-transplant-advanced-hf/52528448) |
| YAP101 (AAV9-SAV1-shRNA) | SALVADOR-HF, first-in-human cardiac *regeneration* gene therapy, started June 2026; pig EF +14% [3](https://www.jfmartinlab.com/blog/our-gene-therapy-work-featured-in-nature) [1](https://www.dongascience.com/en/news/78158) |
| AB-1002 / GenePHIT | n = 173, **non-ischaemic** cardiomyopathy only; results H1 2027 [2](https://www.bayer.com/media/en-us/askbio-presents-baseline-characteristics-from-phase-2-heart-failure-gene-therapy-trial-of-umiposgene-parvec-ab-1002-at-european-society-of-cardiology-congress/) |
| Genotype-guided anti-inflammatory therapy | CANTOS exploratory substudy: TET2-CHIP carriers (n = 103) HR ≈ 0.38; phase 2a NLRP3 trials under way [2](https://www.frontiersin.org/journals/cardiovascular-medicine/articles/10.3389/fcvm.2026.1796328/full) [5](https://pmc.ncbi.nlm.nih.gov/articles/PMC13431735/) |
| CAR-Treg against oxLDL | 2026 proof of concept: less foam-cell formation in vitro, less plaque in mice — atherosclerosis as an autoimmune disease [3](https://academic.oup.com/cardiovascres/advance-article/doi/10.1093/cvr/cvag020/8471361) |
| PCSK9 vaccine | Mice LDL-C −29%; **macaques showed antibodies but no lipid change** [3](https://www.drugtargetreview.com/peptide-vaccine-targeting-pcsk9-proves-effective-in-cardiovascular-disease-models/2135854.article) |
| Coronary sinus reducer | REMEDY-PILOT (TCT 2026, n = 54): no better than sham in ANOCA/CMD [3](https://www.healio.com/news/cardiology/20261104/coronary-sinus-reducer-no-better-than-sham-procedure-in-angina-without-obstruction) |
| Retatrutide | Phase 2: weight −28.7%, SBP −14 mmHg; no outcomes data [3](https://www.peptidejournal.org/research/2025-2026-peptide-research-highlights) |

The catalogue is now **47 therapies, 48 trial anchors, 19 populations**.

---

## 2. New interlocks this sweep forced me to build

Adding these candidates exposed four ways the simulator could have started
lying. Each is now a tested refusal:

1. **Animal data can never become human benefit.** The PCSK9 vaccine and CAR-Treg
   look attractive on paper. A `MECHANISTIC` tier is now barred from biomarker
   transport in *every* mode. Simulated incremental benefit: **0.0 events per 1000
   in strict, surrogate and exploratory**.
2. **An active-comparator hazard ratio is not a placebo-relative effect.**
   SURPASS-CVOT's 0.92 is against dulaglutide, which is itself cardioprotective. A
   new `RANDOMIZED_ACTIVE_COMPARATOR` tier is refused with that exact reason.
3. **Mendelian randomisation is not a treatment effect.** The lifelong genetic
   contrast (RR 0.46 per mmol/L) is used *only* by the lifetime scenario model and
   is excluded from event credit and from per-unit pathway lookup.
4. **Imaging never touches the event model.** The new ASTEROID/GLAGOV/FLOWPROMOTE
   anchors initially shadowed the statin's CTT anchor and silently zeroed statin
   benefit — the test suite caught it. Plaque-imaging anchors are now filtered out
   of credit selection entirely and feed the regression model instead.

One legitimate inversion is now documented rather than hidden: obicetrapib scores
**higher in surrogate mode than in exploratory mode**. In chronic coronary
syndrome, surrogate mode transports the CTT LDL relationship to it; exploratory
mode admits its own neutral BROADWAY event analysis, which is held at the null and
displaces the transport (no double counting). The simulation emits
`neutral_anchor_supersedes_biomarker_transport` so the reason is auditable.

---

## 3. Defining "cure", then searching for it

`src/ihd_solutions/cure.py` states five criteria and scores all 47 therapies from
the encoded evidence. A criterion counts as met at ≥ 0.8. Criteria combine by
**maximum, never by summation** — two half-effective regeneration procedures do
not add up to a restored heart.

| Criterion | Best score any therapy achieves | Why it is capped |
| --- | --- | --- |
| Durable driver control | **1.00** | one-time editing (VERVE-102, CTX310) — the only modality that clears the bar |
| Burden reversal | **0.60** | regression is partial and calcified volume rises |
| Myocardial restoration | **0.50** | BioVAT gives +3.9 LVEF points, not normal myocardium |
| Risk normalisation | **0.12** single agent; **0.74** best stack | treated risk never reaches a low-risk reference |
| Randomised proof | **1.00** | statins, ezetimibe, PCSK9 antibodies, colchicine, aspirin, BP |

The combination search returns `satisfied: false` with binding constraints
`burden_reversal`, `myocardial_restoration`, `risk_normalisation`. Even at a
relaxed threshold of 0.5 it still fails on risk normalisation. The highest
cure-readiness index in the entire catalogue is **0.36** (evolocumab); the one-time
editors score 0.26 because they clear durability but have no outcomes proof.

This is a statement about the encoded evidence, not a proof that a cure is
impossible.

---

## 4. How close does plaque regression actually get?

This is the part of "reversal" that is real. FLOWPROMOTE, published 30 October
2026, randomised statin-naive patients and imaged the whole coronary tree:
rosuvastatin + ezetimibe reduced total plaque volume 21.1% versus 13.6%, and
noncalcified plaque 28.2% versus 20.0%, with soft plaque regressing in 51 of 52
patients at an achieved LDL-C of 42.5 mg/dL — while **calcified volume did not
fall** [5](https://www.dralo.net/blog/reverse-soft-plaque)
[4](https://x.com/CMEINDIA1/status/2117270752689082643). GLAGOV and ASTEROID show
the same pattern: about two-thirds of patients regress, and the proportion rises
to 81.2% in those already below 70 mg/dL
[3](https://deo-medical.com/blog/plaque-regression/).

`regression.py` fits those published points rather than assuming a shape. The fit
puts the switch from progression to regression at **74 mg/dL**, independently
matching the 70–80 mg/dL range in the literature. Simulated over 2 years on top of
statin + aspirin in chronic coronary syndrome:

| Add-on | Achieved LDL-C | Lipid-rich pool removed | Total burden removed |
| --- | --- | --- | --- |
| VERVE-102 | 19.7 mg/dL | 43.1% | 21.5% |
| Enlicitide | 22.5 mg/dL | 42.2% | 21.1% |
| Evolocumab | 22.7 mg/dL | 42.1% | 21.1% |
| Obicetrapib | 34.5 mg/dL | 38.4% | 19.2% |

Note the gap between the columns: **the dangerous fraction shrinks fast, the
artery does not clean out.** Push the model to 30 years at LDL-C 20 mg/dL and the
lipid pool is >90% gone while total burden falls less than 70%, because dense
calcium remains. `eradicated` is `False` by construction, and every projection
past two years is flagged as an extrapolation.

---

## 5. The strongest cure-adjacent lever: start decades earlier

Mendelian randomisation is the reason this is worth modelling at all. Carriers of
*PCSK9* loss-of-function variants have ~28% lower LDL-C and **88% lower coronary
risk** [1](https://www.cardiologywiki.com/wiki/concepts/pcsk-9-inhibitors/), and a
1 mmol/L lower LDL-C sustained from early life associates with ~54% lower CHD risk
— far more than the 22% a five-year statin trial delivers for the same absolute
difference [3](https://www.curingheartdisease.com/the-genetic-architecture-of-cardiovascular-protection/).
The difference is cumulative exposure, not potency.

So the honest question for a one-time editor is not "how much does it lower
LDL-C?" but "how much of a lifetime does it cover?" Scaling the genetic contrast
by the share of the adult exposure window covered:

| Start age | Window covered | Lifetime risk (vs 20%) | Reduction | Events avoided /1000 |
| --- | --- | --- | --- | --- |
| 30 | 83% | 4.2% | 79% | 158 |
| 40 | 67% | 5.7% | 71% | 143 |
| 50 | 50% | 7.9% | 60% | 121 |
| 60 | 33% | 10.9% | 46% | 91 |
| 70 | 17% | 14.8% | 26% | 52 |

(VERVE-102's 62% LDL-C reduction; an editing-plus-statin stack at 81% reaches 87%
lower lifetime risk from age 30.) **This is a genetic scenario, not a trial
result** — it assumes the reduction holds for the whole window and ignores
procedural risk, off-target editing risk and competing mortality. Every row is
labelled `mendelian_randomisation_scenario` in the output, and tests enforce that
label.

---

## 6. What is left over after the best available stack

Five-year risk under the best strict-evidence regimen found, against a low-risk
reference:

| Population | Untreated | Best stack | Reference | Excess risk remaining |
| --- | --- | --- | --- | --- |
| Chronic coronary syndrome | 13.9% | 5.1% | 2.0% | **26%** |
| Secondary post-MI | 20.1% | 8.2% | 2.0% | 34% |
| HeFH primary prevention | 9.5% | 4.5% | 2.0% | 34% |
| Post-MI after PCI | 22.1% | 12.3% | 2.0% | 51% |
| ASCVD + CKD + inflammation | 25.9% | 17.1% | 2.0% | 63% |
| Ischaemic CM with heart failure | 33.0% | 28.0% | 2.0% | **84%** |

Even the best case leaves a quarter of the excess risk in place. In ischaemic
cardiomyopathy, four-fifths remains — which is precisely why regeneration is the
missing criterion rather than an optional extra.

---

## 7. What would actually have to happen

A concrete, falsifiable agenda, in the order the ceilings bind:

1. **Risk normalisation.** Show that a treated cohort's event rate equals a
   low-risk reference, not merely a lower rate than placebo. Nothing currently
   claims this.
2. **Burden reversal.** Serial imaging showing reduction of *total* and calcified
   burden, not just the noncalcified fraction. The new CCTA programmes with
   olpasiran and lepodisiran are the nearest test
   [1](https://elevatedcholesterol.com/lpa-lowering-drugs.html).
3. **Myocardial restoration.** Controlled trials of durable, arrhythmia-free
   remuscularisation with normalised function. BioVAT is the strongest start but
   is single-arm, n = 20, and requires immunosuppression
   [3](https://www.springermedicine.com/chronic-heart-failure/heart-surgery/biovat-bioengineered-heart-muscle-tissue-transplant-advanced-hf/52528448).
4. **Outcomes proof for one-time editing**, plus off-target assays that
   double-strand-break methods cannot provide. Until then the editors are
   durability without proof.
5. **Antigen-specific tolerance.** If atherosclerosis behaves as an autoimmune
   disease, CAR-Treg or ApoB tolerogenic vaccination could address the driver
   rather than the substrate — currently mouse-only
   [1](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2026.1903006/full).
6. **Genotype-guided inflammation.** A prospective TET2-CHIP-stratified
   randomised trial would convert the CANTOS post-hoc signal into evidence
   [5](https://pmc.ncbi.nlm.nih.gov/articles/PMC13431735/).

Trials that could move rows in this registry within ~24 months: OCEAN(a)-Outcomes
(Dec 2026), PREVAIL (Nov 2026 completion), ORION-4 (early 2027), CORALreef
Outcomes, GenePHIT (H1 2027), and the phase-2 readouts for VERVE-102 and CTX310.

---

## 8. Reproducing this

```bash
python run_solutions.py simulate                      # 33 artefacts, ~10 s
cat results/solutions/reports/cure_search_report.txt   # the verdict, with caveats
```

New outputs: `cure_criteria_scores.csv`, `cure_combination_search.json` (+ a
relaxed-threshold variant), `lifetime_exposure_scenarios.csv`,
`regression_calibration.csv`, `regression_calibration_fit.json`,
`plaque_imaging_evidence.csv`, `plaque_regression_grid.csv`, `residual_risk.csv`,
and `reports/cure_search_report.txt`.

New tests: `tests/test_regression.py` and `tests/test_cure.py`, plus integrity
checks that the committed cure outputs still refuse a cure claim. The suite is
**425 tests, 98.5% coverage**, and the cross-product scientific guard now sweeps
19 populations × 47 therapies × 3 modes.

**Non-claims.** No cure is demonstrated, simulated or implied. Animal results
never produce human event credit. Active-comparator and Mendelian-randomisation
contrasts are never applied as placebo-relative treatment effects. Nothing here
is a treatment recommendation or a patient-level prediction.
