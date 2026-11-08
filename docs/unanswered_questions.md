# What this study cannot answer

**A companion to `scientific_audit.md`.** That document asks *what is missing from the
code*. This one asks the harder follow-up: **of the questions the study needs answered,
which are unanswerable without a laboratory — and what does leaving each one open cost?**

```bash
.venv/bin/python run_solutions.py unanswered     # -> results/unanswered/ (6 artefacts)
```

**13 gaps registered. 6 can be closed with arithmetic on this machine. 7 cannot be
closed by any computation at all.** No result in this repository is improved by running
it. Every gap below remains open afterwards.

---

## The one-paragraph answer

The laboratory is not the only thing missing. Of the seven gaps that no computation can
close, only one needs a wet lab; two need a GPU and a dataset; and **four are pure
arithmetic that nobody had done yet**. One of those four is a defect in this study's own
viability table. If the question is *"what should be built next?"*, the answer is not a
simulation — it is **a held-out dataset**, because that is the only item on the list that
would move the repository up its own evidence ladder, and it needs no GPU and no patients.

---

## Four kinds of gap, because they fail differently

| class | meaning | count |
|---|---|---|
| **in principle** | no computation settles it; only a randomised human comparison can | 1 |
| **in practice** | computable, but needs a resource this environment lacks | 3 |
| **analytic** | closable here with arithmetic, **and not yet done** — real omissions | 6 |
| **epistemic** | a number cannot settle it at all; a judgement about what numbers mean | 3 |

The distinction matters because the critique treats all six of its problems as engineering
tasks with a compute budget. Three of them are not.

---

## The gap that decides everything: no outcome data

```
target engagement  →  mechanism  →  surrogate  →  outcome
     MEASURED         INFERRED      MEASURED     NOT MEASURED
   LDL −78 mg/dL    base editing   LVEF +3.9 pp   (uncontrolled)
```

Simulating the links on **either side** of a missing link more precisely does not supply
the missing link. It supplies a more precise version of the assumption that the link
exists. That is the honest boundary of the whole project, and it is not a GPU shortage.

The surrogate cannot stand in either: **STICH** found LVEF change unrelated to death, and
the surrogate meta-analysis gives R² = 0.59 against mortality — **41% of outcome variance
is outside the model by construction.**

## The four arithmetic gaps — these are mine, not the environment's

**1. Neither patch trial could have detected its own result.** With the registered SD of
5 pp:

| study | n | observed | minimum detectable | |
|---|---|---|---|---|
| BioVAT-HF | 20 | +3.9 pp | **4.43 pp** | below its own floor |
| AATS LB8 | 8 | +4.3 pp | **7.00 pp** | below its own floor |
| REVIVED OMT arm | 700 | +3.4 pp | 0.75 pp | properly powered |

Both uncontrolled trials reported effects **smaller than the smallest difference they
could have found.** They cannot separate a real benefit from no benefit, and no re-analysis
changes that. What *would* resolve it:

| effect to resolve | patients per arm | total |
|---|---|---|
| 3.9 pp | 26 | 52 |
| 3.4 pp | 34 | 68 |
| 2.0 pp | 99 | 198 |
| 1.0 pp | 393 | 786 |

This is the most useful number the repository can produce, and it took until now to state
it. It tells anyone planning the next trial exactly what it must be.

**2. Selection is not effect.** The repo refuses to report a stratified gain (population
gain 0%) but never had the mechanism. Simulated: selection on a predictor of outcome
re-orders patients and does not move any of them. In a correctly controlled comparison —
selected-treated against selected-untreated — the effect returns to the population
average. The single-arm readout instead inflates it by **7.6×**.

**3. Composite adverse-event risk has never been reported.** ATLAS-5 stacks three
exposures. At 5.5% per exposure, `1 − (1 − r)ⁿ` = **15.6%** — not the 5.5% quoted. A
reviewer does this in their head; it should be printed before they do.

**4. The winner of a 400-candidate grid is upward biased.** Selection bias inflates the
winner by **≈4.6 pp** at measurement SD 2 pp and true spread 8 pp — larger than the effect
being sought. The held-out split (`SPLIT_RETENTION = 0.9220`) already pays for this; the
magnitude was never stated.

## An open defect in this study's own table

`revival.TRANSMURALITY_RECOVERY`:

| scar band | recovery |
|---|---|
| no LGE | 0.91 |
| 1–50% | 0.44 |
| 51–75% | 0.43 |
| >75% | 0.08 |

The table is monotone, so a casual check passes. The defect is the **size of the steps**.
Burden crosses the mid-wall boundary — where the viability literature says function is
largely lost — and recovery moves **0.44 → 0.43**. A segment with 1% scar and one with 75%
are assigned almost identical recovery, while zero scar scores 0.91. The mapping is a cliff
between 0% and 1%, then flat across the entire range that matters.

**It is not fixed here.** Repair needs the segment counts behind the pooled estimate, and
inventing a corrected curve is the practice this repository forbids. The direction of
repair is stated (LOWER the 51–75% band); the new values are not. This table **flatters the
intervention** wherever the eligible population is described as having "recoverable"
substrate.

## What a GPU would and would not buy

| gap | class | what it needs |
|---|---|---|
| genome-wide off-target | in practice | reference genome + editor-specific model. The legacy screen covers **~0.4%** of sites reachable at 3 mismatches, and two hand-written lists of seven agree on only **four** |
| PCSK9 binding energetics | in practice | multi-replica GPU MD with free-energy convergence |
| 3D plaque mechanics | in practice | CFD/FSI solver + meshed coronary geometry |

A GPU closes these three. It closes **none** of the seven that no computation can close.

## What would actually move the needle

**External validation.** It needs no GPU, no imaging, no patients — only a published
dataset not used to fit any registered constant. It is the only item on the list that moves
the repository from `convergence` to `external`, and it is `NOT DONE`. Splitting data
already used to fit the constants would be internal validation wearing the wrong label, so
it must be a genuinely new table.

What it would still **not** do: external validation of a surrogate does not make the
surrogate an outcome.

## The prior no dataset can supply

"Should the programme continue?" needs a prior over counterfactual outcomes, and there is
no registry of those. Swept across 64 prior combinations with the cost of continuing
included, **the recommendation flips inside a plausible range**. That is the finding: the
decision is prior-dominated, so the honest output is the range with its conditions, never a
single recommendation.

Without a cost term, the decision degenerates to "more evidence is never bad" and cannot
flip — which would make the sweep decorative rather than informative.

## The uncomfortable one

28 patients have been implanted across the two patch trials. One reports 20 of 20 with an
adverse event; the other 0 of 8. Neither is an incidence estimate, because **single-arm
trials do not produce comparable denominators**. The number the repository should be able
to give a patient does not exist here, and pooling would be a statistic about the trials,
not about the therapy.

---

## Code map

| symbol | purpose |
|---|---|
| `gaps.gap_register()` | 13 gaps, classified, with costs where a cost exists |
| `gaps.gap_summary()` | counts by class and the honest headline |
| `gaps.transmurality_reviewer_check()` | the open defect, detected not papered over |
| `gaps.trial_resolution_floor()` | minimum detectable difference for every published trial |
| `gaps.required_trial_size()` | patients per arm for a given target effect |
| `gaps.response_guided_selection_bias()` | why selection cannot produce a population gain |
| `gaps.composite_adverse_event_probability()` | 1 − (1 − r)ⁿ across stacked exposures |
| `gaps.crispr_offtarget_enumerable_count()` | how large a real screen would be |
| `gaps.pick_best_of_k_measurement()` | the winner's curse, quantified |
| `gaps.external_benchmark_gap()` | the cheapest genuine upgrade, and why it is not done |
| `gaps.prior_sensitivity()` | the prior-dominated decision, swept |
| `gaps.write_gap_report()` | 6 artefacts |

Tests: `tests/test_gaps.py` — **52 tests**, `gaps.py` at **100%** line coverage.
