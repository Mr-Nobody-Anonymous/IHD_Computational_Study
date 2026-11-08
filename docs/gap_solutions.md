# Solutions to the six open gaps

> **THESE ARE EXPERIMENTAL DESIGNS AND POWER CALCULATIONS, NOT RESULTS.** Sample sizes
> and detectable differences are arithmetic under stated assumptions; if an assumption
> is wrong the design is wrong. Efficacy is nowhere claimed. External findings belong to
> other people and are tagged by design and species. **No cure is declared.**

Reproduce:

```bash
python run_solutions.py gap-solutions --draws 800
```

Outputs: 10 artefacts in `results/gap_solutions/` (8 tables, 1 report, `DISCLAIMER.txt`).
Tests: `tests/test_gap_solutions.py` (56). Module: `src/ihd_solutions/gap_solutions.py`.

---

## The headline: searching each gap properly overturned half of them

I previously listed six things "nobody has run". Three of those claims were **wrong or
overstated**.

| Gap | Prior claim | Verdict after search |
|---|---|---|
| 1 | Nothing regrows heart muscle in a human | **FALSE** — proven histologically |
| 2 | Response-guided retry never tested | **Stands** |
| 3 | Viability signal needs a reporter never imaged in humans | **Overstated** — imaged in patients already |
| 4 | Repeat dosing awaits CATO | **Stands**, but CATO may be underpowered |
| 5 | Clearance needs 2.18× beyond-precedent potency | **Partly solvable** — reduced to 1.87× |
| 6 | PCSK9-editing durability is assumed | **Now measured** to 18 months |

---

## Gap 1 — remuscularisation has already happened in a human

**I was wrong.** [Jebran et al., *Nature* 2025;639:503–511](https://doi.org/10.1038/s41586-024-08463-0)
examined a human heart explanted three months after engineered-heart-muscle (EHM)
implantation and found **retained, vascularised, viable graft tissue**. That is direct
human histological proof of remuscularisation.

[BioVAT-HF, *NEJM* 2026](https://www.nejm.org/doi/abs/10.1056/NEJMoa2513525) reports the
clinical series: 20 patients, 16 at the maximal dose of 20 EHM units (≈800 million
iPSC-derived cardiomyocytes plus stromal cells).

| Measure | 3 months | Later |
|---|---|---|
| Target wall thickness | **+4.5 mm** | +2.9 mm at 12 months (n=5) |
| LVEF | **+3.9 pp** | **+6.9 pp** at mean 17 months (n=8) |
| KCCQ-OSS | +6.7 points | +15.0 points at 12 months |

No graft-related arrhythmia or tumour. But: open-label, no control, surrogate endpoints,
lifelong immunosuppression, and 196 adverse events (57 serious) in 20 patients.

### So the gap is magnitude, not existence

The repository already carried the 3.9 pp figure. Against the +20 pp criterion:

| Anchor | Required multiplier | EHM units if mass-scaled | Cardiomyocytes |
|---|---|---|---|
| 3-month effect (+3.9 pp) | **5.13×** | 103 | 4.10 billion |
| 17-month effect (+6.9 pp) | **2.90×** | 58 | 2.32 billion |

**Maturation alone buys 1.77× and closes 43.5% of the mass requirement for free.**

### The discordance that changes the answer

Over the same interval, wall thickness **fell 1.6 mm** while LVEF **rose 3.0 pp**.
Function improved as mass declined. That points at **graft maturation and mechanical
entrainment**, not bulk — consistent with the epicardial grafting mechanism the Nature
paper describes.

If that reading is right, mass scaling is the wrong lever and the unit counts above are
an upper bound. **Proposed protocol:** randomise maximal-dose EHM recipients to standard
versus maturation-optimised grafts (longer in-vitro conditioning, paced load training),
active-controlled, followed to 24 months. **Falsifier:** if LVEF gain stops rising after
12 months, maturation is not the lever and the 5.13× mass requirement stands.

---

## Gap 2 — response-guided retry: no new molecule, just an entry rule

Still genuinely untested. But the selection signal is already published: in
Vrtovec/REMEDIUM, patients who gained ≥5 pp from a second dose had **lower LVEF
immediately before it** (34.9% vs 40.0%, p=0.02).

That is available in any clinic today. The design is the existing repeat-dose protocol
with an entry rule at the second dose — re-dose only below a threshold, and randomise
the threshold.

| | Value |
|---|---|
| Baseline responder rate | 10/30 = **33.3%** |
| Best swept enrichment | **2.34×** at LVEF ≤ 36.0 pp |
| Enriched responder rate | **77.9%** |
| Fraction passing the screen | 27.5% |

**Caveat:** the between-patient SD is unreported, so the whole table is swept across
SD 3–5 pp. This is a range, not a promise, and the screening cost is real — roughly
2.6 patients screened out per responder enrolled at the tightest cutoff.

**Falsifier:** if enrichment does not appear prospectively, the Vrtovec subgroup was
chance and fixed-interval dosing is as good.

---

## Gap 3 — the reporter already works in humans; put it in the graft

**Overstated.** Human NIS (sodium-iodide symporter) has already been imaged in patients:

- MV-NIS in relapsed myeloma: **¹²³I scans positive in 8 of 32 patients**
- MV-NIS in ovarian cancer: NIS expression confirmed by ¹²³I SPECT/CT in **3 patients**,
  positive across **repeat cycles** ([Cancer Res 2015](https://cancerres.aacrjournals.org/content/75/1/22))

hNIS is **human-derived and expected non-immunogenic**, and in a direct cardiac
comparison gave **higher signal and contrast than HSV1-tk** (n=34 animals).

| Reporter | Viability-specific | Human imaging | Cardiac graft |
|---|---|---|---|
| **human NIS** | yes | **yes** | **none — this is the gap** |
| HSV1-tk | yes | yes | rodent only; lower contrast |
| Iron-oxide MRI | **no** | yes | signal survives cell death |

**Solution:** the graft is already an engineered allogeneic iPSC line made under GMP, so
adding a stable hNIS cassette is a **manufacturing change, not a new platform**. Image
with ¹²³I SPECT/CT or ¹⁸F-tetrafluoroborate PET. Both components are individually
precedented in humans; the combination in a heart is not.

**Blockers:** thyroid, stomach and salivary glands express NIS natively, so background
needs perchlorate blockade; reporter silencing over time; a second genetic modification
adds regulatory burden.

**Falsifier:** if signal does not track function in humans as it does in rodents
(ρ=0.76–0.80), the feedback loop has nothing to close on.

---

## Gap 4 — the trial meant to answer this may not be able to

This is the most actionable new finding. [CATO](https://clinicaltrials.gov/study/NCT06145035)
randomises 60 patients 1:1:1, so **20 per arm**, with **two** active-versus-control
comparisons (Bonferroni-adjusted here).

| Assumed SD of ΔLVEF | Min. detectable difference | Power for 3 pp | Power for 5 pp | Arm size needed for 3 pp |
|---|---|---|---|---|
| 5 pp | 4.87 pp | 36.5% | 82.1% | 53 |
| **6 pp** | **5.85 pp** | **25.5%** | **65.3%** | **77** |
| 7 pp | 6.82 pp | 18.8% | 50.7% | 104 |
| 8 pp | 7.80 pp | 14.6% | 39.6% | 136 |

**At a plausible SD of 6 pp, CATO has only 25.5% power for a 3 pp increment** — the order
of effect the human repeat-dose literature suggests. It would need 77 per arm (231 total).

**So a null CATO cannot distinguish "no effect" from "an effect too small to see here."**

**Solution:** read it out as planned, but **pre-register that interpretation now**, pool
prospectively with the other repeat-dose trials, and power the successor using the gap-2
enrichment rule rather than repeating an underpowered design.

**Falsifier:** if the true between-patient SD is at the low end, CATO is better powered
than this suggests and the concern is overstated.

---

## Gap 5 — buy the potency requirement down with redundancy

ATLAS-2's ablation found titration worth +72.8 pp, more than any other lever — which
implies the potency requirement is partly substitutable by architecture. Swept rather
than assumed: 6 caps × 4 node counts × 3 gate counts = 72 configurations.

**Cheapest robust point with redundancy inside precedent:** cap **2.40×**, 3 node paths,
3 gate paths → **P = 89.25%**, precedent gap **1.87×** (down from 2.18×).

### An accounting trap I had to fix

My first pass reported cap 2.00× at a 1.56× gap. But that point uses **4 node paths**,
and `MAX_PRECEDENTED_PATHS = 3` — so the redundancy is **itself beyond precedent**.
Paying for a potency breach with an unflagged redundancy breach is not a saving. The
selection now requires both axes to be precedented, which is why the honest answer is
1.87× and not 1.56×.

I also had to stop selecting the *fewest paths* at the cheapest cap: that picked a point
at 82.88%, within 2 SE of the bar — the same "solve exactly on the threshold" failure
this repository already retracted. Selection now takes the most headroom at the cheapest
cap.

**Falsifier:** if no combination below 2.80× is robust with precedented redundancy,
potency is not substitutable and the 2.18× gap stands.

---

## Gap 6 — durability is measured; now power the outcomes trial

[VERVE-102 Heart-2, *NEJM* 2026](https://www.nejm.org/doi/full/10.1056/NEJMoa2601283):
at 1.0 mg/kg, **−88% PCSK9**, **−62% LDL-C**, **−78 mg/dL absolute** (128 → 51 mg/dL),
durable to **18 months** with 15 participants beyond a year. Since hepatocytes turn over
every 200–300 days, persistence through that turnover is real evidence the edit is
stable. No dose-limiting toxicity.

**Durability is no longer an assumption.** What is missing is outcomes — and that trial
can now be powered on a measured effect.

78 mg/dL = **2.02 mmol/L**. Applying the CTT meta-analytic slope log-linearly:

| Assumption | Assumed HR | Events required | Patients at 3%/yr over 3 yr |
|---|---|---|---|
| Full CTT slope | 0.606 | **125** | ~1,389 |
| Three-quarters slope | 0.687 | 223 | ~2,478 |
| Half slope | 0.778 | 500 | ~5,556 |

**On the full slope this is a ~125-event trial, not a mega-trial.** The half-slope arm is
there so the design survives being wrong.

**This is design arithmetic, not a claim.** The CTT slope is imported from statins;
whether editing translates to events on that slope is *exactly what the trial exists to
test*. Asserting the benefit here would be the surrogate-to-outcome leap this study has
refused throughout.

---

## What is actually left

| Gap | Remaining requirement |
|---|---|
| 1 | ~**2.9×** more graft effect than demonstrated; maturation may supply part free |
| 2 | An entry rule nobody has tried, using a signal already in the clinic |
| 3 | A reporter cassette already imaged in patients, never put in a heart |
| 4 | A readout that may be underpowered, plus a successor design if so |
| 5 | A residual **1.87×** potency gap, down from 2.18× |
| 6 | An outcomes trial needing on the order of **125 events** |

Three of these need no new biology at all: gap 2 is an entry rule, gap 3 is a
manufacturing change, gap 6 is a trial design. That is a materially more tractable
picture than "nothing regrows heart muscle."

**None of this is a cure.** Every number above is design arithmetic under stated
assumptions, and every external finding belongs to someone else.

---

## See also

- [`atlas4_precedent_search.md`](atlas4_precedent_search.md) — the precedent registry and exhaustive design grid
- [`atlas3_graft.md`](atlas3_graft.md) — the superseded best-of-k graft model
- [`FULL_REPORT.md`](FULL_REPORT.md) — complete inventory and honesty ledger
- [`scientific_methodology_and_limitations.md`](scientific_methodology_and_limitations.md)
