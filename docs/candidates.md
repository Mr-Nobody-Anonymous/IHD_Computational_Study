# The best candidate for IHD, ranked

**The question:** what is the best candidate for ischaemic heart disease, combining what the
literature says with what this repository's own simulations found?

**The short answer:** deep **LDL-C lowering with agents that are already approved**. It is
not the answer this repository was built to produce. It is what the evidence supports.

This document is generated from `src/ihd_solutions/candidates.py`; run
`ihd-solutions candidates` to reproduce it.

---

## 1. The rule has to come before the answer

A list of candidates is not a ranking. Without a stated criterion, "best" means "most
interesting to the author". This repository already has a criterion, and it is the transfer
gate ([`docs/simulation_frontier.md`](simulation_frontier.md)):

- **item 1** — the endpoint must be a patient outcome, not a surrogate;
- **item 6** — the comparison must be randomised.

So candidates are ranked by the strength of their **outcome** evidence:

| rank | tier | candidates |
|---|---|---|
| 1 | randomised outcome evidence, positive | 2 |
| 2 | randomised outcome evidence, neutral or refuted | 1 |
| 3 | randomised evidence on a surrogate only | 2 |
| 4 | uncontrolled evidence on a surrogate | 1 |
| 5 | mechanism or magnitude claim with no clinical outcome data | 1 |

Two invariants are enforced in code, because both were violated by the first draft of the
module and both change the answer:

1. **A refutation is never counted as a benefit.** The first draft stored ZEUS as
   `effect = 1 − HR = 0.01`, which reads as a 1% win. It now stores the measured hazard
   ratio, 0.99, with `direction="neutral"`.
2. **Magnitudes are never compared across units.** A 3.9 percentage-point LVEF gain and a
   22% relative risk reduction are both "effects"; printing them side by side invites a
   comparison that means nothing. Every candidate carries an `effect_unit`, and the rendered
   string includes it.

---

## 2. The ranking

### 1 — Deep LDL-C lowering to guideline target with approved agents

- **Evidence:** CTT meta-analysis of the statin trials, extended by ezetimibe and the PCSK9
  programmes; ~170,000 patients.
- **Effect:** **22% relative reduction in major vascular events per 1 mmol/L (38.7 mg/dL)
  LDL-C reduction**. This is the most reproduced quantitative relationship in cardiovascular
  prevention.
- **What it licenses:** a quantitative, additive prediction — the same LDL-C reduction buys
  the same proportional risk reduction regardless of how it was achieved.
- **What it does not license:** reversal. It slows a process; it does not remove scar,
  restore myocardium, or eliminate residual risk (roughly two thirds of events in treated
  patients still occur).

### 2 — Low-dose colchicine 0.5 mg daily

- **Evidence:** COLCOT, LoDoCo2, CLEAR SYNERGY; Cochrane 2025 review; ~21,425 patients.
- **Effect:** about **9 fewer MIs and 8 fewer strokes per 1,000 treated**; MI −29%,
  stroke −37%, urgent revascularisation −33%. **No effect on cardiovascular mortality.**
- **Guidelines:** ESC 2024 chronic coronary syndrome **Class IIa, Level A**; ACC/AHA
  Class IIb. A downgrade is possible after the neutral CLEAR SYNERGY result in ACS.
- **What it licenses:** a small, real reduction in non-fatal events, in a cheap generic drug.
- **What it does not license:** a mortality claim — there is none — and it is **not** an
  IL-6 inhibitor, which matters for the next entry.

### 3 — IL-6 ligand inhibition (ziltivekimab) — **refuted**

- **Evidence:** **ZEUS, NCT05021835**, n = **6,385**, reported 31 July 2026, published in
  *NEJM*: ASCVD plus chronic kidney disease plus hsCRP ≥ 2 mg/L.
- **Result:** three-point MACE **HR 0.99 (95% CI 0.88–1.11)** despite **hsCRP −64%** and
  **IL-6 −76%**. Serious infections **20% vs 18%**.
- **Also:** **HERMES and ATHENA were terminated for futility** in September 2026; only
  ARTEMIS (post-AMI) continues.
- **What it licenses:** a refutation, at scale, of the proposition that suppressing IL-6 and
  CRP in high-risk ASCVD reduces events. Target engagement was demonstrated and benefit was
  not.
- **What it does not license:** the claim that inflammation is irrelevant. The population had
  advanced disease and the intervention was downstream of the inflammasome.

**This entry is the one that costs this repository something.** Its inflammatory cascade runs
oxLDL → NF-κB → NLRP3 → IL-1β → IL-6 → CRP, and ZEUS intervened at the last step. The model
would have predicted the 64% CRP fall correctly. It has nothing to say about events, and the
trial now says the CRP fall was not sufficient.

### 4 — Oral PCSK9 inhibitor (enlicitide)

- **Evidence:** CORALreef Lipids, n = 2,909: LDL-C **−56%** placebo-adjusted at week 24,
  sustained ~−48% at week 52, ApoB −50%. **FDA-approved 16 July 2026.**
- **Pending:** **CORALreef Outcomes**, n = 14,500, has enrolled and has not reported.
- **What it licenses:** a large LDL-C reduction delivered orally, without an injection.
- **What it does not license:** an outcome claim of its own — it inherits CTT's class effect
  by argument, not by its own events.

### 5 — Base editing (VERVE-102)

- **Evidence:** Heart-2, NCT06164730, phase 1b, **open-label, n = 35**: LDL-C **−62%** at
  1.0 mg/kg, PCSK9 −88%, durable to 18 months.
- **What it licenses:** proof that a one-time edit can produce a durable LDL-C reduction.
- **What it does not license:** a magnitude claim — 35 uncontrolled patients, and an approved
  oral PCSK9 inhibitor now achieves a comparable reduction **without an irreversible edit**.

### 6 — iPSC-derived cardiomyocyte patch (BioVAT / AATS)

- **Evidence:** BioVAT-HF, n = 20, LVEF **+3.9 pp (90% CI 0.9–6.8)**, uncontrolled; AATS
  late-breaker LB8, n = 8, +4.3 pp, p = 0.156.
- **What it licenses:** a signal worth testing in a randomised trial with an outcome endpoint.
- **What it does not license:** an effect claim. Both estimates sit **below the +6 pp
  surrogate threshold** at which an LVEF change would imply a mortality benefit, and 20
  uncontrolled patients cannot separate an effect from regression to the mean.

### 7 — This repository's combined architecture (cell patch + base editing + OMT)

- **Evidence:** **none.** It has never been administered together in a controlled setting.
- **Own corrected estimate:** 57.56%, with `robust_set_exists=False`.
- **What it licenses:** a hypothesis and an ordered experiment list.
- **What it does not license:** the arithmetic of its parts. The repository measured 15
  pairwise interactions and refused 14.08 pp of double-counted benefit; the wider re-search
  found no robust composition.

**It ranks last.** That is the rule working, not the rule being unkind: a construct with no
clinical data cannot outrank an approved drug with 170,000 patients behind it. Any ranking
that put it first would be measuring novelty, not evidence.

---

## 3. What this does to the repository's own work

The consequence is specific and it is not flattering.

| previously used as | status now |
|---|---|
| `calibration.py`'s CRP target, as a proxy for treatment benefit | **invalid** — ZEUS moved CRP without moving events |
| the same CRP target, as a demonstration of structural non-identifiability | **still valid**, and now better motivated |
| the inflammatory arm as "mechanism of benefit" | demoted to "**marker dynamics**" |

The distinction to keep is that **a marker can be strongly prognostic and still be causally
irrelevant to treatment benefit**. ZEUS separates those two roles for CRP. Nothing in
`calibration.py` needs to change; what needs to change is what any downstream result is
allowed to claim. Any sentence of the form "this lowers CRP, therefore it should help" is now
unsupported.

The correct generalisation is the one the transfer gate already encoded: **a mechanism that
moves a marker is not thereby a mechanism that moves an outcome.**

---

## 4. What the ranking does not say

- It is **not a treatment recommendation** for any individual. Dosing, eligibility and
  interaction are clinical questions and are not answered by a ranking of trial evidence.
- **Nothing here is a cure.** No listed candidate reverses established coronary disease or
  restores infarcted myocardium. The top entry slows progression and reduces events.
- The ranking is **by endpoint strength, not by size of benefit**. A small effect on events
  (colchicine) outranks a large effect on a surrogate (base editing) — deliberately.
- It is **time-stamped**. ZEUS reported in July 2026 and HERMES/ATHENA were stopped in
  September 2026; a positive ARTEMIS or CORALreef Outcomes result would move entries.

---

## 5. Reproducing

```bash
ihd-solutions candidates                  # writes results/candidates/
pytest tests/test_candidates.py           # 29 tests, 100% of the module
```

Artefacts: `candidate_report.txt`, `DISCLAIMER.txt`, and four tables under `tables/`
(ranking, inflammation refutation, CRP surrogate status, best-candidate summary).

Every numeric constant used here is registered in `src/ihd_solutions/provenance.py`; an
unregistered constant fails the build. Registry coverage is 114/114 with 0 values
disagreeing.
