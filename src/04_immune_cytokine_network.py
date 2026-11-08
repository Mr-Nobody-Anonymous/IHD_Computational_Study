"""
================================================================================
NON-USABLE - LEGACY DEMONSTRATION SCRIPT. KEPT FOR REPRODUCIBILITY ONLY.
================================================================================
This script is OUTSIDE the tested package (`src/ihd_solutions/`) and outside the
test suite. It has no validation, no provenance records and no uncertainty
propagation. Its numbers must NOT be quoted as findings and must NOT be fed into
the validated pipeline.

See `docs/scientific_audit.md` for what each script can and cannot support.
================================================================================

What it can support: a reproduction target. The cascade is exactly linear in its
14 rates, and 0 of 14 are jointly identifiable (see
`ihd_solutions.calibration`). It is used in the scientific audit as a HARD TEST
CASE for the identifiability machinery, not as a model of inflammation.
What it CANNOT support: any calibrated statement about cytokine dynamics. The
rates are assumptions; six of them share a single degree of freedom.
================================================================================
SCRIPT 4: Vascular Immune & Cytokine Network Simulation
================================================================================
Models the 180-day dynamics of the vascular inflammatory signaling cascade:
oxLDL -> TLR4/CD36 -> NF-κB -> NLRP3 Inflammasome -> Caspase-1 -> IL-1β -> IL-6 -> hs-CRP -> MMPs

Simulates:
1. Untreated High Risk
2. High-Intensity Statin
3. Canakinumab (Anti-IL-1β, CANTOS trial)
4. Colchicine (NLRP3 / Tubulin Inhibitor, COLCOT trial)
5. Curative Dual Strategy (PCSK9 Silencer + Inflammasome Blockade)

Outputs:
- results/figures/cytokine_dynamics.png
- results/tables/cytokine_network_data.csv
================================================================================
"""

# ------------------------------------------------------------------------------
# This script is legacy. Shout if someone runs it, so its output is never mistaken
# for a validated result. `python -W ignore src/0N_*.py` silences it deliberately.
# ------------------------------------------------------------------------------
import warnings as _warnings

_warnings.warn(
    "NON-USABLE LEGACY SCRIPT: untested, unvalidated, no provenance. "
    "Numbers are NOT findings. See docs/scientific_audit.md.",
    UserWarning,
    stacklevel=2,
)


import os, sys
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
FIG_DIR  = os.path.join(REPO_DIR, "results", "figures")
TAB_DIR  = os.path.join(REPO_DIR, "results", "tables")
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TAB_DIR, exist_ok=True)

print("=" * 75)
print("  SCRIPT 4: VASCULAR IMMUNE & CYTOKINE NETWORK SIMULATION")
print("=" * 75)

DAYS = 180
dt = 0.25
time_steps = np.arange(0, DAYS + dt, dt)
N = len(time_steps)

INTERVENTIONS = [
    {"name": "Untreated High Risk", "oxldl_factor": 1.0, "nlrp3_inh": 0.0, "il1b_inh": 0.0, "color": "#ef4444", "style": "-"},
    {"name": "High-Intensity Statin", "oxldl_factor": 0.50, "nlrp3_inh": 0.15, "il1b_inh": 0.0, "color": "#f59e0b", "style": "--"},
    {"name": "Canakinumab (Anti-IL-1β, CANTOS)", "oxldl_factor": 1.0, "nlrp3_inh": 0.0, "il1b_inh": 0.85, "color": "#8b5cf6", "style": "-."},
    {"name": "Colchicine (NLRP3 Inhibitor)", "oxldl_factor": 1.0, "nlrp3_inh": 0.75, "il1b_inh": 0.10, "color": "#06b6d4", "style": ":"},
    {"name": "Curative Combo (PCSK9 silencer + NLRP3i)", "oxldl_factor": 0.20, "nlrp3_inh": 0.80, "il1b_inh": 0.50, "color": "#10b981", "style": "-"}
]

k_nfkb_act   = 0.08
k_nfkb_deg   = 0.05
k_nlrp3_syn  = 0.12
k_nlrp3_deg  = 0.08
k_casp1_act  = 0.15
k_casp1_deg  = 0.10
k_il1b_prod  = 0.25
k_il1b_deg   = 0.18
k_il6_prod   = 0.35
k_il6_deg    = 0.20
k_crp_prod   = 0.40
k_crp_deg    = 0.15
k_mmp_prod   = 0.18
k_mmp_deg    = 0.12

simulation_results = {}

for drug in INTERVENTIONS:
    name = drug["name"]
    oxldl = drug["oxldl_factor"]
    nlrp3_block = drug["nlrp3_inh"]
    il1b_block  = drug["il1b_inh"]

    nfkb = np.zeros(N)
    nlrp3 = np.zeros(N)
    casp1 = np.zeros(N)
    il1b = np.zeros(N)
    il6  = np.zeros(N)
    crp  = np.zeros(N)
    mmp  = np.zeros(N)

    nfkb[0]  = 1.0
    nlrp3[0] = 1.0
    casp1[0] = 1.0
    il1b[0]  = 14.5
    il6[0]   = 8.2
    crp[0]   = 3.8
    mmp[0]   = 100.0

    treatment_start_step = int(7.0 / dt)

    for t in range(N - 1):
        if t < treatment_start_step:
            curr_oxldl = 1.0
            curr_nlrp3 = 0.0
            curr_il1b  = 0.0
        else:
            curr_oxldl = oxldl
            curr_nlrp3 = nlrp3_block
            curr_il1b  = il1b_block

        d_nfkb = (k_nfkb_act * curr_oxldl - k_nfkb_deg * nfkb[t]) * dt
        nfkb[t + 1] = max(0.05, nfkb[t] + d_nfkb)

        d_nlrp3 = (k_nlrp3_syn * nfkb[t + 1] * (1.0 - curr_nlrp3) - k_nlrp3_deg * nlrp3[t]) * dt
        nlrp3[t + 1] = max(0.02, nlrp3[t] + d_nlrp3)

        d_casp1 = (k_casp1_act * nlrp3[t + 1] - k_casp1_deg * casp1[t]) * dt
        casp1[t + 1] = max(0.02, casp1[t] + d_casp1)

        d_il1b = (k_il1b_prod * casp1[t + 1] * 14.5 * (1.0 - curr_il1b) - k_il1b_deg * il1b[t]) * dt
        il1b[t + 1] = max(0.5, il1b[t] + d_il1b)

        d_il6 = (k_il6_prod * (il1b[t + 1] / 14.5) * 8.2 - k_il6_deg * il6[t]) * dt
        il6[t + 1] = max(0.5, il6[t] + d_il6)

        d_crp = (k_crp_prod * (il6[t + 1] / 8.2) * 3.8 - k_crp_deg * crp[t]) * dt
        crp[t + 1] = max(0.2, crp[t] + d_crp)

        d_mmp = (k_mmp_prod * 0.5 * (il1b[t + 1]/14.5 + il6[t + 1]/8.2) * 100.0 - k_mmp_deg * mmp[t]) * dt
        mmp[t + 1] = max(5.0, mmp[t] + d_mmp)

    simulation_results[name] = {
        "time": time_steps,
        "il1b": il1b,
        "il6": il6,
        "crp": crp,
        "mmp": mmp,
        "color": drug["color"],
        "style": drug["style"]
    }
    print(f"  ✓ {name:<40} -> Day 180 hs-CRP: {crp[-1]:.2f} mg/L | MMP: {mmp[-1]:.1f}%")

fig, axes = plt.subplots(2, 2, figsize=(14, 10), dpi=300)
plt.subplots_adjust(hspace=0.28, wspace=0.25)

ax = axes[0, 0]
for name, data in simulation_results.items():
    ax.plot(data["time"], data["il1b"], label=name, color=data["color"], linestyle=data["style"], linewidth=2.5)
ax.axvline(7.0, color="gray", linestyle=":", alpha=0.7, label="Initiation (Day 7)")
ax.set_title("A. Active Interleukin-1β (IL-1β)", fontsize=12, fontweight="bold")
ax.set_xlabel("Time (Days)")
ax.set_ylabel("IL-1β (pg/mL)")
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8.5, loc="upper right")

ax = axes[0, 1]
for name, data in simulation_results.items():
    ax.plot(data["time"], data["il6"], label=name, color=data["color"], linestyle=data["style"], linewidth=2.5)
ax.set_title("B. Interleukin-6 (IL-6)", fontsize=12, fontweight="bold")
ax.set_xlabel("Time (Days)")
ax.set_ylabel("IL-6 (pg/mL)")
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8.5, loc="upper right")

ax = axes[1, 0]
for name, data in simulation_results.items():
    ax.plot(data["time"], data["crp"], label=name, color=data["color"], linestyle=data["style"], linewidth=2.5)
ax.axhline(2.0, color="darkgreen", linestyle="--", linewidth=1.5, label="Target: hs-CRP < 2.0 mg/L")
ax.set_title("C. High-Sensitivity C-Reactive Protein (hs-CRP)", fontsize=12, fontweight="bold")
ax.set_xlabel("Time (Days)")
ax.set_ylabel("hs-CRP (mg/L)")
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8.5, loc="upper right")

ax = axes[1, 1]
for name, data in simulation_results.items():
    ax.plot(data["time"], data["mmp"], label=name, color=data["color"], linestyle=data["style"], linewidth=2.5)
ax.set_title("D. Plaque-Degrading MMP Activity", fontsize=12, fontweight="bold")
ax.set_xlabel("Time (Days)")
ax.set_ylabel("MMP Activity (%)")
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8.5, loc="upper right")

plot_path = os.path.join(FIG_DIR, "cytokine_dynamics.png")
plt.savefig(plot_path, bbox_inches="tight")
plt.close()

csv_path = os.path.join(TAB_DIR, "cytokine_network_data.csv")
with open(csv_path, "w", encoding="utf-8") as f:
    f.write("Day,Regimen,IL1b_pg_mL,IL6_pg_mL,hsCRP_mg_L,MMP_Activity_pct\n")
    for name, data in simulation_results.items():
        for d_idx in range(0, N, int(15.0 / dt)):
            d = data["time"][d_idx]
            f.write(f"{d:.1f},\"{name}\",{data['il1b'][d_idx]:.2f},{data['il6'][d_idx]:.2f},"
                    f"{data['crp'][d_idx]:.2f},{data['mmp'][d_idx]:.1f}\n")

print(f"✓ Saved plot to: {plot_path}")
print(f"✓ Saved CSV to:  {csv_path}")
print("✅ Script 4 complete!")
