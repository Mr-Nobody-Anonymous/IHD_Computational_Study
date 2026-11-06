"""
================================================================================
SCRIPT 3: Coronary Plaque Biomechanics & Vulnerability Dynamics
================================================================================
Models the 60-month progression of human coronary atherogenesis, fibrous cap
thinning, pulsatile tensile hoop stress, and rupture tipping points under:
1. Untreated High Risk (LDL 190 mg/dL)
2. High-Intensity Statin (LDL 95 mg/dL)
3. VERVE-102 Base Editing (LDL 38 mg/dL)
4. Dual Curative Strategy (LDL 38 mg/dL + Inflammasome Blockade)

Outputs:
- results/figures/plaque_simulation_curves.png
- results/tables/plaque_biomechanics_timeline.csv
================================================================================
"""

import os, sys, math
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
print("  SCRIPT 3: CORONARY ARTERY PLAQUE BIOMECHANICS SIMULATION")
print("=" * 75)

MONTHS = 60
dt = 0.5
time_points = np.arange(0, MONTHS + dt, dt)
N = len(time_points)

LUMEN_RADIUS_INIT = 1.5   # mm
BLOOD_PRESSURE_PEAK = 140  # mmHg
P_KPA = BLOOD_PRESSURE_PEAK * 0.133322  # ~18.66 kPa

CAP_RUPTURE_THRESHOLD_UM = 65.0       # µm
STRESS_RUPTURE_THRESHOLD_KPA = 300.0  # kPa

REGIMENS = [
    {"name": "Untreated (High Risk)", "ldl": 190.0, "mmp_activity": 1.0, "color": "#ef4444", "style": "-"},
    {"name": "High-Intensity Statin", "ldl": 95.0, "mmp_activity": 0.70, "color": "#f59e0b", "style": "--"},
    {"name": "VERVE-102 (CRISPR Base Edit)", "ldl": 38.0, "mmp_activity": 0.45, "color": "#3b82f6", "style": "-."},
    {"name": "Dual Cure (VERVE-102 + Anti-IL-1β)", "ldl": 38.0, "mmp_activity": 0.15, "color": "#10b981", "style": "-"}
]

results_by_regimen = {}

for reg in REGIMENS:
    name = reg["name"]
    ldl = reg["ldl"]
    mmp = reg["mmp_activity"]

    cap_thickness = np.zeros(N)
    cap_thickness[0] = 180.0
    necrotic_core = np.zeros(N)
    necrotic_core[0] = 0.40
    stenosis_pct = np.zeros(N)
    stenosis_pct[0] = 32.0
    hoop_stress = np.zeros(N)

    rupture_occurred = False
    rupture_month = None

    for t in range(N - 1):
        influx_rate = 0.008 * (ldl / 100.0)**1.4
        efflux_rate = 0.012 * max(0, (60.0 - ldl) / 60.0)
        necrotic_core[t + 1] = max(0.05, necrotic_core[t] + (influx_rate - efflux_rate) * dt)

        synthesis_rate = 1.2 / (1.0 + (ldl / 120.0)**1.2)
        degradation_rate = 2.4 * mmp * (1.0 + (ldl / 150.0))
        cap_thickness[t + 1] = max(15.0, cap_thickness[t] + (synthesis_rate - degradation_rate) * dt)

        total_plaque_area = necrotic_core[t + 1] * 2.2
        vessel_area = math.pi * (LUMEN_RADIUS_INIT)**2
        stenosis_pct[t + 1] = min(95.0, (total_plaque_area / vessel_area) * 100.0)

        r_lumen = LUMEN_RADIUS_INIT * math.sqrt(max(0.1, 1.0 - stenosis_pct[t + 1] / 100.0))
        t_cap_mm = cap_thickness[t + 1] / 1000.0
        hoop_stress[t + 1] = (P_KPA * r_lumen) / t_cap_mm

        if not rupture_occurred:
            if cap_thickness[t + 1] <= CAP_RUPTURE_THRESHOLD_UM and hoop_stress[t + 1] >= STRESS_RUPTURE_THRESHOLD_KPA:
                rupture_occurred = True
                rupture_month = time_points[t + 1]

    hoop_stress[0] = (P_KPA * LUMEN_RADIUS_INIT) / (cap_thickness[0] / 1000.0)

    results_by_regimen[name] = {
        "time": time_points,
        "cap_thickness": cap_thickness,
        "necrotic_core": necrotic_core,
        "stenosis_pct": stenosis_pct,
        "hoop_stress": hoop_stress,
        "ruptured": rupture_occurred,
        "rupture_month": rupture_month,
        "color": reg["color"],
        "style": reg["style"]
    }

    if rupture_occurred:
        print(f"  ⚠ {name:<35} -> RUPTURE at Month {rupture_month:.1f}")
    else:
        print(f"  ✓ {name:<35} -> STABLE (Min Cap: {cap_thickness.min():.1f} µm)")

# Plotting
fig, axes = plt.subplots(2, 2, figsize=(14, 10), dpi=300)
plt.subplots_adjust(hspace=0.28, wspace=0.25)

ax = axes[0, 0]
for name, data in results_by_regimen.items():
    ax.plot(data["time"], data["cap_thickness"], label=name, color=data["color"], linestyle=data["style"], linewidth=2.5)
ax.axhline(CAP_RUPTURE_THRESHOLD_UM, color="darkred", linestyle=":", linewidth=1.8, label="Vulnerability (<65 µm)")
ax.set_title("A. Fibrous Cap Thickness Remodeling", fontsize=12, fontweight="bold")
ax.set_xlabel("Time (Months)")
ax.set_ylabel("Cap Thickness (µm)")
ax.set_ylim(0, 220)
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8.5, loc="upper right")

ax = axes[0, 1]
for name, data in results_by_regimen.items():
    ax.plot(data["time"], data["hoop_stress"], label=name, color=data["color"], linestyle=data["style"], linewidth=2.5)
ax.axhline(STRESS_RUPTURE_THRESHOLD_KPA, color="darkred", linestyle=":", linewidth=1.8, label="Rupture Ultimate Stress (>300 kPa)")
ax.set_title("B. Peak Tensile Hoop Stress", fontsize=12, fontweight="bold")
ax.set_xlabel("Time (Months)")
ax.set_ylabel("Tensile Stress (kPa)")
ax.set_ylim(50, 450)
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8.5, loc="upper left")

ax = axes[1, 0]
for name, data in results_by_regimen.items():
    ax.plot(data["time"], data["necrotic_core"], label=name, color=data["color"], linestyle=data["style"], linewidth=2.5)
ax.set_title("C. Necrotic Lipid Core Volume", fontsize=12, fontweight="bold")
ax.set_xlabel("Time (Months)")
ax.set_ylabel("Necrotic Core Area (mm²)")
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8.5, loc="upper left")

ax = axes[1, 1]
for name, data in results_by_regimen.items():
    ax.plot(data["time"], data["stenosis_pct"], label=name, color=data["color"], linestyle=data["style"], linewidth=2.5)
ax.axhline(70.0, color="purple", linestyle="--", alpha=0.7, label="Hemodynamic Significance (70%)")
ax.set_title("D. Lumen Stenosis Percentage", fontsize=12, fontweight="bold")
ax.set_xlabel("Time (Months)")
ax.set_ylabel("Stenosis (%)")
ax.set_ylim(20, 85)
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8.5, loc="upper left")

plot_path = os.path.join(FIG_DIR, "plaque_simulation_curves.png")
plt.savefig(plot_path, bbox_inches="tight")
plt.close()

# CSV Output
csv_path = os.path.join(TAB_DIR, "plaque_biomechanics_timeline.csv")
with open(csv_path, "w", encoding="utf-8") as f:
    f.write("Month,Regimen,CapThickness_um,NecroticCore_mm2,Stenosis_pct,TensileStress_kPa\n")
    for name, data in results_by_regimen.items():
        for m_idx in range(0, N, int(6 / dt)):
            m = data["time"][m_idx]
            f.write(f"{m:.1f},\"{name}\",{data['cap_thickness'][m_idx]:.2f},{data['necrotic_core'][m_idx]:.3f},"
                    f"{data['stenosis_pct'][m_idx]:.2f},{data['hoop_stress'][m_idx]:.2f}\n")

print(f"✓ Saved plot to: {plot_path}")
print(f"✓ Saved CSV to:  {csv_path}")
print("✅ Script 3 complete!")
