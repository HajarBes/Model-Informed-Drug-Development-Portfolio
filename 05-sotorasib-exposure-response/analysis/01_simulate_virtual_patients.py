#!/usr/bin/env python3
"""
01_simulate_virtual_patients.py — Virtual Patient PK Simulation
================================================================
Simulates N=500 virtual patients per dose level (180, 240, 360, 720, 960 mg).

Approach: Direct exposure sampling from log-normal distributions anchored to
published geometric mean AUC and Cmax at each dose level, with correlation
structure from shared CL variability. This is the standard E-R approach when
working from published popPK summary statistics rather than running a full
popPK model.

Saturable bioavailability is captured empirically: the published data shows
that GM AUC barely increases from 180 mg to 960 mg (saturable absorption).
The virtual population reproduces this flat dose-exposure pattern.

Outputs:
  - outputs/tables/virtual_patients.csv          (N=2500 rows)
  - outputs/tables/pk_summary_by_dose.csv
  - outputs/tables/pk_anchor_vs_simulated.csv    (calibration table)
  - figures/dose_exposure_boxplot.png
  - figures/saturable_bioavailability.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ── Import project setup ─────────────────────────────────────────────
sys.path.insert(0, str(Path(__file__).resolve().parent))
import importlib.util as _ilu

_setup_path = Path(__file__).resolve().parent / "00_setup.py"
_spec = _ilu.spec_from_file_location("setup_00", _setup_path)
_mod = _ilu.module_from_spec(_spec)
sys.modules["setup_00"] = _mod
_spec.loader.exec_module(_mod)

PROJECT_ROOT = _mod.PROJECT_ROOT
DIR_FIGURES = _mod.DIR_FIGURES
DIR_TABLES = _mod.DIR_TABLES
PK = _mod.PK
CLINICAL = _mod.CLINICAL
DOSE_LEVELS = _mod.DOSE_LEVELS
DOSE_COLORS = _mod.DOSE_COLORS
COLORS = _mod.COLORS
SEED_VPOP = _mod.SEED_VPOP
MW = _mod.MW
get_logger = _mod.get_logger
anchor_vs_simulated_table = _mod.anchor_vs_simulated_table
apply_publication_theme = _mod.apply_publication_theme

apply_publication_theme()
log = get_logger("01_vpop")

# ── Configuration ────────────────────────────────────────────────────
N_PATIENTS = 500        # per dose level
RNG = np.random.default_rng(SEED_VPOP)

# ── Published PK Anchors ─────────────────────────────────────────────
# Day 8 AUC/Cmax are the only cross-dose comparison (FDA review).
# 960 mg steady-state data available from FDA label.
# We scale Day 8 values → approximate SS using the 960 mg ratio.
DAY8_AUC = {180: 31.7, 360: 38.9, 720: 42.1, 960: 32.4}
DAY8_CMAX = {180: 6.44, 360: 6.31, 720: 5.45, 960: 5.39}
SS_AUC_960 = 65.3    # hr·µg/mL (FDA label)
SS_CMAX_960 = 7.50   # µg/mL (FDA label)

_auc_scale = SS_AUC_960 / DAY8_AUC[960]     # ≈ 2.015
_cmax_scale = SS_CMAX_960 / DAY8_CMAX[960]   # ≈ 1.391

# Target SS geometric mean AUC and Cmax at each dose
TARGET_AUC_SS = {}
TARGET_CMAX_SS = {}
for dose in [180, 360, 720, 960]:
    TARGET_AUC_SS[dose] = DAY8_AUC[dose] * _auc_scale
    TARGET_CMAX_SS[dose] = DAY8_CMAX[dose] * _cmax_scale
# 240 mg: interpolate (no Day 8 data at 240 mg)
TARGET_AUC_SS[240] = (TARGET_AUC_SS[180] + TARGET_AUC_SS[360]) / 2.0
TARGET_CMAX_SS[240] = (TARGET_CMAX_SS[180] + TARGET_CMAX_SS[360]) / 2.0

log.info("Target SS exposure by dose:")
for d in DOSE_LEVELS:
    log.info(f"  {d} mg: GM AUC = {TARGET_AUC_SS[d]:.1f} hr·µg/mL, "
             f"GM Cmax = {TARGET_CMAX_SS[d]:.2f} µg/mL")

# ── IIV Parameters ───────────────────────────────────────────────────
# Published CV%: CL/F = 76%, V/F = 135%
# For exposure metrics (AUC ∝ 1/CL, Cmax ∝ 1/V):
#   omega_auc ≈ omega_cl = sqrt(log(1 + 0.76²)) = 0.67
#   omega_cmax depends on both CL and V → use moderate CV ~60%
OMEGA_AUC = np.sqrt(np.log(1 + PK.cl_f_cv**2))   # ≈ 0.67
OMEGA_CMAX = 0.55   # Moderate IIV for Cmax (AUC and Cmax are correlated)
RHO_AUC_CMAX = 0.65  # Correlation: shared CL variability drives both

log.info(f"IIV: omega_AUC = {OMEGA_AUC:.3f}, omega_Cmax = {OMEGA_CMAX:.3f}, "
         f"rho = {RHO_AUC_CMAX:.2f}")


# ── Simulate Virtual Patients ────────────────────────────────────────
log.info(f"Simulating {N_PATIENTS} virtual patients × {len(DOSE_LEVELS)} doses...")

# Covariance matrix for bivariate log-normal (AUC, Cmax)
cov_matrix = np.array([
    [OMEGA_AUC**2, RHO_AUC_CMAX * OMEGA_AUC * OMEGA_CMAX],
    [RHO_AUC_CMAX * OMEGA_AUC * OMEGA_CMAX, OMEGA_CMAX**2],
])

records = []
for dose in DOSE_LEVELS:
    gm_auc = TARGET_AUC_SS[dose]
    gm_cmax = TARGET_CMAX_SS[dose]

    # Draw correlated log-normal deviates
    etas = RNG.multivariate_normal([0, 0], cov_matrix, size=N_PATIENTS)

    for i in range(N_PATIENTS):
        auc_i = gm_auc * np.exp(etas[i, 0])
        cmax_i = gm_cmax * np.exp(etas[i, 1])

        # Ctrough: approximate from AUC and t½
        # Ctrough ≈ AUC * ke * exp(-ke*tau) with ke from individual AUC
        # Use: ke_eff = dose_absorbed / (AUC * V_eff)
        # Simpler: Ctrough ≈ AUC/24 * exp(-0.693/5 * 18) → low for t½=5hr
        ke_approx = 0.693 / PK.t_half
        ctrough_i = cmax_i * np.exp(-ke_approx * 23.0)  # ~23 hr post-Tmax

        # Covariates (demographics)
        wt = max(40.0, RNG.normal(80, 15))
        ecog = int(RNG.choice([0, 1, 2], p=[0.35, 0.50, 0.15]))
        prior_cpi = int(RNG.choice([0, 1], p=[0.60, 0.40]))
        albumin = max(2.0, RNG.normal(3.8, 0.5))

        records.append({
            "patient_id": f"D{dose}_P{i+1:04d}",
            "dose_mg": dose,
            "auc_ss": round(auc_i, 3),
            "cmax_ss": round(cmax_i, 3),
            "ctrough_ss": round(max(ctrough_i, 0.001), 5),
            "weight_kg": round(wt, 1),
            "ecog": ecog,
            "prior_cpi": prior_cpi,
            "albumin": round(albumin, 2),
        })

df = pd.DataFrame(records)
log.info(f"Generated {len(df)} virtual patients")

# ── Summary Statistics ───────────────────────────────────────────────
summary_rows = []
for dose in DOSE_LEVELS:
    sub = df[df["dose_mg"] == dose]
    gm_auc = np.exp(np.mean(np.log(sub["auc_ss"])))
    gm_cmax = np.exp(np.mean(np.log(sub["cmax_ss"])))
    cv_auc = np.sqrt(np.exp(np.var(np.log(sub["auc_ss"]))) - 1) * 100
    summary_rows.append({
        "dose_mg": dose,
        "n": len(sub),
        "gm_auc_ss": round(gm_auc, 2),
        "gm_cmax_ss": round(gm_cmax, 2),
        "median_auc_ss": round(sub["auc_ss"].median(), 2),
        "cv_auc_pct": round(cv_auc, 1),
    })

df_summary = pd.DataFrame(summary_rows)
log.info("\nPK Summary by Dose (simulated virtual patients):")
log.info("\n" + df_summary.to_string(index=False))

# ── Anchor vs Simulated Table (PK) ──────────────────────────────────
sim_auc = {row["dose_mg"]: row["gm_auc_ss"] for _, row in df_summary.iterrows()}
sim_cmax = {row["dose_mg"]: row["gm_cmax_ss"] for _, row in df_summary.iterrows()}

anchor_vs_simulated_table(
    label="PK Calibration — AUCss",
    dose_levels=DOSE_LEVELS,
    published_values=TARGET_AUC_SS,
    simulated_values=sim_auc,
    metric_name="GM AUCss",
    units="hr·µg/mL",
    logger=log,
)

anchor_vs_simulated_table(
    label="PK Calibration — Cmax,ss",
    dose_levels=DOSE_LEVELS,
    published_values=TARGET_CMAX_SS,
    simulated_values=sim_cmax,
    metric_name="GM Cmax,ss",
    units="µg/mL",
    logger=log,
)

# Save anchor table as CSV
anchor_rows = []
for dose in DOSE_LEVELS:
    pub_auc = TARGET_AUC_SS[dose]
    pub_cmax = TARGET_CMAX_SS[dose]
    anchor_rows.append({
        "dose_mg": dose,
        "published_auc": round(pub_auc, 2),
        "simulated_auc": sim_auc.get(dose, np.nan),
        "auc_ratio": round(sim_auc.get(dose, 0) / pub_auc, 3) if pub_auc > 0 else np.nan,
        "published_cmax": round(pub_cmax, 2),
        "simulated_cmax": sim_cmax.get(dose, np.nan),
        "cmax_ratio": round(sim_cmax.get(dose, 0) / pub_cmax, 3) if pub_cmax > 0 else np.nan,
    })
pd.DataFrame(anchor_rows).to_csv(DIR_TABLES / "pk_anchor_vs_simulated.csv", index=False)

# ── Save Outputs ─────────────────────────────────────────────────────
df.to_csv(DIR_TABLES / "virtual_patients.csv", index=False)
df_summary.to_csv(DIR_TABLES / "pk_summary_by_dose.csv", index=False)
log.info(f"Saved virtual_patients.csv ({len(df)} rows)")
log.info("Saved pk_summary_by_dose.csv")
log.info("Saved pk_anchor_vs_simulated.csv")

# ── Figure 1: Dose-Exposure Boxplot ──────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# AUC boxplot
bp_data_auc = [df[df["dose_mg"] == d]["auc_ss"].values for d in DOSE_LEVELS]
bp1 = axes[0].boxplot(
    bp_data_auc, positions=range(len(DOSE_LEVELS)), widths=0.6,
    patch_artist=True, showfliers=False, medianprops=dict(color="black", linewidth=1.5),
)
for patch, dose in zip(bp1["boxes"], DOSE_LEVELS):
    patch.set_facecolor(DOSE_COLORS[dose])
    patch.set_alpha(0.7)
# Overlay published anchors
for i, dose in enumerate(DOSE_LEVELS):
    axes[0].scatter(i, TARGET_AUC_SS[dose], marker="D", color="red", s=80, zorder=5,
                    edgecolors="black", linewidths=0.8,
                    label="Published anchor" if i == 0 else "")
axes[0].set_xticks(range(len(DOSE_LEVELS)))
axes[0].set_xticklabels([str(d) for d in DOSE_LEVELS])
axes[0].set_xlabel("Dose (mg)")
axes[0].set_ylabel("AUCss (hr·µg/mL)")
axes[0].set_title("Steady-State AUC by Dose\n(model-implied, virtual patients, N=500/dose)")
axes[0].legend(loc="upper right", framealpha=0.9)

# Cmax boxplot
bp_data_cmax = [df[df["dose_mg"] == d]["cmax_ss"].values for d in DOSE_LEVELS]
bp2 = axes[1].boxplot(
    bp_data_cmax, positions=range(len(DOSE_LEVELS)), widths=0.6,
    patch_artist=True, showfliers=False, medianprops=dict(color="black", linewidth=1.5),
)
for patch, dose in zip(bp2["boxes"], DOSE_LEVELS):
    patch.set_facecolor(DOSE_COLORS[dose])
    patch.set_alpha(0.7)
for i, dose in enumerate(DOSE_LEVELS):
    axes[1].scatter(i, TARGET_CMAX_SS[dose], marker="D", color="red", s=80, zorder=5,
                    edgecolors="black", linewidths=0.8,
                    label="Published anchor" if i == 0 else "")
axes[1].set_xticks(range(len(DOSE_LEVELS)))
axes[1].set_xticklabels([str(d) for d in DOSE_LEVELS])
axes[1].set_xlabel("Dose (mg)")
axes[1].set_ylabel("Cmax,ss (µg/mL)")
axes[1].set_title("Steady-State Cmax by Dose\n(model-implied, virtual patients, N=500/dose)")
axes[1].legend(loc="upper right", framealpha=0.9)

fig.suptitle("Sotorasib: Less-Than-Dose-Proportional Exposure (Saturable Absorption)",
             fontsize=14, y=1.02)
plt.tight_layout()
fig.savefig(DIR_FIGURES / "dose_exposure_boxplot.png")
plt.close(fig)
log.info("Saved dose_exposure_boxplot.png")

# ── Figure 2: Saturable Bioavailability Curve ────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Panel A: Relative bioavailability vs dose (conceptual)
doses_cont = np.linspace(10, 1200, 500)
# Empirical relative F from targets: F_rel ∝ AUC/dose
f_ref = TARGET_AUC_SS[180] / 180  # reference at lowest dose
f_rel_empirical = {d: TARGET_AUC_SS[d] / d / f_ref for d in DOSE_LEVELS}

# Fit D50 to these empirical F_rel values
from scipy.optimize import minimize_scalar

def _f_rel_model(dose, d50):
    return d50 / (d50 + dose)

def _obj_frel(d50):
    err = 0
    f0 = _f_rel_model(180, d50)
    for dose in DOSE_LEVELS:
        pred_ratio = _f_rel_model(dose, d50) / f0
        target_ratio = f_rel_empirical[dose]
        err += (pred_ratio - target_ratio)**2
    return err

res = minimize_scalar(_obj_frel, bounds=(5, 500), method="bounded")
D50_PLOT = res.x

f_curve = [_f_rel_model(d, D50_PLOT) / _f_rel_model(180, D50_PLOT) for d in doses_cont]
axes[0].plot(doses_cont, f_curve, color=COLORS["primary"], linewidth=2)
for dose in DOSE_LEVELS:
    axes[0].scatter(dose, f_rel_empirical[dose], color=DOSE_COLORS[dose], s=100,
                    zorder=5, edgecolors="black", linewidths=0.8)
    axes[0].annotate(f"{dose} mg", (dose, f_rel_empirical[dose]),
                     textcoords="offset points", xytext=(8, 5), fontsize=9)
axes[0].set_xlabel("Dose (mg)")
axes[0].set_ylabel("Relative Bioavailability (F/F_180)")
axes[0].set_title("Saturable Absorption\nF(dose) Decreases with Increasing Dose")
axes[0].set_ylim(0, 1.15)

# Panel B: AUC vs dose — flat plateau
auc_cont = [TARGET_AUC_SS[180] * (_f_rel_model(d, D50_PLOT) / _f_rel_model(180, D50_PLOT)) * (d / 180)
            for d in doses_cont]
axes[1].plot(doses_cont, auc_cont, color=COLORS["primary"], linewidth=2,
             label="Model (saturable F)", alpha=0.7)
# Linear dose-proportional reference
auc_linear = [TARGET_AUC_SS[180] * d / 180 for d in doses_cont]
axes[1].plot(doses_cont, auc_linear, "--", color=COLORS["light_gray"], linewidth=1.5,
             label="Dose-proportional (reference)", alpha=0.8)
for dose in DOSE_LEVELS:
    axes[1].scatter(dose, TARGET_AUC_SS[dose], marker="D", color="red", s=80,
                    zorder=5, edgecolors="black", linewidths=0.8,
                    label="Published anchor" if dose == 180 else "")
    axes[1].scatter(dose, sim_auc[dose], marker="o", color=DOSE_COLORS[dose], s=60,
                    zorder=4, edgecolors="black", linewidths=0.5,
                    label="Simulated GM" if dose == 180 else "")
axes[1].set_xlabel("Dose (mg)")
axes[1].set_ylabel("AUCss (hr·µg/mL)")
axes[1].set_title("Dose-Exposure Curve: AUCss Plateaus\nDue to Saturable Absorption")
axes[1].legend(loc="upper left", framealpha=0.9, fontsize=9)

fig.suptitle("Sotorasib: Mechanism of Flat Dose-Exposure Relationship",
             fontsize=14, y=1.02)
plt.tight_layout()
fig.savefig(DIR_FIGURES / "saturable_bioavailability.png")
plt.close(fig)
log.info("Saved saturable_bioavailability.png")

log.info("Script 01 complete.")
