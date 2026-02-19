#!/usr/bin/env python3
# =============================================================================
# 06_formulation_scenarios.py — Particle Size Impact + Virtual BE
# =============================================================================
#
# Evaluates 3 formulations (micronized, reference, coarser) across a virtual
# population with LHS-sampled physiological variability. Computes geometric
# mean ratios with 90% CI for Cmax and AUC.
#
# Disclaimer: Illustrative virtual BE exercise, not a regulatory submission.
#
# Author: Hajar Besbassi

from __future__ import annotations

import os
import sys
import copy
import time
import importlib.util

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import norm as sp_norm

# --- Import setup and ACAT model ---
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "setup_00", os.path.join(_SCRIPT_DIR, "00_setup.py"))
_setup = importlib.util.module_from_spec(_spec)
sys.modules["setup_00"] = _setup
_spec.loader.exec_module(_setup)

sys.path.insert(0, _SCRIPT_DIR)
from acat_model_02 import (
    simulate, compute_pk_metrics, dissolution_rate, saturation_solubility,
    DrugParams, GISegment, FASTED_GI,
)

DIR_FIGURES = _setup.DIR_FIGURES
DIR_TABLES = _setup.DIR_TABLES
SEED_VPOP = _setup.SEED_VPOP
set_pub_style = _setup.set_pub_style
setup_logging = _setup.setup_logging

log = setup_logging("06_formulation")

# =============================================================================
# Formulation definitions
# =============================================================================
FORMULATIONS = {
    "micronized": {"radius_um": 10.0, "label": "Micronized (10 µm)", "color": "#2166ac"},
    "reference":  {"radius_um": 25.0, "label": "Reference (25 µm)",  "color": "#4daf4a"},
    "coarser":    {"radius_um": 50.0, "label": "Coarser (50 µm)",    "color": "#d62728"},
}

# =============================================================================
# In vitro dissolution profiles
# =============================================================================
log.info("Computing in vitro dissolution profiles...")

drug_base = DrugParams()
diss_rows = []
t_diss = np.arange(0, 120.1, 1.0)  # 0-120 min

for form_name, form in FORMULATIONS.items():
    radius_cm = form["radius_um"] * 1e-4
    # Simulate dissolution in pH 6.8 phosphate buffer (fasted intestinal)
    C_sat = saturation_solubility(drug_base.S0_mg_mL, drug_base.pKa, 6.8, 5.0)
    S_mg = drug_base.dose_mg
    D_mg = 0.0
    V_mL = 900.0  # USP dissolution vessel

    for t_min in t_diss:
        # Record current state BEFORE advancing dissolution
        diss_rows.append({
            "formulation": form_name,
            "time_min": t_min,
            "pct_dissolved": D_mg / drug_base.dose_mg * 100,
        })

        # Advance dissolution (Euler forward, 1 min step)
        C_diss = D_mg / V_mL
        rate = dissolution_rate(S_mg, C_diss, C_sat,
                                drug_base.Deff_cm2_s, drug_base.density_g_cm3,
                                radius_cm, drug_base.z_dissolution)
        dt_h = 1.0 / 60.0  # 1 min in hours
        dissolved = rate * dt_h
        dissolved = min(dissolved, S_mg)
        S_mg -= dissolved
        D_mg += dissolved

diss_df = pd.DataFrame(diss_rows)
diss_df.to_csv(os.path.join(DIR_TABLES, "dissolution_profiles.csv"), index=False)
log.info("  Saved dissolution_profiles.csv")

# =============================================================================
# Virtual population — LHS sampling
# =============================================================================
N_VPOP = 200

log.info(f"Generating virtual population (N={N_VPOP} per formulation)...")

VPOP_PARAMS = {
    "Fg":       {"cv": 0.25},
    "Fh":       {"cv": 0.25},
    "CL_L_h":   {"cv": 0.30},
    "Vc_L":     {"cv": 0.20},
    "Peff_cm_s": {"cv": 0.30},
}

GI_VPOP_PARAMS = {
    "gastric_transit_h": {"nominal": 0.25, "cv": 0.50},
    "bile_factor_scale":  {"nominal": 1.0,  "cv": 0.40},
}


def lhs_sample(n: int, n_params: int, rng: np.random.Generator) -> np.ndarray:
    result = np.zeros((n, n_params))
    for j in range(n_params):
        perm = rng.permutation(n)
        for i in range(n):
            result[perm[i], j] = (i + rng.uniform()) / n
    return result


def generate_vpop(n: int, seed: int):
    """Generate virtual patients with varied physiology."""
    rng = np.random.default_rng(seed)
    drug_params = list(VPOP_PARAMS.keys())
    gi_params = list(GI_VPOP_PARAMS.keys())
    all_names = drug_params + gi_params
    U = lhs_sample(n, len(all_names), rng)

    patients = []
    for i in range(n):
        drug = DrugParams()
        gi = [copy.deepcopy(seg) for seg in FASTED_GI]

        for j, pname in enumerate(drug_params):
            nominal = getattr(drug, pname)
            cv = VPOP_PARAMS[pname]["cv"]
            sigma = np.sqrt(np.log(1 + cv ** 2))
            mu = np.log(nominal) - 0.5 * sigma ** 2
            val = np.exp(mu + sigma * sp_norm.ppf(U[i, j]))
            setattr(drug, pname, float(val))

        # GI parameters
        j_base = len(drug_params)
        # Gastric transit
        nom_gt = GI_VPOP_PARAMS["gastric_transit_h"]["nominal"]
        cv_gt = GI_VPOP_PARAMS["gastric_transit_h"]["cv"]
        sigma_gt = np.sqrt(np.log(1 + cv_gt ** 2))
        mu_gt = np.log(nom_gt) - 0.5 * sigma_gt ** 2
        gt_val = np.exp(mu_gt + sigma_gt * sp_norm.ppf(U[i, j_base]))
        gi[0].transit_time_h = float(gt_val)

        # Bile factor scale
        nom_bf = GI_VPOP_PARAMS["bile_factor_scale"]["nominal"]
        cv_bf = GI_VPOP_PARAMS["bile_factor_scale"]["cv"]
        sigma_bf = np.sqrt(np.log(1 + cv_bf ** 2))
        mu_bf = np.log(nom_bf) - 0.5 * sigma_bf ** 2
        bf_scale = np.exp(mu_bf + sigma_bf * sp_norm.ppf(U[i, j_base + 1]))
        for seg in gi:
            seg.bile_factor *= float(bf_scale)

        patients.append((drug, gi))
    return patients


patients = generate_vpop(N_VPOP, SEED_VPOP)

# =============================================================================
# Simulate each formulation across virtual population
# =============================================================================
log.info("Running virtual BE simulations...")
all_rows = []

for form_name, form in FORMULATIONS.items():
    t0 = time.time()
    log.info(f"  Formulation: {form['label']}")

    for i, (drug, gi) in enumerate(patients):
        drug_mod = copy.deepcopy(drug)
        drug_mod.particle_radius_um = form["radius_um"]

        try:
            sim = simulate(drug_mod, gi, t_end_h=24.0, dt_h=0.05)
            met = compute_pk_metrics(sim)
        except Exception:
            met = {"Cmax_ngmL": np.nan, "AUC_ngmL_h": np.nan,
                   "Tmax_h": np.nan, "Fa": np.nan}

        all_rows.append({
            "patient": i + 1,
            "formulation": form_name,
            "Cmax_ngmL": met["Cmax_ngmL"],
            "AUC_ngmL_h": met["AUC_ngmL_h"],
            "Tmax_h": met["Tmax_h"],
            "Fa": met["Fa"],
        })

        if (i + 1) % 50 == 0:
            log.info(f"    {i+1}/{N_VPOP} patients done")

    elapsed = time.time() - t0
    log.info(f"    Done in {elapsed:.1f}s")

vpop_df = pd.DataFrame(all_rows)

# =============================================================================
# Compute geometric mean ratios (Test vs Reference)
# =============================================================================
log.info("Computing virtual BE statistics...")

ref_df = vpop_df[vpop_df["formulation"] == "reference"].set_index("patient")
be_rows = []

for form_name in ["micronized", "coarser"]:
    test_df = vpop_df[vpop_df["formulation"] == form_name].set_index("patient")

    for metric in ["Cmax_ngmL", "AUC_ngmL_h"]:
        ref_vals = ref_df[metric].values
        test_vals = test_df[metric].values

        # Drop any NaN pairs
        mask = ~(np.isnan(ref_vals) | np.isnan(test_vals)) & (ref_vals > 0) & (test_vals > 0)
        ref_log = np.log(ref_vals[mask])
        test_log = np.log(test_vals[mask])
        diff_log = test_log - ref_log

        gmr = np.exp(np.mean(diff_log))
        se = np.std(diff_log, ddof=1) / np.sqrt(len(diff_log))
        ci_lo = np.exp(np.mean(diff_log) - 1.645 * se)
        ci_hi = np.exp(np.mean(diff_log) + 1.645 * se)

        within_be = 0.80 <= ci_lo and ci_hi <= 1.25
        metric_label = "Cmax" if "Cmax" in metric else "AUC"

        be_rows.append({
            "test_formulation": form_name,
            "metric": metric_label,
            "GMR": round(gmr, 4),
            "CI90_lower": round(ci_lo, 4),
            "CI90_upper": round(ci_hi, 4),
            "within_80_125": within_be,
            "n_subjects": int(mask.sum()),
        })
        log.info(f"  {form_name} vs reference — {metric_label}: "
                 f"GMR={gmr:.3f} ({ci_lo:.3f}-{ci_hi:.3f}) "
                 f"{'PASS' if within_be else 'FAIL'}")

be_df = pd.DataFrame(be_rows)
be_df.to_csv(os.path.join(DIR_TABLES, "virtual_be_summary.csv"), index=False)
log.info("  Saved virtual_be_summary.csv")

# =============================================================================
# Figure 6: Dissolution-Exposure Linkage
# =============================================================================
set_pub_style()
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

# Dissolution profiles
for form_name, form in FORMULATIONS.items():
    sub = diss_df[diss_df["formulation"] == form_name]
    ax1.plot(sub["time_min"], sub["pct_dissolved"], "-",
             color=form["color"], linewidth=2.0, label=form["label"])

ax1.axhline(85, color="gray", linestyle="--", linewidth=0.8, alpha=0.7)
ax1.text(115, 86, "85%", fontsize=8, color="gray")
ax1.set_xlabel("Time (min)")
ax1.set_ylabel("% Dissolved")
ax1.set_title("In Vitro Dissolution (pH 6.8, 900 mL)")
ax1.set_xlim(0, 120)
ax1.set_ylim(0, 105)
ax1.legend(fontsize=8)

# PK profiles (population median)
for form_name, form in FORMULATIONS.items():
    drug_med = DrugParams(particle_radius_um=form["radius_um"])
    sim = simulate(drug_med, FASTED_GI, t_end_h=24.0, dt_h=0.05)
    ax2.plot(sim["time_h"], sim["Cp_ngmL"], "-",
             color=form["color"], linewidth=2.0, label=form["label"])

ax2.set_xlabel("Time (h)")
ax2.set_ylabel("Plasma concentration (ng/mL)")
ax2.set_title("Predicted In Vivo PK (Fasted, 10 mg)")
ax2.set_xlim(0, 24)
ax2.set_ylim(bottom=0)
ax2.legend(fontsize=8)

fig.suptitle("Dissolution-Exposure Linkage — Particle Size Effect on Felodipine",
             fontsize=13, fontweight="bold", y=1.02)
fig.tight_layout()
fig.savefig(os.path.join(DIR_FIGURES, "dissolution_exposure_linkage.png"), dpi=300)
plt.close(fig)
log.info("  Saved dissolution_exposure_linkage.png")

# =============================================================================
# Figure 7: Formulation Comparability (virtual BE forest plot)
# =============================================================================
fig, ax = plt.subplots(figsize=(8, 4))

y_positions = []
y_labels = []
y_idx = 0
colors_be = {"micronized": "#2166ac", "coarser": "#d62728"}

for form_name in ["micronized", "coarser"]:
    for metric in ["Cmax", "AUC"]:
        row = be_df[(be_df["test_formulation"] == form_name) &
                    (be_df["metric"] == metric)]
        if len(row) == 0:
            continue
        row = row.iloc[0]

        gmr = row["GMR"]
        ci_lo = row["CI90_lower"]
        ci_hi = row["CI90_upper"]

        ax.plot(gmr, y_idx, "o", color=colors_be[form_name], markersize=8)
        ax.plot([ci_lo, ci_hi], [y_idx, y_idx], "-",
                color=colors_be[form_name], linewidth=2.0)

        label = f"{FORMULATIONS[form_name]['label']} — {metric}"
        y_labels.append(label)
        y_positions.append(y_idx)
        y_idx += 1

# BE limits
ax.axvspan(0.80, 1.25, alpha=0.1, color="green")
ax.axvline(0.80, color="green", linestyle="--", linewidth=0.8, alpha=0.7)
ax.axvline(1.25, color="green", linestyle="--", linewidth=0.8, alpha=0.7)
ax.axvline(1.0, color="black", linestyle="-", linewidth=0.5)

ax.set_yticks(y_positions)
ax.set_yticklabels(y_labels, fontsize=9)
ax.set_xlabel("Geometric Mean Ratio (90% CI)")
ax.set_title("Virtual Bioequivalence Assessment",
             fontsize=11, fontweight="bold")
ax.set_xlim(0.5, 1.6)
ax.invert_yaxis()

fig.tight_layout()
fig.savefig(os.path.join(DIR_FIGURES, "formulation_comparability.png"), dpi=300)
plt.close(fig)
log.info("  Saved formulation_comparability.png")

log.info("Formulation scenarios and virtual BE complete.")
