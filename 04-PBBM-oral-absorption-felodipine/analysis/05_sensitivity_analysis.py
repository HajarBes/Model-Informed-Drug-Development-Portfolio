#!/usr/bin/env python3
# =============================================================================
# 05_sensitivity_analysis.py — Morris Screening (Global Sensitivity)
# =============================================================================
#
# Screens 15 parameters using Morris method (SALib) with targets AUC and Cmax.
# Same pattern as 01-KRAS-G12C-qsp-model/analysis/vpop_sensitivity.py.
#
# Author: Hajar Besbassi

from __future__ import annotations

import os
import sys
import time
import importlib.util

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from SALib.sample import morris as morris_sample
from SALib.analyze import morris as morris_analyze

# --- Import setup and ACAT model ---
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "setup_00", os.path.join(_SCRIPT_DIR, "00_setup.py"))
_setup = importlib.util.module_from_spec(_spec)
sys.modules["setup_00"] = _setup
_spec.loader.exec_module(_setup)

sys.path.insert(0, _SCRIPT_DIR)
from acat_model_02 import (
    simulate, compute_pk_metrics,
    DrugParams, GISegment, FASTED_GI,
)

DIR_FIGURES = _setup.DIR_FIGURES
DIR_TABLES = _setup.DIR_TABLES
SEED_SENSITIVITY = _setup.SEED_SENSITIVITY
set_pub_style = _setup.set_pub_style
setup_logging = _setup.setup_logging

log = setup_logging("05_sensitivity")

# =============================================================================
# Parameter definitions for screening
# =============================================================================
PARAM_DEFS = {
    "S0_ug_mL":           {"nominal": 0.50,   "low": 0.1,   "high": 2.0,    "label": "Intrinsic\nsolubility (S0)"},
    "pKa":                {"nominal": 5.07,   "low": 4.0,   "high": 6.0,    "label": "pKa"},
    "Peff_cm_s":          {"nominal": 5e-4,   "low": 1e-4,  "high": 1e-3,   "label": "Permeability\n(Peff)"},
    "particle_radius_um": {"nominal": 25.0,   "low": 5.0,   "high": 75.0,   "label": "Particle\nradius"},
    "Fg":                 {"nominal": 0.50,   "low": 0.30,  "high": 0.80,   "label": "Gut avail.\n(Fg)"},
    "Fh":                 {"nominal": 0.50,   "low": 0.30,  "high": 0.80,   "label": "Hepatic avail.\n(Fh)"},
    "CL_L_h":             {"nominal": 70.0,   "low": 35.0,  "high": 140.0,  "label": "Clearance\n(CL)"},
    "Vc_L":               {"nominal": 90.0,   "low": 45.0,  "high": 180.0,  "label": "Central vol.\n(Vc)"},
    "density_g_cm3":      {"nominal": 1.3,    "low": 1.0,   "high": 1.6,    "label": "Particle\ndensity"},
    "Deff_cm2_s":         {"nominal": 5e-6,   "low": 1e-6,  "high": 1e-5,   "label": "Diffusion\ncoeff. (Deff)"},
    "dose_mg":            {"nominal": 10.0,   "low": 2.5,   "high": 20.0,   "label": "Dose"},
}

# GI physiology parameters (modify segment 0 for gastric, segments 2-3 for jejunal)
GI_PARAM_DEFS = {
    "gastric_transit_h":  {"nominal": 0.25,  "low": 0.10,  "high": 2.0,   "label": "Gastric\ntransit"},
    "jejunal_transit_h":  {"nominal": 0.75,  "low": 0.30,  "high": 1.50,  "label": "Jejunal\ntransit"},
    "bile_factor_duo":    {"nominal": 2.0,   "low": 1.0,   "high": 15.0,  "label": "Bile factor\n(duodenum)"},
    "gastric_pH":         {"nominal": 1.7,   "low": 1.0,   "high": 5.0,   "label": "Gastric pH"},
}

ALL_PARAMS = {**PARAM_DEFS, **GI_PARAM_DEFS}
param_names = list(ALL_PARAMS.keys())
bounds = [[ALL_PARAMS[p]["low"], ALL_PARAMS[p]["high"]] for p in param_names]

# =============================================================================
# Morris sampling
# =============================================================================
log.info(f"Morris screening: {len(param_names)} parameters")

problem = {
    "num_vars": len(param_names),
    "names": param_names,
    "bounds": bounds,
}

N_TRAJ = 30
X = morris_sample.sample(problem, N=N_TRAJ, seed=SEED_SENSITIVITY)
n_runs = X.shape[0]
log.info(f"  {n_runs} model evaluations ({N_TRAJ} trajectories x {len(param_names)+1} levels)")

# =============================================================================
# Run simulations
# =============================================================================
Y_auc = np.zeros(n_runs)
Y_cmax = np.zeros(n_runs)

import copy

log.info("  Running simulations...")
t0 = time.time()

for i in range(n_runs):
    # Build drug params from sample
    drug = DrugParams()
    gi = [copy.deepcopy(seg) for seg in FASTED_GI]

    for j, pname in enumerate(param_names):
        val = float(X[i, j])
        if pname in PARAM_DEFS:
            setattr(drug, pname, val)
        elif pname == "gastric_transit_h":
            gi[0].transit_time_h = val
        elif pname == "jejunal_transit_h":
            gi[2].transit_time_h = val
            gi[3].transit_time_h = val
        elif pname == "bile_factor_duo":
            gi[1].bile_factor = val
        elif pname == "gastric_pH":
            gi[0].pH = val

    try:
        sim = simulate(drug, gi, t_end_h=24.0, dt_h=0.05)
        met = compute_pk_metrics(sim)
        Y_auc[i] = met["AUC_ngmL_h"]
        Y_cmax[i] = met["Cmax_ngmL"]
    except Exception:
        Y_auc[i] = 0.0
        Y_cmax[i] = 0.0

    if (i + 1) % 50 == 0:
        elapsed = time.time() - t0
        log.info(f"    {i+1}/{n_runs} done ({elapsed:.1f}s)")

elapsed = time.time() - t0
log.info(f"  All {n_runs} evaluations done in {elapsed:.1f}s")

# =============================================================================
# Morris analysis
# =============================================================================
log.info("Analyzing Morris results...")

results_rows = []
for outcome_name, Y in [("AUC", Y_auc), ("Cmax", Y_cmax)]:
    Y_clean = np.nan_to_num(Y, nan=0.0)
    Si = morris_analyze.analyze(problem, X, Y_clean, seed=SEED_SENSITIVITY)

    for j, pname in enumerate(param_names):
        results_rows.append({
            "outcome": outcome_name,
            "parameter": pname,
            "label": ALL_PARAMS[pname]["label"].replace("\n", " "),
            "mu_star": float(Si["mu_star"][j]),
            "mu": float(Si["mu"][j]),
            "sigma": float(Si["sigma"][j]),
        })

sa_df = pd.DataFrame(results_rows)

# Save separate CSVs for AUC and Cmax
for outcome in ["AUC", "Cmax"]:
    subset = sa_df[sa_df["outcome"] == outcome].sort_values("mu_star", ascending=False)
    subset.to_csv(os.path.join(DIR_TABLES, f"sensitivity_morris_{outcome.lower()}.csv"),
                  index=False)

log.info("  Saved sensitivity_morris_auc.csv and sensitivity_morris_cmax.csv")

# =============================================================================
# Figure 5: Tornado plot
# =============================================================================
set_pub_style()
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

for ax, outcome, title, color in [
    (ax1, "AUC", "AUC Sensitivity", "#377eb8"),
    (ax2, "Cmax", "Cmax Sensitivity", "#e41a1c"),
]:
    sub = sa_df[sa_df["outcome"] == outcome].sort_values("mu_star", ascending=True)
    y_pos = range(len(sub))
    labels = [ALL_PARAMS[p]["label"] for p in sub["parameter"]]

    ax.barh(y_pos, sub["mu_star"], xerr=sub["sigma"],
            color=color, alpha=0.85, edgecolor="white",
            capsize=3, height=0.6)
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Morris mu* (mean |elementary effect|)", fontsize=9)
    ax.set_title(title, fontsize=11, fontweight="bold")

fig.suptitle("Global Sensitivity Analysis — Felodipine Oral Absorption (Morris Screening)",
             fontsize=13, fontweight="bold", y=1.02)
fig.tight_layout()
fig.savefig(os.path.join(DIR_FIGURES, "tornado_sensitivity.png"), dpi=300)
plt.close(fig)
log.info("  Saved tornado_sensitivity.png")

log.info("Sensitivity analysis complete.")
