#!/usr/bin/env python3
# =============================================================================
# 04_simulate_fed.py — Fed-State PK Simulation + Food Effect Analysis
# =============================================================================
#
# Simulates fed-state PK for 10 mg IR felodipine using modified GI physiology,
# compares with fasted-state prediction, and decomposes the food effect by
# toggling individual mechanistic knobs.
#
# Author: Hajar Besbassi

from __future__ import annotations

import os
import sys
import copy
import importlib.util

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

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
    DrugParams, GISegment, FASTED_GI, FED_GI,
)

DIR_FIGURES = _setup.DIR_FIGURES
DIR_TABLES = _setup.DIR_TABLES
DIR_OBSERVED_EXTRACTED = _setup.DIR_OBSERVED_EXTRACTED
set_pub_style = _setup.set_pub_style
setup_logging = _setup.setup_logging

log = setup_logging("04_simulate_fed")

# =============================================================================
# Simulate both states
# =============================================================================
drug = DrugParams()

log.info("Simulating fasted-state PK...")
sim_fasted = simulate(drug, FASTED_GI, t_end_h=24.0, dt_h=0.01)
met_fasted = compute_pk_metrics(sim_fasted)

log.info("Simulating fed-state PK...")
sim_fed = simulate(drug, FED_GI, t_end_h=24.0, dt_h=0.01)
met_fed = compute_pk_metrics(sim_fed)

log.info(f"  Fasted: Cmax={met_fasted['Cmax_ngmL']:.2f}, Tmax={met_fasted['Tmax_h']:.2f}, AUC={met_fasted['AUC_ngmL_h']:.2f}")
log.info(f"  Fed:    Cmax={met_fed['Cmax_ngmL']:.2f}, Tmax={met_fed['Tmax_h']:.2f}, AUC={met_fed['AUC_ngmL_h']:.2f}")

# Food effect ratios
cmax_ratio = met_fed["Cmax_ngmL"] / met_fasted["Cmax_ngmL"]
auc_ratio = met_fed["AUC_ngmL_h"] / met_fasted["AUC_ngmL_h"]
tmax_shift = met_fed["Tmax_h"] - met_fasted["Tmax_h"]

log.info(f"\n  Food Effect (IR formulation):")
log.info(f"    Cmax ratio (fed/fasted) = {cmax_ratio:.2f}")
log.info(f"    AUC ratio (fed/fasted)  = {auc_ratio:.2f}")
log.info(f"    Tmax shift              = +{tmax_shift:.1f} h")
log.info(f"  Fed gastric pH (3.0) and emptying time (0.65 h) calibrated to")
log.info(f"  Bratel 1989 IR data. Bile salt enhancement dominates over delayed")
log.info(f"  emptying, yielding a net Cmax increase. The +60% ER food effect is")
log.info(f"  distinct (matrix-controlled release eliminates gastric pH effect).")

# Save summary
food_summary = pd.DataFrame([{
    "parameter": "Cmax_ngmL",
    "fasted": round(met_fasted["Cmax_ngmL"], 2),
    "fed": round(met_fed["Cmax_ngmL"], 2),
    "ratio_fed_fasted": round(cmax_ratio, 3),
    "literature_ratio": "IR: ~1.31 (Bratel 1989); ER: ~1.6",
}, {
    "parameter": "AUC_ngmL_h",
    "fasted": round(met_fasted["AUC_ngmL_h"], 2),
    "fed": round(met_fed["AUC_ngmL_h"], 2),
    "ratio_fed_fasted": round(auc_ratio, 3),
    "literature_ratio": "~1.0",
}, {
    "parameter": "Tmax_h",
    "fasted": round(met_fasted["Tmax_h"], 2),
    "fed": round(met_fed["Tmax_h"], 2),
    "ratio_fed_fasted": round(met_fed["Tmax_h"] / met_fasted["Tmax_h"], 3),
    "literature_ratio": "increased",
}, {
    "parameter": "Fa",
    "fasted": round(met_fasted["Fa"], 3),
    "fed": round(met_fed["Fa"], 3),
    "ratio_fed_fasted": round(met_fed["Fa"] / met_fasted["Fa"], 3),
    "literature_ratio": "increased",
}])
food_summary.to_csv(os.path.join(DIR_TABLES, "food_effect_summary.csv"), index=False)
log.info("  Saved food_effect_summary.csv")

# =============================================================================
# Food effect decomposition — toggle each knob independently
# =============================================================================
log.info("\nDecomposing food effect by mechanism...")

def make_gi_with_single_change(change_name):
    """Create GI parameter set with only one fed-state change applied."""
    gi = [copy.deepcopy(seg) for seg in FASTED_GI]

    if change_name == "gastric_emptying":
        gi[0].transit_time_h = FED_GI[0].transit_time_h  # 0.25 -> 0.65 h
        gi[0].volume_mL = FED_GI[0].volume_mL            # 250 -> 500 mL
    elif change_name == "bile_salts":
        for i in range(len(gi)):
            gi[i].bile_factor = FED_GI[i].bile_factor
    elif change_name == "gastric_pH":
        gi[0].pH = FED_GI[0].pH  # 1.7 -> 3.0
    elif change_name == "all_fed":
        gi = [copy.deepcopy(seg) for seg in FED_GI]

    return gi


decomp_rows = []
knobs = ["gastric_emptying", "bile_salts", "gastric_pH", "all_fed"]
knob_labels = {
    "gastric_emptying": "Delayed Gastric Emptying",
    "bile_salts": "Increased Bile Salts",
    "gastric_pH": "Elevated Gastric pH",
    "all_fed": "All Fed Changes",
}

for knob in knobs:
    gi_mod = make_gi_with_single_change(knob)
    sim_mod = simulate(drug, gi_mod, t_end_h=24.0, dt_h=0.01)
    met_mod = compute_pk_metrics(sim_mod)

    decomp_rows.append({
        "mechanism": knob,
        "label": knob_labels[knob],
        "Cmax_ngmL": round(met_mod["Cmax_ngmL"], 2),
        "Cmax_change_pct": round(
            (met_mod["Cmax_ngmL"] / met_fasted["Cmax_ngmL"] - 1) * 100, 1),
        "AUC_ngmL_h": round(met_mod["AUC_ngmL_h"], 2),
        "AUC_change_pct": round(
            (met_mod["AUC_ngmL_h"] / met_fasted["AUC_ngmL_h"] - 1) * 100, 1),
        "Tmax_h": round(met_mod["Tmax_h"], 2),
    })
    log.info(f"  {knob_labels[knob]}: Cmax={met_mod['Cmax_ngmL']:.2f} "
             f"({decomp_rows[-1]['Cmax_change_pct']:+.1f}%), "
             f"AUC={met_mod['AUC_ngmL_h']:.2f} "
             f"({decomp_rows[-1]['AUC_change_pct']:+.1f}%)")

decomp_df = pd.DataFrame(decomp_rows)
decomp_df.to_csv(os.path.join(DIR_TABLES, "food_effect_decomposition.csv"),
                 index=False)
log.info("  Saved food_effect_decomposition.csv")

# =============================================================================
# Figure 3: Fasted vs Fed PK overlay
# =============================================================================
set_pub_style()
fig, ax = plt.subplots(figsize=(8, 5))

ax.plot(sim_fasted["time_h"], sim_fasted["Cp_ngmL"], "b-", linewidth=2.0,
        label=f"Fasted (Cmax={met_fasted['Cmax_ngmL']:.1f})")
ax.plot(sim_fed["time_h"], sim_fed["Cp_ngmL"], "r-", linewidth=2.0,
        label=f"Fed (Cmax={met_fed['Cmax_ngmL']:.1f})")

# Load observed fed data if available
obs_path = os.path.join(DIR_OBSERVED_EXTRACTED, "mean_profiles_felodipine_po.csv")
if os.path.exists(obs_path):
    obs_df = pd.read_csv(obs_path)
    fed_obs = obs_df[(obs_df["prandial"] == "fed") & (obs_df["is_control"] == True)]
    if len(fed_obs) > 0:
        for study_id, grp in fed_obs.groupby("study_id"):
            ax.plot(grp["time_h"], grp["conc_ngml"], "rs", markersize=4,
                    alpha=0.6, label=f"Obs fed: {grp['study_label'].iloc[0][:25]}")

ax.set_xlabel("Time (h)")
ax.set_ylabel("Plasma concentration (ng/mL)")
ax.set_title("Food Effect on Felodipine PK: Fasted vs Fed (10 mg IR)")

# Annotation box
textstr = (f"Cmax ratio: {cmax_ratio:.2f}\n"
           f"AUC ratio: {auc_ratio:.2f}\n"
           f"Tmax shift: +{tmax_shift:.1f} h")
ax.text(0.98, 0.95, textstr, transform=ax.transAxes, fontsize=9,
        verticalalignment="top", horizontalalignment="right",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="wheat", alpha=0.8))

ax.set_xlim(0, 24)
ax.set_ylim(bottom=0)
ax.legend(fontsize=8, loc="upper right", bbox_to_anchor=(0.75, 1.0))
fig.tight_layout()
fig.savefig(os.path.join(DIR_FIGURES, "pk_fasted_vs_fed.png"), dpi=300)
plt.close(fig)
log.info("  Saved pk_fasted_vs_fed.png")

# =============================================================================
# Figure 4: Food effect decomposition
# =============================================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Exclude "all_fed" from decomposition bars (it's the combined)
decomp_single = decomp_df[decomp_df["mechanism"] != "all_fed"]

colors_bar = ["#e41a1c", "#377eb8", "#4daf4a"]

# Cmax changes
bars1 = ax1.barh(decomp_single["label"], decomp_single["Cmax_change_pct"],
                 color=colors_bar, edgecolor="white", height=0.5)
ax1.axvline(0, color="black", linewidth=0.5)
ax1.set_xlabel("Change in Cmax (%)")
ax1.set_title("Cmax Decomposition")
for bar, val in zip(bars1, decomp_single["Cmax_change_pct"]):
    ax1.text(bar.get_width() + 1, bar.get_y() + bar.get_height() / 2,
             f"{val:+.1f}%", va="center", fontsize=9)

# AUC changes
bars2 = ax2.barh(decomp_single["label"], decomp_single["AUC_change_pct"],
                 color=colors_bar, edgecolor="white", height=0.5)
ax2.axvline(0, color="black", linewidth=0.5)
ax2.set_xlabel("Change in AUC (%)")
ax2.set_title("AUC Decomposition")
for bar, val in zip(bars2, decomp_single["AUC_change_pct"]):
    ax2.text(bar.get_width() + 1, bar.get_y() + bar.get_height() / 2,
             f"{val:+.1f}%", va="center", fontsize=9)

fig.suptitle("Food Effect Decomposition — Individual Mechanism Contributions",
             fontsize=13, fontweight="bold", y=1.02)
fig.tight_layout()
fig.savefig(os.path.join(DIR_FIGURES, "food_effect_decomposition.png"), dpi=300)
plt.close(fig)
log.info("  Saved food_effect_decomposition.png")

log.info("Fed-state simulation and food effect analysis complete.")
