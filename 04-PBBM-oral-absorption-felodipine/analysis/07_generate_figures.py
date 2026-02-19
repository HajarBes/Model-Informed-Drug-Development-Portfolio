#!/usr/bin/env python3
# =============================================================================
# 07_generate_figures.py — Consolidated Publication Figures
# =============================================================================
#
# Regenerates all 8 figures with consistent styling. Some figures are also
# produced by earlier scripts; this script consolidates and adds the
# pH-solubility profile (Figure 8).
#
# Author: Hajar Besbassi

from __future__ import annotations

import os
import sys
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
    simulate, compute_pk_metrics, segment_absorption_breakdown,
    saturation_solubility,
    DrugParams, GISegment, FASTED_GI, FED_GI,
)

DIR_FIGURES = _setup.DIR_FIGURES
DIR_TABLES = _setup.DIR_TABLES
DIR_OBSERVED_EXTRACTED = _setup.DIR_OBSERVED_EXTRACTED
set_pub_style = _setup.set_pub_style
setup_logging = _setup.setup_logging

log = setup_logging("07_figures")

set_pub_style()

drug = DrugParams()

# =============================================================================
# Figure 8: pH-Dependent Solubility Profile (Henderson-Hasselbalch)
# =============================================================================
log.info("Generating Figure 8: pH-dependent solubility profile...")

pH_range = np.linspace(1.0, 8.0, 200)
sol_no_bile = [saturation_solubility(drug.S0_mg_mL, drug.pKa, pH, 1.0) * 1000
               for pH in pH_range]  # convert to µg/mL
sol_bile_2x = [saturation_solubility(drug.S0_mg_mL, drug.pKa, pH, 2.0) * 1000
               for pH in pH_range]
sol_bile_15x = [saturation_solubility(drug.S0_mg_mL, drug.pKa, pH, 15.0) * 1000
                for pH in pH_range]

fig, ax = plt.subplots(figsize=(8, 5))
ax.semilogy(pH_range, sol_no_bile, "k-", linewidth=2.0, label="No bile salts (1x)")
ax.semilogy(pH_range, sol_bile_2x, "b-", linewidth=2.0, label="Fasted bile (2x)")
ax.semilogy(pH_range, sol_bile_15x, "r-", linewidth=2.0, label="Fed bile (15x)")

# Mark GI segment pH values
gi_pHs = {"Stomach": 1.7, "Duodenum": 6.0, "Jejunum": 6.5, "Ileum": 7.1, "Colon": 6.5}
for name, pH_val in gi_pHs.items():
    sol_val = saturation_solubility(drug.S0_mg_mL, drug.pKa, pH_val, 2.0) * 1000
    ax.plot(pH_val, sol_val, "^", color="#e41a1c", markersize=8, zorder=5)
    ax.annotate(name, (pH_val, sol_val), textcoords="offset points",
                xytext=(5, 8), fontsize=7, color="#e41a1c")

ax.axhline(drug.S0_ug_mL, color="gray", linestyle=":", linewidth=0.8)
ax.text(7.5, drug.S0_ug_mL * 1.2, f"S₀ = {drug.S0_ug_mL} µg/mL",
        fontsize=8, color="gray", ha="right")

ax.axvline(drug.pKa, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)
ax.text(drug.pKa + 0.1, ax.get_ylim()[1] * 0.5,
        f"pKa = {drug.pKa}", fontsize=8, color="gray", rotation=90, va="center")

ax.set_xlabel("pH")
ax.set_ylabel("Saturation Solubility (µg/mL)")
ax.set_title("Felodipine pH-Dependent Solubility (Henderson-Hasselbalch + Bile Enhancement)")
ax.legend(fontsize=9)
ax.set_xlim(1, 8)
fig.tight_layout()
fig.savefig(os.path.join(DIR_FIGURES, "drug_solubility_profile.png"), dpi=300)
plt.close(fig)
log.info("  Saved drug_solubility_profile.png")

# =============================================================================
# Re-generate Figures 1-2 with consistent styling (fasted PK + absorption)
# =============================================================================
log.info("Regenerating Figures 1-2 (fasted PK overlay + absorption)...")

sim_fasted = simulate(drug, FASTED_GI, t_end_h=24.0, dt_h=0.01)
met_fasted = compute_pk_metrics(sim_fasted)

# Figure 1: Fasted PK overlay
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(sim_fasted["time_h"], sim_fasted["Cp_ngmL"], "b-", linewidth=2.0,
        label="ACAT predicted (10 mg IR)", zorder=5)

obs_path = os.path.join(DIR_OBSERVED_EXTRACTED, "mean_profiles_felodipine_po.csv")
qual_arms_path = os.path.join(DIR_TABLES, "qualification_arms.csv")
if os.path.exists(obs_path):
    obs_df = pd.read_csv(obs_path)
    qual_df = obs_df[
        (obs_df["is_control"] == True) &
        (obs_df["prandial"].isin(["fasted", "unknown"]))
    ]
    # Load dose info for dose normalization
    dose_map = {}
    if os.path.exists(qual_arms_path):
        qa_tmp = pd.read_csv(qual_arms_path)
        for _, r in qa_tmp.iterrows():
            try:
                dose_map[r["study_id"]] = float(r["dose_mg"])
            except (ValueError, TypeError):
                dose_map[r["study_id"]] = 10.0
    if len(qual_df) > 0:
        colors = plt.cm.Set2(np.linspace(0, 1, min(qual_df["study_id"].nunique(), 8)))
        for idx, (study_id, grp) in enumerate(qual_df.groupby("study_id")):
            label_str = str(grp["study_label"].iloc[0])[:30]
            dose = dose_map.get(study_id, 10.0)
            dn = 10.0 / dose if dose > 0 else 1.0
            # Time-shift multi-dose studies to relative time (start from 0)
            t_rel = grp["time_h"] - grp["time_h"].min()
            ax.errorbar(t_rel, grp["conc_ngml"] * dn,
                        yerr=grp["sd_ngml"].replace("", 0).astype(float) * dn,
                        fmt="o", markersize=4, capsize=2,
                        color=colors[idx % len(colors)], alpha=0.7,
                        label=label_str)

ax.set_xlabel("Time (h)")
ax.set_ylabel("Plasma concentration (ng/mL)")
ax.set_title("Fasted-State PK: ACAT Predicted vs Observed\n"
             "Felodipine 10 mg IR (observed dose-normalized to 10 mg)")
ax.set_xlim(0, 24)
ax.set_ylim(bottom=0)
ax.legend(fontsize=7, loc="upper right", framealpha=0.9)
fig.tight_layout()
fig.savefig(os.path.join(DIR_FIGURES, "pk_overlay_fasted.png"), dpi=300)
plt.close(fig)
log.info("  Saved pk_overlay_fasted.png")

# Figure 2: Absorption by segment
seg_df = segment_absorption_breakdown(drug, FASTED_GI, sim_fasted)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
colors_seg = ["#999999", "#e41a1c", "#377eb8", "#4daf4a",
              "#984ea3", "#ff7f00", "#a65628"]
bars = ax1.bar(seg_df["segment"], seg_df["pct_dose"], color=colors_seg,
               edgecolor="white", linewidth=0.5)
ax1.set_ylabel("Absorption (% of dose)")
ax1.set_xlabel("GI Segment")
ax1.set_title("Regional Absorption — Fasted State")
for bar, val in zip(bars, seg_df["pct_dose"]):
    if val > 1:
        ax1.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                 f"{val:.1f}%", ha="center", va="bottom", fontsize=8)

t_end_idx = len(sim_fasted) - 1
seg_names = [s.name for s in FASTED_GI]
solid_at_end = [float(sim_fasted[f"S_{i}"].iloc[t_end_idx]) for i in range(7)]
diss_at_end = [float(sim_fasted[f"D_{i}"].iloc[t_end_idx]) for i in range(7)]

x = np.arange(7)
ax2.bar(x, solid_at_end, label="Undissolved solid", color="#d62728", alpha=0.8)
ax2.bar(x, diss_at_end, bottom=solid_at_end, label="Dissolved",
        color="#2ca02c", alpha=0.8)
ax2.set_xticks(x)
ax2.set_xticklabels(seg_names, rotation=45, ha="right")
ax2.set_ylabel("Drug remaining (mg)")
ax2.set_xlabel("GI Segment")
ax2.set_title("Drug Mass in GI Tract at t=24 h")
ax2.legend()

fig.suptitle("Felodipine ACAT Model — Absorption Characterization (Fasted)",
             fontsize=13, fontweight="bold", y=1.02)
fig.tight_layout()
fig.savefig(os.path.join(DIR_FIGURES, "absorption_by_segment_fasted.png"), dpi=300)
plt.close(fig)
log.info("  Saved absorption_by_segment_fasted.png")

# =============================================================================
# Re-generate Figures 3-4 (food effect)
# =============================================================================
log.info("Regenerating Figures 3-4 (food effect)...")

sim_fed = simulate(drug, FED_GI, t_end_h=24.0, dt_h=0.01)
met_fed = compute_pk_metrics(sim_fed)
cmax_ratio = met_fed["Cmax_ngmL"] / met_fasted["Cmax_ngmL"]
auc_ratio = met_fed["AUC_ngmL_h"] / met_fasted["AUC_ngmL_h"]
tmax_shift = met_fed["Tmax_h"] - met_fasted["Tmax_h"]

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(sim_fasted["time_h"], sim_fasted["Cp_ngmL"], "b-", linewidth=2.0,
        label=f"Fasted (Cmax={met_fasted['Cmax_ngmL']:.1f})")
ax.plot(sim_fed["time_h"], sim_fed["Cp_ngmL"], "r-", linewidth=2.0,
        label=f"Fed (Cmax={met_fed['Cmax_ngmL']:.1f})")

textstr = (f"Cmax ratio: {cmax_ratio:.2f}\n"
           f"AUC ratio: {auc_ratio:.2f}\n"
           f"Tmax shift: +{tmax_shift:.1f} h")
ax.text(0.98, 0.95, textstr, transform=ax.transAxes, fontsize=9,
        verticalalignment="top", horizontalalignment="right",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="wheat", alpha=0.8))

ax.set_xlabel("Time (h)")
ax.set_ylabel("Plasma concentration (ng/mL)")
ax.set_title("Food Effect on Felodipine PK: Fasted vs Fed (10 mg IR)")
ax.set_xlim(0, 24)
ax.set_ylim(bottom=0)
ax.legend(fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(DIR_FIGURES, "pk_fasted_vs_fed.png"), dpi=300)
plt.close(fig)
log.info("  Saved pk_fasted_vs_fed.png")

# Figure 4 — decomposition (from saved CSV)
decomp_path = os.path.join(DIR_TABLES, "food_effect_decomposition.csv")
if os.path.exists(decomp_path):
    decomp_df = pd.read_csv(decomp_path)
    decomp_single = decomp_df[decomp_df["mechanism"] != "all_fed"]
    colors_bar = ["#e41a1c", "#377eb8", "#4daf4a"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    bars1 = ax1.barh(decomp_single["label"], decomp_single["Cmax_change_pct"],
                     color=colors_bar, edgecolor="white", height=0.5)
    ax1.axvline(0, color="black", linewidth=0.5)
    ax1.set_xlabel("Change in Cmax (%)")
    ax1.set_title("Cmax Decomposition")
    for bar, val in zip(bars1, decomp_single["Cmax_change_pct"]):
        ax1.text(bar.get_width() + 1, bar.get_y() + bar.get_height() / 2,
                 f"{val:+.1f}%", va="center", fontsize=9)

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

# =============================================================================
# Summary
# =============================================================================
fig_list = [
    "pk_overlay_fasted.png",
    "absorption_by_segment_fasted.png",
    "pk_fasted_vs_fed.png",
    "food_effect_decomposition.png",
    "tornado_sensitivity.png",
    "dissolution_exposure_linkage.png",
    "formulation_comparability.png",
    "drug_solubility_profile.png",
    "fg_variability.png",
    "forest_plot_auc_ratio.png",
]

log.info("\nFigure inventory:")
for i, fname in enumerate(fig_list, 1):
    path = os.path.join(DIR_FIGURES, fname)
    exists = os.path.exists(path)
    status = "OK" if exists else "MISSING"
    log.info(f"  [{status}] Fig {i}: {fname}")

log.info("\nAll figures generated.")
