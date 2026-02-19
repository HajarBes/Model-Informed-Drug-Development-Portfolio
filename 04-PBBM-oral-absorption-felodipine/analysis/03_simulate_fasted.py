#!/usr/bin/env python3
# =============================================================================
# 03_simulate_fasted.py — Fasted-State PK Simulation + Qualification
# =============================================================================
#
# Simulates fasted-state PK for 10 mg IR felodipine using the ACAT model
# and compares predicted profiles against observed data from the OSP database.
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
    DrugParams, FASTED_GI,
)

DIR_FIGURES = _setup.DIR_FIGURES
DIR_TABLES = _setup.DIR_TABLES
DIR_OBSERVED_EXTRACTED = _setup.DIR_OBSERVED_EXTRACTED
compute_auc = _setup.compute_auc
compute_cmax_tmax = _setup.compute_cmax_tmax
set_pub_style = _setup.set_pub_style
setup_logging = _setup.setup_logging

log = setup_logging("03_simulate_fasted")

# =============================================================================
# Simulate fasted-state PK
# =============================================================================
log.info("Simulating fasted-state PK (10 mg IR felodipine)...")
drug = DrugParams()
sim_df = simulate(drug, FASTED_GI, t_end_h=24.0, dt_h=0.01)
metrics = compute_pk_metrics(sim_df)

log.info(f"  Predicted Cmax = {metrics['Cmax_ngmL']:.2f} ng/mL")
log.info(f"  Predicted Tmax = {metrics['Tmax_h']:.2f} h")
log.info(f"  Predicted AUC  = {metrics['AUC_ngmL_h']:.2f} ng/mL*h")
log.info(f"  Predicted Fa   = {metrics['Fa']:.3f}")

# Save predicted metrics
metrics_df = pd.DataFrame([metrics])
metrics_df.to_csv(os.path.join(DIR_TABLES, "fasted_pk_metrics.csv"), index=False)
log.info(f"  Saved fasted_pk_metrics.csv")

# =============================================================================
# Mass balance table (QA deliverable)
# =============================================================================
last = sim_df.iloc[-1]
S_total = last["S_total_mg"]
D_total = last["D_total_mg"]
cum_absorbed = last["cum_absorbed_mg"]
fecal_loss = drug.dose_mg - S_total - D_total - cum_absorbed

mb_rows = [
    {"component": "Dose administered", "mass_mg": round(drug.dose_mg, 4),
     "pct_dose": 100.0},
    {"component": "Absorbed through intestinal wall", "mass_mg": round(cum_absorbed, 4),
     "pct_dose": round(cum_absorbed / drug.dose_mg * 100, 2)},
    {"component": "Remaining undissolved (GI)", "mass_mg": round(S_total, 4),
     "pct_dose": round(S_total / drug.dose_mg * 100, 2)},
    {"component": "Remaining dissolved (GI)", "mass_mg": round(D_total, 4),
     "pct_dose": round(D_total / drug.dose_mg * 100, 2)},
    {"component": "Fecal transit loss", "mass_mg": round(max(fecal_loss, 0), 4),
     "pct_dose": round(max(fecal_loss, 0) / drug.dose_mg * 100, 2)},
]
mb_df = pd.DataFrame(mb_rows)
mb_df.to_csv(os.path.join(DIR_TABLES, "mass_balance_fasted.csv"), index=False)
log.info(f"  Saved mass_balance_fasted.csv")
log.info(f"    Absorbed: {cum_absorbed:.4f} mg ({cum_absorbed/drug.dose_mg*100:.2f}%)")
log.info(f"    GI remaining: {S_total + D_total:.4f} mg")
log.info(f"    Fecal loss: {max(fecal_loss,0):.4f} mg ({max(fecal_loss,0)/drug.dose_mg*100:.2f}%)")

# =============================================================================
# Bioavailability breakdown table
# =============================================================================
Fa_pred = metrics["Fa"]
Fg_pred = drug.Fg
Fh_pred = drug.Fh
Foral_pred = Fa_pred * Fg_pred * Fh_pred

ba_rows = [
    {"parameter": "Fa (fraction absorbed)", "predicted": round(Fa_pred, 4),
     "literature": "~1.0 (BCS II, high Peff)"},
    {"parameter": "Fg (gut wall availability)", "predicted": round(Fg_pred, 2),
     "literature": "0.50 (Lundahl 1997)"},
    {"parameter": "Fh (hepatic availability)", "predicted": round(Fh_pred, 2),
     "literature": "0.50 (Lundahl 1997)"},
    {"parameter": "Foral (overall bioavailability)", "predicted": round(Foral_pred, 4),
     "literature": "0.15-0.22 (Lundahl 1997)"},
]
ba_df = pd.DataFrame(ba_rows)
ba_df.to_csv(os.path.join(DIR_TABLES, "bioavailability_breakdown.csv"), index=False)
log.info(f"  Saved bioavailability_breakdown.csv")
log.info(f"    Fa={Fa_pred:.4f}, Fg={Fg_pred:.2f}, Fh={Fh_pred:.2f}, "
         f"Foral={Foral_pred:.4f} (lit: 0.15-0.22)")

# =============================================================================
# Load observed data + qualification set classification
# =============================================================================
obs_path = os.path.join(DIR_OBSERVED_EXTRACTED, "mean_profiles_felodipine_po.csv")
qual_arms_path = os.path.join(DIR_TABLES, "qualification_arms.csv")

if os.path.exists(obs_path):
    obs_df = pd.read_csv(obs_path)
    # Filter: control arms, fasted or unknown prandial state
    qual_df = obs_df[
        (obs_df["is_control"] == True) &
        (obs_df["prandial"].isin(["fasted", "unknown"]))
    ].copy()
    log.info(f"  Loaded {len(qual_df)} fasted/control data points for qualification")

    # Load qualification set classification
    qual_set_ids = set()      # in_qualification_set — used for 2-fold metrics
    ir_study_ids = set()
    er_study_ids = set()
    dose_by_study = {}        # study_id -> dose_mg for dose normalization
    if os.path.exists(qual_arms_path):
        qa = pd.read_csv(qual_arms_path)
        qual_set_ids = set(qa[qa["in_qualification_set"] == True]["study_id"])
        ir_study_ids = set(qa[qa["formulation"].isin(["IR", "solution"])]["study_id"])
        er_study_ids = set(qa[qa["formulation"] == "ER"]["study_id"])
        for _, row in qa.iterrows():
            try:
                dose_by_study[row["study_id"]] = float(row["dose_mg"])
            except (ValueError, TypeError):
                dose_by_study[row["study_id"]] = 10.0
        log.info(f"  Qualification set (IR/solution, fasted): {len(qual_set_ids)} arms")
        log.info(f"  ER arms (informational overlay): {len(er_study_ids)}")
else:
    qual_df = pd.DataFrame()
    qual_set_ids = set()
    ir_study_ids = set()
    er_study_ids = set()
    dose_by_study = {}
    log.info("  No observed data file found — skipping overlay")

# =============================================================================
# Per-study qualification
# =============================================================================
# Dose-normalization strategy (Option B):
# All predictions are at 10 mg. For studies dosed at other strengths (e.g. 5 mg),
# dose-normalize observed profiles UP to 10 mg equivalent assuming linear PK.
# obs_DN = obs × (10 / dose_actual).  pred stays at 10 mg.
REFERENCE_DOSE = 10.0

qual_rows = []
if len(qual_df) > 0:
    studies = qual_df.groupby("study_id")

    for study_id, grp in studies:
        obs_t_raw = grp["time_h"].values
        obs_c = grp["conc_ngml"].values

        if len(obs_t_raw) < 3:
            continue

        # Multi-dose time offset: shift observed time to relative (start from 0)
        # e.g. Blychert 1990 has times 120-132h → shift to 0-12h
        t_offset = obs_t_raw.min()
        obs_t = obs_t_raw - t_offset

        # Study dose
        dose_val = dose_by_study.get(study_id, 10.0)

        # Dose-normalize observed to 10 mg (linear PK assumption)
        dn_factor = REFERENCE_DOSE / dose_val if dose_val > 0 else 1.0
        obs_c_dn = obs_c * dn_factor

        # Observed metrics (dose-normalized, relative time)
        obs_cmax_dn, obs_tmax = compute_cmax_tmax(obs_t, obs_c_dn)
        obs_auc_dn = compute_auc(obs_t, obs_c_dn)

        # Predicted at observed time points (interpolate on relative time)
        pred_c = np.interp(obs_t, sim_df["time_h"].values,
                           sim_df["Cp_ngmL"].values)
        pred_cmax = metrics["Cmax_ngmL"]
        # AUC over the SAME time window as observed (not full 0-24 h)
        pred_auc = compute_auc(obs_t, pred_c)

        # Metrics (pred at 10mg vs obs dose-normalized to 10mg)
        rmse = np.sqrt(np.mean((pred_c - obs_c_dn) ** 2))
        if obs_cmax_dn > 0:
            mpe_pct = (pred_cmax - obs_cmax_dn) / obs_cmax_dn * 100
        else:
            mpe_pct = float("nan")

        cmax_ratio = pred_cmax / obs_cmax_dn if obs_cmax_dn > 0 else float("nan")
        auc_ratio = pred_auc / obs_auc_dn if obs_auc_dn > 0 else float("nan")

        label = grp["study_label"].iloc[0]
        # Classify formulation
        if study_id in ir_study_ids:
            form_class = "IR/solution"
        elif study_id in er_study_ids:
            form_class = "ER"
        else:
            form_class = "unknown"

        in_qual = study_id in qual_set_ids

        qual_rows.append({
            "study_id": study_id,
            "study_label": label,
            "formulation": form_class,
            "dose_mg": dose_val,
            "dn_factor": round(dn_factor, 2),
            "n_points": len(obs_t),
            "obs_Cmax_DN": round(obs_cmax_dn, 2),
            "pred_Cmax": round(pred_cmax, 2),
            "Cmax_ratio": round(cmax_ratio, 3),
            "obs_AUC_DN": round(obs_auc_dn, 2),
            "pred_AUC": round(pred_auc, 2),
            "AUC_ratio": round(auc_ratio, 3),
            "RMSE": round(rmse, 2),
            "MPE_pct": round(mpe_pct, 1),
            "in_qualification_set": in_qual,
            "within_2fold": 0.5 <= auc_ratio <= 2.0,
        })

    qual_summary = pd.DataFrame(qual_rows)
    qual_summary.to_csv(os.path.join(DIR_TABLES, "fasted_qualification.csv"),
                        index=False)
    # Report metrics using ONLY in_qualification_set arms
    qs = qual_summary[qual_summary["in_qualification_set"] == True]
    n_qs_within = qs["within_2fold"].sum() if len(qs) > 0 else 0
    n_qs_total = len(qs)
    # Informational: ER overlay
    er_qual = qual_summary[qual_summary["formulation"] == "ER"]
    n_er_within = er_qual["within_2fold"].sum() if len(er_qual) > 0 else 0
    n_er_total = len(er_qual)
    log.info(f"  Qualification set (IR/solution, fasted): "
             f"{n_qs_within}/{n_qs_total} within 2-fold AUC ratio")
    log.info(f"  ER overlay (informational):              "
             f"{n_er_within}/{n_er_total} within 2-fold")
    log.info(f"  Dose normalization: all observed profiles scaled to {REFERENCE_DOSE} mg "
             f"(linear PK)")
else:
    qual_summary = pd.DataFrame()

# =============================================================================
# Figure 1: PK overlay (fasted)
# =============================================================================
set_pub_style()
fig, ax = plt.subplots(figsize=(8, 5))

ax.plot(sim_df["time_h"], sim_df["Cp_ngmL"], "b-", linewidth=2.0,
        label="ACAT predicted (10 mg)", zorder=5)

if len(qual_df) > 0:
    # Separate qualification set vs informational overlay
    qs_obs = qual_df[qual_df["study_id"].isin(qual_set_ids)]
    er_obs = qual_df[qual_df["study_id"].isin(er_study_ids)]
    other_obs = qual_df[~qual_df["study_id"].isin(qual_set_ids | er_study_ids)]

    # Helpers: dose-normalize + time-shift (multi-dose offset removal)
    def dn_conc(grp):
        dose = dose_by_study.get(grp["study_id"].iloc[0], 10.0)
        return grp["conc_ngml"] * (REFERENCE_DOSE / dose)

    def dn_sd(grp):
        dose = dose_by_study.get(grp["study_id"].iloc[0], 10.0)
        return grp["sd_ngml"].replace("", 0).astype(float) * (REFERENCE_DOSE / dose)

    def rel_time(grp):
        return grp["time_h"] - grp["time_h"].min()

    # Plot qualification set studies with solid markers
    if len(qs_obs) > 0:
        qs_colors = plt.cm.tab10(np.linspace(0, 0.5, qs_obs["study_id"].nunique()))
        for idx, (study_id, grp) in enumerate(qs_obs.groupby("study_id")):
            dose = dose_by_study.get(study_id, 10.0)
            dose_tag = f" [{dose:.0f}mg DN]" if dose != 10.0 else ""
            label_str = str(grp["study_label"].iloc[0])[:25] + f" [IR]{dose_tag}"
            ax.errorbar(rel_time(grp), dn_conc(grp),
                        yerr=dn_sd(grp),
                        fmt="o", markersize=5, capsize=2,
                        color=qs_colors[idx], alpha=0.9,
                        label=label_str)

    # Plot ER studies (informational overlay) with open markers — dose-normalized
    if len(er_obs) > 0:
        er_colors = plt.cm.Set2(np.linspace(0, 1, min(er_obs["study_id"].nunique(), 8)))
        for idx, (study_id, grp) in enumerate(er_obs.groupby("study_id")):
            label_str = str(grp["study_label"].iloc[0])[:25] + " [ER]"
            ax.errorbar(rel_time(grp), dn_conc(grp),
                        yerr=dn_sd(grp),
                        fmt="s", markersize=3, capsize=1,
                        color=er_colors[idx % len(er_colors)], alpha=0.4,
                        label=label_str)

    # Plot unclassified — dose-normalized
    if len(other_obs) > 0:
        for study_id, grp in other_obs.groupby("study_id"):
            label_str = str(grp["study_label"].iloc[0])[:25]
            ax.errorbar(rel_time(grp), dn_conc(grp),
                        yerr=dn_sd(grp),
                        fmt="d", markersize=3, capsize=1,
                        color="gray", alpha=0.4,
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

# =============================================================================
# Figure 2: Absorption by segment
# =============================================================================
seg_df = segment_absorption_breakdown(drug, FASTED_GI, sim_df)
seg_df.to_csv(os.path.join(DIR_TABLES, "absorption_by_segment_fasted.csv"),
              index=False)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Bar chart of absorption fraction
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

# Cumulative solid + dissolved along GI tract
t_end_idx = len(sim_df) - 1
seg_names = [s.name for s in FASTED_GI]
solid_at_end = [float(sim_df[f"S_{i}"].iloc[t_end_idx]) for i in range(7)]
diss_at_end = [float(sim_df[f"D_{i}"].iloc[t_end_idx]) for i in range(7)]

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
# Figure: Dissolved concentration vs C_sat over time (per segment)
# =============================================================================
log.info("Generating dissolved concentration vs C_sat plot...")

# Select key absorbing segments (skip stomach — no absorption)
plot_segments = [
    (1, "Duodenum"),
    (2, "Jejunum 1"),
    (3, "Jejunum 2"),
    (4, "Ileum 1"),
    (6, "Colon"),
]

fig, axes = plt.subplots(len(plot_segments), 1, figsize=(10, 3.0 * len(plot_segments)),
                         sharex=True)

sim_t = sim_df["time_h"].values

for ax_i, (seg_idx, seg_label) in enumerate(plot_segments):
    seg = FASTED_GI[seg_idx]
    D_i = sim_df[f"D_{seg_idx}"].values  # dissolved drug mass (mg)
    C_diss = D_i / seg.volume_mL * 1000.0  # mg/mL -> µg/mL

    C_sat = saturation_solubility(drug.S0_mg_mL, drug.pKa,
                                   seg.pH, seg.bile_factor) * 1000.0  # µg/mL

    ax = axes[ax_i]
    ax.plot(sim_t, C_diss, "b-", linewidth=1.5,
            label=f"C_diss (dissolved)")
    ax.axhline(C_sat, color="red", linestyle="--", linewidth=1.2,
               label=f"C_sat = {C_sat:.2f} µg/mL")

    # Shade supersaturation region (C_diss > C_sat)
    mask_super = C_diss > C_sat
    if mask_super.any():
        ax.fill_between(sim_t, C_sat, C_diss, where=mask_super,
                         alpha=0.15, color="red", label="Supersaturated")

    # Shade absorption-sink region (C_diss < C_sat)
    mask_under = C_diss <= C_sat
    if mask_under.any():
        ax.fill_between(sim_t, C_diss, C_sat, where=mask_under,
                         alpha=0.08, color="green", label="Dissolution driving force")

    ax.set_ylabel("Conc (µg/mL)")
    ax.set_title(f"{seg_label}  (pH={seg.pH}, bile={seg.bile_factor}x, "
                 f"V={seg.volume_mL:.0f} mL)", fontsize=10, fontweight="bold")
    ax.legend(fontsize=7, loc="upper right")
    ax.set_ylim(bottom=0)

axes[-1].set_xlabel("Time (h)")
axes[-1].set_xlim(0, 6)  # Zoom to first 6 h (most dissolution action)

fig.suptitle("Dissolved Concentration vs Saturation Solubility (C_sat) — Fasted State\n"
             "Green = dissolution driving force (C_sat − C_diss); "
             "Red = supersaturation (precipitation)",
             fontsize=11, fontweight="bold", y=1.01)
fig.tight_layout()
fig.savefig(os.path.join(DIR_FIGURES, "cdiss_vs_csat.png"), dpi=300)
plt.close(fig)
log.info("  Saved cdiss_vs_csat.png")

# =============================================================================
# Fg variability simulation
# =============================================================================
log.info("Running Fg variability simulation...")

np.random.seed(_setup.SEED_VPOP)
N_FG = 50
fg_values = np.clip(np.random.normal(loc=0.50, scale=0.10, size=N_FG), 0.15, 0.85)
fg_results = []

for fg_i in fg_values:
    drug_i = DrugParams()
    drug_i.Fg = fg_i
    sim_i = simulate(drug_i, FASTED_GI, t_end_h=24.0, dt_h=0.05)
    met_i = compute_pk_metrics(sim_i)
    fg_results.append({
        "Fg": round(fg_i, 4),
        "Cmax_ngmL": round(met_i["Cmax_ngmL"], 2),
        "AUC_ngmL_h": round(met_i["AUC_ngmL_h"], 2),
        "Foral": round(met_i["Fa"] * fg_i * drug_i.Fh, 4),
    })

fg_var_df = pd.DataFrame(fg_results)
fg_var_df.to_csv(os.path.join(DIR_TABLES, "fg_variability.csv"), index=False)
log.info(f"  Saved fg_variability.csv ({N_FG} samples)")
log.info(f"    Cmax range: {fg_var_df['Cmax_ngmL'].min():.1f} – "
         f"{fg_var_df['Cmax_ngmL'].max():.1f} ng/mL")
log.info(f"    AUC range:  {fg_var_df['AUC_ngmL_h'].min():.1f} – "
         f"{fg_var_df['AUC_ngmL_h'].max():.1f} ng/mL*h")

# Figure: Fg vs Cmax/AUC scatter
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
ax1.scatter(fg_var_df["Fg"], fg_var_df["Cmax_ngmL"], alpha=0.7, s=30, c="#377eb8")
ax1.axvline(0.50, color="red", linestyle="--", linewidth=0.8, label="Fg = 0.50 (ref)")
ax1.set_xlabel("Fg (gut wall availability)")
ax1.set_ylabel("Cmax (ng/mL)")
ax1.set_title("Impact of Fg Variability on Cmax")
ax1.legend(fontsize=8)

ax2.scatter(fg_var_df["Fg"], fg_var_df["AUC_ngmL_h"], alpha=0.7, s=30, c="#e41a1c")
ax2.axvline(0.50, color="red", linestyle="--", linewidth=0.8, label="Fg = 0.50 (ref)")
ax2.set_xlabel("Fg (gut wall availability)")
ax2.set_ylabel("AUC (ng/mL*h)")
ax2.set_title("Impact of Fg Variability on AUC")
ax2.legend(fontsize=8)

fig.suptitle("Fg Variability Simulation (N=50, Normal distribution, CV=20%)",
             fontsize=12, fontweight="bold", y=1.02)
fig.tight_layout()
fig.savefig(os.path.join(DIR_FIGURES, "fg_variability.png"), dpi=300)
plt.close(fig)
log.info("  Saved fg_variability.png")

# =============================================================================
# Figure: Cross-study forest plot of pred/obs AUC ratio
# =============================================================================
if len(qual_rows) > 0:
    log.info("Generating cross-study forest plot...")

    # Sort: qualification set first, then ER
    qdf = qual_summary.copy()
    qdf["sort_key"] = qdf["in_qualification_set"].apply(lambda x: 0 if x else 1)
    qdf = qdf.sort_values(["sort_key", "AUC_ratio"]).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(8, max(6, len(qdf) * 0.35)))
    y_pos = np.arange(len(qdf))
    labels = [f"{r['study_label']} [{r['formulation']}]" for _, r in qdf.iterrows()]

    # Color by qualification status
    colors_fp = []
    for _, r in qdf.iterrows():
        if r["in_qualification_set"]:
            colors_fp.append("#377eb8" if r["within_2fold"] else "#e41a1c")
        else:
            colors_fp.append("#999999")

    ax.barh(y_pos, qdf["AUC_ratio"], height=0.6, color=colors_fp,
            edgecolor="white", linewidth=0.5, alpha=0.85)

    # 2-fold acceptance bounds
    ax.axvline(0.5, color="red", linestyle="--", linewidth=1.0, alpha=0.7)
    ax.axvline(2.0, color="red", linestyle="--", linewidth=1.0, alpha=0.7,
               label="2-fold bounds")
    ax.axvline(1.0, color="black", linestyle="-", linewidth=0.8, alpha=0.5,
               label="Unity")

    # Shade acceptance region
    ax.axvspan(0.5, 2.0, alpha=0.06, color="green")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Predicted / Observed AUC Ratio")
    ax.set_title("Cross-Study Forest Plot: Pred/Obs AUC Ratio\n"
                 "Blue = qualification set (IR), Gray = ER (informational)")

    # Add ratio text
    for i, (_, r) in enumerate(qdf.iterrows()):
        ax.text(r["AUC_ratio"] + 0.05, i, f"{r['AUC_ratio']:.2f}",
                va="center", fontsize=7)

    ax.set_xlim(0, max(qdf["AUC_ratio"].max() * 1.15, 2.5))
    ax.legend(fontsize=8, loc="lower right")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(os.path.join(DIR_FIGURES, "forest_plot_auc_ratio.png"), dpi=300)
    plt.close(fig)
    log.info("  Saved forest_plot_auc_ratio.png")

log.info("Fasted-state simulation complete.")
