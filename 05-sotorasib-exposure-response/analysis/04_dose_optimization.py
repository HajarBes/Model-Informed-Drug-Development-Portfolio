#!/usr/bin/env python3
"""
04_dose_optimization.py — Dose Optimization & Benefit-Risk
===========================================================
Side-by-side comparison: 240 mg vs 960 mg
  - Probability of target attainment
  - Therapeutic index: P(ORR) / P(Grade3+)
  - NNT vs NNH
  - Net clinical benefit visualization

Regulatory framing: connects to FDA PMR, ODAC vote, JCO editorials.

Outputs:
  - outputs/tables/dose_optimization_summary.csv
  - outputs/tables/benefit_risk_comparison.csv
  - figures/dose_optimization_dashboard.png
  - figures/benefit_risk_overlay.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import importlib.util as _ilu

_setup_path = Path(__file__).resolve().parent / "00_setup.py"
_spec = _ilu.spec_from_file_location("setup_00", _setup_path)
_mod = _ilu.module_from_spec(_spec)
sys.modules["setup_00"] = _mod
_spec.loader.exec_module(_mod)

DIR_FIGURES = _mod.DIR_FIGURES
DIR_TABLES = _mod.DIR_TABLES
DOSE_LEVELS = _mod.DOSE_LEVELS
DOSE_COLORS = _mod.DOSE_COLORS
COLORS = _mod.COLORS
CLINICAL = _mod.CLINICAL
get_logger = _mod.get_logger
apply_publication_theme = _mod.apply_publication_theme

apply_publication_theme()
log = get_logger("04_dose_opt")

# ── Load Data ────────────────────────────────────────────────────────
df = pd.read_csv(DIR_TABLES / "virtual_patients.csv")
log.info(f"Loaded {len(df)} virtual patients")

# Load efficacy results (need responder column)
# Re-simulate response using the same model from 02
from scipy.special import expit as _expit

# Efficacy model parameters (must match 02_exposure_efficacy.py)
E0 = 0.12
EMAX = 0.25
EC50 = 30.0
GAMMA = 1.2

def p_response(auc):
    if auc <= 0:
        return E0
    return E0 + EMAX * auc**GAMMA / (EC50**GAMMA + auc**GAMMA)

# Safety model parameters (must match 03_exposure_safety.py)
ALPHA = -3.5
BETA_CMAX = 0.08
BETA_CPI = 1.8

def p_hepatotox(cmax, prior_cpi):
    logit_p = ALPHA + BETA_CMAX * cmax + BETA_CPI * prior_cpi
    return _expit(logit_p)

# Compute probabilities
RNG = np.random.default_rng(20240801)
df["p_response"] = df["auc_ss"].apply(p_response)
df["responder"] = RNG.binomial(1, df["p_response"].values)
df["p_hepatotox"] = df.apply(
    lambda r: p_hepatotox(r["cmax_ss"], r["prior_cpi"]), axis=1
)
df["hepatotox"] = RNG.binomial(1, df["p_hepatotox"].values)

# ── Dose-Level Summaries ─────────────────────────────────────────────
log.info("Computing dose-level benefit-risk summaries...")

summary_rows = []
for dose in DOSE_LEVELS:
    sub = df[df["dose_mg"] == dose]
    orr = sub["responder"].mean()
    hepatotox = sub["hepatotox"].mean()
    gm_auc = np.exp(np.mean(np.log(sub["auc_ss"])))
    gm_cmax = np.exp(np.mean(np.log(sub["cmax_ss"])))

    # Therapeutic index = P(ORR) / P(hepatotox)
    ti = orr / hepatotox if hepatotox > 0 else np.inf
    # NNT = 1 / ORR (number needed to treat for one responder)
    nnt = 1 / orr if orr > 0 else np.inf
    # NNH = 1 / hepatotox (number needed to harm)
    nnh = 1 / hepatotox if hepatotox > 0 else np.inf
    # Net clinical benefit = ORR - hepatotox rate
    ncb = orr - hepatotox

    summary_rows.append({
        "dose_mg": dose,
        "n": len(sub),
        "gm_auc_ss": round(gm_auc, 1),
        "gm_cmax_ss": round(gm_cmax, 2),
        "orr_pct": round(orr * 100, 1),
        "hepatotox_pct": round(hepatotox * 100, 1),
        "therapeutic_index": round(ti, 2),
        "nnt": round(nnt, 1),
        "nnh": round(nnh, 1),
        "net_clinical_benefit_pct": round(ncb * 100, 1),
    })

df_summary = pd.DataFrame(summary_rows)
df_summary.to_csv(DIR_TABLES / "dose_optimization_summary.csv", index=False)
log.info("\nDose Optimization Summary:")
log.info("\n" + df_summary.to_string(index=False))

# ── Head-to-Head: 240 mg vs 960 mg ──────────────────────────────────
log.info("\n240 mg vs 960 mg comparison:")
row_240 = df_summary[df_summary["dose_mg"] == 240].iloc[0]
row_960 = df_summary[df_summary["dose_mg"] == 960].iloc[0]

comparison_data = {
    "metric": [
        "GM AUCss (hr*ug/mL)", "GM Cmax,ss (ug/mL)",
        "ORR (%)", "Grade 3+ Hepatotox (%)",
        "Therapeutic Index", "NNT", "NNH",
        "Net Clinical Benefit (%)",
        "Published ORR (%)", "Published PFS (months)",
    ],
    "dose_240mg": [
        row_240["gm_auc_ss"], row_240["gm_cmax_ss"],
        row_240["orr_pct"], row_240["hepatotox_pct"],
        row_240["therapeutic_index"], row_240["nnt"], row_240["nnh"],
        row_240["net_clinical_benefit_pct"],
        24.8, 5.6,
    ],
    "dose_960mg": [
        row_960["gm_auc_ss"], row_960["gm_cmax_ss"],
        row_960["orr_pct"], row_960["hepatotox_pct"],
        row_960["therapeutic_index"], row_960["nnt"], row_960["nnh"],
        row_960["net_clinical_benefit_pct"],
        32.7, 5.4,
    ],
}
df_comparison = pd.DataFrame(comparison_data)
df_comparison.to_csv(DIR_TABLES / "benefit_risk_comparison.csv", index=False)
log.info("\n" + df_comparison.to_string(index=False))

# ── Figure 6: Dose Optimization Dashboard ────────────────────────────
fig = plt.figure(figsize=(14, 10))
gs = gridspec.GridSpec(2, 3, hspace=0.35, wspace=0.3)

# Panel A: ORR by dose
ax_a = fig.add_subplot(gs[0, 0])
orr_vals = [df_summary[df_summary["dose_mg"] == d]["orr_pct"].values[0] for d in DOSE_LEVELS]
ax_a.bar([str(d) for d in DOSE_LEVELS], orr_vals,
         color=[DOSE_COLORS[d] for d in DOSE_LEVELS], alpha=0.8,
         edgecolor="black", linewidth=0.5)
ax_a.set_ylabel("ORR (%)")
ax_a.set_xlabel("Dose (mg)")
ax_a.set_title("A. Efficacy (ORR)")
ax_a.set_ylim(0, 55)

# Panel B: Hepatotox by dose
ax_b = fig.add_subplot(gs[0, 1])
hep_vals = [df_summary[df_summary["dose_mg"] == d]["hepatotox_pct"].values[0] for d in DOSE_LEVELS]
ax_b.bar([str(d) for d in DOSE_LEVELS], hep_vals,
         color=[DOSE_COLORS[d] for d in DOSE_LEVELS], alpha=0.8,
         edgecolor="black", linewidth=0.5)
ax_b.set_ylabel("Grade 3+ Hepatotox (%)")
ax_b.set_xlabel("Dose (mg)")
ax_b.set_title("B. Safety (Hepatotoxicity)")
ax_b.set_ylim(0, 30)

# Panel C: Therapeutic Index
ax_c = fig.add_subplot(gs[0, 2])
ti_vals = [df_summary[df_summary["dose_mg"] == d]["therapeutic_index"].values[0] for d in DOSE_LEVELS]
ax_c.bar([str(d) for d in DOSE_LEVELS], ti_vals,
         color=[DOSE_COLORS[d] for d in DOSE_LEVELS], alpha=0.8,
         edgecolor="black", linewidth=0.5)
ax_c.set_ylabel("Therapeutic Index\n(ORR / Hepatotox Rate)")
ax_c.set_xlabel("Dose (mg)")
ax_c.set_title("C. Therapeutic Index")
ax_c.axhline(y=1.0, color="gray", linestyle="--", alpha=0.5)

# Panel D: Net Clinical Benefit
ax_d = fig.add_subplot(gs[1, 0])
ncb_vals = [df_summary[df_summary["dose_mg"] == d]["net_clinical_benefit_pct"].values[0] for d in DOSE_LEVELS]
ax_d.bar([str(d) for d in DOSE_LEVELS], ncb_vals,
         color=[DOSE_COLORS[d] for d in DOSE_LEVELS], alpha=0.8,
         edgecolor="black", linewidth=0.5)
ax_d.set_ylabel("Net Clinical Benefit (%)\n(ORR - Hepatotox)")
ax_d.set_xlabel("Dose (mg)")
ax_d.set_title("D. Net Clinical Benefit")
ax_d.axhline(y=0, color="gray", linestyle="--", alpha=0.5)

# Panel E: NNT vs NNH
ax_e = fig.add_subplot(gs[1, 1])
nnt_vals = [df_summary[df_summary["dose_mg"] == d]["nnt"].values[0] for d in DOSE_LEVELS]
nnh_vals = [df_summary[df_summary["dose_mg"] == d]["nnh"].values[0] for d in DOSE_LEVELS]
x = np.arange(len(DOSE_LEVELS))
w = 0.35
ax_e.bar(x - w/2, nnt_vals, w, label="NNT (efficacy)", color=COLORS["primary"], alpha=0.8)
ax_e.bar(x + w/2, nnh_vals, w, label="NNH (safety)", color=COLORS["secondary"], alpha=0.8)
ax_e.set_xticks(x)
ax_e.set_xticklabels([str(d) for d in DOSE_LEVELS])
ax_e.set_ylabel("Number of Patients")
ax_e.set_xlabel("Dose (mg)")
ax_e.set_title("E. NNT vs NNH")
ax_e.legend(loc="upper right", framealpha=0.9, fontsize=9)

# Panel F: Head-to-head summary (text)
ax_f = fig.add_subplot(gs[1, 2])
ax_f.axis("off")
text_lines = [
    "240 mg vs 960 mg",
    "─────────────────────",
    f"ORR:      {row_240['orr_pct']:.0f}% vs {row_960['orr_pct']:.0f}%",
    f"Hepatotox: {row_240['hepatotox_pct']:.0f}% vs {row_960['hepatotox_pct']:.0f}%",
    f"TI:       {row_240['therapeutic_index']:.1f} vs {row_960['therapeutic_index']:.1f}",
    f"PFS:      5.6 vs 5.4 mo (published)",
    "",
    "Interpretation:",
    "Overlapping benefit-risk;",
    "dose-optimization study",
    "justified. Subgroup",
    "effects (CPI, Cmax)",
    "may differentiate.",
]
ax_f.text(0.1, 0.95, "\n".join(text_lines), transform=ax_f.transAxes,
          fontsize=11, verticalalignment="top", fontfamily="monospace",
          bbox=dict(boxstyle="round,pad=0.5", facecolor="#F0F0F0", alpha=0.8))
ax_f.set_title("F. Summary")

fig.suptitle("Sotorasib Dose Optimization: Benefit-Risk Across Dose Levels",
             fontsize=15, y=1.01)
fig.savefig(DIR_FIGURES / "dose_optimization_dashboard.png")
plt.close(fig)
log.info("Saved dose_optimization_dashboard.png")

# ── Figure 7: Benefit-Risk Overlay ───────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 7))

from scipy.stats import binom as _binom

# Simulated points with binomial CI whiskers (ORR from N=500)
N_SIM = 500
# Offsets to repel labels (dose: (dx, dy))
_label_offsets = {180: (10, 8), 240: (-55, 10), 360: (10, -18),
                  720: (-55, -18), 960: (10, 8)}

for dose in DOSE_LEVELS:
    row = df_summary[df_summary["dose_mg"] == dose].iloc[0]
    orr_frac = row["orr_pct"] / 100
    hep_frac = row["hepatotox_pct"] / 100

    # Binomial 95% CI for ORR (y-axis)
    orr_lo = _binom.ppf(0.025, N_SIM, orr_frac) / N_SIM * 100
    orr_hi = _binom.ppf(0.975, N_SIM, orr_frac) / N_SIM * 100
    # Binomial 95% CI for hepatotox (x-axis)
    hep_lo = _binom.ppf(0.025, N_SIM, hep_frac) / N_SIM * 100
    hep_hi = _binom.ppf(0.975, N_SIM, hep_frac) / N_SIM * 100

    ax.errorbar(row["hepatotox_pct"], row["orr_pct"],
                xerr=[[row["hepatotox_pct"] - hep_lo], [hep_hi - row["hepatotox_pct"]]],
                yerr=[[row["orr_pct"] - orr_lo], [orr_hi - row["orr_pct"]]],
                fmt="none", ecolor=DOSE_COLORS[dose], capsize=4, linewidth=1.2, alpha=0.6)
    ax.scatter(row["hepatotox_pct"], row["orr_pct"],
               color=DOSE_COLORS[dose], s=250, zorder=5,
               edgecolors="black", linewidths=1.5)

    # Label only 240 and 960 prominently; others smaller
    dx, dy = _label_offsets[dose]
    if dose in (240, 960):
        ax.annotate(f"{dose} mg", (row["hepatotox_pct"], row["orr_pct"]),
                    textcoords="offset points", xytext=(dx, dy), fontsize=11,
                    fontweight="bold")
    else:
        ax.annotate(f"{dose}", (row["hepatotox_pct"], row["orr_pct"]),
                    textcoords="offset points", xytext=(dx, dy), fontsize=9,
                    color="gray")

# Published comparison points (diamonds) — only 240 and 960
ax.scatter(14.0, 24.8, marker="D", color=DOSE_COLORS[240], s=150, zorder=4,
           edgecolors="black", linewidths=1)
ax.annotate("240 mg (published)", (14.0, 24.8),
            textcoords="offset points", xytext=(-70, -18), fontsize=9,
            arrowprops=dict(arrowstyle="->", color="gray", lw=0.8))
ax.scatter(16.0, 32.7, marker="D", color=DOSE_COLORS[960], s=150, zorder=4,
           edgecolors="black", linewidths=1)
ax.annotate("960 mg (published)", (16.0, 32.7),
            textcoords="offset points", xytext=(12, -22), fontsize=9,
            arrowprops=dict(arrowstyle="->", color="gray", lw=0.8))

# Ideal region annotation
ax.annotate("Ideal:\nHigh efficacy,\nlow toxicity",
            xy=(4, 48), fontsize=10, color=COLORS["tertiary"],
            fontweight="bold", ha="center",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

# Diagonal reference lines (equal net clinical benefit)
for ncb in [0, 10, 20]:
    x_line = np.linspace(0, 30, 100)
    y_line = x_line + ncb
    ax.plot(x_line, y_line, "--", color="gray", alpha=0.3, linewidth=0.8)
    if ncb > 0:
        ax.text(28, 28 + ncb, f"NCB={ncb}%", fontsize=8, color="gray", alpha=0.6)

ax.set_xlabel("Grade 3+ Hepatotoxicity Rate (%)", fontsize=12)
ax.set_ylabel("Objective Response Rate (%)", fontsize=12)
ax.set_title("Sotorasib Benefit-Risk: Efficacy vs Safety by Dose\n"
             "(model-implied, virtual patients; published anchors = diamonds; "
             "error bars = binomial 95% CI)",
             fontsize=11)
ax.set_xlim(0, 28)
ax.set_ylim(15, 55)

from matplotlib.lines import Line2D as _L2D
_legend_els = [
    _L2D([0], [0], marker="o", color="gray", label="Simulated (virtual patients)",
         markeredgecolor="black", markersize=8, linestyle="None"),
    _L2D([0], [0], marker="D", color="gray", label="Published (trial data)",
         markeredgecolor="black", markersize=8, linestyle="None"),
]
ax.legend(handles=_legend_els, loc="lower right", framealpha=0.9, fontsize=9)

plt.tight_layout()
fig.savefig(DIR_FIGURES / "benefit_risk_overlay.png")
plt.close(fig)
log.info("Saved benefit_risk_overlay.png")

log.info("Script 04 complete.")
