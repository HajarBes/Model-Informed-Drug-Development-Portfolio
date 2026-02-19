#!/usr/bin/env python3
"""
02_exposure_efficacy.py — Exposure-Efficacy Analysis (ORR + PFS)
=================================================================
ORR Analysis:
  - Simulates binary response (responder/non-responder) per virtual patient
  - Emax logistic model: P(response) = E0 + Emax × AUC^γ / (EC50^γ + AUC^γ)
  - Calibrated to reproduce published ORR at each dose level
  - Exposure quartile analysis within 960 mg arm

PFS (Descriptive Only):
  - Reports published median PFS by dose level
  - No KM simulation (insufficient basis without IPD)

Key message: E-R appears flat because EXPOSURE is flat, not because drug is inactive.

Outputs:
  - outputs/tables/orr_by_exposure_quartile.csv
  - outputs/tables/pfs_descriptive.csv
  - outputs/tables/logistic_regression_orr.csv
  - outputs/tables/efficacy_anchor_vs_simulated.csv
  - figures/er_efficacy_panel.png
  - figures/dose_vs_exposure_vs_response.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
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
SEED_ER = _mod.SEED_ER
get_logger = _mod.get_logger
anchor_vs_simulated_table = _mod.anchor_vs_simulated_table
apply_publication_theme = _mod.apply_publication_theme

apply_publication_theme()
log = get_logger("02_efficacy")
RNG = np.random.default_rng(SEED_ER)

# ── Load Virtual Patients ────────────────────────────────────────────
df = pd.read_csv(DIR_TABLES / "virtual_patients.csv")
log.info(f"Loaded {len(df)} virtual patients")

# ── Published ORR Targets ────────────────────────────────────────────
# Using Phase 1 dose-escalation data for 180/360/720/960 and dose comparison for 240
ORR_TARGETS = {
    180: 0.33,   # Phase 1, N=3
    240: 0.248,  # Dose comparison, N=105
    360: 0.25,   # Phase 1, N=16
    720: 0.50,   # Phase 1, N=6
    960: 0.36,   # Phase 2 BICR, N=124 (use larger study)
}

# ── Emax Logistic Model for ORR ─────────────────────────────────────
# P(response) = E0 + Emax * AUC^gamma / (EC50^gamma + AUC^gamma)
#
# Because exposure is flat, we calibrate a SHALLOW E-R:
# E0 = baseline response (~0.20), Emax = plateau gain (~0.20)
# EC50 = low (well below observed AUC range) → already on plateau
# gamma = 1 (Hill coefficient)
#
# This produces ORR ~0.25-0.40 across the flat AUC range.

E0 = 0.12       # Baseline response probability
EMAX = 0.25     # Maximum additional response
EC50 = 30.0     # hr*ug/mL — well below typical AUC (~65-85)
GAMMA = 1.2     # Hill coefficient


def p_response(auc: float) -> float:
    """Probability of objective response given AUCss."""
    if auc <= 0:
        return E0
    return E0 + EMAX * auc**GAMMA / (EC50**GAMMA + auc**GAMMA)


# ── Simulate Binary Response ─────────────────────────────────────────
log.info("Simulating binary ORR using Emax logistic model...")
df["p_response"] = df["auc_ss"].apply(p_response)
df["responder"] = RNG.binomial(1, df["p_response"].values)

# Compute simulated ORR by dose
sim_orr = {}
for dose in DOSE_LEVELS:
    sub = df[df["dose_mg"] == dose]
    sim_orr[dose] = sub["responder"].mean()
    log.info(f"  {dose} mg: simulated ORR = {sim_orr[dose]:.3f} "
             f"(target = {ORR_TARGETS[dose]:.3f})")

# ── Anchor vs Simulated Table (Efficacy) ─────────────────────────────
pub_orr_pct = {d: v * 100 for d, v in ORR_TARGETS.items()}
sim_orr_pct = {d: v * 100 for d, v in sim_orr.items()}

anchor_vs_simulated_table(
    label="Efficacy Calibration — ORR",
    dose_levels=DOSE_LEVELS,
    published_values=pub_orr_pct,
    simulated_values=sim_orr_pct,
    metric_name="ORR",
    units="%",
    logger=log,
)

# Save anchor table
anchor_rows = []
for dose in DOSE_LEVELS:
    anchor_rows.append({
        "dose_mg": dose,
        "published_orr_pct": round(ORR_TARGETS[dose] * 100, 1),
        "simulated_orr_pct": round(sim_orr[dose] * 100, 1),
        "ratio": round(sim_orr[dose] / ORR_TARGETS[dose], 3)
            if ORR_TARGETS[dose] > 0 else np.nan,
    })
pd.DataFrame(anchor_rows).to_csv(DIR_TABLES / "efficacy_anchor_vs_simulated.csv", index=False)

# ── Logistic Regression: ORR ~ AUCss ────────────────────────────────
log.info("Fitting logistic regression: ORR ~ AUCss...")
try:
    import statsmodels.api as sm
    X = sm.add_constant(df["auc_ss"])
    model = sm.Logit(df["responder"], X).fit(disp=0)
    log.info(f"  AUC coefficient = {model.params.iloc[1]:.5f} "
             f"(p = {model.pvalues.iloc[1]:.4f})")
    log.info(f"  OR per 10 hr*ug/mL = {np.exp(model.params.iloc[1] * 10):.3f}")

    logistic_results = pd.DataFrame({
        "parameter": model.params.index.tolist(),
        "coefficient": model.params.values.round(5),
        "std_err": model.bse.values.round(5),
        "z_value": model.tvalues.values.round(3),
        "p_value": model.pvalues.values.round(4),
        "or": np.exp(model.params.values).round(4),
    })
    logistic_results.to_csv(DIR_TABLES / "logistic_regression_orr.csv", index=False)
    _has_statsmodels = True
except ImportError:
    log.warning("statsmodels not available; skipping logistic regression")
    _has_statsmodels = False

# ── Exposure Quartile Analysis (960 mg) ──────────────────────────────
log.info("Exposure quartile analysis within 960 mg arm...")
df_960 = df[df["dose_mg"] == 960].copy()
df_960["auc_quartile"] = pd.qcut(df_960["auc_ss"], 4, labels=["Q1", "Q2", "Q3", "Q4"])

quartile_rows = []
for q in ["Q1", "Q2", "Q3", "Q4"]:
    sub = df_960[df_960["auc_quartile"] == q]
    quartile_rows.append({
        "quartile": q,
        "n": len(sub),
        "auc_median": round(sub["auc_ss"].median(), 1),
        "auc_range": f"{sub['auc_ss'].min():.1f}-{sub['auc_ss'].max():.1f}",
        "orr_pct": round(sub["responder"].mean() * 100, 1),
        "n_responders": int(sub["responder"].sum()),
    })
df_quartile = pd.DataFrame(quartile_rows)
df_quartile.to_csv(DIR_TABLES / "orr_by_exposure_quartile.csv", index=False)
log.info("  Quartile ORR:")
for _, row in df_quartile.iterrows():
    log.info(f"    {row['quartile']}: ORR = {row['orr_pct']}% "
             f"(median AUC = {row['auc_median']})")

# ── PFS Descriptive Table ────────────────────────────────────────────
# Published medians only — no simulated time-to-event
pfs_data = [
    {"dose_mg": 960, "median_pfs_months": 5.4, "source": "Dose comparison (EJC 2024)", "n": 104},
    {"dose_mg": 240, "median_pfs_months": 5.6, "source": "Dose comparison (EJC 2024)", "n": 105},
]
# Implied hazard ratio
hr_implied = (np.log(2) / 5.4) / (np.log(2) / 5.6)
pfs_data[0]["implied_hr_vs_240mg"] = round(hr_implied, 3)
pfs_data[1]["implied_hr_vs_240mg"] = 1.0
df_pfs = pd.DataFrame(pfs_data)
df_pfs.to_csv(DIR_TABLES / "pfs_descriptive.csv", index=False)
log.info(f"  PFS: 960 mg = 5.4 mo, 240 mg = 5.6 mo, implied HR = {hr_implied:.3f}")

# ── Figure 3: E-R Efficacy Panel (3 panels) ─────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Panel A: Dose → ORR with N labels and binomial 95% CI
from scipy.stats import binom

# Published sample sizes for CI computation
N_PUBLISHED = {180: 3, 240: 105, 360: 16, 720: 6, 960: 124}

doses_plot = sorted(ORR_TARGETS.keys())
orr_pub = [ORR_TARGETS[d] * 100 for d in doses_plot]
orr_sim = [sim_orr[d] * 100 for d in doses_plot]
x = np.arange(len(doses_plot))
w = 0.35

# Binomial exact 95% CI for published ORR
pub_ci_lo, pub_ci_hi = [], []
for d in doses_plot:
    n_d = N_PUBLISHED[d]
    k_d = round(ORR_TARGETS[d] * n_d)
    lo = binom.ppf(0.025, n_d, ORR_TARGETS[d]) / n_d * 100 if n_d > 0 else 0
    hi = binom.ppf(0.975, n_d, ORR_TARGETS[d]) / n_d * 100 if n_d > 0 else 0
    pub_ci_lo.append(lo)
    pub_ci_hi.append(hi)
pub_yerr_lo = [o - lo for o, lo in zip(orr_pub, pub_ci_lo)]
pub_yerr_hi = [hi - o for o, hi in zip(orr_pub, pub_ci_hi)]

bars_pub = axes[0].bar(x - w/2, orr_pub, w, label="Published", color=COLORS["primary"], alpha=0.8)
axes[0].errorbar(x - w/2, orr_pub, yerr=[pub_yerr_lo, pub_yerr_hi],
                 fmt="none", ecolor="black", capsize=3, linewidth=1)
axes[0].bar(x + w/2, orr_sim, w, label="Simulated", color=COLORS["secondary"], alpha=0.8)

# Add N labels above published bars
for i, d in enumerate(doses_plot):
    axes[0].text(i - w/2, orr_pub[i] + pub_yerr_hi[i] + 1,
                 f"N={N_PUBLISHED[d]}", ha="center", va="bottom", fontsize=8, color="gray")

axes[0].set_xticks(x)
axes[0].set_xticklabels([str(d) for d in doses_plot])
axes[0].set_xlabel("Dose (mg)")
axes[0].set_ylabel("ORR (%)")
axes[0].set_title("A. Dose vs ORR\n(error bars: binomial 95% CI)")
axes[0].legend(loc="upper left", framealpha=0.9)
axes[0].set_ylim(0, 85)

# Panel B: AUC quartile → ORR (960 mg, model-implied)
q_labels = df_quartile["quartile"].tolist()
q_orr = df_quartile["orr_pct"].tolist()
q_colors = [COLORS["primary"], COLORS["tertiary"], COLORS["quaternary"], COLORS["secondary"]]
bars = axes[1].bar(q_labels, q_orr, color=q_colors, alpha=0.8, edgecolor="black", linewidth=0.5)
axes[1].set_xlabel("AUCss Quartile (960 mg)")
axes[1].set_ylabel("ORR (%)")
axes[1].set_title("B. ORR by AUC Quartile — 960 mg\n(model-implied, virtual patients)")
axes[1].set_ylim(0, 65)
# Annotate counts
for bar, row in zip(bars, df_quartile.itertuples()):
    axes[1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                 f"n={row.n}", ha="center", va="bottom", fontsize=9)

# Panel C: PFS by dose (descriptive, published medians)
pfs_doses = [240, 960]
pfs_vals = [5.6, 5.4]
pfs_colors = [DOSE_COLORS[240], DOSE_COLORS[960]]
axes[2].bar([str(d) for d in pfs_doses], pfs_vals, color=pfs_colors, alpha=0.8,
            edgecolor="black", linewidth=0.5)
axes[2].set_xlabel("Dose (mg)")
axes[2].set_ylabel("Median PFS (months)")
axes[2].set_title("C. Median PFS by Dose\n(published, descriptive)")
axes[2].set_ylim(0, 8)
for i, (d, v) in enumerate(zip(pfs_doses, pfs_vals)):
    axes[2].text(i, v + 0.1, f"{v} mo", ha="center", va="bottom", fontsize=10)
axes[2].axhline(y=5.5, color=COLORS["light_gray"], linestyle="--", alpha=0.5)

fig.suptitle("Sotorasib Exposure-Efficacy: Flat E-R Driven by Flat Exposure",
             fontsize=14, y=1.02)
plt.tight_layout()
fig.savefig(DIR_FIGURES / "er_efficacy_panel.png")
plt.close(fig)
log.info("Saved er_efficacy_panel.png")

# ── Figure 4: Dose → Exposure → Response Chain ──────────────────────
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Panel A: Dose → AUC (flat)
gm_auc_by_dose = df.groupby("dose_mg")["auc_ss"].apply(
    lambda x: np.exp(np.mean(np.log(x)))).to_dict()
axes[0].plot(DOSE_LEVELS, [gm_auc_by_dose[d] for d in DOSE_LEVELS],
             "o-", color=COLORS["primary"], markersize=10, linewidth=2)
axes[0].fill_between(DOSE_LEVELS,
                     [df[df["dose_mg"]==d]["auc_ss"].quantile(0.25) for d in DOSE_LEVELS],
                     [df[df["dose_mg"]==d]["auc_ss"].quantile(0.75) for d in DOSE_LEVELS],
                     alpha=0.2, color=COLORS["primary"])
axes[0].set_xlabel("Dose (mg)")
axes[0].set_ylabel("AUCss (hr*ug/mL)")
axes[0].set_title("A. Dose vs Exposure\n(FLAT: saturable absorption)")
axes[0].annotate("Saturable\nabsorption", xy=(600, 70), fontsize=11,
                 color=COLORS["secondary"], fontweight="bold")

# Panel B: AUC → ORR (shallow)
auc_grid = np.linspace(5, 250, 200)
p_grid = [p_response(a) * 100 for a in auc_grid]
axes[1].plot(auc_grid, p_grid, color=COLORS["secondary"], linewidth=2)
# Overlay dose-level points
for dose in DOSE_LEVELS:
    axes[1].scatter(gm_auc_by_dose[dose], sim_orr[dose] * 100,
                    color=DOSE_COLORS[dose], s=120, zorder=5,
                    edgecolors="black", linewidths=0.8)
    axes[1].annotate(f"{dose} mg", (gm_auc_by_dose[dose], sim_orr[dose]*100),
                     textcoords="offset points", xytext=(8, 5), fontsize=9)
# Show the narrow AUC window
auc_min = min(gm_auc_by_dose.values())
auc_max = max(gm_auc_by_dose.values())
axes[1].axvspan(auc_min - 5, auc_max + 5, alpha=0.1, color=COLORS["quaternary"])
axes[1].annotate("Observed AUC\nrange (all doses)", xy=(75, 15), fontsize=9,
                 color=COLORS["dark"], ha="center")
axes[1].set_xlabel("AUCss (hr*ug/mL)")
axes[1].set_ylabel("ORR (%)")
axes[1].set_title("B. Exposure vs Response\n(shallow E-R on plateau)")

# Panel C: Dose → ORR (the combination = flat)
axes[2].bar([str(d) for d in DOSE_LEVELS],
            [sim_orr[d] * 100 for d in DOSE_LEVELS],
            color=[DOSE_COLORS[d] for d in DOSE_LEVELS],
            alpha=0.8, edgecolor="black", linewidth=0.5)
axes[2].set_xlabel("Dose (mg)")
axes[2].set_ylabel("ORR (%)")
axes[2].set_title("C. Dose vs Response\n(flat: exposure doesn't change)")
axes[2].set_ylim(0, 60)

fig.suptitle("Sotorasib: Why Dose-Response Is Flat (Exposure Is the Missing Link)",
             fontsize=14, y=1.02)
plt.tight_layout()
fig.savefig(DIR_FIGURES / "dose_vs_exposure_vs_response.png")
plt.close(fig)
log.info("Saved dose_vs_exposure_vs_response.png")

log.info("Script 02 complete.")
