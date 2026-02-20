#!/usr/bin/env python3
"""
KRAS G12C QSP Model — Integrated PK → Target Engagement → Signaling → Tumor

Modules:
  M1: Oral 1-compartment PK (sotorasib 960 mg QD)
  M2: KRAS GDP/GTP cycling + covalent drug binding
  M3: Receptor feedback (R) + pathway output (E, ERK proxy)
  M4: Tumor growth driven by E(t)
  M5: CRC vs NSCLC parameterization (feedback gain G_fb)
  M6: Combination therapy (EGFR blockade reduces receptor drive)
  M6b: Adaptive resistance state Z(t) + saturating kill

State variables (9):
  y = [A_gut, A_central, KRAS_GDP, KRAS_GTP, KRAS_drug, R, E, T, Z]

Author: Hajar Besbassi
"""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


# ─────────────────────────────────────────
# Parameters
# ─────────────────────────────────────────

@dataclass
class Params:
    """All model parameters with literature-based defaults."""

    # ── PK (M1) ──
    dose_mg: float = 960.0        # sotorasib dose
    tau_h: float = 24.0           # dosing interval (QD)
    n_doses: int = 999            # effectively unlimited for long sims
    ka: float = 1.0               # absorption rate (/h)
    t_half_h: float = 5.0         # elimination half-life (h)
    V: float = 211.0              # Vss/F (L), FDA label
    F: float = 0.73               # bioavailability, FDA label
    MW: float = 560.6             # molecular weight (g/mol)

    # ── KRAS cycling (M2) ──
    k_hyd: float = 0.2            # GTPase hydrolysis rate (/h)
    k_GEF_basal: float = 0.025    # basal GEF exchange rate (/h)
    k_bind: float = 0.15          # covalent drug binding rate (/(µM·h)) — tuned for ~60-80% TE at Cmax
    k_deg: float = 0.029          # protein degradation (/h), t1/2 ~ 24h
    KRAS_total_ss: float = 1.0    # total KRAS at steady state (µM)

    # ── Feedback / signaling (M3) ──
    tau_R: float = 36.0           # receptor drive time constant (h) — RTK upregulation ~1.5 days
    tau_E: float = 6.0            # pathway output time constant (h) — ERK signaling
    R_basal: float = 1.0          # basal receptor drive (normalized)
    G_fb: float = 1.2             # feedback gain (NSCLC default; CRC = 3.5)
    E_ss: float = 1.0             # steady-state pathway output (normalized)

    # ── Tumor (M4) ──
    rho_grow: float = 0.0005      # baseline growth rate (/h) ~ doubling ~58 days (clinical timescale)
    rho_kill: float = 0.0045      # max kill rate (/h) — calibrated for ORR with saturating kill + resistance
    n_hill: float = 2.0           # Hill coefficient for E-driven growth
    T_max: float = 5.0            # carrying capacity (normalized)
    T_0: float = 1.0              # initial tumor volume (normalized)

    # ── Combination (M6) ──
    f_block: float = 0.0          # EGFR blockade fraction [0,1]. 0=mono, 0.55=combo

    # ── Resistance / saturation (M6b) ──
    k_Z_up: float = 0.004         # resistance induction rate (/h) when E is suppressed
    k_Z_down: float = 0.0002      # resistance decay rate (/h) — slow reversibility
    alpha_Z: float = 4.0          # resistance strength: kill / (1 + alpha_Z * Z)
    beta_Z: float = 0.8           # bypass signaling: Z restores E independently of KRAS-GTP
    EC50_kill: float = 0.25       # half-max E suppression for saturating kill (dimensionless)
    Z_0: float = 0.0              # initial resistance level [0,1]; >0 for pre-existing resistance

    @property
    def ke(self) -> float:
        return np.log(2) / self.t_half_h

    @property
    def k_syn(self) -> float:
        """Protein synthesis rate to maintain steady-state KRAS."""
        return self.k_deg * self.KRAS_total_ss

    def kras_ss(self) -> Tuple[float, float, float]:
        """Compute KRAS steady-state (no drug): GDP, GTP, drug=0."""
        # At SS from dKRAS_GTP=0: k_GEF * R * GDP = (k_hyd + k_deg) * GTP
        # So GTP/GDP = k_GEF * R / (k_hyd + k_deg)
        # And GDP + GTP = KRAS_total (from synthesis/degradation balance)
        ratio = self.k_GEF_basal * self.R_basal / (self.k_hyd + self.k_deg)
        GDP_ss = self.KRAS_total_ss / (1.0 + ratio)
        GTP_ss = self.KRAS_total_ss - GDP_ss
        return GDP_ss, GTP_ss, 0.0


# ─────────────────────────────────────────
# Preset scenarios
# ─────────────────────────────────────────

def params_nsclc_mono() -> Params:
    return Params(G_fb=1.2, f_block=0.0)

def params_crc_mono() -> Params:
    return Params(G_fb=3.5, f_block=0.0, rho_grow=0.0007)

def params_crc_combo() -> Params:
    return Params(G_fb=3.5, f_block=0.55, rho_grow=0.0007)

def params_vehicle() -> Params:
    """No drug — just tumor growth."""
    p = Params(G_fb=0.5, f_block=0.0)
    p.dose_mg = 0.0
    return p


SCENARIOS = {
    "vehicle": params_vehicle,
    "nsclc_mono": params_nsclc_mono,
    "crc_mono": params_crc_mono,
    "crc_combo": params_crc_combo,
}


# ─────────────────────────────────────────
# ODE system
# ─────────────────────────────────────────

def rhs(t: float, y: np.ndarray, p: Params, dose_times: np.ndarray,
        KRAS_GTP_ss: float) -> np.ndarray:
    """Right-hand side of the integrated ODE system (9 states)."""
    A_gut, A_cent, KRAS_GDP, KRAS_GTP, KRAS_drug, R, E, T, Z = y

    # Clamp non-negative
    A_gut = max(A_gut, 0.0)
    A_cent = max(A_cent, 0.0)
    KRAS_GDP = max(KRAS_GDP, 1e-12)
    KRAS_GTP = max(KRAS_GTP, 1e-12)
    KRAS_drug = max(KRAS_drug, 0.0)
    R = max(R, 1e-6)
    E = max(E, 1e-6)
    T = max(T, 1e-6)
    Z = max(min(Z, 1.0), 0.0)  # clamp Z ∈ [0, 1]

    # ── PK ──
    dA_gut = -p.ka * A_gut
    dA_cent = p.ka * A_gut - p.ke * A_cent
    C_mgL = A_cent / p.V
    C_uM = C_mgL * 1000.0 / p.MW  # mg/L → µM

    # ── KRAS cycling (M2) ──
    exchange = p.k_GEF_basal * R * KRAS_GDP
    hydrolysis = p.k_hyd * KRAS_GTP
    binding = p.k_bind * C_uM * KRAS_GDP

    dKRAS_GDP = hydrolysis - exchange - binding + p.k_syn - p.k_deg * KRAS_GDP
    dKRAS_GTP = exchange - hydrolysis - p.k_deg * KRAS_GTP
    dKRAS_drug = binding - p.k_deg * KRAS_drug

    # ── Feedback / signaling (M3) ──
    E_ratio = E / p.E_ss
    feedback_release = p.G_fb * max(0.0, 1.0 - E_ratio)
    R_target = (p.R_basal + feedback_release) * (1.0 - p.f_block)
    dR = (1.0 / p.tau_R) * (R_target - R)

    KRAS_GTP_ratio = KRAS_GTP / max(KRAS_GTP_ss, 1e-12)
    # Bypass signaling: Z restores pathway output via alternative pathways
    E_drive = R * KRAS_GTP_ratio + p.beta_Z * Z
    dE = (1.0 / p.tau_E) * (E_drive - E)

    # ── Resistance program (M6b) ──
    # Z accumulates when drug is bound to target (TE > 0 → selective pressure)
    # Decays slowly — partially irreversible on clinical timescales
    KRAS_total = KRAS_GDP + KRAS_GTP + KRAS_drug
    TE_local = KRAS_drug / max(KRAS_total, 1e-12)  # target engagement [0,1]
    dZ = p.k_Z_up * TE_local * (1.0 - Z) - p.k_Z_down * Z

    # ── Tumor (M4) with saturating kill + resistance ──
    growth = p.rho_grow * (E_ratio ** p.n_hill) * T * (1.0 - T / p.T_max)

    # Saturating kill: drug effect saturates at high E suppression
    # Resistance Z reduces kill effectiveness
    E_drug = max(0.0, 1.0 - E_ratio)  # drug-induced suppression [0,1]
    kill_sat = E_drug / (p.EC50_kill + E_drug)  # saturating drug effect
    kill_resist = 1.0 / (1.0 + p.alpha_Z * Z)  # resistance attenuation
    kill = p.rho_kill * kill_sat * kill_resist * T
    dT = growth - kill

    return np.array([dA_gut, dA_cent, dKRAS_GDP, dKRAS_GTP, dKRAS_drug, dR, dE, dT, dZ])


# ─────────────────────────────────────────
# Simulation engine
# ─────────────────────────────────────────

def simulate(p: Params, t_days: int = 180, dt_h: float = 0.5) -> pd.DataFrame:
    """
    Simulate the full QSP model with repeated dosing.

    Uses solve_ivp between dose events for accuracy, with bolus additions at dose times.
    """
    t_end_h = t_days * 24.0
    dose_times_h = np.arange(0.0, t_end_h, p.tau_h)

    # Initial conditions
    GDP_ss, GTP_ss, _ = p.kras_ss()
    y0 = np.array([
        0.0,        # A_gut
        0.0,        # A_central
        GDP_ss,     # KRAS_GDP
        GTP_ss,     # KRAS_GTP
        0.0,        # KRAS_drug
        p.R_basal * (1.0 - p.f_block),  # R (adjusted for combo)
        p.E_ss,     # E
        p.T_0,      # T
        p.Z_0,      # Z (resistance program)
    ])

    # Collect results
    all_t = []
    all_y = []

    y_current = y0.copy()
    t_current = 0.0

    # Create time segments between doses
    event_times = list(dose_times_h)
    if event_times[-1] < t_end_h:
        event_times.append(t_end_h)

    for i in range(len(event_times)):
        t_dose = event_times[i]

        # Apply dose at this time
        if t_dose in dose_times_h or (i < len(dose_times_h) and abs(t_dose - dose_times_h[i]) < 0.01):
            if i < len(dose_times_h):
                y_current[0] += p.F * p.dose_mg  # bolus into gut

        # Determine next event
        if i + 1 < len(event_times):
            t_next = event_times[i + 1]
        else:
            t_next = t_end_h

        if t_next <= t_dose:
            continue

        # Solve ODE from t_dose to t_next
        t_span = (t_dose, t_next)
        n_pts = max(int((t_next - t_dose) / dt_h), 2)
        t_eval = np.linspace(t_dose, t_next, n_pts, endpoint=False)

        sol = solve_ivp(
            lambda t, y: rhs(t, y, p, dose_times_h, GTP_ss),
            t_span,
            y_current,
            t_eval=t_eval,
            method="RK45",
            rtol=1e-6,
            atol=1e-9,
            max_step=1.0,
        )

        if not sol.success:
            print(f"  Warning: solver failed at t={t_dose:.1f}h: {sol.message}")
            break

        all_t.append(sol.t)
        all_y.append(sol.y)

        # Update state for next segment
        y_current = sol.y[:, -1].copy()

    # Concatenate
    t_arr = np.concatenate(all_t)
    y_arr = np.concatenate(all_y, axis=1)

    # Derived quantities
    C_mgL = y_arr[1, :] / p.V
    C_uM = C_mgL * 1000.0 / p.MW
    KRAS_total = y_arr[2, :] + y_arr[3, :] + y_arr[4, :]
    TE = np.where(KRAS_total > 1e-12, y_arr[4, :] / KRAS_total, 0.0)
    KRAS_GTP_norm = y_arr[3, :] / max(GTP_ss, 1e-12)
    E_norm = y_arr[6, :] / p.E_ss
    T_change_pct = (y_arr[7, :] / p.T_0 - 1.0) * 100.0

    df = pd.DataFrame({
        "time_h": t_arr,
        "time_d": t_arr / 24.0,
        "A_gut": y_arr[0, :],
        "A_central": y_arr[1, :],
        "C_mgL": C_mgL,
        "C_uM": C_uM,
        "KRAS_GDP": y_arr[2, :],
        "KRAS_GTP": y_arr[3, :],
        "KRAS_drug": y_arr[4, :],
        "TE": TE,
        "KRAS_GTP_norm": KRAS_GTP_norm,
        "R": y_arr[5, :],
        "E": y_arr[6, :],
        "E_norm": E_norm,
        "T": y_arr[7, :],
        "T_change_pct": T_change_pct,
        "Z": y_arr[8, :],
    })
    return df


# ─────────────────────────────────────────
# Plotting
# ─────────────────────────────────────────

COLORS = {
    "vehicle": "#888888",
    "nsclc_mono": "#377eb8",
    "crc_mono": "#e41a1c",
    "crc_combo": "#4daf4a",
}

LABELS = {
    "vehicle": "Vehicle (no drug)",
    "nsclc_mono": "NSCLC — sotorasib mono",
    "crc_mono": "CRC — sotorasib mono",
    "crc_combo": "CRC — sotorasib + panitumumab",
}


def fig1_pk_profile(results: Dict[str, pd.DataFrame], outdir: str) -> None:
    """Fig 1: PK concentration profile (first 7 days)."""
    df = results["nsclc_mono"] if "nsclc_mono" in results else next(iter(results.values()))
    sub = df[df["time_d"] <= 7.0]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True,
                                    gridspec_kw={"height_ratios": [2, 1], "hspace": 0.08})

    ax1.plot(sub["time_d"], sub["C_uM"], color="#377eb8", linewidth=1.5)
    ax1.set_ylabel("Concentration (µM)", fontsize=10)
    ax1.set_title("Sotorasib 960 mg QD — Plasma Concentration", fontsize=12, fontweight="bold")
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # Dose markers
    for d in range(7):
        ax1.axvline(d, color="#cccccc", linewidth=0.5, linestyle="--", alpha=0.5)

    ax2.plot(sub["time_d"], sub["C_mgL"], color="#e41a1c", linewidth=1.5)
    ax2.set_ylabel("Concentration (mg/L)", fontsize=10)
    ax2.set_xlabel("Time (days)", fontsize=10)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig1_pk_profile.png"), dpi=200)
    plt.close(fig)
    print("  Saved fig1_pk_profile.png")


def fig2_target_engagement(results: Dict[str, pd.DataFrame], outdir: str) -> None:
    """Fig 2: Target engagement and KRAS_GTP suppression (7 days)."""
    df = results["nsclc_mono"] if "nsclc_mono" in results else next(iter(results.values()))
    sub = df[df["time_d"] <= 7.0]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True,
                                    gridspec_kw={"height_ratios": [1, 1], "hspace": 0.12})

    ax1.plot(sub["time_d"], sub["TE"] * 100, color="#984ea3", linewidth=1.5)
    ax1.set_ylabel("Target Engagement (%)", fontsize=10)
    ax1.set_title("KRAS G12C Target Engagement & GTP Suppression", fontsize=12, fontweight="bold")
    ax1.set_ylim(-5, 105)
    ax1.axhline(80, color="#999", linewidth=0.5, linestyle="--", alpha=0.5)
    ax1.text(0.1, 82, "80% TE", fontsize=8, color="#666")
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    ax2.plot(sub["time_d"], sub["KRAS_GTP_norm"] * 100, color="#ff7f00", linewidth=1.5)
    ax2.set_ylabel("KRAS-GTP (% of baseline)", fontsize=10)
    ax2.set_xlabel("Time (days)", fontsize=10)
    ax2.axhline(100, color="#999", linewidth=0.5, linestyle="--", alpha=0.5)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig2_target_engagement.png"), dpi=200)
    plt.close(fig)
    print("  Saved fig2_target_engagement.png")


def fig3_adaptive_rebound(results: Dict[str, pd.DataFrame], outdir: str) -> None:
    """Fig 3: ERK proxy rebound — NSCLC vs CRC."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    # Panel A: E(t) — pathway output
    ax = axes[0]
    for key in ["nsclc_mono", "crc_mono"]:
        if key in results:
            df = results[key]
            sub = df[df["time_d"] <= 14.0]
            ax.plot(sub["time_d"], sub["E_norm"] * 100, color=COLORS[key],
                    linewidth=1.5, label=LABELS[key])
    ax.axhline(100, color="#999", linewidth=0.5, linestyle="--")
    ax.axhline(75, color="#999", linewidth=0.5, linestyle=":", alpha=0.5)
    ax.text(13.5, 77, "75%", fontsize=8, color="#666", ha="right")
    ax.set_xlabel("Time (days)", fontsize=10)
    ax.set_ylabel("Pathway output E (% of baseline)", fontsize=10)
    ax.set_title("A. ERK Rebound", fontsize=11, fontweight="bold")
    ax.legend(fontsize=8, frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Panel B: R(t) — receptor drive
    ax = axes[1]
    for key in ["nsclc_mono", "crc_mono"]:
        if key in results:
            df = results[key]
            sub = df[df["time_d"] <= 14.0]
            ax.plot(sub["time_d"], sub["R"], color=COLORS[key],
                    linewidth=1.5, label=LABELS[key])
    ax.axhline(1.0, color="#999", linewidth=0.5, linestyle="--")
    ax.set_xlabel("Time (days)", fontsize=10)
    ax.set_ylabel("Receptor drive R", fontsize=10)
    ax.set_title("B. Receptor Feedback", fontsize=11, fontweight="bold")
    ax.legend(fontsize=8, frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Panel C: KRAS_GTP normalized
    ax = axes[2]
    for key in ["nsclc_mono", "crc_mono"]:
        if key in results:
            df = results[key]
            sub = df[df["time_d"] <= 14.0]
            ax.plot(sub["time_d"], sub["KRAS_GTP_norm"] * 100, color=COLORS[key],
                    linewidth=1.5, label=LABELS[key])
    ax.axhline(100, color="#999", linewidth=0.5, linestyle="--")
    ax.set_xlabel("Time (days)", fontsize=10)
    ax.set_ylabel("KRAS-GTP (% baseline)", fontsize=10)
    ax.set_title("C. KRAS-GTP Suppression", fontsize=11, fontweight="bold")
    ax.legend(fontsize=8, frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.suptitle("Adaptive Resistance: NSCLC vs CRC", fontsize=13, fontweight="bold", y=1.03)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig3_adaptive_rebound.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("  Saved fig3_adaptive_rebound.png")


def fig4_tumor_growth(results: Dict[str, pd.DataFrame], outdir: str) -> None:
    """Fig 4: Tumor growth curves — vehicle, NSCLC mono, CRC mono."""
    fig, ax = plt.subplots(figsize=(8, 5))

    for key in ["vehicle", "nsclc_mono", "crc_mono"]:
        if key in results:
            df = results[key]
            sub = df[df["time_d"] <= 180]
            ax.plot(sub["time_d"], sub["T"], color=COLORS[key],
                    linewidth=2, label=LABELS[key])

    ax.axhline(results.get("nsclc_mono", results.get("crc_mono", pd.DataFrame({"T_0": [10]}))).iloc[0].get("T", 10),
               color="#999", linewidth=0.5, linestyle="--", alpha=0.3)
    ax.set_xlabel("Time (days)", fontsize=11)
    ax.set_ylabel("Tumor volume (cm³)", fontsize=11)
    ax.set_title("Tumor Growth: Vehicle vs Monotherapy", fontsize=13, fontweight="bold")
    ax.legend(fontsize=9, frameon=False, loc="upper left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig4_tumor_growth.png"), dpi=200)
    plt.close(fig)
    print("  Saved fig4_tumor_growth.png")


def fig5_crc_vs_nsclc(results: Dict[str, pd.DataFrame], outdir: str) -> None:
    """Fig 5: Side-by-side CRC vs NSCLC comparison."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel A: Tumor change (%)
    ax = axes[0]
    for key in ["nsclc_mono", "crc_mono"]:
        if key in results:
            df = results[key]
            sub = df[df["time_d"] <= 180]
            ax.plot(sub["time_d"], sub["T_change_pct"], color=COLORS[key],
                    linewidth=2, label=LABELS[key])
    ax.axhline(0, color="#999", linewidth=0.5, linestyle="--")
    ax.axhline(-30, color="#666", linewidth=0.5, linestyle=":", alpha=0.5)
    ax.text(175, -28, "RECIST PR (-30%)", fontsize=8, color="#666", ha="right")
    ax.set_xlabel("Time (days)", fontsize=10)
    ax.set_ylabel("Tumor change from baseline (%)", fontsize=10)
    ax.set_title("A. Tumor Response", fontsize=11, fontweight="bold")
    ax.legend(fontsize=9, frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Panel B: ERK rebound overlay
    ax = axes[1]
    for key in ["nsclc_mono", "crc_mono"]:
        if key in results:
            df = results[key]
            sub = df[df["time_d"] <= 30]
            ax.plot(sub["time_d"], sub["E_norm"] * 100, color=COLORS[key],
                    linewidth=2, label=LABELS[key])
    ax.axhline(100, color="#999", linewidth=0.5, linestyle="--")
    ax.set_xlabel("Time (days)", fontsize=10)
    ax.set_ylabel("Pathway output E (% baseline)", fontsize=10)
    ax.set_title("B. Pathway Rebound (first 30 days)", fontsize=11, fontweight="bold")
    ax.legend(fontsize=9, frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.suptitle("CRC vs NSCLC: Differential Response to KRAS G12C Inhibition",
                 fontsize=13, fontweight="bold", y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig5_crc_vs_nsclc.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("  Saved fig5_crc_vs_nsclc.png")


def fig6_combination_benefit(results: Dict[str, pd.DataFrame], outdir: str) -> None:
    """Fig 6: CRC mono vs combo — tumor + rebound suppression."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel A: Tumor
    ax = axes[0]
    for key in ["crc_mono", "crc_combo"]:
        if key in results:
            df = results[key]
            sub = df[df["time_d"] <= 180]
            ax.plot(sub["time_d"], sub["T_change_pct"], color=COLORS[key],
                    linewidth=2, label=LABELS[key])
    ax.axhline(0, color="#999", linewidth=0.5, linestyle="--")
    ax.axhline(-30, color="#666", linewidth=0.5, linestyle=":", alpha=0.5)
    ax.text(175, -28, "RECIST PR (-30%)", fontsize=8, color="#666", ha="right")
    ax.set_xlabel("Time (days)", fontsize=10)
    ax.set_ylabel("Tumor change from baseline (%)", fontsize=10)
    ax.set_title("A. Tumor Response: Mono vs Combo", fontsize=11, fontweight="bold")
    ax.legend(fontsize=9, frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Panel B: Rebound suppression
    ax = axes[1]
    for key in ["crc_mono", "crc_combo"]:
        if key in results:
            df = results[key]
            sub = df[df["time_d"] <= 30]
            ax.plot(sub["time_d"], sub["E_norm"] * 100, color=COLORS[key],
                    linewidth=2, label=LABELS[key])
    ax.axhline(100, color="#999", linewidth=0.5, linestyle="--")
    ax.set_xlabel("Time (days)", fontsize=10)
    ax.set_ylabel("Pathway output E (% baseline)", fontsize=10)
    ax.set_title("B. Rebound Suppression by EGFR Blockade", fontsize=11, fontweight="bold")
    ax.legend(fontsize=9, frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.suptitle("Combination Benefit: Sotorasib + Panitumumab in CRC",
                 fontsize=13, fontweight="bold", y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig6_combination_benefit.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("  Saved fig6_combination_benefit.png")


def fig7_model_summary(results: Dict[str, pd.DataFrame], outdir: str) -> None:
    """Fig 7: 4-panel summary dashboard."""
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    # A: PK (7 days)
    ax = axes[0, 0]
    df = results["nsclc_mono"] if "nsclc_mono" in results else next(iter(results.values()))
    sub = df[df["time_d"] <= 7]
    ax.plot(sub["time_d"], sub["C_uM"], color="#377eb8", linewidth=1.5)
    ax.set_xlabel("Time (days)", fontsize=9)
    ax.set_ylabel("Concentration (µM)", fontsize=9)
    ax.set_title("A. Sotorasib PK (960 mg QD)", fontsize=10, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # B: Target Engagement (7 days)
    ax = axes[0, 1]
    sub = df[df["time_d"] <= 7]
    ax.plot(sub["time_d"], sub["TE"] * 100, color="#984ea3", linewidth=1.5)
    ax.set_xlabel("Time (days)", fontsize=9)
    ax.set_ylabel("Target Engagement (%)", fontsize=9)
    ax.set_title("B. KRAS G12C Target Engagement", fontsize=10, fontweight="bold")
    ax.set_ylim(-5, 105)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # C: ERK rebound (30 days, all scenarios)
    ax = axes[1, 0]
    for key in ["nsclc_mono", "crc_mono", "crc_combo"]:
        if key in results:
            df2 = results[key]
            sub2 = df2[df2["time_d"] <= 30]
            ax.plot(sub2["time_d"], sub2["E_norm"] * 100, color=COLORS[key],
                    linewidth=1.5, label=LABELS[key])
    ax.axhline(100, color="#999", linewidth=0.5, linestyle="--")
    ax.set_xlabel("Time (days)", fontsize=9)
    ax.set_ylabel("Pathway output (% baseline)", fontsize=9)
    ax.set_title("C. Adaptive Rebound", fontsize=10, fontweight="bold")
    ax.legend(fontsize=7, frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # D: Tumor curves (180 days, all scenarios)
    ax = axes[1, 1]
    for key in ["vehicle", "nsclc_mono", "crc_mono", "crc_combo"]:
        if key in results:
            df2 = results[key]
            sub2 = df2[df2["time_d"] <= 180]
            ax.plot(sub2["time_d"], sub2["T_change_pct"], color=COLORS[key],
                    linewidth=1.5, label=LABELS[key])
    ax.axhline(0, color="#999", linewidth=0.5, linestyle="--")
    ax.axhline(-30, color="#666", linewidth=0.5, linestyle=":", alpha=0.3)
    ax.set_xlabel("Time (days)", fontsize=9)
    ax.set_ylabel("Tumor change (%)", fontsize=9)
    ax.set_title("D. Tumor Response", fontsize=10, fontweight="bold")
    ax.legend(fontsize=7, frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.suptitle("KRAS G12C QSP Model — Summary Dashboard", fontsize=14, fontweight="bold", y=1.01)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig7_model_summary.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("  Saved fig7_model_summary.png")


# ─────────────────────────────────────────
# Validation metrics
# ─────────────────────────────────────────

def print_validation(results: Dict[str, pd.DataFrame]) -> None:
    """Print key validation metrics."""
    print("\n" + "="*60)
    print("VALIDATION METRICS")
    print("="*60)

    for name, df in results.items():
        print(f"\n--- {LABELS.get(name, name)} ---")

        # PK metrics (at steady state, day 5+)
        ss = df[df["time_d"] >= 5.0]
        if len(ss) > 0 and "C_uM" in ss.columns:
            cmax = ss["C_uM"].max()
            cmin = ss["C_uM"].min()
            print(f"  PK: Cmax = {cmax:.2f} µM, Ctrough = {cmin:.4f} µM")

        # TE metrics
        if "TE" in df.columns and len(ss) > 0:
            te_max = ss["TE"].max() * 100
            te_min = ss["TE"].min() * 100
            print(f"  TE: max = {te_max:.1f}%, min = {te_min:.1f}%")

        # ERK rebound at 72h
        t72 = df[(df["time_h"] >= 70) & (df["time_h"] <= 74)]
        if len(t72) > 0:
            e72 = t72["E_norm"].mean() * 100
            print(f"  ERK at 72h: {e72:.1f}% of baseline")

        # Tumor metrics
        if "T_change_pct" in df.columns:
            min_change = df["T_change_pct"].min()
            t_nadir = df.loc[df["T_change_pct"].idxmin(), "time_d"]
            final_change = df.iloc[-1]["T_change_pct"]
            print(f"  Tumor: best change = {min_change:.1f}% at day {t_nadir:.0f}")
            print(f"  Tumor: final change = {final_change:.1f}% at day {df.iloc[-1]['time_d']:.0f}")
            if min_change < -30:
                print(f"  >> RECIST PR achieved (< -30%)")
            elif min_change < 0:
                print(f"  >> Stable disease (shrinkage but > -30%)")
            else:
                print(f"  >> Progressive disease")

            # Nadir → regrowth detection
            if t_nadir < df.iloc[-1]["time_d"] - 5 and final_change > min_change + 5:
                print(f"  >> Nadir at day {t_nadir:.0f}, then regrowth to {final_change:.1f}%")

        # Resistance state
        if "Z" in df.columns:
            z_final = df.iloc[-1]["Z"]
            z_max = df["Z"].max()
            print(f"  Resistance Z: final = {z_final:.3f}, max = {z_max:.3f}")


# ─────────────────────────────────────────
# CLI
# ─────────────────────────────────────────

def parse_args():
    ap = argparse.ArgumentParser(description="KRAS G12C QSP Model — integrated simulation")
    ap.add_argument("--scenario", type=str, default="all",
                    choices=["vehicle", "nsclc_mono", "crc_mono", "crc_combo", "all"],
                    help="Simulation scenario (default: all)")
    ap.add_argument("--t_days", type=int, default=180,
                    help="Simulation duration in days (default: 180)")
    ap.add_argument("--outdir", type=str, default="figures",
                    help="Output directory for figures")
    ap.add_argument("--csv_dir", type=str, default="outputs",
                    help="Output directory for simulation CSVs")
    ap.add_argument("--save_csv", action="store_true",
                    help="Also save simulation results as CSV")
    return ap.parse_args()


def main():
    args = parse_args()
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    # Determine which scenarios to run
    if args.scenario == "all":
        scenario_names = ["vehicle", "nsclc_mono", "crc_mono", "crc_combo"]
    else:
        scenario_names = [args.scenario]

    # Run simulations
    results: Dict[str, pd.DataFrame] = {}
    for name in scenario_names:
        p = SCENARIOS[name]()
        print(f"Simulating: {LABELS[name]} (t={args.t_days} days)...")
        df = simulate(p, t_days=args.t_days)
        results[name] = df

        if args.save_csv:
            os.makedirs(args.csv_dir, exist_ok=True)
            csv_path = os.path.join(args.csv_dir, f"sim_{name}.csv")
            df.to_csv(csv_path, index=False)
            print(f"  Saved {csv_path}")

    # Validation
    print_validation(results)

    # Generate figures
    print(f"\nGenerating figures in: {outdir}")
    fig1_pk_profile(results, outdir)
    fig2_target_engagement(results, outdir)

    if "nsclc_mono" in results and "crc_mono" in results:
        fig3_adaptive_rebound(results, outdir)
        fig5_crc_vs_nsclc(results, outdir)

    if "vehicle" in results:
        fig4_tumor_growth(results, outdir)

    if "crc_mono" in results and "crc_combo" in results:
        fig6_combination_benefit(results, outdir)

    if len(results) >= 3:
        fig7_model_summary(results, outdir)

    print(f"\nDone. All outputs in: {outdir}")


if __name__ == "__main__":
    main()
