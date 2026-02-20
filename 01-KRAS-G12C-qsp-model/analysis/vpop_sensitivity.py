#!/usr/bin/env python3
"""
M7: Virtual Population (Vpop) Sampling + Waterfall Plot
M8: Morris Global Sensitivity Analysis + Tornado Plot

Generates:
  - fig8_waterfall_vpop.png   (M7)
  - fig9_tornado_sensitivity.png (M8)
  - outputs/vpop_results.csv
  - outputs/sensitivity_morris.csv

Usage:
  python analysis/vpop_sensitivity.py                     # both M7 + M8
  python analysis/vpop_sensitivity.py --vpop-only         # M7 only
  python analysis/vpop_sensitivity.py --sensitivity-only  # M8 only
  python analysis/vpop_sensitivity.py --n_vpop 1000       # larger Vpop

Author: Hajar Besbassi
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.stats import norm as sp_norm, beta as sp_beta

# Import model components from main script
sys.path.insert(0, os.path.dirname(__file__))
from kras_qsp import (
    Params, rhs,
    params_nsclc_mono, params_crc_mono, params_crc_combo,
    COLORS, LABELS,
)


# ─────────────────────────────────────────
# Fast simulation (relaxed tolerances for screening)
# ─────────────────────────────────────────

def simulate_fast(p: Params, t_days: int = 90) -> dict:
    """
    Lightweight simulation returning only summary metrics.
    Uses LSODA, relaxed tolerances, coarse output for speed.
    """
    t_end_h = t_days * 24.0
    dose_times_h = np.arange(0.0, t_end_h, p.tau_h)

    GDP_ss, GTP_ss, _ = p.kras_ss()
    y0 = np.array([
        0.0, 0.0, GDP_ss, GTP_ss, 0.0,
        p.R_basal * (1.0 - p.f_block), p.E_ss, p.T_0, p.Z_0,
    ])

    y_current = y0.copy()
    T_0 = p.T_0

    # Track metrics
    best_change = 0.0
    t_best = 0.0
    T_nadir = T_0
    erk_72h = np.nan
    final_change = 0.0
    change_wk9 = np.nan   # tumor change at week 9 (first RECIST assessment)
    ttp_d = np.nan         # time to progression (T > T_nadir * 1.2 after nadir)

    event_times = list(dose_times_h)
    if event_times[-1] < t_end_h:
        event_times.append(t_end_h)

    for i in range(len(event_times)):
        t_dose = event_times[i]

        if i < len(dose_times_h) and abs(t_dose - dose_times_h[min(i, len(dose_times_h)-1)]) < 0.01:
            if i < len(dose_times_h):
                y_current[0] += p.F * p.dose_mg

        if i + 1 < len(event_times):
            t_next = event_times[i + 1]
        else:
            t_next = t_end_h

        if t_next <= t_dose:
            continue

        t_span = (t_dose, t_next)
        # Only 3 eval points per dose interval (coarse)
        t_eval = np.linspace(t_dose, t_next, 4, endpoint=True)

        sol = solve_ivp(
            lambda t, y: rhs(t, y, p, dose_times_h, GTP_ss),
            t_span, y_current, t_eval=t_eval,
            method="LSODA", rtol=1e-4, atol=1e-6, max_step=4.0,
        )

        if not sol.success:
            break

        # Extract metrics from this segment
        for k in range(sol.t.shape[0]):
            t_h = sol.t[k]
            T_val = sol.y[7, k]
            E_val = sol.y[6, k]
            change = (T_val / T_0 - 1.0) * 100.0

            if T_val < T_nadir:
                T_nadir = T_val
            if change < best_change:
                best_change = change
                t_best = t_h / 24.0

            # TTP: RECIST-like progression (20% increase from nadir + absolute increase)
            # Buffer 21 days past nadir to avoid catching PK oscillations
            if (np.isnan(ttp_d) and t_h / 24.0 > t_best + 21
                    and T_val > T_nadir * 1.2
                    and (T_val - T_nadir) > 0.05 * T_0):
                ttp_d = t_h / 24.0

            # ERK at 72h
            if 70 <= t_h <= 74:
                erk_72h = (E_val / p.E_ss) * 100.0

            # Tumor change at week 9 (~1512h, first RECIST assessment)
            if 1504 <= t_h <= 1520:
                change_wk9 = change

        y_current = sol.y[:, -1].copy()

    final_T = y_current[7]
    final_change = (final_T / T_0 - 1.0) * 100.0

    return {
        "best_change_pct": best_change,
        "final_change_pct": final_change,
        "change_wk9_pct": change_wk9 if not np.isnan(change_wk9) else final_change,
        "erk_rebound_72h": erk_72h,
        "time_nadir_d": t_best,
        "ttp_d": ttp_d,
    }


# ─────────────────────────────────────────
# Parameter sampling
# ─────────────────────────────────────────

VPOP_PARAMS = {
    "G_fb":        {"cv": 0.60},
    "k_GEF_basal": {"cv": 0.30},
    "k_bind":      {"cv": 0.35},
    "rho_grow":    {"cv": 0.40},
    "rho_kill":    {"cv": 0.50},
    "tau_R":       {"cv": 0.25},
    "k_Z_up":      {"cv": 0.70},
    "alpha_Z":     {"cv": 0.70},
    "beta_Z":      {"cv": 0.65},
}

# Z_0 (pre-existing resistance) sampled from Beta(a, b) * scale
# Encodes co-mutations (STK11, KEAP1), prior therapy, intrinsic heterogeneity
# ~40% of KRAS-mutant tumors have co-occurring resistance-associated alterations
Z0_BETA_A = 0.5    # heavy mass near 0 (drug-naive patients)
Z0_BETA_B = 1.5    # broader tail (prior-therapy/co-mutation patients)
Z0_SCALE = 0.6     # max pre-existing resistance = 0.6


def lhs_sample(n: int, n_params: int, rng: np.random.Generator) -> np.ndarray:
    """Latin Hypercube Sampling."""
    result = np.zeros((n, n_params))
    for j in range(n_params):
        perm = rng.permutation(n)
        for i in range(n):
            result[perm[i], j] = (i + rng.uniform()) / n
    return result


def sample_vpop(base_params: Params, n: int, seed: int = 42) -> List[Params]:
    """Generate n virtual patients by LHS around base parameters.

    Samples log-normal distributions for pharmacological parameters,
    Beta distribution for pre-existing resistance (Z_0), and
    log-normal for f_block variability (combo only).
    """
    rng = np.random.default_rng(seed)
    param_names = list(VPOP_PARAMS.keys())

    # LHS for main parameters + 1 extra dim for Z_0 + 1 for f_block
    n_extra = 2  # Z_0 and f_block
    U = lhs_sample(n, len(param_names) + n_extra, rng)

    patients = []
    for i in range(n):
        p = Params(
            G_fb=base_params.G_fb,
            f_block=base_params.f_block,
            rho_grow=base_params.rho_grow,
        )
        for j, name in enumerate(param_names):
            nominal = getattr(base_params, name)
            cv = VPOP_PARAMS[name]["cv"]
            sigma = np.sqrt(np.log(1 + cv**2))
            mu = np.log(nominal) - 0.5 * sigma**2
            val = np.exp(mu + sigma * sp_norm.ppf(U[i, j]))
            setattr(p, name, float(val))

        # Z_0: pre-existing resistance from Beta distribution
        u_z0 = U[i, len(param_names)]
        z0_val = sp_beta.ppf(u_z0, Z0_BETA_A, Z0_BETA_B) * Z0_SCALE
        p.Z_0 = float(np.clip(z0_val, 0.0, 0.7))

        # f_block: vary panitumumab exposure for combo scenarios
        if base_params.f_block > 0:
            u_fb = U[i, len(param_names) + 1]
            cv_fblock = 0.25
            sigma_fb = np.sqrt(np.log(1 + cv_fblock**2))
            mu_fb = np.log(base_params.f_block) - 0.5 * sigma_fb**2
            p.f_block = float(np.clip(
                np.exp(mu_fb + sigma_fb * sp_norm.ppf(u_fb)), 0.0, 0.85
            ))

        patients.append(p)
    return patients


# ─────────────────────────────────────────
# M7: Virtual Population
# ─────────────────────────────────────────

def run_m7(n_vpop: int, t_days: int, outdir: str, csv_dir: str) -> pd.DataFrame:
    """Run M7: Virtual population for all treatment scenarios."""
    print(f"\n{'='*60}")
    print(f"M7: VIRTUAL POPULATION (N={n_vpop} per scenario)")
    print(f"{'='*60}")

    scenarios = {
        "nsclc_mono": params_nsclc_mono,
        "crc_mono": params_crc_mono,
        "crc_combo": params_crc_combo,
    }

    all_rows = []
    for scn_name, scn_fn in scenarios.items():
        t0 = time.time()
        print(f"\n  Running {scn_name}...")
        base = scn_fn()
        patients = sample_vpop(base, n_vpop)

        for i, p in enumerate(patients):
            metrics = simulate_fast(p, t_days=t_days)
            row = {
                "patient": i + 1,
                "scenario": scn_name,
                **metrics,
                "G_fb": p.G_fb,
                "k_GEF_basal": p.k_GEF_basal,
                "k_bind": p.k_bind,
                "rho_grow": p.rho_grow,
                "rho_kill": p.rho_kill,
                "tau_R": p.tau_R,
                "k_Z_up": p.k_Z_up,
                "alpha_Z": p.alpha_Z,
                "beta_Z": p.beta_Z,
                "Z_0": p.Z_0,
                "f_block": p.f_block,
            }
            all_rows.append(row)

            if (i + 1) % 100 == 0:
                print(f"    {i+1}/{n_vpop} patients done")

        elapsed = time.time() - t0
        sub = [r for r in all_rows if r["scenario"] == scn_name]
        orr_nadir = sum(1 for r in sub if r["best_change_pct"] < -30) / len(sub) * 100
        orr_wk9 = sum(1 for r in sub if r["change_wk9_pct"] < -30) / len(sub) * 100
        med = np.median([r["best_change_pct"] for r in sub])
        # PFS: use TTP if progressed, else censored at t_days
        pfs_vals = [r["ttp_d"] if not np.isnan(r["ttp_d"]) else t_days for r in sub]
        mpfs = np.median(pfs_vals) / 30.0  # convert to months
        n_progressed = sum(1 for r in sub if not np.isnan(r["ttp_d"]))
        print(f"    Done in {elapsed:.1f}s — ORR(wk9): {orr_wk9:.0f}%, median best: {med:.1f}%, mPFS: {mpfs:.1f} mo ({n_progressed}/{len(sub)} progressed)")

    vpop_df = pd.DataFrame(all_rows)

    os.makedirs(csv_dir, exist_ok=True)
    vpop_df.to_csv(os.path.join(csv_dir, "vpop_results.csv"), index=False)
    print(f"\n  Saved {csv_dir}/vpop_results.csv")

    fig8_waterfall(vpop_df, outdir)
    return vpop_df


def fig8_waterfall(vpop_df: pd.DataFrame, outdir: str) -> None:
    """Fig 8: Waterfall plot of best tumor change per virtual patient."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=True)

    scenario_order = ["nsclc_mono", "crc_mono", "crc_combo"]
    titles = {
        "nsclc_mono": "NSCLC — Sotorasib Mono",
        "crc_mono": "CRC — Sotorasib Mono",
        "crc_combo": "CRC — Sotorasib + Panitumumab",
    }

    for idx, scn in enumerate(scenario_order):
        ax = axes[idx]
        sub = vpop_df[vpop_df["scenario"] == scn].copy()
        sub = sub.sort_values("best_change_pct", ascending=True).reset_index(drop=True)
        n = len(sub)

        colors = []
        for val in sub["best_change_pct"]:
            if val < -30:
                colors.append("#2166ac")    # PR
            elif val < 0:
                colors.append("#92c5de")    # SD shrinkage
            elif val < 20:
                colors.append("#fddbc7")    # SD growth
            else:
                colors.append("#b2182b")    # PD

        ax.bar(range(n), sub["best_change_pct"], color=colors, width=1.0, edgecolor="none")
        ax.axhline(-30, color="black", linewidth=0.8, linestyle="--", alpha=0.7)
        ax.axhline(20, color="black", linewidth=0.8, linestyle=":", alpha=0.5)
        ax.axhline(0, color="black", linewidth=0.3)

        orr_wk9 = (sub["change_wk9_pct"] < -30).mean() * 100 if "change_wk9_pct" in sub.columns else 0
        # PFS: use TTP if progressed, else censored at sim end
        pfs_vals = sub["ttp_d"].fillna(180.0)  # censor at 180 days
        mpfs_mo = np.median(pfs_vals) / 30.0
        ax.text(0.95, 0.95, f"ORR(wk9): {orr_wk9:.0f}%\nmPFS: {mpfs_mo:.1f} mo",
                transform=ax.transAxes, fontsize=10, fontweight="bold",
                ha="right", va="top",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                         edgecolor="gray", alpha=0.9))

        ax.set_xlabel("Virtual patients (ranked)", fontsize=10)
        if idx == 0:
            ax.set_ylabel("Best tumor change from baseline (%)", fontsize=10)
        ax.set_title(titles[scn], fontsize=11, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xlim(-0.5, n - 0.5)

    ymin = vpop_df["best_change_pct"].min()
    ymax = vpop_df["best_change_pct"].max()
    for ax in axes:
        ax.set_ylim(max(ymin - 5, -100), min(ymax + 10, 150))

    fig.suptitle("Virtual Population: Best Tumor Change (Waterfall Plot)",
                 fontsize=14, fontweight="bold", y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig8_waterfall_vpop.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved fig8_waterfall_vpop.png")


# ─────────────────────────────────────────
# M8: Morris Sensitivity Analysis
# ─────────────────────────────────────────

def run_m8(outdir: str, csv_dir: str, t_days: int = 90) -> pd.DataFrame:
    """Run M8: Morris global sensitivity analysis."""
    print(f"\n{'='*60}")
    print(f"M8: MORRIS GLOBAL SENSITIVITY ANALYSIS")
    print(f"{'='*60}")

    from SALib.sample import morris as morris_sample
    from SALib.analyze import morris as morris_analyze

    param_names = list(VPOP_PARAMS.keys())
    scenarios = {
        "nsclc_mono": params_nsclc_mono,
        "crc_mono": params_crc_mono,
        "crc_combo": params_crc_combo,
    }

    all_sa_rows = []

    for scn_name, scn_fn in scenarios.items():
        print(f"\n  Scenario: {scn_name}")
        base = scn_fn()

        bounds = []
        for pname in param_names:
            nominal = getattr(base, pname)
            bounds.append([nominal * 0.5, nominal * 2.0])

        problem = {
            "num_vars": len(param_names),
            "names": param_names,
            "bounds": bounds,
        }

        n_traj = 20
        X = morris_sample.sample(problem, N=n_traj, seed=123)
        n_runs = X.shape[0]
        print(f"    Morris samples: {n_runs} runs ({n_traj} trajectories)")

        Y_tumor = np.zeros(n_runs)
        Y_erk = np.zeros(n_runs)

        t0 = time.time()
        for i in range(n_runs):
            p = Params(G_fb=base.G_fb, f_block=base.f_block, rho_grow=base.rho_grow)
            for j, pname in enumerate(param_names):
                setattr(p, pname, float(X[i, j]))

            metrics = simulate_fast(p, t_days=t_days)
            Y_tumor[i] = metrics["best_change_pct"]
            Y_erk[i] = metrics["erk_rebound_72h"] if not np.isnan(metrics["erk_rebound_72h"]) else 0.0

            if (i + 1) % 50 == 0:
                print(f"      {i+1}/{n_runs} runs done")

        elapsed = time.time() - t0
        print(f"    Done in {elapsed:.1f}s")

        for outcome_name, Y in [("best_tumor_change", Y_tumor), ("erk_rebound_72h", Y_erk)]:
            Y_clean = np.nan_to_num(Y, nan=0.0)
            Si = morris_analyze.analyze(problem, X, Y_clean, seed=123)

            for j, pname in enumerate(param_names):
                all_sa_rows.append({
                    "scenario": scn_name,
                    "outcome": outcome_name,
                    "parameter": pname,
                    "mu_star": float(Si["mu_star"][j]),
                    "mu": float(Si["mu"][j]),
                    "sigma": float(Si["sigma"][j]),
                })

    sa_df = pd.DataFrame(all_sa_rows)

    os.makedirs(csv_dir, exist_ok=True)
    sa_df.to_csv(os.path.join(csv_dir, "sensitivity_morris.csv"), index=False)
    print(f"\n  Saved {csv_dir}/sensitivity_morris.csv")

    fig9_tornado(sa_df, outdir)
    return sa_df


PARAM_DISPLAY = {
    "G_fb": "Feedback gain\n(G_fb)",
    "k_GEF_basal": "GEF exchange\n(k_GEF)",
    "k_bind": "Drug binding\n(k_bind)",
    "rho_grow": "Tumor growth\n(rho_grow)",
    "rho_kill": "Drug kill rate\n(rho_kill)",
    "tau_R": "Receptor tau\n(tau_R)",
    "k_Z_up": "Resistance rate\n(k_Z_up)",
    "alpha_Z": "Resistance strength\n(alpha_Z)",
    "beta_Z": "Bypass signaling\n(beta_Z)",
}


def fig9_tornado(sa_df: pd.DataFrame, outdir: str) -> None:
    """Fig 9: Tornado plot of Morris mu* for tumor change across scenarios."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=True)

    scenario_order = ["nsclc_mono", "crc_mono", "crc_combo"]
    titles = {
        "nsclc_mono": "NSCLC Mono",
        "crc_mono": "CRC Mono",
        "crc_combo": "CRC Combo",
    }
    bar_colors = {
        "nsclc_mono": "#377eb8",
        "crc_mono": "#e41a1c",
        "crc_combo": "#4daf4a",
    }

    for idx, scn in enumerate(scenario_order):
        ax = axes[idx]
        sub = sa_df[(sa_df["scenario"] == scn) & (sa_df["outcome"] == "best_tumor_change")]
        sub = sub.sort_values("mu_star", ascending=True)

        y_pos = range(len(sub))
        labels = [PARAM_DISPLAY.get(p, p) for p in sub["parameter"]]

        ax.barh(y_pos, sub["mu_star"], xerr=sub["sigma"],
                color=bar_colors[scn], alpha=0.85, edgecolor="white",
                capsize=3, height=0.6)
        ax.set_yticks(list(y_pos))
        ax.set_yticklabels(labels, fontsize=9)
        ax.set_xlabel("Morris mu* (mean |elementary effect|)", fontsize=9)
        ax.set_title(titles[scn], fontsize=11, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.suptitle("Global Sensitivity Analysis: Best Tumor Change (Morris Screening)",
                 fontsize=14, fontweight="bold", y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig9_tornado_sensitivity.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved fig9_tornado_sensitivity.png")


# ─────────────────────────────────────────
# CLI
# ─────────────────────────────────────────

def parse_args():
    ap = argparse.ArgumentParser(
        description="M7: Virtual Population + M8: Morris Sensitivity Analysis"
    )
    ap.add_argument("--n_vpop", type=int, default=200,
                    help="Number of virtual patients per scenario (default: 200)")
    ap.add_argument("--t_days", type=int, default=90,
                    help="Simulation duration in days (default: 90)")
    ap.add_argument("--outdir", type=str, default="figures",
                    help="Output directory for figures")
    ap.add_argument("--csv_dir", type=str, default="outputs",
                    help="Output directory for CSVs")
    ap.add_argument("--vpop-only", action="store_true",
                    help="Run only M7 (Vpop)")
    ap.add_argument("--sensitivity-only", action="store_true",
                    help="Run only M8 (Morris)")
    ap.add_argument("--seed", type=int, default=42,
                    help="Random seed (default: 42)")
    return ap.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    run_vpop = not args.sensitivity_only
    run_sa = not args.vpop_only

    t_total = time.time()

    if run_vpop:
        run_m7(args.n_vpop, args.t_days, args.outdir, args.csv_dir)

    if run_sa:
        run_m8(args.outdir, args.csv_dir, args.t_days)

    elapsed = time.time() - t_total
    print(f"\n{'='*60}")
    print(f"All done in {elapsed:.0f}s")
    print(f"Figures: {args.outdir}/fig8_waterfall_vpop.png, fig9_tornado_sensitivity.png")
    print(f"Data: {args.csv_dir}/vpop_results.csv, sensitivity_morris.csv")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
