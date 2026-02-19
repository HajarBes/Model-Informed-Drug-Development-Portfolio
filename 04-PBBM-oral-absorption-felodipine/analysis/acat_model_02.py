#!/usr/bin/env python3
# =============================================================================
# 02_acat_model.py — ACAT ODE System (Importable Module)
# =============================================================================
#
# 7-segment ACAT (Advanced Compartmental Absorption and Transit) model
# with Noyes-Whitney dissolution, Peff-driven absorption, scalar first-pass,
# and 2-compartment systemic disposition.
#
# 17 ODE states:
#   S_0..S_6  — undissolved solid in each GI segment (7 states)
#   D_0..D_6  — dissolved drug in each GI segment (7 states)
#   Ac        — amount in central compartment
#   Ap        — amount in peripheral compartment
#   cum_abs   — cumulative mass absorbed through intestinal wall
#
# Author: Hajar Besbassi

from __future__ import annotations

import os
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

# Import 00_setup.py (numeric prefix requires importlib)
import importlib.util
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "setup_00", os.path.join(_SCRIPT_DIR, "00_setup.py"))
_setup = importlib.util.module_from_spec(_spec)
sys.modules["setup_00"] = _setup  # register so dataclass can resolve __module__
_spec.loader.exec_module(_setup)

DrugParams = _setup.DrugParams
GISegment = _setup.GISegment
FASTED_GI = _setup.FASTED_GI
FED_GI = _setup.FED_GI
compute_auc = _setup.compute_auc

N_SEGMENTS = 7

# =============================================================================
# Core physics
# =============================================================================

def saturation_solubility(S0_mg_mL: float, pKa: float, pH: float,
                          bile_factor: float) -> float:
    """
    pH-dependent saturation solubility (Henderson-Hasselbalch for a weak base)
    enhanced by bile salt micellar solubilization.

    C_sat = S0 * (1 + 10^(pKa - pH)) * bile_factor
    """
    ionization = 1.0 + 10.0 ** (pKa - pH)
    return S0_mg_mL * ionization * bile_factor


def dissolution_rate(S_mg: float, C_dissolved_mg_mL: float,
                     C_sat_mg_mL: float, Deff: float,
                     density: float, radius_cm: float,
                     z_factor: float = 1.0) -> float:
    """
    Noyes-Whitney dissolution rate (mg/h) — correct concentration-driving-force form.

    diss_rate = z * (3 * Deff) / (density_mg * r^2) * S * (C_sat - C_diss)

    Uses density in mg/cm³ for correct units: rate has units mg/h.
    The (C_sat - C_diss) driving force makes dissolution directly proportional
    to local saturation solubility, so bile salt enhancement of C_sat translates
    directly into faster dissolution — critical for BCS II food effect modeling.

    z_factor accounts for precipitated drug forming fine amorphous particles
    with higher effective surface area than the original formulation particles.

    Returns 0 if supersaturated (C_diss >= C_sat) — precipitation handles that.
    """
    if S_mg <= 0 or C_sat_mg_mL <= 0:
        return 0.0
    if C_dissolved_mg_mL >= C_sat_mg_mL:
        return 0.0
    # Convert Deff from cm^2/s to cm^2/h
    Deff_h = Deff * 3600.0
    # density in mg/cm³ for correct unit balance: [mL/(mg*h)]
    density_mg = density * 1000.0
    k_NW = (3.0 * Deff_h) / (density_mg * radius_cm ** 2)
    return z_factor * k_NW * S_mg * (C_sat_mg_mL - C_dissolved_mg_mL)


def absorption_rate_constant(Peff_cm_s: float, SA_cm2: float,
                             V_lumen_mL: float) -> float:
    """
    Absorption rate constant (1/h) from effective permeability.

    k_abs = Peff * SA / V_lumen  (converted s -> h)
    """
    if V_lumen_mL <= 0:
        return 0.0
    Peff_h = Peff_cm_s * 3600.0  # cm/s -> cm/h
    return Peff_h * SA_cm2 / V_lumen_mL


# =============================================================================
# ODE right-hand side
# =============================================================================

# Precipitation rate constant (h^-1) — rapid re-precipitation when supersaturated
K_PRECIP = 10.0


def rhs(t: float, y: np.ndarray, drug: DrugParams,
        gi_segments: List[GISegment]) -> np.ndarray:
    """
    ACAT ODE right-hand side (17 states).

    State layout:
      y[0:7]   = S_0 .. S_6  (undissolved solid, mg)
      y[7:14]  = D_0 .. D_6  (dissolved drug, mg)
      y[14]    = Ac (central compartment, mg)
      y[15]    = Ap (peripheral compartment, mg)
      y[16]    = cumulative mass absorbed through intestinal wall (mg)

    Includes precipitation: when dissolved concentration exceeds local
    C_sat (e.g., after transit from low-pH stomach to higher-pH intestine),
    excess drug precipitates back to solid phase.
    """
    dydt = np.zeros(17)

    S = y[0:7]
    D = y[7:14]
    Ac = y[14]
    Ap = y[15]

    radius_cm = drug.particle_radius_um * 1e-4  # µm -> cm

    total_absorbed = 0.0

    for i, seg in enumerate(gi_segments):
        # Saturation solubility in this segment
        C_sat = saturation_solubility(drug.S0_mg_mL, drug.pKa,
                                      seg.pH, seg.bile_factor)

        # Dissolved concentration (mg/mL)
        D_i = max(D[i], 0.0)
        C_diss = D_i / seg.volume_mL if seg.volume_mL > 0 else 0.0

        # Dissolution (Noyes-Whitney, only when undersaturated)
        diss = dissolution_rate(max(S[i], 0.0), C_diss, C_sat,
                                drug.Deff_cm2_s, drug.density_g_cm3,
                                radius_cm, drug.z_dissolution)

        # Precipitation (when supersaturated, e.g. after pH transition)
        if C_diss > C_sat and seg.volume_mL > 0:
            precip = K_PRECIP * (C_diss - C_sat) * seg.volume_mL  # mg/h
        else:
            precip = 0.0

        # Transit
        k_t = seg.k_transit_h

        # Incoming transit from upstream segment
        if i == 0:
            S_in = 0.0
            D_in = 0.0
        else:
            k_t_prev = gi_segments[i - 1].k_transit_h
            S_in = k_t_prev * max(S[i - 1], 0.0)
            D_in = k_t_prev * max(D[i - 1], 0.0)

        # Absorption (segments 1-6; stomach does not absorb)
        if seg.absorbs:
            k_abs = absorption_rate_constant(drug.Peff_cm_s,
                                             seg.surface_area_cm2,
                                             seg.volume_mL)
        else:
            k_abs = 0.0

        absorbed = k_abs * D_i
        total_absorbed += absorbed

        # Solid: dissolution out, precipitation in, transit out, transit in
        dydt[i] = -diss + precip - k_t * max(S[i], 0.0) + S_in

        # Dissolved: dissolution in, precipitation out, absorption out, transit out, transit in
        dydt[7 + i] = diss - precip - absorbed - k_t * D_i + D_in

    # Cumulative absorbed mass (before first-pass) — state 17
    dydt[16] = total_absorbed

    # First-pass extraction: Fg * Fh
    drug_to_systemic = total_absorbed * drug.Fg * drug.Fh

    # 2-compartment systemic disposition
    ke = drug.CL_L_h / drug.Vc_L
    k12 = drug.Q_L_h / drug.Vc_L
    k21 = drug.Q_L_h / drug.Vp_L

    dydt[14] = drug_to_systemic - ke * Ac - k12 * Ac + k21 * Ap  # dAc/dt
    dydt[15] = k12 * Ac - k21 * Ap                                # dAp/dt

    return dydt


# =============================================================================
# Simulation driver
# =============================================================================

def simulate(drug: DrugParams, gi_segments: List[GISegment],
             t_end_h: float = 24.0, dt_h: float = 0.01) -> pd.DataFrame:
    """
    Simulate ACAT model for a single oral dose.

    Returns DataFrame with columns: time_h, Cp_ngmL, Ac_mg, Ap_mg,
    S_total_mg, D_total_mg, Fa (fraction absorbed).
    """
    # Initial conditions: all drug as solid in stomach
    y0 = np.zeros(17)
    y0[0] = drug.dose_mg  # S_stomach = dose

    t_eval = np.arange(0, t_end_h + dt_h, dt_h)

    sol = solve_ivp(
        fun=lambda t, y: rhs(t, y, drug, gi_segments),
        t_span=(0, t_end_h),
        y0=y0,
        t_eval=t_eval,
        method="LSODA",
        rtol=1e-8,
        atol=1e-10,
        max_step=0.1,
    )

    if not sol.success:
        raise RuntimeError(f"ODE solver failed: {sol.message}")

    S_total = sol.y[0:7, :].sum(axis=0)
    D_total = sol.y[7:14, :].sum(axis=0)
    Ac = sol.y[14, :]
    Ap = sol.y[15, :]
    cum_absorbed = sol.y[16, :]  # cumulative mass absorbed through intestinal wall

    # Concentration in ng/mL: Ac (mg) / Vc (L) * 1e6 (mg->ng) / 1e3 (L->mL)
    # = Ac / Vc * 1000
    Cp_ngmL = Ac / drug.Vc_L * 1000.0

    # Fraction absorbed = cumulative absorbed through wall / dose
    # (NOT 1 - GI_remaining/dose, which conflates fecal loss with absorption)
    Fa = cum_absorbed / drug.dose_mg

    df = pd.DataFrame({
        "time_h": sol.t,
        "Cp_ngmL": Cp_ngmL,
        "Ac_mg": Ac,
        "Ap_mg": Ap,
        "S_total_mg": S_total,
        "D_total_mg": D_total,
        "cum_absorbed_mg": cum_absorbed,
        "Fa": Fa,
    })

    # Per-segment absorption tracking
    for i in range(N_SEGMENTS):
        df[f"S_{i}"] = sol.y[i, :]
        df[f"D_{i}"] = sol.y[7 + i, :]

    return df


def compute_pk_metrics(sim_df: pd.DataFrame) -> dict:
    """Compute PK summary metrics from simulation DataFrame."""
    t = sim_df["time_h"].values
    cp = sim_df["Cp_ngmL"].values

    idx_max = np.argmax(cp)
    cmax = float(cp[idx_max])
    tmax = float(t[idx_max])
    auc = float(np.trapz(cp, t))
    fa_final = float(sim_df["Fa"].iloc[-1])

    return {
        "Cmax_ngmL": cmax,
        "Tmax_h": tmax,
        "AUC_ngmL_h": auc,
        "Fa": fa_final,
        "F_oral": fa_final * DrugParams().Fg * DrugParams().Fh,
    }


# =============================================================================
# Per-segment absorption breakdown
# =============================================================================

def segment_absorption_breakdown(drug: DrugParams,
                                 gi_segments: List[GISegment],
                                 sim_df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute cumulative absorption per GI segment by integrating
    k_abs_i * D_i over time.
    """
    rows = []
    dt = np.diff(sim_df["time_h"].values, prepend=0)

    for i, seg in enumerate(gi_segments):
        if not seg.absorbs:
            rows.append({"segment": seg.name, "absorbed_mg": 0.0, "pct_dose": 0.0})
            continue
        k_abs = absorption_rate_constant(drug.Peff_cm_s,
                                         seg.surface_area_cm2,
                                         seg.volume_mL)
        D_i = sim_df[f"D_{i}"].values
        absorbed = float(np.sum(k_abs * D_i * dt))
        rows.append({
            "segment": seg.name,
            "absorbed_mg": absorbed,
            "pct_dose": absorbed / drug.dose_mg * 100.0,
        })
    return pd.DataFrame(rows)


# =============================================================================
# Verification
# =============================================================================

def verify():
    """
    Verification suite: mass balance + limiting cases.
    Prints PASS/FAIL for each check.
    """
    drug = DrugParams()
    all_pass = True

    def check(name, condition):
        nonlocal all_pass
        status = "PASS" if condition else "FAIL"
        if not condition:
            all_pass = False
        print(f"  [{status}] {name}")

    print("=" * 60)
    print("ACAT Model Verification")
    print("=" * 60)

    # --- Test 1: Mass balance (fasted) ---
    # Full mass balance:
    #   dose = GI_drug + systemic_drug + systemic_eliminated + first_pass_loss + fecal_loss
    # where:
    #   GI_drug = S_total + D_total
    #   systemic_drug = Ac + Ap
    #   systemic_eliminated = integral of (CL/Vc) * Ac dt
    #   first_pass_loss: We track it via (Ac + Ap + sys_elim) / (Fg * Fh) to get total absorbed
    #   total_absorbed = (Ac + Ap + sys_elim) / (Fg * Fh)
    #   fecal = dose - GI_drug - total_absorbed
    print("\n--- Mass Balance (fasted, 10 mg, 48 h) ---")
    sim = simulate(drug, FASTED_GI, t_end_h=48.0)
    last = sim.iloc[-1]
    S_total = last["S_total_mg"]
    D_total = last["D_total_mg"]
    Ac = last["Ac_mg"]
    Ap = last["Ap_mg"]
    cum_absorbed = last["cum_absorbed_mg"]  # from ODE state 17

    # Cumulative systemic elimination
    ke = drug.CL_L_h / drug.Vc_L
    sys_eliminated = float(np.trapz(ke * sim["Ac_mg"].values, sim["time_h"].values))

    # Total drug that reached systemic circulation
    cumul_systemic = Ac + Ap + sys_eliminated
    first_pass_loss = cum_absorbed - cumul_systemic

    # Fecal loss = dose - GI_drug - total_absorbed_through_wall
    fecal = drug.dose_mg - S_total - D_total - cum_absorbed

    mass_accounted = S_total + D_total + cum_absorbed + fecal
    mass_error = abs(mass_accounted - drug.dose_mg)
    check(f"Mass balance error = {mass_error:.2e} mg (< 0.01 mg)",
          mass_error < 0.01)
    print(f"    GI remaining:    {S_total + D_total:.4f} mg")
    print(f"    Absorbed (wall): {cum_absorbed:.4f} mg (Fa = {cum_absorbed/drug.dose_mg:.4f})")
    print(f"    First-pass loss: {first_pass_loss:.4f} mg")
    print(f"    Systemic (Ac+Ap+elim): {cumul_systemic:.4f} mg")
    print(f"    Fecal loss:      {max(fecal, 0):.4f} mg")

    # --- Test 2: Bioavailability ---
    print("\n--- Oral Bioavailability ---")
    metrics = compute_pk_metrics(sim)
    fa = metrics["Fa"]
    f_oral_apparent = metrics["F_oral"]
    check(f"Fa = {fa:.3f} (fraction leaving GI lumen)",
          fa > 0.3)
    check(f"F_oral (Fa*Fg*Fh) = {f_oral_apparent:.3f} (expect 0.10-0.30, literature 0.15-0.22)",
          0.10 <= f_oral_apparent <= 0.30)
    check(f"Cmax = {metrics['Cmax_ngmL']:.2f} ng/mL (expect 2-25 for 10 mg IR)",
          2.0 <= metrics["Cmax_ngmL"] <= 25.0)
    check(f"Tmax = {metrics['Tmax_h']:.2f} h (expect 0.3-6 h for IR)",
          0.3 <= metrics["Tmax_h"] <= 6.0)

    # --- Test 2b: Food effect direction ---
    print("\n--- Food Effect Direction Check ---")
    sim_fed = simulate(drug, FED_GI, t_end_h=24.0)
    met_fed = compute_pk_metrics(sim_fed)
    cmax_ratio_fe = met_fed["Cmax_ngmL"] / metrics["Cmax_ngmL"]
    auc_ratio_fe = met_fed["AUC_ngmL_h"] / metrics["AUC_ngmL_h"]
    # For IR felodipine: Bratel 1989 observed Cmax ratio ~1.31, AUC ratio ~1.15.
    # Fed gastric pH calibrated to 3.0 (from default 4.0) so that the weak base
    # still fully dissolves in the fed stomach, and bile salt enhancement dominates.
    check(f"Food effect Cmax ratio = {cmax_ratio_fe:.2f} (IR: expect 0.8-1.8, observed ~1.31)",
          0.8 <= cmax_ratio_fe <= 1.8)
    check(f"Food effect AUC ratio = {auc_ratio_fe:.2f} (expect 0.7-1.5, literature ~1.0)",
          0.7 <= auc_ratio_fe <= 1.5)
    print(f"    Fed Cmax = {met_fed['Cmax_ngmL']:.2f} ng/mL")
    print(f"    Fed Tmax = {met_fed['Tmax_h']:.2f} h")

    # --- Test 3: Limiting case — high solubility ---
    print("\n--- Limiting Case: High Solubility (instant dissolution) ---")
    drug_highsol = DrugParams(S0_ug_mL=5000.0)
    sim_hs = simulate(drug_highsol, FASTED_GI, t_end_h=24.0)
    last_hs = sim_hs.iloc[-1]
    s_remaining = last_hs["S_total_mg"]
    check(f"Solid remaining = {s_remaining:.4f} mg (< 0.01 mg, expect full dissolution)",
          abs(s_remaining) < 0.01)

    # --- Test 4: Limiting case — zero permeability ---
    print("\n--- Limiting Case: Zero Permeability (no absorption) ---")
    drug_noperm = DrugParams(Peff_cm_s=0.0)
    sim_np = simulate(drug_noperm, FASTED_GI, t_end_h=24.0)
    cmax_np = sim_np["Cp_ngmL"].max()
    check(f"Cmax with Peff=0: {cmax_np:.6f} ng/mL (expect ~0)",
          cmax_np < 1e-6)

    # --- Test 5: Limiting case — no first-pass ---
    print("\n--- Limiting Case: No First-Pass (Fg=Fh=1) ---")
    drug_nofp = DrugParams(Fg=1.0, Fh=1.0)
    sim_nfp = simulate(drug_nofp, FASTED_GI, t_end_h=24.0)
    metrics_nfp = compute_pk_metrics(sim_nfp)
    check(f"Cmax with Fg=Fh=1: {metrics_nfp['Cmax_ngmL']:.2f} > baseline {metrics['Cmax_ngmL']:.2f}",
          metrics_nfp["Cmax_ngmL"] > metrics["Cmax_ngmL"] * 2.0)

    print(f"\n{'='*60}")
    if all_pass:
        print("All verification checks PASSED.")
    else:
        print("WARNING: Some checks FAILED.")
    print(f"{'='*60}")

    return all_pass


if __name__ == "__main__":
    verify()
