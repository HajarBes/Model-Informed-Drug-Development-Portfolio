#!/usr/bin/env python3
"""
Minimal oral 1-compartment PK model for sotorasib (portfolio v0.1).

Model (amounts):
  dA_gut/dt  = -ka * A_gut + dosing_events(t)
  dA_central/dt =  ka * A_gut - ke * A_central

Concentration:
  C(t) = A_central / V

We pick ke from half-life: ke = ln(2)/t_half.
We choose ka so Tmax ~ 1 h (roughly consistent with label), but ka is adjustable.
"""

import argparse
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


@dataclass
class PKParams:
    dose_mg: float = 960.0
    tau_h: float = 24.0
    n_doses: int = 10
    t_half_h: float = 5.0
    ka_per_h: float = 1.0  # adjustable
    V_L: float = 50.0      # placeholder scaling (portfolio); later calibrate if needed
    F: float = 1.0         # bioavailability placeholder


def simulate_pk(p: PKParams, dt: float = 0.01):
    ke = math.log(2) / p.t_half_h
    t_end = p.tau_h * p.n_doses
    t = np.arange(0.0, t_end + dt, dt)

    A_gut = np.zeros_like(t)
    A_cent = np.zeros_like(t)

    # dosing times: 0, tau, 2*tau, ...
    dose_times = set(int(round(x / dt)) for x in np.arange(0.0, p.tau_h * p.n_doses, p.tau_h))

    for i in range(1, len(t)):
        # dosing into gut as amount (mg)
        if i in dose_times:
            A_gut[i-1] += p.F * p.dose_mg

        dA_gut = -p.ka_per_h * A_gut[i-1]
        dA_cent = p.ka_per_h * A_gut[i-1] - ke * A_cent[i-1]

        A_gut[i] = A_gut[i-1] + dt * dA_gut
        A_cent[i] = A_cent[i-1] + dt * dA_cent

    C = A_cent / p.V_L  # mg/L (units depend on V)
    out = pd.DataFrame({
        "time_h": t,
        "A_gut_mg": A_gut,
        "A_central_mg": A_cent,
        "C_mg_per_L": C,
    })
    return out, ke


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dose_mg", type=float, default=960.0)
    ap.add_argument("--tau_h", type=float, default=24.0)
    ap.add_argument("--n_doses", type=int, default=10)
    ap.add_argument("--t_half_h", type=float, default=5.0)
    ap.add_argument("--ka", type=float, default=1.0)
    ap.add_argument("--V", type=float, default=50.0)
    ap.add_argument("--out_csv", type=str, default="data/processed/pk_sotorasib_onecomp.csv")
    ap.add_argument("--out_fig", type=str, default="figures/pk_sotorasib_onecomp.png")
    args = ap.parse_args()

    p = PKParams(
        dose_mg=args.dose_mg,
        tau_h=args.tau_h,
        n_doses=args.n_doses,
        t_half_h=args.t_half_h,
        ka_per_h=args.ka,
        V_L=args.V,
    )

    df, ke = simulate_pk(p)

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)

    out_fig = Path(args.out_fig)
    out_fig.parent.mkdir(parents=True, exist_ok=True)

    # plot last dosing interval concentration (steady-ish state)
    t_end = p.tau_h * p.n_doses
    t0 = t_end - p.tau_h
    sub = df[(df["time_h"] >= t0) & (df["time_h"] <= t_end)].copy()

    plt.figure()
    plt.plot(sub["time_h"] - t0, sub["C_mg_per_L"])
    plt.xlabel("Time since last dose (h)")
    plt.ylabel("Concentration (mg/L)")
    plt.title(f"Oral 1-comp PK (ke={ke:.3f}/h, t1/2={p.t_half_h} h, ka={p.ka_per_h}/h)")
    plt.tight_layout()
    plt.savefig(out_fig, dpi=200)
    plt.close()

    print(f"[OK] ke = {ke:.4f} 1/h (from t1/2={p.t_half_h} h)")
    print(f"[OK] wrote {out_csv}")
    print(f"[OK] wrote {out_fig}")


if __name__ == "__main__":
    main()
