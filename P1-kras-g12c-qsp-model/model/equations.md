# Model Equations — KRAS G12C QSP Model

## State Variables (9)

| # | Variable | Description | Units | IC |
|---|----------|-------------|-------|-----|
| 1 | A_gut | Drug in gut compartment | mg | 0 (dose added as bolus) |
| 2 | A_central | Drug in central compartment | mg | 0 |
| 3 | KRAS_GDP | Inactive (GDP-bound) KRAS | µM | ~0.90 (from SS) |
| 4 | KRAS_GTP | Active (GTP-bound) KRAS | µM | ~0.10 (from SS) |
| 5 | KRAS_drug | Covalently inhibited KRAS | µM | 0 |
| 6 | R | Receptor drive (RTK proxy) | dimensionless | 1.0 |
| 7 | E | Pathway output (ERK proxy) | dimensionless | 1.0 |
| 8 | T | Tumor volume | normalized | 1.0 |
| 9 | Z | Adaptive resistance program | dimensionless [0,1] | 0.0 (Z_0 > 0 in Vpop for pre-existing resistance) |

## Derived quantities

```
C(t) = A_central / V              [mg/L]
C_µM(t) = C(t) × 1000 / MW       [µM]
TE(t) = KRAS_drug / (KRAS_GDP + KRAS_GTP + KRAS_drug)
```

## ODE System

### Module 1: PK (oral 1-compartment)

```
dA_gut/dt     = -ka × A_gut + dose_bolus(t)
dA_central/dt =  ka × A_gut - ke × A_central
```
where `ke = ln(2) / t_half` and `dose_bolus(t) = F × dose_mg` at each dosing time.

### Module 2: KRAS GDP/GTP cycling + covalent drug binding

```
dKRAS_GDP/dt  = k_hyd × KRAS_GTP - k_GEF × R × KRAS_GDP - k_bind × C_µM × KRAS_GDP + k_syn - k_deg × KRAS_GDP
dKRAS_GTP/dt  = k_GEF × R × KRAS_GDP - k_hyd × KRAS_GTP - k_deg × KRAS_GTP
dKRAS_drug/dt = k_bind × C_µM × KRAS_GDP - k_deg × KRAS_drug
```

### Module 3: Receptor feedback + pathway output (with bypass signaling)

```
R_target = (R_basal + G_fb × max(0, 1 - E/E_ss)) × (1 - f_block)
dR/dt = (1/tau_R) × (R_target - R)

E_drive = R × KRAS_GTP/KRAS_GTP_ss + beta_Z × Z
dE/dt = (1/tau_E) × (E_drive - E)
```

Note: `beta_Z × Z` represents bypass signaling — alternative pathway activation (e.g. NRAS, BRAF, PI3K) that
restores ERK-like output independently of KRAS-GTP, driven by the resistance state Z.

### Module 4: Tumor growth (saturating kill + resistance)

```
growth = rho_grow × (E/E_ss)^n_hill × T × (1 - T/T_max)

E_drug = max(0, 1 - E/E_ss)                    # drug-induced E suppression [0,1]
kill_sat = E_drug / (EC50_kill + E_drug)         # saturating drug effect
kill_resist = 1 / (1 + alpha_Z × Z)             # resistance attenuation
kill = rho_kill × kill_sat × kill_resist × T

dT/dt = growth - kill
```

### Module 6b: Adaptive resistance program

```
TE(t) = KRAS_drug / (KRAS_GDP + KRAS_GTP + KRAS_drug)    # target engagement [0,1]
dZ/dt = k_Z_up × TE × (1 - Z) - k_Z_down × Z
```

Z accumulates when drug is bound to target (TE > 0), reflecting selective pressure for resistance
(e.g. KRAS amplification, bypass RTK activation, MAPK reactivation via parallel pathways).
Decay is slow (`k_Z_down << k_Z_up`), making resistance partially irreversible on clinical timescales.

## Parameters

### PK
| Symbol | Value | Units | Source |
|--------|-------|-------|--------|
| dose_mg | 960 | mg | FDA label |
| tau_h | 24 | h | QD dosing |
| ka | 1.0 | /h | Tmax ~1-2h |
| t_half | 5.0 | h | FDA label |
| V | 211 | L | Vss/F, FDA label |
| F | 0.73 | — | FDA label |
| MW | 560.6 | g/mol | PubChem |

### KRAS cycling
| Symbol | Value | Units | Source |
|--------|-------|-------|--------|
| k_hyd | 0.2 | /h | Stites 2018 (G12C GTPase) |
| k_GEF | 0.025 | /h | Tuned for 11% GTP at SS |
| k_bind | 0.15 | /(µM·h) | Tuned for 77-92% TE |
| k_deg | 0.029 | /h | t1/2_protein = 24h |
| k_syn | 0.029 | µM/h | k_deg × KRAS_total |

### Signaling / feedback
| Symbol | NSCLC | CRC | Units | Source |
|--------|-------|-----|-------|--------|
| tau_R | 36 | 36 | h | RTK upregulation timescale |
| tau_E | 6 | 6 | h | ERK signaling dynamics |
| R_basal | 1.0 | 1.0 | — | Normalized |
| G_fb | 1.2 | 3.5 | — | Key lineage parameter |
| E_ss | 1.0 | 1.0 | — | Normalized |

### Tumor
| Symbol | NSCLC | CRC | Units | Source |
|--------|-------|-----|-------|--------|
| rho_grow | 0.0005 | 0.0007 | /h | Calibrated: DT ~58d (NSCLC), ~41d (CRC) |
| rho_kill | 0.0045 | 0.0045 | /h | Calibrated with saturating kill + resistance for nadir→regrowth |
| n_hill | 2 | 2 | — | Moderate nonlinearity |
| T_max | 5.0 | 5.0 | normalized | Carrying capacity (normalized scale) |
| T_0 | 1.0 | 1.0 | normalized | Initial tumor (normalized scale) |
| EC50_kill | 0.25 | 0.25 | — | Half-max E suppression for saturating kill |

### Combination
| Symbol | Value | Units | Source |
|--------|-------|-------|--------|
| f_block | 0.55 | — | EGFR blockade (panitumumab); tuned for CRC combo ORR |

### Resistance / Adaptation (M6b)
| Symbol | Value | Units | Source |
|--------|-------|-------|--------|
| k_Z_up | 0.004 | /h | Resistance induction rate; TE-driven |
| k_Z_down | 0.0002 | /h | Resistance decay rate; slow (partially irreversible) |
| alpha_Z | 4.0 | — | Resistance strength: kill / (1 + alpha_Z × Z) |
| beta_Z | 0.8 | — | Bypass signaling: Z restores E independently of KRAS-GTP |

## Steady State (no drug)

```
From dKRAS_GTP/dt = 0:  k_GEF × R × GDP_ss = (k_hyd + k_deg) × GTP_ss
=> GTP/GDP ratio = k_GEF × R_basal / (k_hyd + k_deg) = 0.025 / (0.2 + 0.029) ≈ 0.109
=> KRAS_GTP_ss / KRAS_total ≈ 9.8%
=> Z_ss = 0 (no drug → no TE → no resistance)
```

## Outputs / Endpoints

| Endpoint | Definition |
|----------|-----------|
| Target engagement | TE(t) = KRAS_drug / KRAS_total |
| ERK rebound | E(t) / E_ss at t = 72h |
| Best tumor change | min(T(t)/T_0 - 1) × 100% |
| RECIST response | Best change < -30% |
| Time to progression | First t where T(t) > T(t_nadir) × 1.2 |
| Nadir time | t at min(T(t)) |
| Resistance level | Z(t_final) — extent of adaptive resistance at end of simulation |
| Final tumor change | T(t_final)/T_0 - 1 — demonstrates regrowth past nadir |
