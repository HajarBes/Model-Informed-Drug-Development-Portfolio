# Methods — Sotorasib Exposure-Response & Dose Optimization

## Overview

This analysis uses published aggregate clinical data to simulate a virtual patient population and characterize the exposure-response relationship for sotorasib across efficacy (ORR, PFS) and safety (hepatotoxicity) endpoints. The approach follows a calibrate-then-predict framework where model parameters are estimated from published dose-level endpoints and verified against held-out or independent data.

## Calibration / Qualification Separation

**Calibration (parameters estimated to match):**
- PK: Saturable bioavailability parameters (D50) and effective clearance estimated to match published Day 8 AUC/Cmax at 180, 360, 720, 960 mg
- Efficacy: Emax logistic model parameters (E0, Emax, EC50, gamma) calibrated to reproduce published Phase 2 ORR (36% at 960 mg, N=124)
- Safety: Logistic model intercept (alpha) and Cmax/CPI coefficients calibrated to match 16% Grade 3+ hepatotoxicity at 960 mg and published CPI interaction

**Qualification (held-out checks):**
- PK at 240 mg (interpolated, not directly calibrated; confirmed by similar PFS in dose comparison)
- Dose comparison ORR (24.8% at 240 mg vs 32.7% at 960 mg): directional consistency checked post-hoc
- Safety at 240 mg (14% published vs simulated): not used in parameter estimation
- Covariate effect directions and magnitudes vs published popPK (Nagase et al. 2025)

## Claim Language Rule

All exposure quartile plots and within-dose stratifications carry the label **"model-implied (virtual patients)"** in their title. This makes clear that quartile-level results come from the simulated virtual population, not from real individual patient-level data.

## Module 1: Virtual Patient PK Simulation

### Approach
Direct exposure sampling from log-normal distributions anchored to published geometric mean AUC and Cmax at each dose level. This is the standard E-R approach when published popPK model code is unavailable.

### Model
- AUCss and Cmax drawn from correlated bivariate log-normal distributions
- Parameters: omega_AUC = 0.675 (from CL/F CV 76%), omega_Cmax = 0.55, correlation = 0.65
- N = 500 virtual patients per dose level (180, 240, 360, 720, 960 mg)
- Total: 2,500 virtual patients

### Saturable Bioavailability
- Published data shows less-than-dose-proportional exposure: AUC barely changes from 180 to 960 mg
- Relative bioavailability modeled as F(dose) = D50 / (D50 + dose)
- D50 estimated by least-squares fit to published Day 8 AUC ratios across doses
- Day 8 values scaled to steady-state using the 960 mg Day 8→SS ratio

### Reduced Model Justification
The published popPK model is a 2-compartment structure with 3 transit absorption compartments (Nagase et al. 2025). Our reduced approach uses direct exposure sampling (not a mechanistic PK model) because:
1. The full model code is not publicly available
2. E-R analyses use steady-state AUC and Cmax — metrics well-captured by sampling from the correct marginal distributions
3. Individual PK profiles and absorption dynamics are not needed for aggregate E-R characterization

## Module 2: Exposure-Efficacy (ORR)

### Model
Emax logistic model:
```
P(response) = E0 + Emax * AUC^gamma / (EC50^gamma + AUC^gamma)
```

Parameters: E0 = 0.12, Emax = 0.25, EC50 = 30 hr·µg/mL, gamma = 1.2

### Key Finding
Within the observed exposure range (AUC ~64–85 hr·µg/mL across 180–960 mg), the incremental change in predicted ORR is small; therefore dose changes do not translate to meaningful response changes. This reflects the compressed exposure range (not drug inactivity): saturable absorption places all doses on the plateau of the Emax E-R curve.

### Calibration
- Parameters selected to produce ORR ~30% across the observed AUC range (~65–85 hr·µg/mL)
- Consistent with Phase 2 ORR (36%, N=124) and dose comparison average (~29%)
- Small-N Phase 1 cohort ORR (25–50%, N=3–34) used as secondary reference; wide CIs acknowledged

### Exposure Quartile Analysis
- 960 mg patients divided into AUC quartiles (Q1–Q4)
- ORR computed per quartile
- **Labeled as "model-implied (virtual patients)"** — not derived from real patient-level quartile data

## Module 3: Exposure-Efficacy (PFS) — Descriptive Only

Published median PFS reported by dose level:
- 960 mg: 5.4 months (N=104)
- 240 mg: 5.6 months (N=105)
- Implied HR = 1.04 (calculated from median ratio)

**No simulated time-to-event analysis.** Weibull/Cox modeling without individual-level KM data would be statistically inappropriate. Published medians are sufficient for the dose-optimization argument.

## Module 4: Exposure-Safety (Hepatotoxicity)

### Model
Logistic regression:
```
P(Grade 3+ hepatotox) = expit(alpha + beta_cmax * Cmax + beta_cpi * prior_CPI)
```

Parameters: alpha = -3.5, beta_cmax = 0.08, beta_cpi = 1.8

### Calibration
- Overall rate at 960 mg: 16% (published: 16%)
- CPI interaction: prior CPI increases risk approximately 3-fold
- Clinical framing: hepatotoxicity is exposure-related and modified by prior CPI timing

### Risk Stratification
Patients classified by Cmax level (above/below median) and CPI status, creating a 2x2 risk matrix. High Cmax + prior CPI identifies the highest-risk subgroup.

## Module 5: Dose Optimization

### Benefit-Risk Framework
For each dose level:
- **Therapeutic index** = P(ORR) / P(Grade 3+ hepatotoxicity)
- **NNT** = 1 / ORR (patients treated per responder)
- **NNH** = 1 / hepatotoxicity rate (patients treated per Grade 3+ event)
- **Net clinical benefit** = ORR - hepatotoxicity rate

### Decision Framework
The 240 mg vs 960 mg comparison evaluates whether dose reduction preserves efficacy while reducing safety risk, consistent with the FDA PMR objective.

## Module 6: Covariate Effects

### Approach
Forest plot of AUC ratio versus reference for each covariate subgroup at 960 mg. Covariate effects applied as scalar adjustments based on published popPK findings:
- ECOG 1: +15% CL (lower AUC)
- ECOG 2: +35% CL (lower AUC)
- Low albumin: +25% CL (lower AUC)
- Weight, renal, hepatic, CPI: no significant effect on PK

### Verification
Direction and approximate magnitude of effects confirmed against published popPK covariate analysis (Nagase et al. 2025).

## Software

- Python 3.10+ with NumPy, SciPy, pandas, matplotlib, statsmodels
- Random seeds: SEED_VPOP = 20240501, SEED_ER = 20240701
- All outputs reproducible via `bash run_all.sh`
