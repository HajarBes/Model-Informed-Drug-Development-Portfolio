# Static DDI Risk Assessment Framework

**Mechanistic static model for drug-drug interaction screening per FDA 2020 guidance**

## Overview

This project implements a complete static DDI risk assessment framework following the FDA 2020 In Vitro DDI Guidance. It applies the mechanistic static model equations to predict DDI magnitude (AUCR) and classify interaction severity, with sensitivity analysis to characterize prediction robustness.

Two cases demonstrate the framework:

| Case | Perpetrator | Victim | Key Finding |
|------|------------|--------|-------------|
| **A (Benchmark)** | Ketoconazole 400 mg | Midazolam | AUCR = 27.6 (Strong), observed = 11.2, validates framework |
| **B (Signature)** | Sotorasib 960 mg | General assessment | Dual TDI + induction on CYP3A4, 6/9 transporters flagged |

## Quick Start

```bash
cd 03-static-DDI-framework
bash run_all.sh
```

## Project Structure

```
03-static-DDI-framework/
├── README.md
├── run_all.sh                          # One-command pipeline execution
├── data/
│   ├── assumptions.md                  # Parameter documentation with sources
│   ├── inputs_midazolam_case.csv       # Case A input parameters
│   ├── inputs_sotorasib_case.csv       # Case B input parameters
│   └── processed/                      # Intermediate RDS files (auto-generated)
├── analysis/
│   ├── 00_setup.R                      # Packages, paths, DDI computation functions
│   ├── 01_compute_static_ddi.R         # CYP/transporter screening + AUCR prediction
│   ├── 02_sensitivity_uncertainty.R    # Tornado, heatmap, Monte Carlo, Igut check
│   ├── 03_make_figures.R              # 5 publication-quality figures
│   └── 04_generate_tables.R           # Formatted summary tables
├── figures/                            # Output PNGs (auto-generated)
├── outputs/
│   ├── tables/                         # CSV results (auto-generated)
│   └── logs/                           # Run logs + Igut sanity check
└── report/
    ├── REPORT.md                       # Comprehensive cross-case report
    ├── CASE_A_midazolam.md             # Benchmark case report
    ├── CASE_B_sotorasib.md             # Signature case report
    └── METHODS_APPENDIX.md             # Full mathematical derivations
```

## Methods

### DDI Screening (per FDA 2020 Guidance)

| Mechanism | Equation | Threshold |
|-----------|----------|-----------|
| Reversible CYP inhibition | R1 = 1 + [I]max,u / Ki | R1 >= 1.02 |
| Intestinal CYP inhibition | R1,gut = 1 + [I]gut / Ki | R1,gut >= 11 |
| Time-dependent inhibition | TDI factor = kdeg / (kdeg + kinact*[I]u/(KI+[I]u)) | TDI factor < 0.5 |
| CYP induction | R3 = 1 / (1 + d*Emax*[I]u/(EC50+[I]u)) | R3 <= 0.8 |
| Transporter (systemic) | [I]max,u / IC50 | >= 0.1 |
| Transporter (intestinal) | [I]gut / IC50 | >= 10 |

### AUCR Prediction

Total AUCR = AUCR_hepatic x AUCR_gut

- AUCR_hepatic = 1 / (fm * (1/R1) * TDI_factor + (1-fm))
- AUCR_gut = 1 / (fg * (1/R1,gut) + (1-fg))

### Classification (FDA)

| AUCR | Classification |
|------|---------------|
| >= 5 | Strong |
| >= 2 | Moderate |
| >= 1.25 | Weak |
| < 1.25 | No interaction |

## Key Results

### Case A: Midazolam + Ketoconazole
- Predicted AUCR = 27.6 (Strong), Observed = 11.2 (Strong)
- Pred/Obs ratio = 2.46 (expected conservative overprediction)
- Monte Carlo: 100% of 10,000 samples classify as Strong
- Most sensitive to fg and fm; robust to inhibitor parameter uncertainty

### Case B: Sotorasib
- CYP3A4: Both TDI (TDI factor = 0.027) and induction (R3 = 0.23) flagged
- CYP2C8: R1 = 2.93 (clinically confirmed inhibition)
- Transporters: 6/9 flagged (P-gp, BCRP, OATP1B1/1B3, MATE1/2-K)
- Static model predicts net inhibition; clinical data shows net induction (AUCR = 0.48)
- Demonstrates the need for PBPK to resolve opposing CYP3A4 mechanisms

## Figures

| # | Figure | Description |
|---|--------|-------------|
| 1 | `midazolam_tornado.png` | Tornado sensitivity: parameter impact on AUCR |
| 2 | `midazolam_heatmap_Ki_vs_Imaxu.png` | Decision boundary: Ki vs [I]max,u |
| 3 | `midazolam_pred_obs_band.png` | Monte Carlo AUCR distribution vs observed |
| 4 | `sotorasib_mechanism_dashboard.png` | 3-panel: CYP inhibition, induction, transporters |
| 5 | `sotorasib_net_effect_scenarios.png` | Net CYP3A4 effect: TDI vs induction scenarios |

## References

- FDA (2020). In Vitro Drug Interaction Studies — Cytochrome P450 Enzyme- and Transporter-Mediated Drug Interactions. Guidance for Industry.
- Obach RS et al. (2007). Mechanism-based inactivation of human cytochrome P450 enzymes and the prediction of drug-drug interactions. Drug Metab Dispos.
- Huang SM et al. (2007). New era in drug interaction evaluation. Clin Pharmacol Ther.
- Yang J et al. (2008). Cytochrome P450 turnover. Drug Metab Dispos.
- LUMAKRAS (sotorasib) NDA 214665, FDA Clinical Pharmacology Review.

## Requirements

- R >= 4.1
- Packages: tidyverse, readr, dplyr, ggplot2, scales, patchwork (auto-installed)
