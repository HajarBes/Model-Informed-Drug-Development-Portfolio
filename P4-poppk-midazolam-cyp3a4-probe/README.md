# P4: Population PK of Oral Midazolam

**Industry-aligned PopPK workflow | nlmixr2 + rxode2 | OSP PBPK-calibrated**

Population pharmacokinetic analysis of oral midazolam (7.5 mg, single dose)
in 120 healthy adults. The dataset is simulated with parameters calibrated to
the [OSP Midazolam PBPK model](https://github.com/Open-Systems-Pharmacology/Midazolam-Model)
(validated against 40+ clinical studies), producing publication-realistic
concentration-time profiles. No real patient data are used.

The workflow follows the structure of the FDA PopPK Guidance (2022): model
development, diagnostics, covariate evaluation, simulation-based dose
assessment, and a CYP3A4 inhibition bridge to DDI risk (P10).

> **Start here:** [`report/EXECUTIVE_SUMMARY.md`](report/EXECUTIVE_SUMMARY.md)
> then [`report/REPORT_SUBMISSION_STYLE.md`](report/REPORT_SUBMISSION_STYLE.md)

---

## At a Glance

| | |
|---|---|
| **Model** | 2-compartment, first-order absorption, SAEM (dAIC = 467 vs 1-comp) |
| **CL/F** | 52 L/h estimated at 70 kg (true: 50) — consistent with CL ~25 L/h, F ~0.5 |
| **BSV on CL** | ~30% CV (eta shrinkage: 46%) |
| **Covariate** | Allometric WT on CL (^0.75) and V (^1.0); dOFV = 109, Vc BSV reduced 52% |
| **Weight effect** | Median shift: 50 kg +33%, 90 kg -14%, 100 kg -22% vs 70 kg reference |
| **DDI bridge** | 50% CL reduction (strong CYP3A4 inhibitor) yields 2.0x AUC fold-change |

---

## Why Midazolam?

Midazolam is the standard sensitive CYP3A4 probe substrate. Its PK is
dominated by CYP3A4-mediated first-pass metabolism in both the gut wall
and liver (oral F < 50%). This makes it the reference compound for DDI
studies and connects this analysis directly to P10 (Static DDI Framework):
the baseline PK characterized here is the same PK that perpetrator drugs
alter in a DDI scenario.

---

## Diagnostic Figures

Every figure answers a specific decision question:

| Figure | Question |
|--------|----------|
| [`gof_4panel.png`](figures/gof_4panel.png) | Does the model describe the data? |
| [`vpc.png`](figures/vpc.png) | Does the model capture population variability? |
| [`eta_distributions.png`](figures/eta_distributions.png) | Are random effects well-behaved? |
| [`individual_fits.png`](figures/individual_fits.png) | Does the model work per-subject? |
| [`forest_covariate_auc.png`](figures/forest_covariate_auc.png) | Are covariate effects clinically relevant? |
| [`wt_auc_ratio.png`](figures/wt_auc_ratio.png) | Is dose adjustment needed by weight? |
| [`exposure_dose_auc.png`](figures/exposure_dose_auc.png) | How does AUC vary by dose and weight? |
| [`cyp3a4_inhibition_fold_change.png`](figures/cyp3a4_inhibition_fold_change.png) | What AUC change under CYP3A4 inhibition? |

---

## Weight-Exposure Assessment (Verified)

The 0.80-1.25 band is used as a **clinical relevance heuristic**, not
formal bioequivalence criteria. Two metrics separate the weight effect
from individual variability (see [`aucr_definition.md`](outputs/tables/aucr_definition.md)):

**Stratum-level median shift** (isolates the weight effect):

| Weight | Median Shift | 90% Bootstrap CI | Within 0.80-1.25? |
|--------|-------------|------------------|-------------------|
| 50 kg | 1.33 | [1.29, 1.36] | No (above) |
| 70 kg | 1.00 | [0.98, 1.02] | Reference |
| 90 kg | 0.86 | [0.84, 0.88] | Yes |
| 100 kg | 0.78 | [0.76, 0.80] | No (below) |

Individual pharmacokinetic variability (BSV ~30% CV on CL) substantially
exceeds the weight-related shift at all strata. Even at 70 kg, 45% of
subjects fall outside the 0.80-1.25 band due to variability alone.

---

## Biology of the Model

**CL/F = 52 L/h** is apparent oral clearance. It reflects systemic
CL (~25 L/h) divided by oral bioavailability (~0.5). The model cannot
separate CL from F using oral data alone.

**Why 2-compartment?** Midazolam distributes rapidly into well-perfused
tissues (central) and more slowly into muscle/fat (peripheral). The
1-compartment model cannot capture this biphasic decline.

**What the random effects mean:**
- **eta(CL):** CYP3A4 metabolic capacity — enzyme expression, liver size, CYP3A5 genotype
- **eta(Ka):** Absorption rate — gastric emptying, intestinal motility
- **eta(Vc):** Distribution volume — body composition, protein binding

---

## Portfolio Connection

| Project | Role |
|---------|------|
| **P4 (this)** | Characterizes midazolam baseline PK |
| **P10 (Static DDI)** | Evaluates how CYP3A4 inhibitors alter midazolam AUC |
| **Future PBPK** | OSP model bridges PopPK to mechanistic modeling |

The CYP3A4 inhibition scenario (script 07) demonstrates the logic chain:
baseline PopPK quantifies CL/F, CYP3A4 inhibition reduces CL/F, and
the resulting ~2x AUC fold-change is consistent with P10's mechanistic
static DDI predictions.

---

## Quick Start

```bash
cd P4-poppk-midazolam-cyp3a4-probe
bash run_all.sh    # 7 scripts, ~10 min total
```

Step 1 (simulation) needs only base R + deSolve. Steps 2-7 require
nlmixr2/rxode2.

## Project Structure

```
P4-poppk-midazolam-cyp3a4-probe/
├── run_all.sh                              # Single reproducible entrypoint
├── analysis/
│   ├── 00_setup.R                          # Packages, paths, true parameters
│   ├── 01_simulate_trial_dataset.R         # 120-subject OSP-calibrated trial
│   ├── 02_fit_models.R                     # 1-comp vs 2-comp SAEM
│   ├── 03_model_diagnostics.R              # GOF, VPC, etas, individual fits
│   ├── 04_covariates.R                     # Allometric WT + forest plot
│   ├── 05_simulation_decision.R            # Dose scenarios + recommendation
│   ├── 06_verify_aucr.R                    # AUCR verification + bootstrap CI
│   └── 07_cyp3a4_inhibition_bridge.R       # CYP3A4 sensitivity → P10 bridge
├── data/
│   ├── simulated/                          # NONMEM-format dataset + dictionary
│   └── sources/                            # OSP calibration notes + citations
├── figures/                                # 10 decision-oriented figures
├── outputs/tables/                         # Estimates, verified AUCR, DDI summary
└── report/
    ├── EXECUTIVE_SUMMARY.md                # 1-page decision summary
    ├── REPORT_SUBMISSION_STYLE.md          # 10-section guidance-aligned report
    ├── METHODS.md                          # Equations, estimation, software
    └── LIMITATIONS_AND_SCOPE.md            # Scope + what a real analysis needs
```

## Traceability

| Item | Detail |
|------|--------|
| R | 4.5.0 |
| nlmixr2 / rxode2 | 5.0.0 / 5.0.1 |
| Random seeds | Documented per script (42, 20240201, 20240315) |
| Scope | Exploratory methodological demonstration; not GxP-validated |

**File provenance** — each script's outputs are documented in
[`report/METHODS.md`](report/METHODS.md) and the provenance table below:

| Script | Key Outputs |
|--------|-------------|
| `01_simulate_trial_dataset.R` | `analysis_dataset.csv`, `qc_spaghetti_log.png` |
| `02_fit_models.R` | `model_comparison.csv`, `parameter_estimates_base.csv` |
| `03_model_diagnostics.R` | `gof_4panel.png`, `vpc.png`, `eta_distributions.png`, `eta_shrinkage.csv` |
| `04_covariates.R` | `forest_covariate_auc.png`, `bsv_comparison.csv` |
| `05_simulation_decision.R` | `exposure_dose_auc.png`, `wt_auc_ratio.png` |
| `06_verify_aucr.R` | `aucr_stratum_shift_summary.csv`, `aucr_definition.md` |
| `07_cyp3a4_inhibition_bridge.R` | `cyp3a4_inhibition_fold_change.png` |

## References

- OSP Midazolam Model v2.0. [github.com/Open-Systems-Pharmacology/Midazolam-Model](https://github.com/Open-Systems-Pharmacology/Midazolam-Model)
- Hanke N et al. CPT Pharmacometrics Syst Pharmacol 2018;7(10):647-659
- Anderson BJ, Holford NHG. Annu Rev Pharmacol Toxicol 2008;48:303-32
- FDA (2022). Population Pharmacokinetics. Guidance for Industry.
