# Report — Sotorasib Exposure-Response Analysis
### Structured Along FDA E-R Guidance Sections

---

## 1. Objectives

**Primary:** Characterize the exposure-response relationship for sotorasib across efficacy (ORR, PFS) and safety (Grade 3+ hepatotoxicity) endpoints, and evaluate whether 960 mg QD is the optimal dose compared to 240 mg QD.

**Secondary:** Identify patient-level factors (covariates) that drive exposure variability and assess whether exposure-based risk stratification can improve the therapeutic index.

## 2. Data Sources & Methodological Transparency

All analyses are based on published clinical data from peer-reviewed sources and FDA regulatory reviews. This project demonstrates an end-to-end E-R workflow using public, aggregate clinical results and a virtual population calibrated to published PK summaries (no IPD). All within-arm stratifications are model-implied for illustration.

| Data Source | Type | Reference |
|-------------|------|-----------|
| Population PK parameters | Summary statistics | Nagase et al., AAPS J 2025; FDA label |
| Efficacy by dose (ORR) | Aggregate dose-level | CodeBreaK 100 Phase 1/2 |
| PFS by dose | Published medians | Dose comparison, EJC 2024 |
| Safety (hepatotoxicity) | Aggregate rates | FDA clinical review; pooled safety |
| Dose comparison (240 vs 960) | Published endpoints | EJC 2024 |

## 3. Population PK Summary

### Published Model
2-compartment model with 3 transit absorption compartments (Nagase et al. 2025). Saturable absorption produces less-than-dose-proportional exposure across 180–960 mg.

### Reduced Model for E-R Analysis
Direct exposure sampling from log-normal distributions anchored to published geometric means. Justification: E-R analyses require steady-state exposure metrics (AUCss, Cmax,ss), which are well-captured by sampling from correct marginal distributions without requiring full PK profile simulation.

### PK Calibration (Anchor vs Simulated)

| Dose (mg) | Published GM AUC | Simulated GM AUC | Ratio |
|-----------|-----------------|------------------|-------|
| 180 | 63.9 | ~67 | 1.05 |
| 240 | 71.1* | ~73 | 1.03 |
| 360 | 78.4 | ~81 | 1.03 |
| 720 | 84.8 | ~86 | 1.01 |
| 960 | 65.3 | ~65 | 1.00 |

*240 mg interpolated from 180 and 360 mg data. All values in hr·µg/mL.

## 4. Exposure Metrics

| Metric | Units | 960 mg SS | Source |
|--------|-------|-----------|--------|
| AUC0-24h | hr·µg/mL | 65.3 | FDA label |
| Cmax | µg/mL | 7.50 | FDA label |
| Cmax (µM) | µM | 13.4 | Converted (MW 560.6) |
| CL/F | L/hr | 26.2 (CV 76%) | PopPK |
| t½ | hr | 5.0 | FDA label |

## 5. Exposure-Efficacy: ORR

### Dose-Level ORR

| Dose (mg) | Published ORR | Simulated ORR | N (published) |
|-----------|--------------|---------------|---------------|
| 180 | 33% | ~29% | 3 |
| 240 | 24.8% | ~33% | 105 |
| 360 | 25% | ~33% | 16 |
| 720 | 50% | ~30% | 6 |
| 960 | 36% | ~31% | 124 |

### Key Finding
Within the observed exposure range (AUC ~64–85 hr·µg/mL across 180–960 mg), the incremental change in predicted ORR is small; therefore dose changes do not translate to meaningful response changes. This reflects the compressed exposure range (not drug inactivity): saturable absorption places all doses on the plateau of the Emax E-R curve, limiting the achievable exposure gradient.

### Exposure Quartile Analysis (960 mg, model-implied, virtual patients)
Shallow positive trend: Q1 ORR ~22%, Q4 ORR ~46%. The narrow AUC range within a dose level limits the observable E-R gradient.

## 6. Exposure-Efficacy: PFS (Descriptive)

| Dose | Median PFS | N | Source |
|------|-----------|---|--------|
| 960 mg | 5.4 months | 104 | Dose comparison |
| 240 mg | 5.6 months | 105 | Dose comparison |

Implied HR = 1.04. No simulated time-to-event analysis performed (insufficient basis without IPD).

## 7. Exposure-Safety: Hepatotoxicity

### Safety Calibration (Anchor vs Simulated)

| Dose (mg) | Published Rate | Simulated Rate | Ratio |
|-----------|---------------|----------------|-------|
| 240 | 14% | ~17% | 1.20 |
| 960 | 16% | ~16% | 1.00 |

### CPI Interaction
Prior checkpoint inhibitor use increases hepatotoxicity risk approximately 3-fold. Published data shows extreme time-dependence: 75% rate within 30 days of CPI vs 0% beyond 90 days.

### Risk Stratification (model-implied, virtual patients)
The 2x2 matrix (low/high Cmax × CPI−/CPI+) identifies a high-risk subgroup (High Cmax + Prior CPI) with hepatotoxicity rates substantially above the population average.

## 8. Dose Optimization

### 240 mg vs 960 mg Benefit-Risk Summary

| Metric | 240 mg | 960 mg |
|--------|--------|--------|
| GM AUCss (hr·µg/mL) | ~73 | ~65 |
| ORR (simulated) | ~30% | ~31% |
| ORR (published) | 24.8% | 32.7% |
| PFS (published) | 5.6 mo | 5.4 mo |
| Grade 3+ Hepatotox | ~17% | ~16% |
| Therapeutic Index | ~1.8 | ~2.0 |

### Interpretation
Given saturable absorption, higher dose does not guarantee higher exposure; the benefit-risk profiles at 240 mg and 960 mg overlap substantially. A dose-optimization study is justified, and benefit-risk may depend on subgroups (e.g., patients with prior CPI and high Cmax). This analysis does not definitively conclude that one dose is superior — it demonstrates that the E-R data support formal evaluation of lower doses.

## 9. Covariate Effects on Exposure

### Significant Covariates
- **ECOG 1/2:** Increased CL → reduced AUC (ratio 0.84/0.78 vs ECOG 0)
- **Low albumin (≤3.5 g/dL):** Increased CL → reduced AUC (ratio 0.81)

### Not Significant
Body weight, renal function, hepatic function (mild), prior CPI — all within 0.80–1.25 bioequivalence bounds.

## 10. Discussion & Regulatory Implications

### Summary
Sotorasib's compressed dose-response is a pharmacokinetic phenomenon, not a pharmacodynamic one. Saturable absorption creates a ceiling on achievable exposure, placing all dose levels (180–960 mg) in a narrow region of the E-R curve. The hepatotoxicity signal, while modest, is Cmax-related and substantially modified by prior CPI timing. Given this pharmacology, a dose-optimization study is justified, and benefit-risk may depend on subgroups (e.g., CPI timing + high Cmax).

### Regulatory Alignment
This analysis is consistent with:
1. **FDA's dose optimization PMR** for formal 240 mg vs 960 mg comparison
2. **ODAC's 10-2 vote** recommending further evaluation of lower doses
3. **ICH E4 principles** that dose-response should inform dose selection
4. **Project Optimus** goals of selecting the optimal (not maximal) dose

### Limitations
See LIMITATIONS_AND_SCOPE.md for full discussion. Key: no IPD, reduced PK model, aggregate anchoring, safety model is clinical (not mechanistic), PFS descriptive only, methodological demonstration (not regulatory filing).
