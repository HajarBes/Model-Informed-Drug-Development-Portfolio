# Assumptions — Project 05: Sotorasib Exposure-Response

## PK Simulation Assumptions

1. **Direct exposure sampling:** Virtual patients are simulated by drawing AUCss and Cmax from log-normal distributions centered on published geometric means at each dose level, with published IIV (CV% for CL/F and V/F). This is the standard E-R approach when full popPK model code is unavailable.

2. **Day 8 → Steady-state scaling:** Cross-dose PK comparison uses Day 8 data (only available multi-dose timepoint). Scaled to steady state using the 960 mg Day 8→SS ratio (factor ≈ 2.0). This assumes similar accumulation dynamics across doses.

3. **240 mg PK interpolation:** No published Day 8 PK data at 240 mg. Exposure targets interpolated from 180 mg and 360 mg data. The dose comparison study confirms similar PFS at 240 and 960 mg, consistent with similar exposure.

4. **IIV magnitude:** AUC IIV based on published CL/F CV (76%). Cmax IIV set at 55% (moderate, reflecting both CL and V variability). AUC-Cmax correlation set at 0.65 (shared CL contribution).

## Efficacy Model Assumptions

5. **Emax logistic model:** ORR modeled as P(response) = E0 + Emax × AUC^γ / (EC50^γ + AUC^γ). Parameters chosen to produce ~30% ORR across the flat AUC range.

6. **Flat ORR by design:** Because exposure is approximately flat across 180–960 mg, the model correctly predicts flat ORR. Published dose-level ORR variability (25–50%) is driven by small sample sizes (N=3 to N=34) and not captured by the model.

7. **Phase 2 as primary anchor:** Simulated ORR calibrated to the Phase 2 estimate (36%, N=124) and dose comparison average (~29%). Small-N Phase 1 cohorts used as secondary reference only.

## Safety Model Assumptions

8. **Logistic safety model:** P(Grade 3+ hepatotox) = expit(α + β_cmax × Cmax + β_cpi × CPI). Parameters calibrated to 16% at 960 mg and ~14% at 240 mg.

9. **CPI interaction modeled as binary:** Prior checkpoint inhibitor use modeled as 0/1. Published data shows a timing effect (≤30 days = 75%, >90 days = 0%), but our model uses a single binary flag. The simulated CPI effect captures the average risk elevation.

10. **No mechanistic claims:** Hepatotoxicity is modeled as exposure-related and CPI-modified based on clinical observations. No assertion about underlying mechanism (covalent reactivity, immune-mediated hepatitis, etc.).

## Covariate Analysis Assumptions

11. **Published effects applied directly:** ECOG and albumin effects on CL are applied as scalar adjustments (ECOG 1: +15% CL, ECOG 2: +35% CL, low albumin: +25% CL). These magnitudes are approximate based on published popPK.

12. **Non-significant covariates:** Weight, renal function, prior CPI, and hepatic function are confirmed non-significant on PK based on published popPK results. Forest plot shows these covariates with ratio ≈ 1.0.

## General Assumptions

13. **No individual patient data:** All analyses use published aggregate endpoints. Within-dose exposure quartile stratifications are model-implied from the virtual population.

14. **PFS descriptive only:** Published median PFS reported without simulated time-to-event analysis. Insufficient basis for Weibull/Cox modeling without individual KM data.

15. **Dose comparison not powered:** 240 mg vs 960 mg results are numerically descriptive, not statistically hypothesis-tested in our analysis.
