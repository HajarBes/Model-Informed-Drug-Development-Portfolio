# Limitations and Scope — Sotorasib Exposure-Response & Dose Optimization

## Scope

This is a **methodological demonstration** of exposure-response analysis and dose optimization using published clinical data for sotorasib. It is designed as a portfolio project to illustrate quantitative pharmacology skills in the context of a real and ongoing regulatory debate. It is not a regulatory filing and does not use proprietary data.

## Limitations

### 1. No Individual Patient-Level Data (IPD)
Virtual patients are simulated from published popPK summary statistics (geometric means, CVs). All within-dose exposure quartile stratifications are **model-implied from the virtual population**, not derived from real patient-level quartile data. This limits the ability to capture individual-level exposure-outcome associations and patient-specific covariate effects.

### 2. Reduced PK Approach
The published popPK model is a 2-compartment structure with 3 transit absorption compartments (Nagase et al. 2025). Our approach uses direct exposure sampling from log-normal distributions (not a mechanistic PK model). This is sufficient for steady-state AUC/Cmax metrics used in E-R analysis but does not capture full absorption dynamics, time-course profiles, or individual PK trajectories.

### 3. Aggregate Anchoring
E-R relationships are calibrated to published dose-level aggregate endpoints (ORR, Grade 3+ rates), not individual exposure-outcome pairs. This means:
- The simulated E-R gradient within a dose level is model-implied, not data-driven
- True individual-level E-R may differ from the aggregate relationship (ecological fallacy)
- The Emax model shape is assumed, not estimated from individual data

### 4. Safety Model Is Clinical, Not Mechanistic
Hepatotoxicity is modeled as exposure-related (Cmax coefficient) and modified by prior checkpoint inhibitor use (binary CPI flag). No mechanistic assertion is made about the underlying biology:
- Covalent reactivity of the acrylamide warhead with off-target cysteines
- Immune-mediated hepatitis from CPI priming
- Direct hepatocellular toxicity
- Idiosyncratic drug reaction

The model captures clinical risk factors, not mechanisms.

### 5. PFS Descriptive Only
Published median PFS values are reported (5.4 mo at 960 mg, 5.6 mo at 240 mg). No simulated time-to-event analysis is performed because:
- Individual-level KM data are not available
- Weibull/Cox modeling from aggregate medians alone would be under-determined
- The similar median PFS values are sufficient to support the dose-optimization conclusion

### 6. Dose Comparison Study Not Powered
The 240 mg vs 960 mg dose comparison (N=104 vs N=105) was not powered as a formal non-inferiority or superiority trial. ORR and PFS differences are numerically descriptive; statistical hypothesis testing was not conducted in this analysis.

### 7. Small-N Phase 1 Cohorts
Published ORR at 180 mg (N=3), 360 mg (N=16), and 720 mg (N=6) have wide 95% confidence intervals. The Emax model cannot and should not reproduce this sampling noise. Simulated ORR is calibrated to the larger Phase 2 and dose comparison datasets.

### 8. Covariate Effects Are Applied, Not Estimated
ECOG and albumin effects on clearance are applied as scalar adjustments based on published popPK results, not estimated de novo from the virtual population. The forest plot is a visualization of published findings, not an independent covariate analysis.

### 9. CPI Timing Not Modeled
The published CPI-hepatotoxicity interaction shows dramatic timing dependence (75% within 30 days, 0% beyond 90 days). Our model uses a binary CPI flag (yes/no), which captures average risk elevation but not the temporal gradient. A more refined analysis would require individual-level CPI timing data.

### 10. Methodological Demonstration
This project demonstrates E-R methodology for portfolio purposes. Findings are directionally consistent with published analyses and regulatory advisory committee recommendations but should not be interpreted as a definitive clinical recommendation. The analysis supports the rationale for a dose-optimization study — it does not conclude that one dose is superior to another.
