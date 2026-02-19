# Executive Summary — Sotorasib Exposure-Response & Dose Optimization

## Decision Question

**Does sotorasib exhibit a clinically meaningful exposure-response relationship for efficacy (ORR, PFS) or safety (hepatotoxicity), and do the data support 960 mg QD as the optimal dose over 240 mg?**

## Key Findings

### 1. Exposure Is Flat Across 180–960 mg (Saturable Absorption)
Sotorasib displays saturable oral absorption: geometric mean AUCss ranges from approximately 64 to 85 hr·µg/mL across the entire 180–960 mg dose range — a less than 1.4-fold change over a 5.3-fold dose range. Cmax,ss similarly plateaus (~7.5–9.0 µg/mL). This is the root cause of the flat dose-response relationship.

### 2. Compressed E-R Range — Dose Changes Don't Change Response
Within the observed exposure range (AUC ~64–85 hr·µg/mL across 180–960 mg), the incremental change in predicted ORR is small; therefore dose changes do not translate to meaningful response changes. The Emax logistic model produces ORR of ~30% across all dose levels, consistent with the Phase 2 BICR estimate (36%, N=124) and dose comparison results (24.8–32.7%). This is driven by saturable absorption compressing the exposure range, not by drug inactivity — all doses place patients on the plateau of the E-R curve.

### 3. Hepatotoxicity Is Cmax-Related and CPI-Modified
Grade 3+ hepatotoxicity rates (~16% at 960 mg, ~14% at 240 mg) show a modest Cmax-driven component modified substantially by prior checkpoint inhibitor (CPI) use. Patients with prior CPI show hepatotoxicity rates approximately 3x higher than CPI-naive patients. This interaction is clinically actionable for risk stratification.

### 4. Dose-Optimization Decision Remains Open
- Similar ORR (~25–33% published; ~30% simulated at both doses)
- Near-identical PFS (5.6 vs 5.4 months)
- Comparable hepatotoxicity rates (14% vs 16%)
- Therapeutic index is similar across doses

## Interpretation

Given saturable absorption, higher dose does not guarantee higher exposure; a dose-optimization study is justified, and benefit-risk may depend on subgroups (e.g., CPI timing + high Cmax). This analysis supports the FDA's dose optimization post-marketing requirement (PMR) and aligns with the ODAC 10-2 vote recommending further evaluation of lower doses.

Sotorasib is an active agent — the question is not whether the drug works, but whether 960 mg QD provides clinically meaningful benefit over lower doses that achieve similar exposure.

## Methodological Note

This project demonstrates an end-to-end E-R workflow using public, aggregate clinical results and a virtual population calibrated to published PK summaries (no IPD). All within-arm stratifications are model-implied for illustration. See METHODS.md and LIMITATIONS_AND_SCOPE.md for details.
