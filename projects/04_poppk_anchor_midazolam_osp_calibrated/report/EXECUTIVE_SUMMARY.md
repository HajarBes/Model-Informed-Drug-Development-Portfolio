# Executive Summary — P4 PopPK Analysis of Oral Midazolam

## Decision Question

Does body weight require dose adjustment for oral midazolam 7.5 mg in
healthy adults, and how does baseline PK characterization inform CYP3A4
drug-drug interaction risk?

## Evidence

A 2-compartment PopPK model with first-order absorption was fitted to
a simulated dataset (N=120, publication-calibrated to the OSP Midazolam
PBPK model) using SAEM estimation in nlmixr2/rxode2.

**Model selection:** 2-compartment selected over 1-compartment (dAIC = 467),
consistent with known biphasic midazolam disposition.

**Key parameter estimates (70 kg reference):**
- CL/F = 52 L/h (consistent with CL ~25 L/h, F ~0.5)
- Vc/F = 47 L, Vp/F = 63 L
- BSV on CL: ~30% CV

**Covariate model:** Allometric WT scaling (CL ~ WT^0.75, V ~ WT^1.0)
improved fit by dOFV = 109 and reduced Vc BSV by 52%.

**Verified weight-exposure analysis (7.5 mg, stratum median shift vs 70 kg):**

| Weight | Median Shift | 90% Bootstrap CI | Within 0.80-1.25? |
|--------|-------------|------------------|-------------------|
| 50 kg | 1.33 | [1.29, 1.36] | No (above) |
| 70 kg | 1.00 | [0.98, 1.02] | Yes (reference) |
| 90 kg | 0.86 | [0.84, 0.88] | Yes |
| 100 kg | 0.78 | [0.76, 0.80] | No (below) |

These shifts isolate the weight effect on typical exposure. Individual
variability from BSV in CL (~30% CV) is substantially larger: even at
70 kg, 45% of subjects fall outside 0.80-1.25 due to pharmacokinetic
variability alone.

**CYP3A4 sensitivity:** 50% CL reduction → 2.0x AUC fold-change,
consistent with P10 static DDI predictions for strong CYP3A4 inhibitors.

## Conclusion

Allometric WT scaling is statistically and physiologically justified.
Median weight-exposure shifts are modest but marginally exceed the
0.80-1.25 clinical relevance heuristic at the extremes (50 kg, 100+ kg).
For midazolam's sedation indication (titration-to-effect), fixed 7.5 mg
dosing is generally acceptable. In contexts with narrow therapeutic
windows, weight-based adjustment may warrant consideration.

The baseline PopPK characterization directly informs DDI risk: CYP3A4
inhibition-mediated AUC changes are predictable from CL/F and fm, linking
this analysis to P10's mechanistic static DDI framework.

## Scope Statement

This is a **methodological demonstration** using simulated data calibrated
to published PK evidence. It is not a clinical study analysis or regulatory
filing. See `LIMITATIONS_AND_SCOPE.md`.
