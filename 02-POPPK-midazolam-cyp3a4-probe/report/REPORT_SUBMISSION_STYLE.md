# Population Pharmacokinetic Analysis of Oral Midazolam
## Guidance-Aligned Analysis Report

---

## 1. Context of Use

This PopPK analysis characterizes the population pharmacokinetics of oral
midazolam following single-dose administration (7.5 mg tablet) in healthy
adults. The analysis addresses the following decision questions:

- Structural model selection (1-compartment vs 2-compartment)
- Quantification of covariate effects on apparent clearance and volume
- Simulation-based exposure evaluation across body weight strata
- Sensitivity analysis linking baseline PK to CYP3A4 inhibition risk

The analysis uses a publication-calibrated simulated dataset anchored to
the OSP Midazolam PBPK model (validated across 40+ studies). The workflow
follows the structure and conventions of the FDA PopPK Guidance (2022).

## 2. Data

### 2.1 Dataset

| Attribute | Value |
|-----------|-------|
| Source | Simulated, calibrated to OSP Midazolam PBPK (v2.0) |
| Subjects | 120 (40 rich + 80 sparse sampling) |
| Dose | 7.5 mg oral midazolam (single dose) |
| Observations | ~700 total (above LLOQ) |
| LLOQ | 0.5 ng/mL |
| BLQ handling | M1 (excluded with MDV=1) |

### 2.2 Calibration Anchor

True simulation parameters produce concentration-time profiles consistent
with clinical PK from the OSP evaluation studies:

| Study | Dose / Route | Key Feature |
|-------|-------------|-------------|
| Hohmann 2015 | 3 mg PO | Microdose CYP3A activity |
| Link 2008 | 7.5 mg PO | Plasma PK profiles |
| Gorski 2003 | PO | Age/sex effects on CYP3A |
| van Dyk 2018 | PO | Inter-individual variability |

Expected midazolam oral PK: Cmax ~50-150 ng/mL, Tmax ~0.5-1 h,
terminal t1/2 ~2-4 h, substantial CYP3A4-mediated first-pass effect (F < 0.5).

### 2.3 Why Midazolam?

Midazolam is the standard CYP3A4 probe substrate because:
- CYP3A4 mediates >90% of its oxidative metabolism (fm ~0.93)
- Both hepatic and intestinal CYP3A4 contribute to first-pass extraction
- Oral bioavailability (~40-50%) is directly modulated by CYP3A4 activity
- UGT1A4-mediated glucuronidation is a minor secondary pathway

This makes midazolam PK highly sensitive to CYP3A4 modulation and
an ideal candidate for DDI risk assessment (see Section 8).

### 2.4 Demographics

| Covariate | Distribution | Range |
|-----------|-------------|-------|
| Body weight | Lognormal (sex-stratified) | 40-150 kg |
| Sex | 55% male / 45% female | Binary |
| Age | Uniform | 22-72 years |

## 3. Model

### 3.1 Structural Model

A 2-compartment model with first-order absorption was selected based on
AIC comparison with a 1-compartment alternative (dAIC = 467). This is
consistent with the known biphasic disposition of midazolam: rapid
distribution into well-perfused central tissues followed by slower
equilibration with peripheral compartments.

Parameters (all apparent, incorporating oral bioavailability F):
- **Ka**: first-order absorption rate constant (1/h)
- **CL/F**: apparent oral clearance (L/h) — CL ~25 L/h divided by F ~0.5
- **Vc/F**: apparent central volume of distribution (L)
- **Q/F**: apparent intercompartmental clearance (L/h)
- **Vp/F**: apparent peripheral volume of distribution (L)

**Note on CL/F:** Apparent clearance from oral data cannot separate
systemic clearance from bioavailability. CL/F ~50 L/h is consistent
with CL ~25 L/h and F ~0.5, but other combinations are possible.
IV data would be needed to resolve this.

### 3.2 Statistical Model

**Between-subject variability (BSV):** Lognormal on Ka, CL/F, Vc/F:

    P_i = theta * exp(eta_i),  eta_i ~ N(0, omega^2)

**Residual error:** Combined proportional + additive:

    DV = IPRED * (1 + eps_prop) + eps_add

The combined model accommodates concentration-dependent assay error at
high concentrations (proportional component) and fixed assay noise near
the LLOQ (additive component, ~0.5 ng/mL).

### 3.3 Covariate Model

Allometric scaling with fixed exponents (Anderson & Holford 2008):
- CL/F, Q/F: (WT/70 kg)^0.75
- Vc/F, Vp/F: (WT/70 kg)^1.0

This is a biologically motivated parameterization. The 0.75 exponent on
clearance reflects the known scaling of metabolic rate with body mass
(CYP3A4 expression correlates with liver size). The 1.0 exponent on
volume reflects direct proportionality between distribution volume and
body mass.

## 4. Estimation

### 4.1 Method

SAEM (Stochastic Approximation Expectation-Maximization) via nlmixr2:
- Burn-in: 200 iterations
- Expectation-maximization: 100 iterations (300 total)
- Convergence assessed by parameter trace stability in final 50 iterations
- -2LL computed by Gaussian quadrature (nnodes=3, nsd=1.6)

### 4.2 Model Comparison

| Model | OFV | AIC | dAIC | Decision |
|-------|-----|-----|------|----------|
| 1-compartment | 3497.8 | 4800.4 | +467 | Rejected |
| 2-compartment | 3027.2 | 4333.7 | ref | Selected |
| + Allometric WT | 2918.4 | 4224.9 | -109 | Final model |

### 4.3 Parameter Estimates (Base 2-Compartment Model)

| Parameter | Estimate | RSE (%) | Unit | Biological Interpretation |
|-----------|----------|---------|------|---------------------------|
| Ka | 2.26 (exp(0.815)) | 14.5% | 1/h | Absorption rate |
| CL/F | 52.3 (exp(3.948)) | 0.8% | L/h | Apparent clearance (70 kg) |
| Vc/F | 46.6 (exp(3.842)) | 2.2% | L | Central volume (70 kg) |
| Q/F | 16.8 (exp(2.821)) | 1.9% | L/h | Intercompartmental CL |
| Vp/F | 62.5 (exp(4.136)) | 1.9% | L | Peripheral volume |
| sigma prop | 0.249 | — | SD | Proportional error (~25%) |
| sigma add | 0.494 | — | ng/mL | Additive error |

Parameter estimates are close to the simulation truth (CL/F true = 50,
estimated = 52.3; Vc/F true = 45, estimated = 46.6), confirming
adequate recovery by SAEM.

### 4.4 Shrinkage Summary

| Random Effect | Omega (var) | Eta SD | Shrinkage | Interpretation |
|---------------|-------------|--------|-----------|----------------|
| eta(CL) | 0.090 | 0.278 | 46.3% | High — limited individual CL information from sparse design |
| eta(Ka) | 0.267 | 0.418 | -39.7% | Negative (inflation) — Ka poorly identified individually |
| eta(Vc) | 0.028 | 0.111 | 33.6% | Borderline — typical for volume in sparse data |

**Interpretation thresholds:** Shrinkage < 20% indicates good individual
estimation. Shrinkage 20-40% is common and acceptable for population-level
inference. Shrinkage > 40% indicates individual parameter estimates are
pulled toward the population mean and should not be used for individual
predictions. Negative shrinkage (eta inflation) suggests the random effect
may be absorbing model misspecification.

The observed shrinkage pattern is consistent with a study design where
80/120 subjects have only 4 observations — insufficient to precisely
estimate individual Ka or CL. Population-level conclusions (typical
values, BSV magnitude, covariate effects) remain valid.

## 5. Diagnostics

### 5.1 Goodness-of-Fit

**Figure: `gof_4panel.png`**

- DV vs PRED: Adequate agreement with no systematic bias
- DV vs IPRED: Tighter correlation confirms BSV captures individual differences
- IRES vs TIME: Random scatter around zero, no time-dependent trends
- IRES vs PRED: No concentration-dependent bias

### 5.2 Random Effects

**Figure: `eta_distributions.png`**

Eta distributions for CL and Ka are approximately normal and centered
near zero, supporting the lognormal BSV assumption. The Vc eta
distribution is narrower, consistent with lower BSV on volume.

### 5.3 Visual Predictive Check

**Figure: `vpc.png`**

The VPC overlays observed percentiles (5th, 50th, 95th) on the 90%
prediction interval from 500 model simulations. In the 0-8 h window
(where both rich and sparse sampling arms contribute data), observed
percentiles fall within simulated prediction intervals, confirming
adequate model performance.

**Late-time instability note:** At t > 8 h the lower (5th percentile)
prediction interval widens sharply on the log scale. This is expected:
at 12 h only the rich-arm subjects (40/120) contribute observations,
and simulated concentrations approach the LLOQ (0.5 ng/mL). The log
scale amplifies small absolute differences near zero. This is a display
artifact of sparse terminal-phase data on a log axis, not evidence of
model misspecification. A prediction-corrected VPC (pcVPC) or
stratification by sampling density would mitigate this visual effect
but was not implemented here.

### 5.4 Individual Fits

**Figure: `individual_fits.png`**

Selected individual concentration-time profiles with PRED (population)
and IPRED (individual) overlays. IPRED tracks individual data more
closely than PRED, confirming that BSV improves individual-level
predictions for subjects with sufficient data.

## 6. Covariate Analysis

### 6.1 Rationale

Body weight is the primary covariate for midazolam clearance and volume:
- CYP3A4 expression correlates with liver size, which scales with body mass
- Volume of distribution scales directly with body mass
- Allometric exponents (0.75/1.0) are fixed based on physiological scaling
  theory (Anderson & Holford 2008), not empirically estimated

### 6.2 Results

The allometric WT model improved fit substantially (dOFV = 109, p < 0.001)
with meaningful BSV reduction:

| Parameter | Base Omega | Covariate Omega | Reduction |
|-----------|-----------|-----------------|-----------|
| eta(Ka) | 0.268 | 0.268 | -0.2% |
| eta(CL) | 0.090 | 0.084 | 6.6% |
| eta(Vc) | 0.028 | 0.013 | 52.4% |

The 52% reduction in Vc BSV confirms that body weight explains a
substantial fraction of volume variability.

**Figure: `forest_covariate_auc.png`** — Forest plot showing AUC ratio
by weight stratum and sex, with the 0.80-1.25 clinical relevance band.

### 6.3 Clinical Relevance Assessment

The 0.80-1.25 band is used here as a **clinical relevance heuristic**
for evaluating whether weight-related exposure differences are likely to
be meaningful. This is not a formal bioequivalence assessment.

## 7. Weight-Exposure Analysis (Verified)

Two metrics are reported to separate the weight effect from individual
variability. Full definitions: `outputs/tables/aucr_definition.md`.

### 7.1 Metric 1: Stratum-Level Median Shift (Isolates Weight Effect)

**Definition:** AUCR_shift = median(AUC at weight W) / median(AUC at 70 kg)

| Stratum | Median Shift | Theoretical | 90% Bootstrap CI | Within 0.80-1.25? |
|---------|-------------|-------------|------------------|-------------------|
| 50 kg | **1.33** | 1.29 | [1.29, 1.36] | No (above) |
| 70 kg (ref) | 1.00 | 1.00 | [0.98, 1.02] | Yes |
| 90 kg | 0.86 | 0.83 | [0.84, 0.88] | Yes |
| 100 kg | **0.78** | 0.77 | [0.76, 0.80] | No (below) |
| 110 kg | **0.73** | 0.71 | [0.71, 0.75] | No (below) |

**Source:** `outputs/tables/aucr_stratum_shift_summary.csv` (script 06,
N=2000 bootstrap replicates)

The narrow bootstrap CIs confirm these shifts are precisely estimated.
They match the theoretical allometric prediction AUCR = (70/WT)^0.75,
confirming internal consistency.

### 7.2 Metric 2: Individual-Level Spread (Weight + BSV Combined)

**Definition:** AUCR_i = AUC_i / median(AUC at 70 kg) for each individual

| Stratum | Median AUCR | 90% PI (individual) | % Outside 0.80-1.25 |
|---------|-------------|---------------------|---------------------|
| 50 kg | 1.33 | [0.84, 2.20] | 62.8% |
| 70 kg (ref) | 1.00 | [0.61, 1.61] | 45.0% |
| 90 kg | 0.86 | [0.53, 1.39] | 52.5% |
| 100 kg | 0.78 | [0.48, 1.26] | 58.7% |
| 110 kg | 0.73 | [0.43, 1.19] | 64.5% |

**Source:** `outputs/tables/aucr_weight_summary_verified.csv`

**Why 45% of the 70 kg reference group falls outside 0.80-1.25:**
This is not a weight effect. Even at the reference weight, individual
AUC varies because of BSV in CL (~30% CV). The lognormal distribution
of CL produces a spread in individual AUCs where roughly half of
subjects fall outside the 0.80-1.25 band around the population median.

### 7.3 Interpretation

The weight effect on typical exposure (Metric 1) is modest but real:
allometric 0.75 scaling on CL shifts median AUC by +33% at 50 kg and
-22% at 100 kg relative to the 70 kg reference.

Individual pharmacokinetic variability (Metric 2) substantially
exceeds the weight-related shift at all strata. The 90% PI for an
individual subject spans 2-3 fold regardless of weight, driven
primarily by BSV in CL.

### 7.3 Dose Evaluation Statement

For the 70-90 kg range, median exposure differences from the 70 kg
reference are within the 0.80-1.25 heuristic band. For patients at
50 kg or 100+ kg, the allometric model predicts median AUC shifts of
~30% and ~22% respectively. Whether these shifts warrant dose
modification depends on the therapeutic index of the specific clinical
context. For midazolam's sedation indication, where titration to effect
is standard practice, fixed 7.5 mg dosing is generally acceptable
across this weight range.

**Figure: `wt_auc_ratio.png`**

### 7.4 Dose-Proportionality

**Figures: `exposure_dose_auc.png`, `exposure_dose_cmax.png`**

AUC and Cmax scale proportionally with dose across 5, 7.5, and 15 mg
scenarios, confirming dose-linear PK within this range. Within-dose
variability (driven by BSV in CL) is consistent across weight strata.

## 8. CYP3A4 Inhibition Sensitivity (Portfolio Bridge)

### 8.1 Scenario

As a conceptual bridge to DDI risk assessment (P10 Static DDI Framework),
CL/F was reduced by 50% to approximate strong CYP3A4 inhibition. This
reflects midazolam's high fm(CYP3A4) ~0.93 under near-complete
enzyme inhibition.

### 8.2 Results

| Metric | Value |
|--------|-------|
| Median AUC fold-change | 1.99x |
| 90% PI | [1.92, 2.03] |
| % subjects with FC > 2x | 29.2% |

**Figure: `cyp3a4_inhibition_fold_change.png`**

The tight fold-change distribution (narrow 90% PI) indicates that
the AUC ratio under inhibition is primarily determined by the degree
of CL reduction, with BSV contributing minimal additional variability
to the fold-change itself.

### 8.3 Connection to P10

This sensitivity analysis demonstrates the logic chain:
1. **Baseline PopPK** (this project) quantifies CL/F and its variability
2. **CYP3A4 inhibition** reduces CL/F proportionally to fm and inhibitor potency
3. **Static DDI assessment** (P10) formalizes this using the mechanistic
   static model (R1, R1gut equations from FDA 2020 DDI Guidance)

The ~2x fold-change from 50% CL reduction is consistent with the AUCR
predictions from P10's mechanistic static model for moderate-to-strong
CYP3A4 inhibitors.

## 9. Limitations

1. **Simulated data:** This analysis uses simulated data calibrated to
   published evidence, not real clinical data. See
   `report/LIMITATIONS_AND_SCOPE.md` for full scope statement.

2. **Single dose only:** Multiple-dose PK, steady-state, and time-dependent
   CYP3A4 effects are not evaluated.

3. **Healthy adults only:** Hepatic/renal impairment, pediatric, and
   geriatric populations are not simulated.

4. **Shrinkage:** High eta shrinkage on CL (46%) reflects the sparse
   study design. Population-level inferences are valid; individual
   empirical Bayes estimates should be interpreted with caution.

5. **Fixed allometric exponents:** The 0.75/1.0 exponents are biologically
   justified but not empirically estimated from the data. In a clinical
   analysis, both fixed and estimated approaches should be compared.

6. **No bootstrap or external validation:** Parameter confidence intervals
   are from asymptotic covariance only. A clinical-grade analysis would
   include 500-1000 bootstrap replicates.

## 10. External Literature Qualification

External structural qualification assessment: evaluates consistency of
model-predicted central tendency (typical prediction; ETAs = 0) against
digitized study-level mean +/- SD profiles. This assessment does not
evaluate individual-level predictive performance.

Because the OSP source data are aggregated (study-level mean +/- SD)
without individual concentrations or covariates, this module supports
structural/typical-profile benchmarking only (not IIV/covariate
validation).

### 10.1 Data Source

12 control/baseline arms from published midazolam PO studies (N = 10-65,
doses 1-15 mg, 10-21 timepoints per profile), extracted from the OSP
Database for Observed Data (Open-Systems-Pharmacology, 2024). All data
are aggregated (study-level means +/- SD); no individual-level
concentrations or subject-level covariates are available. See
`data/literature_osp_midazolam/qualification_set.md` for inclusion
criteria and study details.

### 10.2 Method

Typical profiles were generated at each study's actual dose using fitted
THETAs (Ka = 2.26 h^-1, CL/F = 52.3 L/h, Vc/F = 46.6 L, Q/F = 16.8
L/h, Vp/F = 62.5 L). Population prediction intervals (90%) were
computed from 500 simulated subjects per dose using the simulated PopPK
model's IIV (TRUE_OMEGA), not literature-derived variability.

### 10.3 Metrics

Primary metrics: Unweighted RMSE and AUC ratio (predicted / observed).
The 0.5-2.0 AUC ratio band is used as a conventional benchmarking
reference (guidance-aligned heuristic), not a formal acceptance criterion.

Secondary: Weighted RMSE (1/SD^2 weights) as sensitivity analysis.
Digitized SD values may reflect digitization uncertainty in addition to
biological variability.

### 10.4 Results

**Figures:**
- `lit_overlay_mean_profiles.png` — All studies overlaid with typical curves
- `lit_dose_stratified_overlays.png` — Dose-faceted with 90% PI bands
- `lit_external_vpc_like.png` — External PI overlay (not pcVPC)
- `lit_qualification_summary.png` — Per-study AUC ratio bar chart

**Table:** `outputs/tables/literature_qualification_summary.csv`

### 10.5 Interpretation

Agreement is strongest around 7.5 mg (the probe dose used in the
simulated PopPK dataset), where Mueller 2009 achieves an AUC ratio of
1.10 and Ahonen 1995 reaches 1.28. The 15 mg studies (Olkkola 1993,
Backman 1996, Bornemann 1986) show AUC ratios of 1.24-1.50, consistent
with the model's dose-linear structure.

The main deviation occurs at Hohmann 2015 (3 mg, AUC ratio 2.39), the
only study outside the 2-fold band. This likely reflects a combination
of the low dose (where digitization transfer error has proportionally
larger impact on the absolute concentration values), the study's mixed
microdose/standard design (0.003 mg and 3 mg arms in the same
protocol), and potential differences in assay sensitivity or
formulation between this study and the studies used to calibrate the
OSP PBPK model parameters.

The 1 mg dose group (Wiesinger, Zahner, Chattopadhyay) shows systematic
overprediction (AUC ratios 1.42-1.96). At these low doses, absolute
concentrations are small (Cmax ~5-10 ng/mL), amplifying the relative
impact of digitization transfer error and assay sensitivity differences.

Overall, 11/12 studies fall within the 2-fold benchmarking band, and the
model captures the dose-concentration relationship across a 15-fold dose
range (1-15 mg). This supports structural plausibility of the fitted
PopPK model parameters while preserving the simulated dataset as the
primary basis for the full PopPK workflow demonstration (IIV, covariates,
simulation-based decision analysis).

### 10.6 Limitations

OSP profiles are digitized from published figures; SD values reflect
digitization uncertainty in addition to biological variability. No
covariate matching is possible (no individual WT, SEX, AGE). This
qualification assesses structural model consistency only — it does
not constitute external validation of population-level variability
or covariate effects. PI bands shown in overlay figures are derived
from the simulated PopPK model's IIV parameters, not from
literature-derived variability estimates.

## 11. Conclusions

A 2-compartment PopPK model with first-order absorption adequately describes
the pharmacokinetics of oral midazolam 7.5 mg in a simulated healthy adult
population. Key findings:

- Parameter estimates recover simulation truth within 5-10%, confirming
  SAEM estimation adequacy
- Allometric WT scaling improves fit (dOFV = 109) with 52% reduction in Vc BSV
- Weight-related exposure shifts are modest (median AUCR 0.73-1.33 across
  50-110 kg) but individual variability from BSV dominates the weight effect
- A 50% CL reduction (CYP3A4 inhibition scenario) produces ~2x AUC
  fold-change, consistent with P10 static DDI predictions

This analysis demonstrates an industry-aligned PopPK workflow: data
preparation, model development, diagnostics, covariate evaluation,
and simulation-based exposure assessment — structured to reflect
the reporting conventions of the FDA PopPK Guidance (2022).

---
*P4 PopPK Analysis | OSP-Calibrated Midazolam*
*Submission-style structure | Methodological demonstration*
