# Methods — PopPK Analysis of Oral Midazolam

## 1. Data Generation

A publication-calibrated simulated trial dataset was generated to
demonstrate an industry-aligned PopPK workflow. True PK parameters were
anchored to the Open Systems Pharmacology (OSP) Midazolam PBPK model
(v2.0), validated against 40+ clinical studies spanning IV and oral
administration across a 4-log dose range (0.001-40 mg).

### Model Hierarchy: PBPK to PopPK

The PBPK model provides the mechanistic foundation:
- CYP3A4-mediated hepatic and intestinal metabolism (kcat = 8.76/min, Km = 4.0 uM)
- UGT1A4 secondary glucuronidation pathway
- Rodgers-Rowland tissue partition coefficients (logP = 2.90, fu = 0.031)
- Bioavailability < 50% due to combined first-pass extraction

The PopPK model approximates this mechanistic system as a 2-compartment
disposition model with first-order absorption. Apparent parameters (CL/F, V/F)
implicitly incorporate bioavailability. This approximation is standard practice
for oral drug PopPK analysis when the primary goal is describing population
variability and informing dosing decisions.

### Study Design

| Feature | Specification |
|---------|--------------|
| Drug | Midazolam 7.5 mg oral tablet |
| Population | Healthy adults, N = 120 |
| Rich sampling arm (n=40) | 0.25, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12 h |
| Sparse sampling arm (n=80) | 0.5, 2, 4, 8 h |
| Covariates | WT (continuous), SEX (binary), AGE (continuous) |
| LLOQ | 0.5 ng/mL |
| BLQ handling | M1 method (excluded with MDV=1) |

## 2. Structural Model

Two candidate models were evaluated:

**Model 1: 1-compartment with first-order absorption**

    dA_depot/dt   = -Ka * A_depot
    dA_central/dt =  Ka * A_depot - (CL/F)/V * A_central
    Cp = A_central / V

**Model 2: 2-compartment with first-order absorption** (selected)

    dA_depot/dt      = -Ka * A_depot
    dA_central/dt    =  Ka * A_depot
                       - (CL/F)/(Vc/F) * A_central
                       - (Q/F)/(Vc/F) * A_central
                       + (Q/F)/(Vp/F) * A_peripheral
    dA_peripheral/dt =  (Q/F)/(Vc/F) * A_central
                       - (Q/F)/(Vp/F) * A_peripheral
    Cp = A_central / (Vc/F)

Where:
- A_depot: amount in absorption compartment (ug)
- A_central: amount in central compartment (ug)
- A_peripheral: amount in peripheral compartment (ug)
- Cp: plasma concentration (ug/L = ng/mL)
- Ka: first-order absorption rate constant (1/h)
- CL/F: apparent oral clearance (L/h)
- Vc/F: apparent central volume (L)
- Q/F: apparent intercompartmental clearance (L/h)
- Vp/F: apparent peripheral volume (L)

Initial conditions: A_depot(0) = Dose, A_central(0) = 0, A_peripheral(0) = 0.

Model selection was based on AIC comparison (dAIC = 467 favoring 2-comp).
The 2-compartment model captures the biphasic concentration-time profile
of midazolam: rapid initial distribution (central ↔ peripheral exchange)
followed by a terminal elimination phase dominated by CYP3A4 metabolism.

## 3. Statistical Model

### Between-Subject Variability (BSV)

Log-normal distribution on Ka, CL/F, and Vc/F:

    P_i = theta * exp(eta_i),    eta_i ~ N(0, omega^2)

Where theta is the typical (population) value and eta_i is the individual
deviation. CV% ≈ sqrt(omega^2) * 100 for omega^2 < 0.5.

| Parameter | omega^2 (variance) | Approximate CV |
|-----------|--------------------|----------------|
| Ka | 0.267 | ~52% |
| CL/F | 0.090 | ~30% |
| Vc/F | 0.028 | ~17% |

### Residual Error Model

Combined proportional + additive:

    DV_ij = IPRED_ij * (1 + eps_prop) + eps_add

Where:
- eps_prop ~ N(0, sigma_prop^2), sigma_prop ≈ 0.25 (25% proportional)
- eps_add ~ N(0, sigma_add^2), sigma_add ≈ 0.5 ng/mL

**Justification:** The proportional component captures concentration-dependent
measurement variability (e.g., chromatographic peak integration error scales
with signal magnitude). The additive component captures fixed-noise sources
(e.g., extraction recovery variability, instrument baseline noise) that
dominate near the LLOQ. A purely proportional model would underpredict
error at low concentrations; a purely additive model would overpredict
error at high concentrations.

## 4. Estimation

### 4.1 Algorithm

SAEM (Stochastic Approximation Expectation-Maximization) via nlmixr2 v5.0.0.

| Setting | Value |
|---------|-------|
| Burn-in iterations | 200 |
| EM iterations | 100 |
| Total iterations | 300 |
| -2LL method | Gaussian quadrature (nnodes=3, nsd=1.6) |

### 4.2 Convergence Assessment

Parameter traces were inspected for stability in the final 50 EM iterations.
All fixed-effect parameters (log-transformed) showed stable oscillation
around their final values, with no drift or divergence.

The 2-compartment base model converged in ~2 minutes. The allometric
covariate model converged in ~3.5 minutes.

### 4.3 Parameter Precision

Relative standard errors (RSE) were computed from the asymptotic
covariance matrix for fixed effects:

| Parameter | Estimate | SE | RSE (%) |
|-----------|----------|-----|---------|
| log(Ka) | 0.815 | 0.118 | 14.5% |
| log(CL/F) | 3.948 | 0.031 | 0.8% |
| log(Vc/F) | 3.842 | 0.086 | 2.2% |
| log(Q/F) | 2.821 | 0.055 | 1.9% |
| log(Vp/F) | 4.136 | 0.080 | 1.9% |

RSE < 30% for all fixed effects, indicating adequate parameter precision.
Ka has the highest RSE (14.5%), reflecting the inherent difficulty of
estimating absorption rate from sparse sampling designs.

Note: Asymptotic SEs for variance components (omega, sigma) were not
available from the SAEM covariance matrix. A bootstrap analysis would
be needed for confidence intervals on BSV parameters.

### 4.4 Eta Shrinkage

| Random Effect | Shrinkage | Threshold | Assessment |
|---------------|-----------|-----------|------------|
| eta(CL) | 46.3% | >40% = high | Individual CL estimates pulled toward population mean |
| eta(Ka) | -39.7% | <0% = inflation | Ka poorly identified individually; eta may absorb misspecification |
| eta(Vc) | 33.6% | 20-40% = moderate | Acceptable for population inference |

**Epsilon shrinkage** was not separately computed but is expected to
be moderate given the combined residual error model.

**Impact:** High CL shrinkage means individual empirical Bayes estimates
(IPRED, individual CL) are biased toward the population mean. This does
not invalidate population-level parameter estimates, BSV magnitudes, or
covariate effects — only individual-level predictions should be interpreted
cautiously.

## 5. Covariate Analysis

### Allometric Scaling (Biologically Motivated)

    CL/F = theta_CL * (WT/70)^0.75
    Vc/F = theta_Vc * (WT/70)^1.0
    Q/F  = theta_Q  * (WT/70)^0.75
    Vp/F = theta_Vp * (WT/70)^1.0

Allometric exponents (0.75 for clearance, 1.0 for volume) are fixed based on
established physiological scaling theory (Anderson and Holford 2008). This is
not an empirical fit but a mechanistically justified parameterization:
- CYP3A4 hepatic metabolism scales with liver size, which scales allometrically
  with body weight (exponent ~0.75)
- Distribution volume scales approximately linearly with body mass (exponent ~1.0)

### Evaluation Criteria

1. **Statistical:** dOFV = 108.8 (p < 0.001 for chi-sq, 0 additional df
   since exponents are fixed)
2. **BSV reduction:** 52% reduction in omega(Vc), 6.6% in omega(CL)
3. **Clinical relevance:** Median AUCR evaluated against 0.80-1.25 heuristic band

## 6. Decision Simulations

1000 virtual subjects per scenario, simulated from final model parameters
with BSV sampled from estimated omega matrix.

### Dose-proportionality scenarios
- 5 mg, 7.5 mg, 15 mg across weight strata (50, 70, 100 kg)
- Exposure metrics: AUC(0-inf) by trapezoidal rule, Cmax

### Weight-exposure verification
- 5 strata: 50, 70, 90, 100, 110 kg at 7.5 mg
- AUCR computed relative to 70 kg median
- 90% prediction intervals and % outside 0.80-1.25 band

### CYP3A4 inhibition sensitivity
- CL/F reduced by 50% (strong inhibitor approximation)
- AUC fold-change distribution (inhibited vs baseline)
- Conceptual bridge to P10 static DDI framework

## 7. Software

| Tool | Version | Purpose |
|------|---------|---------|
| R | 4.5.0 | Analysis environment |
| nlmixr2 | 5.0.0 | PopPK estimation (SAEM) |
| rxode2 | 5.0.1 | ODE solver for nlmixr2 |
| deSolve | 1.41 | Simulation ODE solver |
| tidyverse | 2.0.0 | Data manipulation and visualization |
| ggplot2 | 3.5.2 | Figures |
| patchwork | 1.3.0 | Figure composition |

Platform: x86_64-apple-darwin20, macOS.
