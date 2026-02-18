# Methods Appendix — Mathematical Derivations and Assumptions

## A1. Concentration Calculations

### A1.1 Unbound Systemic Concentration

The unbound maximum plasma concentration is the primary metric for systemic DDI screening:

```
Cmax_uM = (Cmax_total [ug/mL] / MW [g/mol]) * 1000
[I]max,u = Cmax_uM * fu
```

Where:
- Cmax_total: Total (bound + unbound) peak plasma concentration
- MW: Molecular weight of the perpetrator
- fu: Fraction unbound in plasma

**Case A (Ketoconazole):** [I]max,u = (8.07/531.43)*1000 * 0.01 = 0.152 uM
**Case B (Sotorasib):** [I]max,u = (7.87/560.6)*1000 * 0.11 = 1.544 uM

### A1.2 Intestinal Concentration ([I]gut)

The FDA assumes complete dissolution of the oral dose in 250 mL of gastrointestinal fluid:

```
dose_ug = dose_mg * 1000
[I]gut_ugml = dose_ug / 250 [mL]
[I]gut_uM = ([I]gut_ugml / MW) * 1000
```

Equivalently (verified via dual-method unit check):

```
dose_umol = (dose_mg / MW) * 1000
[I]gut_uM = dose_umol / 0.250 [L]
```

Both methods are verified to agree within floating-point precision (<1e-6 uM) in `outputs/logs/igut_sanity.txt`.

**Rationale:** This is deliberately conservative. In vivo, incomplete dissolution, GI transit, and absorption mean the actual luminal concentration is lower. The purpose is high-sensitivity screening.

## A2. CYP Inhibition Screening

### A2.1 Reversible Inhibition (R1)

The basic screening ratio for competitive/reversible inhibition:

```
R1 = 1 + [I]max,u / Ki
```

- **Threshold:** R1 >= 1.02 triggers further evaluation
- **Interpretation:** R1 represents the fold-increase in the victim's Km,apparent when the inhibitor occupies the enzyme active site
- **Conservatism:** Uses unbound Cmax (worst-case systemic exposure)

### A2.2 Intestinal R1 (R1,gut)

For gut wall metabolism:

```
R1,gut = 1 + [I]gut / Ki
```

- **Threshold:** R1,gut >= 11 triggers further evaluation
- **Note:** The higher threshold (11 vs 1.02) partially compensates for the extreme conservatism of [I]gut

### A2.3 Relationship Between R1 and AUCR

For a single-enzyme inhibition without TDI or gut contribution:

```
AUCR = 1 / (fm/R1 + (1-fm))
```

This simplifies to: AUCR approaches 1/((1-fm)) when R1 >> 1 (complete inhibition), capped by the fraction metabolized by alternative pathways.

## A3. Time-Dependent Inhibition (TDI)

### A3.1 Mechanism

TDI (mechanism-based inhibition) involves covalent modification of the CYP enzyme, leading to irreversible loss of catalytic activity. The enzyme can only be replaced through de novo synthesis.

### A3.2 Inactivation Rate

The apparent first-order inactivation rate constant:

```
lambda = kinact * [I]max,u / (KI + [I]max,u)
```

Where:
- kinact: Maximum inactivation rate constant (min-1)
- KI: Inhibitor concentration producing half-maximal inactivation (uM)
- [I]max,u: Unbound systemic concentration (uM)

### A3.3 Steady-State Enzyme Fraction

At steady state, the fraction of functional enzyme remaining is:

```
TDI_factor = kdeg / (kdeg + lambda)
```

Where kdeg is the natural degradation rate constant of the CYP enzyme.

**Interpretation:**
- TDI_factor = 1.0: No TDI (enzyme fully functional)
- TDI_factor = 0.5: 50% of enzyme inactivated
- TDI_factor = 0.05: 95% of enzyme inactivated (severe)

**Case A:** TDI_factor = 0.00032 / (0.00032 + 0.048*0.152/(0.86+0.152)) = 0.062 (94% loss)
**Case B:** TDI_factor = 0.00032 / (0.00032 + 0.04*1.544/(3.5+1.544)) = 0.027 (97% loss)

### A3.4 Enzyme Degradation Rate (kdeg)

```
kdeg = ln(2) / t1/2_enzyme
```

For hepatic CYP3A4: t1/2 ~36 hours -> kdeg = 0.693/36/60 = 0.00032 min-1

This is a critical parameter because it determines the rate of enzyme recovery. Longer enzyme half-life means slower recovery from TDI and greater DDI magnitude.

## A4. CYP Induction

### A4.1 R3 Equation

The induction screening ratio incorporates the maximal induction effect and the inhibitor concentration relative to the EC50:

```
R3 = 1 / (1 + d * Emax * [I]max,u / (EC50 + [I]max,u))
```

Where:
- Emax: Maximum fold-increase in mRNA/enzyme activity (from in vitro)
- EC50: Concentration producing half-maximal induction (uM)
- d: Scaling factor for in vitro-to-in vivo translation (default = 1)

**Threshold:** R3 <= 0.8 indicates clinically relevant induction potential.

### A4.2 The d-Scaling Factor

The d parameter is the most uncertain element of the induction assessment:

- d = 1.0: Full in vitro Emax translates to in vivo (conservative)
- d < 1.0: Empirical calibration suggests in vivo effect is attenuated
- d is compound-specific and ideally calibrated against clinical data

For sotorasib with d=1.0: R3 = 1/(1 + 1*8.1*1.544/(1.5+1.544)) = 0.226

### A4.3 Why Induction Cannot Simply Be Combined with TDI in Static Models

Induction (increased enzyme synthesis) and TDI (irreversible enzyme destruction) act on the same enzyme pool but through fundamentally different mechanisms:

- **TDI** is instantaneous (relative to enzyme turnover) and proportional to current drug concentration
- **Induction** requires transcription, translation, and protein maturation — onset takes 2-5 days
- **At steady state**, TDI maximally destroys enzyme while induction maximally produces new enzyme

The static model applies both simultaneously at [I]max,u, but in reality:
- Early in treatment: TDI may dominate (fast onset, induction not yet established)
- At steady state: Induction may dominate if transcriptional upregulation exceeds TDI-mediated loss
- The balance is time-dependent and concentration-profile-dependent

This is why PBPK modeling (which tracks enzyme pool dynamics over time) is essential for compounds with dual TDI + induction.

## A5. AUCR Prediction

### A5.1 Hepatic Component

```
CL_int_ratio = (1/R1) * TDI_factor
AUCR_hepatic = 1 / (fm * CL_int_ratio + (1-fm))
```

The fm term captures the fraction of total clearance mediated by the affected enzyme. When fm is high (e.g., 0.94 for midazolam/CYP3A4), even modest inhibition produces large AUCR changes. When fm is low, alternative pathways compensate.

### A5.2 Gut Wall Component

```
CL_gut_ratio = 1 / R1,gut
AUCR_gut = 1 / (fg * CL_gut_ratio + (1-fg))
```

Where fg = fraction of oral dose metabolized in the gut wall (= 1 - Fg, where Fg is gut availability).

### A5.3 Total AUCR

```
AUCR_total = AUCR_hepatic * AUCR_gut
```

This multiplicative assumption treats hepatic and gut contributions as independent and sequential, which is an approximation. In reality, changes in gut wall metabolism affect the amount of drug reaching the liver.

### A5.4 Worked Example: Case A

```
R1 = 1 + 0.152/0.015 = 11.12
R1,gut = 1 + 3010.7/0.015 = 200,716
TDI_factor = 0.00032 / (0.00032 + 0.048*0.152/(0.86+0.152)) = 0.0618
CL_int_ratio = (1/11.12) * 0.0618 = 0.00556
AUCR_h = 1 / (0.94*0.00556 + 0.06) = 1 / 0.0652 = 15.72
CL_gut_ratio = 1/200716 = 4.98e-6
AUCR_g = 1 / (0.43*4.98e-6 + 0.57) = 1/0.570 = 1.754
AUCR_total = 15.72 * 1.754 = 27.59
```

## A6. DDI Classification

Per FDA 2020 guidance:

| AUCR Range | Classification |
|------------|---------------|
| >= 5.0 | Strong inhibitor |
| 2.0 to <5.0 | Moderate inhibitor |
| 1.25 to <2.0 | Weak inhibitor |
| < 1.25 | No clinically relevant interaction |

For induction (AUCR < 1):

| AUCR Range | Classification |
|------------|---------------|
| < 0.5 | Strong inducer |
| 0.5 to <0.8 | Moderate inducer |
| 0.8 to <1.0 | Weak inducer |

## A7. Transporter Screening

### A7.1 Concentration Selection

| Transporter Location | Relevant Concentration | Rationale |
|---------------------|----------------------|-----------|
| Intestinal efflux (P-gp, BCRP) | [I]gut | Pre-absorption luminal concentration |
| Hepatic uptake (OATP1B1/1B3) | [I]max,u | Portal/systemic unbound concentration |
| Renal (OAT1/3, OCT2, MATE1/2-K) | [I]max,u | Systemic unbound concentration |

### A7.2 Decision Criteria

```
Ratio = I_relevant / IC50
```

| Transporter Type | Threshold | Action if Exceeded |
|-----------------|-----------|-------------------|
| P-gp, BCRP | Ratio >= 10 | Clinical DDI study recommended |
| OATP1B1/1B3 | Ratio >= 0.1 | Clinical DDI study recommended |
| OAT, OCT, MATE | Ratio >= 0.1 | Clinical DDI study recommended |

## A8. Sensitivity Analysis Methods

### A8.1 Tornado (One-at-a-Time)

For each parameter p_i:
1. Hold all other parameters at base values
2. Compute AUCR at p_i,low and p_i,high
3. Record delta_low = AUCR(p_i,low) - AUCR_base and delta_high = AUCR(p_i,high) - AUCR_base
4. Rank parameters by max(|delta_low|, |delta_high|)

Parameter ranges used for Case A:

| Parameter | Base | Low | High | Rationale |
|-----------|------|-----|------|-----------|
| Ki (uM) | 0.015 | 0.0015 | 0.15 | 10-fold range around measured value |
| [I]max,u (uM) | 0.152 | 0.076 | 0.228 | +/- 50% |
| fm | 0.94 | 0.50 | 0.95 | Physiological range |
| fg | 0.43 | 0.00 | 0.70 | Zero to high gut extraction |
| kinact (1/min) | 0.048 | 0.024 | 0.072 | +/- 50% |
| KI (uM) | 0.86 | 0.43 | 1.29 | +/- 50% |
| kdeg (1/min) | 0.00032 | 0.00016 | 0.00048 | +/- 50% |

### A8.2 Monte Carlo

Simultaneous variation of all parameters using probability distributions:

| Parameter | Distribution | Parameters | Rationale |
|-----------|-------------|------------|-----------|
| Ki | Lognormal | mu=ln(0.015), sigma=0.3 | Positive-valued, right-skewed |
| [I]max,u | Lognormal | mu=ln(0.152), sigma=0.3 | PK variability |
| fm | Normal (truncated) | mu=0.94, sd=0.05, [0.30, 0.99] | Bounded fraction |
| fg | Normal (truncated) | mu=0.43, sd=0.10, [0.00, 0.80] | Bounded fraction |
| kinact | Lognormal | mu=ln(0.048), sigma=0.3 | Positive-valued |
| KI | Lognormal | mu=ln(0.86), sigma=0.3 | Positive-valued |
| kdeg | Lognormal | mu=ln(0.00032), sigma=0.3 | Positive-valued |

N = 10,000 samples with seed = 42 for reproducibility.

### A8.3 Scenario Analysis (Sotorasib)

For the net effect analysis, the hepatic CL ratio combines TDI and induction:

```
CL_ratio_net = (1/R1) * TDI_factor * induction_fold
```

Where:
```
induction_fold = 1 + d * Emax * [I]max,u / (EC50 + [I]max,u)
```

The induction_fold increases CL (decreases AUC), while (1/R1)*TDI_factor decreases CL (increases AUC). The net direction depends on which effect dominates.

---
*P10 Static DDI Risk Assessment Framework — Methods Appendix*
