# Static DDI Risk Assessment Framework — Comprehensive Report

**FDA 2020 Mechanistic Static Model Implementation**

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [Regulatory Context](#2-regulatory-context)
3. [Methods](#3-methods)
4. [Case A: Midazolam + Ketoconazole (Benchmark)](#4-case-a-midazolam--ketoconazole)
5. [Case B: Sotorasib (Multi-Mechanism Perpetrator)](#5-case-b-sotorasib)
6. [Cross-Case Comparison](#6-cross-case-comparison)
7. [Sensitivity & Uncertainty Analysis](#7-sensitivity--uncertainty-analysis)
8. [Limitations of the Static Approach](#8-limitations-of-the-static-approach)
9. [Conclusions & Regulatory Implications](#9-conclusions--regulatory-implications)
10. [References](#10-references)

---

## 1. Introduction

Drug-drug interactions (DDIs) mediated by cytochrome P450 (CYP) enzymes and drug transporters represent a major safety concern in drug development. The FDA requires sponsors to evaluate DDI potential using a tiered approach: in vitro screening, static mechanistic modeling, and (when warranted) dynamic PBPK modeling or clinical studies.

This project implements the **mechanistic static model** described in the FDA 2020 guidance "In Vitro Drug Interaction Studies — Cytochrome P450 Enzyme- and Transporter-Mediated Drug Interactions." The framework:

- Computes all required screening ratios (R1, R1,gut, R3) and transporter metrics
- Predicts AUCR using the mechanistic static equation incorporating reversible inhibition, time-dependent inhibition (TDI), and gut wall metabolism
- Classifies DDI severity per FDA thresholds
- Provides sensitivity and uncertainty analysis to characterize prediction robustness

Two cases are analyzed to demonstrate different facets of the framework:

| Case | Purpose | Key Challenge |
|------|---------|---------------|
| A: Ketoconazole + Midazolam | Benchmark validation | Verify framework against a well-characterized strong inhibitor |
| B: Sotorasib | Real-world complexity | Resolve opposing TDI and induction on CYP3A4; multi-transporter assessment |

## 2. Regulatory Context

### 2.1 FDA Decision Framework

The FDA 2020 guidance establishes a stepwise decision process:

1. **Basic screening** (R1, R1,gut): Identifies CYP enzymes where inhibition may be clinically relevant
2. **Mechanistic static model**: Integrates reversible inhibition, TDI, and substrate-specific parameters (fm, fg) to predict AUCR
3. **PBPK modeling**: Required when static models flag risk but quantitative prediction is needed, or when opposing mechanisms exist
4. **Clinical DDI study**: Gold standard for confirming or refuting DDI magnitude

### 2.2 Classification Thresholds

| AUCR | Classification | Regulatory Action |
|------|---------------|-------------------|
| >= 5.0 | Strong inhibitor | Contraindication or dose adjustment usually needed |
| 2.0 - 5.0 | Moderate inhibitor | Dose adjustment may be needed |
| 1.25 - 2.0 | Weak inhibitor | Monitor; labeling language |
| < 1.25 | No interaction | No action needed |

For induction: AUCR < 0.5 (strong), 0.5-0.8 (moderate), 0.8-1.0 (weak).

### 2.3 Transporter Thresholds

| Transporter Type | Metric | Threshold |
|-----------------|--------|-----------|
| Intestinal efflux (P-gp, BCRP) | [I]gut / IC50 | >= 10 |
| Hepatic uptake (OATP1B1/3) | [I]max,u / IC50 | >= 0.1 |
| Renal (OAT, OCT, MATE) | [I]max,u / IC50 | >= 0.1 |

## 3. Methods

### 3.1 Concentration Metrics

Two perpetrator concentrations drive the screening:

**Unbound systemic concentration ([I]max,u):**
$$[I]_{max,u} = C_{max,total} \times f_u / MW \times 1000$$

**Intestinal concentration ([I]gut):**
$$[I]_{gut} = \frac{Dose_{mg} \times 1000}{250 \text{ mL} \times MW} \times 1000$$

The 250 mL GI volume is the FDA-mandated conservative assumption. See `outputs/logs/igut_sanity.txt` for dual-method unit verification.

### 3.2 Screening Equations

**Reversible inhibition:**
- R1 = 1 + [I]max,u / Ki (threshold: R1 >= 1.02)
- R1,gut = 1 + [I]gut / Ki (threshold: R1,gut >= 11)

**Time-dependent inhibition:**
- lambda = kinact * [I]max,u / (KI + [I]max,u)
- TDI_factor = kdeg / (kdeg + lambda)

**Induction:**
- R3 = 1 / (1 + d * Emax * [I]max,u / (EC50 + [I]max,u))
- Threshold: R3 <= 0.8

### 3.3 AUCR Prediction

**Hepatic:**
$$AUCR_h = \frac{1}{f_m \times \frac{1}{R1} \times TDI_{factor} + (1 - f_m)}$$

**Gut wall:**
$$AUCR_g = \frac{1}{f_g \times \frac{1}{R1_{gut}} + (1 - f_g)}$$

**Total:** AUCR = AUCR_h x AUCR_g

### 3.4 Sensitivity Analysis Methods

- **Tornado (one-at-a-time):** Each parameter varied independently while others held at base values
- **Heatmap grid:** Ki and [I]max,u co-varied on a 100x100 log-scale grid
- **Monte Carlo:** 10,000 samples with lognormal distributions (CV ~30%) for all key parameters
- **Scenario analysis:** Sotorasib net effect evaluated across induction scaling factors

Full mathematical details in [Methods Appendix](METHODS_APPENDIX.md).

## 4. Case A: Midazolam + Ketoconazole

### 4.1 Rationale

Ketoconazole is the FDA-recommended index strong CYP3A4 inhibitor. Its DDI with midazolam (AUCR = 11.2) is one of the best-characterized interactions in clinical pharmacology. This case serves as a positive control to validate the framework.

### 4.2 Key Results

| Metric | Predicted | Observed | Match |
|--------|-----------|----------|-------|
| AUCR (total) | 27.6 | 11.2 | Classification matches (Strong) |
| AUCR (hepatic) | 15.7 | — | Dominant contributor |
| AUCR (gut) | 1.76 | — | Modest gut contribution |
| Pred/Obs ratio | 2.46 | — | Expected overprediction |

### 4.3 Sensitivity Summary

- **Most sensitive:** fg (gut extraction) and fm (fraction metabolized)
- **Least sensitive:** KI, kdeg (TDI kinetic parameters)
- **Monte Carlo:** Median = 28.0, 90% CI [11.5, 145.4], 100% classified Strong
- **Robustness:** Classification is insensitive to all parameter uncertainty

### 4.4 Interpretation

The 2.46x overprediction is a well-known property of the static model — it uses peak concentrations and assumes steady-state maximal inhibition. PBPK models for this pair predict AUCR = 10-13, much closer to the observed 11.2. The static model's purpose is screening with high sensitivity, not quantitative prediction.

See [Case A Report](CASE_A_midazolam.md) for full details.

## 5. Case B: Sotorasib

### 5.1 Rationale

Sotorasib (LUMAKRAS) represents a modern drug development challenge: a KRAS G12C inhibitor with opposing CYP3A4 mechanisms (TDI + induction), clinically relevant CYP2C8 inhibition, and broad transporter inhibition. This case demonstrates the framework's ability to identify complex risk profiles and recognize its own limitations.

### 5.2 CYP Profile

| Enzyme | Mechanism | R1 | Flag | Clinical Relevance |
|--------|-----------|-----|------|-------------------|
| CYP3A4 | Reversible + TDI | 1.15 | Yes | Net inducer in vivo |
| CYP2C8 | Reversible | 2.93 | Yes | Confirmed in label |
| CYP2C9 | Reversible | 1.03 | Marginal | Not clinically relevant |
| CYP2C19 | Reversible | 1.03 | Marginal | Not clinically relevant |

### 5.3 The TDI vs. Induction Problem

| Mechanism | In Vitro Signal | Static Prediction | Clinical Reality |
|-----------|----------------|-------------------|-----------------|
| TDI (CYP3A4) | TDI factor = 0.027 | Strong inhibition | — |
| Induction (CYP3A4) | R3 = 0.23 | Strong induction | — |
| **Net effect** | — | **Inhibition dominant** | **Induction dominant (AUCR = 0.48)** |

The static model predicts net inhibition across all tested d-values (0.2-1.0), while the clinical result is a 52% AUC decrease (net induction). This discrepancy is expected: the static model cannot account for the different timescales of TDI (rapid, concentration-dependent) and induction (slow, transcriptional).

### 5.4 Transporter Profile

6 of 9 screened transporters are flagged: P-gp, BCRP (gut efflux), OATP1B1/1B3 (hepatic uptake), MATE1/2-K (renal). Clinical studies confirmed P-gp inhibition and modest OATP effects.

See [Case B Report](CASE_B_sotorasib.md) for full details.

## 6. Cross-Case Comparison

| Feature | Case A (Ketoconazole) | Case B (Sotorasib) |
|---------|----------------------|---------------------|
| Primary CYP mechanism | Potent reversible + TDI | Moderate TDI + strong induction |
| Direction | Unambiguous inhibition | Ambiguous (opposing mechanisms) |
| Static model accuracy | Correct classification, 2.5x overprediction | Incorrect direction (predicts inhibition, observed induction) |
| Sensitivity drivers | Victim parameters (fm, fg) | Induction scaling factor (d) |
| PBPK needed? | Optional (for quantification) | **Essential** (to resolve net direction) |
| Transporter risk | None flagged | 6/9 flagged |
| Regulatory complexity | Low (standard index inhibitor) | High (multi-mechanism, multi-target) |

This comparison highlights the central lesson: **the static model is a screening tool, not a universal predictor**. It excels at Case A (single dominant mechanism) and correctly flags all risks in Case B, but cannot resolve the quantitative net effect when opposing mechanisms compete.

## 7. Sensitivity & Uncertainty Analysis

### 7.1 Tornado Results (Case A)

The one-at-a-time analysis reveals that victim parameters (fm, fg) dominate the prediction uncertainty, while inhibitor potency parameters have modest impact. This is because ketoconazole's Ki = 0.015 uM is so far below [I]max,u = 0.152 uM that even 10x changes in Ki do not substantially alter the near-complete inhibition.

### 7.2 Monte Carlo (Case A)

With 10,000 samples and ~30% CV on all parameters:
- The entire distribution classifies as Strong (100%)
- The observed AUCR (11.2) falls near the 5th percentile
- This confirms both the robustness of the classification and the conservative bias

### 7.3 Scenario Analysis (Case B)

The d-scaling sensitivity analysis is the most informative output for sotorasib:
- At d = 1.0 (conservative): AUCR = 6.0 (still predicts inhibition)
- At d = 0.5 (moderate): AUCR = 8.1 (still predicts inhibition)
- At d = 0.2 (empirical): AUCR = 10.2 (still predicts inhibition)
- Clinical reality: AUCR = 0.48 (induction dominant)

No tested d-value produces the observed result, indicating the static framework's fundamental inability to capture the in vivo interplay. The effective d for sotorasib would need to exceed the mathematically possible range in the static model, confirming that dynamic (PBPK) modeling is required.

## 8. Limitations of the Static Approach

1. **Peak concentration assumption:** Uses Cmax rather than time-averaged concentrations, leading to systematic overprediction of inhibition magnitude
2. **No temporal dynamics:** Cannot model the time-dependent onset and offset of induction vs. the rapid kinetics of TDI
3. **No parallel pathway competition:** Assumes independent hepatic and gut contributions without accounting for shifts in metabolic routing
4. **Linear scaling of induction:** The d-parameter is a crude approximation of the in vitro-to-in vivo translation
5. **No auto-inhibition/auto-induction:** Cannot account for the perpetrator's effect on its own metabolism
6. **Single dose level:** Does not capture dose-dependent shifts in the balance of opposing mechanisms

These limitations are acknowledged in the FDA guidance, which positions the static model as a screening tier, not a replacement for PBPK or clinical studies.

## 9. Conclusions & Regulatory Implications

### 9.1 Framework Validation

The static DDI framework is validated by:
- Correct Strong classification for the ketoconazole-midazolam benchmark (AUCR = 27.6 vs. observed 11.2)
- Appropriate identification of all risk signals for sotorasib (TDI, induction, CYP2C8, 6 transporters)
- Transparent sensitivity analysis demonstrating classification robustness (Case A) and model limitations (Case B)

### 9.2 When to Escalate

| Scenario | Action |
|----------|--------|
| Single mechanism, clear flag | Static model sufficient for screening; PBPK optional |
| Opposing mechanisms (TDI + induction) | **PBPK essential** to resolve net direction |
| Quantitative prediction needed | PBPK recommended (static models overpredict) |
| Multiple flagged transporters | Clinical studies for highest-risk substrates |
| Regulatory submission | Static screening + PBPK + clinical confirmation |

### 9.3 Portfolio Significance

This project demonstrates:
- Fluency with FDA 2020 DDI guidance methodology
- Implementation of the complete mechanistic static model in R
- Ability to identify when static models are sufficient vs. when PBPK is needed
- Sensitivity analysis skills (tornado, Monte Carlo, scenario analysis)
- Scientific communication through regulatory-style reports and publication-quality figures

## 10. References

1. FDA (2020). *In Vitro Drug Interaction Studies — Cytochrome P450 Enzyme- and Transporter-Mediated Drug Interactions.* Guidance for Industry. January 2020.
2. Obach RS, Walsky RL, Venkatakrishnan K, et al. (2007). Mechanism-based inactivation of human cytochrome P450 enzymes and the prediction of drug-drug interactions. *Drug Metab Dispos*, 35(2):246-255.
3. Huang SM, Temple R, Throckmorton DC, Lesko LJ (2007). Drug interaction studies: study design, data analysis, and implications for dosing and labeling. *Clin Pharmacol Ther*, 81(2):298-304.
4. Yang J, Liao M, Shou M, et al. (2008). Cytochrome P450 turnover: regulation of synthesis and degradation. *Drug Metab Dispos*, 36(10):2060-2072.
5. Galetin A, Burt H, Gibbons L, Houston JB (2006). Prediction of time-dependent CYP3A4 drug-drug interactions: impact of enzyme degradation, parallel elimination pathways, and intestinal inhibition. *Drug Metab Dispos*, 34(1):166-175.
6. LUMAKRAS (sotorasib) NDA 214665, FDA Clinical Pharmacology Review, 2021.
7. Fahmi OA, Maurer TS, Kish M, et al. (2008). A combined model for predicting CYP3A4 clinical net drug-drug interaction based on CYP3A4 inhibition, inactivation, and induction determined in vitro. *Drug Metab Dispos*, 36(8):1698-1708.

---
*P10 Static DDI Risk Assessment Framework*
*FDA 2020 Mechanistic Static Model*
