# Case B: Sotorasib — Multi-Mechanism DDI Perpetrator Assessment

## 1. Executive Summary

Sotorasib (LUMAKRAS, 960 mg QD) is a KRAS G12C inhibitor with a **complex DDI perpetrator profile**: simultaneous CYP3A4 time-dependent inhibition (TDI) and CYP3A4 induction, strong CYP2C8 reversible inhibition, and inhibition of multiple transporters (P-gp, BCRP, OATP1B1/1B3, MATE1/2-K).

The key regulatory challenge is that the **net CYP3A4 effect** (TDI vs. induction) **cannot be resolved by the static model**. The static framework correctly identifies all risk signals and demonstrates that the net prediction is highly sensitive to the induction scaling factor (d). Clinically, sotorasib is a **net inducer** of CYP3A4 (52% decrease in midazolam AUC), but the static model with default parameters predicts net inhibition — a known limitation that mandates PBPK modeling.

**Regulatory recommendation:** PBPK modeling is essential to resolve the opposing CYP3A4 mechanisms. Clinical DDI studies are warranted for CYP2C8 substrates and flagged transporters.

---

## 2. Perpetrator Characterization

| Parameter | Value | Source |
|-----------|-------|--------|
| Drug | Sotorasib (LUMAKRAS) | KRAS G12C inhibitor |
| Indication | NSCLC with KRAS G12C mutation | FDA approved 2021 |
| Dose | 960 mg QD | LUMAKRAS label |
| MW | 560.6 g/mol | PubChem |
| Cmax,total | 7.87 ug/mL | LUMAKRAS label (geometric mean) |
| fu | 0.11 (11%) | LUMAKRAS label (89% protein bound) |
| **[I]max,u** | **1.544 uM** | Computed |
| **[I]gut** | **6,850 uM** | Computed |

### [I]gut Sanity Check
- Method 1: 960,000 ug / 250 mL = 3,840 ug/mL -> (3,840/560.6)*1000 = **6,850 uM**
- Method 2: 960/560.6*1000 = 1,712 umol -> 1,712/0.250 L = **6,850 uM**
- Unit consistency: **PASS**
- [I]gut / [I]max,u ratio = **4,436x**

## 3. CYP Inhibition Screening

| Enzyme | Ki (uM) | R1 | R1 Flag | Interpretation |
|--------|---------|-----|---------|----------------|
| **CYP3A4** | 10.0 | 1.154 | **FLAGGED** | Weak reversible + TDI |
| **CYP2C8** | 0.8 | 2.930 | **FLAGGED** | Clinically relevant (confirmed in label) |
| CYP2C9 | 50.0 | 1.031 | **FLAGGED** | Marginal |
| CYP2C19 | 50.0 | 1.031 | **FLAGGED** | Marginal |
| CYP2D6 | 50.0 | 1.031 | **FLAGGED** | Marginal |
| CYP1A2 | 50.0 | 1.031 | **FLAGGED** | Marginal |

CYP2C8 shows the strongest reversible inhibition (R1 = 2.93), consistent with the LUMAKRAS label that identifies clinically relevant CYP2C8 inhibition (IC50 ~1.6 uM). CYP3A4 reversible inhibition is modest (R1 = 1.15), but this is compounded by TDI.

### CYP3A4 Time-Dependent Inhibition

| Parameter | Value | Source |
|-----------|-------|--------|
| kinact | 0.04 min-1 | FDA NDA 214665 review |
| KI | 3.5 uM | FDA NDA 214665 review |
| kdeg | 0.00032 min-1 | Yang et al. 2008 |
| **TDI factor** | **0.0268** | Severe inactivation (~97% CYP3A4 loss) |

## 4. CYP3A4 Induction Screening

| Parameter | Value | Source |
|-----------|-------|--------|
| Emax | 8.1-fold | FDA NDA review (hepatocyte data) |
| EC50 | 1.5 uM | Estimated from dose-response |
| d (scaling) | 1.0 (default) | Conservative — no empirical calibration |
| **R3** | **0.226** | Well below 0.8 threshold |
| **R3 Flag** | **FLAGGED** | Clinically relevant induction |

The R3 value of 0.226 indicates that sotorasib at [I]max,u = 1.54 uM produces substantial CYP3A4 induction in vitro. The induction signal is strong and unambiguous.

## 5. The Net Effect Problem: TDI vs. Induction

This is the central analytical challenge of the sotorasib DDI profile. Both TDI and induction are independently flagged for CYP3A4, and they act in **opposing directions**:
- TDI reduces CYP3A4 activity (increases victim AUC)
- Induction increases CYP3A4 expression (decreases victim AUC)

### 5.1 Scenario Analysis

Using a hypothetical CYP3A4 substrate (fm = 0.94, fg = 0.43), the net AUCR depends critically on the induction scaling factor (d):

| Scenario | AUCR | Direction | Recommendation |
|----------|------|-----------|----------------|
| TDI only (no induction) | **12.39** | Inhibition | Strong inhibition -> clinical study |
| Induction only (no TDI) | **0.24** | Induction | Strong induction -> clinical study |
| Combined, d=1.0 (conservative) | **6.03** | Inhibition | Strong inhibition -> clinical study |
| Combined, d=0.5 (moderate) | **8.11** | Inhibition | Strong inhibition -> clinical study |
| Combined, d=0.2 (empirical) | **10.23** | Inhibition | Strong inhibition -> clinical study |
| **Observed clinical** | **0.48** | **Induction** | **Net induction confirmed** |

### 5.2 Critical Insight

The static model predicts **net inhibition** across all d-values tested (d = 0.2 to 1.0), yet the clinical result shows **net induction** (AUCR = 0.48, i.e., 52% decrease in midazolam AUC). This discrepancy reveals a fundamental limitation:

1. **The static model cannot resolve opposing mechanisms.** It applies peak concentrations to both TDI and induction simultaneously, but in vivo these processes operate on different timescales:
   - TDI is concentration-dependent and maximal at Cmax
   - Induction is a transcriptional process with onset over days and persistence beyond drug elimination

2. **The induction scaling factor d is poorly constrained.** The default d=1 assumes a direct translation from in vitro Emax to in vivo effect. In reality, the effective in vivo induction for sotorasib exceeds what d=1 captures, because:
   - Sotorasib accumulates in hepatocytes
   - The in vitro system may underestimate true induction potential
   - The 960 mg dose produces sustained high hepatocyte concentrations

3. **This is precisely the scenario where PBPK is essential.** The FDA 2020 guidance explicitly identifies dual TDI/induction as requiring dynamic modeling to resolve the net effect. Published PBPK analyses of sotorasib (e.g., FDA NDA 214665 clinical pharmacology review) correctly predict the net induction.

## 6. Transporter Screening

| Transporter | IC50 (uM) | Concentration | Ratio | Threshold | Flag |
|-------------|-----------|---------------|-------|-----------|------|
| **P-gp** | 3.0 | [I]gut = 6,850 | 2,283 | 10 | **FLAGGED** |
| **BCRP** | 0.54 | [I]gut = 6,850 | 12,685 | 10 | **FLAGGED** |
| **OATP1B1** | 0.29 | [I]max,u = 1.54 | 5.32 | 0.1 | **FLAGGED** |
| **OATP1B3** | 1.1 | [I]max,u = 1.54 | 1.40 | 0.1 | **FLAGGED** |
| **MATE1** | 1.0 | [I]max,u = 1.54 | 1.54 | 0.1 | **FLAGGED** |
| **MATE2-K** | 3.0 | [I]max,u = 1.54 | 0.51 | 0.1 | **FLAGGED** |
| OAT1 | 50.0 | [I]max,u = 1.54 | 0.031 | 0.1 | No |
| OAT3 | 50.0 | [I]max,u = 1.54 | 0.031 | 0.1 | No |
| OCT2 | 50.0 | [I]max,u = 1.54 | 0.031 | 0.1 | No |

**6 of 9** transporters flagged. Key clinical implications:
- **P-gp/BCRP inhibition (gut):** May increase oral bioavailability of P-gp/BCRP substrates (e.g., digoxin). The LUMAKRAS label includes warnings for P-gp substrates.
- **OATP1B1/1B3 inhibition (hepatic):** May increase AUC of OATP substrates (e.g., statins). The label reports a 1.48-fold increase in rosuvastatin AUC.
- **MATE1/2-K inhibition (renal):** May decrease renal secretion of cation substrates (e.g., metformin).

## 7. Integrated Risk Summary

| Mechanism | Status | Clinical Relevance |
|-----------|--------|-------------------|
| CYP3A4 reversible inhibition | Flagged (R1 = 1.15) | Weak |
| CYP3A4 TDI | **Severe** (TDI factor = 0.027) | Major — but offset by induction |
| CYP3A4 induction | **Strong** (R3 = 0.23) | Dominant in vivo |
| CYP2C8 inhibition | **Flagged** (R1 = 2.93) | Clinically confirmed |
| Transporter inhibition | 6/9 flagged | P-gp, BCRP, OATP1B1/3, MATE1/2-K |

### Regulatory Recommendations

1. **CYP3A4 net effect:** PBPK modeling is essential. The static model cannot resolve TDI vs. induction. Clinical data confirm net induction (AUCR = 0.48).
2. **CYP2C8 substrates:** Clinical DDI study or PBPK recommended. The R1 = 2.93 and label-confirmed inhibition warrant dose adjustment guidance for co-administered CYP2C8 substrates.
3. **P-gp/BCRP substrates:** Clinical DDI study warranted. The LUMAKRAS label already includes this recommendation.
4. **OATP substrates:** Monitor; the rosuvastatin clinical study showed a modest 1.48-fold increase.
5. **Renal transporters:** Monitor for co-administered MATE substrates.

## 8. Why This Case Matters for a Portfolio

This case demonstrates several competencies valued in pharma MIDD roles:

1. **Multi-mechanism complexity:** Real-world DDI assessment is rarely single-mechanism. Sotorasib exemplifies the modern challenge of compounds that simultaneously inhibit and induce the same enzyme.

2. **Knowing the limits of the tool:** The static model correctly flags all signals but cannot resolve the net direction. Recognizing this limitation and recommending PBPK is the appropriate scientific response — not forcing a single AUCR number.

3. **Regulatory alignment:** The analysis follows the exact framework of FDA 2020 DDI guidance Figure 7 (inhibition) and Figure 8 (induction), demonstrating fluency with regulatory methodology.

4. **Clinical validation:** Comparing predictions against the observed 52% AUC decrease provides a reality check that strengthens the analysis.

5. **Transporter completeness:** Screening all 9 transporters with appropriate concentration metrics ([I]gut for efflux, [I]max,u for uptake/renal) shows thoroughness expected in regulatory submissions.

## 9. Figures

| Figure | File | Description |
|--------|------|-------------|
| Mechanism dashboard | `figures/sotorasib_mechanism_dashboard.png` | 3-panel overview: CYP R1, R3 induction, transporter ratios |
| Net effect scenarios | `figures/sotorasib_net_effect_scenarios.png` | AUCR predictions across d-scaling scenarios vs. clinical observation |

## 10. Conclusion

Sotorasib is a paradigmatic example of a multi-mechanism DDI perpetrator where static screening correctly identifies all risk signals but cannot resolve the net CYP3A4 effect. The framework demonstrates:
- Comprehensive screening across CYP inhibition (6 enzymes), TDI, induction, and transporters (9)
- Transparent scenario analysis showing the sensitivity to induction scaling
- Clear identification of where PBPK modeling is needed to complement static screening
- Integration of clinical validation data to contextualize model predictions

This case showcases the analytical maturity to use static models as intended — as screening tools within a broader MIDD decision framework — rather than as standalone quantitative predictors.

---
*P10 Static DDI Risk Assessment Framework*
*FDA 2020 DDI Guidance | Mechanistic Static Model*
*Clinical reference: LUMAKRAS (sotorasib) NDA 214665*
