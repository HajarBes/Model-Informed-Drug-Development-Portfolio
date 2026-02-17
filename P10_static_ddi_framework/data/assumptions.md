# Parameter Sources and Assumptions

> Every parameter used in this framework is documented here with value, units, source, and notes.

---

## Case A: Midazolam + Ketoconazole (Benchmark)

### Perpetrator: Ketoconazole

| Parameter | Value | Units | Source | Notes |
|-----------|-------|-------|--------|-------|
| Dose (oral) | 400 | mg | FDA DDI Guidance (2020) | Standard DDI study dose |
| Molecular weight | 531.43 | g/mol | PubChem CID 456201 | |
| Cmax,total | 8.07 | ug/mL | Huang et al. (2004) Clin Pharmacol Ther | After 400 mg single oral dose |
| fu,plasma | 0.01 | — | FDA DDI Guidance (2020) | ~99% protein bound |
| **Derived: [I]max,u** | **0.152** | **uM** | Calculated | = (Cmax/MW) x 1000 x fu |
| **Derived: [I]gut** | **3012** | **uM** | Calculated | = (Dose_ug / 250 mL) / MW x 1000 |

### CYP Inhibition Parameters (Ketoconazole)

| CYP | Ki (uM) | Source | kinact (1/min) | KI (uM) | Source (TDI) |
|-----|---------|--------|----------------|---------|-------------|
| CYP3A4 | 0.015 | FDA DDI Guidance (2020), Table 7 | 0.048 | 0.86 | Mayhew et al. (2000) Drug Metab Dispos |
| CYP2C9 | 10.0 | Niwa et al. (2005) | — | — | — |
| CYP2C19 | 3.4 | Niwa et al. (2005) | — | — | — |
| CYP2D6 | >100 | FDA DDI Guidance (2020) | — | — | Negligible |
| CYP1A2 | >100 | FDA DDI Guidance (2020) | — | — | Negligible |

### Victim: Midazolam

| Parameter | Value | Units | Source | Notes |
|-----------|-------|-------|--------|-------|
| fm,CYP3A4 (hepatic) | 0.94 | — | Gorski et al. (2003) Clin Pharmacol Ther | |
| Fg (fraction escaping gut) | 0.57 | — | Thummel et al. (1996) Clin Pharmacol Ther | fg = 1 - Fg = 0.43 |

### Enzyme Degradation Rate

| Parameter | Value | Units | Source | Notes |
|-----------|-------|-------|--------|-------|
| kdeg,CYP3A4 (hepatic) | 0.00032 | 1/min | Yang et al. (2008) | t1/2 ~ 36 h |
| kdeg,CYP3A4 (gut) | 0.00048 | 1/min | Galetin et al. (2008) | t1/2 ~ 24 h |

### Observed Clinical DDI (Validation)

| Metric | Observed Value | Source |
|--------|---------------|--------|
| Midazolam AUC ratio | ~10-15x | Olkkola et al. (1994) Clin Pharmacol Ther |
| Classification | Strong inhibitor | FDA DDI Guidance |

---

## Case B: Sotorasib (LUMAKRAS) DDI Assessment

### Perpetrator: Sotorasib

| Parameter | Value | Units | Source | Notes |
|-----------|-------|-------|--------|-------|
| Dose (oral) | 960 | mg | LUMAKRAS FDA Label (2021) | Dose used in DDI studies |
| Molecular weight | 560.6 | g/mol | PubChem CID 137278711 | |
| Cmax,total | 7.87 | ug/mL | LUMAKRAS FDA Label, Clinical Pharmacology | Geometric mean at 960 mg QD |
| fu,plasma | 0.11 | — | LUMAKRAS FDA Label | 89% protein bound |
| **Derived: [I]max,u** | **1.54** | **uM** | Calculated | = (7.87/560.6) x 1000 x 0.11 |
| **Derived: [I]gut** | **6850** | **uM** | Calculated | = (960000/250) / 560.6 x 1000 |

### CYP Inhibition Parameters (Sotorasib as Perpetrator)

| CYP | Ki_rev (uM) | Source | kinact (1/min) | KI (uM) | Source (TDI) | Notes |
|-----|-------------|--------|----------------|---------|-------------|-------|
| CYP3A4 | 10.0 | Estimated from label | 0.04 | 3.5 | Estimated from FDA review | TDI confirmed in label |
| CYP2C8 | 0.8 | Estimated (IC50 ~1.6 uM) | — | — | — | Label confirms clinically relevant inhibition |
| CYP2C9 | >50 | LUMAKRAS Label | — | — | — | Not clinically relevant |
| CYP2C19 | >50 | LUMAKRAS Label | — | — | — | Not clinically relevant |
| CYP2D6 | >50 | LUMAKRAS Label | — | — | — | Not clinically relevant |
| CYP1A2 | >50 | LUMAKRAS Label | — | — | — | Not clinically relevant |

**Note on CYP3A4 TDI parameters:** The kinact and KI values are estimated from the FDA clinical pharmacology review (NDA 214665) and published characterization. Exact values should be verified against the original review document. Sensitivity analysis is performed over a plausible range.

### CYP Induction Parameters (Sotorasib)

| CYP | Emax (fold) | EC50 (uM) | Source | Notes |
|-----|-------------|-----------|--------|-------|
| CYP3A4 | 8.1 | 1.5 | Estimated from FDA review | ~8x induction at 3 uM in hepatocytes; positive control rifampin ~20x |

**Note:** Emax and EC50 are estimated from reported fold-induction data in cultured human hepatocytes. The scaling factor d = 1 is used (conservative).

### Transporter Inhibition Parameters (Sotorasib)

| Transporter | IC50 (uM) | Screening conc. | Source | Notes |
|-------------|-----------|-----------------|--------|-------|
| P-gp | 3.0 | [I]gut | LUMAKRAS Label | Intestinal efflux transporter |
| BCRP | 0.54 | [I]gut | LUMAKRAS Label | Intestinal efflux transporter |
| OATP1B1 | 0.29 | [I]max,u | LUMAKRAS Label | Hepatic uptake transporter |
| OATP1B3 | 1.1 | [I]max,u | LUMAKRAS Label | Hepatic uptake transporter |
| MATE1 | 1.0 | [I]max,u | LUMAKRAS Label | Renal transporter |
| MATE2-K | 3.0 | [I]max,u | LUMAKRAS Label | Renal transporter |
| OAT1 | >50 | — | LUMAKRAS Label | Not clinically relevant |
| OAT3 | >50 | — | LUMAKRAS Label | Not clinically relevant |
| OCT2 | >50 | — | LUMAKRAS Label | Not clinically relevant |

### Observed Clinical DDI (Validation)

| Interaction | Observed | Source |
|-------------|----------|--------|
| Sotorasib on midazolam (CYP3A4) | AUC decreased 52% | LUMAKRAS Label |
| Net CYP3A4 effect | Weak inducer (induction > TDI) | LUMAKRAS Label |
| Sotorasib + ketoconazole (as victim) | AUC increased 2.1x | LUMAKRAS Label |

---

## Transporter Screening Thresholds (FDA 2020)

| Transporter | Relevant [I] | Threshold | Action if exceeded |
|-------------|-------------|-----------|-------------------|
| P-gp (intestinal) | [I]gut / IC50 | >= 10 | Clinical study recommended |
| BCRP (intestinal) | [I]gut / IC50 | >= 10 | Clinical study recommended |
| OATP1B1/1B3 (hepatic) | [I]max,u / IC50 | >= 0.1 | Clinical study recommended |
| OCT/OAT/MATE (renal) | [I]max,u / IC50 | >= 0.1 | Clinical study recommended |

---

## Key Assumptions

1. [I]gut calculated as dose (ug) / 250 mL, per FDA guidance assumption for intestinal lumen volume.
2. fu,plasma used as proxy for fu,blood unless fu,blood is specifically available.
3. Ki estimated as IC50/2 where only IC50 is reported (competitive inhibition assumption).
4. kdeg for CYP3A4 hepatic = 0.00032 min-1 (t1/2 = 36 h) used for all cases.
5. Scaling factor d = 1 for induction (conservative; no empirical calibration applied).
6. Sensitivity analysis ranges: parameters varied +/- 50% from baseline unless otherwise specified.
7. Monte Carlo: lognormal distributions with CV = 30% for most parameters.

---

## Primary References

1. FDA Guidance for Industry: In Vitro Drug Interaction Studies (2020)
2. FDA Guidance for Industry: Clinical Drug Interaction Studies (2020)
3. LUMAKRAS (sotorasib) FDA Label, NDA 214665 (2021, updated)
4. Huang SM et al. (2004) Clin Pharmacol Ther
5. Gorski JC et al. (2003) Clin Pharmacol Ther
6. Mayhew BS et al. (2000) Drug Metab Dispos
7. Yang J et al. (2008) Curr Drug Metab
