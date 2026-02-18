# OSP Database Audit: Midazolam Oral PK Data Candidates

**Source:** [Open-Systems-Pharmacology/Database-for-observed-data](https://github.com/Open-Systems-Pharmacology/Database-for-observed-data)
**File:** `ObsDataPK_OSP.xlsx` (downloaded 2026-02-18)
**Filter:** Compound = Midazolam, Route = PO, Species = Human

---

## Summary

| Metric | Count |
|--------|-------|
| Total entries (midazolam PO) | 208 |
| Unique studies | 65 |
| Studies with control/baseline arms | ~50 |
| Studies with **individual-level** data | 2 |
| Studies with **aggregated** (mean ± SD) data | ~63 |
| Studies with plasma concentration-time profiles | ~60 |

---

## Critical Finding

**Nearly all data in the OSP database is aggregated (mean ± SD).** Only 2 studies
are flagged as "Individual" data type (Heizmann 1983, Greenblatt 1986), and even
these contain digitized individual curves from published figures — not raw
bioanalytical data with covariates.

**Implication for PopPK:**
- Structural parameters (Ka, CL/F, Vc/F, Q/F, Vp/F) CAN be estimated from
  mean profiles (naive pooled or two-stage approach)
- Between-subject variability (BSV/IIV) CANNOT be estimated from mean ± SD alone
- Covariate analysis (WT, SEX, AGE) is NOT possible — no individual covariates
- Full PopPK workflow (NONMEM/nlmixr2 with random effects) requires
  individual-level data that this database does not contain

---

## Grade Definitions

| Grade | Meaning | What it supports |
|-------|---------|------------------|
| **A** | Individual-level data (multiple IDs, covariates) | Full PopPK with IIV + covariate analysis |
| **B** | Mean ± SD profiles, usable for structural fit | Parameter estimation, model structure selection, no IIV |
| **C** | Not suitable (DDI-only, urine, no profiles, too sparse) | Cannot use |

---

## Grade A Candidates: Suitable for PopPK

### NONE

No study in the OSP database provides individual-level concentration-time data
with subject-level covariates (WT, SEX, AGE) for oral midazolam in healthy adults.

The two studies flagged as "Individual" data type are:

#### Heizmann 1983 (closest to Grade A)
- **Ref:** https://pubmed.ncbi.nlm.nih.gov/6138080/
- **Data type:** Individual (digitized per-subject curves from publication)
- **N:** 6 subjects, each dosed at 10, 20, and 40 mg PO (crossover)
- **Timepoints:** 13 per curve, range 0.25–12 h
- **Compartment:** Plasma
- **Covariates available:** NONE (no WT, SEX, AGE in database)
- **Limitation:** 6 subjects is too few for meaningful BSV estimation. No
  covariates. Individual curves are digitized from figures, not raw data.
- **GRADE: B** (useful for structural fit with dose-linearity, not for IIV)

#### Greenblatt 1986
- **Ref:** https://pubmed.ncbi.nlm.nih.gov/2935051/
- **Data type:** Individual
- **N:** 6 healthy males
- **Dose:** 15 mg PO
- **Timepoints:** 12 per curve, range 0.5–10 h
- **Compartment:** Plasma
- **Covariates:** NONE in database
- **Context:** Cimetidine DDI study — control arm usable
- **GRADE: B** (N=6, no covariates, insufficient for PopPK)

---

## Grade B Candidates: Suitable for Structural Re-fit (No IIV)

Ranked by utility (N, sampling richness, dose, relevance):

### Tier 1 — Best for structural parameter estimation

| Study | N | Dose (mg) | Timepoints | Range (h) | Fasted | Notes |
|-------|---|-----------|------------|-----------|--------|-------|
| **Wiesinger 2020** | 65 | 1 | 14 | 0.5–24 | — | Largest N; control arm of rifampicin DDI; all female |
| **Gorski 2003** | 52 | 4/6 | 13 | 0.08–10 | Fasted | Large N; **whole blood** (not plasma) |
| **Darwish 2008** | 24 | 5 | 15 | 0–24 | Fasted | Standalone; multiethnic (8 white, 16 black) |
| **Lutz 2018a** | 20 | 2 | 9 | 0.25–12 | Fasted | 2 control cohorts; 50% female |
| **Mueller 2009** | 20 | 7.5 | 16 | 0.08–12 | Fasted | Rich sampling; pre-SJW baseline |
| **Zahner 2019** | 20 | 1 | 18 | 0–72 | — | Very rich (18 tp); Day 1 control |
| **Bornemann 1986** | 18 | 15 | 13–15 | 0.25–12 | Both | 4 food conditions; good for Ka |
| **Gorski 1998** | 16 | 4 | 11 | 0.25–8 | Fasted | **Whole blood**; control arm |
| **Hohmann 2015** | 16 | 3 | 13 | 0.25–24 | Fasted | Two dose levels (0.003 and 3 mg) |
| **Kharasch 2011** | 12 | 3 | 21 | 0–9 | — | Very rich (21 tp); control arm |
| **Chattopadhyay 2018** | 11 | 1 | 15 | 0.5–24 | — | Control arm; female only |

### Tier 2 — Standard DDI control arms (7.5 or 15 mg)

| Study | N | Dose (mg) | Timepoints | Range (h) | Notes |
|-------|---|-----------|------------|-----------|-------|
| **Olkkola 1993** | 12 | 15 | 12 | 0–18 | Finland; 9F/3M |
| **Olkkola 1996** | 12 | 7.5 | 11 | 0–17 | Finland; 2 control arms |
| **Yeates 1996** | 12 | 15 | 10 | 0.5–24 | Clarithromycin DDI control |
| **Ahonen 1995** | 12 | 7.5 | 11 | 0–17 | Finland; itraconazole DDI control |
| **Zimmermann 1996** | 12 | 15 | 9 | 0.5–12 | Erythromycin DDI control |
| **Backman 1996** | 10 | 15 | 10 | 0.5–10 | WT/AGE in paper (not in DB) |
| **Saari 2006** | 10 | 7.5 | 9 | 0–8 | Finland; voriconazole DDI control |
| **Greenblat 1984** | 10 | 10 | 12–13 | 0.17–24 | "Typical" profiles by age/sex/obesity subgroups |
| **Olkkola 1994** | 9 | 7.5 | 11 | 0–17 | Finland; ketoconazole DDI control |

### Tier 3 — Usable but limited

| Study | N | Dose (mg) | Timepoints | Range (h) | Limitation |
|-------|---|-----------|------------|-----------|------------|
| Backman 1994 | 9 | 15 | 10 | 0.5–17 | All female |
| Allonen 1981 | 6 | 15 | 12 | 0.17–7 | N=6, short range |
| Smith 1981 | 6 | 10 | 11–13 | 0.08–8 | Tablet vs solution |
| Salonen 1986 | 6 | 15 | 10 | 0–5 | Very short range |
| Fee 1987 | 8 | 15 | 14 | 0.25–9 | N=8 |
| Markert 2013 | 11 | 3 | 4 | 2–4 | Only 4 timepoints |

---

## Grade C: Not Suitable

| Study | Reason |
|-------|--------|
| Reitman 2011 | DDI arms only (post-rifampicin recovery), no baseline |
| Hyland 2009 | Urine data only (1 timepoint at 24h) |
| Thummel 1996 | Urine data only |
| Carls 2014 | Microdose (0.003 mg) with erythromycin, no standalone control |
| Björkhem-Bergman 2013 | No concentration-time profiles in database |
| Wang 2005 | No concentration-time profiles in database |
| Scholz 2021 | Urine + sparse plasma; SJW induction only |
| Gurley 2002 | Post-induction only, no baseline profile |
| Gurley 2005 | Post-induction only, no baseline profile (elderly) |
| Lutz 2018b | Carbamazepine arm only, no control profile |
| Martinez 1999 | Cimetidine arm only, no control profile |

---

## Recommended Strategy for Literature-Derived Upgrade

### What the OSP database gives us:
- **~35 usable mean ± SD concentration-time profiles** across control/baseline arms
- Dose range: 0.075–40 mg (mostly 1–15 mg)
- Rich sampling in many studies (9–21 timepoints per profile)
- Plasma (most studies) and whole blood (Gorski 1998, 2003)

### What it does NOT give us:
- Individual-level data (ID-level concentration-time records)
- Subject-level covariates (WT, SEX, AGE per subject)
- Raw bioanalytical data (all digitized from publications)

### Recommendation:

**Option 1: Multi-study pooled structural fit (RECOMMENDED)**
Pool control/baseline mean profiles from 10–15 best studies spanning
multiple doses (1–15 mg). Fit structural 2-compartment model to mean
data using NLS or naive pooled approach. This replaces the simulated
dataset's structural parameters with literature-calibrated values.
**Limitation:** No IIV, no covariate analysis.

**Option 2: Two-stage approach**
Extract published PK parameter estimates (AUC, CL/F, Cmax) from the
PK-Parameter sheet for ~30+ studies. Use these to inform prior
distributions for BSV. Combined with mean profile fit, this provides
literature-grounded structural + variability parameters.
**Limitation:** Indirect IIV estimation, still no individual covariates.

**Option 3: Keep simulated data, cite OSP calibration explicitly**
The current approach (OSP PBPK-calibrated simulation) remains the most
honest way to demonstrate a **complete PopPK workflow** including IIV,
covariate analysis, and simulation-based decision-making. Literature
mean profiles cannot support these analyses.
**Advantage:** Full workflow demonstration intact.
**Enhancement:** Add a validation step comparing model-predicted mean
profile against the best OSP literature profiles (external visual
predictive check).

### Bottom line:
The OSP database is a goldmine for **structural parameter validation**
but contains **no individual-level data suitable for population PK modeling
with random effects.** A hybrid approach (Option 3 + structural validation
against literature means) is the strongest defensible position.
