# OSP Midazolam PBPK Model — Calibration Notes

Source: https://github.com/Open-Systems-Pharmacology/Midazolam-Model (v2.0, 2025)

## Model Hierarchy: PBPK to PopPK

The Open Systems Pharmacology (OSP) Midazolam PBPK model is a whole-body,
physiologically-based model validated across 40+ clinical studies and a
4-log dose range (0.001-40 mg, IV and oral). It mechanistically represents
CYP3A4-mediated hepatic and intestinal metabolism, UGT1A4 glucuronidation,
and renal elimination.

This PopPK project does NOT reproduce the PBPK model. Instead, it uses the
PBPK-validated parameter space as a **calibration anchor**: the "true"
population PK parameters used to simulate the trial dataset are constrained
to produce concentration-time profiles consistent with the PBPK evaluation
studies. This ensures the simulated data are publication-calibrated rather
than arbitrary.

## Key PBPK Parameters Used for Calibration

| Parameter | PBPK Value | PopPK Calibration | Rationale |
|-----------|-----------|-------------------|-----------|
| CL (hepatic) | ~25-30 L/h (CYP3A4 kcat=8.76/min, Km=4.0 uM) | CL/F = 50 L/h (70 kg) | F ~0.4-0.5 for oral tablet |
| Vss | ~70-100 L (Rodgers-Rowland) | Vc/F=45 L, Vp/F=55 L | Partitioned into central+peripheral |
| Ka | Intestinal perm = 1.55e-4 cm/min, dissolution T50=0.01 min | Ka = 2.5 /h | Effective first-order approximation |
| F | ~0.3-0.5 (combined hepatic + gut CYP3A4 first-pass) | Absorbed into CL/F, V/F | Oral model estimates apparent parameters |
| fu | 0.031 (97% albumin-bound) | Not modeled explicitly | Absorbed into apparent CL/F |

## Primary Enzyme: CYP3A4

- Hepatic CYP3A4: kcat = 8.76 /min, Km = 4.0 uM (Meyer 2012 expression)
- Intestinal CYP3A4: same kinetics, tissue-specific expression
- UGT1A4: secondary pathway (kcat = 3.59 /min, Km = 37.8 uM)
- Renal: GFR fraction = 0.64 (minor, <1% of dose)

## Physicochemical Properties

| Property | Value |
|----------|-------|
| MW | 325.78 g/mol |
| logP | 2.90 |
| pKa (base) | 6.2 |
| pKa (acid) | 10.95 |
| fu plasma | 0.031 |
| Aqueous solubility | 0.13 mg/mL (pH 5.0) |

## Population Variability Anchoring

The OSP model evaluates PK across diverse populations (European, Korean),
age groups, and dose levels. Published PopPK analyses of midazolam report:

- CL/F CV: 30-45%
- V/F CV: 25-40%
- Ka CV: 50-80%
- WT effect on CL: allometric exponent ~0.75 (established)
- WT effect on V: allometric exponent ~1.0 (established)

These variability estimates inform the BSV (omega) values used in the
simulated trial dataset.
