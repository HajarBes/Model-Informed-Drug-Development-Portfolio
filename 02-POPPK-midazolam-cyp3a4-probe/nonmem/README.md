# NONMEM Control Stream — Midazolam 2-Compartment PopPK

This directory contains a NONMEM 7.5-compatible control stream that mirrors
the nlmixr2 model structure: 2-compartment oral, first-order absorption,
allometric WT scaling, and combined residual error. Results may differ
slightly from nlmixr2 due to estimator implementation differences.

## Files

| File | Description |
|------|-------------|
| `midazolam_2cmt_saem.ctl` | NONMEM control stream (SAEM + IMP + COV) |
| `export_for_nonmem.R` | R script to create NM-ready dataset from `analysis_dataset.csv` |
| `data_nm.csv` | NM-ready dataset (created by `export_for_nonmem.R`) |

## How to Run

```bash
# 1. Create NM-ready dataset
Rscript export_for_nonmem.R

# 2. Run NONMEM
nmfe75 midazolam_2cmt_saem.ctl midazolam_2cmt_saem.lst
```

## Parameter Mapping

| Parameter | nlmixr2 Name | NONMEM | Initial Value | Interpretation |
|-----------|-------------|--------|---------------|----------------|
| Ka | `lka` | THETA(1) | 0.8153 | log(Ka) → Ka ≈ 2.26 /h |
| CL/F | `lcl` | THETA(2) | 3.9476 | log(CL/F) → CL/F ≈ 52 L/h at 70 kg |
| Vc/F | `lvc` | THETA(3) | 3.8417 | log(Vc/F) → Vc/F ≈ 46 L at 70 kg |
| Q/F | `lq` | THETA(4) | 2.8205 | log(Q/F) → Q/F ≈ 17 L/h at 70 kg |
| Vp/F | `lvp` | THETA(5) | 4.1356 | log(Vp/F) → Vp/F ≈ 63 L at 70 kg |
| BSV Ka | `eta.ka` | ETA(1) / OMEGA(1,1) | 0.0895 | CV ≈ 30% |
| BSV CL | `eta.cl` | ETA(2) / OMEGA(2,2) | 0.267 | CV ≈ 52% |
| BSV Vc | `eta.vc` | ETA(3) / OMEGA(3,3) | 0.0279 | CV ≈ 17% |
| Prop. error | `prop.sd` | EPS(1) / SIGMA(1,1) | 0.0617 | SD ≈ 0.25 (~25%) |
| Add. error | `add.sd` | EPS(2) / SIGMA(2,2) | 0.2437 | SD ≈ 0.49 ng/mL |

## Estimation Strategy

1. **SAEM** (METHOD=SAEM): stochastic approximation EM, 200 burn-in + 100 iterations
2. **IMP** (METHOD=IMP, EONLY=1): importance sampling evaluation of -2LL at SAEM estimates
3. **$COV**: R-matrix covariance for standard errors

## Data Preparation Notes

`export_for_nonmem.R` transforms the analysis CSV for NONMEM compatibility:
- Drops `IPRED` column (not a NONMEM input)
- Encodes `ARM`: RICH → 1, SPARSE → 2
- Replaces `NA` with `.` (NONMEM missing-value convention)
- Writes space-delimited format with an `@`-prefixed header line (e.g., `@ID TIME AMT ...`) so `IGNORE=@` in the `.ctl` skips it

## Notes

- Initial estimates (THETAs, OMEGAs, SIGMAs) are taken from nlmixr2 final estimates
- Allometric exponents are fixed at theoretical values (0.75 for clearances, 1.0 for volumes)
- IIV is estimated on Ka, CL, and V2 only (Q and V3 have no random effects)
- Commented-out `$SIM` and `$EST METHOD=1` (FOCEI) blocks are included for convenience
