$PROBLEM Midazolam 2-Cpt PopPK | Oral Dosing | SAEM + Allometric WT
; -----------------------------------------------------------------------
; Structural:  2-compartment, first-order absorption (ADVAN4 TRANS4)
; Covariates:  Allometric WT on CL/Q (^0.75) and V2/V3 (^1.0)
; IIV:         Log-normal on Ka, CL, V2 (diagonal OMEGA)
; Residual:    Combined proportional + additive
; Estimation:  SAEM (burn-in 200, EM 100) → IMP evaluation → COV
;
; nlmixr2 counterpart: analysis/02_fit_models.R (fit_2comp)
; Dataset:     120 subjects, 720 observations, simulated
; NM-ready data created by: export_for_nonmem.R
; -----------------------------------------------------------------------

$INPUT ID TIME AMT EVID CMT MDV DV WT SEX AGE DROP DROP
;      columns:                                  ARM  BLQ

$DATA data_nm.csv
      IGNORE=@           ; skip header row

$SUBROUTINE ADVAN4 TRANS4
; TRANS4 parameterization: CL, V2, Q, V3, KA

; -----------------------------------------------------------------------
; $PK: Fixed effects (log-scale) + allometric scaling + IIV
; -----------------------------------------------------------------------
$PK

  ; --- Typical values with allometric scaling (ref = 70 kg) ---
  TVKA = EXP(THETA(1))
  TVCL = EXP(THETA(2)) * (WT/70)**0.75
  TVV2 = EXP(THETA(3)) * (WT/70)**1.0
  TVQ  = EXP(THETA(4)) * (WT/70)**0.75
  TVV3 = EXP(THETA(5)) * (WT/70)**1.0

  ; --- Individual values (log-normal IIV) ---
  KA = TVKA * EXP(ETA(1))
  CL = TVCL * EXP(ETA(2))
  V2 = TVV2 * EXP(ETA(3))
  Q  = TVQ                  ; no IIV on Q
  V3 = TVV3                 ; no IIV on V3

  S2 = V2                   ; scaling: F = A(2)/S2 → concentration

; -----------------------------------------------------------------------
; $ERROR: Combined proportional + additive
; -----------------------------------------------------------------------
$ERROR

  IPRED = F
  Y     = IPRED * (1 + EPS(1)) + EPS(2)

; -----------------------------------------------------------------------
; $THETA: Initial estimates (log-scale, from nlmixr2 final estimates)
; -----------------------------------------------------------------------
$THETA
  0.8153         ; 1  log(Ka)   → Ka   ≈ 2.26 /h
  3.9476         ; 2  log(CL/F) → CL/F ≈ 52   L/h  (at 70 kg)
  3.8417         ; 3  log(Vc/F) → Vc/F ≈ 46   L    (at 70 kg)
  2.8205         ; 4  log(Q/F)  → Q/F  ≈ 17   L/h  (at 70 kg)
  4.1356         ; 5  log(Vp/F) → Vp/F ≈ 63   L    (at 70 kg)

; -----------------------------------------------------------------------
; $OMEGA: Between-subject variability (diagonal, variances)
; -----------------------------------------------------------------------
$OMEGA
  0.0895         ; ETA(1) on Ka  (CV ≈ 30%)
  0.267          ; ETA(2) on CL  (CV ≈ 52%)
  0.0279         ; ETA(3) on V2  (CV ≈ 17%)

; -----------------------------------------------------------------------
; $SIGMA: Residual error variances
; -----------------------------------------------------------------------
$SIGMA
  0.0617         ; 1  Proportional variance (SD ≈ 0.2485, ~25%)
  0.2437         ; 2  Additive variance      (SD ≈ 0.494 ng/mL)

; -----------------------------------------------------------------------
; Estimation: SAEM → IMP objective function evaluation → Covariance
; -----------------------------------------------------------------------
$EST METHOD=SAEM NBURN=200 NITER=100 SEED=12345 PRINT=10 NOABORT
;    ISAMPLE=2 CTYPE=3  ; version-dependent — uncomment for NONMEM 7.4+

$EST METHOD=IMP EONLY=1 ISAMPLE=1000 NITER=5 PRINT=1
; IMP EONLY=1: evaluates -2LL at SAEM final estimates (no re-estimation)

$COV MATRIX=R UNCONDITIONAL PRINT=E

; -----------------------------------------------------------------------
; Standard diagnostic tables
; -----------------------------------------------------------------------
$TABLE ID TIME DV PRED IPRED CWRES EWRES
       ETA1 ETA2 ETA3 CL V2 KA WT
       NOPRINT ONEHEADER FILE=sdtab001

$TABLE ID CL V2 KA Q V3 ETA1 ETA2 ETA3 WT SEX AGE
       NOPRINT ONEHEADER FIRSTONLY FILE=patab001

; -----------------------------------------------------------------------
; Simulation (commented out — uncomment for VPC or predictive checks)
; -----------------------------------------------------------------------
; $SIM (12345) ONLYSIM NSUB=500

; -----------------------------------------------------------------------
; Alternative estimation: FOCEI (uncomment to compare with SAEM)
; -----------------------------------------------------------------------
; $EST METHOD=1 INTER MAXEVAL=9999 PRINT=5 NSIG=3 SIGL=9
; $COV MATRIX=R UNCONDITIONAL PRINT=E
