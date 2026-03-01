# Derivation Rules - NONMEM PK Analysis Dataset

## Study: CDISCPILOT01

## Time Variables

### TIME (hours from first dose)
- Reference time = datetime of first dosing record (EVID=1) per subject
- TIME = (record_datetime - reference_time) in hours
- All records sorted by USUBJID, DATETIME, desc(EVID) before derivation
- Dose records at same timepoint as observations are placed first

### TAD (time after most recent dose)
- For dose records (EVID=1): TAD = 0
- For observation records (EVID=0): TAD = TIME - TIME_of_last_dose
- Pre-dose observations: TAD = NA (flagged for review)

## Concentration Variables

### DV (dependent variable)
- Source: PC.PCSTRESN (standardized numeric result)
- Unit: ug/mL (standardized from PC.PCSTRESU)
- For dose records (EVID=1): DV = NA (displayed as "." in NONMEM output)

### BLQ Handling
- LLOQ = 0.01 ug/mL (study-specific, from PCLLOQ, defined in config)
- BLQ flag: 1 if PCSTRESN < LLOQ or PCSTRESN is NA (PCORRES = "<BLQ"), 0 otherwise
- BLQ rule applied: LLOQ/2 (DV set to 0.005 ug/mL for BLQ records)
- MDV remains 0 for BLQ records (included in analysis with LLOQ/2 imputation)

## Dosing Variables

### AMT (dose amount)
- Source: EX.EXDOSE
- Unit: mg (from EX.EXDOSU)
- For observation records (EVID=0): AMT = 0

### CUMDOSE (cumulative dose)
- Running sum of AMT for EVID=1 records within each subject
- Carried forward to observation records

## Event Flags

### EVID (event identifier)
- 1 = dosing event (from EX domain)
- 0 = observation event (from PC domain)

### MDV (missing dependent variable)
- 1 for dose records (EVID=1)
- 1 for observations with missing concentration
- 0 for observations with valid DV (including BLQ with LLOQ/2 imputation)

### CMT (compartment)
- 1 = dosing compartment (oral input)
- 2 = observation compartment (plasma)

## Covariate Derivations

### Baseline Covariates (time-invariant)
- AGE: from DM.AGE (years)
- SEX: recoded from DM.SEX (0=Female, 1=Male)
- RACE: recoded from DM.RACE (1=White, 2=Black, 3=Asian, 4=Other)

### Baseline Vitals (from VS domain, VSBLFL="Y")
- WT: VSTESTCD="WEIGHT", unit kg
- HT: VSTESTCD="HEIGHT", unit cm
- BMI: derived as WT/(HT/100)^2

### Baseline Labs (from LB domain, LBBLFL="Y")
- CREAT: serum creatinine (if available)
- ALT, AST, BILI: hepatic function markers (if available)

## Unit Harmonization

| Variable | Target Unit | Source | Conversion |
|----------|-------------|--------|------------|
| DV       | ug/mL       | PC.PCSTRESU | Verify match; convert if needed |
| AMT      | mg          | EX.EXDOSU   | Verify match; convert if needed |
| TIME/TAD | hours       | Datetime diff | difftime(..., units="hours") |
| WT       | kg          | VS.VSSTRESU  | Verify match |
| HT       | cm          | VS.VSSTRESU  | Verify match |
