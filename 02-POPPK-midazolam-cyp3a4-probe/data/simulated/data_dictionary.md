# Data Dictionary — analysis_dataset.csv

## Dataset Description

Publication-calibrated simulated trial dataset for oral midazolam 7.5 mg.
True PK parameters are anchored to the OSP Midazolam PBPK model evaluation
studies. This is NOT real clinical data.

## Study Design

- Single oral dose: midazolam 7.5 mg tablet
- N = 120 subjects (40 rich sampling + 80 sparse sampling)
- Population: healthy adults (age 22-72 y, weight 40-150 kg)
- LLOQ: 0.5 ng/mL
- BLQ handling: M1 method (excluded as MDV=1)

## Column Definitions

| Column | Type | Description |
|--------|------|-------------|
| ID | Integer | Subject identifier (1-120) |
| TIME | Numeric | Time since dose (hours) |
| AMT | Numeric | Dose amount (micrograms); 7500 for dosing records, 0 for observations |
| EVID | Integer | Event ID: 1 = dosing, 0 = observation |
| CMT | Integer | Compartment: 1 = depot (dose), 2 = central (observation) |
| MDV | Integer | Missing dependent variable: 1 = no DV (dose record or BLQ), 0 = DV present |
| DV | Numeric | Dependent variable: observed concentration (ng/mL); NA for dosing/BLQ records |
| IPRED | Numeric | Individual prediction without residual error (ng/mL); for QC only |
| WT | Numeric | Body weight (kg) |
| SEX | Integer | Sex: 0 = female, 1 = male |
| AGE | Numeric | Age (years) |
| ARM | Character | Sampling arm: "RICH" (10 samples) or "SPARSE" (4 samples) |
| BLQ | Integer | Below limit of quantification flag: 1 = BLQ, 0 = above LLOQ |

## Units

| Quantity | Unit |
|----------|------|
| Dose (AMT) | micrograms (ug) |
| Concentration (DV) | ng/mL (= ug/L) |
| Time | hours |
| Weight | kg |
| Age | years |

## True Simulation Parameters

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| Ka | 2.5 | 1/h | First-order absorption rate |
| CL/F | 50 | L/h | Apparent clearance (70 kg reference) |
| Vc/F | 45 | L | Apparent central volume (70 kg reference) |
| Q/F | 15 | L/h | Intercompartmental clearance |
| Vp/F | 55 | L | Peripheral volume |
| omega_Ka | 0.36 | variance | BSV on Ka (CV ~60%) |
| omega_CL | 0.09 | variance | BSV on CL (CV ~30%) |
| omega_Vc | 0.04 | variance | BSV on Vc (CV ~20%) |
| sigma_prop | 0.20 | SD | Proportional residual error (20%) |
| sigma_add | 0.5 | ng/mL | Additive residual error |
| LLOQ | 0.5 | ng/mL | Lower limit of quantification |

## Allometric Scaling (Built into Simulation)

- CL/F: (WT/70)^0.75
- Vc/F, Vp/F: (WT/70)^1.0
- Q/F: (WT/70)^0.75

## Sampling Schedules

| Arm | N | Samples (h post-dose) |
|-----|---|----------------------|
| RICH | 40 | 0.25, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12 |
| SPARSE | 80 | 0.5, 2, 4, 8 |
