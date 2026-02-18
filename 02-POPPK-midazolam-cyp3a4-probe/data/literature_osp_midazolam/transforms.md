# Unit Conversion and Data Transforms

## Concentration Unit Normalization

All concentrations are normalized to ng/mL using the following factors:

| Source Unit | Conversion Factor | Formula |
|------------|-------------------|----------|
| mg/L | 1000.0 | conc * 1000.0 |
| ng/L | 0.001 | conc * 0.001 |
| ng/mL | 1.0 | conc * 1.0 |
| ng/ml | 1.0 | conc * 1.0 |
| nmol/L | 0.3258 | conc_nmol/L * MW / 1000 (MW = 325.78) |
| nmol/l | 0.3258 | conc_nmol/L * MW / 1000 (MW = 325.78) |
| pg/mL | 0.001 | conc * 0.001 |
| ug/L | 1.0 | conc * 1.0 |
| µg/L | 1.0 | conc * 1.0 |

MW (Midazolam) = 325.78 g/mol

## Variability Conversion

| Source VarType | Conversion to SD |
|---------------|------------------|
| arith. SD | Used directly (with unit conversion) |
| arith. SEM | SD = SEM * sqrt(N), where N from Studies sheet |
| geom. SD | Not converted (geometric SD not directly comparable) |
| 95th CI | Not converted |

## Dose Parsing

- Dose strings parsed from Studies sheet column 10
- Parenthetical notes stripped: '15 (actually 7.5)' → 15
- mg/kg doses excluded (cannot convert without individual body weight)
- Dose unit 'µg' converted to mg by dividing by 1000
