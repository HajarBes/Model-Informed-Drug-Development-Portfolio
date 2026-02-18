# Data Provenance

## Source

- **Database:** Open-Systems-Pharmacology/Database-for-observed-data
- **File:** ObsDataPK_OSP.xlsx
- **Download date:** 2026-02-18
- **URL:** https://github.com/Open-Systems-Pharmacology/Database-for-observed-data

## Extraction Pipeline

Script: `analysis/08_extract_osp_profiles.py`

### Filtering Steps

| Step | Filter | Rows Remaining |
|------|--------|----------------|
| 0 | All PK-Profiles rows | 31098 |
| 1 | Study ID in midazolam PO human studies | 1952 |
| 2 | Compartment = Plasma | 1842 |
| 3 | Valid concentration units (ng/mL convertible) | 1791 |
| 4a | Mean profiles (excluding Heizmann individual) | 1653 |
| 4b | Heizmann 1983 individual curves | 138 |

### Summary

- Unique studies with extracted profiles: 164
- Total mean profile data points: 1653
- Heizmann 1983 individual data points: 138
- Excluded non-plasma compartments: Serum, Urine, Whole Blood
