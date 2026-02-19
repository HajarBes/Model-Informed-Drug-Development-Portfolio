# Data Provenance — Felodipine PK Profiles

## Source

- **Database:** Open-Systems-Pharmacology/Database-for-observed-data
- **File:** ObsDataPK_OSP.xlsx
- **Extraction date:** 2026-02-19
- **URL:** https://github.com/Open-Systems-Pharmacology/Database-for-observed-data

## Extraction Pipeline

Script: `analysis/01_extract_osp_felodipine.py`

### Filtering Steps

| Step | Filter | Rows Remaining |
|------|--------|----------------|
| 0 | All PK-Profiles rows | 31098 |
| 1 | Study ID in felodipine PO human studies | 422 |
| 2 | Compartment = Plasma | 408 |
| 3 | Valid concentration units (ng/mL convertible) | 408 |
| 4 | Control/baseline arms only | 384 |

### Summary

- Unique studies with extracted profiles: 31
- Control arm studies (qualification set): 29
- Total data points: 408
- Control arm data points: 384
- Excluded compartments: Serum

### DDI Arms Excluded

Perpetrator co-administration arms were excluded from the qualification set:
- Study 16062: Bailey 1996 — Erythromycin DDI, treatment group
- Study 16063: Jalava 1997 — Itraconazole DDI, treatment group
