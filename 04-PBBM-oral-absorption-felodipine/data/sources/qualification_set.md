# Qualification Set — Felodipine PBBM Model

## Positioning Statement

All observed data are **arm-level mean +/- SD profiles** (digitized) from the
OSP Database for Observed Data. This work supports **arm-level / central-tendency
qualification only**. No claims are made about individual-level predictive
performance, inter-individual variability (IIV) validation, or covariate
validation.

## Selection Criteria

| Criterion | Rule |
|-----------|------|
| Compound | Felodipine |
| Route | Oral (PO) |
| Species | Human |
| Compartment | Plasma only (serum excluded) |
| Arm type | Control/baseline only |
| Min time points | >= 3 per profile |

## Exclusion Criteria

| Criterion | Rule |
|-----------|------|
| DDI perpetrator | Grapefruit juice, itraconazole, erythromycin co-admin |
| Non-plasma | Serum, urine, whole blood compartments |
| mg/kg dosing | Excluded unless convertible with documented assumptions |
| Unclear arms | Arms with ambiguous grouping labels |

## Two-Level Classification

Each arm is classified at two levels:

1. **in_database** — Extracted from OSP, control arm, >=3 time points.
   These arms appear in the data set and can be plotted.

2. **in_qualification_set** — Used for pred/obs 2-fold AUC ratio metrics.
   Requires: in_database AND IR/solution formulation AND fasted/unknown prandial state.
   ER arms are NOT in the qualification set because the ACAT model simulates IR dissolution.

- Arms in database: 29
- Arms in qualification set: 9
- ER arms (in database, not in qualification set): 18
- Excluded from database: 2

## Qualification Set (used for 2-fold metrics)

| Study ID | Study | Grouping | Form. | Dose | N | Points | Prandial | Notes |
|----------|-------|----------|-------|------|---|--------|----------|-------|
| 16002 | Edgar 1987 | 10mg-PO solution | solution | 10.0 mg | 12 | 13 | unknown | — |
| 16003 | Edgar 1987 | 10mg-PO IR tablet | IR | 10.0 mg | 12 | 15 | unknown | — |
| 16007 | Blychert 1990 | Study I-10mg-PO solution | solution | 10.0 mg | 18 | 14 | unknown | multi-dose |
| 16008 | Blychert 1990 | Study I-10mg-PO IR tablet-trial | IR | 10.0 mg | 18 | 14 | unknown | multi-dose |
| 16009 | Blychert 1990 | Study I-10mg-PO IR tablet-market | IR | 10.0 mg | 18 | 14 | unknown | multi-dose |
| 16010 | Blychert 1990 | Study II-10mg-PO IR tablet-market | IR | 10.0 mg | 15 | 14 | unknown | multi-dose |
| 16019 | Smith 1987 | 10mg-PO IR tablet | IR | 10.0 mg | 8 | 10 | fasted | multi-dose |
| 16027 | Bailey 1993 | 5mg-PO IR tablet | IR | 5.0 mg | 9 | 11 | fasted | — |
| 16029 | Edgar 1992 | 5mg-PO IR tablet | IR | 5.0 mg | 9 | 11 | fasted | — |

## In Database but NOT in Qualification Set

| Study ID | Study | Form. | Dose | Prandial | Reason not in qualification set |
|----------|-------|-------|------|----------|--------------------------------|
| 16000 | Lundahl 1997 | ER | 10.0 mg | fasted | ER formulation |
| 16004 | Edgar 1987 | ER | 10.0 mg | unknown | ER formulation |
| 16005 | Bailey 1996 | ER | 10.0 mg | fasted | ER formulation |
| 16006 | Jalava 1997 | unknown | 5.0 mg | fasted | unknown formulation |
| 16011 | Blychert 1990 | ER | 10.0 mg | unknown | ER formulation |
| 16012 | Blychert 1990 | ER | 10.0 mg | unknown | ER formulation |
| 16014 | Bailey 2003 | ER | 10.0 mg | fasted | ER formulation |
| 16015 | Bailey 2000 | ER | 10.0 mg | fasted | ER formulation |
| 16016 | Bailey 1998 | ER | 10.0 mg | fasted | ER formulation |
| 16017 | Guo 2007  | ER | 10.0 mg | fasted | ER formulation |
| 16018 | Dresser 2002 | ER | 10.0 mg | fasted | ER formulation |
| 16020 | Madsen 1996 | ER | 10.0 mg | fasted | ER formulation |
| 16021 | Lundahl 1998 | ER | 10.0 mg | fasted | ER formulation |
| 16022 | Dresser 2017 | ER | 10.0 mg | fasted | ER formulation |
| 16023 | Bailey 1995 | ER | 10.0 mg | fasted | ER formulation |
| 16024 | Aberg 1997 | ER | 10.0 mg | fasted | ER formulation |
| 16025 | Dresser 2000 | ER | 5.0 mg | fasted | ER formulation |
| 16026 | Dresser 2000 | ER | nan mg | fasted | ER formulation |
| 16028 | Goosen 2004 | ER | 5.0 mg | fasted | ER formulation |
| 16030 | Bratel 1989 | IR | 5.0 mg | fed | fed prandial state |

## Excluded from Database

| Study ID | Study | Grouping | Reason |
|----------|-------|----------|--------|
| 16062 | Bailey 1996 | Erythromycin DDI, treatment group | DDI perpetrator co-administration arm; ER formulation — not directly comparable to IR ACAT model |
| 16063 | Jalava 1997 | Itraconazole DDI, treatment group | DDI perpetrator co-administration arm; Multi-dose study (time offset > 48 h) |

## Summary

- Total arms in OSP: 31
- In database (extracted control arms): 29
- In qualification set (IR/solution, fasted): 9
- Excluded from database: 2
- Data points in database: 384
- Data points in qualification set: 116
