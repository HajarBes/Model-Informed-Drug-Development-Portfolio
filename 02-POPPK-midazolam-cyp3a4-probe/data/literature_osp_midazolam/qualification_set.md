# Qualification Set: OSP Literature Profiles for Structural Qualification

## Inclusion Criteria

1. **Route:** Oral (PO)
2. **Compartment:** Plasma only (whole blood studies excluded: Gorski 1998, 2003)
3. **Arm:** Control or baseline arm only (no DDI treatment arms)
4. **Timepoints:** >= 10 per profile, spanning both absorption and elimination phases
5. **Dose:** Fixed mg dose (mg/kg excluded)
6. **Species:** Human
7. **Data type:** Aggregated (mean +/- SD or SEM)
8. **Source grading:** Grade B Tier 1 or Tier 2 from `study_candidates.md`

## Selected Profiles (N = 12)

| ID | Study | N | Dose (mg) | Timepoints | Time Range (h) | Fasted/Fed | Dose Group |
|----|-------|---|-----------|------------|----------------|------------|------------|
| 203 | Wiesinger 2020 | 65 | 1 | 14 | 0.5-24 | -- | low |
| 16647 | Zahner 2019 | 20 | 1 | 14 | 0-72 | -- | low |
| 1361 | Chattopadhyay 2018 | 11 | 1 | 15 | 0.5-24 | -- | low |
| 264 | Hohmann 2015 | 16 | 3 | 13 | 0.25-24 | Fasted | low-mid |
| 303 | Kharasch 2011 | 12 | 3 | 21 | 0-9 | -- | low-mid |
| 118 | Darwish 2008 | 24 | 5 | 14 | 0-24 | Fasted | mid |
| 16500 | Mueller 2009 | 20 | 7.5 | 15 | 0.08-12 | Fasted | mid |
| 374 | Olkkola 1996 | 12 | 7.5 | 11 | 0-17 | Fasted | mid |
| 49 | Ahonen 1995 | 12 | 7.5 | 11 | 0-17 | Fasted | mid |
| 107 | Bornemann 1986 | 18 | 15 | 14 | 0.25-12 | Fasted | high |
| 365 | Olkkola 1993 | 12 | 15 | 12 | 0-18 | Fasted | high |
| 53 | Backman 1996 | 10 | 15 | 10 | 0.5-10 | Fasted | high |

## Dose Groups

| Group | Doses (mg) | Studies | Rationale |
|-------|-----------|---------|-----------|
| low | 1 | Wiesinger, Zahner, Chattopadhyay | Microdose CYP3A4 phenotyping range |
| low-mid | 3 | Hohmann, Kharasch | Low therapeutic range |
| mid | 5, 7.5 | Darwish, Mueller, Olkkola 1996, Ahonen | Standard probe dose range |
| high | 15 | Bornemann, Olkkola 1993, Backman | Classic midazolam PO dose |

## Exclusions

| Study | Reason |
|-------|--------|
| Gorski 1998, 2003 | Whole blood (not plasma) |
| Bornemann 1986 ID 104, 106 | Fed conditions (model has no food effect) |
| Bornemann 1986 ID 105 | "1 h before a meal" — ambiguous fasting state |
| All DDI treatment arms | Only control/baseline arms qualify |
| mg/kg dose studies | Cannot convert without individual body weight |

## Notes

- All data are aggregated (study-level mean +/- SD or SEM), digitized from published figures
- No individual-level concentrations or subject-level covariates are available
- This qualification set spans a 15-fold dose range (1-15 mg) with N = 10-65 per study
- The 7.5 mg dose used in the simulated PopPK dataset is directly represented by 3 studies
