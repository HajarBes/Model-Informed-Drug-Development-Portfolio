# Regulatory-Style Report (Portfolio): KRAS G12C QSP Model

## 1. Context of Use (COU)
**Purpose:** Mechanistically evaluate adaptive resistance to KRAS G12C inhibition and explain differential clinical response between NSCLC and CRC, and the rationale for EGFR blockade combination in CRC.

**Intended audience:** QSP / MIDD reviewers (portfolio demonstration).

**Key questions:**
1) Why does KRAS-G12C monotherapy show stronger activity in NSCLC than CRC?
2) Can EGFR blockade reduce rebound signaling and improve tumor control in CRC?

## 2. Model Overview
- PK: oral 1-compartment
- Target engagement: KRAS(GDP) covalent binding
- Signaling: receptor-driven feedback -> ERK proxy rebound
- Tumor: growth modulated by pathway output

See: model/equations.md

## 3. Data Sources & Traceability
Primary local sources:
- Sotorasib label: data/raw/pk/sotorasib_LUMAKRAS_FDA_label_2025.pdf
- Panitumumab label/PI: data/raw/pk/panitumumab_VECTIBIX_FDA_label.pdf ; data/raw/pk/panitumumab_VECTIBIX_PI.pdf
- CodeBreaK 100 (NSCLC): data/raw/clinical/NEJM_CodeBreaK100_sotorasib_NSCLC_2021.pdf
- CodeBreaK 300 (mCRC combo): data/raw/clinical/NEJM_CodeBreaK300_sotorasib_panitumumab_mCRC_2023.pdf
- Review context: data/raw/references/Riedl_2026_KRAS_landscape_Cancer_Cell.pdf

### Traceability table (fill)
| Quantity | Value | Units | Source PDF | Where used |
|---|---:|---|---|---|
| Dose (sotorasib) | 960 | mg QD | LUMAKRAS label | PK module |
| Half-life | ~5 | h | LUMAKRAS label | PK ke |
| ORR NSCLC mono | 37.1 | % | CodeBreaK 100 | tumor calibration target |
| mPFS NSCLC mono | 6.8 | months | CodeBreaK 100 | tumor timescale check |
| ORR mCRC combo | 26.4 | % | CodeBreaK 300 | combo calibration target |
| mPFS mCRC combo | 5.6 | months | CodeBreaK 300 | tumor timescale check |
| Panitumumab dose | 6 | mg/kg q2w | VECTIBIX PI | combo module |

## 4. Verification (Did we implement what we intended?)
- Unit checks (parameters, time units)
- Scenario reproducibility (one command generates figs + CSV)
- Sanity checks (steady state PK, TE within [0,100])

## 5. Validation (Fit-for-purpose checks)
**Qualitative validation targets:**
- NSCLC: smaller rebound, better tumor control vs CRC mono
- CRC: stronger rebound, weaker tumor control
- CRC combo: rebound suppression and improved tumor control

**Quantitative anchors (current):**
- ORR/PFS summary targets from clinical trials (no individual patient data)

## 6. Sensitivity & Uncertainty Plan (to add)
- Global sensitivity on rebound amplitude and tumor shrinkage
- Virtual population sampling of key parameters

## 7. Assumptions & Limitations
- ERK represented as a proxy output
- WT-RAS and parallel pathways not explicit
- No biomarker time-course fitting yet (pERK etc.)
- CRC monotherapy paper not included yet (gap)

## 8. Readiness & Gaps
**Ready:** mechanistic prototype demonstrates hypothesis.
**Gaps to reach “strong portfolio”:**
- Add CRC monotherapy clinical anchor
- Add sensitivity analysis + VPop
- Optional: digitize biomarker curves (pERK) from literature
