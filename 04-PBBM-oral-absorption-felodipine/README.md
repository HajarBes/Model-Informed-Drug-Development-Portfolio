# PBBM Oral Absorption Model — Felodipine (BCS Class II)

A physiologically based biopharmaceutics model (PBBM) demonstrating mechanistic oral absorption modeling: formulation + dissolution + physiology → absorption → systemic PK.

## Decision Question

**Can a mechanistic ACAT-based oral absorption model, parameterized from published biopharmaceutic properties, predict observed felodipine plasma profiles and mechanistically decompose the food effect into competing physiological drivers?**

## Model Overview

- **Drug:** Felodipine (BCS Class II — low solubility, high permeability)
- **Architecture:** 7-segment ACAT model with 17 ODE states
- **Dissolution:** Noyes-Whitney with pH-dependent (Henderson-Hasselbalch) and bile-enhanced solubility
- **Absorption:** Peff-driven across intestinal segments
- **First-pass:** Scalar Fg × Fh (Lundahl 1997)
- **Disposition:** 2-compartment linear PK (Edgar 1987, Blychert 1990)
- **Food effect:** 3 mechanistic knobs (gastric emptying, bile salts, gastric pH)

## Data

- **408 data points** from **31 studies** in the OSP Database for Observed Data
- 29 arms included (10 IR/solution + 18 ER); 2 DDI treatment arms excluded
- Published physicochemical properties from PubChem, DrugBank, and indexed literature

## Key Results

| Analysis | Finding |
|----------|---------|
| Fasted qualification | 9/9 IR studies within 2-fold AUC ratio; Cmax = 8.95 ng/mL |

Model qualification achieved 9/9 IR arms within 2-fold AUC ratio under mean-level comparison. No parameter was tuned per study.
| Food effect (IR) | Cmax ratio 1.31, AUC ratio 1.03 — bile enhancement dominates |
| Mechanistic decomposition | Bile salts +60%, gastric emptying −16%, gastric pH −8% (three competing effects) |
| Dominant drivers | CL, Fg, Fh, dose, Peff are top sensitivity drivers |
| Formulation impact | Micronized (10 µm): Cmax GMR 1.26; Coarser (50 µm): AUC GMR 0.83 |

## Figures

| # | Figure | Description |
|---|--------|-------------|
| 1 | `pk_overlay_fasted.png` | Predicted vs observed fasted-state PK profiles |
| 2 | `absorption_by_segment_fasted.png` | Regional absorption along GI tract |
| 3 | `pk_fasted_vs_fed.png` | Food effect on PK (fasted vs fed overlay) |
| 4 | `food_effect_decomposition.png` | Individual mechanism contributions to food effect |
| 5 | `tornado_sensitivity.png` | Morris global sensitivity analysis (AUC + Cmax) |
| 6 | `dissolution_exposure_linkage.png` | In vitro dissolution → in vivo PK |
| 7 | `formulation_comparability.png` | Virtual BE forest plot (3 particle sizes) |
| 8 | `drug_solubility_profile.png` | pH-dependent solubility characterization |

## Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Ensure OSP data is available (from project 02)
ls ../02-POPPK-midazolam-cyp3a4-probe/data/literature_osp/ObsDataPK_OSP.xlsx

# Run full pipeline
bash run_all.sh
```

## Project Structure

```
analysis/
  00_setup.py                     — Constants, GI physiology, drug params, helpers
  01_extract_osp_felodipine.py    — OSP data extraction (openpyxl)
  01b_qualify_and_food_effect.py  — Qualification set audit + observed food effect
  acat_model_02.py                — ACAT ODE system (importable module)
  03_simulate_fasted.py           — Fasted PK simulation + qualification
  04_simulate_fed.py              — Fed PK + food effect analysis
  05_sensitivity_analysis.py      — Morris screening (SALib)
  06_formulation_scenarios.py     — Particle size impact + virtual BE
  07_generate_figures.py          — Consolidated publication figures

data/
  observed/extracted/             — Extracted CSV profiles
  sources/                        — Citations, provenance, assumptions

figures/                          — 8 publication-quality PNGs
outputs/tables/                   — CSV summaries
report/                           — Regulatory-style documentation
```

## Tools

**Python 3.10+** with NumPy, SciPy, Matplotlib, Pandas, SALib, openpyxl

## Domain

Biopharmaceutics, Oral Absorption, PBBM, BCS Class II, Food Effect, Dissolution Modeling

## Author

Hajar Besbassi
