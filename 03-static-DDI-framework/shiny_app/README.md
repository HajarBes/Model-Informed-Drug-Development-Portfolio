# DDI Static Screening Calculator — Interactive Shiny App

Interactive tool for mechanistic static DDI risk assessment per FDA 2020 In Vitro DDI Guidance.

## Quick Start

```r
# From the project root (03-static-DDI-framework/)
shiny::runApp("shiny_app")

# Or from any directory
shiny::runApp("path/to/03-static-DDI-framework/shiny_app")
```

## Requirements

- R >= 4.4.0
- Packages: `shiny`, `bslib`, `DT`, `plotly`, `readr`, `dplyr`, `tidyr`, `purrr`, `tidyverse`, `scales`

Install all at once:
```r
install.packages(c("shiny", "bslib", "DT", "plotly", "tidyverse", "scales"))
```

## Features

### Presets
Load the ketoconazole benchmark (Case A) or sotorasib (Case B) from validated CSV input files. Fields auto-populate and lock. Click "Unlock" to customize.

### Tab 1: Screening Results (Core)
- CYP enzyme screening table with R1, R1,gut, TDI factor, and flag status (color-coded)
- Induction screening: R3 value and flag
- AUCR breakdown: hepatic, gut, and total contributions
- FDA classification badge (Strong / Moderate / Weak / No interaction)
- Regulatory recommendation text

### Tab 2: Sensitivity
- Tornado plot: one-at-a-time parameter variation showing AUCR impact
- Ki vs [I]max,u heatmap with compound position marked (white X)

### Tab 3: Monte Carlo
- Configurable N (100-10,000) and CV (10-60%)
- Button-triggered simulation (prevents UI freeze)
- AUCR distribution histogram colored by FDA classification
- Summary: median, 90% CI, percentage in each classification

### Tab 4: Transporter Dashboard
- Horizontal bar chart of [I]/IC50 ratio vs threshold per transporter
- Color-coded: green (below threshold), red (flagged)
- Supports all 9 transporters (P-gp, BCRP, OATP1B1/1B3, MATE1, MATE2-K, OAT1/3, OCT2)

### Download
Export a CSV report of all screening results via the Download button.

## Architecture

Single-file `app.R` that sources computation functions from `analysis/00_setup.R`. All DDI equations (R1, R1,gut, TDI factor, R3, AUCR) are reused from the validated analysis pipeline — no duplicated logic.

## Regulatory Basis

- FDA In Vitro Drug Interaction Studies Guidance (2020)
- ICH M12 Drug Interaction Studies (2024)
- Screening thresholds: R1 >= 1.02, R1,gut >= 11, TDI factor < 0.5, R3 <= 0.8
- AUCR classification: >= 5 (Strong), >= 2 (Moderate), >= 1.25 (Weak), < 1.25 (No interaction)
