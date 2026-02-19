# Executive Summary

## PBBM Oral Absorption Model — Felodipine (BCS Class II)

### Objective

To build and qualify a physiologically based biopharmaceutics model (PBBM) that predicts felodipine oral absorption, explores the food effect through mechanistic decomposition (with limited benchmarking against one observed IR study), and quantifies the impact of formulation particle size on bioavailability.

### Approach

A 7-segment ACAT (Advanced Compartmental Absorption and Transit) model was constructed using published biopharmaceutic properties of felodipine — a BCS Class II drug with low solubility and high permeability. The model couples Noyes-Whitney dissolution kinetics with Peff-driven intestinal absorption, scalar first-pass metabolism (Fg × Fh), and 2-compartment systemic disposition. All parameters were sourced from peer-reviewed literature without fitting to observed plasma data.

### Key Findings

1. **Fasted-state qualification:** The ACAT model predicts fasted-state felodipine PK profiles consistent with published clinical data from the OSP database. For IR/solution arms (model-comparable), 9/9 studies fall within 2-fold predicted/observed AUC ratio.

2. **Food effect mechanistic exploration (limited benchmarking):** Fed-state gastric pH (3.0) and emptying time (0.65 h) were calibrated to match the single available observed IR food effect study (Bratel 1989, Cmax ratio ~1.31). The model decomposes the food effect into three competing mechanisms: (a) bile salt enhancement (+60% Cmax), (b) delayed gastric emptying (~−16%), and (c) reduced gastric ionization (~−8%). For the **IR formulation**, the net effect is a Cmax increase (ratio ≈ 1.31) with preserved AUC. This is a mechanistic exploration calibrated to one study, not a validated food-effect prediction. It is mechanistically distinct from the +60% food effect documented for the **ER formulation** (Plendil).

3. **Sensitivity analysis:** Morris global screening of 15 parameters identified clearance, gut wall availability (Fg), hepatic availability (Fh), dose, and permeability as the top drivers of both Cmax and AUC — consistent with the rapid absorption kinetics of this highly permeable weak base.

4. **Formulation scenarios:** Virtual bioequivalence analysis across a 200-subject population demonstrated that particle size reduction from 25 µm (reference) to 10 µm (micronized) increases Cmax (GMR 1.26, failing 80-125% criterion), while coarser particles (50 µm) decrease AUC (GMR 0.83, within 80-125%), illustrating the dissolution-absorption linkage for BCS Class II compounds.

### Data

- **408 data points** from **31 studies** in the OSP Database for Observed Data
- 29 arms in database (control arms); 2 DDI treatment arms excluded
- Qualification set (IR/solution, fasted): used for pred/obs 2-fold metrics
- ER arms in database but not in qualification set (informational overlay)
- Observed profiles dose-normalized to 10 mg reference (linear PK assumption)
- Published biopharmaceutic properties from PubChem, DrugBank, and indexed literature

### Regulatory Context

This analysis aligns with the FDA Draft Guidance on PBBM (2023) and EMA PBPK reflection paper (2018). The model demonstrates the mechanistic linkage between formulation properties (particle size, dissolution), GI physiology, and systemic exposure — a framework applicable to biowaivers, food effect exploration, and clinically relevant dissolution specifications.

### Limitations

This is a methodological demonstration, not a regulatory submission. Key simplifications include monodisperse particles (no size distribution), scalar first-pass extraction (no CYP3A4 saturation), and bulk food effect knobs. See `report/LIMITATIONS_AND_SCOPE.md` for complete discussion.
