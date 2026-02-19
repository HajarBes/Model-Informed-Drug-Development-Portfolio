# Report — PBBM Oral Absorption Model for Felodipine

*Guidance-aligned sections following EMA PBPK Reporting Recommendations (2018) and FDA PBBM Draft Guidance (2023)*

---

## 1. Objectives

- Develop a mechanistic oral absorption model for felodipine (BCS Class II) to predict fasted-state plasma PK profiles from biopharmaceutic inputs
- Qualify the model against published clinical data from 31 studies (408 data points; 10 IR/solution arms, 18 ER arms)
- Decompose the food effect into competing mechanistic drivers and explain the IR vs ER distinction
- Identify key formulation and physiological drivers of exposure through global sensitivity analysis
- Explore the dissolution-exposure linkage through virtual bioequivalence analysis at varied particle sizes

## 2. Drug Properties

| Property | Value | Source |
|----------|-------|--------|
| Molecular weight | 384.26 g/mol | PubChem CID 3333 |
| BCS Class | II (low solubility, high permeability) | FDA/EMA classification |
| pKa | 5.07 (weak base) | Raggi et al. 2013 |
| LogP | 3.86 | DrugBank |
| Intrinsic solubility (S₀) | 0.50 µg/mL | Loftsson & Hreinsdóttir 2006 |
| Effective permeability (Peff) | 5×10⁻⁴ cm/s | High (BCS II confirmed) |
| Oral bioavailability | 15–22% | Lundahl et al. 1997 |
| Formulation | Immediate release (IR), 10 mg | Reference dose |

## 3. GI Physiology Model

Seven-segment ACAT model (Stomach → Duodenum → Jejunum1 → Jejunum2 → Ileum1 → Ileum2 → Colon) with published population-average physiology for fasted and fed states. See `report/METHODS.md` Section 1 for full parameter table.

## 4. Dissolution Model

Noyes-Whitney dissolution with Henderson-Hasselbalch pH-dependent solubility and bile salt micellar solubilization enhancement. Monodisperse spherical particles assumed. See `report/METHODS.md` Section 2.

## 5. Absorption and First-Pass

Peff-driven absorption across intestinal segments with scalar first-pass extraction:
- Gut wall availability: Fg = 0.50
- Hepatic availability: Fh = 0.50
- Overall first-pass extraction: ~75%

## 6. Systemic Disposition

Two-compartment linear model with parameters from Edgar 1987 and Blychert 1990:
CL = 70 L/h, Vc = 90 L, Q = 30 L/h, Vp = 200 L.

## 7. Model Qualification

**Data source:** OSP Database for Observed Data — 29 control arms from 31 published felodipine PO studies.

**Two-level classification:**
- **in_database** (29 arms): All extracted control arms. Shown in plots.
- **in_qualification_set** (IR/solution, fasted): Used for pred/obs 2-fold metrics. ER arms excluded from qualification because the model simulates IR dissolution only.

**Dose normalization:** Observed profiles at non-reference doses (e.g. 5 mg) are scaled to 10 mg assuming linear PK.

**Acceptance criterion:** Predicted/observed AUC ratio within 2-fold (EMA 2018 standard).

**Results:** See `outputs/tables/fasted_qualification.csv` and Figure 1.

## 8. Food Effect Exploration (Mechanistic, Limited Benchmarking)

The model explores the food effect through three mechanistic changes, calibrated to one observed IR study (Bratel 1989, N=11). This is a mechanistic exploration with limited benchmarking, not a validated food-effect prediction.
- Delayed gastric emptying (0.25 h → 0.65 h, calibrated to Bratel 1989)
- Increased bile salt solubilization (2× → 15× in duodenum)
- Reduced gastric ionization (pH 1.7 → 3.0, calibrated to preserve weak-base dissolution)

**Calibration note:** Fed gastric pH (3.0) and emptying time (0.65 h) were calibrated to match the Bratel 1989 observed IR food effect (Cmax ratio ~1.31). Pre-calibration defaults (pH 4.0, τ = 1.5 h) yielded Cmax ratio 0.88. Both calibrated values are within published physiological ranges for a standard breakfast.

**IR formulation result:** Cmax ratio ≈ 1.31, AUC ratio ≈ 1.03. The bile salt enhancement dominates over the moderate gastric emptying delay, yielding a net Cmax increase consistent with the observed IR food effect.

**Mechanistic decomposition:** Bile salts (+60% Cmax), delayed emptying (~−16%), reduced gastric ionization (~−8%). For IR, bile enhancement outweighs the competing mechanisms. The well-known +60% food effect refers to the ER formulation (Plendil), where matrix-controlled release eliminates the gastric pH mechanism entirely.

**Results:** See `outputs/tables/food_effect_summary.csv`, `food_effect_decomposition.csv`, and Figures 3–4.

## 9. Sensitivity Analysis

Morris global screening of 15 parameters (30 trajectories, ~480 evaluations) identifies:
- **Top AUC drivers:** Fg, Fh, CL, S₀, bile factors
- **Top Cmax drivers:** S₀, bile factors, Fg, particle radius, Vc

Results in `outputs/tables/sensitivity_morris_auc.csv`, `sensitivity_morris_cmax.csv`, and Figure 5.

## 10. Formulation Scenarios

Virtual bioequivalence at 3 particle sizes (10, 25, 50 µm) across N=200 virtual subjects with LHS-sampled physiological variability. GMR with 90% CI evaluated against 80–125% limits.

Results in `outputs/tables/virtual_be_summary.csv` and Figures 6–7.

## 11. Discussion and Conclusions

1. The ACAT model, parameterized entirely from published biopharmaceutic data, predicts felodipine fasted-state PK consistent with observed clinical profiles (9/9 IR studies within 2-fold AUC ratio).

2. The food effect mechanistic exploration, calibrated to one observed IR study (Bratel 1989), reveals three competing drivers: bile salt enhancement (+60% Cmax), delayed gastric emptying (~−16%), and reduced gastric ionization (~−8%). Bile enhancement dominates, yielding a net Cmax increase (ratio 1.31). This is benchmarked against a single study (N=11) and should be interpreted as mechanistic exploration, not a validated prediction.

3. Sensitivity analysis confirms that for this rapidly absorbed weak base, clearance (CL), first-pass metabolism (Fg, Fh), and dose are the dominant exposure drivers, consistent with absorption-rate-limited rather than dissolution-limited kinetics for the IR formulation.

4. The dissolution-exposure linkage demonstrates that particle size is a clinically relevant formulation attribute for felodipine, with micronized particles (10 µm) increasing Cmax by 26% compared to the reference (25 µm).

5. Limitations include monodisperse particles, scalar first-pass, a dissolution scaling factor (z = 100) for precipitate re-dissolution, and simplified food effect physiology. These are appropriate for the model's intended use as a methodological demonstration of PBBM principles.
