# Limitations and Scope

## Model Limitations

### 1. Simplified ACAT: 7 Segments, Monodisperse Particles

The model uses a 7-segment ACAT framework with monodisperse spherical particles. Real formulations have particle size distributions that could alter dissolution kinetics. The constant particle radius assumption (no shrinking core) is conservative — actual dissolution may be faster as particles shrink.

**Impact:** May slightly underpredict dissolution rate for polydisperse formulations.

### 2. Scalar First-Pass (Fg, Fh): No CYP3A4 Saturation

First-pass metabolism is modeled as fixed fractions (Fg = 0.50, Fh = 0.50) independent of drug concentration. In reality, CYP3A4 saturation could increase bioavailability at higher doses or when absorption rate increases (as in the fed state). This simplification is reasonable at the 10 mg dose level but may not extend to higher doses.

**Impact:** Food effect on F_oral may be slightly underestimated if CYP3A4 saturation occurs in the fed state.

### 3. Observed Data Are Aggregated (Mean ± SD)

All observed data from the OSP database are published group-level summary statistics (arithmetic mean ± SD or SEM). Individual-level pharmacokinetic data are not available. This prevents individual-level validation and limits assessment of between-subject variability.

**Impact:** Model qualification compares predicted population-typical profile against study means; individual variability is only explored through virtual population simulations.

### 4. Food Effect: Three Bulk Physiological Knobs (Two Calibrated)

The food effect is modeled by modifying three GI parameters (gastric emptying, bile salts, gastric pH) rather than through a detailed meal digestion model. Fed gastric pH (3.0) and emptying time (0.65 h) were **calibrated** to match the Bratel 1989 IR observed food effect (Cmax ratio ~1.31); bile factors were estimated from FaSSIF/FeSSIF literature. The calibrated values are within published physiological ranges for a standard breakfast. The model omits:
- Luminal fluid volume dynamics
- Viscosity and motility pattern changes
- Fat-mediated lymphatic absorption (minor for felodipine)
- Meal composition-dependent effects

**Impact:** The model matches the single available paired IR food effect study (Bratel 1989) but the calibration is based on one observed data point. Generalization to other meal types or doses is uncertain.

### 5. Virtual BE Is Illustrative

The virtual bioequivalence analysis uses simplified physiological variability (log-normal LHS) and does not include:
- Crossover study design simulation
- Period effects or sequence effects
- Intra-subject variability
- Regulatory-compliant statistical methods (e.g., ANOVA on log-transformed data with proper degrees of freedom)

**Impact:** Results are qualitatively informative for dissolution-exposure linkage but should not be used for regulatory filing.

### 6. Permeability Assumed Uniform

Effective permeability (Peff) is applied uniformly across all absorbing segments. In reality, transporter expression and surface morphology vary along the GI tract. However, for high-permeability BCS Class II drugs like felodipine, absorption is dissolution-limited, making this assumption reasonable.

**Impact:** Negligible for felodipine; would matter more for BCS Class IV or actively transported drugs.

### 7. Precipitation and Dissolution Scaling

The model includes a first-order precipitation term (k_precip = 10 h⁻¹) for when dissolved concentration exceeds local C_sat (e.g., stomach-dissolved drug entering higher-pH intestine). A dissolution scaling factor (z = 100) accounts for the fine amorphous precipitate formed during pH transition having higher effective surface area than the original formulation particles. The exact precipitate morphology (amorphous vs. nanocrystalline) is not tracked.

**Impact:** The z-factor is a calibration parameter; sensitivity analysis can quantify its effect on PK predictions.

### 8. IR Model vs ER Food Effect Literature

The ACAT model simulates **IR (immediate release)** dissolution only. The well-documented +60% Cmax food effect for felodipine refers to the **ER formulation (Plendil)**, not IR. For IR, the calibrated model predicts a +31% Cmax increase (matching Bratel 1989), driven by bile salt enhancement outweighing the moderate gastric emptying delay. ER formulations, where the matrix controls release rate independently of gastric pH, are expected to show a larger positive food effect driven by bile salt enhancement alone.

**Impact:** The food effect predictions are specific to the IR formulation. ER comparison would require a modified release module.

## Scope Statement

This project is a **methodological demonstration** of PBBM principles for portfolio purposes. It is:

- **NOT** a regulatory filing or submission-ready analysis
- **NOT** validated against individual-level clinical data
- **NOT** intended for clinical dose recommendations

The model demonstrates the author's capability to:
1. Build mechanistic oral absorption models from first principles
2. Connect formulation properties to systemic exposure through dissolution-absorption coupling
3. Apply global sensitivity analysis to identify critical quality attributes
4. Generate publication-quality figures and regulatory-style documentation
5. Work with real published clinical data (OSP database)
