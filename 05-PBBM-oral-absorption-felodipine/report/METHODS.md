# Methods

## 1. ACAT Model Structure

The Advanced Compartmental Absorption and Transit (ACAT) model divides the gastrointestinal tract into 7 segments:

| Segment | pH (fasted) | Transit time (h) | Surface area (cm²) | Volume (mL) | Bile factor |
|---------|-------------|-------------------|---------------------|-------------|-------------|
| Stomach | 1.7 | 0.25 | 0 (no absorption) | 250 | 1.0 |
| Duodenum | 6.0 | 0.25 | 150 | 50 | 2.0 |
| Jejunum 1 | 6.4 | 0.75 | 5400 | 100 | 1.5 |
| Jejunum 2 | 6.6 | 0.75 | 5400 | 100 | 1.2 |
| Ileum 1 | 7.0 | 1.00 | 3600 | 100 | 1.0 |
| Ileum 2 | 7.2 | 1.00 | 3600 | 100 | 1.0 |
| Colon | 6.5 | 18.0 | 900 | 500 | 1.0 |

Each segment tracks two drug states (undissolved solid, dissolved drug), yielding 14 GI states plus 2 systemic states (central and peripheral compartments) plus 1 cumulative absorption state = **17 ODEs total**. The 17th state tracks cumulative mass absorbed through the intestinal wall, providing a rigorous Fa calculation (see `data/sources/Fa_definition.md`).

## 2. Dissolution Model (Noyes-Whitney)

pH-dependent saturation solubility (Henderson-Hasselbalch for a weak base):

```
C_sat_i = S₀ × (1 + 10^(pKa - pH_i)) × bile_factor_i
```

Dissolution rate per segment (concentration-driving-force Noyes-Whitney):

```
diss_rate_i = z × (3 × D_eff) / (ρ_mg × r²) × S_i × (C_sat_i - C_diss_i)
```

Where:
- S₀ = 0.50 µg/mL (intrinsic solubility)
- pKa = 5.07 (weak base)
- D_eff = 5×10⁻⁶ cm²/s (diffusion coefficient)
- ρ_mg = 1300 mg/cm³ (particle density in mg/cm³ for correct unit balance)
- r = 25 µm (reference particle radius)
- z = 100 (dissolution scaling factor — accounts for fine amorphous precipitate formed when stomach-dissolved drug encounters higher-pH intestinal fluid)

The (C_sat − C_diss) driving force makes dissolution rate directly proportional to local saturation solubility. This correctly captures the bile salt effect: higher bile factors increase C_sat, directly accelerating dissolution rate — essential for modeling the food effect mechanism in BCS Class II drugs.

**Precipitation:** When C_diss > C_sat (supersaturation, e.g., after transit from stomach to intestine), drug precipitates at rate k_precip × (C_diss − C_sat) × V_segment, with k_precip = 10 h⁻¹.

## 3. Absorption Model

Absorption is driven by effective permeability:

```
k_abs_i = P_eff × SA_i / V_lumen_i
```

Where P_eff = 5×10⁻⁴ cm/s (high permeability, BCS Class II confirmed).

## 4. First-Pass Metabolism

Scalar extraction model:

```
drug_to_systemic = Σ(k_abs_i × D_i) × Fg × Fh
```

Where Fg = 0.50 (gut wall availability) and Fh = 0.50 (hepatic availability), per Lundahl et al. 1997.

## 5. Systemic Disposition

Two-compartment linear model:

```
dAc/dt = drug_to_systemic - (CL/Vc)×Ac - (Q/Vc)×Ac + (Q/Vp)×Ap
dAp/dt = (Q/Vc)×Ac - (Q/Vp)×Ap
Cp (ng/mL) = Ac / Vc × 1000
```

Parameters: CL = 70 L/h, Vc = 90 L, Q = 30 L/h, Vp = 200 L (Edgar 1987, Blychert 1990).

## 6. Food Effect Model

Three mechanistic knobs modify GI physiology in the fed state:

| Parameter | Fasted | Fed | Source | Mechanism |
|-----------|--------|-----|--------|-----------|
| Gastric emptying τ | 0.25 h | 0.65 h | **Calibrated** | Standard breakfast emptying (Bratel 1989 target) |
| Bile salt factor (duodenum) | 2× | 15× | Estimated | FaSSIF → FeSSIF bile secretion |
| Bile salt factor (jejunum) | 1.5× | 12× | Estimated | Micellar solubilization enhancement |
| Gastric pH | 1.7 | 3.0 | **Calibrated** | Preserves weak-base dissolution in stomach |

**Calibration note:** Fed gastric pH and emptying time were calibrated to reproduce the Bratel 1989 observed IR food effect (Cmax ratio ~1.31). Pre-calibration defaults (pH 4.0, τ = 1.5 h) yielded Cmax ratio 0.88, inconsistent with observed data. The calibrated values (pH 3.0, τ = 0.65 h) are within published physiological ranges for a standard breakfast: pH 3.0 preserves weak-base dissolution in the stomach (ionization = 118×, capacity = 29.6 mg > 10 mg dose), and 0.65 h is a standard-meal gastric emptying time (high-fat-meal default is 1.0–2.0 h).

**IR vs ER food effect:** The mechanistic decomposition reveals three competing effects: bile salts (+60% Cmax), delayed emptying (~−16%), and reduced gastric ionization (~−8%). For the IR formulation, bile salt enhancement dominates, yielding a net Cmax increase (~31%) with preserved AUC. This contrasts with the well-known +60% food effect for the ER formulation (Plendil), where matrix-controlled release eliminates the gastric pH mechanism and bile effects on dissolution dominate even more.

## 7. Model Qualification

Qualification against OSP database published mean PK profiles:

- **Acceptance criterion:** Predicted/observed AUC ratio within 2-fold (0.50–2.00)
- **Two-level classification:**
  - **in_database** — All extracted control arms (29 arms). These appear in plots.
  - **in_qualification_set** — IR/solution, fasted/unknown prandial state only. These are used for the X/Y within 2-fold metric. ER arms are excluded from qualification metrics because the IR ACAT model does not simulate extended-release dissolution.
- **DDI arms excluded:** Erythromycin and itraconazole perpetrator co-administration arms (2 arms)

**Dose normalization:** Studies dosed at strengths other than the 10 mg reference (e.g. 5 mg) are dose-normalized to 10 mg assuming linear PK:

```
C_obs_DN(t) = C_obs(t) × (D_ref / D_actual)
AUC_obs_DN  = AUC_obs  × (D_ref / D_actual)
```

Where D_ref = 10 mg and D_actual is the administered dose (e.g. 5 mg → DN factor = 2.0). The model predicts at 10 mg; observed profiles are scaled to the same reference. Predicted AUC is computed over the same time window as observed data points (not the full 0–24 h simulation) to ensure a fair comparison. This approach is appropriate for felodipine at therapeutic doses where PK is linear (no CYP3A4 saturation).

**Bioavailability breakdown:** See `outputs/tables/bioavailability_breakdown.csv`.

| Parameter | Predicted | Literature |
|-----------|-----------|------------|
| Fa (fraction absorbed) | 0.9715 | ~1.0 (BCS II, high Peff) |
| Fg (gut wall availability) | 0.50 | 0.50 (Lundahl 1997) |
| Fh (hepatic availability) | 0.50 | 0.50 (Lundahl 1997) |
| Foral (overall bioavailability) | 0.2429 | 0.15–0.22 (Lundahl 1997) |

The predicted Foral (24.3%) is slightly above the upper bound of the literature range (22%), consistent with the model's high Fa (~97%) combined with scalar first-pass extraction. The ~3% gap reflects the simplified first-pass model (constant Fg, Fh) vs. in vivo variability.

**Mass balance (fasted, 10 mg, 48 h):** See `outputs/tables/mass_balance_fasted.csv`.

| Component | Mass (mg) | % Dose |
|-----------|-----------|--------|
| Absorbed through intestinal wall | 9.715 | 97.15% |
| Remaining undissolved (GI) | ~0 | ~0% |
| Remaining dissolved (GI) | ~0 | ~0% |
| Fecal transit loss | 0.285 | 2.85% |
| **Total** | **10.000** | **100%** |

Mass balance closes to < 0.01 mg over 48 h simulation. The 17th ODE state tracks cumulative absorbed mass, providing a rigorous Fa = 0.9715 that matches the per-segment absorption sum (see `data/sources/Fa_definition.md`).

## 8. Sensitivity Analysis

Morris method (global screening) via SALib:
- **Parameters:** 15 (physicochemical + GI physiology + PK)
- **Trajectories:** 30 (yielding ~480 model evaluations)
- **Targets:** AUC and Cmax
- **Output:** μ* (mean absolute elementary effect) and σ (interaction/nonlinearity indicator)

## 9. Formulation Scenarios

Three particle size formulations evaluated:
- **Micronized:** r = 10 µm
- **Reference:** r = 25 µm
- **Coarser:** r = 50 µm

Virtual bioequivalence assessed across N=200 virtual subjects with LHS-sampled physiological variability. Geometric mean ratio (GMR) with 90% CI compared against 80–125% acceptance limits.

## 10. Numerical Methods

- **ODE solver:** scipy.integrate.solve_ivp with LSODA method
- **Tolerances:** rtol = 1×10⁻⁸, atol = 1×10⁻¹⁰
- **Time step:** dt = 0.01 h (36 seconds) for PK simulations; 0.05 h for screening
- **Mass balance verification:** Tolerance < 0.01 mg over 48 h simulation
