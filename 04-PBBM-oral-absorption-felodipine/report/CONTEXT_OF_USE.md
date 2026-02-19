# Context of Use

## Decision Question

**Primary:** Can a mechanistic ACAT-based oral absorption model, parameterized from published biopharmaceutic properties, predict observed felodipine plasma profiles and decompose the IR food effect into competing mechanistic drivers?

**Secondary:** Which formulation and physiological parameters are the dominant drivers of Cmax and AUC?

## Intended Application

This PBBM model is designed to:

1. **Explore food effect mechanisms** through mechanistic decomposition (bile salts, gastric emptying, gastric pH), calibrated against one observed IR study (Bratel 1989) — limited benchmarking, not a validated food-effect model
2. **Qualify against published clinical PK data** using the OSP database
3. **Explore dissolution-absorption linkage** for BCS Class II compounds
4. **Demonstrate formulation impact** through virtual bioequivalence analysis with varied particle size

## Data Foundation

| Source | Description | N |
|--------|-------------|---|
| OSP Database for Observed Data | Published felodipine mean ± SD plasma profiles | 408 data points, 31 studies |
| Published literature | Biopharmaceutic properties (MW, pKa, LogP, S₀, Peff) | Indexed references |
| FDA PLENDIL label | Food effect documentation | NDA 019834 |
| Edgar 1987, Blychert 1990 | Systemic PK parameters (CL, Vc, Q, Vp) | Clinical studies |

## Model Boundaries

| Component | Scope | Out of Scope |
|-----------|-------|-------------|
| Dissolution | Noyes-Whitney, pH-dependent, bile-enhanced | Supersaturation, amorphous, controlled-release |
| Absorption | Peff-driven, 7-segment ACAT | Active transport, efflux, regional Peff variation |
| First-pass | Scalar Fg × Fh | CYP3A4 concentration-dependent metabolism, DDI |
| Disposition | 2-compartment linear | Nonlinear PK, metabolite formation |
| Population | Log-normal LHS for virtual BE | Covariates (age, sex, genotype) |
| Food effect | 3 bulk physiological knobs | Detailed meal composition model |

## Regulatory Alignment

- **FDA Draft Guidance:** Physiologically Based Biopharmaceutics Modeling (PBBM), October 2023
- **EMA Guideline:** Qualification and Reporting of PBPK Modelling and Simulation (2018)
- **ICH M9:** Biopharmaceutics Classification System-Based Biowaivers (2019)

## Risk Assessment

| Risk | Mitigation |
|------|-----------|
| Model overpredicts absorption due to simplified dissolution | Mass balance verification; qualification against observed data |
| Scalar first-pass misses dose-dependent nonlinearity | Verified at 10 mg (well below CYP3A4 saturation); noted as limitation |
| Virtual BE uses simplified variability model | Clearly labeled as illustrative; not for regulatory filing |
| Food effect driven by bile salt factors may be overparameterized | Decomposition analysis separates individual contributions |
