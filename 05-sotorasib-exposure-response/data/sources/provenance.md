# Data Provenance — Project 05: Sotorasib Exposure-Response

## Data Flow

```
Published Clinical Data (peer-reviewed + FDA review)
        │
        ▼
Published Aggregate Endpoints ──────────────────────┐
  • GM AUC/Cmax by dose (FDA review, Day 8)         │
  • ORR by dose (CodeBreaK 100 Phase 1/2)           │
  • Safety rates (pooled + dose comparison)          │
  • PopPK summary stats (Nagase et al. 2025)         │
        │                                            │
        ▼                                            ▼
Virtual Patient Simulation              Calibration Targets
  • N=500 per dose level                • PK: GM AUC/Cmax
  • Log-normal exposure sampling        • Efficacy: ORR by dose
  • Covariate assignment                • Safety: Grade 3+ rate
        │                                            │
        ▼                                            │
E-R Analysis & Figures ◄─── Anchor vs Simulated ─────┘
  • Logistic regression (ORR, safety)
  • Quartile analyses (model-implied)
  • Benefit-risk comparison
  • Covariate forest plot
```

## Source-to-Output Mapping

| Output | Primary Data Source | Calibration Target |
|--------|--------------------|--------------------|
| virtual_patients.csv | PopPK params (Nagase 2025, FDA label) | GM AUC/Cmax by dose |
| orr_by_exposure_quartile.csv | CodeBreaK ORR by dose | ORR ~30% (Phase 2) |
| safety_by_exposure_quartile.csv | Pooled safety rates (FDA review) | 16% Grade 3+ (960 mg) |
| dose_optimization_summary.csv | Dose comparison (EJC 2024) | Published ORR/PFS/OS |
| covariate_effects.csv | PopPK covariate analysis (Nagase 2025) | Direction/magnitude |

## Transparency Statement

This project uses **NO individual patient-level data**. All analyses are built from published aggregate endpoints and published popPK summary statistics. Virtual patients are simulated to illustrate E-R concepts and reproduce published dose-level summaries.

Any within-dose stratification (e.g., exposure quartile plots) is **model-implied from the virtual population** and labeled as such — not derived from real patient-level quartile data.

## Reproducibility

All outputs are fully reproducible from the published data tables in `data/published/` and the analysis scripts in `analysis/`. Random seeds are fixed (SEED_VPOP=20240501, SEED_ER=20240701).
