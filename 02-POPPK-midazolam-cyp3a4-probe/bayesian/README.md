# Bayesian PopPK Estimation (Stan/cmdstanr) with Informative Priors

This module re-estimates the 2-compartment oral midazolam model using
full Bayesian MCMC (Stan/cmdstanr) with informative priors derived from
the OSP Midazolam PBPK model. The structural model mirrors the
nlmixr2 and NONMEM implementations, creating a three-tool comparison:

> **nlmixr2 SAEM | NONMEM SAEM | Stan MCMC**

The same 2-compartment model estimated three ways, demonstrating fluency
in the emerging Bayesian pharmacometrics ecosystem.

---

## Files

| File | Description |
|------|-------------|
| `midazolam_2cmt_poppk.stan` | Stan model: 2-cmt oral PopPK with informative priors |
| `01_fit_bayesian.R` | Data preparation, Stan compilation, MCMC sampling |
| `02_bayesian_diagnostics.R` | Convergence diagnostics, PPC, prior sensitivity |
| `03_compare_methods.R` | Three-way comparison + LOO-CV + structural check |
| `README.md` | This file |

---

## How to Run

### Prerequisites

- **R** (>= 4.2) with packages: `cmdstanr`, `posterior`, `bayesplot`, `loo`, `tidyverse`, `patchwork`
- **CmdStan** (>= 2.33): Install via `cmdstanr::install_cmdstan()`

```r
# Install cmdstanr (if needed)
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)
install_cmdstan()   # downloads and builds CmdStan
```

### Run Sequence

```bash
cd 02-POPPK-midazolam-cyp3a4-probe

# Step 1: Fit Bayesian model (requires CmdStan)
Rscript bayesian/01_fit_bayesian.R

# Step 2: Diagnostics + prior sensitivity
Rscript bayesian/02_bayesian_diagnostics.R

# Step 3: Three-way comparison + LOO-CV
Rscript bayesian/03_compare_methods.R
```

**Note:** Step 1 requires a compiled Stan model and runs 4 MCMC chains.
Runtime depends on hardware (typically 10-30 minutes on modern CPUs).

---

## Prior Specification

All prior centers are derived from `analysis/00_setup.R` TRUE_PARAMS, which
encode the OSP Midazolam PBPK model (Hanke et al. CPT:PSP 2018, qualified
against 40+ clinical studies).

| Parameter | Prior | OSP Source | 95% Prior Range |
|-----------|-------|------------|-----------------|
| lka (log Ka) | Normal(log(2.5), 0.5) | TRUE_PARAMS\$ka = 2.5 /h | 0.9--6.8 /h |
| lcl (log CL/F) | Normal(log(50), 0.3) | TRUE_PARAMS\$cl_f = 50 L/h | 27--91 L/h |
| lv1 (log Vc/F) | Normal(log(45), 0.3) | TRUE_PARAMS\$vc_f = 45 L | 25--82 L |
| lq (log Q/F) | Normal(log(15), 0.5) | TRUE_PARAMS\$q_f = 15 L/h | 5.5--41 L/h |
| lv2 (log Vp/F) | Normal(log(55), 0.5) | TRUE_PARAMS\$vp_f = 55 L | 20--150 L |
| omega_ka | half-Normal(0, 1.0) | Weakly informative | -- |
| omega_cl | half-Normal(0, 0.5) | TRUE_OMEGA\$cl ~ CV 30% | -- |
| omega_v1 | half-Normal(0, 0.5) | TRUE_OMEGA\$vc ~ CV 20% | -- |
| sigma_prop | half-Normal(0, 0.5) | TRUE_SIGMA\$prop = 0.20 | -- |
| sigma_add | half-Normal(0, 2.0) | TRUE_SIGMA\$add = 0.5 ng/mL | -- |

**Prior width rationale:**
- CL/F and Vc/F: SD = 0.3 (narrower) -- OSP CL is well-validated across multiple studies
- Ka, Q/F, Vp/F: SD = 0.5 (wider) -- less precisely identified from oral data
- All omega/sigma priors: weakly informative half-normal -- prevents pathological shrinkage to zero

---

## Parameter Mapping

| Parameter | nlmixr2 | NONMEM | Stan | Scale Note |
|-----------|---------|--------|------|------------|
| Ka | `lka` / `eta.ka` | THETA(1) / ETA(1) | `lka` / `z_ka` | All log-scale fixed effects |
| CL/F | `lcl` / `eta.cl` | THETA(2) / ETA(2) | `lcl` / `z_cl` | All log-scale fixed effects |
| Vc/F | `lv1` / `eta.v1` | THETA(3) / ETA(3) | `lv1` / `z_v1` | All log-scale fixed effects |
| Q/F | `lq` | THETA(4) | `lq` | No random effect |
| Vp/F | `lv2` | THETA(5) | `lv2` | No random effect |
| BSV Ka | omega^2 (variance) | OMEGA(1,1) (variance) | `omega_ka` (SD) | **nlmixr2/NONMEM = variance; Stan = SD** |
| BSV CL | omega^2 (variance) | OMEGA(2,2) (variance) | `omega_cl` (SD) | **nlmixr2/NONMEM = variance; Stan = SD** |
| BSV Vc | omega^2 (variance) | OMEGA(3,3) (variance) | `omega_v1` (SD) | **nlmixr2/NONMEM = variance; Stan = SD** |
| Prop err | `prop.err` (SD) | SIGMA(1,1) (variance) | `sigma_prop` (SD) | nlmixr2/Stan = SD; NONMEM = variance |
| Add err | `add.err` (SD) | SIGMA(2,2) (variance) | `sigma_add` (SD) | nlmixr2/Stan = SD; NONMEM = variance |

**Conversion:** SD = sqrt(variance). For example, TRUE_OMEGA\$cl = 0.09 (variance) -> omega_cl = 0.30 (SD).

### ETA Mapping

BSV (between-subject variability) is estimated on **Ka, CL, and Vc only** --
Q and Vp are fixed effects only (no random effects). This is consistent
across all three tools (nlmixr2, NONMEM, Stan).

---

## ODE System

Three-state system (depot, central, peripheral) for 2-compartment oral:

```
dA_depot/dt      = -Ka * A_depot
dA_central/dt    =  Ka * A_depot - ke * A_central - k12 * A_central + k21 * A_peripheral
dA_peripheral/dt =  k12 * A_central - k21 * A_peripheral
```

Where micro-constants: `ke = CL/V1`, `k12 = Q/V1`, `k21 = Q/V2`.

### Allometric Scaling

```
CL_i = exp(lcl + omega_cl * z_cl) * (WT_i / 70)^0.75
V1_i = exp(lv1 + omega_v1 * z_v1) * (WT_i / 70)^1.0
Q_i  = exp(lq)  * (WT_i / 70)^0.75    (no BSV)
V2_i = exp(lv2) * (WT_i / 70)^1.0     (no BSV)
Ka_i = exp(lka + omega_ka * z_ka)      (no WT scaling on Ka)
```

### Residual Error

Combined proportional + additive:

```
SD_j = sqrt((sigma_prop * IPRED_j)^2 + sigma_add^2)
DV_j ~ Normal(IPRED_j, SD_j)
```

This matches nlmixr2's `cp ~ prop(prop.err) + add(add.err)` and NONMEM's
`Y = F*(1+EPS(1)) + EPS(2)`.

---

## Output Files

| File | Description |
|------|-------------|
| `fit_bayesian.rds` | CmdStanMCMC fit object (full posterior draws) |
| `bayesian_summary.csv` | Parameter summary: median, mean, sd, q5, q95, Rhat, ESS |
| `bayesian_data_map.rds` | Stan data list + subject info + time alignment |
| `bayesian_diagnostics.csv` | Prior sensitivity comparison summary |
| `three_way_comparison.csv` | Parameter estimates: OSP true / nlmixr2 / NONMEM / Bayesian |
| `bsv_comparison_bayesian.csv` | BSV comparison (variance vs SD with conversion) |
| `loo_summary.csv` | LOO-CV results: elpd_loo, p_loo, LOOIC, Pareto-k counts |
| `bayes_trace_plots.png` | MCMC trace plots (10 parameters, 4 chains) |
| `bayes_rhat_ess.png` | Rhat + ESS convergence panel |
| `bayes_posterior_densities.png` | Natural-scale posteriors with OSP true values |
| `bayes_prior_vs_posterior.png` | Prior vs posterior overlay (data informativeness) |
| `bayes_ppc.png` | Posterior predictive check: density + time-course 90% PI |
| `bayes_bsv_posterior.png` | BSV omega SD posterior histograms |
| `bayes_prior_sensitivity.png` | Informative vs weakly informative posteriors |
| `comparison_forest.png` | Forest-style three-way parameter comparison |
| `comparison_pk_profile.png` | Typical PK profile overlay (3 methods + OSP, semi-log) |
| `comparison_ratio.png` | Parameter ratio plot (estimate/true, 0.80-1.25 band) |

---

## LOO-CV Results

Leave-one-out cross-validation via Pareto-smoothed importance sampling
(PSIS-LOO; Vehtari et al. 2017) provides a measure of out-of-sample
predictive performance.

- **elpd_loo**: expected log pointwise predictive density
- **p_loo**: effective number of parameters
- **LOOIC**: LOO information criterion (-2 * elpd_loo)

Pareto-k diagnostics flag observations where the LOO approximation may
be unreliable (k > 0.7). Results are saved in `loo_summary.csv`.

---

## Prior Sensitivity Results

A mini-experiment compares posteriors under two prior specifications:

1. **Informative** (default): priors centered on OSP PBPK values with moderate SDs
2. **Weakly informative**: same centers, SDs doubled (0.3 -> 0.6, 0.5 -> 1.0)

If posteriors are similar under both specifications, the data dominates the
prior information. If posteriors differ, the prior is influential and should
be carefully justified. Results are saved in `bayes_prior_sensitivity.png`
and `bayesian_diagnostics.csv`.

---

## Technical Notes

### Non-Centered Parameterization

Individual parameters use non-centered parameterization (z-scores):
```
Ka_i = exp(lka + omega_ka * z_ka_i),  z_ka_i ~ N(0,1)
```
This avoids the "funnel" geometry that causes divergent transitions in
centered parameterizations, especially when omega is small relative to
the number of subjects.

### Analytical Solution (No ODE Solver)

The 2-compartment oral model uses a closed-form tri-exponential analytical
solution rather than a numerical ODE solver. This is equivalent to
NONMEM's ADVAN4/TRANS4 and dramatically faster for MCMC sampling
(~100x vs `ode_rk45`). The analytical solution is exact for linear
compartmental systems and avoids numerical integration failures during
exploration of extreme parameter regions in warmup.

### adapt_delta = 0.90

The target acceptance rate is set to 0.90 (default: 0.80) to reduce
divergent transitions. Higher values increase computation time but
improve sampling reliability for hierarchical models.

### IPRED Floor

`fmax(A_central / V1, 1e-12)` prevents zero or negative predicted
concentrations during MCMC exploration of extreme parameter regions,
which would produce NaN log-likelihood values.

### Torsten Note

For multi-dose studies or complex PK models, the
[Torsten Stan library](https://github.com/metrumresearchgroup/Torsten)
(Margossian & Gillespie 2018) provides specialized PK ODE solvers and
analytical solutions (PMXSolve). This single-dose model uses Stan's
built-in `ode_rk45`; for production multi-dose PopPK, Torsten's
event-driven solver would be preferred.

---

> **Bayesian Workflow Checklist**
>
> 1. Prior specification (informative + sensitivity analysis)
> 2. Non-centered hierarchical parameterization
> 3. Convergence diagnostics (Rhat < 1.01, ESS > 400, zero divergences)
> 4. Posterior predictive checks (distribution + time-course)
> 5. Model comparison (LOO-PSIS + Pareto-k diagnostics)

---

## References

- Carpenter B et al. Stan: A probabilistic programming language. J Stat Softw 2017;76(1)
- Margossian CC, Gillespie WR. Stan functions for pharmacometrics modeling. J Pharmacokinet Pharmacodyn 2018;45(Suppl 1):S81
- Gabry J et al. Visualization in Bayesian workflow. J R Stat Soc A 2019;182(2):389-402
- Hanke N et al. PBPK models for CYP3A4 and P-gp substrates. CPT Pharmacometrics Syst Pharmacol 2018;7(10):647-659
- Vehtari A, Gelman A, Gabry J. Practical Bayesian model evaluation using LOO-CV and WAIC. Stat Comput 2017;27(5):1413-1432
