// =============================================================================
// midazolam_2cmt_poppk.stan
// 2-Compartment Oral PopPK Model for Midazolam (CYP3A4 Probe)
//
// Mirrors the nlmixr2/NONMEM structural model:
//   - 2-compartment with first-order absorption (analytical solution)
//   - BSV on Ka, CL, Vc only (no BSV on Q or Vp)
//   - Allometric scaling: CL/Q ~ (WT/70)^0.75, Vc/Vp ~ (WT/70)^1.0
//   - Combined proportional + additive residual error
//   - Non-centered parameterization for efficient MCMC sampling
//
// Uses closed-form analytical solution (tri-exponential) instead of ODE
// solver for performance. Equivalent to NONMEM ADVAN4/TRANS4.
//
// Priors centered on OSP Midazolam PBPK (Hanke et al. CPT:PSP 2018)
// =============================================================================

functions {
  // Analytical solution for 2-compartment oral model
  // Returns central compartment concentration at time t
  // Equivalent to depot -> central -> peripheral with first-order absorption
  real two_cmt_oral_conc(real t, real dose, real Ka,
                         real ke, real k12, real k21, real V1) {
    // Hybrid rate constants (eigenvalues of disposition matrix)
    real sum_k = ke + k12 + k21;
    real disc = sqrt(fmax(square(sum_k) - 4.0 * ke * k21, 1e-20));
    real alpha = 0.5 * (sum_k + disc);
    real beta  = 0.5 * (sum_k - disc);

    // Coefficients for tri-exponential solution
    // C(t) = dose*Ka/V1 * [A*exp(-alpha*t) + B*exp(-beta*t) + C*exp(-Ka*t)]
    real A_coef = (k21 - alpha) / ((beta - alpha) * (Ka - alpha));
    real B_coef = (k21 - beta)  / ((alpha - beta) * (Ka - beta));
    real C_coef = (k21 - Ka)    / ((alpha - Ka) * (beta - Ka));

    real conc = (dose * Ka / V1) *
                (A_coef * exp(-alpha * t) +
                 B_coef * exp(-beta * t) +
                 C_coef * exp(-Ka * t));

    return fmax(conc, 1e-12);  // IPRED floor
  }
}

data {
  int<lower=1> N_subj;                        // number of subjects
  int<lower=1> N_obs;                         // total observations
  array[N_obs] real<lower=0> DV;              // observed concentrations (ng/mL)
  array[N_obs] real<lower=0> TIME;            // observation times (h, aligned to dose at t=0)
  array[N_subj] real<lower=0> DOSE;           // dose per subject (ug)
  array[N_subj] real<lower=0> WT;             // body weight (kg)
  array[N_subj] int<lower=1> start_idx;       // first obs index per subject
  array[N_subj] int<lower=1> end_idx;         // last obs index per subject
  array[N_subj] int<lower=1> n_obs_per_subj;  // number of obs per subject
}

transformed data {
  real WT_REF = 70.0;  // reference body weight (kg)
}

parameters {
  // Population fixed effects (log-scale)
  real lka;   // log(Ka)   — absorption rate constant (1/h)
  real lcl;   // log(CL/F) — apparent clearance (L/h) at 70 kg
  real lv1;   // log(Vc/F) — apparent central volume (L) at 70 kg
  real lq;    // log(Q/F)  — apparent intercompartmental clearance (L/h) at 70 kg
  real lv2;   // log(Vp/F) — apparent peripheral volume (L) at 70 kg

  // BSV standard deviations (half-normal priors enforce positivity)
  real<lower=0> omega_ka;   // SD of eta on Ka
  real<lower=0> omega_cl;   // SD of eta on CL
  real<lower=0> omega_v1;   // SD of eta on Vc

  // Individual z-scores (non-centered parameterization)
  array[N_subj] real z_ka;
  array[N_subj] real z_cl;
  array[N_subj] real z_v1;

  // Residual error SDs
  real<lower=0> sigma_prop;  // proportional error SD
  real<lower=0> sigma_add;   // additive error SD (ng/mL)
}

transformed parameters {
  // Individual PK parameters
  array[N_subj] real<lower=0> Ka_ind;
  array[N_subj] real<lower=0> CL_ind;
  array[N_subj] real<lower=0> V1_ind;
  array[N_subj] real<lower=0> Q_ind;
  array[N_subj] real<lower=0> V2_ind;

  // Individual predicted concentrations
  array[N_obs] real<lower=0> IPRED;

  for (i in 1:N_subj) {
    // Non-centered parameterization with allometric scaling
    real wt_ratio = WT[i] / WT_REF;
    Ka_ind[i] = exp(lka + omega_ka * z_ka[i]);
    CL_ind[i] = exp(lcl + omega_cl * z_cl[i]) * pow(wt_ratio, 0.75);
    V1_ind[i] = exp(lv1 + omega_v1 * z_v1[i]) * pow(wt_ratio, 1.0);
    Q_ind[i]  = exp(lq) * pow(wt_ratio, 0.75);   // no BSV
    V2_ind[i] = exp(lv2) * pow(wt_ratio, 1.0);    // no BSV

    // Micro-rate constants
    real ke  = CL_ind[i] / V1_ind[i];
    real k12 = Q_ind[i]  / V1_ind[i];
    real k21 = Q_ind[i]  / V2_ind[i];

    // Analytical solution for each observation time
    for (j in 1:n_obs_per_subj[i]) {
      int idx = start_idx[i] + j - 1;
      IPRED[idx] = two_cmt_oral_conc(TIME[idx], DOSE[i], Ka_ind[i],
                                      ke, k12, k21, V1_ind[i]);
    }
  }
}

model {
  // --- Priors (centered on OSP Midazolam PBPK, Hanke et al. 2018) ---

  // Population fixed effects (log-scale)
  // lka: log(2.5) = 0.916, SD=0.5 -> 95% prior range: ~0.9-6.8 /h
  lka ~ normal(log(2.5), 0.5);
  // lcl: log(50) = 3.912, SD=0.3 -> 95% prior range: ~27-91 L/h
  lcl ~ normal(log(50), 0.3);
  // lv1: log(45) = 3.807, SD=0.3 -> 95% prior range: ~25-82 L
  lv1 ~ normal(log(45), 0.3);
  // lq: log(15) = 2.708, SD=0.5 -> 95% prior range: ~5.5-41 L/h
  lq  ~ normal(log(15), 0.5);
  // lv2: log(55) = 4.007, SD=0.5 -> 95% prior range: ~20-150 L
  lv2 ~ normal(log(55), 0.5);

  // BSV SDs — weakly informative half-normal (lower bound enforced by <lower=0>)
  omega_ka ~ normal(0, 1.0);   // wide: Ka variability less precisely known
  omega_cl ~ normal(0, 0.5);   // TRUE_OMEGA$cl = 0.09 -> SD = 0.30
  omega_v1 ~ normal(0, 0.5);   // TRUE_OMEGA$vc = 0.04 -> SD = 0.20

  // Residual error SDs — weakly informative half-normal
  sigma_prop ~ normal(0, 0.5);  // TRUE_SIGMA$prop = 0.20
  sigma_add  ~ normal(0, 2.0);  // TRUE_SIGMA$add = 0.5 ng/mL

  // Individual z-scores (non-centered)
  z_ka ~ std_normal();
  z_cl ~ std_normal();
  z_v1 ~ std_normal();

  // --- Likelihood (combined proportional + additive residual error) ---
  for (j in 1:N_obs) {
    real sd_j = sqrt(square(sigma_prop * IPRED[j]) + square(sigma_add));
    DV[j] ~ normal(IPRED[j], sd_j);
  }
}

generated quantities {
  // Population-level parameters (natural scale, at reference weight)
  real Ka_pop = exp(lka);
  real CL_pop = exp(lcl);
  real V1_pop = exp(lv1);
  real Q_pop  = exp(lq);
  real V2_pop = exp(lv2);

  // Posterior predictive replications (for PPC)
  array[N_obs] real dv_rep;
  for (j in 1:N_obs) {
    real sd_j = sqrt(square(sigma_prop * IPRED[j]) + square(sigma_add));
    dv_rep[j] = normal_rng(IPRED[j], sd_j);
  }

  // Pointwise log-likelihood (for LOO-CV)
  array[N_obs] real log_lik;
  for (j in 1:N_obs) {
    real sd_j = sqrt(square(sigma_prop * IPRED[j]) + square(sigma_add));
    log_lik[j] = normal_lpdf(DV[j] | IPRED[j], sd_j);
  }
}
