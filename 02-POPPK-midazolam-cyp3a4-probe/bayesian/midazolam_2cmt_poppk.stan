// =============================================================================
// midazolam_2cmt_poppk.stan
// 2-Compartment Oral PopPK Model for Midazolam (CYP3A4 Probe)
//
// Mirrors the nlmixr2/NONMEM structural model:
//   - 3 ODE states: depot, central, peripheral
//   - BSV on Ka, CL, Vc only (no BSV on Q or Vp)
//   - Allometric scaling: CL/Q ~ (WT/70)^0.75, Vc/Vp ~ (WT/70)^1.0
//   - Combined proportional + additive residual error
//   - Non-centered parameterization for efficient MCMC sampling
//
// Priors centered on OSP Midazolam PBPK (Hanke et al. CPT:PSP 2018)
// =============================================================================

functions {
  // 2-compartment oral ODE system
  // y[1] = depot, y[2] = central, y[3] = peripheral
  // theta[1] = Ka, theta[2] = ke, theta[3] = k12, theta[4] = k21
  vector ode_2cmt(real t, vector y, array[] real theta,
                  array[] real x_r, array[] int x_i) {
    vector[3] dydt;
    real Ka  = theta[1];
    real ke  = theta[2];
    real k12 = theta[3];
    real k21 = theta[4];

    dydt[1] = -Ka * y[1];                                        // depot
    dydt[2] =  Ka * y[1] - ke * y[2] - k12 * y[2] + k21 * y[3]; // central
    dydt[3] =  k12 * y[2] - k21 * y[3];                          // peripheral

    return dydt;
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
  array[0] real x_r;   // no real auxiliary data for ODE
  array[0] int x_i;    // no integer auxiliary data for ODE
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
    array[4] real theta = {Ka_ind[i], ke, k12, k21};

    // Initial condition: full dose in depot
    vector[3] y0 = [DOSE[i], 0.0, 0.0]';

    // Collect observation times for this subject
    int n_i = n_obs_per_subj[i];
    array[n_i] real times_i;
    for (j in 1:n_i) {
      times_i[j] = TIME[start_idx[i] + j - 1];
    }

    // Solve ODE (batched per subject)
    array[n_i] vector[3] sol = ode_rk45(ode_2cmt, y0, 0.0, times_i,
                                         theta, x_r, x_i);

    // Extract concentrations with IPRED floor
    for (j in 1:n_i) {
      int idx = start_idx[i] + j - 1;
      IPRED[idx] = fmax(sol[j][2] / V1_ind[i], 1e-12);
    }
  }
}

model {
  // --- Priors (centered on OSP Midazolam PBPK, Hanke et al. 2018) ---

  // Population fixed effects (log-scale)
  // lka: log(2.5) = 0.916, SD=0.5 → 95% prior range: ~0.9–6.8 /h
  lka ~ normal(log(2.5), 0.5);
  // lcl: log(50) = 3.912, SD=0.3 → 95% prior range: ~27–91 L/h
  lcl ~ normal(log(50), 0.3);
  // lv1: log(45) = 3.807, SD=0.3 → 95% prior range: ~25–82 L
  lv1 ~ normal(log(45), 0.3);
  // lq: log(15) = 2.708, SD=0.5 → 95% prior range: ~5.5–41 L/h
  lq  ~ normal(log(15), 0.5);
  // lv2: log(55) = 4.007, SD=0.5 → 95% prior range: ~20–150 L
  lv2 ~ normal(log(55), 0.5);

  // BSV SDs — weakly informative half-normal (lower bound enforced by <lower=0>)
  omega_ka ~ normal(0, 1.0);   // wide: Ka variability less precisely known
  omega_cl ~ normal(0, 0.5);   // TRUE_OMEGA$cl = 0.09 → SD = 0.30
  omega_v1 ~ normal(0, 0.5);   // TRUE_OMEGA$vc = 0.04 → SD = 0.20

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
