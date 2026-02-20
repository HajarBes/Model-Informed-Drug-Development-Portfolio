# =============================================================================
# 01_fit_bayesian.R — Bayesian PopPK Estimation (Stan/cmdstanr)
# 2-Compartment Oral Midazolam with Informative OSP Priors
# =============================================================================
# Mirrors the nlmixr2/NONMEM structural model:
#   - 3 ODE states: depot, central, peripheral
#   - BSV on Ka, CL, Vc only
#   - Allometric scaling on CL/Q (^0.75) and Vc/Vp (^1.0)
#   - Combined proportional + additive residual error
#
# Requires: cmdstanr, posterior, bayesplot, dplyr, readr
# =============================================================================

# --- Setup -------------------------------------------------------------------
script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(script_dir) || !file.exists(file.path(script_dir, "..", "analysis", "00_setup.R"))) {
  script_dir <- file.path(getwd(), "bayesian")
}
source(file.path(script_dir, "..", "analysis", "00_setup.R"))

# Bayesian-specific packages
for (pkg in c("cmdstanr", "posterior", "bayesplot")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(pkg, " is required. Install with: install.packages('", pkg, "')")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

DIR_BAYESIAN <- file.path(PROJECT_ROOT, "bayesian")
dir.create(DIR_BAYESIAN, showWarnings = FALSE, recursive = TRUE)

log_msg("=== Starting Bayesian PopPK Estimation ===")

# --- Load and prepare data ---------------------------------------------------
dat_raw <- read_csv(file.path(DIR_SIM, "analysis_dataset.csv"), show_col_types = FALSE)

# Observations: EVID=0, MDV=0, DV > 0
dat_obs <- dat_raw %>%

  filter(EVID == 0, MDV == 0, DV > 0) %>%
  select(ID, TIME, DV, WT, SEX, AGE)

# Dosing records: extract dose and weight per subject
dat_dose <- dat_raw %>%
  filter(EVID == 1) %>%
  select(ID, TIME_DOSE = TIME, AMT, WT) %>%
  distinct(ID, .keep_all = TRUE)

log_msg("Observations: ", nrow(dat_obs), " from ", n_distinct(dat_obs$ID), " subjects")
log_msg("Dose records: ", nrow(dat_dose), " subjects")

# --- Dose time alignment -----------------------------------------------------
# Subtract each subject's dosing time from observation times so dose is at t=0
dat_obs <- dat_obs %>%
  left_join(dat_dose %>% select(ID, TIME_DOSE, AMT), by = "ID") %>%
  mutate(
    TIME_ALIGNED = TIME - TIME_DOSE,
    DOSE = AMT
  ) %>%
  filter(TIME_ALIGNED > 0) %>%   # exclude time <= 0 (pre-dose)
  arrange(ID, TIME_ALIGNED)

log_msg("After time alignment: ", nrow(dat_obs), " observations")

# --- Build Stan data list ----------------------------------------------------
subj_ids <- sort(unique(dat_obs$ID))
N_subj <- length(subj_ids)
N_obs <- nrow(dat_obs)

# Per-subject indices (1-based, for batched ODE solving)
subj_info <- dat_obs %>%
  mutate(row_idx = row_number()) %>%
  group_by(ID) %>%
  summarise(
    start_idx = min(row_idx),
    end_idx = max(row_idx),
    n_obs = n(),
    .groups = "drop"
  ) %>%
  left_join(dat_dose %>% select(ID, AMT, WT), by = "ID")

stan_data <- list(
  N_subj        = N_subj,
  N_obs         = N_obs,
  DV            = dat_obs$DV,
  TIME          = dat_obs$TIME_ALIGNED,
  DOSE          = subj_info$AMT,
  WT            = subj_info$WT,
  start_idx     = subj_info$start_idx,
  end_idx       = subj_info$end_idx,
  n_obs_per_subj = subj_info$n_obs
)

log_msg("Stan data: ", N_subj, " subjects, ", N_obs, " observations")
log_msg("Obs per subject range: ", min(subj_info$n_obs), "-", max(subj_info$n_obs))

# Save data mapping for diagnostics
bayesian_data_map <- list(
  stan_data = stan_data,
  subj_info = subj_info,
  dat_obs = dat_obs,
  time_alignment = "Observation times shifted by subtracting dose time (TIME_DOSE) per subject"
)
saveRDS(bayesian_data_map, file.path(DIR_BAYESIAN, "bayesian_data_map.rds"))
log_msg("Data mapping saved to bayesian_data_map.rds")

# --- Compile Stan model ------------------------------------------------------
stan_file <- file.path(DIR_BAYESIAN, "midazolam_2cmt_poppk.stan")
log_msg("Compiling Stan model: ", stan_file)
mod <- cmdstan_model(stan_file)
log_msg("Stan model compiled successfully")

# --- Initial values -----------------------------------------------------------
# Centered on nlmixr2 final estimates with small perturbation
nlmixr2_est <- list(
  lka = 0.8153, lcl = 3.9476, lv1 = 3.8417, lq = 2.8205, lv2 = 4.1356
)

init_fun <- function() {
  list(
    lka  = nlmixr2_est$lka + rnorm(1, 0, 0.05),
    lcl  = nlmixr2_est$lcl + rnorm(1, 0, 0.05),
    lv1  = nlmixr2_est$lv1 + rnorm(1, 0, 0.05),
    lq   = nlmixr2_est$lq  + rnorm(1, 0, 0.05),
    lv2  = nlmixr2_est$lv2 + rnorm(1, 0, 0.05),
    omega_ka = abs(0.60 + rnorm(1, 0, 0.05)),
    omega_cl = abs(0.30 + rnorm(1, 0, 0.05)),
    omega_v1 = abs(0.20 + rnorm(1, 0, 0.05)),
    z_ka = rnorm(N_subj, 0, 0.1),
    z_cl = rnorm(N_subj, 0, 0.1),
    z_v1 = rnorm(N_subj, 0, 0.1),
    sigma_prop = abs(0.25 + rnorm(1, 0, 0.02)),
    sigma_add  = abs(0.50 + rnorm(1, 0, 0.05))
  )
}

# --- MCMC sampling ------------------------------------------------------------
log_msg("Starting MCMC sampling: 4 chains, 1000 warmup + 1000 sampling")
log_msg("  adapt_delta = 0.90, max_treedepth = 12, seed = 20240315")

t_start <- Sys.time()

fit <- mod$sample(
  data            = stan_data,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  adapt_delta     = 0.90,
  max_treedepth   = 12,
  seed            = 20240315,
  init            = init_fun,
  refresh         = 100
)

t_elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
log_msg("MCMC completed in ", round(t_elapsed, 1), " minutes")

# --- Quick convergence summary ------------------------------------------------
pop_pars <- c("lka", "lcl", "lv1", "lq", "lv2",
              "omega_ka", "omega_cl", "omega_v1",
              "sigma_prop", "sigma_add")

summ <- fit$summary(variables = pop_pars)
log_msg("Convergence summary:")
log_msg("  Max Rhat: ", round(max(summ$rhat, na.rm = TRUE), 4))
log_msg("  Min ESS bulk: ", round(min(summ$ess_bulk, na.rm = TRUE), 0))
log_msg("  Min ESS tail: ", round(min(summ$ess_tail, na.rm = TRUE), 0))

# Divergent transitions
diag_summ <- fit$diagnostic_summary()
n_divergent <- sum(diag_summ$num_divergent)
log_msg("  Divergent transitions: ", n_divergent)
if (n_divergent > 0) {
  log_msg("  WARNING: ", n_divergent, " divergent transitions detected. ",
          "Consider increasing adapt_delta or reparameterizing.")
}

# Max treedepth
n_max_td <- sum(diag_summ$num_max_treedepth)
if (n_max_td > 0) {
  log_msg("  Max treedepth hits: ", n_max_td)
}

# --- Save outputs -------------------------------------------------------------

# Save CmdStanMCMC fit object
fit$save_object(file.path(DIR_BAYESIAN, "fit_bayesian.rds"))
log_msg("Fit object saved to fit_bayesian.rds")

# Summary table
summ_full <- fit$summary(variables = c(pop_pars,
                                        "Ka_pop", "CL_pop", "V1_pop", "Q_pop", "V2_pop"))
write_csv(summ_full, file.path(DIR_BAYESIAN, "bayesian_summary.csv"))
log_msg("Summary table saved to bayesian_summary.csv")

# Print key estimates
log_msg("Population parameter estimates (median [95% CrI]):")
for (p in c("Ka_pop", "CL_pop", "V1_pop", "Q_pop", "V2_pop")) {
  row <- summ_full %>% filter(variable == p)
  log_msg(sprintf("  %s: %.2f [%.2f, %.2f]", p, row$median, row$q5, row$q95))
}

log_msg("BSV SDs (median [95% CrI]):")
for (p in c("omega_ka", "omega_cl", "omega_v1")) {
  row <- summ_full %>% filter(variable == p)
  log_msg(sprintf("  %s: %.3f [%.3f, %.3f]", p, row$median, row$q5, row$q95))
}

log_msg("Residual error SDs (median [95% CrI]):")
for (p in c("sigma_prop", "sigma_add")) {
  row <- summ_full %>% filter(variable == p)
  log_msg(sprintf("  %s: %.3f [%.3f, %.3f]", p, row$median, row$q5, row$q95))
}

log_msg("=== Bayesian PopPK Estimation Complete ===")
log_msg("Runtime: ", round(t_elapsed, 1), " minutes")
log_msg("Next: Run 02_bayesian_diagnostics.R for convergence and posterior checks")
