# =============================================================================
# 02_bayesian_diagnostics.R — Bayesian Diagnostics + Prior Sensitivity
# =============================================================================
# Generates 6 diagnostic figures + prior sensitivity comparison:
#   1. Trace plots (10 population parameters, 4 chains)
#   2. Rhat + ESS panel (convergence diagnostics)
#   3. Posterior densities (natural-scale with OSP true values)
#   4. Prior vs posterior overlay (data informativeness)
#   5. Posterior predictive check (PPC: density + time-course)
#   6. BSV posterior (omega SD histograms)
#   7. Prior sensitivity: informative vs weakly informative
#
# Requires: cmdstanr, posterior, bayesplot, ggplot2, dplyr, patchwork
# =============================================================================

# --- Setup -------------------------------------------------------------------
script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(script_dir) || !file.exists(file.path(script_dir, "..", "analysis", "00_setup.R"))) {
  script_dir <- file.path(getwd(), "bayesian")
}
source(file.path(script_dir, "..", "analysis", "00_setup.R"))

for (pkg in c("cmdstanr", "posterior", "bayesplot")) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

DIR_BAYESIAN <- file.path(PROJECT_ROOT, "bayesian")

log_msg("=== Starting Bayesian Diagnostics ===")

# --- Load fit and data -------------------------------------------------------
fit <- readRDS(file.path(DIR_BAYESIAN, "fit_bayesian.rds"))
data_map <- readRDS(file.path(DIR_BAYESIAN, "bayesian_data_map.rds"))
dat_obs <- data_map$dat_obs

pop_pars <- c("lka", "lcl", "lv1", "lq", "lv2",
              "omega_ka", "omega_cl", "omega_v1",
              "sigma_prop", "sigma_add")

nat_pars <- c("Ka_pop", "CL_pop", "V1_pop", "Q_pop", "V2_pop")

draws <- fit$draws(format = "df")

# =============================================================================
# FIGURE 1: Trace Plots
# =============================================================================
log_msg("Generating trace plots...")

p_trace <- mcmc_trace(fit$draws(), pars = pop_pars,
                      facet_args = list(ncol = 2, strip.position = "left")) +
  ggtitle("MCMC Trace Plots — Population Parameters (4 Chains)") +
  theme_poppk

ggsave(file.path(DIR_BAYESIAN, "bayes_trace_plots.png"), p_trace,
       width = 12, height = 14, dpi = 150)
log_msg("  Saved bayes_trace_plots.png")

# =============================================================================
# FIGURE 2: Rhat + ESS Panel
# =============================================================================
log_msg("Generating Rhat + ESS panel...")

summ <- fit$summary(variables = pop_pars)

p_rhat <- ggplot(summ, aes(x = reorder(variable, rhat), y = rhat)) +
  geom_point(size = 3, color = COL_PRIMARY) +
  geom_hline(yintercept = 1.01, linetype = "dashed", color = COL_SECONDARY) +
  coord_flip() +
  labs(x = NULL, y = expression(hat(R)),
       title = "Convergence: Rhat (target < 1.01)") +
  theme_poppk

p_ess <- ggplot(summ, aes(x = reorder(variable, ess_bulk), y = ess_bulk)) +
  geom_point(size = 3, color = COL_PRIMARY) +
  geom_hline(yintercept = 400, linetype = "dashed", color = COL_SECONDARY) +
  coord_flip() +
  labs(x = NULL, y = "ESS (bulk)",
       title = "Effective Sample Size (target > 400)") +
  theme_poppk

p_convergence <- p_rhat + p_ess +
  plot_annotation(title = "MCMC Convergence Diagnostics",
                  theme = theme(plot.title = element_text(face = "bold", size = 14)))

ggsave(file.path(DIR_BAYESIAN, "bayes_rhat_ess.png"), p_convergence,
       width = 12, height = 6, dpi = 150)
log_msg("  Saved bayes_rhat_ess.png")

# =============================================================================
# FIGURE 3: Posterior Densities (Natural Scale) with OSP True Values
# =============================================================================
log_msg("Generating posterior density plots...")

true_vals <- data.frame(
  variable = c("Ka_pop", "CL_pop", "V1_pop", "Q_pop", "V2_pop"),
  true_val = c(TRUE_PARAMS$ka, TRUE_PARAMS$cl_f, TRUE_PARAMS$vc_f,
               TRUE_PARAMS$q_f, TRUE_PARAMS$vp_f),
  label = c("Ka (1/h)", "CL/F (L/h)", "Vc/F (L)", "Q/F (L/h)", "Vp/F (L)")
)

nat_draws <- draws %>%
  select(all_of(nat_pars)) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  left_join(true_vals, by = "variable")

p_posterior <- ggplot(nat_draws, aes(x = value)) +
  geom_density(fill = COL_PRIMARY, alpha = 0.4, color = COL_PRIMARY) +
  geom_vline(aes(xintercept = true_val), linetype = "dashed",
             color = COL_SECONDARY, linewidth = 0.8) +
  facet_wrap(~ label, scales = "free", ncol = 3) +
  labs(x = "Parameter Value", y = "Density",
       title = "Posterior Distributions (Natural Scale)",
       subtitle = "Dashed red = OSP true value (Hanke et al. 2018)") +
  theme_poppk

ggsave(file.path(DIR_BAYESIAN, "bayes_posterior_densities.png"), p_posterior,
       width = 12, height = 8, dpi = 150)
log_msg("  Saved bayes_posterior_densities.png")

# =============================================================================
# FIGURE 4: Prior vs Posterior Overlay
# =============================================================================
log_msg("Generating prior vs posterior overlay...")

prior_spec <- data.frame(
  variable = c("lka", "lcl", "lv1", "lq", "lv2"),
  prior_mean = c(log(2.5), log(50), log(45), log(15), log(55)),
  prior_sd = c(0.5, 0.3, 0.3, 0.5, 0.5),
  label = c("log(Ka)", "log(CL/F)", "log(Vc/F)", "log(Q/F)", "log(Vp/F)")
)

log_draws <- draws %>%
  select(all_of(c("lka", "lcl", "lv1", "lq", "lv2"))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  left_join(prior_spec, by = "variable")

# Generate prior samples for overlay
set.seed(42)
prior_samples <- do.call(rbind, lapply(1:nrow(prior_spec), function(i) {
  data.frame(
    variable = prior_spec$variable[i],
    value = rnorm(4000, prior_spec$prior_mean[i], prior_spec$prior_sd[i]),
    label = prior_spec$label[i],
    source = "Prior"
  )
}))

log_draws$source <- "Posterior"
prior_samples$prior_mean <- NA
prior_samples$prior_sd <- NA

combined <- bind_rows(
  log_draws %>% select(variable, value, label, source),
  prior_samples %>% select(variable, value, label, source)
)

p_prior_post <- ggplot(combined, aes(x = value, fill = source, color = source)) +
  geom_density(alpha = 0.3) +
  facet_wrap(~ label, scales = "free", ncol = 3) +
  scale_fill_manual(values = c("Prior" = "grey60", "Posterior" = COL_PRIMARY)) +
  scale_color_manual(values = c("Prior" = "grey40", "Posterior" = COL_PRIMARY)) +
  labs(x = "Parameter Value (log scale)", y = "Density",
       title = "Prior vs Posterior Distributions",
       subtitle = "Grey = informative prior (OSP), Blue = posterior after MCMC",
       fill = NULL, color = NULL) +
  theme_poppk

ggsave(file.path(DIR_BAYESIAN, "bayes_prior_vs_posterior.png"), p_prior_post,
       width = 12, height = 8, dpi = 150)
log_msg("  Saved bayes_prior_vs_posterior.png")

# =============================================================================
# FIGURE 5: Posterior Predictive Check (PPC)
# =============================================================================
log_msg("Generating posterior predictive check...")

# 5a: Density overlay
dv_rep <- fit$draws("dv_rep", format = "matrix")
n_rep <- min(200, nrow(dv_rep))
rep_idx <- sample(1:nrow(dv_rep), n_rep)
dv_rep_sub <- dv_rep[rep_idx, ]

p_ppc_dens <- ppc_dens_overlay(y = dat_obs$DV, yrep = dv_rep_sub, alpha = 0.2) +
  labs(title = "PPC: Density Overlay (200 replications)",
       subtitle = "Dark line = observed; light lines = posterior predictive") +
  theme_poppk +
  coord_cartesian(xlim = c(0, quantile(dat_obs$DV, 0.99) * 1.5))

# 5b: Time-course PPC with 90% prediction interval
dv_rep_all <- fit$draws("dv_rep", format = "matrix")

ppc_summary <- data.frame(
  TIME = dat_obs$TIME_ALIGNED,
  DV = dat_obs$DV,
  median_rep = apply(dv_rep_all, 2, median),
  lo = apply(dv_rep_all, 2, quantile, probs = 0.05),
  hi = apply(dv_rep_all, 2, quantile, probs = 0.95)
)

p_ppc_time <- ggplot(ppc_summary, aes(x = TIME)) +
  geom_ribbon(aes(ymin = pmax(lo, 0), ymax = hi), fill = COL_PRIMARY, alpha = 0.2) +
  geom_line(aes(y = median_rep), color = COL_PRIMARY, linewidth = 0.6) +
  geom_point(aes(y = DV), color = COL_OBSERVED, size = 0.8, alpha = 0.5) +
  labs(x = "Time (h)", y = "Concentration (ng/mL)",
       title = "PPC: Time-Course with 90% Prediction Interval",
       subtitle = "Points = observed; ribbon = 90% PI; line = median prediction") +
  theme_poppk

p_ppc <- p_ppc_dens / p_ppc_time +
  plot_annotation(title = "Posterior Predictive Checks",
                  theme = theme(plot.title = element_text(face = "bold", size = 14)))

ggsave(file.path(DIR_BAYESIAN, "bayes_ppc.png"), p_ppc,
       width = 10, height = 10, dpi = 150)
log_msg("  Saved bayes_ppc.png")

# =============================================================================
# FIGURE 6: BSV Posterior (Omega SD Histograms)
# =============================================================================
log_msg("Generating BSV posterior plots...")

omega_draws <- draws %>%
  select(omega_ka, omega_cl, omega_v1) %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
  mutate(label = case_when(
    parameter == "omega_ka" ~ "omega[Ka] (SD)",
    parameter == "omega_cl" ~ "omega[CL] (SD)",
    parameter == "omega_v1" ~ "omega[Vc] (SD)"
  ))

# True BSV SDs (sqrt of TRUE_OMEGA variances)
true_omega_sd <- data.frame(
  parameter = c("omega_ka", "omega_cl", "omega_v1"),
  true_sd = c(sqrt(TRUE_OMEGA$ka), sqrt(TRUE_OMEGA$cl), sqrt(TRUE_OMEGA$vc)),
  label = c("omega[Ka] (SD)", "omega[CL] (SD)", "omega[Vc] (SD)")
)

p_bsv <- ggplot(omega_draws, aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40,
                 fill = COL_PRIMARY, alpha = 0.4, color = "white") +
  geom_density(color = COL_PRIMARY, linewidth = 0.8) +
  geom_vline(data = true_omega_sd, aes(xintercept = true_sd),
             linetype = "dashed", color = COL_SECONDARY, linewidth = 0.8) +
  facet_wrap(~ label, scales = "free_x", ncol = 3) +
  labs(x = "Standard Deviation", y = "Density",
       title = "BSV Posterior Distributions (SD Scale)",
       subtitle = "Dashed red = true simulation value (sqrt of TRUE_OMEGA variance))") +
  theme_poppk

ggsave(file.path(DIR_BAYESIAN, "bayes_bsv_posterior.png"), p_bsv,
       width = 12, height = 5, dpi = 150)
log_msg("  Saved bayes_bsv_posterior.png")

# =============================================================================
# PRIOR SENSITIVITY MINI-EXPERIMENT
# =============================================================================
log_msg("Starting prior sensitivity analysis...")
log_msg("  Refitting with weakly informative priors (SDs doubled)")

# Read and modify the Stan model for weakly informative priors
stan_code <- readLines(file.path(DIR_BAYESIAN, "midazolam_2cmt_poppk.stan"))

# Replace prior SDs: 0.3 -> 0.6, 0.5 -> 1.0 for fixed effects
# Replace omega priors: 0.5 -> 1.0, 1.0 -> 2.0
# Replace sigma priors: 0.5 -> 1.0, 2.0 -> 4.0
weak_code <- stan_code
weak_code <- gsub(
  "lka ~ normal(log(2.5), 0.5);",
  "lka ~ normal(log(2.5), 1.0);",
  weak_code, fixed = TRUE
)
weak_code <- gsub(
  "lcl ~ normal(log(50), 0.3);",
  "lcl ~ normal(log(50), 0.6);",
  weak_code, fixed = TRUE
)
weak_code <- gsub(
  "lv1 ~ normal(log(45), 0.3);",
  "lv1 ~ normal(log(45), 0.6);",
  weak_code, fixed = TRUE
)
weak_code <- gsub(
  "lq  ~ normal(log(15), 0.5);",
  "lq  ~ normal(log(15), 1.0);",
  weak_code, fixed = TRUE
)
weak_code <- gsub(
  "lv2 ~ normal(log(55), 0.5);",
  "lv2 ~ normal(log(55), 1.0);",
  weak_code, fixed = TRUE
)
weak_code <- gsub(
  "omega_ka ~ normal(0, 1.0);",
  "omega_ka ~ normal(0, 2.0);",
  weak_code, fixed = TRUE
)
weak_code <- gsub(
  "omega_cl ~ normal(0, 0.5);",
  "omega_cl ~ normal(0, 1.0);",
  weak_code, fixed = TRUE
)
weak_code <- gsub(
  "omega_v1 ~ normal(0, 0.5);",
  "omega_v1 ~ normal(0, 1.0);",
  weak_code, fixed = TRUE
)
weak_code <- gsub(
  "sigma_prop ~ normal(0, 0.5);",
  "sigma_prop ~ normal(0, 1.0);",
  weak_code, fixed = TRUE
)
weak_code <- gsub(
  "sigma_add  ~ normal(0, 2.0);",
  "sigma_add  ~ normal(0, 4.0);",
  weak_code, fixed = TRUE
)

weak_stan_file <- file.path(DIR_BAYESIAN, "midazolam_2cmt_weak_priors.stan")
writeLines(weak_code, weak_stan_file)

mod_weak <- cmdstan_model(weak_stan_file)
log_msg("  Weakly informative Stan model compiled")

# Define init function for sensitivity refit
N_subj_sens <- data_map$stan_data$N_subj
init_fun_sens <- function() {
  list(
    lka  = 0.8153 + rnorm(1, 0, 0.05),
    lcl  = 3.9476 + rnorm(1, 0, 0.05),
    lv1  = 3.8417 + rnorm(1, 0, 0.05),
    lq   = 2.8205 + rnorm(1, 0, 0.05),
    lv2  = 4.1356 + rnorm(1, 0, 0.05),
    omega_ka = abs(0.60 + rnorm(1, 0, 0.05)),
    omega_cl = abs(0.30 + rnorm(1, 0, 0.05)),
    omega_v1 = abs(0.20 + rnorm(1, 0, 0.05)),
    z_ka = rnorm(N_subj_sens, 0, 0.1),
    z_cl = rnorm(N_subj_sens, 0, 0.1),
    z_v1 = rnorm(N_subj_sens, 0, 0.1),
    sigma_prop = abs(0.25 + rnorm(1, 0, 0.02)),
    sigma_add  = abs(0.50 + rnorm(1, 0, 0.05))
  )
}

# Refit with weakly informative priors
fit_weak <- mod_weak$sample(
  data            = data_map$stan_data,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 1500,
  iter_sampling   = 1500,
  adapt_delta     = 0.95,
  max_treedepth   = 12,
  seed            = 20240315,
  init            = init_fun_sens,
  refresh         = 200
)

log_msg("  Weakly informative fit completed")

# --- FIGURE 7: Prior Sensitivity Comparison ---
draws_inform <- fit$draws(format = "df") %>%
  select(all_of(c("lka", "lcl", "lv1", "lq", "lv2"))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  mutate(prior_type = "Informative")

draws_weak <- fit_weak$draws(format = "df") %>%
  select(all_of(c("lka", "lcl", "lv1", "lq", "lv2"))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  mutate(prior_type = "Weakly Informative")

sens_combined <- bind_rows(draws_inform, draws_weak) %>%
  mutate(label = case_when(
    variable == "lka" ~ "log(Ka)",
    variable == "lcl" ~ "log(CL/F)",
    variable == "lv1" ~ "log(Vc/F)",
    variable == "lq"  ~ "log(Q/F)",
    variable == "lv2" ~ "log(Vp/F)"
  ))

p_sensitivity <- ggplot(sens_combined, aes(x = value, fill = prior_type, color = prior_type)) +
  geom_density(alpha = 0.3) +
  facet_wrap(~ label, scales = "free", ncol = 3) +
  scale_fill_manual(values = c("Informative" = COL_PRIMARY, "Weakly Informative" = COL_SECONDARY)) +
  scale_color_manual(values = c("Informative" = COL_PRIMARY, "Weakly Informative" = COL_SECONDARY)) +
  labs(x = "Parameter Value (log scale)", y = "Density",
       title = "Prior Sensitivity Analysis",
       subtitle = "Blue = informative OSP priors; Red = weakly informative (SDs doubled)",
       fill = "Prior Type", color = "Prior Type") +
  theme_poppk

ggsave(file.path(DIR_BAYESIAN, "bayes_prior_sensitivity.png"), p_sensitivity,
       width = 12, height = 8, dpi = 150)
log_msg("  Saved bayes_prior_sensitivity.png")

# Quantify sensitivity
sens_summary <- sens_combined %>%
  group_by(variable, label, prior_type) %>%
  summarise(
    median = median(value),
    q025 = quantile(value, 0.025),
    q975 = quantile(value, 0.975),
    sd = sd(value),
    .groups = "drop"
  ) %>%
  arrange(variable, prior_type)

write_csv(sens_summary, file.path(DIR_BAYESIAN, "bayesian_diagnostics.csv"))
log_msg("  Prior sensitivity summary saved to bayesian_diagnostics.csv")

# Interpret sensitivity
for (v in unique(sens_summary$variable)) {
  inform <- sens_summary %>% filter(variable == v, prior_type == "Informative")
  weak <- sens_summary %>% filter(variable == v, prior_type == "Weakly Informative")
  sd_ratio <- weak$sd / inform$sd
  med_diff <- abs(weak$median - inform$median)
  log_msg(sprintf("  %s: SD ratio (weak/inform) = %.2f, median diff = %.4f",
                  v, sd_ratio, med_diff))
}

# Clean up temporary Stan file
unlink(weak_stan_file)

log_msg("=== Bayesian Diagnostics Complete ===")
log_msg("Figures saved to: ", DIR_BAYESIAN)
log_msg("Next: Run 03_compare_methods.R for three-way comparison")
