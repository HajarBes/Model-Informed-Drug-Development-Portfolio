# =============================================================================
# 03_compare_methods.R — Three-Way Method Comparison + LOO-CV
# nlmixr2 SAEM | NONMEM SAEM | Stan MCMC
# =============================================================================
# Creates:
#   - Three-way parameter comparison table
#   - BSV comparison table (variance vs SD, with conversion note)
#   - LOO-CV with Pareto-k diagnostics
#   - Structural consistency check (typical profile overlay)
#   - Figure 1: Forest-style parameter comparison
#   - Figure 2: Typical PK profile overlay (3 methods + OSP true, semi-log)
#   - Figure 3: Parameter ratio plot (estimate/true)
#
# Requires: cmdstanr, posterior, loo, dplyr, ggplot2, patchwork, deSolve
# =============================================================================

# --- Setup -------------------------------------------------------------------
script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(script_dir) || !file.exists(file.path(script_dir, "..", "analysis", "00_setup.R"))) {
  script_dir <- file.path(getwd(), "bayesian")
}
source(file.path(script_dir, "..", "analysis", "00_setup.R"))

for (pkg in c("cmdstanr", "posterior", "loo")) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

DIR_BAYESIAN <- file.path(PROJECT_ROOT, "bayesian")

log_msg("=== Starting Three-Way Method Comparison ===")

# --- Load fits ---------------------------------------------------------------
fit_bayes <- readRDS(file.path(DIR_BAYESIAN, "fit_bayesian.rds"))
fit_nlmixr2 <- tryCatch(
  readRDS(file.path(DIR_TABLES, "fit_base.rds")),
  error = function(e) { log_msg("nlmixr2 fit not found"); NULL }
)

# --- Extract parameter estimates ---------------------------------------------

# Bayesian: median [5th, 95th] from posterior
bayes_summ <- fit_bayes$summary(
  variables = c("Ka_pop", "CL_pop", "V1_pop", "Q_pop", "V2_pop",
                "omega_ka", "omega_cl", "omega_v1",
                "sigma_prop", "sigma_add")
)

# nlmixr2 estimates (from parameter_estimates_base.csv)
nlmixr2_est <- read_csv(file.path(DIR_TABLES, "parameter_estimates_base.csv"),
                         show_col_types = FALSE)

# NONMEM estimates (from control stream THETAs, OMEGAs, SIGMAs)
nonmem_est <- list(
  lka = 0.8153, lcl = 3.9476, lv1 = 3.8417, lq = 2.8205, lv2 = 4.1356,
  Ka  = exp(0.8153), CL = exp(3.9476), V1 = exp(3.8417),
  Q   = exp(2.8205), V2 = exp(4.1356),
  omega_ka_var = 0.0895, omega_cl_var = 0.267, omega_v1_var = 0.0279,
  sigma_prop_var = 0.0617, sigma_add_var = 0.2437
)

# nlmixr2 natural-scale values
nlmixr2_natural <- list(
  Ka = exp(nlmixr2_est$Estimate[nlmixr2_est$Parameter == "lka"]),
  CL = exp(nlmixr2_est$Estimate[nlmixr2_est$Parameter == "lcl"]),
  V1 = exp(nlmixr2_est$Estimate[nlmixr2_est$Parameter == "lv1"]),
  Q  = exp(nlmixr2_est$Estimate[nlmixr2_est$Parameter == "lq"]),
  V2 = exp(nlmixr2_est$Estimate[nlmixr2_est$Parameter == "lv2"])
)

# --- Three-Way Comparison Table -----------------------------------------------
log_msg("Building three-way comparison table...")

get_bayes_str <- function(par) {
  row <- bayes_summ %>% filter(variable == par)
  sprintf("%.2f [%.2f, %.2f]", row$median, row$q5, row$q95)
}

comparison <- tibble(
  Parameter = c("Ka (1/h)", "CL/F (L/h)", "Vc/F (L)", "Q/F (L/h)", "Vp/F (L)"),
  OSP_True = c(TRUE_PARAMS$ka, TRUE_PARAMS$cl_f, TRUE_PARAMS$vc_f,
               TRUE_PARAMS$q_f, TRUE_PARAMS$vp_f),
  nlmixr2_SAEM = c(nlmixr2_natural$Ka, nlmixr2_natural$CL, nlmixr2_natural$V1,
                    nlmixr2_natural$Q, nlmixr2_natural$V2),
  NONMEM_SAEM = c(nonmem_est$Ka, nonmem_est$CL, nonmem_est$V1,
                  nonmem_est$Q, nonmem_est$V2),
  Bayesian_Median_95CrI = c(
    get_bayes_str("Ka_pop"), get_bayes_str("CL_pop"), get_bayes_str("V1_pop"),
    get_bayes_str("Q_pop"), get_bayes_str("V2_pop")
  )
)

write_csv(comparison, file.path(DIR_BAYESIAN, "three_way_comparison.csv"))
log_msg("  Saved three_way_comparison.csv")

# Print comparison
log_msg("Three-way parameter comparison:")
for (i in seq_len(nrow(comparison))) {
  r <- comparison[i, ]
  log_msg(sprintf("  %s: True=%.1f | nlmixr2=%.2f | NONMEM=%.2f | Bayes=%s",
                  r$Parameter, r$OSP_True, r$nlmixr2_SAEM, r$NONMEM_SAEM,
                  r$Bayesian_Median_95CrI))
}

# --- BSV Comparison Table -----------------------------------------------------
log_msg("Building BSV comparison table...")

bayes_omega <- bayes_summ %>%
  filter(variable %in% c("omega_ka", "omega_cl", "omega_v1"))

bsv_comparison <- tibble(
  Parameter = c("BSV Ka", "BSV CL", "BSV Vc"),
  OSP_True_Variance = c(TRUE_OMEGA$ka, TRUE_OMEGA$cl, TRUE_OMEGA$vc),
  OSP_True_SD = c(sqrt(TRUE_OMEGA$ka), sqrt(TRUE_OMEGA$cl), sqrt(TRUE_OMEGA$vc)),
  nlmixr2_Variance = c(
    # nlmixr2 stores omega as variance in the fit object
    if (!is.null(fit_nlmixr2)) fit_nlmixr2$omega["eta.ka", "eta.ka"] else NA,
    if (!is.null(fit_nlmixr2)) fit_nlmixr2$omega["eta.cl", "eta.cl"] else NA,
    if (!is.null(fit_nlmixr2)) fit_nlmixr2$omega["eta.v1", "eta.v1"] else NA
  ),
  NONMEM_Variance = c(nonmem_est$omega_ka_var, nonmem_est$omega_cl_var, nonmem_est$omega_v1_var),
  Bayesian_SD_Median = bayes_omega$median,
  Bayesian_SD_Q5 = bayes_omega$q5,
  Bayesian_SD_Q95 = bayes_omega$q95,
  Note = "nlmixr2/NONMEM report variance (omega^2); Stan estimates SD (omega). Convert: SD = sqrt(variance)."
)

write_csv(bsv_comparison, file.path(DIR_BAYESIAN, "bsv_comparison_bayesian.csv"))
log_msg("  Saved bsv_comparison_bayesian.csv")

# --- LOO-CV (Leave-One-Out Cross-Validation) ----------------------------------
log_msg("Computing LOO-CV...")

log_lik <- fit_bayes$draws("log_lik", format = "matrix")
loo_result <- loo::loo(log_lik)

log_msg("LOO-CV results:")
log_msg(sprintf("  elpd_loo = %.1f (SE = %.1f)", loo_result$estimates["elpd_loo", "Estimate"],
                loo_result$estimates["elpd_loo", "SE"]))
log_msg(sprintf("  p_loo    = %.1f (SE = %.1f)", loo_result$estimates["p_loo", "Estimate"],
                loo_result$estimates["p_loo", "SE"]))
log_msg(sprintf("  looic    = %.1f (SE = %.1f)", loo_result$estimates["looic", "Estimate"],
                loo_result$estimates["looic", "SE"]))

# Pareto-k diagnostics
pareto_k <- loo_result$diagnostics$pareto_k
k_summary <- data.frame(
  category = c("good (k < 0.5)", "ok (0.5 < k < 0.7)", "bad (k > 0.7)"),
  count = c(sum(pareto_k < 0.5), sum(pareto_k >= 0.5 & pareto_k < 0.7), sum(pareto_k >= 0.7)),
  pct = round(c(sum(pareto_k < 0.5), sum(pareto_k >= 0.5 & pareto_k < 0.7),
                sum(pareto_k >= 0.7)) / length(pareto_k) * 100, 1)
)

log_msg("Pareto-k diagnostics:")
for (i in seq_len(nrow(k_summary))) {
  log_msg(sprintf("  %s: %d (%.1f%%)", k_summary$category[i],
                  k_summary$count[i], k_summary$pct[i]))
}

if (any(pareto_k >= 0.7)) {
  log_msg("  CAUTION: ", sum(pareto_k >= 0.7),
          " observations with Pareto-k > 0.7 — LOO estimates may be unreliable for these points.")
}

loo_summary <- tibble(
  metric = c("elpd_loo", "p_loo", "looic",
             "n_pareto_k_good", "n_pareto_k_ok", "n_pareto_k_bad"),
  estimate = c(
    loo_result$estimates["elpd_loo", "Estimate"],
    loo_result$estimates["p_loo", "Estimate"],
    loo_result$estimates["looic", "Estimate"],
    k_summary$count[1], k_summary$count[2], k_summary$count[3]
  ),
  se = c(
    loo_result$estimates["elpd_loo", "SE"],
    loo_result$estimates["p_loo", "SE"],
    loo_result$estimates["looic", "SE"],
    NA, NA, NA
  )
)

write_csv(loo_summary, file.path(DIR_BAYESIAN, "loo_summary.csv"))
log_msg("  Saved loo_summary.csv")

# --- Structural Consistency Check ---------------------------------------------
log_msg("Running structural consistency check (typical profile overlay)...")

# Simulate typical PK profile from each method's estimates
dose_ug <- 7500  # 7.5 mg
times <- seq(0.1, 12, by = 0.1)

# OSP true
pk_true <- simulate_subject_pk(dose_ug, times,
                                ka = TRUE_PARAMS$ka, cl = TRUE_PARAMS$cl_f,
                                v1 = TRUE_PARAMS$vc_f, q = TRUE_PARAMS$q_f,
                                v2 = TRUE_PARAMS$vp_f)
pk_true$Method <- "OSP True"

# nlmixr2
pk_nlmixr2 <- simulate_subject_pk(dose_ug, times,
                                    ka = nlmixr2_natural$Ka, cl = nlmixr2_natural$CL,
                                    v1 = nlmixr2_natural$V1, q = nlmixr2_natural$Q,
                                    v2 = nlmixr2_natural$V2)
pk_nlmixr2$Method <- "nlmixr2 SAEM"

# NONMEM
pk_nonmem <- simulate_subject_pk(dose_ug, times,
                                  ka = nonmem_est$Ka, cl = nonmem_est$CL,
                                  v1 = nonmem_est$V1, q = nonmem_est$Q,
                                  v2 = nonmem_est$V2)
pk_nonmem$Method <- "NONMEM SAEM"

# Bayesian (posterior median)
bayes_med <- bayes_summ %>%
  select(variable, median) %>%
  deframe()

pk_bayes <- simulate_subject_pk(dose_ug, times,
                                 ka = bayes_med["Ka_pop"], cl = bayes_med["CL_pop"],
                                 v1 = bayes_med["V1_pop"], q = bayes_med["Q_pop"],
                                 v2 = bayes_med["V2_pop"])
pk_bayes$Method <- "Bayesian (Stan)"

pk_all <- bind_rows(pk_true, pk_nlmixr2, pk_nonmem, pk_bayes)
pk_all$Method <- factor(pk_all$Method,
                        levels = c("OSP True", "nlmixr2 SAEM", "NONMEM SAEM", "Bayesian (Stan)"))

# =============================================================================
# FIGURE 1: Forest-Style Parameter Comparison
# =============================================================================
log_msg("Generating forest-style parameter comparison figure...")

bayes_draws <- fit_bayes$draws(format = "df")

forest_data <- tibble(
  Parameter = rep(c("Ka (1/h)", "CL/F (L/h)", "Vc/F (L)", "Q/F (L/h)", "Vp/F (L)"), 3),
  Method = rep(c("nlmixr2 SAEM", "NONMEM SAEM", "Bayesian (Stan)"), each = 5),
  Estimate = c(
    nlmixr2_natural$Ka, nlmixr2_natural$CL, nlmixr2_natural$V1,
    nlmixr2_natural$Q, nlmixr2_natural$V2,
    nonmem_est$Ka, nonmem_est$CL, nonmem_est$V1, nonmem_est$Q, nonmem_est$V2,
    bayes_med["Ka_pop"], bayes_med["CL_pop"], bayes_med["V1_pop"],
    bayes_med["Q_pop"], bayes_med["V2_pop"]
  )
)

# Add Bayesian CrI
bayes_ci <- bayes_summ %>%
  filter(variable %in% c("Ka_pop", "CL_pop", "V1_pop", "Q_pop", "V2_pop")) %>%
  mutate(Parameter = c("Ka (1/h)", "CL/F (L/h)", "Vc/F (L)", "Q/F (L/h)", "Vp/F (L)"))

forest_data <- forest_data %>%
  left_join(
    bayes_ci %>%
      select(Parameter, q5, q95) %>%
      mutate(Method = "Bayesian (Stan)"),
    by = c("Parameter", "Method")
  )

true_ref <- tibble(
  Parameter = c("Ka (1/h)", "CL/F (L/h)", "Vc/F (L)", "Q/F (L/h)", "Vp/F (L)"),
  true_val = c(TRUE_PARAMS$ka, TRUE_PARAMS$cl_f, TRUE_PARAMS$vc_f,
               TRUE_PARAMS$q_f, TRUE_PARAMS$vp_f)
)

method_colors <- c("nlmixr2 SAEM" = COL_PRIMARY,
                    "NONMEM SAEM" = COL_SECONDARY,
                    "Bayesian (Stan)" = "#4DAF4A")

p_forest <- ggplot(forest_data, aes(x = Estimate, y = Method, color = Method)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = q5, xmax = q95), height = 0.2, na.rm = TRUE) +
  geom_vline(data = true_ref, aes(xintercept = true_val),
             linetype = "dashed", color = "grey40", linewidth = 0.6) +
  facet_wrap(~ Parameter, scales = "free_x", ncol = 3) +
  scale_color_manual(values = method_colors) +
  labs(x = "Parameter Value", y = NULL,
       title = "Three-Way Parameter Comparison",
       subtitle = "Dashed line = OSP true value | Error bars = Bayesian 90% CrI",
       color = NULL) +
  theme_poppk +
  theme(legend.position = "top")

ggsave(file.path(DIR_BAYESIAN, "comparison_forest.png"), p_forest,
       width = 12, height = 8, dpi = 150)
log_msg("  Saved comparison_forest.png")

# =============================================================================
# FIGURE 2: Typical PK Profile Overlay (Semi-Log)
# =============================================================================
log_msg("Generating typical PK profile overlay...")

method_colors_full <- c("OSP True" = "grey40",
                         "nlmixr2 SAEM" = COL_PRIMARY,
                         "NONMEM SAEM" = COL_SECONDARY,
                         "Bayesian (Stan)" = "#4DAF4A")
method_linetypes <- c("OSP True" = "dashed",
                       "nlmixr2 SAEM" = "solid",
                       "NONMEM SAEM" = "solid",
                       "Bayesian (Stan)" = "solid")

p_profile <- ggplot(pk_all, aes(x = time, y = CP, color = Method, linetype = Method)) +
  geom_line(linewidth = 0.9) +
  scale_y_log10(labels = scales::label_number()) +
  scale_color_manual(values = method_colors_full) +
  scale_linetype_manual(values = method_linetypes) +
  labs(x = "Time (h)", y = "Concentration (ng/mL, log scale)",
       title = "Typical PK Profile — Three-Way Comparison (7.5 mg, 70 kg)",
       subtitle = "Same 2-compartment model estimated by three methods",
       color = NULL, linetype = NULL) +
  annotation_logticks(sides = "l") +
  theme_poppk +
  theme(legend.position = "top")

ggsave(file.path(DIR_BAYESIAN, "comparison_pk_profile.png"), p_profile,
       width = 10, height = 6, dpi = 150)
log_msg("  Saved comparison_pk_profile.png")

# =============================================================================
# FIGURE 3: Parameter Ratio Plot (Estimate / True)
# =============================================================================
log_msg("Generating parameter ratio plot...")

ratio_data <- tibble(
  Parameter = rep(c("Ka", "CL/F", "Vc/F", "Q/F", "Vp/F"), 3),
  Method = rep(c("nlmixr2 SAEM", "NONMEM SAEM", "Bayesian (Stan)"), each = 5),
  True = rep(c(TRUE_PARAMS$ka, TRUE_PARAMS$cl_f, TRUE_PARAMS$vc_f,
               TRUE_PARAMS$q_f, TRUE_PARAMS$vp_f), 3),
  Estimate = c(
    nlmixr2_natural$Ka, nlmixr2_natural$CL, nlmixr2_natural$V1,
    nlmixr2_natural$Q, nlmixr2_natural$V2,
    nonmem_est$Ka, nonmem_est$CL, nonmem_est$V1, nonmem_est$Q, nonmem_est$V2,
    bayes_med["Ka_pop"], bayes_med["CL_pop"], bayes_med["V1_pop"],
    bayes_med["Q_pop"], bayes_med["V2_pop"]
  )
) %>%
  mutate(Ratio = Estimate / True)

p_ratio <- ggplot(ratio_data, aes(x = Parameter, y = Ratio, color = Method, shape = Method)) +
  geom_hline(yintercept = 1.0, linetype = "solid", color = "grey40") +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0.80, ymax = 1.25),
            fill = "grey90", alpha = 0.01, inherit.aes = FALSE) +
  geom_hline(yintercept = c(0.80, 1.25), linetype = "dashed", color = "grey60") +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = method_colors) +
  labs(x = NULL, y = "Estimate / True",
       title = "Parameter Accuracy — Ratio to OSP True Values",
       subtitle = "Shaded band = 0.80-1.25 (clinical relevance heuristic)",
       color = NULL, shape = NULL) +
  coord_cartesian(ylim = c(0.5, 1.8)) +
  theme_poppk +
  theme(legend.position = "top")

ggsave(file.path(DIR_BAYESIAN, "comparison_ratio.png"), p_ratio,
       width = 10, height = 6, dpi = 150)
log_msg("  Saved comparison_ratio.png")

# --- Summary -----------------------------------------------------------------
log_msg("=== Three-Way Method Comparison Complete ===")
log_msg("Output files:")
log_msg("  Tables: three_way_comparison.csv, bsv_comparison_bayesian.csv, loo_summary.csv")
log_msg("  Figures: comparison_forest.png, comparison_pk_profile.png, comparison_ratio.png")
