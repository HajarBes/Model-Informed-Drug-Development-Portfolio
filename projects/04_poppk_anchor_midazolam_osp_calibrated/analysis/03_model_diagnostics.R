# =============================================================================
# 03_model_diagnostics.R — Production-grade GOF, VPC, eta diagnostics
# =============================================================================
# Each figure answers a specific model qualification question:
#   GOF:  Does the model describe the observed data?
#   CWRES: Is there systematic bias in predictions?
#   Etas:  Are random effects well-estimated and normally distributed?
#   VPC:   Does the model capture population variability?
# =============================================================================

script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(script_dir) || !file.exists(file.path(script_dir, "00_setup.R"))) {
  script_dir <- file.path(getwd(), "analysis")
}
source(file.path(script_dir, "00_setup.R"))

log_msg("=== Starting Model Diagnostics ===")

# Load base model fit
fit <- readRDS(file.path(DIR_TABLES, "fit_base.rds"))

if (is.null(fit)) {
  log_msg("No model fit available — skipping diagnostics")
  quit(save = "no")
}

# Extract diagnostic data
diag_data <- as.data.frame(fit) %>% as_tibble()

SAVE_W <- 8
SAVE_H <- 6

# =============================================================================
# FIGURE: GOODNESS-OF-FIT (4-PANEL)
# =============================================================================

# Panel 1: DV vs PRED
p1 <- ggplot(diag_data, aes(x = PRED, y = DV)) +
  geom_point(alpha = 0.3, size = 1, color = COL_PRIMARY) +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "black") +
  scale_x_log10() + scale_y_log10() +
  labs(title = "DV vs Population Prediction", x = "PRED (ng/mL)", y = "DV (ng/mL)")

# Panel 2: DV vs IPRED
p2 <- ggplot(diag_data, aes(x = IPRED, y = DV)) +
  geom_point(alpha = 0.3, size = 1, color = COL_IPRED) +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "black") +
  scale_x_log10() + scale_y_log10() +
  labs(title = "DV vs Individual Prediction", x = "IPRED (ng/mL)", y = "DV (ng/mL)")

# Panel 3: CWRES vs TIME
cwres_col <- if ("CWRES" %in% names(diag_data)) "CWRES" else "IRES"
p3 <- ggplot(diag_data, aes(x = TIME, y = .data[[cwres_col]])) +
  geom_point(alpha = 0.3, size = 1, color = COL_OBSERVED) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "grey50") +
  geom_smooth(method = "loess", se = FALSE, color = COL_SECONDARY, linewidth = 0.8) +
  labs(title = paste0(cwres_col, " vs Time"), x = "Time (h)", y = cwres_col)

# Panel 4: CWRES vs PRED
p4 <- ggplot(diag_data, aes(x = PRED, y = .data[[cwres_col]])) +
  geom_point(alpha = 0.3, size = 1, color = COL_OBSERVED) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "grey50") +
  geom_smooth(method = "loess", se = FALSE, color = COL_SECONDARY, linewidth = 0.8) +
  labs(title = paste0(cwres_col, " vs PRED"), x = "PRED (ng/mL)", y = cwres_col)

fig_gof <- (p1 + p2) / (p3 + p4) +
  plot_annotation(
    title = "Goodness-of-Fit: Base PopPK Model",
    subtitle = "Oral midazolam 7.5 mg | 2-compartment with first-order absorption",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )
  )

ggsave(file.path(DIR_FIGURES, "gof_4panel.png"), fig_gof,
       width = 10, height = 8, dpi = 300, bg = "white")
log_msg("Saved: gof_4panel.png")

# =============================================================================
# FIGURE: ETA DISTRIBUTIONS + SHRINKAGE
# =============================================================================

# Extract etas
eta_cols <- grep("^eta\\.", names(diag_data), value = TRUE)
if (length(eta_cols) == 0) eta_cols <- grep("^ETA", names(diag_data), value = TRUE)

if (length(eta_cols) > 0) {
  eta_long <- diag_data %>%
    select(ID, all_of(eta_cols)) %>%
    distinct() %>%
    pivot_longer(-ID, names_to = "eta", values_to = "value")

  # Compute shrinkage
  shrinkage <- eta_long %>%
    group_by(eta) %>%
    summarise(
      mean = round(mean(value), 3),
      sd   = round(sd(value), 3),
      shrinkage_pct = round((1 - sd(value) / sd(value[1])) * 100, 1),
      .groups = "drop"
    )

  # Compute eta-shrinkage properly: 1 - SD(eta_i) / omega
  # Use variance from fit if available
  omega_diag <- tryCatch(diag(fit$omega), error = function(e) NULL)
  if (!is.null(omega_diag)) {
    omega_sd <- sqrt(omega_diag[seq_along(eta_cols)])
    shrinkage$omega_sd <- omega_sd
    shrinkage$shrinkage_pct <- round((1 - shrinkage$sd / omega_sd) * 100, 1)
  }

  fig_etas <- ggplot(eta_long, aes(x = value)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30,
                   fill = COL_PRIMARY, alpha = 0.6, color = "white") +
    geom_density(color = COL_SECONDARY, linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
    facet_wrap(~eta, scales = "free") +
    labs(
      title = "Random Effect (eta) Distributions",
      subtitle = "Dashed line = 0 (population typical value)",
      x = "Eta value", y = "Density"
    )

  ggsave(file.path(DIR_FIGURES, "eta_distributions.png"), fig_etas,
         width = SAVE_W, height = SAVE_H, dpi = 300, bg = "white")
  log_msg("Saved: eta_distributions.png")

  write_csv(shrinkage, file.path(DIR_TABLES, "eta_shrinkage.csv"))
  log_msg("Eta shrinkage computed:")
  for (i in seq_len(nrow(shrinkage))) {
    log_msg(sprintf("  %s: SD=%.3f, shrinkage=%.1f%%",
                    shrinkage$eta[i], shrinkage$sd[i], shrinkage$shrinkage_pct[i]))
  }
}

# =============================================================================
# FIGURE: VISUAL PREDICTIVE CHECK (VPC)
# =============================================================================

# Simulate from fitted model for VPC
log_msg("Generating VPC simulations...")

vpc_available <- requireNamespace("vpc", quietly = TRUE)

if (vpc_available && !is.null(fit)) {
  library(vpc)

  tryCatch({
    vpc_data <- vpcSim(fit, n = 500)

    vpc_fig <- vpcPlot(fit, n = 500,
                       show = list(obs_dv = TRUE, obs_ci = TRUE,
                                   pi = TRUE, pi_ci = TRUE),
                       log_y = TRUE) +
      labs(
        title = "Visual Predictive Check (VPC)",
        subtitle = "500 simulations | Shaded = 90% prediction interval | Lines = observed percentiles",
        x = "Time (h)", y = "Concentration (ng/mL)",
        caption = paste0(
          "Note: Lower PI widens at t > 8 h on log scale. This reflects sparse ",
          "observations (only rich-arm subjects sampled at 12 h)\n",
          "and simulated concentrations approaching the LLOQ (0.5 ng/mL). ",
          "Log-scale amplifies variability near zero. Population-level ",
          "model adequacy\nis best assessed in the 0-8 h window where both ",
          "sampling arms contribute data."
        )
      ) +
      theme(plot.caption = element_text(size = 7, color = "grey40", hjust = 0))

    ggsave(file.path(DIR_FIGURES, "vpc.png"), vpc_fig,
           width = SAVE_W, height = SAVE_H + 0.5, dpi = 300, bg = "white")
    log_msg("Saved: vpc.png")
  }, error = function(e) {
    log_msg("VPC generation failed: ", conditionMessage(e))
    log_msg("Generating manual VPC from simulated data...")

    # Manual VPC: simulate and overlay percentiles
    sim_list <- list()
    for (rep in 1:200) {
      sim_df <- as.data.frame(simulate(fit)) %>%
        select(TIME, sim) %>%
        mutate(REP = rep)
      sim_list[[rep]] <- sim_df
    }
    sim_all <- bind_rows(sim_list)

    sim_pctiles <- sim_all %>%
      group_by(TIME) %>%
      summarise(
        p5  = quantile(sim, 0.05, na.rm = TRUE),
        p50 = quantile(sim, 0.50, na.rm = TRUE),
        p95 = quantile(sim, 0.95, na.rm = TRUE),
        .groups = "drop"
      )

    obs_pctiles <- diag_data %>%
      group_by(TIME) %>%
      summarise(
        obs_p50 = median(DV, na.rm = TRUE),
        .groups = "drop"
      )

    vpc_manual <- ggplot() +
      geom_ribbon(data = sim_pctiles, aes(x = TIME, ymin = p5, ymax = p95),
                  fill = COL_PRIMARY, alpha = 0.2) +
      geom_line(data = sim_pctiles, aes(x = TIME, y = p50),
                color = COL_PRIMARY, linewidth = 0.8) +
      geom_point(data = obs_pctiles, aes(x = TIME, y = obs_p50),
                 color = COL_SECONDARY, size = 2) +
      scale_y_log10() +
      labs(
        title = "Visual Predictive Check (Manual)",
        subtitle = "Blue band = 90% PI from 200 simulations | Red = observed median",
        x = "Time (h)", y = "Concentration (ng/mL)"
      )

    ggsave(file.path(DIR_FIGURES, "vpc.png"), vpc_manual,
           width = SAVE_W, height = SAVE_H, dpi = 300, bg = "white")
    log_msg("Saved: vpc.png (manual)")
  })
} else {
  log_msg("vpc package not available — skipping VPC")
}

# =============================================================================
# FIGURE: INDIVIDUAL FITS (SELECTED SUBJECTS)
# =============================================================================

# Show 9 representative subjects
selected_ids <- sort(sample(unique(diag_data$ID), min(9, n_distinct(diag_data$ID))))

ind_data <- diag_data %>%
  filter(ID %in% selected_ids)

fig_ind <- ggplot(ind_data, aes(x = TIME)) +
  geom_point(aes(y = DV), color = COL_OBSERVED, size = 1.5) +
  geom_line(aes(y = IPRED), color = COL_IPRED, linewidth = 0.8) +
  geom_line(aes(y = PRED), color = COL_PRED, linewidth = 0.6, linetype = "dashed") +
  facet_wrap(~ID, ncol = 3, scales = "free_y") +
  labs(
    title = "Individual Fits (Selected Subjects)",
    subtitle = "Points = observed | Solid = IPRED | Dashed = PRED",
    x = "Time (h)", y = "Concentration (ng/mL)"
  )

ggsave(file.path(DIR_FIGURES, "individual_fits.png"), fig_ind,
       width = 10, height = 8, dpi = 300, bg = "white")
log_msg("Saved: individual_fits.png")

log_msg("=== Model Diagnostics Complete ===")
