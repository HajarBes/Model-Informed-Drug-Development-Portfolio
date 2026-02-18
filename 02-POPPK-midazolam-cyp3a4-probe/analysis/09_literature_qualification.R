# =============================================================================
# 09_literature_qualification.R
# External Literature Qualification: overlay model predictions on OSP profiles
# =============================================================================
#
# This script overlays PopPK model-predicted typical concentration-time curves
# on digitized mean +/- SD profiles from the OSP observed-data database.
# It assesses structural consistency only — not population variability or
# covariate effects.
#
# Inputs:
#   data/literature_osp_midazolam/extracted/mean_profiles_midazolam_po.csv
#   data/literature_osp_midazolam/qualification_set_ids.txt
#   outputs/tables/parameter_estimates_base.csv
#
# Outputs:
#   figures/lit_overlay_mean_profiles.png
#   figures/lit_dose_stratified_overlays.png
#   figures/lit_external_vpc_like.png
#   figures/lit_qualification_summary.png
#   outputs/tables/literature_qualification_summary.csv

# --- Setup ------------------------------------------------------------------
source(file.path("analysis", "00_setup.R"))

log_msg("09_literature_qualification.R started.")

# --- Load fitted THETAs from SAEM estimation --------------------------------
param_file <- file.path(DIR_TABLES, "parameter_estimates_base.csv")
param_df   <- read_csv(param_file, show_col_types = FALSE)

get_param <- function(name) {
  val <- param_df$Estimate[param_df$Parameter == name]
  if (length(val) == 0) stop(paste("Parameter not found:", name))
  val
}

# Exponentiate log-transformed parameters
FITTED <- list(
  ka   = exp(get_param("lka")),   # 2.26 h^-1

cl_f = exp(get_param("lcl")),   # 52.3 L/h
  vc_f = exp(get_param("lv1")),   # 46.6 L
  q_f  = exp(get_param("lq")),    # 16.8 L/h
  vp_f = exp(get_param("lv2"))    # 62.5 L
)

log_msg(sprintf("Fitted THETAs: Ka=%.2f, CL/F=%.1f, Vc/F=%.1f, Q/F=%.1f, Vp/F=%.1f",
                FITTED$ka, FITTED$cl_f, FITTED$vc_f, FITTED$q_f, FITTED$vp_f))

# --- Load qualification set -------------------------------------------------
DIR_LIT <- file.path(DIR_DATA, "literature_osp_midazolam")
qual_ids <- as.integer(readLines(file.path(DIR_LIT, "qualification_set_ids.txt")))
qual_ids <- qual_ids[!is.na(qual_ids)]
log_msg(sprintf("Qualification set: %d studies", length(qual_ids)))

profiles <- read_csv(file.path(DIR_LIT, "extracted", "mean_profiles_midazolam_po.csv"),
                     show_col_types = FALSE) %>%
  filter(study_id %in% qual_ids)

log_msg(sprintf("Loaded %d data points from qualification set", nrow(profiles)))

# --- Assign dose groups -----------------------------------------------------
profiles <- profiles %>%
  mutate(
    dose_group = case_when(
      dose_mg <= 1   ~ "1 mg",
      dose_mg <= 3   ~ "3 mg",
      dose_mg <= 7.5 ~ "5-7.5 mg",
      TRUE           ~ "15 mg"
    ),
    dose_group = factor(dose_group, levels = c("1 mg", "3 mg", "5-7.5 mg", "15 mg"))
  )

# Unique study label for display
profiles <- profiles %>%
  mutate(study_display = paste0(study_label, " (", dose_mg, " mg, N=", n, ")"))

# --- Generate typical curves per study --------------------------------------
log_msg("Generating typical predicted curves per study...")

# Get max observed time per study
study_info <- profiles %>%
  group_by(study_id, study_label, dose_mg, dose_group, study_display) %>%
  summarise(
    t_max_obs = max(time_h),
    n_tp = n(),
    n_subj = first(n),
    .groups = "drop"
  )

# Simulate typical curve for each study (per-study time grid)
typical_curves <- map_dfr(seq_len(nrow(study_info)), function(i) {
  si <- study_info[i, ]
  dose_ug <- si$dose_mg * 1000
  t_max <- si$t_max_obs * 1.1  # 10% buffer
  times <- seq(0.01, t_max, by = 0.1)

  pred <- simulate_subject_pk(
    dose_ug = dose_ug,
    times   = times,
    ka = FITTED$ka,
    cl = FITTED$cl_f,
    v1 = FITTED$vc_f,
    q  = FITTED$q_f,
    v2 = FITTED$vp_f
  )

  tibble(
    study_id      = si$study_id,
    study_label   = si$study_label,
    dose_mg       = si$dose_mg,
    dose_group    = si$dose_group,
    study_display = si$study_display,
    time_h        = pred$time,
    pred_ngml     = pred$CP
  )
})

log_msg(sprintf("Generated %d prediction points across %d studies",
                nrow(typical_curves), nrow(study_info)))

# --- Generate population PI bands per dose group ----------------------------
log_msg("Generating population PI bands (500 subjects/dose group)...")

set.seed(20240301)
N_SIM <- 500

dose_groups <- study_info %>%
  group_by(dose_group) %>%
  summarise(
    dose_mg  = first(dose_mg),
    t_max    = max(t_max_obs) * 1.1,
    .groups  = "drop"
  )

pi_bands <- map_dfr(seq_len(nrow(dose_groups)), function(g) {
  dg <- dose_groups[g, ]
  dose_ug <- dg$dose_mg * 1000
  times <- seq(0.01, dg$t_max, by = 0.1)

  # Simulate N_SIM subjects with IIV from TRUE_OMEGA
  sim_matrix <- matrix(NA_real_, nrow = length(times), ncol = N_SIM)
  for (j in seq_len(N_SIM)) {
    eta_ka <- rnorm(1, 0, sqrt(TRUE_OMEGA$ka))
    eta_cl <- rnorm(1, 0, sqrt(TRUE_OMEGA$cl))
    eta_vc <- rnorm(1, 0, sqrt(TRUE_OMEGA$vc))

    ka_j <- FITTED$ka   * exp(eta_ka)
    cl_j <- FITTED$cl_f * exp(eta_cl)
    vc_j <- FITTED$vc_f * exp(eta_vc)

    pred_j <- simulate_subject_pk(
      dose_ug = dose_ug,
      times   = times,
      ka = ka_j,
      cl = cl_j,
      v1 = vc_j,
      q  = FITTED$q_f,
      v2 = FITTED$vp_f
    )
    sim_matrix[, j] <- pred_j$CP
  }

  tibble(
    dose_group = dg$dose_group,
    dose_mg    = dg$dose_mg,
    time_h     = times,
    p05        = apply(sim_matrix, 1, quantile, probs = 0.05),
    p50        = apply(sim_matrix, 1, quantile, probs = 0.50),
    p95        = apply(sim_matrix, 1, quantile, probs = 0.95)
  )
})

log_msg(sprintf("PI bands: %d rows across %d dose groups", nrow(pi_bands), nrow(dose_groups)))

# --- Compute per-study qualification metrics --------------------------------
log_msg("Computing per-study qualification metrics...")

# For each study, interpolate predicted curve at observed time points
metrics <- map_dfr(seq_len(nrow(study_info)), function(i) {
  si <- study_info[i, ]

  obs <- profiles %>%
    filter(study_id == si$study_id, time_h > 0, conc_ngml > 0)
  if (nrow(obs) < 3) return(NULL)

  # Get typical curve for this study
  tc <- typical_curves %>% filter(study_id == si$study_id)

  # Interpolate predicted values at observed time points
  pred_at_obs <- approx(tc$time_h, tc$pred_ngml, xout = obs$time_h, rule = 2)$y

  # Unweighted RMSE (primary)
  resid <- pred_at_obs - obs$conc_ngml
  rmse <- sqrt(mean(resid^2))

  # MPE%
  mpe_pct <- mean(resid / obs$conc_ngml) * 100

  # AUC ratio (trapezoidal, observed time grid only)
  auc_pred <- compute_auc(obs$time_h, pred_at_obs)
  auc_obs  <- compute_auc(obs$time_h, obs$conc_ngml)
  auc_ratio <- if (!is.na(auc_pred) && !is.na(auc_obs) && auc_obs > 0) {
    auc_pred / auc_obs
  } else { NA_real_ }

  # Cmax ratio
  cmax_pred <- max(pred_at_obs, na.rm = TRUE)
  cmax_obs  <- max(obs$conc_ngml, na.rm = TRUE)
  cmax_ratio <- if (cmax_obs > 0) cmax_pred / cmax_obs else NA_real_

  # Weighted RMSE (secondary — sensitivity analysis)
  sd_vals <- obs$sd_ngml
  has_sd <- !is.na(sd_vals) & sd_vals > 0
  wrmse <- if (sum(has_sd) >= 3) {
    w <- 1 / sd_vals[has_sd]^2
    sqrt(sum(w * resid[has_sd]^2) / sum(w))
  } else { NA_real_ }

  tibble(
    study_id     = si$study_id,
    study_label  = si$study_label,
    dose_mg      = si$dose_mg,
    n            = si$n_subj,
    n_timepoints = nrow(obs),
    rmse_ngml    = round(rmse, 2),
    auc_ratio    = round(auc_ratio, 3),
    cmax_ratio   = round(cmax_ratio, 3),
    mpe_pct      = round(mpe_pct, 1),
    wrmse_ngml   = if (!is.na(wrmse)) round(wrmse, 2) else NA_real_,
    within_2fold = !is.na(auc_ratio) & auc_ratio >= 0.5 & auc_ratio <= 2.0
  )
})

log_msg("Qualification metrics computed:")
print(metrics)

# Write summary table
metrics_path <- file.path(DIR_TABLES, "literature_qualification_summary.csv")
write_csv(metrics, metrics_path)
log_msg(sprintf("Wrote %s", metrics_path))

# =============================================================================
# FIGURES
# =============================================================================

# --- Palettes ---------------------------------------------------------------
dose_colors <- c("1 mg" = "#1B9E77", "3 mg" = "#D95F02",
                 "5-7.5 mg" = "#7570B3", "15 mg" = "#E7298A")

# 12 distinct shapes for individual studies
study_labels <- sort(unique(profiles$study_label))
study_shapes <- setNames(c(16, 17, 15, 18, 3, 4, 8, 0, 1, 2, 5, 6)[seq_along(study_labels)],
                         study_labels)

# --- Figure 1: Overlay all mean profiles -----------------------------------
log_msg("Generating Figure 1: lit_overlay_mean_profiles.png")

p1 <- ggplot() +
  # Observed means as points
  geom_point(data = profiles %>% filter(conc_ngml > 0),
             aes(x = time_h, y = conc_ngml, color = dose_group,
                 shape = study_label),
             size = 1.8, alpha = 0.7) +
  # Error bars where SD available
  geom_errorbar(data = profiles %>% filter(conc_ngml > 0, !is.na(sd_ngml), sd_ngml > 0),
                aes(x = time_h,
                    ymin = pmax(conc_ngml - sd_ngml, 0.01),
                    ymax = conc_ngml + sd_ngml,
                    color = dose_group),
                width = 0, alpha = 0.3, linewidth = 0.3) +
  # Typical predicted curves
  geom_line(data = typical_curves,
            aes(x = time_h, y = pred_ngml, color = dose_group,
                group = study_id),
            linewidth = 0.7, alpha = 0.8) +
  scale_y_log10(labels = scales::label_number()) +
  scale_color_manual(values = dose_colors, name = "Dose") +
  scale_shape_manual(values = study_shapes) +
  coord_cartesian(ylim = c(0.1, NA)) +
  labs(
    title = "External Literature Qualification: Model vs OSP Mean Profiles",
    subtitle = "Lines = PopPK typical predictions (fitted THETAs); Points = OSP digitized means",
    x = "Time (h)", y = "Concentration (ng/mL, log scale)",
    shape = "Study",
    caption = paste0("Observed profiles are digitized from publications (OSP database). ",
                     "Uncertainty reflects reported SD/SEM\n",
                     "(converted where applicable) and digitization transfer error.")
  ) +
  guides(shape = guide_legend(ncol = 2, override.aes = list(size = 2))) +
  theme(legend.position = "right",
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"),
        plot.caption = element_text(size = 8, face = "italic", hjust = 0))

ggsave(file.path(DIR_FIGURES, "lit_overlay_mean_profiles.png"),
       p1, width = 12, height = 7, dpi = 300, bg = "white")

# --- Figure 2: Dose-stratified overlays with PI bands ----------------------
log_msg("Generating Figure 2: lit_dose_stratified_overlays.png")

obs_for_facet <- profiles %>%
  filter(conc_ngml > 0) %>%
  select(dose_group, study_label, time_h, conc_ngml, sd_ngml)

tc_for_facet <- typical_curves %>%
  select(dose_group, study_id, study_label, time_h, pred_ngml)

pi_for_facet <- pi_bands %>%
  select(dose_group, time_h, p05, p50, p95)

p2 <- ggplot() +
  # PI band (90%)
  geom_ribbon(data = pi_for_facet,
              aes(x = time_h, ymin = pmax(p05, 0.01), ymax = p95),
              fill = "grey80", alpha = 0.5) +
  # PI median
  geom_line(data = pi_for_facet,
            aes(x = time_h, y = p50),
            color = "grey50", linewidth = 0.5, linetype = "dashed") +
  # Typical curves per study
  geom_line(data = tc_for_facet,
            aes(x = time_h, y = pred_ngml, group = study_id),
            color = COL_PRIMARY, linewidth = 0.8) +
  # Observed means + SD
  geom_pointrange(data = obs_for_facet %>% filter(!is.na(sd_ngml), sd_ngml > 0),
                  aes(x = time_h, y = conc_ngml,
                      ymin = pmax(conc_ngml - sd_ngml, 0.01),
                      ymax = conc_ngml + sd_ngml,
                      shape = study_label),
                  size = 0.3, fatten = 2, color = COL_SECONDARY, alpha = 0.7) +
  geom_point(data = obs_for_facet %>% filter(is.na(sd_ngml) | sd_ngml == 0),
             aes(x = time_h, y = conc_ngml, shape = study_label),
             size = 1.5, color = COL_SECONDARY, alpha = 0.7) +
  facet_wrap(~ dose_group, scales = "free", ncol = 2) +
  scale_y_log10(labels = scales::label_number()) +
  scale_shape_manual(values = study_shapes) +
  coord_cartesian(ylim = c(0.1, NA)) +
  labs(
    title = "Dose-Stratified Literature Qualification",
    subtitle = paste0("Blue lines = typical predictions (fitted THETAs); ",
                      "Grey band = 90% PI from simulated PopPK model IIV (TRUE_OMEGA);\n",
                      "Red points = OSP digitized means +/- SD"),
    x = "Time (h)", y = "Concentration (ng/mL, log scale)",
    shape = "Study",
    caption = paste0("PI bands reflect IIV from the simulated PopPK model (TRUE_OMEGA), ",
                     "not literature-derived variability.")
  ) +
  guides(shape = guide_legend(ncol = 2, override.aes = list(size = 2))) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 7),
        plot.caption = element_text(size = 8, face = "italic", hjust = 0))

ggsave(file.path(DIR_FIGURES, "lit_dose_stratified_overlays.png"),
       p2, width = 11, height = 9, dpi = 300, bg = "white")

# --- Figure 3: External VPC-like overlay ------------------------------------
log_msg("Generating Figure 3: lit_external_vpc_like.png")

p3 <- ggplot() +
  # PI band (90%)
  geom_ribbon(data = pi_for_facet,
              aes(x = time_h, ymin = pmax(p05, 0.01), ymax = p95),
              fill = "#D4E6F1", alpha = 0.6) +
  # PI 50th
  geom_line(data = pi_for_facet,
            aes(x = time_h, y = p50),
            color = COL_PRIMARY, linewidth = 0.8) +
  # PI 5th and 95th
  geom_line(data = pi_for_facet,
            aes(x = time_h, y = p05),
            color = COL_PRIMARY, linewidth = 0.4, linetype = "dashed") +
  geom_line(data = pi_for_facet,
            aes(x = time_h, y = p95),
            color = COL_PRIMARY, linewidth = 0.4, linetype = "dashed") +
  # Observed means
  geom_point(data = obs_for_facet,
             aes(x = time_h, y = conc_ngml),
             color = COL_SECONDARY, size = 1.5, alpha = 0.6) +
  # Observed SD bars
  geom_errorbar(data = obs_for_facet %>% filter(!is.na(sd_ngml), sd_ngml > 0),
                aes(x = time_h,
                    ymin = pmax(conc_ngml - sd_ngml, 0.01),
                    ymax = conc_ngml + sd_ngml),
                color = COL_SECONDARY, width = 0, alpha = 0.3, linewidth = 0.3) +
  facet_wrap(~ dose_group, scales = "free", ncol = 2) +
  scale_y_log10(labels = scales::label_number()) +
  coord_cartesian(ylim = c(0.1, NA)) +
  labs(
    title = "External PI Overlay (not pcVPC)",
    subtitle = paste0("Blue = model-predicted PI (5th/50th/95th from simulated PopPK IIV); ",
                      "Red = OSP digitized means +/- SD"),
    x = "Time (h)", y = "Concentration (ng/mL, log scale)",
    caption = paste0("Prediction intervals are derived from the simulated PopPK model ",
                     "IIV (TRUE_OMEGA), not literature-derived variability.\n",
                     "Observed data are study-level means +/- SD, not individual ",
                     "concentrations. This is NOT a formal prediction-corrected VPC.")
  ) +
  theme(plot.caption = element_text(size = 8, face = "italic", hjust = 0))

ggsave(file.path(DIR_FIGURES, "lit_external_vpc_like.png"),
       p3, width = 11, height = 9, dpi = 300, bg = "white")

# --- Figure 4: AUC ratio bar chart -----------------------------------------
log_msg("Generating Figure 4: lit_qualification_summary.png")

metrics_plot <- metrics %>%
  mutate(
    bar_label = paste0(study_label, "\n(", dose_mg, " mg)"),
    bar_label = fct_reorder(bar_label, dose_mg),
    fill_cat  = ifelse(within_2fold, "Within 2-fold", "Outside 2-fold")
  )

p4 <- ggplot(metrics_plot, aes(x = bar_label, y = auc_ratio, fill = fill_cat)) +
  geom_col(width = 0.7, alpha = 0.85) +
  geom_hline(yintercept = 1.0, linewidth = 0.6, color = "black") +
  geom_hline(yintercept = c(0.5, 2.0), linewidth = 0.5,
             linetype = "dashed", color = "grey40") +
  annotate("text", x = 0.5, y = 2.05, label = "2-fold upper",
           hjust = 0, vjust = 0, size = 3, color = "grey40") +
  annotate("text", x = 0.5, y = 0.45, label = "2-fold lower",
           hjust = 0, vjust = 1, size = 3, color = "grey40") +
  scale_fill_manual(values = c("Within 2-fold" = "#2CA02C",
                               "Outside 2-fold" = "#D62728"),
                    name = "") +
  scale_y_continuous(breaks = seq(0, 3, by = 0.5)) +
  coord_cartesian(ylim = c(0, max(metrics$auc_ratio, na.rm = TRUE) * 1.2)) +
  labs(
    title = "Literature Qualification: AUC Ratio (Predicted / Observed)",
    subtitle = "Per-study comparison of model-predicted vs OSP-digitized AUC",
    x = NULL, y = "AUC Ratio (Predicted / Observed)",
    caption = paste0("The 0.5\u20132.0 AUC ratio band is used as a conventional benchmarking ",
                     "reference (guidance-aligned heuristic),\n",
                     "not a formal acceptance criterion. ",
                     "AUC computed by trapezoidal rule on observed time grid.")
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.caption = element_text(size = 8, face = "italic", hjust = 0))

ggsave(file.path(DIR_FIGURES, "lit_qualification_summary.png"),
       p4, width = 10, height = 6, dpi = 300, bg = "white")

# =============================================================================
# OPTIONAL: Heizmann 1983 micro-IIV sensitivity
# =============================================================================
log_msg("Heizmann 1983 individual curve analysis (exploratory)...")

heiz_file <- file.path(DIR_LIT, "extracted", "heizmann_1983_individual.csv")
if (file.exists(heiz_file)) {
  heiz <- read_csv(heiz_file, show_col_types = FALSE)

  if (nrow(heiz) > 0) {
    # Each grouping represents a different dose or individual
    heiz_groups <- unique(heiz$grouping)
    log_msg(sprintf("  Heizmann groupings: %s", paste(heiz_groups, collapse = ", ")))

    # Fit 2-comp model to each group using optim
    heiz_fits <- map_dfr(heiz_groups, function(grp) {
      dat <- heiz %>% filter(grouping == grp, time_h > 0, conc_ngml > 0)
      if (nrow(dat) < 5) return(NULL)

      # Determine dose from grouping string
      dose_mg <- tryCatch({
        # Try to extract dose from grouping
        m <- regmatches(grp, regexpr("[0-9]+\\.?[0-9]*\\s*mg", grp))
        if (length(m) > 0) as.numeric(gsub("\\s*mg", "", m[1]))
        else 15  # default
      }, error = function(e) 15)
      dose_ug <- dose_mg * 1000

      obj <- function(par) {
        ka_i <- exp(par[1])
        cl_i <- exp(par[2])
        vc_i <- exp(par[3])
        pred <- tryCatch({
          simulate_subject_pk(dose_ug, dat$time_h, ka_i, cl_i, vc_i,
                              FITTED$q_f, FITTED$vp_f)
        }, error = function(e) NULL)
        if (is.null(pred)) return(1e6)
        resid <- log(dat$conc_ngml) - log(pmax(pred$CP, 0.001))
        sum(resid^2)
      }

      init <- c(log(FITTED$ka), log(FITTED$cl_f), log(FITTED$vc_f))
      fit <- tryCatch(
        optim(init, obj, method = "Nelder-Mead",
              control = list(maxit = 5000)),
        error = function(e) NULL
      )

      if (is.null(fit) || fit$convergence != 0) return(NULL)

      tibble(
        grouping = grp,
        dose_mg  = dose_mg,
        ka       = exp(fit$par[1]),
        cl_f     = exp(fit$par[2]),
        vc_f     = exp(fit$par[3]),
        obj_val  = fit$value
      )
    })

    if (nrow(heiz_fits) >= 3) {
      log_msg("  Heizmann individual fit results:")
      print(heiz_fits)

      # Empirical BSV (SD of log-parameters)
      bsv <- tibble(
        parameter = c("Ka", "CL/F", "Vc/F"),
        sd_log    = c(sd(log(heiz_fits$ka)),
                      sd(log(heiz_fits$cl_f)),
                      sd(log(heiz_fits$vc_f))),
        cv_pct    = round(sqrt(exp(c(sd(log(heiz_fits$ka)),
                                     sd(log(heiz_fits$cl_f)),
                                     sd(log(heiz_fits$vc_f)))^2) - 1) * 100, 1)
      )
      log_msg("  Exploratory micro-IIV (N=6, insufficient for formal BSV estimation):")
      print(bsv)
    } else {
      log_msg("  Heizmann: fewer than 3 successful fits, skipping micro-IIV")
    }
  }
} else {
  log_msg("  Heizmann individual file not found, skipping")
}

# =============================================================================
log_msg("09_literature_qualification.R complete.")
log_msg(sprintf("Figures: %s/lit_*.png", DIR_FIGURES))
log_msg(sprintf("Table: %s", metrics_path))
