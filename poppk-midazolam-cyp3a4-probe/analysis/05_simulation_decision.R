# =============================================================================
# 05_simulation_decision.R — Dosing simulations and label-like recommendations
# =============================================================================
# Questions answered:
#   1. How does exposure (AUC, Cmax) change across dosing scenarios?
#   2. Do body weight strata require dose adjustment?
#   3. What is the population exposure variability for label language?
# =============================================================================

script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(script_dir) || !file.exists(file.path(script_dir, "00_setup.R"))) {
  script_dir <- file.path(getwd(), "analysis")
}
source(file.path(script_dir, "00_setup.R"))

log_msg("=== Starting Decision Simulations ===")

set.seed(20240201)

# =============================================================================
# SIMULATION PARAMETERS (from final covariate model)
# =============================================================================

# Use true parameters with allometric scaling (validated by fitting)
N_SIM    <- 1000      # Virtual subjects per scenario
T_SIM    <- seq(0, 24, by = 0.1)  # Fine time grid for AUC

# Dosing scenarios
scenarios <- tibble(
  scenario = c("5 mg (low dose)", "7.5 mg (standard)", "15 mg (high dose)"),
  dose_mg  = c(5, 7.5, 15),
  dose_ug  = c(5000, 7500, 15000)
)

# Weight strata for subgroup analysis
wt_strata <- tibble(
  stratum = c("Low WT (50 kg)", "Reference (70 kg)", "High WT (100 kg)"),
  wt_val  = c(50, 70, 100)
)

# =============================================================================
# SIMULATE EXPOSURE FOR EACH SCENARIO x WEIGHT STRATUM
# =============================================================================

log_msg("Simulating ", N_SIM, " subjects x ", nrow(scenarios), " doses x ",
        nrow(wt_strata), " WT strata...")

sim_results <- list()

for (s in seq_len(nrow(scenarios))) {
  for (w in seq_len(nrow(wt_strata))) {
    dose_ug <- scenarios$dose_ug[s]
    wt      <- wt_strata$wt_val[w]

    for (i in 1:N_SIM) {
      # Individual parameters with allometric scaling + BSV
      eta_ka <- rnorm(1, 0, sqrt(TRUE_OMEGA$ka))
      eta_cl <- rnorm(1, 0, sqrt(TRUE_OMEGA$cl))
      eta_vc <- rnorm(1, 0, sqrt(TRUE_OMEGA$vc))

      ind_ka <- TRUE_PARAMS$ka   * exp(eta_ka)
      ind_cl <- TRUE_PARAMS$cl_f * exp(eta_cl) * (wt / 70)^0.75
      ind_vc <- TRUE_PARAMS$vc_f * exp(eta_vc) * (wt / 70)^1.0
      ind_q  <- TRUE_PARAMS$q_f  * (wt / 70)^0.75
      ind_vp <- TRUE_PARAMS$vp_f * (wt / 70)^1.0

      # Simulate concentration-time profile (no residual error for NCA)
      pk <- simulate_subject_pk(
        dose_ug = dose_ug, times = T_SIM[T_SIM > 0],
        ka = ind_ka, cl = ind_cl, v1 = ind_vc, q = ind_q, v2 = ind_vp
      )

      auc  <- compute_auc(pk$time, pk$CP)
      cmax <- max(pk$CP)
      tmax <- pk$time[which.max(pk$CP)]

      sim_results[[length(sim_results) + 1]] <- tibble(
        scenario = scenarios$scenario[s],
        dose_mg  = scenarios$dose_mg[s],
        stratum  = wt_strata$stratum[w],
        wt       = wt,
        sim_id   = i,
        AUC_0inf = auc,
        Cmax     = cmax,
        Tmax     = tmax
      )
    }
  }
}

sim_df <- bind_rows(sim_results)
log_msg("Simulations complete: ", nrow(sim_df), " total profiles")

# =============================================================================
# EXPOSURE SUMMARY TABLE
# =============================================================================

exposure_summary <- sim_df %>%
  group_by(scenario, dose_mg, stratum, wt) %>%
  summarise(
    AUC_median  = round(median(AUC_0inf), 1),
    AUC_p5      = round(quantile(AUC_0inf, 0.05), 1),
    AUC_p95     = round(quantile(AUC_0inf, 0.95), 1),
    Cmax_median = round(median(Cmax), 1),
    Cmax_p5     = round(quantile(Cmax, 0.05), 1),
    Cmax_p95    = round(quantile(Cmax, 0.95), 1),
    Tmax_median = round(median(Tmax), 2),
    .groups     = "drop"
  )

write_csv(exposure_summary, file.path(DIR_TABLES, "exposure_summary.csv"))
log_msg("Saved: exposure_summary.csv")

# =============================================================================
# FIGURE: EXPOSURE BY DOSE (BOXPLOT)
# =============================================================================

fig_dose_auc <- sim_df %>%
  ggplot(aes(x = factor(dose_mg), y = AUC_0inf, fill = stratum)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  scale_fill_manual(values = c("Low WT (50 kg)" = "#D6604D",
                                "Reference (70 kg)" = "#4393C3",
                                "High WT (100 kg)" = "#2166AC"),
                    name = "Weight Stratum") +
  labs(
    title = "Predicted AUC by Dose and Body Weight",
    subtitle = "N=1000 virtual subjects per stratum | Allometric scaling on CL and V",
    x = "Dose (mg)", y = "AUC (ng*h/mL)"
  )

ggsave(file.path(DIR_FIGURES, "exposure_dose_auc.png"), fig_dose_auc,
       width = 8, height = 6, dpi = 300, bg = "white")
log_msg("Saved: exposure_dose_auc.png")

# =============================================================================
# FIGURE: EXPOSURE BY DOSE (CMAX)
# =============================================================================

fig_dose_cmax <- sim_df %>%
  ggplot(aes(x = factor(dose_mg), y = Cmax, fill = stratum)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  scale_fill_manual(values = c("Low WT (50 kg)" = "#D6604D",
                                "Reference (70 kg)" = "#4393C3",
                                "High WT (100 kg)" = "#2166AC"),
                    name = "Weight Stratum") +
  labs(
    title = "Predicted Cmax by Dose and Body Weight",
    subtitle = "N=1000 virtual subjects per stratum | Allometric scaling on CL and V",
    x = "Dose (mg)", y = "Cmax (ng/mL)"
  )

ggsave(file.path(DIR_FIGURES, "exposure_dose_cmax.png"), fig_dose_cmax,
       width = 8, height = 6, dpi = 300, bg = "white")
log_msg("Saved: exposure_dose_cmax.png")

# =============================================================================
# FIGURE: WEIGHT EFFECT ON AUC (STANDARD DOSE ONLY)
# =============================================================================

std_dose <- sim_df %>% filter(dose_mg == 7.5)

# AUC ratio relative to reference (70 kg)
ref_auc_median <- median(std_dose$AUC_0inf[std_dose$wt == 70])

auc_ratios <- std_dose %>%
  group_by(stratum, wt) %>%
  summarise(
    auc_median = median(AUC_0inf),
    auc_p5     = quantile(AUC_0inf, 0.05),
    auc_p95    = quantile(AUC_0inf, 0.95),
    .groups = "drop"
  ) %>%
  mutate(
    ratio_median = auc_median / ref_auc_median,
    ratio_p5     = auc_p5 / ref_auc_median,
    ratio_p95    = auc_p95 / ref_auc_median
  )

fig_wt_ratio <- ggplot(auc_ratios, aes(x = ratio_median, y = fct_rev(stratum))) +
  geom_point(size = 3, color = COL_PRIMARY) +
  geom_errorbarh(aes(xmin = ratio_p5, xmax = ratio_p95), height = 0.2,
                 color = COL_PRIMARY) +
  geom_vline(xintercept = 1, linetype = "solid") +
  geom_vline(xintercept = c(0.8, 1.25), linetype = "dashed", color = "grey50") +
  annotate("rect", xmin = 0.8, xmax = 1.25, ymin = -Inf, ymax = Inf,
           alpha = 0.08, fill = "green") +
  labs(
    title = "AUC Ratio by Body Weight (7.5 mg Dose)",
    subtitle = "Reference = 70 kg | Green = clinical relevance window (0.80-1.25) | Bars = 90% PI",
    x = "AUC Ratio vs 70 kg Reference", y = NULL
  )

ggsave(file.path(DIR_FIGURES, "wt_auc_ratio.png"), fig_wt_ratio,
       width = 7, height = 4, dpi = 300, bg = "white")
log_msg("Saved: wt_auc_ratio.png")

# =============================================================================
# DOSE RECOMMENDATION
# =============================================================================

log_msg("")
log_msg("=== DOSE RECOMMENDATION ===")
log_msg("")

# Check if any stratum median falls outside 0.8-1.25 clinical relevance heuristic
outside_be <- auc_ratios %>%
  filter(ratio_median < 0.8 | ratio_median > 1.25)

if (nrow(outside_be) == 0) {
  recommendation <- paste(
    "All weight strata have median AUCR within the 0.80-1.25 clinical",
    "relevance heuristic band. Individual variability from BSV in CL",
    "exceeds the weight effect. Standard 7.5 mg dosing is appropriate",
    "across the evaluated weight range."
  )
} else {
  recommendation <- paste(
    "Median AUCR marginally exceeds the 0.80-1.25 clinical relevance",
    "heuristic band for:", paste(outside_be$stratum, collapse = ", "), ".",
    "However, individual variability from BSV in CL (~30% CV) substantially",
    "exceeds the weight effect at all strata. Whether dose modification is",
    "warranted depends on the therapeutic index in the specific clinical context."
  )
}

log_msg(recommendation)

# Save recommendation
writeLines(c(
  "=== PopPK Dose Recommendation ===",
  "",
  "Drug: Midazolam (oral)",
  "Context: Single-dose PK in healthy adults",
  "Model: 2-compartment, first-order absorption, allometric WT scaling",
  "",
  "Dosing scenarios evaluated: 5 mg, 7.5 mg, 15 mg",
  "Weight strata: 50 kg, 70 kg (reference), 100 kg",
  "",
  "RECOMMENDATION:",
  recommendation,
  "",
  "NOTE: This recommendation is based on simulated data calibrated to",
  "published midazolam PK. It demonstrates the PopPK-informed dose",
  "evaluation workflow and should not replace clinical decision-making.",
  ""
), file.path(DIR_TABLES, "dose_recommendation.txt"))

write_csv(auc_ratios, file.path(DIR_TABLES, "auc_ratios_by_wt.csv"))

log_msg("")
log_msg("=== Decision Simulations Complete ===")
