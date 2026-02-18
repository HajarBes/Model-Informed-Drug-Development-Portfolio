# =============================================================================
# 07_cyp3a4_inhibition_bridge.R — CYP3A4 inhibition sensitivity scenario
# =============================================================================
# Conceptual bridge between baseline PopPK and DDI risk assessment:
#   - Reduces CL/F by 50% (strong CYP3A4 inhibitor, fm~1.0 scenario)
#   - Simulates AUC distribution fold-change vs uninhibited baseline
#   - Produces 1 figure + summary statistics
#
# This is NOT a PBPK-DDI model. It is a PopPK-based sensitivity analysis
# illustrating how baseline PK informs DDI risk quantification.
# =============================================================================

script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(script_dir) || !file.exists(file.path(script_dir, "00_setup.R"))) {
  script_dir <- file.path(getwd(), "analysis")
}
source(file.path(script_dir, "00_setup.R"))

log_msg("=== CYP3A4 Inhibition Bridge Scenario ===")

set.seed(20240315)

N_SIM  <- 1000
T_SIM  <- seq(0, 24, by = 0.1)
DOSE   <- 7500  # 7.5 mg
WT_REF <- 70

# CYP3A4 inhibition scenario:
# Strong inhibitor reduces CL/F by ~50% (equivalent to AUC fold-change ~2x)
# This reflects midazolam fm(CYP3A4) near 1.0 and strong inhibition
CL_REDUCTION <- 0.50

log_msg("Scenario: Strong CYP3A4 inhibitor (CL reduced by ",
        round(CL_REDUCTION * 100), "%)")
log_msg("N = ", N_SIM, " virtual subjects at 70 kg reference weight")

# --- Simulate baseline and inhibited -----------------------------------------

results <- list()

for (i in 1:N_SIM) {
  eta_ka <- rnorm(1, 0, sqrt(TRUE_OMEGA$ka))
  eta_cl <- rnorm(1, 0, sqrt(TRUE_OMEGA$cl))
  eta_vc <- rnorm(1, 0, sqrt(TRUE_OMEGA$vc))

  ind_ka <- TRUE_PARAMS$ka   * exp(eta_ka)
  ind_cl <- TRUE_PARAMS$cl_f * exp(eta_cl)
  ind_vc <- TRUE_PARAMS$vc_f * exp(eta_vc)
  ind_q  <- TRUE_PARAMS$q_f
  ind_vp <- TRUE_PARAMS$vp_f

  # Baseline (uninhibited)
  pk_base <- simulate_subject_pk(
    dose_ug = DOSE, times = T_SIM[T_SIM > 0],
    ka = ind_ka, cl = ind_cl, v1 = ind_vc, q = ind_q, v2 = ind_vp
  )

  # Inhibited (CL reduced by 50%)
  pk_inhib <- simulate_subject_pk(
    dose_ug = DOSE, times = T_SIM[T_SIM > 0],
    ka = ind_ka, cl = ind_cl * (1 - CL_REDUCTION), v1 = ind_vc,
    q = ind_q, v2 = ind_vp
  )

  auc_base  <- compute_auc(pk_base$time, pk_base$CP)
  auc_inhib <- compute_auc(pk_inhib$time, pk_inhib$CP)

  results[[i]] <- tibble(
    sim_id       = i,
    AUC_baseline = auc_base,
    AUC_inhibited = auc_inhib,
    fold_change  = auc_inhib / auc_base
  )
}

ddi_df <- bind_rows(results)

# --- Summary statistics -------------------------------------------------------

fc_summary <- tibble(
  scenario         = "Strong CYP3A4 inhibition (50% CL reduction)",
  n                = nrow(ddi_df),
  fc_median        = round(median(ddi_df$fold_change), 2),
  fc_mean          = round(mean(ddi_df$fold_change), 2),
  fc_p5            = round(quantile(ddi_df$fold_change, 0.05), 2),
  fc_p95           = round(quantile(ddi_df$fold_change, 0.95), 2),
  auc_base_median  = round(median(ddi_df$AUC_baseline), 1),
  auc_inhib_median = round(median(ddi_df$AUC_inhibited), 1),
  pct_fc_above_2   = round(100 * mean(ddi_df$fold_change > 2), 1),
  pct_fc_above_5   = round(100 * mean(ddi_df$fold_change > 5), 1)
)

log_msg("")
log_msg("Results:")
log_msg("  Median AUC fold-change: ", fc_summary$fc_median, "x")
log_msg("  90% PI: [", fc_summary$fc_p5, ", ", fc_summary$fc_p95, "]")
log_msg("  Baseline median AUC: ", fc_summary$auc_base_median, " ng*h/mL")
log_msg("  Inhibited median AUC: ", fc_summary$auc_inhib_median, " ng*h/mL")
log_msg("  Subjects with FC > 2x: ", fc_summary$pct_fc_above_2, "%")

# --- Figure: fold-change distribution ----------------------------------------

fig_ddi <- ggplot(ddi_df, aes(x = fold_change)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40,
                 fill = COL_PRIMARY, alpha = 0.7, color = "white") +
  geom_density(color = COL_SECONDARY, linewidth = 0.8) +
  geom_vline(xintercept = median(ddi_df$fold_change),
             linetype = "solid", color = "black", linewidth = 0.7) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "grey40") +
  annotate("text", x = median(ddi_df$fold_change), y = Inf,
           label = paste0("Median: ", fc_summary$fc_median, "x"),
           vjust = 2, hjust = -0.1, fontface = "bold", size = 3.5) +
  annotate("text", x = 2, y = Inf,
           label = "2x threshold",
           vjust = 3.5, hjust = -0.1, color = "grey40", size = 3) +
  labs(
    title = "Predicted AUC Fold-Change Under CYP3A4 Inhibition",
    subtitle = paste0(
      "50% CL/F reduction (strong inhibitor) | N=", N_SIM,
      " | Median FC = ", fc_summary$fc_median, "x [90% PI: ",
      fc_summary$fc_p5, "-", fc_summary$fc_p95, "]"
    ),
    x = "AUC Fold-Change (inhibited / baseline)",
    y = "Density",
    caption = paste0(
      "Conceptual sensitivity analysis: CL/F reduced by 50% to approximate ",
      "strong CYP3A4 inhibition.\n",
      "This links baseline PopPK characterization to DDI risk logic ",
      "(see P10 Static DDI Framework for formal assessment)."
    )
  )

ggsave(file.path(DIR_FIGURES, "cyp3a4_inhibition_fold_change.png"), fig_ddi,
       width = 7, height = 5, dpi = 300, bg = "white")
log_msg("Saved: cyp3a4_inhibition_fold_change.png")

write_csv(fc_summary, file.path(DIR_TABLES, "cyp3a4_inhibition_summary.csv"))
log_msg("Saved: cyp3a4_inhibition_summary.csv")

log_msg("")
log_msg("=== CYP3A4 Inhibition Bridge Complete ===")
