# =============================================================================
# 04_generate_tables.R — Generate formatted summary tables for the report
# =============================================================================

script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(script_dir) || !file.exists(file.path(script_dir, "00_setup.R"))) {
  script_dir <- file.path(getwd(), "analysis")
}
source(file.path(script_dir, "00_setup.R"))

log_msg("=== Starting Table Generation ===")

# Load all results
ddi_res <- readRDS(file.path(DIR_PROCESSED, "ddi_results.rds"))
sens    <- readRDS(file.path(DIR_PROCESSED, "sensitivity_results.rds"))

# =============================================================================
# TABLE 1: CROSS-CASE RISK SUMMARY
# =============================================================================

risk_all <- bind_rows(
  ddi_res[["midazolam_ketoconazole"]]$summary,
  ddi_res[["sotorasib_perpetrator"]]$summary
)

write_csv(risk_all, file.path(DIR_TABLES, "summary_all_cases.csv"))
log_msg("Saved: summary_all_cases.csv")

# =============================================================================
# TABLE 2: CONCENTRATION SUMMARY
# =============================================================================

conc_summary <- tibble(
  case = c("midazolam_ketoconazole", "sotorasib_perpetrator"),
  perpetrator = c("Ketoconazole", "Sotorasib"),
  dose_mg = c(ddi_res[["midazolam_ketoconazole"]]$params$dose_mg,
              ddi_res[["sotorasib_perpetrator"]]$params$dose_mg),
  MW = c(ddi_res[["midazolam_ketoconazole"]]$params$mw,
         ddi_res[["sotorasib_perpetrator"]]$params$mw),
  Cmax_total_ugml = c(ddi_res[["midazolam_ketoconazole"]]$params$cmax_total,
                      ddi_res[["sotorasib_perpetrator"]]$params$cmax_total),
  fu = c(ddi_res[["midazolam_ketoconazole"]]$params$fu,
         ddi_res[["sotorasib_perpetrator"]]$params$fu),
  Imax_u_uM = c(ddi_res[["midazolam_ketoconazole"]]$params$Imax_u,
                ddi_res[["sotorasib_perpetrator"]]$params$Imax_u),
  Igut_uM = c(ddi_res[["midazolam_ketoconazole"]]$params$Igut,
              ddi_res[["sotorasib_perpetrator"]]$params$Igut),
  Igut_over_Imaxu = Igut_uM / Imax_u_uM
)

write_csv(conc_summary, file.path(DIR_TABLES, "concentration_summary.csv"))
log_msg("Saved: concentration_summary.csv")

# =============================================================================
# TABLE 3: CYP SCREENING — ALL ENZYMES, ALL CASES
# =============================================================================

cyp_all <- bind_rows(
  ddi_res[["midazolam_ketoconazole"]]$cyp,
  ddi_res[["sotorasib_perpetrator"]]$cyp
)

write_csv(cyp_all, file.path(DIR_TABLES, "cyp_screening_all.csv"))
log_msg("Saved: cyp_screening_all.csv")

# =============================================================================
# TABLE 4: TRANSPORTER SCREENING — ALL CASES
# =============================================================================

trans_all <- bind_rows(
  ddi_res[["midazolam_ketoconazole"]]$transporter,
  ddi_res[["sotorasib_perpetrator"]]$transporter
)

write_csv(trans_all, file.path(DIR_TABLES, "transporter_screening_all.csv"))
log_msg("Saved: transporter_screening_all.csv")

# =============================================================================
# TABLE 5: VALIDATION COMPARISON
# =============================================================================

validation <- tibble(
  case = "midazolam_ketoconazole",
  perpetrator = "Ketoconazole",
  victim = "Midazolam",
  AUCR_predicted = ddi_res[["midazolam_ketoconazole"]]$aucr$AUCR_total,
  AUCR_observed = ddi_res[["midazolam_ketoconazole"]]$params$observed_aucr,
  pred_obs_ratio = AUCR_predicted / AUCR_observed,
  class_predicted = ddi_res[["midazolam_ketoconazole"]]$aucr$classification,
  class_observed = ddi_res[["midazolam_ketoconazole"]]$params$observed_class,
  class_match = class_predicted == class_observed
)

write_csv(validation, file.path(DIR_TABLES, "validation_comparison.csv"))
log_msg("Saved: validation_comparison.csv")

# =============================================================================
# TABLE 6: FULL PARAMETER TABLE FOR REPRODUCIBILITY
# =============================================================================

param_table <- tibble(
  parameter = c("Ki_CYP3A4", "Imax_u", "Igut", "fm_CYP3A4", "fg_CYP3A4",
                "kinact", "KI", "kdeg_hepatic", "Emax_CYP3A4", "EC50_CYP3A4"),
  units = c("uM", "uM", "uM", "fraction", "fraction",
            "1/min", "uM", "1/min", "fold", "uM"),
  ketoconazole = c(0.015,
                   ddi_res[["midazolam_ketoconazole"]]$params$Imax_u,
                   ddi_res[["midazolam_ketoconazole"]]$params$Igut,
                   0.94, 0.43, 0.048, 0.86, 0.00032, NA, NA),
  sotorasib = c(10.0,
                ddi_res[["sotorasib_perpetrator"]]$params$Imax_u,
                ddi_res[["sotorasib_perpetrator"]]$params$Igut,
                NA, NA, 0.04, 3.5, 0.00032, 8.1, 1.5),
  source = c("Obach 2007 / FDA NDA review", "Computed from Cmax*fu/MW",
             "Computed from dose/250mL/MW", "Midazolam literature",
             "Midazolam literature", "Obach 2007 / FDA NDA review",
             "Obach 2007 / FDA NDA review", "Yang 2008",
             "FDA NDA 214665 review", "FDA NDA 214665 review")
)

write_csv(param_table, file.path(DIR_TABLES, "parameter_table.csv"))
log_msg("Saved: parameter_table.csv")

# =============================================================================
# SUMMARY LOG
# =============================================================================

log_msg("=== Table Generation Complete ===")
log_msg("Generated 6 summary tables in: ", DIR_TABLES)
