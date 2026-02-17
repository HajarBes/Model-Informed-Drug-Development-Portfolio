# =============================================================================
# 01_compute_static_ddi.R — Compute static DDI screening ratios and AUCR
# Runs on ALL cases found in data/ (inputs_*_case.csv)
# =============================================================================

# Source setup — works from Rscript or interactive
script_dir <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) NULL
)
if (is.null(script_dir) || !file.exists(file.path(script_dir, "00_setup.R"))) {
  script_dir <- file.path(getwd(), "analysis")
}
source(file.path(script_dir, "00_setup.R"))

log_msg("=== Starting Static DDI Computation ===")

# --- Discover input files ----------------------------------------------------
input_files <- list.files(DIR_DATA, pattern = "^inputs_.*_case\\.csv$",
                          full.names = TRUE)
log_msg("Found ", length(input_files), " case file(s): ",
        paste(basename(input_files), collapse = ", "))

# --- Process each case -------------------------------------------------------
all_results <- list()

for (input_file in input_files) {

  log_msg("--- Processing: ", basename(input_file), " ---")
  params <- parse_inputs(input_file)
  case   <- params$case_name
  log_msg("Case: ", case)
  log_msg("Perpetrator: ", params$perpetrator)
  log_msg("[I]max,u = ", round(params$Imax_u, 4), " uM")
  log_msg("[I]gut   = ", round(params$Igut, 1), " uM")

  # === CYP INHIBITION SCREENING =============================================
  cyp_results <- map_dfr(params$cyp, function(cyp) {
    enz <- cyp$enzyme

    R1     <- compute_R1(params$Imax_u, cyp$ki)
    R1_gut <- compute_R1_gut(params$Igut, cyp$ki)

    tdi_factor <- compute_TDI_factor(
      cyp$kinact, params$Imax_u, cyp$KI, params$kdeg_hepatic
    )

    tibble(
      case          = case,
      enzyme        = enz,
      Ki_uM         = cyp$ki,
      Imax_u_uM     = params$Imax_u,
      Igut_uM       = params$Igut,
      R1            = R1,
      R1_flag       = !is.na(R1) & R1 >= 1.02,
      R1_gut        = R1_gut,
      R1_gut_flag   = !is.na(R1_gut) & R1_gut >= 11,
      has_TDI       = !is.na(cyp$kinact),
      kinact        = cyp$kinact,
      KI_uM         = cyp$KI,
      TDI_factor    = tdi_factor,
      kdeg          = params$kdeg_hepatic
    )
  })

  log_msg("CYP screening complete: ", nrow(cyp_results), " enzymes assessed")
  for (i in seq_len(nrow(cyp_results))) {
    r <- cyp_results[i, ]
    log_msg("  ", r$enzyme, ": R1=", round(r$R1, 3),
            " (flag=", r$R1_flag, ")",
            " R1,gut=", round(r$R1_gut, 1),
            " (flag=", r$R1_gut_flag, ")")
  }

  # === INDUCTION SCREENING ==================================================
  induction_result <- tibble(
    case       = case,
    enzyme     = "CYP3A4",
    has_data   = params$induction$available & !is.na(params$induction$emax),
    Emax       = params$induction$emax,
    EC50_uM    = params$induction$ec50,
    d_scaling  = params$induction$d,
    Imax_u_uM  = params$Imax_u,
    R3         = compute_R3(params$induction$emax, params$Imax_u,
                            params$induction$ec50, params$induction$d),
    R3_flag    = NA
  )
  induction_result$R3_flag <- !is.na(induction_result$R3) &
                                induction_result$R3 <= 0.8

  if (induction_result$has_data) {
    log_msg("Induction: R3=", round(induction_result$R3, 4),
            " (flag=", induction_result$R3_flag, ")")
  } else {
    log_msg("Induction: No data available")
  }

  # === AUCR PREDICTION (CYP3A4 only, if victim fm available) ================
  aucr_result <- tibble(
    case          = case,
    enzyme        = "CYP3A4",
    fm            = params$fm_cyp3a4,
    fg            = params$fg_cyp3a4,
    R1            = cyp_results %>% filter(enzyme == "CYP3A4") %>% pull(R1),
    TDI_factor    = cyp_results %>% filter(enzyme == "CYP3A4") %>% pull(TDI_factor),
    R1_gut        = cyp_results %>% filter(enzyme == "CYP3A4") %>% pull(R1_gut),
    R3            = induction_result$R3,
    AUCR_hepatic  = NA_real_,
    AUCR_gut      = NA_real_,
    AUCR_total    = NA_real_,
    classification = NA_character_
  )

  if (!is.na(params$fm_cyp3a4)) {
    aucr_result$AUCR_hepatic <- compute_AUCR_hepatic(
      params$fm_cyp3a4,
      aucr_result$R1,
      aucr_result$TDI_factor
    )

    fg_val <- ifelse(is.na(params$fg_cyp3a4), 0, params$fg_cyp3a4)
    aucr_result$AUCR_gut <- compute_AUCR_gut(fg_val, aucr_result$R1_gut)

    aucr_result$AUCR_total <- aucr_result$AUCR_hepatic * aucr_result$AUCR_gut
    aucr_result$classification <- classify_ddi(aucr_result$AUCR_total)

    log_msg("AUCR prediction (CYP3A4):")
    log_msg("  Hepatic: ", round(aucr_result$AUCR_hepatic, 2))
    log_msg("  Gut:     ", round(aucr_result$AUCR_gut, 2))
    log_msg("  Total:   ", round(aucr_result$AUCR_total, 2))
    log_msg("  Class:   ", aucr_result$classification)
    if (!is.na(params$observed_aucr)) {
      log_msg("  Observed AUCR: ", params$observed_aucr,
              " (pred/obs = ", round(aucr_result$AUCR_total / params$observed_aucr, 2), ")")
    }
  } else {
    log_msg("AUCR: No victim fm available — screening ratios only")
  }

  # === TRANSPORTER SCREENING ================================================
  trans_results <- map_dfr(names(params$transporters), function(tn) {
    ic50 <- params$transporters[[tn]]

    # Determine relevant concentration and threshold
    if (tn %in% c("P-gp", "BCRP")) {
      I_relevant <- params$Igut
      threshold  <- 10
      conc_type  <- "[I]gut"
    } else {
      I_relevant <- params$Imax_u
      threshold  <- 0.1
      conc_type  <- "[I]max,u"
    }

    result <- screen_transporter(I_relevant, ic50, threshold)

    tibble(
      case         = case,
      transporter  = tn,
      IC50_uM      = ic50,
      I_relevant   = I_relevant,
      conc_type    = conc_type,
      ratio        = result$ratio,
      threshold    = threshold,
      flag         = result$flag
    )
  })

  flagged_trans <- trans_results %>% filter(flag == TRUE)
  log_msg("Transporter screening: ", nrow(trans_results), " assessed, ",
          nrow(flagged_trans), " flagged")
  if (nrow(flagged_trans) > 0) {
    for (i in seq_len(nrow(flagged_trans))) {
      r <- flagged_trans[i, ]
      log_msg("  FLAGGED: ", r$transporter, " ratio=", round(r$ratio, 2),
              " (threshold=", r$threshold, ")")
    }
  }

  # === RISK SUMMARY ==========================================================
  # Determine overall recommendation
  cyp_flagged <- any(cyp_results$R1_flag, na.rm = TRUE) |
                 any(cyp_results$R1_gut_flag, na.rm = TRUE)
  tdi_present <- any(cyp_results$has_TDI & cyp_results$TDI_factor < 0.5,
                     na.rm = TRUE)
  induction_flagged <- isTRUE(induction_result$R3_flag)
  trans_flagged_any <- any(trans_results$flag, na.rm = TRUE)

  # Build recommendation
  recommendations <- c()
  if (cyp_flagged)       recommendations <- c(recommendations, "Clinical DDI study or PBPK for CYP inhibition")
  if (tdi_present)       recommendations <- c(recommendations, "Evaluate TDI contribution via PBPK")
  if (induction_flagged) recommendations <- c(recommendations, "Clinical DDI study or PBPK for CYP induction")
  if (cyp_flagged & induction_flagged) {
    recommendations <- c(recommendations, "PBPK recommended to resolve net inhibition vs. induction")
  }
  if (trans_flagged_any) recommendations <- c(recommendations, "Clinical DDI study for flagged transporters")
  if (length(recommendations) == 0) recommendations <- "Low DDI risk — no further studies needed"

  risk_summary <- tibble(
    case           = case,
    perpetrator    = params$perpetrator,
    victim         = ifelse(is.na(params$victim) || params$victim == "NA",
                           "General assessment", params$victim),
    cyp_flagged    = cyp_flagged,
    tdi_present    = tdi_present,
    induction_flag = induction_flagged,
    trans_flagged  = trans_flagged_any,
    AUCR_pred      = aucr_result$AUCR_total,
    classification = aucr_result$classification,
    recommendation = paste(recommendations, collapse = "; ")
  )

  log_msg("RECOMMENDATION: ", risk_summary$recommendation)

  # === STORE RESULTS ==========================================================
  all_results[[case]] <- list(
    params     = params,
    cyp        = cyp_results,
    induction  = induction_result,
    aucr       = aucr_result,
    transporter = trans_results,
    summary    = risk_summary
  )

  # === SAVE PER-CASE OUTPUTS =================================================
  case_tag <- gsub("[^a-zA-Z0-9]", "_", tolower(case))

  write_csv(cyp_results,
            file.path(DIR_TABLES, paste0("cyp_screening_", case_tag, ".csv")))
  write_csv(trans_results,
            file.path(DIR_TABLES, paste0("transporter_screening_", case_tag, ".csv")))
  write_csv(risk_summary,
            file.path(DIR_TABLES, paste0("risk_summary_", case_tag, ".csv")))

  log_msg("Results saved for case: ", case)
}

# === SAVE COMBINED RESULTS ===================================================
saveRDS(all_results, file.path(DIR_PROCESSED, "ddi_results.rds"))
log_msg("All results saved to data/processed/ddi_results.rds")
log_msg("=== Static DDI Computation Complete ===")
