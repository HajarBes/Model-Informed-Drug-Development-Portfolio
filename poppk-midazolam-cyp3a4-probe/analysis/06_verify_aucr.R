# =============================================================================
# 06_verify_aucr.R — Independent AUCR verification with expanded weight strata
# =============================================================================
# Recomputes AUC ratios vs 70 kg reference for 7.5 mg standard dose.
# Includes: median AUCR, 90% PI, % subjects outside 0.80-1.25 heuristic band.
# Weight strata: 50, 70 (ref), 90, 100, 110 kg.
# =============================================================================

script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(script_dir) || !file.exists(file.path(script_dir, "00_setup.R"))) {
  script_dir <- file.path(getwd(), "analysis")
}
source(file.path(script_dir, "00_setup.R"))

log_msg("=== AUCR Verification (Independent Recomputation) ===")

set.seed(20240201)  # Same seed as 05_simulation_decision.R

N_SIM  <- 1000
T_SIM  <- seq(0, 24, by = 0.1)
DOSE   <- 7500  # 7.5 mg in ug

wt_strata <- c(50, 70, 90, 100, 110)

log_msg("Dose: 7.5 mg | N = ", N_SIM, " per stratum | WT strata: ",
        paste(wt_strata, collapse = ", "), " kg")

# --- Simulate each stratum ---------------------------------------------------

all_results <- list()

for (wt in wt_strata) {
  for (i in 1:N_SIM) {
    eta_ka <- rnorm(1, 0, sqrt(TRUE_OMEGA$ka))
    eta_cl <- rnorm(1, 0, sqrt(TRUE_OMEGA$cl))
    eta_vc <- rnorm(1, 0, sqrt(TRUE_OMEGA$vc))

    ind_ka <- TRUE_PARAMS$ka   * exp(eta_ka)
    ind_cl <- TRUE_PARAMS$cl_f * exp(eta_cl) * (wt / 70)^0.75
    ind_vc <- TRUE_PARAMS$vc_f * exp(eta_vc) * (wt / 70)^1.0
    ind_q  <- TRUE_PARAMS$q_f  * (wt / 70)^0.75
    ind_vp <- TRUE_PARAMS$vp_f * (wt / 70)^1.0

    pk <- simulate_subject_pk(
      dose_ug = DOSE, times = T_SIM[T_SIM > 0],
      ka = ind_ka, cl = ind_cl, v1 = ind_vc, q = ind_q, v2 = ind_vp
    )

    all_results[[length(all_results) + 1]] <- tibble(
      wt     = wt,
      sim_id = i,
      AUC    = compute_auc(pk$time, pk$CP)
    )
  }
}

sim_df <- bind_rows(all_results)

# --- Compute individual-level AUCR vs 70 kg reference -------------------------

ref_aucs <- sim_df %>% filter(wt == 70) %>% pull(AUC)
ref_median <- median(ref_aucs)

log_msg("Reference (70 kg) median AUC: ", round(ref_median, 1), " ng*h/mL")

# Per-stratum summary
aucr_summary <- sim_df %>%
  group_by(wt) %>%
  summarise(
    n              = n(),
    auc_median     = round(median(AUC), 1),
    auc_p5         = round(quantile(AUC, 0.05), 1),
    auc_p95        = round(quantile(AUC, 0.95), 1),
    aucr_median    = round(median(AUC) / ref_median, 3),
    aucr_p5        = round(quantile(AUC, 0.05) / ref_median, 3),
    aucr_p95       = round(quantile(AUC, 0.95) / ref_median, 3),
    pct_below_0.80 = round(100 * mean(AUC / ref_median < 0.80), 1),
    pct_above_1.25 = round(100 * mean(AUC / ref_median > 1.25), 1),
    pct_outside    = round(100 * mean(AUC / ref_median < 0.80 | AUC / ref_median > 1.25), 1),
    .groups = "drop"
  ) %>%
  mutate(
    stratum = case_when(
      wt == 50  ~ "50 kg (low)",
      wt == 70  ~ "70 kg (reference)",
      wt == 90  ~ "90 kg",
      wt == 100 ~ "100 kg",
      wt == 110 ~ "110 kg (high)"
    ),
    # Theoretical AUCR from allometric formula: (70/WT)^0.75
    aucr_theoretical = round((70 / wt)^0.75, 3)
  ) %>%
  select(stratum, wt, n, auc_median, auc_p5, auc_p95,
         aucr_median, aucr_theoretical, aucr_p5, aucr_p95,
         pct_below_0.80, pct_above_1.25, pct_outside)

# --- Report to log ------------------------------------------------------------

log_msg("")
log_msg("AUCR Summary (7.5 mg, vs 70 kg reference):")
log_msg("The 0.80-1.25 band is used as a clinical relevance heuristic,")
log_msg("NOT as formal bioequivalence criteria.")
log_msg("")
log_msg(sprintf("%-20s  Median AUCR  Theoretical  90%% PI            %%Outside",
                "Stratum"))
log_msg(paste(rep("-", 80), collapse = ""))

for (i in seq_len(nrow(aucr_summary))) {
  r <- aucr_summary[i, ]
  log_msg(sprintf("%-20s  %6.3f       %6.3f       [%5.3f, %5.3f]     %4.1f%%",
                  r$stratum, r$aucr_median, r$aucr_theoretical,
                  r$aucr_p5, r$aucr_p95, r$pct_outside))
}

log_msg("")

# --- Classify strata ----------------------------------------------------------

within_band  <- aucr_summary %>% filter(aucr_median >= 0.80, aucr_median <= 1.25)
outside_band <- aucr_summary %>% filter(aucr_median < 0.80 | aucr_median > 1.25)

if (nrow(outside_band) > 0) {
  log_msg("Strata with median AUCR OUTSIDE 0.80-1.25 band:")
  for (i in seq_len(nrow(outside_band))) {
    r <- outside_band[i, ]
    direction <- ifelse(r$aucr_median > 1.25, "above 1.25", "below 0.80")
    log_msg(sprintf("  %s: median AUCR = %.3f (%s)", r$stratum, r$aucr_median, direction))
  }
} else {
  log_msg("All strata have median AUCR within the 0.80-1.25 band.")
}

if (nrow(within_band) > 0) {
  log_msg("Strata with median AUCR WITHIN 0.80-1.25 band:")
  for (i in seq_len(nrow(within_band))) {
    r <- within_band[i, ]
    log_msg(sprintf("  %s: median AUCR = %.3f", r$stratum, r$aucr_median))
  }
}

log_msg("")
log_msg("NOTE: Wide 90% PIs reflect ~30% BSV in CL. Individual variability")
log_msg("substantially exceeds the weight effect for ALL strata.")

# --- Save individual-level table -----------------------------------------------

write_csv(aucr_summary, file.path(DIR_TABLES, "aucr_weight_summary_verified.csv"))
log_msg("")
log_msg("Saved: aucr_weight_summary_verified.csv")

# =============================================================================
# METRIC 2: STRATUM-LEVEL MEDIAN SHIFT (isolates weight effect)
# =============================================================================
# AUCR_shift = median(AUC_W) / median(AUC_70)
# This is a single number per stratum, free of BSV mixing.
# 90% CI via bootstrap of the ratio of medians.

log_msg("")
log_msg("--- Stratum-Level Median Shift (bootstrap CI) ---")

N_BOOT <- 2000
ref_auc_vec <- sim_df %>% filter(wt == 70) %>% pull(AUC)

shift_results <- list()
for (w in wt_strata) {
  w_auc_vec <- sim_df %>% filter(wt == w) %>% pull(AUC)

  # Point estimate
  shift_point <- median(w_auc_vec) / median(ref_auc_vec)

  # Bootstrap CI for the ratio of medians
  boot_ratios <- numeric(N_BOOT)
  for (b in 1:N_BOOT) {
    boot_w   <- sample(w_auc_vec, length(w_auc_vec), replace = TRUE)
    boot_ref <- sample(ref_auc_vec, length(ref_auc_vec), replace = TRUE)
    boot_ratios[b] <- median(boot_w) / median(boot_ref)
  }

  shift_results[[length(shift_results) + 1]] <- tibble(
    stratum          = case_when(
      w == 50  ~ "50 kg (low)",
      w == 70  ~ "70 kg (reference)",
      w == 90  ~ "90 kg",
      w == 100 ~ "100 kg",
      w == 110 ~ "110 kg (high)"
    ),
    wt               = w,
    median_auc       = round(median(w_auc_vec), 1),
    aucr_shift       = round(shift_point, 3),
    aucr_theoretical = round((70 / w)^0.75, 3),
    aucr_boot_lo     = round(quantile(boot_ratios, 0.05), 3),
    aucr_boot_hi     = round(quantile(boot_ratios, 0.95), 3),
    within_band      = shift_point >= 0.80 & shift_point <= 1.25
  )
}

shift_df <- bind_rows(shift_results)

log_msg("")
log_msg("AUCR Stratum Shift (median ratio, 90% bootstrap CI):")
log_msg("  Denominator: median AUC across N=1000 subjects at 70 kg")
log_msg("")
log_msg(sprintf("%-20s  Shift   Theoretical  90%% CI             Within 0.80-1.25?",
                "Stratum"))
log_msg(paste(rep("-", 80), collapse = ""))

for (i in seq_len(nrow(shift_df))) {
  r <- shift_df[i, ]
  flag <- ifelse(r$within_band, "YES", "NO")
  log_msg(sprintf("%-20s  %5.3f   %5.3f        [%5.3f, %5.3f]     %s",
                  r$stratum, r$aucr_shift, r$aucr_theoretical,
                  r$aucr_boot_lo, r$aucr_boot_hi, flag))
}

write_csv(shift_df, file.path(DIR_TABLES, "aucr_stratum_shift_summary.csv"))
log_msg("")
log_msg("Saved: aucr_stratum_shift_summary.csv")
log_msg("")
log_msg("INTERPRETATION:")
log_msg("  The stratum shift isolates the weight effect on TYPICAL exposure.")
log_msg("  The individual-level '% outside' in the verified table reflects")
log_msg("  BSV + weight combined. At 70 kg (reference), 45% of subjects fall")
log_msg("  outside 0.80-1.25 due to BSV in CL alone — this is pharmacokinetic")
log_msg("  variability, not a weight effect.")
log_msg("  See outputs/tables/aucr_definition.md for full documentation.")
log_msg("")
log_msg("=== AUCR Verification Complete ===")
