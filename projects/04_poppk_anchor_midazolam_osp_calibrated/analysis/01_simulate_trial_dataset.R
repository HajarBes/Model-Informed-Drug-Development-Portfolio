# =============================================================================
# 01_simulate_trial_dataset.R — Generate publication-calibrated simulated trial
# =============================================================================
# Produces a NONMEM-format analysis dataset for oral midazolam 7.5 mg.
# "True" PK parameters are calibrated to reproduce profiles consistent with
# OSP Midazolam PBPK evaluation studies (Hohmann 2015, Link 2008, Gorski 2003).
#
# IMPORTANT: This is simulated data. It does not represent real clinical data.
# =============================================================================

script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(script_dir) || !file.exists(file.path(script_dir, "00_setup.R"))) {
  script_dir <- file.path(getwd(), "analysis")
}
source(file.path(script_dir, "00_setup.R"))

log_msg("=== Starting Trial Dataset Simulation ===")

set.seed(20240101)

# =============================================================================
# STUDY DESIGN
# =============================================================================

N_TOTAL    <- 120         # Total subjects
N_RICH     <- 40          # Rich PK sampling arm
N_SPARSE   <- 80          # Sparse clinical sampling arm
DOSE_MG    <- 7.5         # Oral midazolam dose
DOSE_UG    <- DOSE_MG * 1000

# Sampling schedules (hours post-dose)
TIMES_RICH   <- c(0.25, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12)
TIMES_SPARSE <- c(0.5, 2, 4, 8)

log_msg("Design: N=", N_TOTAL, " (", N_RICH, " rich + ", N_SPARSE, " sparse)")
log_msg("Dose: ", DOSE_MG, " mg oral midazolam")
log_msg("Rich schedule: ", paste(TIMES_RICH, collapse=", "), " h")
log_msg("Sparse schedule: ", paste(TIMES_SPARSE, collapse=", "), " h")

# =============================================================================
# DEMOGRAPHICS
# =============================================================================

# Generate realistic covariate distributions
demographics <- tibble(
  ID  = 1:N_TOTAL,
  ARM = c(rep("RICH", N_RICH), rep("SPARSE", N_SPARSE)),
  SEX = sample(c(0L, 1L), N_TOTAL, replace = TRUE, prob = c(0.45, 0.55)),
  AGE = round(runif(N_TOTAL, 22, 72)),
  WT  = round(ifelse(
    SEX == 1,
    rlnorm(N_TOTAL, log(78), 0.18),  # Male: median ~78 kg
    rlnorm(N_TOTAL, log(66), 0.20)   # Female: median ~66 kg
  ), 1)
)

# Constrain to physiological range
demographics <- demographics %>%
  mutate(WT = pmin(pmax(WT, 40), 150))

log_msg("Demographics: median WT=", median(demographics$WT), " kg, ",
        sum(demographics$SEX == 1), " male, median AGE=", median(demographics$AGE))

# =============================================================================
# SIMULATE INDIVIDUAL PK PROFILES
# =============================================================================

log_msg("Simulating individual PK profiles...")

all_records <- list()

for (i in seq_len(N_TOTAL)) {
  subj <- demographics[i, ]

  # Individual parameters: true population + BSV + allometric scaling
  eta_ka <- rnorm(1, 0, sqrt(TRUE_OMEGA$ka))
  eta_cl <- rnorm(1, 0, sqrt(TRUE_OMEGA$cl))
  eta_vc <- rnorm(1, 0, sqrt(TRUE_OMEGA$vc))

  wt_ratio <- subj$WT / TRUE_PARAMS$wt_ref

  ind_ka <- TRUE_PARAMS$ka   * exp(eta_ka)
  ind_cl <- TRUE_PARAMS$cl_f * exp(eta_cl) * wt_ratio^0.75
  ind_vc <- TRUE_PARAMS$vc_f * exp(eta_vc) * wt_ratio^1.0
  ind_q  <- TRUE_PARAMS$q_f  * wt_ratio^0.75
  ind_vp <- TRUE_PARAMS$vp_f * wt_ratio^1.0

  # Select sampling schedule
  obs_times <- if (subj$ARM == "RICH") TIMES_RICH else TIMES_SPARSE

  # Solve ODE
  pk_profile <- simulate_subject_pk(
    dose_ug = DOSE_UG, times = obs_times,
    ka = ind_ka, cl = ind_cl, v1 = ind_vc, q = ind_q, v2 = ind_vp
  )

  # Add residual error: combined proportional + additive
  pk_profile$DV <- pk_profile$CP * (1 + rnorm(nrow(pk_profile), 0, TRUE_SIGMA$prop)) +
                   rnorm(nrow(pk_profile), 0, TRUE_SIGMA$add)
  pk_profile$DV <- pmax(pk_profile$DV, 0)  # Concentrations cannot be negative

  # Create NONMEM-format records
  # Dosing record
  dose_rec <- tibble(
    ID   = subj$ID,
    TIME = 0,
    AMT  = DOSE_UG,
    EVID = 1L,
    CMT  = 1L,
    MDV  = 1L,
    DV   = NA_real_,
    IPRED = NA_real_,
    WT   = subj$WT,
    SEX  = subj$SEX,
    AGE  = subj$AGE,
    ARM  = subj$ARM,
    BLQ  = 0L
  )

  # Observation records
  obs_recs <- tibble(
    ID    = subj$ID,
    TIME  = pk_profile$time,
    AMT   = 0,
    EVID  = 0L,
    CMT   = 2L,
    MDV   = 0L,
    DV    = round(pk_profile$DV, 3),
    IPRED = round(pk_profile$CP, 3),
    WT    = subj$WT,
    SEX   = subj$SEX,
    AGE   = subj$AGE,
    ARM   = subj$ARM,
    BLQ   = 0L
  )

  # Handle BLQ: M1 method (set MDV=1 for observations below LLOQ)
  obs_recs <- obs_recs %>%
    mutate(
      BLQ = as.integer(DV < LLOQ & EVID == 0),
      MDV = ifelse(BLQ == 1, 1L, MDV),
      DV  = ifelse(BLQ == 1, NA_real_, DV)
    )

  all_records[[i]] <- bind_rows(dose_rec, obs_recs)
}

dataset <- bind_rows(all_records) %>%
  arrange(ID, TIME, desc(EVID))

# =============================================================================
# DATASET SUMMARY
# =============================================================================

obs_data <- dataset %>% filter(EVID == 0)
n_obs   <- nrow(obs_data)
n_blq   <- sum(dataset$BLQ == 1, na.rm = TRUE)
n_above <- n_obs - n_blq

log_msg("Dataset generated:")
log_msg("  Total records: ", nrow(dataset))
log_msg("  Dosing records: ", sum(dataset$EVID == 1))
log_msg("  Observation records: ", n_obs)
log_msg("  Above LLOQ: ", n_above, " (", round(n_above/n_obs*100, 1), "%)")
log_msg("  BLQ (M1 excluded): ", n_blq, " (", round(n_blq/n_obs*100, 1), "%)")

conc_summary <- obs_data %>%
  filter(BLQ == 0) %>%
  summarise(
    median_Cmax = median(DV, na.rm = TRUE),
    range_DV    = paste0("[", round(min(DV, na.rm=TRUE), 2), ", ",
                         round(max(DV, na.rm=TRUE), 2), "]")
  )
log_msg("  Concentration range: ", conc_summary$range_DV, " ng/mL")

# =============================================================================
# SAVE
# =============================================================================

write_csv(dataset, file.path(DIR_SIM, "analysis_dataset.csv"))
write_csv(demographics, file.path(DIR_SIM, "demographics.csv"))
log_msg("Saved: data/simulated/analysis_dataset.csv")
log_msg("Saved: data/simulated/demographics.csv")

# Quick concentration-time plot for QC
qc_plot <- obs_data %>%
  filter(BLQ == 0) %>%
  ggplot(aes(x = TIME, y = DV, group = ID)) +
  geom_line(alpha = 0.15, color = COL_PRIMARY) +
  geom_point(alpha = 0.2, size = 0.8, color = COL_PRIMARY) +
  scale_y_log10(limits = c(0.5, 500)) +
  labs(
    title = "Simulated Midazolam PK Profiles (QC)",
    subtitle = sprintf("N=%d | Oral %g mg | Rich + sparse sampling", N_TOTAL, DOSE_MG),
    x = "Time (h)", y = "Concentration (ng/mL)"
  )

ggsave(file.path(DIR_FIGURES, "qc_spaghetti_log.png"), qc_plot,
       width = 8, height = 5, dpi = 300, bg = "white")
log_msg("Saved: figures/qc_spaghetti_log.png")

log_msg("=== Trial Dataset Simulation Complete ===")
