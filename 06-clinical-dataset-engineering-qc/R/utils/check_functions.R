# ==============================================================================
# Utility: check_functions.R - Reusable QC Check Functions
# ==============================================================================

#' Run all QC checks on NONMEM dataset
#'
#' @param data NONMEM-ready data frame
#' @param config Study configuration list
#' @return Data frame with check results (check_name, status, n_issues, details)
run_all_qc_checks <- function(data, config) {

  results <- bind_rows(
    check_missing_values(data),
    check_duplicate_records(data),
    check_time_consistency(data),
    check_dose_records(data),
    check_blq_handling(data, config),
    check_covariate_ranges(data),
    check_subject_completeness(data)
  )

  return(results)
}

#' Check for unexpected missing values
check_missing_values <- function(data) {
  critical_vars <- c("ID", "TIME", "EVID", "MDV", "CMT")
  n_missing <- sapply(critical_vars, function(v) sum(is.na(data[[v]])))

  tibble(
    check_name = paste0("Missing: ", critical_vars),
    status     = if_else(n_missing == 0, "PASS", "FAIL"),
    n_issues   = n_missing,
    details    = if_else(n_missing == 0, "No missing values",
                         paste(n_missing, "records with missing", critical_vars))
  )
}

#' Check for duplicate records
check_duplicate_records <- function(data) {
  n_dup <- data %>%
    group_by(ID, TIME, EVID, CMT) %>%
    filter(n() > 1) %>%
    nrow()

  tibble(
    check_name = "Duplicate records (ID+TIME+EVID+CMT)",
    status     = if_else(n_dup == 0, "PASS", "WARNING"),
    n_issues   = n_dup,
    details    = if_else(n_dup == 0, "No duplicates found",
                         paste(n_dup, "duplicate record combinations"))
  )
}

#' Check time consistency (pre-dose acceptable, flag data errors)
check_time_consistency <- function(data) {
  # Pre-dose samples (TIME slightly < 0) are normal in clinical PK
  predose <- sum(data$TIME < 0 & data$TIME >= -24, na.rm = TRUE)
  neg_error <- sum(data$TIME < -24, na.rm = TRUE)
  tad_gt_time <- sum(data$TAD > data$TIME & data$TIME >= 0, na.rm = TRUE)

  bind_rows(
    tibble(check_name = "Pre-dose samples (TIME < 0, >= -24h)",
           status = if_else(predose == 0, "PASS", "INFO"),
           n_issues = predose,
           details = if_else(predose == 0, "No pre-dose samples",
                             sprintf("%d pre-dose samples (expected, baseline PK)", predose))),
    tibble(check_name = "Negative TIME errors (< -24h)",
           status = if_else(neg_error == 0, "PASS", "FAIL"),
           n_issues = neg_error,
           details = if_else(neg_error == 0, "No TIME errors",
                             paste(neg_error, "records with TIME < -24h (data error)"))),
    tibble(check_name = "TAD > TIME (post-dose only)",
           status = if_else(tad_gt_time == 0, "PASS", "FAIL"),
           n_issues = tad_gt_time,
           details = if_else(tad_gt_time == 0, "All TAD <= TIME for post-dose records",
                             paste(tad_gt_time, "records where TAD exceeds TIME")))
  )
}

#' Check dosing records for consistency
check_dose_records <- function(data) {
  dose_data <- data %>% filter(EVID == 1)
  zero_amt <- sum(dose_data$AMT == 0 | is.na(dose_data$AMT))
  dose_has_dv <- sum(!is.na(dose_data$DV) & dose_data$DV != 0)

  bind_rows(
    tibble(check_name = "Dose records with zero/missing AMT",
           status = if_else(zero_amt == 0, "PASS", "FAIL"),
           n_issues = zero_amt,
           details = if_else(zero_amt == 0, "All dose records have AMT > 0",
                             paste(zero_amt, "dose records with AMT = 0 or NA"))),
    tibble(check_name = "Dose records with non-zero DV",
           status = if_else(dose_has_dv == 0, "PASS", "WARNING"),
           n_issues = dose_has_dv,
           details = if_else(dose_has_dv == 0, "No DV values on dose records",
                             paste(dose_has_dv, "dose records have DV values")))
  )
}

#' Check BLQ handling
check_blq_handling <- function(data, config) {
  obs <- data %>% filter(EVID == 0)
  n_blq <- sum(obs$BLQ == 1, na.rm = TRUE)
  pct_blq <- round(100 * n_blq / nrow(obs), 1)

  # Check BLQ values match rule
  blq_obs <- obs %>% filter(BLQ == 1)
  if (nrow(blq_obs) > 0 && config$blq_rule == "LLOQ/2") {
    expected_dv <- config$lloq / 2
    wrong_dv <- sum(abs(blq_obs$DV - expected_dv) > 1e-6, na.rm = TRUE)
  } else {
    wrong_dv <- 0
  }

  bind_rows(
    tibble(check_name = "BLQ summary",
           status = "INFO",
           n_issues = n_blq,
           details = sprintf("%d BLQ observations (%.1f%% of total)", n_blq, pct_blq)),
    tibble(check_name = "BLQ rule applied correctly",
           status = if_else(wrong_dv == 0, "PASS", "FAIL"),
           n_issues = wrong_dv,
           details = if_else(wrong_dv == 0,
                             sprintf("All BLQ values set to %s", config$blq_rule),
                             paste(wrong_dv, "BLQ records with incorrect DV")))
  )
}

#' Check covariate ranges
check_covariate_ranges <- function(data) {
  subj <- data %>% distinct(ID, .keep_all = TRUE)

  expected_ranges <- list(
    AGE = c(18, 100),
    WT  = c(30, 200),
    HT  = c(120, 220),
    BMI = c(12, 60)
  )

  results <- map_dfr(names(expected_ranges), function(var) {
    if (!var %in% names(subj)) {
      return(tibble(check_name = paste0("Covariate range: ", var),
                    status = "SKIP", n_issues = 0,
                    details = paste(var, "not in dataset")))
    }
    vals <- subj[[var]]
    rng <- expected_ranges[[var]]
    n_na  <- sum(is.na(vals))

    # Handle all-NA columns
    if (all(is.na(vals))) {
      return(tibble(check_name = paste0("Covariate range: ", var),
                    status = "INFO",
                    n_issues = n_na,
                    details = sprintf("All %d values missing (not available in source)", n_na)))
    }

    n_oor <- sum(vals < rng[1] | vals > rng[2], na.rm = TRUE)

    tibble(check_name = paste0("Covariate range: ", var),
           status = if_else(n_oor == 0 & n_na == 0, "PASS",
                            if_else(n_oor > 0, "WARNING", "INFO")),
           n_issues = n_oor + n_na,
           details = sprintf("Range [%.1f-%.1f], %d out-of-range, %d missing",
                             min(vals, na.rm = TRUE), max(vals, na.rm = TRUE),
                             n_oor, n_na))
  })

  return(results)
}

#' Check subject completeness
check_subject_completeness <- function(data) {
  subj_summary <- data %>%
    group_by(ID) %>%
    summarise(
      n_doses = sum(EVID == 1),
      n_obs   = sum(EVID == 0),
      .groups = "drop"
    )

  no_doses <- sum(subj_summary$n_doses == 0)
  no_obs   <- sum(subj_summary$n_obs == 0)

  bind_rows(
    tibble(check_name = "Subjects with no dosing records",
           status = if_else(no_doses == 0, "PASS", "WARNING"),
           n_issues = no_doses,
           details = if_else(no_doses == 0, "All subjects have dose records",
                             paste(no_doses, "subjects without doses"))),
    tibble(check_name = "Subjects with no observations",
           status = if_else(no_obs == 0, "PASS", "WARNING"),
           n_issues = no_obs,
           details = if_else(no_obs == 0, "All subjects have observations",
                             paste(no_obs, "subjects without observations")))
  )
}

#' Perform subject-level spot audit
#'
#' @param data NONMEM dataset
#' @param n_subjects Number of subjects to randomly sample
#' @param seed Random seed for reproducibility
#' @return List with sampled subject IDs and their data
spot_audit <- function(data, n_subjects = 5, seed = 42) {
  set.seed(seed)
  sampled_ids <- sample(unique(data$ID), min(n_subjects, n_distinct(data$ID)))

  audit_data <- data %>%
    filter(ID %in% sampled_ids) %>%
    arrange(ID, TIME) %>%
    select(ID, TIME, TAD, AMT, DV, MDV, EVID, CMT, BLQ)

  cat(sprintf("Spot audit: %d subjects sampled (IDs: %s)\n",
              length(sampled_ids), paste(sampled_ids, collapse = ", ")))

  return(list(
    subject_ids = sampled_ids,
    data = audit_data,
    seed = seed,
    n_checked = length(sampled_ids)
  ))
}
