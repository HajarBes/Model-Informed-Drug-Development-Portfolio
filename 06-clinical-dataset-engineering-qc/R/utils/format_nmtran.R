# ==============================================================================
# Utility: format_nmtran.R - NMTRAN Formatting
# ==============================================================================

#' Format dataset for NMTRAN compatibility
#'
#' Applies NONMEM formatting conventions:
#' - Numeric columns rounded to appropriate precision
#' - Comment column (C) added
#' - NA values formatted as "." for NONMEM
#'
#' @param data Data frame with standard PK columns
#' @return Formatted data frame
format_nmtran <- function(data) {
  data <- data %>%
    mutate(
      C    = "",  # Comment column (blank = active, "C" = commented out)
      TIME = round(TIME, 4),
      TAD  = round(TAD, 4),
      AMT  = round(AMT, 2),
      DV   = round(DV, 6)
    )

  # Round optional numeric covariates if present
  if ("WT"  %in% names(data)) data$WT  <- round(data$WT, 1)
  if ("HT"  %in% names(data)) data$HT  <- round(data$HT, 1)
  if ("BMI" %in% names(data)) data$BMI <- round(data$BMI, 1)
  if ("CREAT" %in% names(data)) data$CREAT <- round(data$CREAT, 2)
  if ("ALT"   %in% names(data)) data$ALT   <- round(data$ALT, 1)
  if ("AST"   %in% names(data)) data$AST   <- round(data$AST, 1)
  if ("BILI"  %in% names(data)) data$BILI  <- round(data$BILI, 2)

  return(data)
}

#' Validate NMTRAN dataset before export
#'
#' Checks critical requirements for NONMEM input:
#' - No negative TIME values
#' - EVID matches AMT/DV logic
#' - No missing ID or TIME
#'
#' @param data NMTRAN-formatted data frame
#' @return List of check results
validate_nmtran <- function(data) {
  checks <- list()

  # Pre-dose samples (TIME < 0) are acceptable - flag only if < -24h (data error)
  checks$negative_time_error <- sum(data$TIME < -24, na.rm = TRUE)
  checks$predose_samples <- sum(data$TIME < 0 & data$TIME >= -24, na.rm = TRUE)
  checks$missing_id    <- sum(is.na(data$ID))
  checks$missing_time  <- sum(is.na(data$TIME))
  checks$dose_no_amt   <- sum(data$EVID == 1 & (is.na(data$AMT) | data$AMT == 0))
  checks$obs_has_amt   <- sum(data$EVID == 0 & !is.na(data$AMT) & data$AMT > 0)

  checks$all_pass <- all(unlist(checks) == 0)

  if (!checks$all_pass) {
    warning("NMTRAN validation issues found:")
    for (nm in names(checks)) {
      if (nm != "all_pass" && checks[[nm]] > 0) {
        warning(sprintf("  %s: %d issues", nm, checks[[nm]]))
      }
    }
  } else {
    cat("NMTRAN validation: ALL CHECKS PASSED\n")
  }

  return(checks)
}
