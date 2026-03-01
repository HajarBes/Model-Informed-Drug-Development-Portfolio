# ==============================================================================
# Utility: derive_tad.R - Time After Most Recent Dose
# ==============================================================================

#' Derive time after most recent dose (TAD)
#'
#' @param datetime POSIXct datetime vector (sorted within subject)
#' @param evid Event ID vector (1 = dose, 0 = observation)
#' @param time Numeric time from first dose (hours)
#' @return Numeric vector of TAD values (hours)
derive_tad <- function(datetime, evid, time) {
  tad <- numeric(length(time))
  last_dose_time <- NA_real_

  for (i in seq_along(time)) {
    if (evid[i] == 1) {
      last_dose_time <- time[i]
      tad[i] <- 0
    } else {
      tad[i] <- if_else(!is.na(last_dose_time), time[i] - last_dose_time, NA_real_)
    }
  }
  return(tad)
}
