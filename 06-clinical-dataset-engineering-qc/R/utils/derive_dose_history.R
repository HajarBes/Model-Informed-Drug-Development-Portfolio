# ==============================================================================
# Utility: derive_dose_history.R - Cumulative Dose and Dose History
# ==============================================================================

#' Derive cumulative dose per subject
#'
#' @param data Data frame with ID, EVID, AMT columns (sorted by ID, TIME)
#' @return Data frame with CUMDOSE column added
derive_cumulative_dose <- function(data) {
  data %>%
    group_by(ID) %>%
    mutate(
      CUMDOSE = cumsum(if_else(EVID == 1 & !is.na(AMT), AMT, 0))
    ) %>%
    ungroup()
}

#' Derive number of prior doses
#'
#' @param data Data frame with ID, EVID columns (sorted by ID, TIME)
#' @return Data frame with NDOSE column added
derive_ndose <- function(data) {
  data %>%
    group_by(ID) %>%
    mutate(
      NDOSE = cumsum(EVID == 1)
    ) %>%
    ungroup()
}
