# ==============================================================================
# Utility: derive_blq.R - BLQ Flag and Handling
# ==============================================================================

#' Flag BLQ observations
#'
#' @param concentration Numeric concentration values
#' @param lloq Lower limit of quantification
#' @return Integer vector (1 = BLQ, 0 = above LLOQ)
derive_blq_flag <- function(concentration, lloq) {
  if_else(!is.na(concentration) & concentration < lloq, 1L, 0L)
}

#' Apply BLQ handling rule to observations
#'
#' @param data Data frame with DV and BLQ columns
#' @param lloq Lower limit of quantification
#' @param rule Character: "LLOQ/2" (set to LLOQ/2), "zero" (set to 0),
#'   "exclude" (set MDV=1)
#' @return Modified data frame
apply_blq_rule <- function(data, lloq, rule = "LLOQ/2") {
  data %>%
    mutate(
      DV = case_when(
        BLQ == 0 ~ DV,
        rule == "LLOQ/2" ~ lloq / 2,
        rule == "zero"   ~ 0,
        rule == "exclude" ~ NA_real_,
        TRUE ~ DV
      ),
      MDV = if_else(rule == "exclude" & BLQ == 1, 1, MDV)
    )
}
