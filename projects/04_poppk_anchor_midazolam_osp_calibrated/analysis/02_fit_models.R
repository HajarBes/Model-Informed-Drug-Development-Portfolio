# =============================================================================
# 02_fit_models.R — Fit 1-comp and 2-comp PopPK models via nlmixr2 SAEM
# =============================================================================
# Candidate models:
#   Model 1: 1-compartment with first-order absorption
#   Model 2: 2-compartment with first-order absorption (expected best)
#
# Both include BSV on CL, V (and Ka), combined residual error.
# Estimation: SAEM (stochastic approximation expectation-maximization).
# =============================================================================

script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(script_dir) || !file.exists(file.path(script_dir, "00_setup.R"))) {
  script_dir <- file.path(getwd(), "analysis")
}
source(file.path(script_dir, "00_setup.R"))

if (!NLMIXR2_AVAILABLE) {
  stop("nlmixr2 is required for model fitting. Install with: install.packages('nlmixr2')")
}

log_msg("=== Starting Model Fitting ===")

# Load analysis dataset
dat_raw <- read_csv(file.path(DIR_SIM, "analysis_dataset.csv"), show_col_types = FALSE)

# Prepare for nlmixr2: keep only records with valid DV or dosing events
dat <- dat_raw %>%
  filter(MDV == 0 | EVID == 1) %>%
  select(ID, TIME, AMT, EVID, CMT, DV, MDV, WT, SEX, AGE) %>%
  mutate(
    DV  = ifelse(EVID == 1, 0, DV),
    AMT = ifelse(is.na(AMT), 0, AMT)
  )

log_msg("Fitting dataset: ", n_distinct(dat$ID), " subjects, ",
        sum(dat$EVID == 0), " observations")

# =============================================================================
# MODEL 1: ONE-COMPARTMENT ORAL
# =============================================================================

one_comp <- function() {
  ini({
    lka  <- log(2.5)    # Ka (1/h)
    lcl  <- log(50)     # CL/F (L/h)
    lv   <- log(60)     # V/F (L)
    eta.ka ~ 0.4
    eta.cl ~ 0.1
    eta.v  ~ 0.05
    prop.err <- 0.2
    add.err  <- 0.5
  })
  model({
    ka <- exp(lka + eta.ka)
    cl <- exp(lcl + eta.cl)
    v  <- exp(lv  + eta.v)
    d/dt(depot)   = -ka * depot
    d/dt(central) = ka * depot - cl/v * central
    cp = central / v
    cp ~ prop(prop.err) + add(add.err)
  })
}

log_msg("Fitting Model 1: 1-compartment oral...")
t1 <- Sys.time()
fit1 <- tryCatch(
  nlmixr2(one_comp, dat, est = "saem",
          saemControl(nBurn = 200, nEm = 100, print = 50)),
  error = function(e) {
    log_msg("Model 1 fitting failed: ", conditionMessage(e))
    NULL
  }
)
t1_elapsed <- as.numeric(difftime(Sys.time(), t1, units = "mins"))

if (!is.null(fit1)) {
  log_msg("Model 1 converged in ", round(t1_elapsed, 1), " min")
  log_msg("  OFV = ", round(fit1$objf, 2))
  log_msg("  AIC = ", round(AIC(fit1), 2))
} else {
  log_msg("Model 1: fitting failed — skipping")
}

# =============================================================================
# MODEL 2: TWO-COMPARTMENT ORAL
# =============================================================================

two_comp <- function() {
  ini({
    lka  <- log(2.5)    # Ka (1/h)
    lcl  <- log(50)     # CL/F (L/h)
    lv1  <- log(45)     # Vc/F (L)
    lq   <- log(15)     # Q/F (L/h)
    lv2  <- log(55)     # Vp/F (L)
    eta.ka  ~ 0.4
    eta.cl  ~ 0.1
    eta.v1  ~ 0.05
    prop.err <- 0.2
    add.err  <- 0.5
  })
  model({
    ka <- exp(lka + eta.ka)
    cl <- exp(lcl + eta.cl)
    v1 <- exp(lv1 + eta.v1)
    q  <- exp(lq)
    v2 <- exp(lv2)
    d/dt(depot)      = -ka * depot
    d/dt(central)    = ka * depot - cl/v1 * central - q/v1 * central + q/v2 * peripheral
    d/dt(peripheral) = q/v1 * central - q/v2 * peripheral
    cp = central / v1
    cp ~ prop(prop.err) + add(add.err)
  })
}

log_msg("Fitting Model 2: 2-compartment oral...")
t2 <- Sys.time()
fit2 <- tryCatch(
  nlmixr2(two_comp, dat, est = "saem",
          saemControl(nBurn = 200, nEm = 100, print = 50)),
  error = function(e) {
    log_msg("Model 2 fitting failed: ", conditionMessage(e))
    NULL
  }
)
t2_elapsed <- as.numeric(difftime(Sys.time(), t2, units = "mins"))

if (!is.null(fit2)) {
  log_msg("Model 2 converged in ", round(t2_elapsed, 1), " min")
  log_msg("  OFV = ", round(fit2$objf, 2))
  log_msg("  AIC = ", round(AIC(fit2), 2))
} else {
  log_msg("Model 2: fitting failed — skipping")
}

# =============================================================================
# MODEL COMPARISON
# =============================================================================

comparison <- tibble(
  Model = c("1-compartment", "2-compartment"),
  Parameters = c(
    ifelse(!is.null(fit1), length(fixef(fit1)), NA),
    ifelse(!is.null(fit2), length(fixef(fit2)), NA)
  ),
  OFV = c(
    ifelse(!is.null(fit1), round(fit1$objf, 2), NA),
    ifelse(!is.null(fit2), round(fit2$objf, 2), NA)
  ),
  AIC = c(
    ifelse(!is.null(fit1), round(AIC(fit1), 2), NA),
    ifelse(!is.null(fit2), round(AIC(fit2), 2), NA)
  ),
  Runtime_min = c(round(t1_elapsed, 1), round(t2_elapsed, 1))
)

log_msg("Model comparison:")
for (i in seq_len(nrow(comparison))) {
  r <- comparison[i, ]
  log_msg(sprintf("  %s: OFV=%.1f, AIC=%.1f (%.1f min)",
                  r$Model, r$OFV, r$AIC, r$Runtime_min))
}

write_csv(comparison, file.path(DIR_TABLES, "model_comparison.csv"))

# Select best model (lowest AIC)
if (!is.null(fit2) && !is.null(fit1)) {
  best_fit <- if (AIC(fit2) < AIC(fit1)) fit2 else fit1
  best_name <- if (AIC(fit2) < AIC(fit1)) "2-compartment" else "1-compartment"
} else if (!is.null(fit2)) {
  best_fit <- fit2
  best_name <- "2-compartment"
} else {
  best_fit <- fit1
  best_name <- "1-compartment"
}

log_msg("Selected base model: ", best_name)

# Save fits
saveRDS(fit1, file.path(DIR_TABLES, "fit_1comp.rds"))
saveRDS(fit2, file.path(DIR_TABLES, "fit_2comp.rds"))
saveRDS(best_fit, file.path(DIR_TABLES, "fit_base.rds"))

# Save parameter estimates
if (!is.null(best_fit)) {
  fe <- fixef(best_fit)
  se_vec <- tryCatch({
    v <- vcov(best_fit)
    fe_names <- names(fe)
    # vcov may contain more params than fixef — subset to matching names
    common <- intersect(fe_names, rownames(v))
    se_all <- rep(NA_real_, length(fe))
    names(se_all) <- fe_names
    se_all[common] <- sqrt(diag(v[common, common, drop = FALSE]))
    round(se_all, 4)
  }, error = function(e) rep(NA_real_, length(fe)))
  param_table <- tibble(
    Parameter = names(fe),
    Estimate  = round(unname(fe), 4),
    SE        = unname(se_vec)
  )
  write_csv(param_table, file.path(DIR_TABLES, "parameter_estimates_base.csv"))
  log_msg("Parameter estimates saved")
}

log_msg("=== Model Fitting Complete ===")
