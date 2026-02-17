# =============================================================================
# 02_sensitivity_uncertainty.R — Tornado, heatmap grid, Monte Carlo, Igut check
# =============================================================================

script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(script_dir) || !file.exists(file.path(script_dir, "00_setup.R"))) {
  script_dir <- file.path(getwd(), "analysis")
}
source(file.path(script_dir, "00_setup.R"))

log_msg("=== Starting Sensitivity & Uncertainty Analysis ===")

# Load results from 01
res <- readRDS(file.path(DIR_PROCESSED, "ddi_results.rds"))

# =============================================================================
# IGUT SANITY CHECK
# =============================================================================

igut_check <- function(params) {
  lines <- c(
    "=== [I]gut Sanity Check ===",
    "",
    "Formula (method 1 - via ug/mL):",
    "  dose_ug  = dose_mg * 1000",
    "  Igut_ugml = dose_ug / V_gut_mL   (FDA assumes V_gut = 250 mL)",
    "  Igut_uM   = (Igut_ugml / MW) * 1000",
    "",
    "Formula (method 2 - via umol/L):",
    "  dose_umol = dose_mg / MW * 1000",
    "  V_gut_L   = 0.250 L",
    "  Igut_uM   = dose_umol / V_gut_L",
    "",
    sprintf("Perpetrator: %s", params$perpetrator),
    sprintf("Dose: %g mg", params$dose_mg),
    sprintf("MW: %g g/mol", params$mw),
    "",
    "--- Method 1 ---",
    sprintf("  dose_ug   = %g ug", params$dose_mg * 1000),
    sprintf("  Igut_ugml = %g / 250 = %g ug/mL",
            params$dose_mg * 1000, params$dose_mg * 1000 / 250),
    sprintf("  Igut_uM   = (%g / %g) * 1000 = %g uM",
            params$dose_mg * 1000 / 250, params$mw,
            (params$dose_mg * 1000 / 250 / params$mw) * 1000),
    "",
    "--- Method 2 ---",
    sprintf("  dose_umol = %g / %g * 1000 = %g umol",
            params$dose_mg, params$mw, params$dose_mg / params$mw * 1000),
    sprintf("  V_gut_L   = 0.250 L"),
    sprintf("  Igut_uM   = %g / 0.250 = %g uM",
            params$dose_mg / params$mw * 1000,
            (params$dose_mg / params$mw * 1000) / 0.250),
    "",
    sprintf("Computed [I]gut = %g uM", params$Igut),
    "",
    "--- Unit consistency assertion ---"
  )

  # Assertion: both methods agree
  m1 <- (params$dose_mg * 1000 / 250 / params$mw) * 1000
  m2 <- (params$dose_mg / params$mw * 1000) / 0.250
  check <- abs(m1 - m2) < 1e-6
  lines <- c(lines,
    sprintf("  Method 1: %g uM", m1),
    sprintf("  Method 2: %g uM", m2),
    sprintf("  Match: %s", ifelse(check, "PASS", "FAIL")),
    "",
    "--- Interpretation ---",
    "[I]gut is intentionally a worst-case screening concentration.",
    "It assumes the entire oral dose dissolves in 250 mL of GI fluid.",
    "This value can be orders of magnitude higher than systemic [I]max,u.",
    sprintf("[I]gut / [I]max,u ratio = %g / %g = %.0fx",
            params$Igut, params$Imax_u, params$Igut / params$Imax_u),
    "This conservative bias is by design (FDA 2020 guidance) to ensure",
    "potential intestinal DDIs are not missed at the screening stage.",
    "Results flagged by [I]gut should be interpreted qualitatively,",
    "not as quantitative predictions of clinical DDI magnitude.",
    ""
  )
  return(lines)
}

# Write Igut sanity check for both cases
igut_lines <- c()
for (nm in names(res)) {
  igut_lines <- c(igut_lines, igut_check(res[[nm]]$params), "")
}
writeLines(igut_lines, file.path(DIR_LOGS, "igut_sanity.txt"))
log_msg("Igut sanity check written to outputs/logs/igut_sanity.txt")

# =============================================================================
# MIDAZOLAM CASE: TORNADO SENSITIVITY
# =============================================================================

mdz <- res[["midazolam_ketoconazole"]]
p   <- mdz$params

# Base AUCR computation function
compute_aucr_from_params <- function(Ki, Imax_u, Igut, fm, fg, kinact, KI, kdeg) {
  R1        <- compute_R1(Imax_u, Ki)
  R1_gut    <- compute_R1_gut(Igut, Ki)
  tdi_fac   <- compute_TDI_factor(kinact, Imax_u, KI, kdeg)
  aucr_h    <- compute_AUCR_hepatic(fm, R1, tdi_fac)
  aucr_g    <- compute_AUCR_gut(fg, R1_gut)
  list(total = aucr_h * aucr_g, hepatic = aucr_h, gut = aucr_g)
}

base <- compute_aucr_from_params(
  Ki = 0.015, Imax_u = p$Imax_u, Igut = p$Igut,
  fm = 0.94, fg = 0.43, kinact = 0.048, KI = 0.86, kdeg = 0.00032
)

log_msg("Midazolam base AUCR: total=", round(base$total, 2),
        " hepatic=", round(base$hepatic, 2), " gut=", round(base$gut, 2))

# Define parameter ranges for tornado
tornado_params <- tribble(
  ~param,   ~label,               ~base,     ~low,      ~high,
  "Ki",     "Ki CYP3A4 (uM)",     0.015,     0.0015,    0.15,
  "Imax_u", "[I]max,u (uM)",      p$Imax_u,  p$Imax_u*0.5, p$Imax_u*1.5,
  "fm",     "fm CYP3A4",          0.94,      0.50,      0.95,
  "fg",     "fg (gut extraction)", 0.43,      0.00,      0.70,
  "kinact", "kinact (1/min)",      0.048,     0.024,     0.072,
  "KI",     "KI TDI (uM)",        0.86,      0.43,      1.29,
  "kdeg",   "kdeg CYP3A4 (1/min)",0.00032,   0.00016,   0.00048
)

# Run one-at-a-time sensitivity
base_aucr <- base$total
base_args <- list(Ki=0.015, Imax_u=p$Imax_u, Igut=p$Igut,
                  fm=0.94, fg=0.43, kinact=0.048, KI=0.86, kdeg=0.00032)

run_tornado_row <- function(param_name, low_val, high_val) {
  args_lo <- base_args
  args_lo[[param_name]] <- low_val
  args_hi <- base_args
  args_hi[[param_name]] <- high_val
  c(aucr_low  = do.call(compute_aucr_from_params, args_lo)$total,
    aucr_high = do.call(compute_aucr_from_params, args_hi)$total)
}

tornado_computed <- pmap_dfr(
  list(tornado_params$param, tornado_params$low, tornado_params$high),
  function(pn, lo, hi) as_tibble_row(run_tornado_row(pn, lo, hi))
)

tornado_data <- bind_cols(tornado_params, tornado_computed) %>%
  mutate(
    aucr_base  = base_aucr,
    delta_low  = aucr_low - aucr_base,
    delta_high = aucr_high - aucr_base,
    max_delta  = pmax(abs(delta_low), abs(delta_high)),
    range_label = sprintf("[%.4g, %.4g]", low, high)
  ) %>%
  arrange(desc(max_delta))

log_msg("Tornado sensitivity computed:")
for (i in seq_len(nrow(tornado_data))) {
  r <- tornado_data[i, ]
  log_msg(sprintf("  %s: AUCR range [%.1f, %.1f] (base %.1f)",
                  r$label, min(r$aucr_low, r$aucr_high),
                  max(r$aucr_low, r$aucr_high), r$aucr_base))
}

# =============================================================================
# MIDAZOLAM: HEATMAP GRID (Ki vs Imax_u)
# =============================================================================

ki_grid    <- 10^seq(log10(0.001), log10(10), length.out = 100)
imax_grid  <- 10^seq(log10(0.01), log10(10), length.out = 100)

heatmap_data <- expand_grid(Ki_uM = ki_grid, Imax_u_uM = imax_grid) %>%
  rowwise() %>%
  mutate(
    R1 = compute_R1(Imax_u_uM, Ki_uM),
    AUCR = compute_aucr_from_params(
      Ki = Ki_uM, Imax_u = Imax_u_uM, Igut = p$Igut,
      fm = 0.94, fg = 0.43, kinact = 0.048, KI = 0.86, kdeg = 0.00032
    )$total,
    classification = classify_ddi(AUCR)
  ) %>%
  ungroup()

log_msg("Heatmap grid computed: ", nrow(heatmap_data), " points")

# =============================================================================
# MIDAZOLAM: MONTE CARLO
# =============================================================================

set.seed(42)
N_MC <- 10000

mc_samples <- tibble(
  Ki     = rlnorm(N_MC, log(0.015), 0.3),
  Imax_u = rlnorm(N_MC, log(p$Imax_u), 0.3),
  fm     = pmin(0.99, pmax(0.3, rnorm(N_MC, 0.94, 0.05))),
  fg     = pmin(0.80, pmax(0.0, rnorm(N_MC, 0.43, 0.10))),
  kinact = rlnorm(N_MC, log(0.048), 0.3),
  KI     = rlnorm(N_MC, log(0.86), 0.3),
  kdeg   = rlnorm(N_MC, log(0.00032), 0.3)
)

mc_results <- mc_samples %>%
  rowwise() %>%
  mutate(
    AUCR = compute_aucr_from_params(Ki, Imax_u, p$Igut,
                                     fm, fg, kinact, KI, kdeg)$total,
    class = classify_ddi(AUCR)
  ) %>%
  ungroup()

mc_summary <- mc_results %>%
  summarise(
    median = median(AUCR),
    p5     = quantile(AUCR, 0.05),
    p95    = quantile(AUCR, 0.95),
    pct_strong   = mean(AUCR >= 5) * 100,
    pct_moderate = mean(AUCR >= 2 & AUCR < 5) * 100,
    pct_weak     = mean(AUCR >= 1.25 & AUCR < 2) * 100,
    pct_none     = mean(AUCR < 1.25) * 100
  )

log_msg("Monte Carlo summary (N=", N_MC, "):")
log_msg("  Median AUCR: ", round(mc_summary$median, 1))
log_msg("  90% CI: [", round(mc_summary$p5, 1), ", ", round(mc_summary$p95, 1), "]")
log_msg("  Strong: ", round(mc_summary$pct_strong, 1), "%")
log_msg("  Moderate: ", round(mc_summary$pct_moderate, 1), "%")
log_msg("  Weak: ", round(mc_summary$pct_weak, 1), "%")
log_msg("  No interaction: ", round(mc_summary$pct_none, 1), "%")

# =============================================================================
# SOTORASIB: NET EFFECT SCENARIOS
# =============================================================================

sot <- res[["sotorasib_perpetrator"]]
sp  <- sot$params

# Hypothetical victim: midazolam-like CYP3A4 substrate (fm=0.94, fg=0.43)
fm_victim <- 0.94
fg_victim <- 0.43

# Define scenarios with varying d (induction scaling factor)
scenarios <- tibble(
  scenario = c(
    "TDI only (no induction)",
    "Induction only (no TDI)",
    "Combined, d=1.0 (conservative)",
    "Combined, d=0.5 (moderate calibration)",
    "Combined, d=0.2 (empirical calibration)",
    "Observed clinical"
  ),
  use_tdi     = c(TRUE,  FALSE, TRUE,  TRUE,  TRUE,  NA),
  use_induct  = c(FALSE, TRUE,  TRUE,  TRUE,  TRUE,  NA),
  d_value     = c(NA,    1.0,   1.0,   0.5,   0.2,   NA),
  observed    = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)
)

compute_net_aucr <- function(use_tdi, use_induct, d_val, sp, fm, fg) {
  Ki_val     <- 10.0
  kinact_val <- ifelse(use_tdi, 0.04, NA)
  KI_val     <- ifelse(use_tdi, 3.5, NA)
  emax_val   <- ifelse(use_induct, 8.1, NA)
  ec50_val   <- ifelse(use_induct, 1.5, NA)

  R1      <- compute_R1(sp$Imax_u, Ki_val)
  tdi_fac <- compute_TDI_factor(kinact_val, sp$Imax_u, KI_val, 0.00032)

  if (use_induct && !is.na(emax_val)) {
    induct_fold <- 1 + d_val * emax_val * sp$Imax_u / (ec50_val + sp$Imax_u)
  } else {
    induct_fold <- 1.0
  }

  # Net CL ratio: inhibition reduces CL, induction increases CL
  cl_ratio_net <- (1 / R1) * tdi_fac * induct_fold
  AUCR <- 1 / (fm * cl_ratio_net + (1 - fm))
  return(AUCR)
}

scenario_results <- scenarios %>%
  rowwise() %>%
  mutate(
    AUCR = if (observed) {
      0.48  # Observed: midazolam AUC decreased 52%
    } else {
      compute_net_aucr(use_tdi, use_induct, d_value, sp, fm_victim, fg_victim)
    },
    direction = case_when(
      AUCR > 1.25 ~ "Inhibition dominant",
      AUCR < 0.80 ~ "Induction dominant",
      TRUE ~ "Mixed / balanced"
    ),
    recommendation = case_when(
      observed ~ "Clinical result confirms net induction",
      AUCR > 5 ~ "Strong inhibition -> clinical study needed",
      AUCR > 2 ~ "Moderate inhibition -> clinical study or PBPK",
      AUCR > 1.25 ~ "Weak inhibition -> evaluate with PBPK",
      AUCR < 0.5 ~ "Strong induction -> clinical study needed",
      AUCR < 0.8 ~ "Moderate induction -> PBPK recommended",
      TRUE ~ "Unclear net effect -> PBPK essential"
    )
  ) %>%
  ungroup()

log_msg("Sotorasib net effect scenarios:")
for (i in seq_len(nrow(scenario_results))) {
  r <- scenario_results[i, ]
  log_msg(sprintf("  %s: AUCR=%.2f (%s)", r$scenario, r$AUCR, r$direction))
}

# =============================================================================
# SAVE ALL SENSITIVITY RESULTS
# =============================================================================

sensitivity_results <- list(
  tornado      = tornado_data,
  heatmap      = heatmap_data,
  mc_samples   = mc_results,
  mc_summary   = mc_summary,
  scenarios    = scenario_results
)

saveRDS(sensitivity_results, file.path(DIR_PROCESSED, "sensitivity_results.rds"))
write_csv(tornado_data, file.path(DIR_TABLES, "tornado_midazolam.csv"))
write_csv(scenario_results, file.path(DIR_TABLES, "net_effect_scenarios_sotorasib.csv"))
write_csv(mc_summary, file.path(DIR_TABLES, "monte_carlo_summary_midazolam.csv"))

log_msg("=== Sensitivity & Uncertainty Analysis Complete ===")
