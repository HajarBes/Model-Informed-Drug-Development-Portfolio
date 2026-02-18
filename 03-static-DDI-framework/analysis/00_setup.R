# =============================================================================
# 00_setup.R — Static DDI Risk Assessment Framework
# Packages, paths, helper functions, DDI computation functions
# =============================================================================

# --- Packages ----------------------------------------------------------------
required_pkgs <- c("tidyverse", "readr", "dplyr", "ggplot2", "scales", "patchwork")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# --- Paths -------------------------------------------------------------------
# Resolve project root (works from Rscript and interactive)
PROJECT_ROOT <- getwd()
if (!file.exists(file.path(PROJECT_ROOT, "data", "inputs_midazolam_case.csv"))) {
  candidates <- c(
    file.path(Sys.getenv("HOME"), "Desktop",
              "Model-Informed-Drug-Development-Portfolio",
              "03-static-DDI-framework"),
    dirname(getwd())
  )
  for (p in candidates) {
    if (file.exists(file.path(p, "data", "inputs_midazolam_case.csv"))) {
      PROJECT_ROOT <- p
      break
    }
  }
}

DIR_DATA      <- file.path(PROJECT_ROOT, "data")
DIR_PROCESSED <- file.path(PROJECT_ROOT, "data", "processed")
DIR_ANALYSIS  <- file.path(PROJECT_ROOT, "analysis")
DIR_FIGURES   <- file.path(PROJECT_ROOT, "figures")
DIR_TABLES    <- file.path(PROJECT_ROOT, "outputs", "tables")
DIR_LOGS      <- file.path(PROJECT_ROOT, "outputs", "logs")

for (d in c(DIR_PROCESSED, DIR_FIGURES, DIR_TABLES, DIR_LOGS)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# --- Logging -----------------------------------------------------------------
LOG_FILE <- file.path(DIR_LOGS, paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(...))
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}

# --- ggplot2 theme -----------------------------------------------------------
theme_ddi <- theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 10, color = "grey40"),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    strip.text       = element_text(face = "bold")
  )
theme_set(theme_ddi)

COLORS_RISK <- c("No interaction" = "#2CA02C", "Weak" = "#FFD700",
                 "Moderate" = "#FF8C00", "Strong" = "#D62728")
COLOR_FLAG  <- "#D62728"
COLOR_PASS  <- "#2CA02C"

# =============================================================================
# DDI COMPUTATION FUNCTIONS
# =============================================================================

#' Convert Cmax (ug/mL) to uM
cmax_to_uM <- function(cmax_ugml, mw) {
  (cmax_ugml / mw) * 1000
}

#' Compute unbound systemic Cmax [I]max,u in uM
compute_Imax_u <- function(cmax_ugml, fu, mw) {
  cmax_uM <- cmax_to_uM(cmax_ugml, mw)
  Imax_u  <- cmax_uM * fu
  return(Imax_u)
}

#' Compute intestinal concentration [I]gut in uM
#' FDA assumption: 250 mL gastrointestinal volume
compute_Igut <- function(dose_mg, mw) {
  dose_ug  <- dose_mg * 1000
  Igut_ugml <- dose_ug / 250
  Igut_uM   <- (Igut_ugml / mw) * 1000
  return(Igut_uM)
}

#' R1: Reversible CYP inhibition screening ratio
#' R1 = 1 + [I]max,u / Ki
#' Threshold: R1 >= 1.02 triggers further evaluation
compute_R1 <- function(Imax_u, Ki) {
  if (is.na(Ki) || Ki <= 0) return(NA_real_)
  1 + Imax_u / Ki
}

#' R1,gut: Intestinal CYP inhibition screening ratio
#' R1,gut = 1 + [I]gut / Ki
#' Threshold: R1,gut >= 11 triggers further evaluation
compute_R1_gut <- function(Igut, Ki) {
  if (is.na(Ki) || Ki <= 0) return(NA_real_)
  1 + Igut / Ki
}

#' TDI factor: fraction of CYP remaining after time-dependent inactivation
#' TDI_factor = kdeg / (kdeg + lambda)
#' where lambda = kinact * [I]max,u / (KI + [I]max,u)
compute_TDI_factor <- function(kinact, Imax_u, KI, kdeg) {
  if (any(is.na(c(kinact, KI, kdeg)))) return(1.0)  # No TDI data → factor = 1
  lambda <- kinact * Imax_u / (KI + Imax_u)
  tdi_factor <- kdeg / (kdeg + lambda)
  return(tdi_factor)
}

#' R3: Induction screening ratio
#' R3 = 1 / (1 + d * Emax * [I]max,u / (EC50 + [I]max,u))
#' Threshold: R3 <= 0.8 indicates clinically relevant induction
compute_R3 <- function(Emax, Imax_u, EC50, d = 1) {
  if (any(is.na(c(Emax, EC50)))) return(1.0)  # No induction data → factor = 1
  R3 <- 1 / (1 + d * Emax * Imax_u / (EC50 + Imax_u))
  return(R3)
}

#' Predicted AUCR from mechanistic static model (hepatic only)
#' AUCR_h = 1 / (fm * CLint_ratio + (1 - fm))
#' CLint_ratio = (1/R1) * TDI_factor
compute_AUCR_hepatic <- function(fm, R1, tdi_factor = 1) {
  if (is.na(fm) || is.na(R1)) return(NA_real_)
  cl_ratio <- (1 / R1) * tdi_factor
  AUCR <- 1 / (fm * cl_ratio + (1 - fm))
  return(AUCR)
}

#' Predicted AUCR contribution from gut wall inhibition
#' AUCR_g = 1 / (fg * (1/R1_gut) + (1 - fg))
#' fg = fraction metabolized in gut wall (= 1 - Fg)
compute_AUCR_gut <- function(fg, R1_gut) {
  if (is.na(fg) || is.na(R1_gut) || fg == 0) return(1.0)
  cl_ratio_gut <- 1 / R1_gut
  AUCR <- 1 / (fg * cl_ratio_gut + (1 - fg))
  return(AUCR)
}

#' Total predicted AUCR = hepatic * gut contributions
compute_AUCR_total <- function(fm, R1, tdi_factor, fg, R1_gut) {
  AUCR_h <- compute_AUCR_hepatic(fm, R1, tdi_factor)
  AUCR_g <- compute_AUCR_gut(fg, R1_gut)
  return(AUCR_h * AUCR_g)
}

#' Classify DDI magnitude based on AUCR (FDA classification)
classify_ddi <- function(AUCR) {
  if (is.na(AUCR)) return("Not assessable")
  if (AUCR >= 5)   return("Strong")
  if (AUCR >= 2)   return("Moderate")
  if (AUCR >= 1.25) return("Weak")
  return("No interaction")
}

#' Screen transporter DDI risk
#' Returns ratio and flag status
screen_transporter <- function(I_relevant, IC50, threshold = 0.1) {
  if (is.na(IC50) || IC50 <= 0 || is.na(I_relevant)) {
    return(list(ratio = NA_real_, flag = NA))
  }
  ratio <- I_relevant / IC50
  flag  <- ratio >= threshold
  return(list(ratio = ratio, flag = flag))
}

# =============================================================================
# INPUT PARSING
# =============================================================================

#' Read and parse an input CSV file into a structured parameter list
parse_inputs <- function(csv_path) {
  df <- read_csv(csv_path, show_col_types = FALSE,
                 col_types = cols(.default = "c"))

  get_val <- function(cat, param, enzyme = NA) {
    rows <- df %>% filter(category == cat, parameter == param)
    if (!is.na(enzyme)) {
      rows <- rows %>% filter(enzyme_transporter == enzyme)
    }
    if (nrow(rows) == 0) return(NA_character_)
    rows$value[1]
  }

  get_num <- function(cat, param, enzyme = NA) {
    v <- get_val(cat, param, enzyme)
    suppressWarnings(as.numeric(v))
  }

  # Build parameter list
  params <- list(
    case_name       = get_val("info", "case_name"),
    perpetrator     = get_val("info", "perpetrator_name"),
    victim          = get_val("info", "victim_name"),
    dose_mg         = get_num("perpetrator", "dose_mg"),
    mw              = get_num("perpetrator", "mw"),
    cmax_total      = get_num("perpetrator", "cmax_total_ugml"),
    fu              = get_num("perpetrator", "fu_plasma"),
    fm_cyp3a4       = get_num("victim", "fm", "CYP3A4"),
    fg_cyp3a4       = get_num("victim", "fg", "CYP3A4"),
    kdeg_hepatic    = get_num("degradation", "kdeg_per_min", "CYP3A4"),
    kdeg_gut        = get_num("degradation", "kdeg_gut_per_min", "CYP3A4"),
    observed_aucr   = get_num("validation", "observed_aucr", "CYP3A4"),
    observed_class  = get_val("validation", "observed_class", "CYP3A4"),
    observed_effect = get_val("validation", "observed_effect", "CYP3A4")
  )

  # Derived concentrations
  params$Imax_u <- compute_Imax_u(params$cmax_total, params$fu, params$mw)
  params$Igut   <- compute_Igut(params$dose_mg, params$mw)

  # CYP inhibition data (all enzymes)
  cyp_rows <- df %>% filter(category == "inhibition") %>%
    select(enzyme = enzyme_transporter, parameter, value)
  cyp_enzymes <- unique(cyp_rows$enzyme)

  params$cyp <- map(cyp_enzymes, function(enz) {
    list(
      enzyme  = enz,
      ki      = get_num("inhibition", "ki_uM", enz),
      kinact  = get_num("inhibition", "kinact_per_min", enz),
      KI      = get_num("inhibition", "KI_uM", enz)
    )
  }) %>% set_names(cyp_enzymes)

  # Induction data
  avail_val <- get_val("induction", "available")
  emax_val  <- get_num("induction", "emax_fold", "CYP3A4")
  ec50_val  <- get_num("induction", "ec50_uM", "CYP3A4")
  d_val     <- get_num("induction", "d_scaling", "CYP3A4")

  has_induction <- !identical(avail_val, "FALSE") & !is.na(emax_val)

  params$induction <- list(
    available = has_induction,
    emax      = emax_val,
    ec50      = ec50_val,
    d         = ifelse(is.na(d_val), 1.0, d_val)
  )

  # Transporter data
  trans_rows <- df %>%
    filter(category == "transporter", parameter == "ic50_uM") %>%
    mutate(value = suppressWarnings(as.numeric(value)))
  params$transporters <- trans_rows %>%
    select(transporter = enzyme_transporter, ic50 = value) %>%
    deframe() %>% as.list()

  return(params)
}

log_msg("00_setup.R loaded successfully.")
