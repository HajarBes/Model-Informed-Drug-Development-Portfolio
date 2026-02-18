# =============================================================================
# 00_setup.R — PopPK Analysis Setup
# Packages, paths, themes, helper functions
# =============================================================================

# --- Packages ----------------------------------------------------------------
core_pkgs <- c("tidyverse", "dplyr", "ggplot2", "readr", "patchwork", "scales")
sim_pkgs  <- c("deSolve")
fit_pkgs  <- c("nlmixr2", "rxode2")

install_if_missing <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("[setup] Installing ", pkg, "...")
      install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

install_if_missing(core_pkgs)
install_if_missing(sim_pkgs)

# nlmixr2/rxode2: attempt installation but don't fail setup if unavailable
NLMIXR2_AVAILABLE <- tryCatch({
  install_if_missing(fit_pkgs)
  TRUE
}, error = function(e) {
  message("[setup] nlmixr2/rxode2 not available: ", conditionMessage(e))
  message("[setup] Simulation scripts will work. Fitting scripts require nlmixr2.")
  FALSE
})

# --- Paths -------------------------------------------------------------------
PROJECT_ROOT <- getwd()
if (!file.exists(file.path(PROJECT_ROOT, "analysis", "00_setup.R"))) {
  candidates <- c(
    file.path(Sys.getenv("HOME"), "Desktop",
              "Model-Informed-Drug-Development-Portfolio",
              "projects", "04_poppk_anchor_midazolam_osp_calibrated"),
    dirname(getwd())
  )
  for (p in candidates) {
    if (file.exists(file.path(p, "analysis", "00_setup.R"))) {
      PROJECT_ROOT <- p
      break
    }
  }
}

DIR_DATA      <- file.path(PROJECT_ROOT, "data")
DIR_SIM       <- file.path(PROJECT_ROOT, "data", "simulated")
DIR_SOURCES   <- file.path(PROJECT_ROOT, "data", "sources")
DIR_ANALYSIS  <- file.path(PROJECT_ROOT, "analysis")
DIR_FIGURES   <- file.path(PROJECT_ROOT, "figures")
DIR_TABLES    <- file.path(PROJECT_ROOT, "outputs", "tables")
DIR_LOGS      <- file.path(PROJECT_ROOT, "outputs", "logs")

for (d in c(DIR_SIM, DIR_FIGURES, DIR_TABLES, DIR_LOGS)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# --- Logging -----------------------------------------------------------------
LOG_FILE <- file.path(DIR_LOGS, paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(...))
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}

# --- Theme -------------------------------------------------------------------
theme_poppk <- theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 10, color = "grey40"),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    strip.text       = element_text(face = "bold")
  )
theme_set(theme_poppk)

COL_PRIMARY   <- "#2166AC"
COL_SECONDARY <- "#B2182B"
COL_OBSERVED  <- "grey30"
COL_PRED      <- "#4393C3"
COL_IPRED     <- "#D6604D"

# --- Constants ---------------------------------------------------------------
# OSP-calibrated "true" parameters for oral midazolam 2-comp model
# These produce PK profiles consistent with OSP PBPK evaluation studies
TRUE_PARAMS <- list(
  ka   = 2.5,     # 1/h — effective first-order absorption rate
  cl_f = 50,      # L/h — apparent clearance (70 kg reference); CL~25 L/h, F~0.5
  vc_f = 45,      # L   — apparent central volume (70 kg)
  q_f  = 15,      # L/h — apparent intercompartmental clearance
  vp_f = 55,      # L   — apparent peripheral volume
  wt_ref = 70     # kg  — reference body weight
)

# BSV (omega^2 on log scale), informed by published PopPK variability
TRUE_OMEGA <- list(
  ka  = 0.36,    # CV ~60%
  cl  = 0.09,    # CV ~30%
  vc  = 0.04     # CV ~20%
)

# Residual error
TRUE_SIGMA <- list(
  prop = 0.20,   # proportional (20%)
  add  = 0.5     # additive (0.5 ng/mL)
)

LLOQ <- 0.5  # ng/mL — lower limit of quantification

# --- Helper functions --------------------------------------------------------

#' Two-compartment oral ODE system
#' State: depot, central, peripheral
#' Concentrations in ng/mL when dose in ug and volumes in L
ode_2comp_oral <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dDepot      <- -ka * Depot
    dCentral    <- ka * Depot - (cl/v1) * Central - (q/v1) * Central + (q/v2) * Peripheral
    dPeripheral <- (q/v1) * Central - (q/v2) * Peripheral
    list(c(dDepot, dCentral, dPeripheral))
  })
}

#' Simulate a single subject PK profile
simulate_subject_pk <- function(dose_ug, times, ka, cl, v1, q, v2) {
  parms <- c(ka = ka, cl = cl, v1 = v1, q = q, v2 = v2)
  state <- c(Depot = dose_ug, Central = 0, Peripheral = 0)
  out   <- deSolve::ode(y = state, times = c(0, times), func = ode_2comp_oral,
                        parms = parms, method = "lsoda")
  out   <- as.data.frame(out)
  # Concentration = Central / V1 (ug/L = ng/mL)
  out$CP <- out$Central / v1
  out[out$time > 0, c("time", "CP")]
}

#' Compute NCA-like AUC by trapezoidal rule
compute_auc <- function(time, conc) {
  n <- length(time)
  if (n < 2) return(NA_real_)
  sum(diff(time) * (conc[-n] + conc[-1]) / 2)
}

log_msg("00_setup.R loaded. nlmixr2 available: ", NLMIXR2_AVAILABLE)
