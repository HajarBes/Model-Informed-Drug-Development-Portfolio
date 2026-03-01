# ==============================================================================
# Project 06: Clinical Dataset Engineering for Pharmacometrics
# Script: 00_setup.R - Package loading, paths, configuration
# ==============================================================================

# ---- Required packages ----
required_packages <- c(
  "tidyverse",        # Data manipulation and visualization
  "haven",            # Read/write SAS XPT files
  "pharmaversesdtm",  # CDISC SDTM pilot data
  "pharmaverseadam",  # CDISC ADaM pilot data
  "lubridate",        # Date handling
  "janitor",          # Clean column names
  "rmarkdown",        # QC report rendering
  "shiny",            # QC dashboard
  "DT",               # Interactive tables in Shiny
  "plotly"            # Interactive plots in Shiny
)

# Install missing packages
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}
invisible(lapply(required_packages, install_if_missing))

# Load packages
invisible(lapply(required_packages, library, character.only = TRUE))

# ---- Project paths ----
proj_root <- here::here("06-clinical-dataset-engineering-qc")
if (!dir.exists(proj_root)) {
  proj_root <- getwd()  # fallback if running from project dir
}

paths <- list(
  sdtm       = file.path(proj_root, "data", "sdtm"),
  adam       = file.path(proj_root, "data", "adam"),
  spec       = file.path(proj_root, "spec"),
  outputs    = file.path(proj_root, "outputs"),
  submission = file.path(proj_root, "submission_package_mock")
)

# Create directories if needed
invisible(lapply(paths, function(p) dir.create(p, recursive = TRUE, showWarnings = FALSE)))

# ---- Study configuration ----
# All study-specific parameters in one place (not hard-coded in scripts)
config <- list(
  study_id   = "CDISCPILOT01",
  drug       = "Xanomeline",
  lloq       = 0.01,             # LLOQ from PCLLOQ (ug/mL)
  blq_rule   = "LLOQ/2",         # BLQ handling: set to LLOQ/2
  dose_unit  = "mg",
  conc_unit  = "ug/mL",
  time_unit  = "hours",
  # Covariates to extract
  vs_tests   = c("WEIGHT", "HEIGHT"),
  lb_tests   = c("CREAT", "ALT", "AST", "BILI"),
  # Nominal PK timepoints (hours post-dose) - study-specific
  nominal_times = c(0, 0.5, 1, 2, 4, 6, 8, 12, 24),
  time_tolerance = 0.5
)

cat("Setup complete.\n")
cat(sprintf("Study: %s (%s)\n", config$study_id, config$drug))
cat(sprintf("LLOQ: %s %s\n", config$lloq, config$conc_unit))
cat(sprintf("Project root: %s\n", proj_root))
