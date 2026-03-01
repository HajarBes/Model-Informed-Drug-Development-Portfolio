# ==============================================================================
# Project 06: Clinical Dataset Engineering for Pharmacometrics
# Script: 06_assemble_submission.R - Assemble Mock e-Submission Package
# ==============================================================================
#
# Creates a submission-ready folder structure following eCTD conventions.
# Renames files per naming conventions, creates define-style metadata,
# and documents the package assembly.
# ==============================================================================

source("R/00_setup.R")

cat("=== Assembling Mock e-Submission Package ===\n\n")

sub_root <- paths$submission

# ---- 1. Define folder structure ----
sub_dirs <- c(
  file.path(sub_root, "datasets"),
  file.path(sub_root, "control_streams"),
  file.path(sub_root, "outputs"),
  file.path(sub_root, "tables"),
  file.path(sub_root, "docs")
)
invisible(lapply(sub_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# ---- 2. Copy and rename datasets ----
cat("Copying datasets with submission naming...\n")

# NONMEM dataset
nm_source <- file.path(paths$outputs, "nm_dataset.csv")
if (file.exists(nm_source)) {
  nm_dest <- file.path(sub_root, "datasets", "pk-nm-input-cdiscpilot01.csv")
  file.copy(nm_source, nm_dest, overwrite = TRUE)
  cat(sprintf("  %s -> %s\n", basename(nm_source), basename(nm_dest)))
}

# ADNCA dataset
adnca_source <- file.path(paths$outputs, "adnca_dataset.csv")
if (file.exists(adnca_source)) {
  adnca_dest <- file.path(sub_root, "datasets", "pk-adnca-cdiscpilot01.csv")
  file.copy(adnca_source, adnca_dest, overwrite = TRUE)
  cat(sprintf("  %s -> %s\n", basename(adnca_source), basename(adnca_dest)))
}

# ---- 3. Create placeholder control stream ----
cat("\nCreating placeholder control stream...\n")

control_content <- paste(
  "$PROBLEM PK Analysis - CDISCPILOT01",
  "$INPUT C ID TIME TAD AMT DV MDV EVID CMT BLQ AGE SEX RACE WT HT BMI ARM CUMDOSE",
  "$DATA pk-nm-input-cdiscpilot01.csv IGNORE=C",
  "",
  "; Placeholder control stream",
  "; Model development performed by PM Leader",
  "",
  sep = "\n"
)

writeLines(control_content,
           file.path(sub_root, "control_streams", "run001-pk-base.ctl"))

# ---- 4. Create define-style metadata ----
cat("Creating define and variable description files...\n")

# Variable names and descriptions (for define linking)
var_descr <- read_csv(file.path(paths$spec, "var_names_descr.csv"),
                      show_col_types = FALSE)
write_csv(var_descr, file.path(sub_root, "docs", "var-names-descr.csv"))

# Dataset-level define
define <- tibble(
  dataset_name = c("pk-nm-input-cdiscpilot01.csv",
                   "pk-adnca-cdiscpilot01.csv"),
  description = c("NONMEM PK Analysis Input Dataset",
                   "NCA Analysis Input Dataset (ADNCA)"),
  location = c("datasets/pk-nm-input-cdiscpilot01.csv",
                "datasets/pk-adnca-cdiscpilot01.csv"),
  structure = c("One record per subject per event",
                "One record per subject per analyte per timepoint"),
  variable_spec = c("docs/var-names-descr.csv",
                     "docs/var-names-descr.csv")
)
write_csv(define, file.path(sub_root, "docs", "define.csv"))

# ---- 5. Create README for submission package ----
readme_content <- paste(
  "# Mock e-Submission Package - CDISCPILOT01",
  "",
  "## Structure",
  "",
  "```",
  "submission_package_mock/",
  "  datasets/          - Analysis input datasets (renamed per convention)",
  "  control_streams/   - NONMEM control files (placeholder)",
  "  outputs/           - Model output files (placeholder)",
  "  tables/            - Output tables (placeholder)",
  "  docs/              - Define, variable descriptions, assembly notes",
  "```",
  "",
  "## Naming Convention",
  "",
  "- Datasets: `{analysis-type}-{dataset-type}-{study-id}.csv`",
  "- Control streams: `run{NNN}-{description}.ctl`",
  "",
  "## Files Included",
  "",
  sprintf("- Assembled by: H. Besbassi"),
  "",
  "## Notes",
  "",
  "This is a MOCK submission package for demonstration purposes.",
  "In production, coordinate with EPOD team for:",
  "- Final folder structure requirements",
  "- Format conversion (CSV to SAS XPT if required)",
  "- Assembly server directory placement",
  "",
  sep = "\n"
)

writeLines(readme_content,
           file.path(sub_root, "docs", "README_submission.md"))

# ---- 6. Package summary ----
cat("\n--- Package Contents ---\n")
all_files <- list.files(sub_root, recursive = TRUE)
for (f in all_files) {
  fsize <- file.size(file.path(sub_root, f))
  cat(sprintf("  %s (%.1f KB)\n", f, fsize / 1024))
}

cat(sprintf("\nSubmission package assembled: %d files in %s\n",
            length(all_files), sub_root))
