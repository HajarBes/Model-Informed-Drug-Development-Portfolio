# ==============================================================================
# Project 06: Clinical Dataset Engineering for Pharmacometrics
# Script: 02_build_nm_dataset.R - Build NONMEM-ready PK analysis dataset
# ==============================================================================
#
# Input:  SDTM domains (DM, EX, PC, VS, LB) from CDISCPILOT01
# Output: NONMEM-ready analysis dataset (nm_dataset.csv)
#
# Study: Xanomeline transdermal patch, Alzheimer's disease
# Drug: Xanomeline (analyte code: XAN)
# LLOQ: 0.01 ug/mL
#
# Key derivations:
#   - TIME: hours from first dose
#   - TAD: time after most recent dose
#   - AMT/DV: dose amounts and observed concentrations
#   - BLQ: below limit of quantification flag and handling
#   - MDV/EVID: missing dependent variable / event ID
#   - Baseline covariates: AGE, SEX, WT, HT, BMI, CREAT, ALT, AST, BILI
# ==============================================================================

source("R/00_setup.R")
source("R/utils/derive_tad.R")
source("R/utils/derive_blq.R")
source("R/utils/derive_dose_history.R")
source("R/utils/format_nmtran.R")

# ---- 1. Load SDTM domains ----
cat("Loading SDTM domains...\n")
dm <- read_csv(file.path(paths$sdtm, "dm.csv"), show_col_types = FALSE)
ex <- read_csv(file.path(paths$sdtm, "ex.csv"), show_col_types = FALSE)
pc <- read_csv(file.path(paths$sdtm, "pc.csv"), show_col_types = FALSE)
vs <- read_csv(file.path(paths$sdtm, "vs.csv"), show_col_types = FALSE)
lb <- read_csv(file.path(paths$sdtm, "lb.csv"), show_col_types = FALSE)

# ---- 2. Build subject-level covariates ----
cat("Deriving baseline covariates...\n")

# Demographics
demo <- dm %>%
  transmute(
    USUBJID,
    AGE   = AGE,
    SEX   = if_else(SEX == "F", 0L, 1L),  # 0=Female, 1=Male
    RACE  = case_when(
      RACE == "WHITE" ~ 1L,
      RACE == "BLACK OR AFRICAN AMERICAN" ~ 2L,
      RACE == "ASIAN" ~ 3L,
      TRUE ~ 4L
    ),
    ARM   = ARMCD
  )

# Baseline vitals (weight, height from VS domain, VSBLFL="Y")
baseline_vs <- vs %>%
  filter(VSBLFL == "Y", VSTESTCD %in% c("WEIGHT", "HEIGHT")) %>%
  select(USUBJID, VSTESTCD, VSSTRESN) %>%
  pivot_wider(names_from = VSTESTCD, values_from = VSSTRESN, values_fn = first)

# Derive WT, HT, BMI - handle missing columns gracefully
baseline_vs <- baseline_vs %>%
  transmute(
    USUBJID,
    WT  = if ("WEIGHT" %in% names(.)) WEIGHT else NA_real_,
    HT  = if ("HEIGHT" %in% names(.)) HEIGHT else NA_real_,
    BMI = if_else(!is.na(WT) & !is.na(HT) & HT > 0,
                  WT / (HT / 100)^2, NA_real_)
  )

# Baseline labs (renal/hepatic function from LB domain, LBBLFL="Y")
baseline_lb <- lb %>%
  filter(LBBLFL == "Y", LBTESTCD %in% c("CREAT", "ALT", "AST", "BILI")) %>%
  select(USUBJID, LBTESTCD, LBSTRESN) %>%
  pivot_wider(names_from = LBTESTCD, values_from = LBSTRESN, values_fn = first)

# Rename LB columns to match spec
lb_names <- intersect(c("CREAT", "ALT", "AST", "BILI"), names(baseline_lb))
baseline_lb <- baseline_lb %>% select(USUBJID, all_of(lb_names))

# Merge covariates
covariates <- demo %>%
  left_join(baseline_vs, by = "USUBJID") %>%
  left_join(baseline_lb, by = "USUBJID")

cat(sprintf("  Covariates derived for %d subjects\n", nrow(covariates)))

# ---- 3. Build dosing records ----
cat("Processing dosing records...\n")

# EX: EXSTDTC is Date class (no time), convert to datetime at midnight
doses <- ex %>%
  filter(EXDOSE > 0) %>%  # Exclude placebo zero-dose records
  transmute(
    USUBJID,
    DATETIME = as.POSIXct(EXSTDTC, tz = "UTC"),  # Date -> POSIXct at midnight
    AMT      = EXDOSE,
    DV       = NA_real_,
    MDV      = 1L,
    EVID     = 1L,
    CMT      = 1L,
    BLQ      = 0L
  )

cat(sprintf("  Dosing records: %d (excluded %d placebo zero-dose)\n",
            nrow(doses), sum(ex$EXDOSE == 0)))

# ---- 4. Build observation records ----
cat("Processing PK observations...\n")

# PC: PCDTC is POSIXct, PCSTRESN is numeric (NA for BLQ)
# Filter PLASMA only (exclude URINE samples to avoid duplicates at same timepoint)
observations <- pc %>%
  filter(PCSPEC == "PLASMA") %>%
  transmute(
    USUBJID,
    DATETIME = PCDTC,  # Already POSIXct
    AMT      = 0,
    DV       = PCSTRESN,
    MDV      = if_else(is.na(PCSTRESN), 1L, 0L),
    EVID     = 0L,
    CMT      = 2L,
    BLQ      = derive_blq_flag(PCSTRESN, config$lloq)
  )

# Handle NAs in PCSTRESN (these are BLQ: PCORRES = "<BLQ")
# For records where PCSTRESN is NA but PCORRES indicates BLQ, flag as BLQ
observations <- observations %>%
  mutate(
    BLQ = if_else(is.na(DV), 1L, BLQ)
  )

# Apply BLQ handling rule
observations <- apply_blq_rule(observations, config$lloq, config$blq_rule)

cat(sprintf("  Observations: %d (BLQ: %d, quantifiable: %d)\n",
            nrow(observations), sum(observations$BLQ == 1),
            sum(observations$BLQ == 0)))

# ---- 5. Combine and derive time variables ----
cat("Combining records and deriving time variables...\n")

nm_data <- bind_rows(
  doses %>% select(USUBJID, DATETIME, AMT, DV, MDV, EVID, CMT, BLQ),
  observations %>% select(USUBJID, DATETIME, AMT, DV, MDV, EVID, CMT, BLQ)
) %>%
  arrange(USUBJID, DATETIME, desc(EVID)) %>%  # Doses before obs at same time
  group_by(USUBJID) %>%
  mutate(
    # Reference time = first dose for this subject
    REFTIME = min(DATETIME[EVID == 1], na.rm = TRUE),
    TIME    = as.numeric(difftime(DATETIME, REFTIME, units = "hours")),
    TAD     = derive_tad(DATETIME, EVID, TIME),
    ID      = cur_group_id()
  ) %>%
  ungroup()

# Remove subjects with no dosing records (can't derive TIME)
no_dose_subjects <- nm_data %>%
  group_by(USUBJID) %>%
  filter(all(EVID == 0)) %>%
  pull(USUBJID) %>%
  unique()

if (length(no_dose_subjects) > 0) {
  cat(sprintf("  Removing %d subjects with no dose records (placebo/screen failure)\n",
              length(no_dose_subjects)))
  nm_data <- nm_data %>% filter(!USUBJID %in% no_dose_subjects)
  # Reassign IDs
  nm_data <- nm_data %>%
    mutate(ID = as.integer(factor(USUBJID)))
}

# ---- 6. Merge covariates ----
cat("Merging covariates...\n")

nm_data <- nm_data %>%
  left_join(covariates, by = "USUBJID")

# ---- 7. Derive cumulative dose ----
nm_data <- derive_cumulative_dose(nm_data)

# ---- 8. Format for NMTRAN ----
cat("Formatting for NMTRAN...\n")

nm_final <- nm_data %>%
  format_nmtran() %>%
  select(
    C,          # Comment flag
    ID,         # Subject ID
    TIME,       # Time from first dose (hours)
    TAD,        # Time after most recent dose (hours)
    AMT,        # Dose amount (mg)
    DV,         # Dependent variable (ug/mL)
    MDV,        # Missing DV flag
    EVID,       # Event ID
    CMT,        # Compartment
    BLQ,        # BLQ flag
    AGE,        # Age (years)
    SEX,        # Sex (0=F, 1=M)
    RACE,       # Race (numeric)
    WT,         # Weight (kg)
    HT,         # Height (cm)
    BMI,        # BMI (kg/m2)
    ARM,        # Treatment arm code
    any_of(c("CREAT", "ALT", "AST", "BILI")),  # Lab covariates
    CUMDOSE,    # Cumulative dose (mg)
    USUBJID     # Keep for traceability
  )

# Validate
validation <- validate_nmtran(nm_final)

# ---- 9. Save ----
write_csv(nm_final, file.path(paths$outputs, "nm_dataset.csv"), na = ".")
cat(sprintf("\nNONMEM dataset saved: %d records, %d subjects, %d variables\n",
            nrow(nm_final), n_distinct(nm_final$ID), ncol(nm_final)))
cat(sprintf("  Dosing records: %d\n", sum(nm_final$EVID == 1)))
cat(sprintf("  Observations: %d (quantifiable: %d, BLQ: %d)\n",
            sum(nm_final$EVID == 0),
            sum(nm_final$EVID == 0 & nm_final$BLQ == 0),
            sum(nm_final$EVID == 0 & nm_final$BLQ == 1)))
cat(sprintf("  Time range: %.1f to %.1f hours\n",
            min(nm_final$TIME, na.rm = TRUE), max(nm_final$TIME, na.rm = TRUE)))
