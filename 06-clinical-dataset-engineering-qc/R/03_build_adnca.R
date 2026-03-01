# ==============================================================================
# Project 06: Clinical Dataset Engineering for Pharmacometrics
# Script: 03_build_adnca.R - Build ADNCA Input Dataset
# ==============================================================================
#
# Generates an NCA-ready analysis dataset (ADNCA) from SDTM PC/EX domains.
# Produces concentration-time profiles per subject suitable for NCA parameter
# derivation (Cmax, Tmax, AUC).
#
# Output: adnca_dataset.csv + adnca_metadata.csv
# ==============================================================================

source("R/00_setup.R")

# ---- 1. Load source data ----
cat("Building ADNCA dataset...\n")

pc <- read_csv(file.path(paths$sdtm, "pc.csv"), show_col_types = FALSE)
ex <- read_csv(file.path(paths$sdtm, "ex.csv"), show_col_types = FALSE)
dm <- read_csv(file.path(paths$sdtm, "dm.csv"), show_col_types = FALSE)

# Load ADSL for subject-level variables
adsl <- read_csv(file.path(paths$adam, "adsl.csv"), show_col_types = FALSE)

# ---- 2. Derive ADNCA variables ----

# Subject-level info from ADSL
subj_info <- adsl %>%
  select(USUBJID, STUDYID, SITEID, SUBJID,
         any_of(c("TRT01P", "TRT01A", "AGE", "AGEU", "SEX", "RACE",
                   "SAFFL", "PKFL")))

# PK concentrations with NCA-relevant derivations
adnca <- pc %>%
  left_join(subj_info, by = c("USUBJID", "STUDYID")) %>%
  mutate(
    # Parse datetime
    ADTM = parse_date_time(PCDTC, orders = c("ymd HMS", "ymd")),
    # Analyte info
    PARAM   = PCTEST,
    PARAMCD = PCTESTCD,
    # Concentration result
    AVAL    = PCSTRESN,
    AVALU   = PCSTRESU,
    # BLQ flag
    ABLFL   = if_else(!is.na(AVAL) & AVAL < config$lloq, "Y", "N"),
    # Analysis flags
    ANL01FL = "Y"  # Include in primary NCA analysis
  )

# Derive nominal and actual time relative to dose
# TODO: Link to specific dose records for TAD derivation
# This requires matching each PC record to the most recent EX record

# ---- 3. Add NCA-specific metadata columns ----
adnca <- adnca %>%
  mutate(
    DTYPE   = NA_character_,  # Derivation type (for derived records)
    ATPT    = VISIT,          # Analysis timepoint (map from visit)
    SRCDOM  = "PC",           # Source domain
    SRCSEQ  = PCSEQ           # Source sequence number
  )

# ---- 4. Select and order ADNCA columns ----
adnca_final <- adnca %>%
  select(
    STUDYID, USUBJID, SUBJID, SITEID,
    any_of(c("TRT01P", "TRT01A")),
    PARAM, PARAMCD,
    AVAL, AVALU,
    ADTM, ATPT,
    ABLFL, ANL01FL, DTYPE,
    SRCDOM, SRCSEQ,
    any_of(c("AGE", "SEX", "RACE"))
  ) %>%
  arrange(USUBJID, PARAMCD, ADTM)

# ---- 5. Generate metadata ----
adnca_metadata <- tibble(
  variable = names(adnca_final),
  label = c(
    "Study Identifier", "Unique Subject Identifier", "Subject Identifier",
    "Site Identifier",
    # Treatment columns (conditional)
    if ("TRT01P" %in% names(adnca_final)) c("Planned Treatment", "Actual Treatment") else character(0),
    "Parameter", "Parameter Code",
    "Analysis Value", "Analysis Value Unit",
    "Analysis Datetime", "Analysis Timepoint",
    "Baseline Flag (BLQ)", "Analysis Flag 01", "Derivation Type",
    "Source Domain", "Source Sequence Number",
    if ("AGE" %in% names(adnca_final)) c("Age", "Sex", "Race") else character(0)
  )[seq_along(names(adnca_final))],
  type = sapply(adnca_final, class) %>% sapply(`[`, 1),
  source = "Derived from SDTM PC + ADSL"
)

# ---- 6. Save outputs ----
write_csv(adnca_final, file.path(paths$outputs, "adnca_dataset.csv"))
write_csv(adnca_metadata, file.path(paths$outputs, "adnca_metadata.csv"))

cat(sprintf("\nADNCA dataset saved: %d records, %d subjects, %d parameters\n",
            nrow(adnca_final), n_distinct(adnca_final$USUBJID),
            n_distinct(adnca_final$PARAMCD)))
cat(sprintf("Metadata saved: %d variables documented\n", nrow(adnca_metadata)))
