# ==============================================================================
# Project 06: Clinical Dataset Engineering for Pharmacometrics
# Script: 01_extract_sdtm.R - Extract SDTM domains from pharmaversesdtm
# ==============================================================================

source("R/00_setup.R")

# ---- Extract SDTM domains ----
cat("Extracting SDTM domains from pharmaversesdtm...\n")

# Demographics
dm <- pharmaversesdtm::dm %>%
  select(STUDYID, USUBJID, SUBJID, SITEID, AGE, AGEU, SEX, RACE, ETHNIC,
         ARM, ARMCD, ACTARM, ACTARMCD, COUNTRY, RFSTDTC, RFENDTC)
write_csv(dm, file.path(paths$sdtm, "dm.csv"))
cat(sprintf("  DM: %d subjects, %d variables\n", nrow(dm), ncol(dm)))

# Exposure (dosing records)
ex <- pharmaversesdtm::ex %>%
  select(STUDYID, USUBJID, EXSEQ, EXTRT, EXDOSE, EXDOSU, EXDOSFRM,
         EXROUTE, EXSTDTC, EXENDTC, VISITNUM, VISIT)
write_csv(ex, file.path(paths$sdtm, "ex.csv"))
cat(sprintf("  EX: %d records, %d variables\n", nrow(ex), ncol(ex)))

# Pharmacokinetic Concentrations
pc <- pharmaversesdtm::pc
write_csv(pc, file.path(paths$sdtm, "pc.csv"))
cat(sprintf("  PC: %d records, %d variables\n", nrow(pc), ncol(pc)))

# Vital Signs (for baseline covariates)
vs <- pharmaversesdtm::vs %>%
  select(STUDYID, USUBJID, VSSEQ, VSTESTCD, VSTEST, VSORRES, VSORRESU,
         VSSTRESC, VSSTRESN, VSSTRESU, VISITNUM, VISIT, VSDTC, VSBLFL)
write_csv(vs, file.path(paths$sdtm, "vs.csv"))
cat(sprintf("  VS: %d records, %d variables\n", nrow(vs), ncol(vs)))

# Laboratory (for renal/hepatic function covariates)
lb <- pharmaversesdtm::lb
write_csv(lb, file.path(paths$sdtm, "lb.csv"))
cat(sprintf("  LB: %d records, %d variables\n", nrow(lb), ncol(lb)))

# ---- Extract ADaM reference datasets ----
cat("\nExtracting ADaM reference datasets from pharmaverseadam...\n")

adsl <- pharmaverseadam::adsl
write_csv(adsl, file.path(paths$adam, "adsl.csv"))
cat(sprintf("  ADSL: %d subjects, %d variables\n", nrow(adsl), ncol(adsl)))

# Check if ADPC exists
if ("adpc" %in% ls(asNamespace("pharmaverseadam"))) {
  adpc <- pharmaverseadam::adpc
  write_csv(adpc, file.path(paths$adam, "adpc.csv"))
  cat(sprintf("  ADPC: %d records, %d variables\n", nrow(adpc), ncol(adpc)))
}

cat("\nExtraction complete. Files saved to data/sdtm/ and data/adam/\n")

# ---- Quick summary ----
cat("\n--- Domain Summary ---\n")
cat(sprintf("Unique subjects (DM): %d\n", n_distinct(dm$USUBJID)))
cat(sprintf("Treatment arms: %s\n", paste(unique(dm$ARM), collapse = ", ")))
cat(sprintf("Dosing records (EX): %d\n", nrow(ex)))
cat(sprintf("PK observations (PC): %d\n", nrow(pc)))
