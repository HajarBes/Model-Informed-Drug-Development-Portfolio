# ============================================================================
# export_for_nonmem.R
# Create a NONMEM-ready dataset from the analysis CSV
#
# Reads:  ../data/simulated/analysis_dataset.csv
# Writes: data_nm.csv (space-delimited, @-prefixed header for IGNORE=@)
#
# Transformations:
#   - Drop IPRED column (not needed by NONMEM)
#   - Encode ARM: RICH → 1, SPARSE → 2
#   - Replace NA with '.' (NONMEM missing value)
#   - Write space-delimited with @-prefixed header line
# ============================================================================

# --- Read source dataset ---
infile <- "../data/simulated/analysis_dataset.csv"
if (!file.exists(infile)) stop("Cannot find input file: ", infile)

cat("Reading:", infile, "\n")
dat <- read.csv(infile, stringsAsFactors = FALSE)

# --- Drop IPRED (not needed by NONMEM) ---
if ("IPRED" %in% names(dat)) dat$IPRED <- NULL

# --- Encode ARM: RICH → 1, SPARSE → 2 ---
if ("ARM" %in% names(dat)) {
  dat$ARM <- ifelse(dat$ARM == "RICH", 1,
             ifelse(dat$ARM == "SPARSE", 2, NA))
}

# --- Safety check: required columns ---
needed <- c("ID", "TIME", "AMT", "EVID", "CMT", "MDV", "DV",
            "WT", "SEX", "AGE")
missing <- setdiff(needed, names(dat))
if (length(missing) > 0) {
  stop("Missing required columns: ", paste(missing, collapse = ", "))
}

# --- Compute summary counts BEFORE character conversion ---
n_subj <- length(unique(dat$ID))
n_obs  <- sum(dat$EVID == 0, na.rm = TRUE)
n_dose <- sum(dat$EVID == 1, na.rm = TRUE)

# --- Convert NA → '.' for NONMEM (character conversion at the end) ---
dat_out <- as.data.frame(lapply(dat, function(x) {
  x <- as.character(x)
  x[is.na(x) | x == "NA"] <- "."
  x
}), stringsAsFactors = FALSE)

# --- Write space-delimited NM-ready dataset with @-prefixed header ---
outfile <- "data_nm.csv"
cat("Writing:", outfile, "\n")

con <- file(outfile, open = "wt")
writeLines(paste0("@", paste(names(dat_out), collapse = " ")), con)
write.table(dat_out, file = con, quote = FALSE, row.names = FALSE,
            col.names = FALSE, sep = " ")
close(con)

# --- Summary ---
cat(sprintf("  Subjects: %d | Observations: %d | Dose records: %d\n",
            n_subj, n_obs, n_dose))
cat("Done.\n")
