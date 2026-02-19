#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────────────
# run_all.sh — Sotorasib Exposure-Response & Dose Optimization Pipeline
# Regenerates all outputs, figures, and tables from scratch.
# Usage:  bash run_all.sh
# ──────────────────────────────────────────────────────────────────────
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_DIR="${SCRIPT_DIR}/analysis"
LOG_DIR="${SCRIPT_DIR}/outputs/logs"

mkdir -p "$LOG_DIR"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOGFILE="${LOG_DIR}/run_${TIMESTAMP}.log"

echo "═══════════════════════════════════════════════════════════════"
echo " Sotorasib Exposure-Response & Dose Optimization Pipeline"
echo " Started: $(date)"
echo " Log: ${LOGFILE}"
echo "═══════════════════════════════════════════════════════════════"

run_step() {
    local step_num="$1"
    local script="$2"
    local label="$3"
    echo ""
    echo "── Step ${step_num}: ${label} ──"
    python "${ANALYSIS_DIR}/${script}" 2>&1 | tee -a "$LOGFILE"
    echo "   ✓ ${label} complete"
}

run_step 0 "00_setup.py"                     "Configuration & sanity checks"
run_step 1 "01_simulate_virtual_patients.py"  "Virtual patient PK simulation"
run_step 2 "02_exposure_efficacy.py"          "Exposure-efficacy analysis"
run_step 3 "03_exposure_safety.py"            "Exposure-safety analysis"
run_step 4 "04_dose_optimization.py"          "Dose optimization & benefit-risk"
run_step 5 "05_covariate_forest.py"           "Covariate effects on exposure"
run_step 6 "06_generate_figures.py"           "Consolidated publication figures"

echo ""
echo "═══════════════════════════════════════════════════════════════"
echo " Pipeline complete: $(date)"
echo " Figures: ${SCRIPT_DIR}/figures/"
echo " Tables:  ${SCRIPT_DIR}/outputs/tables/"
echo " Log:     ${LOGFILE}"
echo "═══════════════════════════════════════════════════════════════"
