#!/bin/bash
# =============================================================================
# run_all.sh — Execute the PopPK analysis pipeline
# =============================================================================
#
# Usage:
#   cd 02-POPPK-midazolam-cyp3a4-probe
#   bash run_all.sh                          # default: simulated PopPK only
#   bash run_all.sh simulated                # same as default
#   bash run_all.sh literature_qualification # OSP literature qualification only
#   bash run_all.sh all                      # full pipeline (simulated + qualification)
#
# Prerequisites:
#   R (>= 4.2) with: tidyverse, deSolve, nlmixr2, rxode2, patchwork, scales
#   python3 with openpyxl (for literature qualification only)
#   (Step 1 works without nlmixr2; steps 2-7 require it)

set -e

cd "$(dirname "$0")"

MODE="${1:-simulated}"

case "$MODE" in
  simulated|literature_qualification|all)
    ;;
  *)
    echo "Usage: bash run_all.sh [simulated|literature_qualification|all]"
    echo ""
    echo "  simulated                 Simulated PopPK pipeline (steps 1-7, default)"
    echo "  literature_qualification  OSP literature qualification (steps 8-9)"
    echo "  all                       Full pipeline (steps 1-9)"
    exit 1
    ;;
esac

echo "============================================="
echo "  PopPK Analysis: Oral Midazolam"
echo "  OSP-Calibrated | nlmixr2 + rxode2"
echo "  Mode: $MODE"
echo "============================================="
echo ""

# --- Simulated PopPK pipeline (steps 1-7) ----------------------------------
if [ "$MODE" = "simulated" ] || [ "$MODE" = "all" ]; then

echo "[1/7] Simulating trial dataset..."
Rscript analysis/01_simulate_trial_dataset.R
echo ""

echo "[2/7] Fitting candidate models (1-comp + 2-comp)..."
Rscript analysis/02_fit_models.R
echo ""

echo "[3/7] Generating model diagnostics..."
Rscript analysis/03_model_diagnostics.R
echo ""

echo "[4/7] Covariate analysis (allometric WT)..."
Rscript analysis/04_covariates.R
echo ""

echo "[5/7] Decision simulations..."
Rscript analysis/05_simulation_decision.R
echo ""

echo "[6/7] AUCR verification (expanded weight strata)..."
Rscript analysis/06_verify_aucr.R
echo ""

echo "[7/7] CYP3A4 inhibition bridge scenario..."
Rscript analysis/07_cyp3a4_inhibition_bridge.R
echo ""

fi

# --- Literature qualification (steps 8-9) -----------------------------------
if [ "$MODE" = "literature_qualification" ] || [ "$MODE" = "all" ]; then

echo "[8/9] Extracting OSP midazolam PO profiles (Python)..."
python3 analysis/08_extract_osp_profiles.py
echo ""

echo "[9/9] Literature qualification analysis (R)..."
Rscript analysis/09_literature_qualification.R
echo ""

fi

# --- Summary ----------------------------------------------------------------
echo "============================================="
echo "  Pipeline complete! (mode: $MODE)"
echo ""
if [ "$MODE" = "simulated" ] || [ "$MODE" = "all" ]; then
echo "  Simulated PopPK outputs:"
echo "    data/simulated/   — Analysis dataset + dictionary"
echo "    figures/           — GOF, VPC, etas, forest, exposure, DDI bridge"
echo "    outputs/tables/    — Parameter estimates, verified AUCR, DDI summary"
echo "    outputs/logs/      — Run logs"
echo "    report/            — Submission-style report + methods + exec summary"
fi
if [ "$MODE" = "literature_qualification" ] || [ "$MODE" = "all" ]; then
echo "  Literature qualification outputs:"
echo "    data/literature_osp_midazolam/ — Extracted profiles + provenance"
echo "    figures/lit_*.png              — 4 qualification figures"
echo "    outputs/tables/                — Qualification summary CSV"
fi
echo "============================================="
