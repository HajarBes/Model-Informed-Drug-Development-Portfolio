#!/bin/bash
# =============================================================================
# run_all.sh — Execute the complete PopPK analysis pipeline
# =============================================================================
#
# Usage:
#   cd 02-POPPK-midazolam-cyp3a4-probe
#   bash run_all.sh
#
# Prerequisites:
#   R (>= 4.2) with: tidyverse, deSolve, nlmixr2, rxode2, patchwork, scales
#   (Step 1 works without nlmixr2; steps 2-7 require it)

set -e

echo "============================================="
echo "  PopPK Analysis: Oral Midazolam"
echo "  OSP-Calibrated | nlmixr2 + rxode2"
echo "============================================="
echo ""

cd "$(dirname "$0")"

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

echo "============================================="
echo "  Pipeline complete!"
echo ""
echo "  Outputs:"
echo "    data/simulated/   — Analysis dataset + dictionary"
echo "    figures/           — GOF, VPC, etas, forest, exposure, DDI bridge"
echo "    outputs/tables/    — Parameter estimates, verified AUCR, DDI summary"
echo "    outputs/logs/      — Run logs"
echo "    report/            — Submission-style report + methods + exec summary"
echo "============================================="
