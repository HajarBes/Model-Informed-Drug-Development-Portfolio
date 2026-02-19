#!/bin/bash
# =============================================================================
# run_all.sh — Execute the PBBM Oral Absorption pipeline for Felodipine
# =============================================================================
#
# Usage:
#   cd 05-PBBM-oral-absorption-felodipine
#   bash run_all.sh
#
# Prerequisites:
#   Python 3.10+ with packages: numpy, scipy, matplotlib, pandas, SALib, openpyxl
#   OSP observed-data xlsx at ../02-POPPK-midazolam-cyp3a4-probe/data/literature_osp/ObsDataPK_OSP.xlsx

set -e

cd "$(dirname "$0")"

echo "============================================="
echo "  PBBM Oral Absorption Model — Felodipine"
echo "  ACAT Dissolution + Absorption | BCS II"
echo "============================================="
echo ""

# --- GATING: Entry criteria must pass before ACAT model work ---
echo "=== ENTRY CRITERIA (A-E) ==="
echo ""

echo "[1/8] Extracting felodipine PK profiles from OSP database..."
python3 analysis/01_extract_osp_felodipine.py
echo ""

echo "[2/8] Qualification set audit + observed food effect..."
python3 analysis/01b_qualify_and_food_effect.py
echo ""

echo "[3/8] Verifying ACAT model (mass balance + limiting cases)..."
python3 analysis/acat_model_02.py
echo ""

# --- MODEL WORK (proceeds only after gating checks pass) ---
echo "=== MODEL SIMULATIONS ==="
echo ""

echo "[4/8] Simulating fasted-state PK + qualification..."
python3 analysis/03_simulate_fasted.py
echo ""

echo "[5/8] Simulating fed-state PK + food effect analysis..."
python3 analysis/04_simulate_fed.py
echo ""

echo "[6/8] Running Morris sensitivity analysis..."
python3 analysis/05_sensitivity_analysis.py
echo ""

echo "[7/8] Running formulation scenarios + virtual BE..."
python3 analysis/06_formulation_scenarios.py
echo ""

echo "[8/8] Generating consolidated figures..."
python3 analysis/07_generate_figures.py
echo ""

echo "============================================="
echo "  Pipeline complete!"
echo ""
echo "  Outputs:"
echo "    figures/          — 9 publication-quality PNGs"
echo "    outputs/tables/   — CSV tables for all results"
echo "    outputs/logs/     — Run logs"
echo "    data/observed/    — Extracted observed PK profiles"
echo "    data/sources/     — Provenance, assumptions, qualification"
echo "    report/           — Regulatory-style documentation"
echo "============================================="
