#!/bin/bash
# =============================================================================
# run_literature_qualification.sh
# Standalone wrapper: extract OSP profiles then run R qualification analysis
# =============================================================================
#
# Usage:
#   bash scripts/run_literature_qualification.sh
#
# Prerequisites:
#   python3 with openpyxl
#   R (>= 4.1) with: tidyverse, deSolve, ggplot2, scales, patchwork

set -e

cd "$(dirname "$0")/.."

echo "============================================="
echo "  External Literature Qualification Module"
echo "  OSP Mean Profile Overlay Analysis"
echo "============================================="
echo ""

# Check prerequisites
if ! command -v python3 &> /dev/null; then
  echo "ERROR: python3 not found. Install Python 3 with openpyxl."
  exit 1
fi

if ! command -v Rscript &> /dev/null; then
  echo "ERROR: Rscript not found. Install R >= 4.1."
  exit 1
fi

echo "[1/2] Extracting OSP midazolam PO profiles (Python)..."
python3 analysis/08_extract_osp_profiles.py
echo ""

echo "[2/2] Running literature qualification analysis (R)..."
Rscript analysis/09_literature_qualification.R
echo ""

echo "============================================="
echo "  Literature qualification complete!"
echo ""
echo "  Outputs:"
echo "    data/literature_osp_midazolam/ — Extracted profiles + provenance"
echo "    figures/lit_*.png              — 4 qualification figures"
echo "    outputs/tables/                — Qualification summary CSV"
echo "============================================="
