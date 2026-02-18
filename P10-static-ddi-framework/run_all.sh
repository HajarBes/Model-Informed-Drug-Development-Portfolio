#!/bin/bash
# =============================================================================
# run_all.sh — Execute the complete Static DDI Risk Assessment pipeline
# =============================================================================
#
# Usage:
#   cd P10-static-ddi-framework
#   bash run_all.sh
#
# Prerequisites:
#   R (>= 4.1) with packages: tidyverse, readr, dplyr, ggplot2, scales, patchwork
#   (packages will be auto-installed if missing)

set -e

echo "============================================="
echo "  Static DDI Risk Assessment Framework"
echo "  FDA 2020 Mechanistic Static Model"
echo "============================================="
echo ""

cd "$(dirname "$0")"

echo "[1/4] Computing static DDI screening ratios and AUCR predictions..."
Rscript analysis/01_compute_static_ddi.R
echo ""

echo "[2/4] Running sensitivity & uncertainty analysis..."
Rscript analysis/02_sensitivity_uncertainty.R
echo ""

echo "[3/4] Generating publication-quality figures..."
Rscript analysis/03_make_figures.R
echo ""

echo "[4/4] Generating formatted summary tables..."
Rscript analysis/04_generate_tables.R
echo ""

echo "============================================="
echo "  Pipeline complete!"
echo ""
echo "  Outputs:"
echo "    figures/          — 5 publication-quality PNGs"
echo "    outputs/tables/   — CSV tables for all results"
echo "    outputs/logs/     — Run logs + Igut sanity check"
echo "    data/processed/   — RDS intermediate data"
echo "    report/           — Case reports (markdown)"
echo "============================================="
