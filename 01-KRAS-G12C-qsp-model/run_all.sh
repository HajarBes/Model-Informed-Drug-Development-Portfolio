#!/bin/bash
# ──────────────────────────────────────────────────────────
# KRAS G12C QSP Model — Full Regeneration Script
# Generates all figures, simulation CSVs, Vpop, and SA outputs
# ──────────────────────────────────────────────────────────
set -e

echo "============================================"
echo "KRAS G12C QSP Model — Full Pipeline"
echo "============================================"

# 1. Base model: all scenarios (270 days)
echo ""
echo "[1/3] Running base model (4 scenarios, 270 days)..."
python analysis/kras_qsp.py --scenario all --t_days 270 \
    --outdir figures --save_csv --csv_dir outputs

# 2. Virtual population (200 patients x 3 scenarios, 180 days for PFS)
echo ""
echo "[2/3] Running virtual population (N=200, 180 days)..."
python analysis/vpop_sensitivity.py --n_vpop 200 --t_days 180 --vpop-only

# 3. Morris global sensitivity analysis (90 days)
echo ""
echo "[3/3] Running Morris sensitivity analysis..."
python analysis/vpop_sensitivity.py --t_days 90 --sensitivity-only

echo ""
echo "============================================"
echo "Done. All outputs regenerated."
echo "  Figures: figures/fig1-fig9"
echo "  CSVs:    outputs/"
echo "============================================"
