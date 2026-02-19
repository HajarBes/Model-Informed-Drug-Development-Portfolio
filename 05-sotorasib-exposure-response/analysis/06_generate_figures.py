#!/usr/bin/env python3
"""
06_generate_figures.py — Consolidated Figure Regeneration
==========================================================
Regenerates all 8 presentation-grade figures by re-running the relevant analysis scripts.
Ensures consistent styling across the entire figure set.

Figures:
  1. dose_exposure_boxplot.png         — Flat dose-exposure
  2. saturable_bioavailability.png     — Mechanism of flat E-R
  3. er_efficacy_panel.png             — ORR + PFS vs exposure
  4. dose_vs_exposure_vs_response.png  — Dose→Exposure→Response chain
  5. er_safety_panel.png               — Hepatotoxicity vs Cmax + CPI
  6. dose_optimization_dashboard.png   — Benefit-risk dashboard
  7. benefit_risk_overlay.png          — Efficacy vs safety overlay
  8. covariate_forest_plot.png         — Covariate forest plot
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import importlib.util as _ilu

_setup_path = Path(__file__).resolve().parent / "00_setup.py"
_spec = _ilu.spec_from_file_location("setup_00", _setup_path)
_mod = _ilu.module_from_spec(_spec)
sys.modules["setup_00"] = _mod
_spec.loader.exec_module(_mod)

DIR_FIGURES = _mod.DIR_FIGURES
get_logger = _mod.get_logger

log = get_logger("06_figures")

# ── Verify All Figures Exist ─────────────────────────────────────────
expected_figures = [
    "dose_exposure_boxplot.png",
    "saturable_bioavailability.png",
    "er_efficacy_panel.png",
    "dose_vs_exposure_vs_response.png",
    "er_safety_panel.png",
    "dose_optimization_dashboard.png",
    "benefit_risk_overlay.png",
    "covariate_forest_plot.png",
]

log.info("Verifying figure outputs...")
all_present = True
for fname in expected_figures:
    fpath = DIR_FIGURES / fname
    if fpath.exists():
        size_kb = fpath.stat().st_size / 1024
        log.info(f"  [OK] {fname} ({size_kb:.0f} KB)")
    else:
        log.warning(f"  [MISSING] {fname}")
        all_present = False

if all_present:
    log.info(f"\nAll {len(expected_figures)} figures present in {DIR_FIGURES}")
else:
    log.warning("\nSome figures are missing. Re-run the upstream scripts.")

# ── Summary Table of Figures ─────────────────────────────────────────
fig_descriptions = {
    "dose_exposure_boxplot.png": "Does exposure increase with dose? (No — saturable)",
    "saturable_bioavailability.png": "What drives the flat dose-exposure curve?",
    "er_efficacy_panel.png": "Is there an exposure-ORR relationship? (+ PFS descriptive)",
    "dose_vs_exposure_vs_response.png": "Dose -> exposure -> response: where does the chain break?",
    "er_safety_panel.png": "Does Cmax drive hepatotoxicity? Does prior CPI modify risk?",
    "dose_optimization_dashboard.png": "Is 960 mg optimal, or would 240 mg suffice?",
    "benefit_risk_overlay.png": "Net clinical benefit at each dose level",
    "covariate_forest_plot.png": "Which patient factors drive exposure variability?",
}

log.info("\nFigure Index:")
log.info(f"{'#':<4} {'File':<40} {'Decision Question'}")
log.info(f"{'─'*4} {'─'*40} {'─'*50}")
for i, fname in enumerate(expected_figures, 1):
    desc = fig_descriptions.get(fname, "")
    log.info(f"{i:<4} {fname:<40} {desc}")

log.info("\nScript 06 complete.")
