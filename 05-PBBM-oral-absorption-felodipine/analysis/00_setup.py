#!/usr/bin/env python3
# =============================================================================
# 00_setup.py — Configuration, Constants, GI Physiology, Drug Parameters
# =============================================================================
#
# Shared configuration for the PBBM oral absorption model (felodipine).
# Imported by all downstream analysis scripts.
#
# Author: Hajar Besbassi

from __future__ import annotations

import os
import sys
import logging
from dataclasses import dataclass, field
from datetime import datetime
from typing import List, Optional

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# =============================================================================
# Path resolution
# =============================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
PORTFOLIO_ROOT = os.path.dirname(PROJECT_ROOT)

DIR_DATA = os.path.join(PROJECT_ROOT, "data")
DIR_OBSERVED = os.path.join(DIR_DATA, "observed")
DIR_OBSERVED_RAW = os.path.join(DIR_OBSERVED, "raw")
DIR_OBSERVED_EXTRACTED = os.path.join(DIR_OBSERVED, "extracted")
DIR_SOURCES = os.path.join(DIR_DATA, "sources")
DIR_FIGURES = os.path.join(PROJECT_ROOT, "figures")
DIR_TABLES = os.path.join(PROJECT_ROOT, "outputs", "tables")
DIR_LOGS = os.path.join(PROJECT_ROOT, "outputs", "logs")

# OSP data file (shared across portfolio)
OSP_XLSX = os.path.join(
    PORTFOLIO_ROOT,
    "02-POPPK-midazolam-cyp3a4-probe",
    "data",
    "literature_osp",
    "ObsDataPK_OSP.xlsx",
)

# Auto-create output directories
for d in [DIR_OBSERVED_RAW, DIR_OBSERVED_EXTRACTED, DIR_SOURCES,
          DIR_FIGURES, DIR_TABLES, DIR_LOGS]:
    os.makedirs(d, exist_ok=True)

# =============================================================================
# Logging infrastructure
# =============================================================================
def setup_logging(script_name: str) -> logging.Logger:
    """Create a logger that writes to both console and timestamped log file."""
    logger = logging.getLogger(script_name)
    logger.setLevel(logging.INFO)

    if logger.handlers:
        return logger

    fmt = logging.Formatter("[%(asctime)s] %(levelname)s — %(message)s",
                            datefmt="%H:%M:%S")

    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(fmt)
    logger.addHandler(ch)

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = os.path.join(DIR_LOGS, f"run_{script_name}_{ts}.log")
    fh = logging.FileHandler(log_path)
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    return logger

# =============================================================================
# Matplotlib publication theme
# =============================================================================
def set_pub_style():
    """Apply publication-quality matplotlib styling."""
    plt.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "figure.figsize": (8, 5),
        "axes.spines.top": False,
        "axes.spines.right": False,
        "lines.linewidth": 1.5,
        "lines.markersize": 5,
        "axes.grid": False,
        "font.family": "sans-serif",
    })

set_pub_style()

# =============================================================================
# Random seed registry
# =============================================================================
SEED_DEFAULT = 42
SEED_VPOP = 20240501
SEED_SENSITIVITY = 20240601

# =============================================================================
# Drug parameters — Felodipine
# =============================================================================
@dataclass
class DrugParams:
    """Felodipine physicochemical and PK parameters."""
    name: str = "Felodipine"
    MW: float = 384.26        # g/mol (PubChem CID 3333)
    pKa: float = 5.07         # weak base (Raggi et al. 2013)
    logP: float = 3.86        # DrugBank
    S0_ug_mL: float = 0.50    # Intrinsic solubility at neutral pH (µg/mL)
    Peff_cm_s: float = 5.0e-4 # Effective permeability (high, BCS II)
    particle_radius_um: float = 25.0  # Reference IR formulation (µm)
    density_g_cm3: float = 1.3        # Particle density
    Deff_cm2_s: float = 5.0e-6        # Effective diffusion coefficient

    # Dissolution scaling factor (accounts for precipitate fineness + GI mixing)
    # Precipitated drug from stomach supersaturation forms fine amorphous particles
    # with much higher effective surface area than the original formulation.
    z_dissolution: float = 100.0

    # First-pass metabolism
    Fg: float = 0.50          # Gut wall availability (Lundahl 1997)
    Fh: float = 0.50          # Hepatic availability (Lundahl 1997)

    # Systemic PK — 2-compartment (Edgar 1987, Blychert 1990)
    CL_L_h: float = 70.0     # Clearance (L/h)
    Vc_L: float = 90.0       # Central volume (L)
    Q_L_h: float = 30.0      # Inter-compartmental clearance (L/h)
    Vp_L: float = 200.0      # Peripheral volume (L)

    # Dosing
    dose_mg: float = 10.0    # Standard IR dose

    @property
    def S0_mg_mL(self) -> float:
        return self.S0_ug_mL / 1000.0

    @property
    def F_oral(self) -> float:
        return self.Fg * self.Fh

# =============================================================================
# GI segment physiology
# =============================================================================
@dataclass
class GISegment:
    """Physiology for a single GI tract segment."""
    name: str
    pH: float
    transit_time_h: float          # Mean residence time (h)
    surface_area_cm2: float        # Absorptive surface area
    volume_mL: float               # Luminal fluid volume
    bile_factor: float = 1.0       # Solubility enhancement from bile salts
    absorbs: bool = True           # Whether absorption occurs

    @property
    def k_transit_h(self) -> float:
        """First-order transit rate constant (1/h)."""
        return 1.0 / self.transit_time_h if self.transit_time_h > 0 else 0.0

# --- Fasted-state GI parameters ---
FASTED_GI = [
    GISegment("Stomach",   pH=1.7, transit_time_h=0.25, surface_area_cm2=0.0,
              volume_mL=250.0, bile_factor=1.0, absorbs=False),
    GISegment("Duodenum",  pH=6.0, transit_time_h=0.25, surface_area_cm2=150.0,
              volume_mL=50.0,  bile_factor=2.0),
    GISegment("Jejunum1",  pH=6.4, transit_time_h=0.75, surface_area_cm2=5400.0,
              volume_mL=100.0, bile_factor=1.5),
    GISegment("Jejunum2",  pH=6.6, transit_time_h=0.75, surface_area_cm2=5400.0,
              volume_mL=100.0, bile_factor=1.2),
    GISegment("Ileum1",    pH=7.0, transit_time_h=1.00, surface_area_cm2=3600.0,
              volume_mL=100.0, bile_factor=1.0),
    GISegment("Ileum2",    pH=7.2, transit_time_h=1.00, surface_area_cm2=3600.0,
              volume_mL=100.0, bile_factor=1.0),
    GISegment("Colon",     pH=6.5, transit_time_h=18.0, surface_area_cm2=900.0,
              volume_mL=500.0, bile_factor=1.0),
]

# --- Fed-state GI parameters ---
FED_GI = [
    # Fed gastric pH and transit calibrated to match Bratel 1989 IR food effect
    # (observed Cmax ratio ~1.31). Standard breakfast: pH 3.0, GE time 0.65 h.
    # Rationale: pH 3.0 preserves weak-base dissolution in stomach; GE 0.65 h
    # is a standard-meal value (heavy meal would be 1.5 h).
    GISegment("Stomach",   pH=3.0, transit_time_h=0.65, surface_area_cm2=0.0,
              volume_mL=500.0, bile_factor=1.0, absorbs=False),
    GISegment("Duodenum",  pH=5.5, transit_time_h=0.25, surface_area_cm2=150.0,
              volume_mL=80.0,  bile_factor=15.0),
    GISegment("Jejunum1",  pH=6.2, transit_time_h=0.75, surface_area_cm2=5400.0,
              volume_mL=150.0, bile_factor=12.0),
    GISegment("Jejunum2",  pH=6.4, transit_time_h=0.75, surface_area_cm2=5400.0,
              volume_mL=150.0, bile_factor=10.0),
    GISegment("Ileum1",    pH=6.8, transit_time_h=1.00, surface_area_cm2=3600.0,
              volume_mL=100.0, bile_factor=4.0),
    GISegment("Ileum2",    pH=7.0, transit_time_h=1.00, surface_area_cm2=3600.0,
              volume_mL=100.0, bile_factor=3.0),
    GISegment("Colon",     pH=6.5, transit_time_h=18.0, surface_area_cm2=900.0,
              volume_mL=500.0, bile_factor=1.0),
]

# =============================================================================
# Helper functions
# =============================================================================
def compute_auc(time_h: np.ndarray, conc: np.ndarray) -> float:
    """Trapezoidal AUC (ng/mL * h)."""
    return float(np.trapz(conc, time_h))

def compute_cmax_tmax(time_h: np.ndarray, conc: np.ndarray):
    """Return (Cmax ng/mL, Tmax h)."""
    idx = np.argmax(conc)
    return float(conc[idx]), float(time_h[idx])
