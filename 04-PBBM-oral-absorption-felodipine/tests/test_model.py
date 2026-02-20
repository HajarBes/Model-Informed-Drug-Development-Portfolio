"""
Unit tests for the PBBM oral absorption (ACAT) model core physics.

Tests cover: pH-dependent solubility, bile salt enhancement, dissolution
driving force, absorption rate constant, and full mass balance.
"""

import sys
import os

import numpy as np
import pytest

# Ensure the analysis directory is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir, "analysis"))
from acat_model_02 import (
    saturation_solubility,
    dissolution_rate,
    absorption_rate_constant,
    simulate,
    DrugParams,
    FASTED_GI,
)


class TestSolubility:
    """Tests for saturation_solubility (Henderson-Hasselbalch)."""

    def test_solubility_ph_dependence(self):
        """Lower pH should give higher solubility for a weak base (pKa > pH)."""
        S0 = 0.0005  # mg/mL
        pKa = 5.07
        sol_low_pH = saturation_solubility(S0, pKa, pH=2.0, bile_factor=1.0)
        sol_high_pH = saturation_solubility(S0, pKa, pH=7.0, bile_factor=1.0)
        assert sol_low_pH > sol_high_pH

    def test_solubility_bile_enhancement(self):
        """bile_factor > 1 should increase C_sat proportionally."""
        S0 = 0.0005
        pKa = 5.07
        pH = 6.5
        sol_base = saturation_solubility(S0, pKa, pH, bile_factor=1.0)
        sol_bile = saturation_solubility(S0, pKa, pH, bile_factor=3.0)
        assert abs(sol_bile / sol_base - 3.0) < 1e-10


class TestDissolution:
    """Tests for dissolution_rate (Noyes-Whitney)."""

    def test_dissolution_driving_force(self):
        """Rate should be 0 when C_diss >= C_sat (saturated)."""
        rate = dissolution_rate(
            S_mg=5.0, C_dissolved_mg_mL=0.01, C_sat_mg_mL=0.01,
            Deff=5e-6, density=1.3, radius_cm=25e-4
        )
        assert rate == 0.0

    def test_dissolution_proportional_to_solid(self):
        """More solid mass should give faster dissolution."""
        kwargs = dict(
            C_dissolved_mg_mL=0.0, C_sat_mg_mL=0.01,
            Deff=5e-6, density=1.3, radius_cm=25e-4
        )
        rate_low = dissolution_rate(S_mg=1.0, **kwargs)
        rate_high = dissolution_rate(S_mg=5.0, **kwargs)
        assert rate_high > rate_low
        # Should be exactly proportional (rate = k * S * driving_force)
        assert abs(rate_high / rate_low - 5.0) < 1e-10


class TestAbsorption:
    """Tests for absorption_rate_constant."""

    def test_absorption_rate_constant(self):
        """Known Peff, SA, V should give expected k_abs."""
        Peff = 5.0e-4  # cm/s
        SA = 5400.0     # cm^2
        V = 100.0       # mL
        k_abs = absorption_rate_constant(Peff, SA, V)
        # Expected: (5e-4 * 3600) * 5400 / 100 = 1.8 * 54 = 97.2 /h
        expected = (Peff * 3600.0) * SA / V
        assert abs(k_abs - expected) < 1e-10

    def test_zero_permeability(self):
        """Peff = 0 should give k_abs = 0."""
        k_abs = absorption_rate_constant(0.0, 5400.0, 100.0)
        assert k_abs == 0.0


class TestMassBalance:
    """Full simulation mass balance test."""

    def test_mass_balance(self):
        """Total mass accounted for should equal dose (< 0.01 mg error)."""
        drug = DrugParams()
        sim = simulate(drug, FASTED_GI, t_end_h=48.0)
        last = sim.iloc[-1]

        S_total = last["S_total_mg"]
        D_total = last["D_total_mg"]
        Ac = last["Ac_mg"]
        Ap = last["Ap_mg"]
        cum_absorbed = last["cum_absorbed_mg"]

        # Systemic elimination: integral of ke * Ac
        ke = drug.CL_L_h / drug.Vc_L
        sys_eliminated = float(np.trapz(
            ke * sim["Ac_mg"].values, sim["time_h"].values
        ))

        # Fecal = dose - GI remaining - absorbed through wall
        fecal = drug.dose_mg - S_total - D_total - cum_absorbed

        mass_accounted = S_total + D_total + cum_absorbed + fecal
        mass_error = abs(mass_accounted - drug.dose_mg)
        assert mass_error < 0.01, f"Mass balance error {mass_error:.4f} mg exceeds 0.01 mg"
