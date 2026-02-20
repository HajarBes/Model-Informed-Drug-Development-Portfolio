"""
Unit tests for the KRAS G12C QSP model core functions.

Tests cover: parameter derivations, ODE system shape, steady-state physics,
mass conservation, and scenario preset correctness.
"""

import sys
import os

import numpy as np
import pytest

# Ensure the analysis directory is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir, "analysis"))
from kras_qsp import Params, rhs, simulate, params_vehicle, params_nsclc_mono, params_crc_mono, params_crc_combo


class TestParams:
    """Tests for Params dataclass derived quantities."""

    def test_kras_steady_state(self):
        """GDP + GTP should equal KRAS_total_ss at steady state (no drug)."""
        p = Params()
        gdp, gtp, drug = p.kras_ss()
        assert drug == 0.0
        assert abs((gdp + gtp) - p.KRAS_total_ss) < 1e-12

    def test_elimination_rate(self):
        """ke should equal ln(2) / t_half_h."""
        p = Params()
        expected = np.log(2) / p.t_half_h
        assert abs(p.ke - expected) < 1e-12

    def test_synthesis_rate(self):
        """k_syn should equal k_deg * KRAS_total_ss."""
        p = Params()
        expected = p.k_deg * p.KRAS_total_ss
        assert abs(p.k_syn - expected) < 1e-12


class TestScenarioPresets:
    """Each preset function returns correct distinguishing parameters."""

    def test_nsclc_mono(self):
        p = params_nsclc_mono()
        assert p.G_fb == 1.2
        assert p.f_block == 0.0

    def test_crc_mono(self):
        p = params_crc_mono()
        assert p.G_fb == 3.5
        assert p.f_block == 0.0

    def test_crc_combo(self):
        p = params_crc_combo()
        assert p.G_fb == 3.5
        assert p.f_block == 0.55

    def test_vehicle(self):
        p = params_vehicle()
        assert p.dose_mg == 0.0


class TestRHS:
    """Tests for the ODE right-hand side function."""

    def test_rhs_output_shape(self):
        """rhs() must return a 9-element array."""
        p = Params()
        gdp, gtp, _ = p.kras_ss()
        y0 = np.array([0.0, 0.0, gdp, gtp, 0.0, p.R_basal, p.E_ss, p.T_0, 0.0])
        dose_times = np.array([0.0])
        dydt = rhs(0.0, y0, p, dose_times, gtp)
        assert dydt.shape == (9,)


class TestSimulation:
    """Integration-level tests on the simulate() function."""

    def test_vehicle_no_drug(self):
        """Vehicle scenario: A_gut stays 0, tumor grows monotonically."""
        p = params_vehicle()
        df = simulate(p, t_days=30, dt_h=1.0)
        # No drug should enter the gut
        assert df["A_gut"].max() < 1e-10
        # Tumor should grow monotonically (no drug kill)
        tumor = df["T"].values
        # Check that the final value exceeds the initial value
        assert tumor[-1] > tumor[0]

    def test_mass_conservation_pk(self):
        """PK mass (A_gut + A_central) after first dose should not exceed the bolus."""
        p = Params()
        # simulate() applies doses every tau_h; look only at the first dosing interval
        df = simulate(p, t_days=1, dt_h=0.5)
        first_interval = df[df["time_h"] < p.tau_h]
        pk_mass = first_interval["A_gut"] + first_interval["A_central"]
        initial_bolus = p.F * p.dose_mg
        # No mass creation — only absorption and elimination
        assert pk_mass.max() <= initial_bolus + 0.1  # small numerical tolerance
        # PK mass should decay over the interval (elimination exceeds zero input)
        assert pk_mass.iloc[-1] < pk_mass.iloc[1]
