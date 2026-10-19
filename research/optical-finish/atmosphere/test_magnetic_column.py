"""Independent analytic limits, not a Vela opacity validation."""
from pathlib import Path
import sys
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from magnetic_column import (Column, Parameters, response_opacities, cyclic_powers,
                             planck_nu, KAPPA_T, MP, KB, H, KEV_ERG)


def test_transverse_polarization_completeness():
    mu = np.linspace(0, 1, 21)
    q = cyclic_powers(mu).reshape(2, len(mu), 3)
    np.testing.assert_allclose(q.sum(-1), 1.0, atol=2e-16)
    np.testing.assert_allclose(q.sum(0)[:, 0], (1 + mu**2) / 2)
    np.testing.assert_allclose(q.sum(0)[:, 2], 1 - mu**2)


def test_low_energy_electron_proton_cancellation():
    # Corrected transverse absorption / parallel absorption tends to E^4
    # below *both* cyclotron energies. Electron-only formulas tend to E^2.
    e = np.array([1e-7, 2e-7, 4e-7])
    a, _ = response_opacities(e, 1e-12, 1e6, 3e12)
    suppression = a[:, 0] / a[:, 2]
    np.testing.assert_allclose(np.diff(np.log(suppression)) / np.diff(np.log(e)), 4, rtol=2e-5)


def test_reject_unvalidated_zero_field_limit():
    # The magnetic Gaunt factor convention and its ground-Landau limit are
    # not the nonmagnetic Kramers g_ff convention. Setting both to unity
    # produced a 0.412 ratio in an attempted B=0 comparison; retain that
    # failed audit in README rather than changing constants to force a match.
    with pytest.raises(ValueError, match='no validated B=0 limit'):
        response_opacities(0.2, 0.1, 1e6, 0)


def test_scattering_conservation_and_detailed_balance():
    column = Column(Parameters(depth_points=10, energy_points=4, angles_per_hemisphere=5))
    _, factors = response_opacities(0.1, 1.0, 1e6, 3e12)
    kernel = 1.5 * KAPPA_T * np.einsum('a,ca,ba,b->cb', factors, column.q, column.q, column.channel_w)
    expected_out = KAPPA_T * np.einsum('a,ca->c', factors, column.q)
    np.testing.assert_allclose(kernel.sum(-1), expected_out, atol=1e-16)
    weighted = kernel * column.channel_w[:, None]
    np.testing.assert_allclose(weighted, weighted.T, atol=1e-16)


def test_absorption_only_isothermal_slab_exact():
    p = Parameters(depth_points=12, energy_points=10, angles_per_hemisphere=4,
                   column_min=1e-8, column_max=100, magnetic_field_g=3e12)
    column = Column(p)
    t = np.full(p.depth_points, 7e5)
    out = column.transfer(t, True, scattering_multiplier=0, bottom_intensity_factor=0)
    optical_depth = np.sum(out['absorption'] * column.volume[None, :, None], axis=1)
    expected = planck_nu(column.nu[:, None], t[0]) / 2 * (
        -np.expm1(-optical_depth / column.channel_mu[None, :])
    )
    np.testing.assert_allclose(out['surface_i_nu_cgs'], expected, rtol=3e-13, atol=0)


def test_integral_transfer_energy_ledger():
    p = Parameters(depth_points=16, energy_points=16, angles_per_hemisphere=3,
                   column_min=1e-8, column_max=10)
    column = Column(p)
    t = np.geomspace(3e5, 3e6, p.depth_points)
    out = column.transfer(t, True)
    # F(bottom)-F(top) = 4 pi ∫rho kappa_abs (J-B/2) dz;
    # coherent scattering only exchanges modes/directions and cancels.
    expected_jump = 4 * np.pi * out['heating'] * column.volume
    actual_jump = np.diff(out['flux'])
    scale = max(np.max(np.abs(actual_jump)), 1e-100)
    np.testing.assert_allclose(actual_jump / scale, expected_jump / scale, atol=3e-13, rtol=3e-11)
    assert np.all(out['surface_i_nu_cgs'] >= 0)


def test_uniform_thermal_reservoir_absorptive_limit():
    p = Parameters(depth_points=10, energy_points=8, angles_per_hemisphere=3)
    column = Column(p)
    t = np.full(p.depth_points, 7e5)
    out = column.transfer(t, True, scattering_multiplier=0)
    expected = np.broadcast_to(planck_nu(column.nu[:, None], t[0]) / 2, out['surface_i_nu_cgs'].shape)
    np.testing.assert_allclose(out['surface_i_nu_cgs'], expected, rtol=4e-15)
