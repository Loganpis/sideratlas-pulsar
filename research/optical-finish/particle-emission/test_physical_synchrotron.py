"""Physics checks; these do not validate Vela injection or atmosphere."""
import numpy as np
import pytest
from scipy.constants import c, e, epsilon_0, m_e
from scipy.integrate import quad

from physical_synchrotron import (
    PairModel, angular_synchrotron_pnu, harmonic_line_power_w_sr,
    binned_harmonic_spectrum, synchrotron_power, column_transport,
    boost_to_magnetic_frame,
    generate_events,
)
from load_camera_events import load_camera_events


@pytest.mark.parametrize('gamma,pitch,nmax', [
    (1.005, .5, 12), (1.05, 1.2, 30), (2., .5, 200),
    (5., .3, 500), (10., .1, 500),
])
def test_exact_harmonics_integrate_to_lienard(gamma, pitch, nmax):
    def angular_sum(cosine):
        _, power = harmonic_line_power_w_sr(gamma, 1., pitch, np.arccos(cosine), np.arange(1, nmax+1))
        return 2*np.pi*power.sum()
    measured = quad(angular_sum, -1, 1, epsabs=1e-100, epsrel=2e-8, limit=200)[0]
    expected = synchrotron_power(gamma, 1., pitch)
    assert measured == pytest.approx(expected, rel=1e-7)


def test_nonrelativistic_first_harmonic_angular_pattern():
    beta = 1e-4
    gamma = 1/np.sqrt(1-beta*beta)
    theta = np.array([0., .3, 1.1, np.pi])
    frequency, power = harmonic_line_power_w_sr(gamma, 1., np.pi/2, theta, 1)
    omega = e/(gamma*m_e)
    expected = e*e*omega*omega*beta*beta/(32*np.pi*np.pi*epsilon_0*c)*(1+np.cos(theta)**2)
    np.testing.assert_allclose(power, expected, rtol=5e-8)
    np.testing.assert_allclose(frequency, omega/(2*np.pi), rtol=1e-12)


def test_schott_source_time_power_matches_direct_lienard_wiechert():
    # Independent time-domain angular radiation expression, including source
    # denominator^5. This catches an accidental observer-time factor in Schott.
    gamma, field, pitch, viewing = 2., 1., .5, .7
    speed = np.sqrt(1-1/gamma**2)
    bp, bt = speed*np.cos(pitch), speed*np.sin(pitch)
    omega = e*field/(gamma*m_e)
    n = np.array([np.sin(viewing), 0., np.cos(viewing)])
    def source_power(phi):
        beta = np.array([bt*np.cos(phi), bt*np.sin(phi), bp])
        acceleration = omega*np.array([-bt*np.sin(phi), bt*np.cos(phi), 0.])
        numerator = np.cross(n, np.cross(n-beta, acceleration))
        return e*e/(16*np.pi*np.pi*epsilon_0*c)*(numerator@numerator)/(1-n@beta)**5
    direct = quad(source_power, 0, 2*np.pi, epsabs=1e-100, epsrel=1e-10)[0]/(2*np.pi)
    _, lines = harmonic_line_power_w_sr(gamma, field, pitch, viewing, np.arange(1, 301))
    np.testing.assert_allclose(lines.sum(), direct, rtol=1e-9)


def test_bin_partition_preserves_exact_line_powers():
    args = (np.array([20.]), np.array([120.]), np.array([.15]),
            np.array([.13]), np.ones(1), np.ones(1))
    coarse, _, _ = binned_harmonic_spectrum(*args, np.linspace(400e-9, 700e-9, 5))
    fine, _, _ = binned_harmonic_spectrum(*args, np.linspace(400e-9, 700e-9, 101))
    np.testing.assert_allclose(coarse.sum(axis=1), fine.sum(axis=1), rtol=1e-14)


def test_numerical_harmonic_truncation_bounds_actual_discarded_power():
    # Deliberately loose tolerance exercises the bound with measurable drops.
    viewing = np.linspace(.02, 3.10, 17)
    gamma = np.full_like(viewing, 5.)
    field = np.ones_like(viewing)
    pitch = np.full_like(viewing, .2)
    edge = np.geomspace(1e-5, .01, 25)
    args = (gamma, field, pitch, viewing, np.ones_like(viewing), np.ones_like(viewing), edge)
    exact, _, _ = binned_harmonic_spectrum(*args, negligible_power_bound_w_sr=0.)
    limited, _, diagnostics = binned_harmonic_spectrum(*args, negligible_power_bound_w_sr=1e-15)
    dropped = exact.sum(axis=1)-limited.sum(axis=1)
    bound = diagnostics['omitted_power_bounds_w_sr']
    assert np.all(dropped <= bound*(1+1e-10)+1e-100)
    assert np.all(dropped >= -1e-100)


def test_line_exact_xyz_and_boost_jacobian():
    gamma, b, pitch, theta = 10., 100., .2, .22
    nu1, power1 = harmonic_line_power_w_sr(gamma, b, pitch, theta, 1)
    wavelength = c/nu1
    edges = np.array([wavelength*.95, wavelength*1.05])
    cie_wave = np.array([wavelength*.9, wavelength*1.1])
    cie_xyz = np.tile([.2, .5, .7], (2, 1))
    bands, xyz, _ = binned_harmonic_spectrum(
        np.array([gamma]), np.array([b]), np.array([pitch]), np.array([theta]),
        np.array([1.]), np.array([1.]), edges, cie_wave, cie_xyz)
    np.testing.assert_allclose(bands[0, 0], power1, rtol=1e-12)
    np.testing.assert_allclose(xyz[0], 683*power1*np.array([.2, .5, .7]), rtol=1e-12)
    boosted, _, _ = binned_harmonic_spectrum(
        np.array([gamma]), np.array([b]), np.array([pitch]), np.array([theta]),
        np.array([2.]), np.array([1.25]), edges/2)
    np.testing.assert_allclose(boosted.sum(), power1*2**3/1.25, rtol=1e-12)


def test_zero_pitch_emits_no_cyclotron_power():
    _, power = harmonic_line_power_w_sr(10., 1., 0., .4, np.arange(1, 10))
    assert np.all(power == 0)


def test_boost_direction_and_frequency():
    n = np.array([0., 0., 1.])
    transformed, d, g = boost_to_magnetic_frame(n, np.array([[0., 0., 0.], [0., 0., .2], [.2, 0., 0.]]))
    np.testing.assert_allclose(np.linalg.norm(transformed, axis=1), 1., atol=1e-14)
    np.testing.assert_allclose(transformed[0], n, atol=1e-14)
    assert d[1] == pytest.approx(np.sqrt(1.2/.8))
    assert d[2] == pytest.approx(np.sqrt(1-.2**2))


def test_cooling_column_energy_and_adiabatic_scale():
    result = column_transport(PairModel(), .12, .85)
    diag = result['diagnostics']
    assert diag['relative_energy_ledger_error'] < 1e-8
    assert diag['remaining_perpendicular_energy_bound_j'] < 2e-8*diag['emitted_energy_j']
    assert diag['max_gyro_radius_over_b_gradient_length'] < 1e-8
    assert diag['max_gyroperiod_over_cooling_time'] < 1e-3
    assert diag['max_quantum_chi'] < 1e-3
    assert result['residence_time_weights_s'].sum() == pytest.approx(diag['duration_s'], rel=1e-10)
    assert diag['final_transverse_momentum'] == pytest.approx(.001, rel=1e-7)


def test_continuum_kernel_angular_recovery_at_high_harmonics():
    gamma, b, pitch = 100., 1., .5
    nuc = 3*e*b*gamma*gamma*np.sin(pitch)/(4*np.pi*m_e)
    nu = nuc
    def f(theta):
        return float(angular_synchrotron_pnu(gamma, b, pitch, theta, nu))*2*np.pi*np.sin(theta)
    measured = quad(f, pitch-.2, pitch+.2, epsabs=1e-100, epsrel=1e-7)[0]
    from scipy.special import kv
    F = quad(lambda u: kv(5/3, u), 1., np.inf)[0]
    expected = np.sqrt(3)*e**3*b*np.sin(pitch)/(4*np.pi*epsilon_0*m_e*c)*F
    assert measured == pytest.approx(expected, rel=2e-5)


def test_refuse_superluminal_camera_boost_and_bad_harmonic():
    with pytest.raises(ValueError):
        boost_to_magnetic_frame([1., 0., 0.], np.array([1., 0., 0.]))
    with pytest.raises(ValueError):
        harmonic_line_power_w_sr(20., 1., .2, .3, 1.5)


def test_compact_camera_reconstructs_geometry_with_measured_color_error(tmp_path):
    model = PairModel(phase_nodes=8, azimuth_nodes=8, cooling_nodes=4, radial_nodes=1, wavelength_nodes=4)
    full = tmp_path/'full.npz'
    compact = tmp_path/'compact.npz'
    generate_events(model, full, tmp_path/'full.json')
    report = generate_events(model, compact, tmp_path/'compact.json', compact_camera=True)
    df, dc = load_camera_events(full), load_camera_events(compact)
    assert len(df['positions_m']) > 0
    for key in ('positions_m', 'source_time_s', 'source_interval_s', 'number_weights', 'pattern_beta'):
        np.testing.assert_allclose(dc[key], df[key], rtol=0, atol=1e-12 if key == 'positions_m' else 0)
    np.testing.assert_allclose(dc['xyz_lm_sr'], df['xyz_lm_sr'], rtol=6e-8, atol=1e-44)
    assert max(report['xyz_float32_storage_relative_energy_error_bound']) < 6e-8
    np.testing.assert_allclose(dc['event_band_power_w_sr'], df['power_bins_w_sr'].sum(axis=1), rtol=0, atol=0)
