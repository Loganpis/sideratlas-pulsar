"""Physical, directional synchrotron column model; not a Vela reconstruction.

SI units throughout. The radiation is incoherent and gyro-averaged. Physical
particle energies and magnetic fields are specified independently of reduced
PIC output. Production photons use the exact classical Schott harmonic
angular spectrum (Marcowith & Malzac 2003, eqs 1,2,6). A high-harmonic continuum
kernel (Lai et al. 2025, eqs 2--4) is retained for benchmarks only. Neither
kernel is an angle-integrated spectrum multiplied by a beaming function.

The spatial model is deliberately narrower than a PIC claim: guiding centres
follow a prescribed near-zone dipole, with adiabatic focusing and classical
synchrotron cooling. Pitch activation, pair supply, dipole geometry and rigid
corotation are prescribed. No cascade, resonant radio absorption, current
closure, kinetic feedback, atmosphere, GR or extinction is calculated here.
"""
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.constants import c, e, epsilon_0, m_e, hbar
from scipy.integrate import solve_ivp
from scipy.special import kv, jv, jvp

B_Q_T = m_e**2*c**2/(e*hbar)
RR_COEFFICIENT = e**4/(6*np.pi*epsilon_0*m_e**3*c**3)


def _positive(value, name, allow_zero=False):
    a = np.asarray(value, float)
    if np.any(~np.isfinite(a)) or np.any(a < 0 if allow_zero else a <= 0):
        raise ValueError(f'{name} must be finite and positive')
    return a


def angular_synchrotron_pnu(gamma, magnetic_field_t, pitch_rad,
                            observer_field_angle_rad, frequency_hz):
    """Gyro-averaged emitted P_nu/Omega, W Hz^-1 sr^-1 per lepton.

    Local magnetic frame, E'=0. Gyrophases are uniformly distributed and the
    population is steady on a gyroperiod. theta=observer-B angle minus pitch.
    The pulse-energy kernel is multiplied ONCE by eB/(2 pi gamma m_e).
    This is power per source-frame time, not the power of a translating single
    perfect helix per compressed observer time. Source-cell translation and
    arrival-time mapping are handled separately.

    This is the high-harmonic small-formation-angle continuum, not exact
    cyclotron harmonics: gamma>=10, gamma*sin(pitch)>=10, chi<=0.1 are enforced.
    The regime cutoffs are guards, not precision guarantees. The gamma^-2
    curvature kernel uses beta=1 at leading ultrarelativistic order.
    """
    g, b, pitch, viewing, nu = np.broadcast_arrays(
        np.asarray(gamma, float), np.asarray(magnetic_field_t, float),
        np.asarray(pitch_rad, float), np.asarray(observer_field_angle_rad, float),
        np.asarray(frequency_hz, float))
    for a in (g, b, pitch, viewing, nu):
        if np.any(~np.isfinite(a)):
            raise ValueError('All radiation inputs must be finite')
    if (np.any(g < 10) or np.any(b <= 0) or np.any(nu <= 0)
            or np.any(pitch <= 0) or np.any(pitch >= np.pi)
            or np.any(viewing < 0) or np.any(viewing > np.pi)):
        raise ValueError('Require gamma>=10, B,nu>0, pitch in (0,pi), viewing in [0,pi]')
    if np.any(g*np.sin(pitch) < 10*(1-1e-12)):
        raise ValueError('High-harmonic continuum requires gamma*sin(pitch)>=10')
    if np.any(g*b*np.sin(pitch)/B_Q_T > .1):
        raise ValueError('Classical synchrotron invalid for chi_e>0.1')
    omega_b = e*b/(g*m_e)
    rho = c/(omega_b*np.sin(pitch))
    omega = 2*np.pi*nu
    theta = viewing-pitch
    a = 1/g**2+theta**2
    xi = omega*rho/(3*c)*a**1.5
    perpendicular = a*a*kv(2/3, xi)**2
    parallel = theta*theta*a*kv(1/3, xi)**2
    pulse_per_omega = e*e/(12*np.pi**3*epsilon_0*c)*(omega*rho/c)**2
    # d omega/d nu=2*pi; omega_b/(2*pi) independent pulses per source second.
    return pulse_per_omega*(perpendicular+parallel)*omega_b


def synchrotron_power(gamma, magnetic_field_t, pitch_rad):
    """Exact classical total power for supplied physical fixed-pitch momenta."""
    g, b, a = np.broadcast_arrays(gamma, magnetic_field_t, pitch_rad)
    return e**4*b*b*(g*g-1)*np.sin(a)**2/(6*np.pi*epsilon_0*m_e*m_e*c)


def harmonic_line_power_w_sr(gamma, magnetic_field_t, pitch_rad,
                             observer_field_angle_rad, harmonic):
    """Exact classical Schott line frequency and SOURCE-time power per sr.

    Eq.1,2,6 Marcowith & Malzac2003, doi:10.1051/0004-6361:20031060,
    converted to SI and integrated analytically across each delta line.
    A=1-beta_parallel*cos(theta), nu_n=n*eB/(2pi*gamma*m_e*A).
    The emitted source-time power has A^-3; a single translating perfect helix
    would acquire another arrival-time compression factor, not included here.
    No artificial linewidth is assigned. Photons/electron intensity is also
    the positron intensity. gamma need only exceed1, not gamma*sin(pitch)>10.
    """
    g, b, a, theta, order = np.broadcast_arrays(
        np.asarray(gamma, float), np.asarray(magnetic_field_t, float),
        np.asarray(pitch_rad, float), np.asarray(observer_field_angle_rad, float),
        np.asarray(harmonic))
    if (np.any(~np.isfinite(g)) or np.any(g <= 1) or np.any(~np.isfinite(b))
            or np.any(b <= 0) or np.any(~np.isfinite(a)) or np.any(a < 0)
            or np.any(a > np.pi) or np.any(~np.isfinite(theta))
            or np.any(theta < 0) or np.any(theta > np.pi)
            or np.any(~np.isfinite(order)) or np.any(order < 1)
            or np.any(order != np.floor(order))):
        raise ValueError('Finite gamma>1,B>0,physical angles,positive integer harmonics required')
    if np.any(g*b*np.sin(a)/B_Q_T > .1):
        raise ValueError('Classical harmonics invalid for chi_e>0.1')
    speed = np.sqrt((g-1)*(g+1))/g
    bp, bt = speed*np.cos(a), speed*np.sin(a)
    co, si = np.cos(theta), np.sin(theta)
    # Stable in the forward direction for an almost longitudinal particle.
    den = (1-bp)+bp*(1-co)
    nub = e*b/(2*np.pi*g*m_e)
    frequency = order*nub/den
    argument = order*bt*si/den
    j = jv(order, argument)
    derivative = jvp(order, argument)
    along = np.divide((co-bp)*j, si, out=np.zeros_like(j), where=np.abs(si) > 1e-14)
    # At either magnetic axis, only n=1 survives; J1(z)/sin(theta)->n bt/(2A).
    axis = np.abs(si) <= 1e-14
    along = np.where(axis & (order == 1), (co-bp)*bt/(2*den), along)
    line_power = e*e/(2*epsilon_0*c)*nub*nub*order*order/den**3*(along*along+bt*bt*derivative*derivative)
    return frequency, line_power


def binned_harmonic_spectrum(gamma, magnetic_field_t, pitch_rad, viewing,
                             doppler, bulk_gamma, wavelength_edges_m,
                             cie_wavelength_m=None, cie_xyz=None,
                             max_harmonic=200_000, negligible_power_bound_w_sr=1e-45):
    """Exact delta-line bin powers, W/sr/lepton, and optional line-exact XYZ.

    Inputs describe N steady comoving magnetic cells. For their local source
    P'_nu', source lab P_nu=D²/Gamma P'_nu'(nu/D); a delta-line's integrated
    power therefore transforms as D³/Gamma. No line is broadened physically.
    Spectral bins are for reporting; optional CIE interpolation is evaluated
    at each individual line wavelength, before grouping it into a bin.
    """
    g, b, a, theta, d, gb = np.broadcast_arrays(gamma, magnetic_field_t, pitch_rad, viewing, doppler, bulk_gamma)
    if g.ndim != 1:
        raise ValueError('One-dimensional cell arrays required')
    edges = _positive(wavelength_edges_m, 'wavelength_edges_m')
    if edges.ndim != 1 or np.any(np.diff(edges) <= 0):
        raise ValueError('Increasing wavelength edges required')
    if np.any(np.asarray(d) <= 0) or np.any(np.asarray(gb) < 1):
        raise ValueError('Positive Doppler and bulk gamma>=1 required')
    speed = np.sqrt((g-1)*(g+1))/g
    bp, bt = speed*np.cos(a), speed*np.sin(a)
    co, si = np.cos(theta), np.sin(theta)
    den = (1-bp)+bp*(1-co)
    fundamental_lab = d*e*b/(2*np.pi*g*m_e*den)
    lower, upper = c/edges[-1], c/edges[0]
    nmin = np.maximum(1, np.ceil(lower/fundamental_lab).astype(np.int64))
    nmax = np.floor(upper/fundamental_lab).astype(np.int64)
    if np.any(nmax > max_harmonic):
        raise ValueError('Harmonic work limit exceeded; increase explicit max_harmonic')
    counts = np.maximum(0, nmax-nmin+1)
    # Certified Cauchy generating-function bound for real z<n:
    # |J_n(z)| <= exp(sqrt(n²-z²)-n*acosh(n/z)). This follows by taking
    # the Laurent coefficient of exp[(z/2)(t-1/t)] on |t|=r, minimizing
    # r^-n exp[z(r-1/r)/2]. For z>=n the bound1 suffices here.
    # J'_n=(J_{n-1}-J_{n+1})/2. The n-1 bound dominates the n+1 bound and
    # decreases with n at fixed q=z/n<1. Using that larger neighbouring
    # bound and sum(n²) therefore bounds the ENTIRE band tail.
    def bessel_bound(order, argument):
        order, argument = np.broadcast_arrays(np.asarray(order, float), np.asarray(argument, float))
        result = np.ones_like(order)
        valid = (order > 0) & (argument > 0) & (argument < order)
        oo, zz = order[valid], argument[valid]
        result[valid] = np.exp(np.sqrt((oo-zz)*(oo+zz))-oo*np.arccosh(oo/zz))
        result[(argument == 0) & (order > 0)] = 0
        return result
    first_argument = nmin*bt*si/den
    jb = bessel_bound(nmin, first_argument)
    jdb = bessel_bound(nmin-1, first_argument)
    ratio = np.divide(co-bp, si, out=np.full_like(si, np.inf), where=np.abs(si) > 1e-14)
    sum_order_squared = (nmax.astype(float)*(nmax+1)*(2*nmax+1)-(nmin-1).astype(float)*nmin*(2*nmin-1))/6
    sum_order_squared = np.maximum(0., sum_order_squared)
    base = e*e/(2*epsilon_0*c)*(e*b/(2*np.pi*g*m_e))**2/den**3
    with np.errstate(invalid='ignore'):
        upper_bound = base*sum_order_squared*(ratio*ratio*jb*jb+bt*bt*jdb*jdb)*d**3/gb
    upper_bound = np.where(np.abs(si) <= 1e-14, np.inf, upper_bound)
    negligible = upper_bound < negligible_power_bound_w_sr
    omitted_bound = np.where(negligible, upper_bound, 0.)
    counts = np.where(negligible, 0, counts)
    # For fixed argument/order ratio<1, J_n(nq) decays at high n. Eliminate
    # only rows already in numerical underflow at the FIRST line; no signal
    # threshold or arbitrary physical beaming cone is introduced.
    q = bt*si/den
    first_j = jv(nmin, nmin*q)
    first_prev = jv(np.maximum(0, nmin-1), nmin*q)
    numerically_nonzero = (first_j != 0) | (first_prev != 0)
    counts = np.where(numerically_nonzero, counts, 0)
    cell = np.repeat(np.arange(len(g)), counts)
    if len(cell):
        starts = np.repeat(np.cumsum(counts)-counts, counts)
        order = np.repeat(nmin, counts)+np.arange(len(cell))-starts
        _, power_prime = harmonic_line_power_w_sr(g[cell], b[cell], a[cell], theta[cell], order)
        wave = c/(fundamental_lab[cell]*order)
        power = power_prime*d[cell]**3/gb[cell]
        band = np.searchsorted(edges, wave, side='right')-1
        band = np.clip(band, 0, len(edges)-2)
    else:
        order, wave, power, band = (np.array([], dtype=float) for _ in range(4))
    spectral = np.zeros((len(g), len(edges)-1))
    if len(cell):
        np.add.at(spectral, (cell, band.astype(int)), power)
    xyz = np.zeros((len(g), 3))
    if cie_wavelength_m is not None:
        if cie_xyz is None:
            raise ValueError('Supply CIE wavelengths and all three matching functions')
        for k in range(3):
            value = np.interp(wave, cie_wavelength_m, cie_xyz[:, k], left=0, right=0)
            np.add.at(xyz[:, k], cell, 683*power*value)
    return spectral, xyz, dict(harmonics_evaluated=len(cell),
                               max_intersecting_harmonic=int(np.max(nmax)),
                               numerically_nonzero_cells=int(np.count_nonzero(spectral.sum(axis=1))),
                               negligible_power_bound_per_cell_w_sr=negligible_power_bound_w_sr,
                               sum_omitted_power_bounds_w_sr_per_lepton_samples=float(omitted_bound.sum()),
                               omitted_power_bounds_w_sr=omitted_bound)


def rotate_y(x, angle):
    x = np.asarray(x)
    co, si = np.cos(angle), np.sin(angle)
    return np.stack((co*x[..., 0]+si*x[..., 2], x[..., 1],
                     -si*x[..., 0]+co*x[..., 2]), axis=-1)


def rotate_z(x, angle):
    x = np.asarray(x)
    co, si = np.cos(angle), np.sin(angle)
    return np.stack((co*x[..., 0]-si*x[..., 1],
                     si*x[..., 0]+co*x[..., 1], x[..., 2]), axis=-1)


def dipole_b_t(position_m, polar_field_t, star_radius_m):
    """Static dipole in its magnetic coordinates; polar field convention."""
    x = np.asarray(position_m, float)
    r = np.linalg.norm(x, axis=-1)
    if np.any(r < star_radius_m):
        raise ValueError('Dipole transport is restricted to outside the stellar surface')
    n = x/r[..., None]
    z = np.zeros_like(x)
    z[..., 2] = 1
    return .5*polar_field_t*(star_radius_m/r)[..., None]**3*(3*n*n[..., 2, None]-z)


def grad_log_dipole_b(position_m):
    x = np.asarray(position_m, float)
    r = np.linalg.norm(x)
    n = x/r
    co = n[2]
    return -3*n/r + 3*co/(1+3*co*co)*(np.array([0., 0., 1.])-co*n)/r


@dataclass(frozen=True)
class PairModel:
    period_s: float = .08937108859
    star_radius_m: float = 13_000.
    polar_field_t: float = 3e8
    inclination_deg: float = 75.
    observer_deg: float = 65.
    pair_multiplicity_per_primary: float = 6000.
    initial_gamma: float = 300.
    initial_pitch_rad: float = .10
    final_transverse_momentum: float = .001
    activation_radius_min_rlc: float = .08
    activation_radius_max_rlc: float = .18
    max_radius_rlc: float = .20
    annulus_inner: float = .80
    annulus_outer: float = .90
    radial_nodes: int = 3
    annulus_nodes: int = 1
    azimuth_nodes: int = 64
    cooling_nodes: int = 32
    phase_nodes: int = 128
    wavelength_nodes: int = 31
    max_cooling_time_s: float = .01
    aberration: bool = True

    def __post_init__(self):
        vals = asdict(self)
        if any(isinstance(v, float) and not np.isfinite(v) for v in vals.values()):
            raise ValueError('Finite parameters required')
        if not (self.period_s > 0 and self.star_radius_m > 0 and self.polar_field_t > 0
                and self.initial_gamma >= 10 and 0 < self.initial_pitch_rad < np.pi/2
                and self.initial_gamma*np.sin(self.initial_pitch_rad) > self.final_transverse_momentum > 0
                and 0 < self.activation_radius_min_rlc <= self.activation_radius_max_rlc < self.max_radius_rlc <= .2
                and self.star_radius_m/self.rlc_m < self.activation_radius_min_rlc
                and 0 < self.annulus_inner < self.annulus_outer <= 1
                and 0 <= self.inclination_deg <= 90 and 0 < self.observer_deg < 180
                and self.pair_multiplicity_per_primary > 0
                and self.max_cooling_time_s > 0):
            raise ValueError('Physical/model domain violation')
        for name in ('radial_nodes', 'annulus_nodes', 'azimuth_nodes', 'cooling_nodes', 'phase_nodes', 'wavelength_nodes'):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < (2 if name in ('wavelength_nodes', 'phase_nodes') else 1):
                raise ValueError(f'Positive integral {name} required')

    @property
    def omega(self):
        return 2*np.pi/self.period_s

    @property
    def rlc_m(self):
        return c/self.omega

    @property
    def observer_direction(self):
        z = np.deg2rad(self.observer_deg)
        return np.array([np.sin(z), 0., np.cos(z)])


def column_transport(model, activation_radius_rlc, footprint_ovc, pole=1,
                     rtol=2e-9, atol=1e-11):
    """Cooling/focusing calculation, independent of reduced-PIC energies.

    The activation layer gives existing cascade pairs a chosen pitch. The
    omitted radio process is not silently solved: both activation location and
    initial momenta are explicit. Momentum is m_e*c units, time is local static
    magnetic-frame seconds. No acceleration E_parallel, inverse Compton,
    radio reheating, curvature-radiation loss or pair creation is supplied.
    Transport continues through low harmonics to an explicit small residual
    perpendicular momentum; the final residual energy is bounded in the report.
    """
    if pole not in (-1, 1):
        raise ValueError('Pole must be +1 or -1')
    r = activation_radius_rlc*model.rlc_m
    re = model.rlc_m/footprint_ovc**2
    theta = np.arcsin(np.sqrt(r/re))
    p = np.sqrt(model.initial_gamma**2-1)
    y0 = np.array([r*np.sin(theta), 0., pole*r*np.cos(theta),
                   p*np.cos(model.initial_pitch_rad), p*np.sin(model.initial_pitch_rad), 0.])

    def rhs(t, y):
        x, par, perp = y[:3], y[3], y[4]
        b = dipole_b_t(x, model.polar_field_t, model.star_radius_m)
        bm = np.linalg.norm(b)
        direction = pole*b/bm
        g = np.sqrt(1+par*par+perp*perp)
        velocity = c*par/g*direction
        dlogb = grad_log_dipole_b(x)@velocity
        # Magnetic-moment focusing conserves total kinetic energy on its own.
        focus_perp = .5*perp*dlogb
        focus_par = -.5*perp*perp/par*dlogb
        a = RR_COEFFICIENT*bm*bm
        loss_par = -a*par*perp*perp/g
        loss_perp = -a*perp*(1+perp*perp)/g
        # Lienard energy loss integrated as a separate ledger variable, in J.
        energy_rate = m_e*c*c*a*perp*perp
        return np.r_[velocity, focus_par+loss_par, focus_perp+loss_perp, energy_rate]

    def cooling_end(t, y):
        return y[4]-model.final_transverse_momentum
    cooling_end.terminal = True
    cooling_end.direction = -1

    def radial_end(t, y):
        return model.max_radius_rlc*model.rlc_m-np.linalg.norm(y[:3])
    radial_end.terminal = True
    radial_end.direction = -1
    sol = solve_ivp(rhs, [0, model.max_cooling_time_s], y0, method='DOP853',
                    rtol=rtol, atol=atol, dense_output=True,
                    events=(cooling_end, radial_end))
    if not sol.success or sol.t[-1] == model.max_cooling_time_s:
        raise RuntimeError('Column failed to reach a declared physical/continuum boundary')
    # Log-time Gaussian quadrature resolves the initial bright, rapid cooling
    # and the long weak low-harmonic tail. Weights remain source seconds.
    nodes, weights = np.polynomial.legendre.leggauss(model.cooling_nodes)
    duration = sol.t[-1]
    initial_power = m_e*c*c*RR_COEFFICIENT*np.linalg.norm(dipole_b_t(y0[:3], model.polar_field_t, model.star_radius_m))**2*y0[4]**2
    tscale = .01*(model.initial_gamma-1)*m_e*c*c/initial_power
    logspan = np.log1p(duration/tscale)
    s = .5*logspan*(nodes+1)
    time = tscale*np.expm1(s)
    dt = tscale*np.exp(s)*.5*logspan*weights
    y = sol.sol(time).T
    b = dipole_b_t(y[:, :3], model.polar_field_t, model.star_radius_m)
    bm = np.linalg.norm(b, axis=1)
    g = np.sqrt(1+y[:, 3]**2+y[:, 4]**2)
    pitch = np.arctan2(y[:, 4], y[:, 3])
    initial_energy = model.initial_gamma*m_e*c*c
    final = sol.y[:, -1]
    final_g = np.sqrt(1+final[3]**2+final[4]**2)
    energy_difference = (model.initial_gamma-final_g)*m_e*c*c
    ledger = float(final[5])
    radius = np.linalg.norm(y[:, :3], axis=1)
    # True gyro radius perpendicular to B (curvature radius is different).
    gyro_radius = m_e*c*y[:, 4]/(e*bm)
    grad = np.stack([grad_log_dipole_b(x) for x in y[:, :3]])
    b_length = 1/np.linalg.norm(grad, axis=1)
    gyroperiod = 2*np.pi*g*m_e/(e*bm)
    cooling_time = g/(RR_COEFFICIENT*bm*bm*y[:, 4]**2)
    direction = pole*b/bm[:, None]
    # Independently estimate omitted static-field-line curvature radiation.
    # This is a diagnostic bound for these trajectories, not added emission.
    h_geometry = 1.
    bplus = dipole_b_t(y[:, :3]+h_geometry*direction, model.polar_field_t, model.star_radius_m)
    bminus = dipole_b_t(y[:, :3]-h_geometry*direction, model.polar_field_t, model.star_radius_m)
    d_direction_ds = pole*(bplus/np.linalg.norm(bplus, axis=1)[:, None]-bminus/np.linalg.norm(bminus, axis=1)[:, None])/(2*h_geometry)
    inv_radius_squared = np.sum(d_direction_ds*d_direction_ds, axis=1)
    curvature_power = e*e*c/(6*np.pi*epsilon_0)*(g*g-1)**2*inv_radius_squared
    synch_power = synchrotron_power(g, bm, pitch)
    # Remaining perpendicular kinetic energy at fixed beta_parallel for a
    # locally uniform static field is gamma_parallel*(sqrt(1+p_perp²)-1).
    parallel_g = final_g/np.sqrt(1+final[4]**2)
    residual_perp = parallel_g*final[4]**2/(np.sqrt(1+final[4]**2)+1)*m_e*c*c
    sampled_bolometric_energy = float(np.sum(synchrotron_power(g, bm, pitch)*dt))
    return dict(position_m=y[:, :3], magnetic_field_t=b, gamma=g,
                pitch_rad=pitch, residence_time_weights_s=dt,
                diagnostics=dict(duration_s=float(duration), integration_steps=len(sol.t),
                    final_gamma=float(final_g), final_transverse_momentum=float(final[4]),
                    final_radius_rlc=float(np.linalg.norm(final[:3])/model.rlc_m),
                    emitted_energy_j=ledger, initial_energy_j=float(initial_energy),
                    kinetic_energy_difference_j=float(energy_difference),
                    relative_energy_ledger_error=float(abs(ledger-energy_difference)/energy_difference),
                    sampled_bolometric_radiated_energy_j=sampled_bolometric_energy,
                    relative_bolometric_residence_quadrature_error=float(abs(sampled_bolometric_energy-ledger)/ledger),
                    max_gyro_radius_over_b_gradient_length=float(np.max(gyro_radius/b_length)),
                    max_gyroperiod_over_cooling_time=float(np.max(gyroperiod/cooling_time)),
                    max_gyroperiod_over_field_variation_time=float(np.max(gyroperiod*c/b_length)),
                    max_static_dipole_curvature_power_over_synchrotron=float(np.max(curvature_power/synch_power)),
                    max_static_dipole_curvature_critical_frequency_hz=float(np.max(3*c*g**3*np.sqrt(inv_radius_squared)/(4*np.pi))),
                    remaining_perpendicular_energy_bound_j=float(residual_perp),
                    residual_perpendicular_energy_over_radiated=float(residual_perp/ledger),
                    min_gamma_sin_pitch=float(np.min(g*np.sin(pitch))),
                    max_quantum_chi=float(np.max(g*bm*np.sin(pitch)/B_Q_T)),
                    largest_displacement_m=float(np.linalg.norm(final[:3]-y0[:3])),
                    gravity_omitted=True, local_column_rotational_motion_omitted=True))


def boost_to_magnetic_frame(observer_direction, bulk_beta):
    """Aberrate a lab photon into a cell moving with bulk_beta, E'=0 there.

    Returns n', Doppler D=nu/nu' and bulk Lorentz factor. This is a local
    inertial boost; it does not turn the prescribed rotating dipole into a
    solution of Maxwell's equations.
    """
    beta = np.asarray(bulk_beta, float)
    n = np.broadcast_to(np.asarray(observer_direction, float), beta.shape)
    b2 = np.sum(beta*beta, axis=-1)
    if np.any(b2 >= 1) or not np.allclose(np.linalg.norm(n, axis=-1), 1, atol=1e-12):
        raise ValueError('Subluminal boost and unit photon direction required')
    g = 1/np.sqrt(1-b2)
    bn = np.sum(beta*n, axis=-1)
    doppler = 1/(g*(1-bn))
    coeff = np.divide((g-1)*bn, b2, out=np.zeros_like(b2), where=b2 > 0)-g
    nprime = (n+coeff[..., None]*beta)/(g*(1-bn))[..., None]
    return nprime, doppler, g


def _cells(model):
    """Generate physical residence-time weighted cells and absolute supply."""
    incl = np.deg2rad(model.inclination_deg)
    azimuth = (np.arange(model.azimuth_nodes)+.5)*2*np.pi/model.azimuth_nodes
    # Activation is an explicit uniform radial probability, not a cascade.
    rx, rw = np.polynomial.legendre.leggauss(model.radial_nodes)
    radii = (model.activation_radius_max_rlc+model.activation_radius_min_rlc)/2 + rx*(model.activation_radius_max_rlc-model.activation_radius_min_rlc)/2
    radial_probability = rw/2
    ux, uw = np.polynomial.legendre.leggauss(model.annulus_nodes)
    # Integrate magnetic flux in u^2, since dipole cap flux is proportional u^2.
    u2 = .5*(model.annulus_inner**2+model.annulus_outer**2)+.5*(model.annulus_outer**2-model.annulus_inner**2)*ux
    uweights = .5*(model.annulus_outer**2-model.annulus_inner**2)*uw
    data = {k: [] for k in ('position_m', 'magnetic_field_t', 'gamma', 'pitch_rad', 'number_weight', 'footprint_vertices_m')}
    diagnostics = []
    injection_sum = 0.
    energy_rate_sum = 0.
    for pole in (1, -1):
        for u_sq, du_sq in zip(u2, uweights):
            u = np.sqrt(u_sq)
            theta_star = np.arcsin(np.sqrt(model.star_radius_m/model.rlc_m)*u)
            foot = np.array([model.star_radius_m*np.sin(theta_star), 0., pole*model.star_radius_m*np.cos(theta_star)])
            bfoot0 = dipole_b_t(foot, model.polar_field_t, model.star_radius_m)
            # Surface area quadrature for u^2=sin^2(theta)*RLC/Rstar.
            area = model.star_radius_m**3/model.rlc_m/(2*np.cos(theta_star))*du_sq*2*np.pi/model.azimuth_nodes
            for rad, probability in zip(radii, radial_probability):
                tr = column_transport(model, float(rad), float(u), pole)
                diagnostics.append(dict(pole=pole, footprint_ovc=float(u), activation_radius_rlc=float(rad), **tr['diagnostics']))
                for phi in azimuth:
                    bfoot = rotate_y(rotate_z(bfoot0, phi), incl)
                    rfoot = rotate_y(rotate_z(foot, phi), incl)
                    outward = pole*bfoot/np.linalg.norm(bfoot)
                    # n_GJ=|rho_GJ|/e, rho_GJ=-2 epsilon0 Omega.B in SI.
                    ngj = 2*epsilon_0*model.omega*abs(bfoot[2])/e
                    primary_rate = ngj*c*max(0., outward@(rfoot/model.star_radius_m))*area
                    # Each pair supplies two leptons; both have the same spectrum.
                    rate = 2*model.pair_multiplicity_per_primary*primary_rate*probability
                    injection_sum += rate
                    energy_rate_sum += rate*tr['diagnostics']['emitted_energy_j']
                    data['position_m'].append(rotate_y(rotate_z(tr['position_m'], phi), incl))
                    data['magnetic_field_t'].append(rotate_y(rotate_z(pole*tr['magnetic_field_t'], phi), incl))
                    data['gamma'].append(tr['gamma'])
                    data['pitch_rad'].append(tr['pitch_rad'])
                    data['number_weight'].append(rate*tr['residence_time_weights_s'])
                    # Finite physical injection-coordinate footprint for
                    # spatial quadrature. Gauss weights partition the radial
                    # and u² interval into equal-weight support cells. This
                    # geometry is not an angular-emission broadening kernel.
                    ri = np.flatnonzero(radii == rad)[0]
                    ui = np.flatnonzero(u2 == u_sq)[0]
                    rwid = radial_probability[ri]*(model.activation_radius_max_rlc-model.activation_radius_min_rlc)
                    uwid = uweights[ui]
                    # Supports use cumulative quadrature weights, not a claim
                    # that Gauss-node values equal every interior value.
                    rlo = model.activation_radius_min_rlc+sum(radial_probability[:ri])*(model.activation_radius_max_rlc-model.activation_radius_min_rlc)
                    ulo = model.annulus_inner**2+sum(uweights[:ui])
                    corners = []
                    radial_offset = np.linalg.norm(tr['position_m'], axis=1)-rad*model.rlc_m
                    for rr in (rlo, rlo+rwid):
                        radius_here = rr*model.rlc_m+radial_offset
                        for uu in (ulo, ulo+uwid):
                            st = np.sqrt(radius_here*uu/model.rlc_m)
                            ct = np.sqrt(1-st*st)
                            for pp in (phi-np.pi/model.azimuth_nodes, phi+np.pi/model.azimuth_nodes):
                                corner = radius_here[:, None]*np.stack((st*np.cos(pp), st*np.sin(pp), pole*ct), axis=1)
                                corners.append(rotate_y(corner, incl))
                    data['footprint_vertices_m'].append(np.stack(corners, axis=1))
    data = {k: np.concatenate(v, axis=0) for k, v in data.items()}
    return data, diagnostics, injection_sum, energy_rate_sum


def generate_events(model, output_npz, output_json, compact_camera=False):
    """Write observer-specific physical source events for a sparse camera.

    NPZ units: positions_m [m], source_time_s [s], source_interval_s [s],
    number_weights [physical leptons], power_bins_w_sr [W sr^-1 per bin]
    PER ONE PHYSICAL LEPTON, wavelength_edges_m [m]. Multiply power by number weight
    and source interval ONCE; deposit energy at (t-n.r/c) modulo period.
    No arrival Jacobian is already applied. Last cell/time axis is not an
    observed optical particle trail. The array is a stationary population
    quadrature repeated over a prescribed rotation, not a time-dependent PIC.
    """
    cells, columns, injection_rate, loss_luminosity = _cells(model)
    edges = np.linspace(400e-9, 700e-9, model.wavelength_nodes+1)
    lam = .5*(edges[1:]+edges[:-1])
    cie_path = Path(__file__).resolve().parents[3]/'data/color/CIE_xyz_1931_2deg.csv'
    cie_data = np.loadtxt(cie_path, delimiter=',')
    cie_wave, cie = cie_data[:, 0]*1e-9, cie_data[:, 1:]
    keys = ('xyz_lm_sr', 'phase_index', 'cell_index', 'event_band_power_w_sr') if compact_camera else ('positions_m', 'source_time_s', 'source_interval_s', 'number_weights', 'power_bins_w_sr', 'xyz_lm_sr', 'phase_index', 'cell_index', 'pattern_beta')
    events = {k: [] for k in keys}
    n = model.observer_direction
    phase_power = []
    phase_count = []
    r = cells['position_m']
    b = cells['magnetic_field_t']
    bm = np.linalg.norm(b, axis=1)
    max_bulk = 0.
    harmonic_diagnostics = []
    phase_spectral_totals = []
    xyz_quantization_error_bound = np.zeros(3)
    exact_xyz_energy = np.zeros(3)
    for index in range(model.phase_nodes):
        # Midpoint source-time quadrature; one complete physical rotation.
        phase = (index+.5)*2*np.pi/model.phase_nodes
        time = phase/model.omega
        xyz = rotate_z(r, phase)
        local_b = rotate_z(b, phase)
        beta = model.omega*np.stack((-xyz[:, 1], xyz[:, 0], np.zeros(len(xyz))), axis=1)/c
        max_bulk = max(max_bulk, float(np.max(np.linalg.norm(beta, axis=1))))
        if model.aberration:
            nprime, doppler, bulk_gamma = boost_to_magnetic_frame(n, beta)
        else:
            nprime = np.broadcast_to(n, beta.shape)
            doppler, bulk_gamma = np.ones(len(xyz)), np.ones(len(xyz))
        view = np.arccos(np.clip(np.sum(nprime*local_b, axis=1)/bm, -1, 1))
        bands, xyz_lm, diag = binned_harmonic_spectrum(cells['gamma'], bm, cells['pitch_rad'], view, doppler, bulk_gamma, edges, cie_wave, cie)
        omitted_bound = diag.pop('omitted_power_bounds_w_sr')
        diag['weighted_omitted_power_bound_w_sr'] = float(omitted_bound@cells['number_weight'])
        harmonic_diagnostics.append(diag)
        keep = np.sum(bands, axis=1) > 0
        weights = cells['number_weight'][keep]
        total = float(np.sum(bands[keep]*weights[:, None]))
        phase_spectral_totals.append(np.sum(bands[keep]*weights[:, None], axis=0))
        phase_power.append(total)
        phase_count.append(int(np.count_nonzero(keep)))
        if compact_camera:
            stored_xyz = xyz_lm[keep].astype(np.float32)
            xyz_quantization_error_bound += np.sum(np.abs(stored_xyz.astype(float)-xyz_lm[keep])*weights[:, None], axis=0)*model.period_s/model.phase_nodes
            exact_xyz_energy += np.sum(xyz_lm[keep]*weights[:, None], axis=0)*model.period_s/model.phase_nodes
            events['xyz_lm_sr'].append(stored_xyz)
            events['event_band_power_w_sr'].append(np.sum(bands[keep], axis=1))
        else:
            events['positions_m'].append(xyz[keep])
            events['source_time_s'].append(np.full(len(weights), time))
            events['source_interval_s'].append(np.full(len(weights), model.period_s/model.phase_nodes))
            events['number_weights'].append(weights)
            events['power_bins_w_sr'].append(bands[keep])
            events['xyz_lm_sr'].append(xyz_lm[keep])
            events['pattern_beta'].append(beta[keep])
        events['phase_index'].append(np.full(len(weights), index, dtype=np.int16))
        events['cell_index'].append(np.flatnonzero(keep).astype(np.int32))
    arrays = {k: np.concatenate(v, axis=0) for k, v in events.items()}
    arrays.update(wavelength_m=lam, wavelength_edges_m=edges, observer_direction=n, period_s=np.array(model.period_s),
                  phase_angles_rad=(np.arange(model.phase_nodes)+.5)*2*np.pi/model.phase_nodes,
                  reference_cell_positions_m=r, cell_footprint_vertices_reference_m=cells['footprint_vertices_m'],
                  reference_cell_number_weights=cells['number_weight'],
                  phase_spectral_power_w_sr=np.asarray(phase_spectral_totals),
                  source_interval_scalar_s=np.array(model.period_s/model.phase_nodes),
                  compact_camera=np.array(compact_camera))
    output_npz = Path(output_npz)
    output_npz.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output_npz, **arrays)
    q = np.asarray(phase_power)
    if compact_camera:
        ci, pi = arrays['cell_index'], arrays['phase_index']
        angle = arrays['phase_angles_rad'][pi]
        pos = rotate_z(r[ci], angle)
        event_band_energy = arrays['event_band_power_w_sr']*cells['number_weight'][ci]*(model.period_s/model.phase_nodes)
        source_time = angle/model.omega
    else:
        pos = arrays['positions_m']
        event_band_energy = np.sum(arrays['power_bins_w_sr'], axis=1)*arrays['number_weights']*arrays['source_interval_s']
        source_time = arrays['source_time_s']
    arrival = np.remainder(source_time-pos@n/c, model.period_s)
    received_energy, _ = np.histogram(arrival, bins=np.linspace(0, model.period_s, model.phase_nodes+1), weights=event_band_energy)
    received = received_energy/(model.period_s/model.phase_nodes)
    report = dict(schema='sideratlas.physical-synchrotron-events.compact.v1' if compact_camera else 'sideratlas.physical-synchrotron-events.v1', model=asdict(model),
        scientific_status='Vela-informed prescribed physical synchrotron model; not fitted or validated Vela, and not PIC-derived emission',
        observing_band='400-700 nm, vacuum, unextincted, directional photon power',
        observer_direction_toward_camera=n.tolist(),
        units={k: v for k, v in [('positions_m','m'), ('source_time_s','s'), ('source_interval_s','s'), ('number_weights','physical leptons'), ('power_bins_w_sr','W sr^-1 per bin per physical lepton'), ('xyz_lm_sr','XYZ luminous power per sr per physical lepton; 683 applied ONCE'), ('wavelength_edges_m','m'), ('pattern_beta','v_pattern/c; rotation of prescribed cell, not lepton velocity')]},
        time_frame='lab source time; rotating steady population, midpoint phase quadrature; not arrival time',
        spectrum_frame='lab photon frequency and solid angle after local magnetic-cell boost',
        number_weight_convention='2 * pair_multiplicity * chosen Goldreich-Julian primary flux * local residence-time quadrature weight; linear incoherent sum',
        periodicity='exactly periodic by construction in prescribed rotating population; source interval is P/Nphase, not cooling residence interval',
        physical_rlc_m=model.rlc_m, cell_quadrature_count=len(r), source_event_count=len(pos),
        compact_camera=compact_camera,
        compact_reconstruction='positions=Rz(phase_angle)*reference_cell_positions[cell_index]; time=phase_angle/Omega; dt=P/Nphase; number_weights=reference_cell_number_weights[cell_index]; beta_pattern=Omega*(-y,x,0)/c' if compact_camera else None,
        xyz_float32_storage_absolute_energy_error_bound_lm_s_sr=xyz_quantization_error_bound.tolist() if compact_camera else None,
        xyz_float32_storage_relative_energy_error_bound=np.divide(xyz_quantization_error_bound, exact_xyz_energy, out=np.zeros(3), where=exact_xyz_energy > 0).tolist() if compact_camera else None,
        retained_events_per_phase=phase_count,
        zeros_discarded='Harmonic rows are discarded only when a Cauchy/Bessel upper bound on TOTAL per-lepton band power is below1e-45 W/sr, or floating-point underflow; the weighted omitted-power bounds are recorded per phase. No angular cone or display-brightness cutoff is introduced.',
        injected_total_leptons_per_local_second=float(injection_rate),
        local_synchrotron_loss_luminosity_w=float(loss_luminosity),
        total_source_frame_band_energy_j_sr=float(event_band_energy.sum()),
        retarded_band_energy_j_sr=float(received_energy.sum()),
        relative_retarded_energy_histogram_error=float(abs(received_energy.sum()-event_band_energy.sum())/event_band_energy.sum()) if event_band_energy.sum() > 0 else 0.,
        positive_optical_band_signal=bool(event_band_energy.sum() > 0),
        source_band_power_per_phase_w_sr=q.tolist(), retarded_band_power_per_phase_w_sr=received.tolist(),
        phases_above_one_percent_peak=int(np.count_nonzero(received > .01*received.max())),
        phases_above_one_ten_thousandth_peak=int(np.count_nonzero(received > .0001*received.max())),
        max_corotation_beta=max_bulk, column_diagnostics=columns,
        harmonic_diagnostics_per_phase=harmonic_diagnostics,
        total_omitted_band_energy_upper_bound_j_sr=float(sum(d['weighted_omitted_power_bound_w_sr'] for d in harmonic_diagnostics)*model.period_s/model.phase_nodes),
        spectral_representation='Exact delta-line power integrated into wavelength bins; XYZ computed at actual individual line wavelengths using piecewise-linear official CIE1931 data, not from the binned density',
        cie_dataset_sha256=hashlib.sha256(cie_path.read_bytes()).hexdigest(),
        footprint_definition='Eight physical injection-coordinate corners, carried along the same dipole geometry and cooling-node radial displacement; finite Gauss supports partitioned by weights. Footprint integration is a spatial quadrature approximation, NOT angular beam broadening.',
        assumptions=[
            'Chosen physical B_p=3e12 G, R=13 km and P=0.08937108859 s; values are model inputs',
            'Prescribed rotating near-zone dipole; local E=0 magnetic cells; no global Maxwell/PIC feedback',
            'Static dipole cap mapping and annulus u=0.8-0.9; both poles and all azimuths',
            'Existing pairs receive chosen gamma and pitch at an assumed uniformly distributed activation radius; no radio heating or cascade is calculated',
            'Pair multiplicity is adopted from Harding2018; their full pair spectrum is NOT reproduced by a chosen monoenergetic input',
            'Guiding-centre focusing plus gyro-averaged classical synchrotron drag; no curvature loss, GR, inverse Compton, optical absorption, extinction or atmosphere',
            'Particle counts use local static-column residence times; moving-cell population is prescribed, not a relativistic continuity solution',
            'Exact Schott harmonics, including low harmonics; stop at p_perp=.001 m_e c with residual perpendicular-energy bound recorded',
            'Uniform random gyro phases, optically thin incoherent emission; positrons/electrons have identical intensity colors',
            'Sampling is explicitly finite; no Vela pulse shape, optical slope or absolute normalization is fit or validated'],
        references=[
            'Harding et al. 2018, doi:10.3847/2041-8213/aaf3b2, https://arxiv.org/abs/1811.11157',
            'Harding & Kalapotharakos 2015, doi:10.1088/0004-637X/811/1/63, https://arxiv.org/abs/1508.06251',
            'Lai et al. 2025, doi:10.1093/mnras/staf1295, https://arxiv.org/abs/2509.04572',
            'Marcowith & Malzac2003, doi:10.1051/0004-6361:20031060, https://arxiv.org/abs/astro-ph/0309107',
            'Mignani et al. 2007, doi:10.1051/0004-6361:20077774, https://arxiv.org/abs/0707.2036'],
        npz_sha256=hashlib.sha256(output_npz.read_bytes()).hexdigest(),
        physical_pic_energy_rescaling=False, pair_cascade_calculated=False,
        radio_absorption_calculated=False, self_consistent_em_feedback=False,
        atmosphere_calculated=False, vela_optical_validation=False, sampling_convergence_established=False)
    Path(output_json).write_text(json.dumps(report, indent=2)+'\n')
    return report
