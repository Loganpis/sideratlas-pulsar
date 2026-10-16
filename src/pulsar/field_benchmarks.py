"""Independent vacuum/source and absorbing-layer verification experiments.

These experiments do not alter or validate the production plasma closure.
The evolved dipole test uses a point source; a separate analytic Deutsch
reference handles the finite conducting sphere. The matched absorber below
is a one-dimensional normal-incidence prototype, not a three-dimensional
PML. Units are c=Omega=1, with static A=m cross r/r**3.
"""
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.special import spherical_jn, spherical_yn

from .plasma import Magnetosphere, PlasmaConfig, cross, curl, divergence, energy


def dipole_moments(t, inclination_deg=65.):
    """Uniformly rotating unit moment and its first two time derivatives."""
    t = np.asarray(t)
    s, c = np.sin(np.deg2rad(inclination_deg)), np.cos(np.deg2rad(inclination_deg))
    m = np.array([s*np.cos(t), s*np.sin(t), np.full_like(t, c)])
    dm = np.array([-s*np.sin(t), s*np.cos(t), np.zeros_like(t)])
    ddm = np.array([-s*np.cos(t), -s*np.sin(t), np.zeros_like(t)])
    return m, dm, ddm


def retarded_dipole(xyz, time, inclination_deg=65., min_radius=None):
    """Return independent analytic A, E, B at component-first positions.

    All moments are evaluated at t-r. A is in the zero-scalar-potential
    gauge, so E=-dA/dt. A nonzero min_radius is only a bounded continuation
    for excluded numerical source cells; it must never be used in validation
    shells. The exact solution is singular at the point source.
    """
    xyz = np.asarray(xyz, dtype=float)
    r = np.sqrt(np.sum(xyz*xyz, axis=0))
    if min_radius is None:
        if np.any(r == 0):
            raise ValueError('Point-dipole field is undefined at the origin')
    else:
        if min_radius <= 0:
            raise ValueError('min_radius must be positive')
        r = np.maximum(r, min_radius)
    n = xyz/r
    m, dm, ddm = dipole_moments(time-r, inclination_deg)
    a = cross(m, n)/r**2 + cross(dm, n)/r
    e = -cross(dm, n)/r**2 - cross(ddm, n)/r
    b = ((3*n*np.sum(n*m, axis=0)-m)/r**3
         + (3*n*np.sum(n*dm, axis=0)-dm)/r**2
         + (n*np.sum(n*ddm, axis=0)-ddm)/r)
    return a, e, b


def point_dipole_luminosity(inclination_deg=65.):
    """Flux for this code's unit-moment normalization and rationalized EM."""
    return float(8*np.pi/3*np.sin(np.deg2rad(inclination_deg))**2)


def deutsch_luminosity(star_radius=.35, inclination_deg=65.):
    """Finite conducting sphere: Petri 2015 Eq. 46a, unit static moment.

    a=R/R_LC; the second term is the induced electric-quadrupole radiation.
    A net electric monopole does not change the radial radiative luminosity.
    """
    a = float(star_radius)
    if not 0 < a < 1:
        raise ValueError('The stellar surface must be inside the light cylinder')
    correction = 1/(1+a*a)+3*a**4/(5*(a**6-3*a**4+36))
    return point_dipole_luminosity(inclination_deg)*correction


def deutsch_fields(xyz, time, star_radius=.35, inclination_deg=65.):
    """Exact exterior E/B of an uncharged, uniformly rotating conducting sphere.

    The magnetic surface normal is 2 m(t).n/R**3 and tangential E satisfies
    perfect corotation. The outgoing magnetic l=1 and electric l=2 spherical
    Hankel modes are the Deutsch solution, equivalent to Petri 2015 Eq. 24.
    Phase convention exp(i(phi-t)), c=Omega=1, B_equatorial*R**3=1.
    This is a reference only: no baseline boundary or snapshot is modified.
    """
    R = float(star_radius)
    xyz = np.asarray(xyz, dtype=float)
    r = np.sqrt(np.sum(xyz*xyz, axis=0))
    if not 0 < R < 1 or np.any(r < R*(1-1e-13)):
        raise ValueError('Deutsch reference is defined on/outside a sub-light-cylinder star')
    ct = xyz[2]/r
    st = np.sqrt(np.maximum(0., 1-ct*ct))
    phi = np.arctan2(xyz[1], xyz[0])
    cp, sp = np.cos(phi), np.sin(phi)
    phase = np.exp(1j*(phi-time))
    sc, cc = np.sin(np.deg2rad(inclination_deg)), np.cos(np.deg2rad(inclination_deg))
    def hankel(l, x):
        return spherical_jn(l, x)+1j*spherical_yn(l, x)
    def radial_derivative(l, x):
        return hankel(l, x)+x*(spherical_jn(l, x, derivative=True)
                              +1j*spherical_yn(l, x, derivative=True))
    h1, h2 = hankel(1, r), hankel(2, r)
    h1R, d2R = hankel(1, R), radial_derivative(2, R)
    d1, d2 = radial_derivative(1, r), radial_derivative(2, r)
    br = 2*(R**3/r**3*cc*ct+R/r*h1/h1R*sc*st*phase)/R**3
    bt = (R**3/r**3*cc*st+(R/r*d1/h1R+R**2*h2/d2R)*sc*ct*phase)/R**3
    bp = (R/r*d1/h1R+R**2*h2/d2R*(2*ct*ct-1))*1j*sc*phase/R**3
    er = (-R**4/r**4*(3*ct*ct-1)*cc+6*sc*st*ct*phase*R/r*h2/d2R)/R**2
    et = (sc*phase*(R/r*d2/d2R*(2*ct*ct-1)-h1/h1R)-R**4/r**4*2*st*ct*cc)/R**2
    ep = (R/r*d2/d2R-h1/h1R)*1j*sc*ct*phase/R**2
    erhat = np.array([st*cp, st*sp, ct])
    ethat = np.array([ct*cp, ct*sp, -st])
    ephat = np.array([-sp, cp, np.zeros_like(r)])
    return np.real(erhat*er+ethat*et+ephat*ep), np.real(erhat*br+ethat*bt+ephat*bp)


def deutsch_spherical_flux(radius, time=0., star_radius=.35, inclination_deg=65.):
    dirs, weights = sphere_quadrature()
    e, b = deutsch_fields(radius*dirs, time, star_radius, inclination_deg)
    return float(np.sum(np.sum(cross(e, b)*dirs, axis=0)*weights)*radius**2)


def finite_conductor_reference(star_radius=.35, inclination_deg=65.):
    """Machine-checkable finite-source reference, surface data and flux."""
    dirs, _ = sphere_quadrature()
    t = .731
    e, b = deutsch_fields(star_radius*dirs, t, star_radius, inclination_deg)
    m, _, _ = dipole_moments(t, inclination_deg)
    br_expected = 2*np.einsum('i,ijk->jk', m, dirs)/star_radius**3
    beta = star_radius*np.array([-dirs[1], dirs[0], np.zeros_like(dirs[0])])
    conductor_residual = e+cross(beta, b)
    tangent = conductor_residual-dirs*np.sum(conductor_residual*dirs, axis=0)
    return dict(star_radius=star_radius, inclination_deg=inclination_deg,
                point_luminosity=point_dipole_luminosity(inclination_deg),
                finite_luminosity=deutsch_luminosity(star_radius, inclination_deg),
                radial_B_surface_relative_error=float(np.linalg.norm(np.sum(b*dirs, axis=0)-br_expected)/np.linalg.norm(br_expected)),
                tangential_E_surface_relative_error=float(np.linalg.norm(tangent)/np.linalg.norm(e)),
                fluxes=[dict(radius=r, phase=t, flux=deutsch_spherical_flux(r, t, star_radius, inclination_deg)) for r in (star_radius, 1., 1.5, 2.5, 10.)],
                reference='https://arxiv.org/pdf/1503.05307; Eqs. 24, 46a',
                caveat='Physical finite sphere reference; baseline mask is not a sharply resolved conductor surface')


def sphere_quadrature(n_mu=32, n_phi=96):
    mu, weights = np.polynomial.legendre.leggauss(n_mu)
    phi = np.arange(n_phi)*2*np.pi/n_phi
    dirs = np.array(np.broadcast_arrays(np.sqrt(1-mu[:, None]**2)*np.cos(phi),
                    np.sqrt(1-mu[:, None]**2)*np.sin(phi), mu[:, None]))
    return dirs, weights[:, None]*2*np.pi/n_phi


def analytic_spherical_flux(radius, time=0., inclination_deg=65.):
    dirs, weights = sphere_quadrature()
    _, e, b = retarded_dipole(radius*dirs, time, inclination_deg)
    return float(np.sum(np.sum(cross(e, b)*dirs, axis=0)*weights)*radius**2)


def grid_spherical_flux(axis, e, b, radii=(1.25, 1.75, 2.25)):
    dirs, weights = sphere_quadrature()
    # Interpolate fields, then form flux; use the identical operation for
    # numerical fields and sampled analytic-reference fields.
    values = []
    ie = RegularGridInterpolator((axis,)*3, np.moveaxis(e, 0, -1))
    ib = RegularGridInterpolator((axis,)*3, np.moveaxis(b, 0, -1))
    for radius in radii:
        pts = np.moveaxis(radius*dirs, 0, -1)
        es, bs = np.moveaxis(ie(pts), -1, 0), np.moveaxis(ib(pts), -1, 0)
        values.append(float(np.sum(np.sum(cross(es, bs)*dirs, axis=0)*weights)*radius**2))
    return values


def analytic_maxwell_residual(n, extent=4., time=.731):
    """Continuum reference sampled on the baseline central-difference mesh.

    The exact continuum fields are used, rather than curl(A) for both sides
    of an identity. Temporal derivatives use independent symmetric samples.
    Fixed physical shells exclude the source and periodic grid seams.
    """
    h = 2*extent/n
    axis = (np.arange(n)+.5)*h-extent
    xyz = np.array(np.meshgrid(axis, axis, axis, indexing='ij'))
    r = np.sqrt(np.sum(xyz*xyz, axis=0))
    mask = (r > 1.0) & (r < 2.4)
    a, e, b = retarded_dipole(xyz, time)
    eps = 1e-5
    _, em, bm = retarded_dipole(xyz, time-eps)
    _, ep, bp = retarded_dipole(xyz, time+eps)
    edt, bdt = (ep-em)/(2*eps), (bp-bm)/(2*eps)
    norm = lambda x: float(np.linalg.norm(x[..., mask]))
    return dict(n=n, h=h, shell=[1., 2.4],
                faraday_relative_L2=norm(bdt+curl(e, h))/norm(bdt),
                ampere_relative_L2=norm(edt-curl(b, h))/norm(edt),
                curlA_relative_L2=norm(curl(a, h)-b)/norm(b),
                divB_relative_L2=norm(divergence(b, h))/norm(b),
                divE_relative_L2=norm(divergence(e, h))/norm(e))


class ExactDipoleBoundary(Magnetosphere):
    """Vacuum benchmark with a compatible A/E source at every RK substage.

    Both the source ball and outer guard cells are exact time-dependent
    Dirichlet data. Only the exterior is evolved. The outer boundary supplies
    an exact known solution: this class DOES NOT test absorbing boundaries.
    Baseline Maxwell operators/RK4 are retained; physical closures, filters,
    masks, spinup and the A sponge are excluded from this controlled test.
    """
    def __init__(self, n=48, source_radius=.65, extent=4.):
        if source_radius < .6:
            raise ValueError('Benchmark source must enclose the bounded continuation')
        config = PlasmaConfig(n=n, half_extent=extent, plasma=False,
                              conductor=False, filter_strength=0., sponge_strength=0.)
        super().__init__(config)
        self.source_radius = source_radius
        self.guard = ((self.r <= source_radius)
                      | (np.max(np.abs(self.xyz), axis=0) > extent-3*self.h))
        self.a, self.e, _ = self.exact(0.)
        self.initial_energy = energy(self.e, self.b, self.h)

    def exact(self, time):
        return retarded_dipole(self.xyz, time, self.config.inclination_deg, min_radius=.2)

    def boundary_rhs(self, a, e, time):
        a, e = a.copy(), e.copy()
        exact_a, exact_e, _ = retarded_dipole(self.xyz[:, self.guard], time,
                                             self.config.inclination_deg, min_radius=.2)
        a[:, self.guard], e[:, self.guard] = exact_a, exact_e
        return self.rhs(a, e)

    def step(self, dt=None):
        dt = self.config.courant*self.h if dt is None else dt
        if not 0 < dt <= .35*self.h:
            raise ValueError('Timestep violates stability bound')
        a, e, t = self.a, self.e, self.time
        ka1, ke1 = self.boundary_rhs(a, e, t)
        ka2, ke2 = self.boundary_rhs(a+.5*dt*ka1, e+.5*dt*ke1, t+.5*dt)
        ka3, ke3 = self.boundary_rhs(a+.5*dt*ka2, e+.5*dt*ke2, t+.5*dt)
        ka4, ke4 = self.boundary_rhs(a+dt*ka3, e+dt*ke3, t+dt)
        self.a = a+dt/6*(ka1+2*ka2+2*ka3+ka4)
        self.e = e+dt/6*(ke1+2*ke2+2*ke3+ke4)
        self.time += dt
        self.step_number += 1
        exact_a, exact_e, _ = retarded_dipole(self.xyz[:, self.guard], self.time,
                                             self.config.inclination_deg, min_radius=.2)
        self.a[:, self.guard] = exact_a
        self.e[:, self.guard] = exact_e

    def benchmark_diagnostics(self):
        _, e_ref, b_ref = self.exact(self.time)
        b = self.b
        mask = (self.r > 1.0) & (self.r < 2.4)
        norm = lambda f: np.linalg.norm(f[:, mask])
        flux = np.asarray(grid_spherical_flux(self.axis, self.e, b))
        analytic_sampled = np.asarray(grid_spherical_flux(self.axis, e_ref, b_ref))
        target = point_dipole_luminosity(self.config.inclination_deg)
        return dict(n=self.config.n, h=self.h, time=self.time,
                    source_radius=self.source_radius, source_radius_in_cells=self.source_radius/self.h,
                    steps=self.step_number,
                    e_relative_L2=float(norm(self.e-e_ref)/norm(e_ref)),
                    b_relative_L2=float(norm(b-b_ref)/norm(b_ref)),
                    divB_relative=float(np.linalg.norm(divergence(b, self.h))*self.h/np.linalg.norm(b)),
                    luminosity_exact=target, flux_radii=[1.25, 1.75, 2.25],
                    flux=flux.tolist(), sampled_analytic_flux=analytic_sampled.tolist(),
                    flux_relative_error=(np.abs(flux/target-1)).tolist(),
                    interpolation_relative_error=(np.abs(analytic_sampled/target-1)).tolist())


def conductor_gauge_seam(n=64, time=2., dt=None):
    """Measure magnetic forcing caused by masking a source gauge correction.

    For a rigid dipole, E_corotation + dA_static/dt is a gradient in the
    continuum. Masking that gradient introduces a nonzero curl at the source
    seam. This is a diagnostic of the actual baseline source, not a claim
    that ideal-conductor E and a rotating B are physically incompatible.
    """
    model = Magnetosphere(PlasmaConfig(n=n, plasma=False))
    dt = model.h*.3 if dt is None else dt
    a = model.star_potential(time)
    b = curl(a, model.h)
    _, omega = model.rotation(time)
    e = -cross(model.beta*omega, b)
    eps = 1e-5
    da_dt = (model.star_potential(time+eps)-model.star_potential(time-eps))/(2*eps)
    gradient = da_dt+e
    masked_curl = curl(model.star_mask[None]*gradient, model.h)
    seam = (model.r > model.config.star_radius-2*model.h) & (model.r < model.config.star_radius+3*model.h)
    return dict(n=n, star_radius_in_cells=model.config.star_radius/model.h,
                seam_relative_B_change_per_step=float(dt*np.linalg.norm(masked_curl[:, seam])/np.linalg.norm(b[:, seam])),
                mask_curl_norm=float(np.linalg.norm(masked_curl[:, seam])),
                interpretation='First-order magnetic correction from masked gauge-reset term; not total step error')


def _d1(a, h):
    return (np.roll(a, -1)-np.roll(a, 1))/(2*h)


def _filter1(a, strength):
    return a-strength/16*(np.roll(a, 2)-4*np.roll(a, 1)+6*a-4*np.roll(a, -1)+np.roll(a, -2))


def wave_packet_absorber(n=128, kind='potential', strength=4., filter_strength=0.,
                         end_time=10., record_stride=4):
    """Normal-incidence wave, baseline grid/domain/sponge and RK4 splitting.

    kind='potential' reproduces damping A_y and E_y, B_z=D_x A_y.
    kind='matched_fields' damps physical E_y and B_z equally. It is matched
    for THIS transverse 1-D problem; multidimensional PML remains separate.
    E=B initially is an exactly rightgoing *discrete* characteristic. The
    independent leftgoing energy (E-B)^2/4 measures generated reflection.
    Returned wrap energy includes leakage through both periodic edge layers.
    """
    if kind not in ('potential', 'matched_fields'):
        raise ValueError('Unknown absorber')
    if n < 32 or n % 2 or strength < 0 or not 0 <= filter_strength <= .3 or end_time <= 0:
        raise ValueError('Invalid packet parameters')
    h, extent, start = 8/n, 4., 2.8
    x = (np.arange(n)+.5)*h-extent
    a = np.exp(-.5*((x+1.)/.45)**2)
    b = _d1(a, h)
    e = b.copy()
    initial = float(np.sum(.5*(e*e+b*b))*h)
    sigma = strength*np.clip((np.abs(x)-start)/(extent-start), 0, 1)**3
    q = a.copy() if kind == 'potential' else b.copy()
    t, step = 0., 0
    history = []
    interior = np.abs(x) < 2.5
    left_interior = interior & (x < 0)
    max_reflected, max_wrapped = 0., 0.
    def rhs(q, e):
        if kind == 'potential':
            return -e, -_d1(_d1(q, h), h)
        return -_d1(e, h), -_d1(q, h)
    while t < end_time:
        dt = min(.3*h, end_time-t)
        k1q, k1e = rhs(q, e)
        k2q, k2e = rhs(q+dt*k1q/2, e+dt*k1e/2)
        k3q, k3e = rhs(q+dt*k2q/2, e+dt*k2e/2)
        k4q, k4e = rhs(q+dt*k3q, e+dt*k3e)
        q += dt/6*(k1q+2*k2q+2*k3q+k4q)
        e += dt/6*(k1e+2*k2e+2*k3e+k4e)
        q, e = _filter1(q, filter_strength), _filter1(e, filter_strength)
        damp = np.exp(-dt*sigma)
        q, e = q*damp, e*damp
        b = _d1(q, h) if kind == 'potential' else q
        t += dt
        step += 1
        right, left = (e+b)/2, (e-b)/2
        reflected = float(np.sum(left[interior]**2)*h/initial)
        wrapped = float(np.sum(right[left_interior]**2)*h/initial)
        if t > 3.:
            max_reflected = max(max_reflected, reflected)
        if t > 6.:
            max_wrapped = max(max_wrapped, wrapped)
        if step % record_stride == 0 or t >= end_time:
            history.append(dict(time=t, reflected_energy_fraction=reflected,
                                left_interior_rightgoing_energy_fraction=wrapped,
                                total_energy_fraction=float(np.sum(e*e+b*b)*.5*h/initial)))
    return dict(n=n, h=h, kind=kind, strength=strength, filter_strength=filter_strength,
                duration=end_time, initial_energy=initial,
                peak_reflected_energy_fraction=max_reflected,
                peak_periodic_wrap_energy_fraction=max_wrapped,
                continuum_two_layer_transmission_energy=float(np.exp(-strength*1.2)),
                history=history)


def gauge_sponge_check(n=64, gauge_offset=1., dt=.1):
    """A spatially constant A is pure gauge; the baseline sponge creates B."""
    h = 8/n
    x = (np.arange(n)+.5)*h-4
    sigma = 4*np.clip((np.abs(x)-2.8)/1.2, 0, 1)**3
    a = np.full(n, gauge_offset)
    after_b = _d1(a*np.exp(-dt*sigma), h)
    return dict(n=n, gauge_offset=gauge_offset, dt=dt,
                initial_field_energy=0.,
                added_magnetic_energy_per_transverse_area=float(.5*np.sum(after_b*after_b)*h),
                max_created_B=float(np.max(np.abs(after_b))),
                matched_fields_added_energy=0.)
