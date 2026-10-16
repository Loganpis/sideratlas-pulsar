"""Prescribed-velocity transport through supplied laboratory E and B fields.

Implements Kalapotharakos, Harding & Kazanas (2014), Eq.12, not corotation
continued outside the light cylinder and not a kinetic/PIC orbit solver.
https://arxiv.org/html/1310.3545v2

Positions are x/R_LC and times are Omega*t; fields(position,time) returns
Cartesian (E,B), in statvolt/cm and gauss. Batch time must broadcast to the
position batch. No field scaling, charge selection, or absolute-value work
is hidden. A fixed signed branch specifies motion along/against local B;
charge_sign independently specifies positron/electron charge. Both are
prescribed populations, not pair creation or current closure.
"""
from dataclasses import asdict, dataclass

import numpy as np

from .radiation import (
    CurvatureTransport, E_STATC, MEC2_ERG, _advance_constant,
    critical_energy_gev, curvature_photon_yield, curvature_power_erg_s,
)


def _vectors(value, name):
    value = np.asarray(value, dtype=float)
    if value.shape[-1:] != (3,) or not np.all(np.isfinite(value)):
        raise ValueError(f"{name} must contain finite Cartesian three-vectors")
    return value


def _signs(value, shape, name):
    value = np.broadcast_to(np.asarray(value, dtype=float), shape)
    if not np.all((value == 1) | (value == -1)):
        raise ValueError(f"{name} must be +1 or -1")
    return value


def prescribed_velocity(electric, magnetic, branch=1):
    """K14 Eq.12 unit laboratory velocity, including electric dominance.

    beta = E cross B/(B^2+E0^2) + branch*sqrt(1-beta_drift^2)*Bhat.
    E0 is the nonnegative invariant electric field defined in K14 Eq.6.
    The branch is independent of charge. A magnetic null is undefined under
    this prescription and raises; it is never repaired by inventing Bhat.
    Near-null/strongly electric regions require physical validity checks.
    """
    e, b = np.broadcast_arrays(_vectors(electric, "E"), _vectors(magnetic, "B"))
    sign = _signs(branch, e.shape[:-1], "branch")
    # Scale locally before invariants, avoiding overflow and cancellation.
    scale = np.maximum(np.max(np.abs(e), axis=-1), np.max(np.abs(b), axis=-1))
    if np.any(scale == 0):
        raise ValueError("Velocity is undefined where both fields vanish")
    e = e/scale[..., None]
    b = b/scale[..., None]
    b2 = np.sum(b*b, axis=-1)
    if np.any(b2 <= np.finfo(float).tiny):
        raise ValueError("K14 velocity is undefined at a magnetic null")
    e2 = np.sum(e*e, axis=-1)
    eb = np.sum(e*b, axis=-1)
    delta = b2-e2
    discriminant = np.hypot(delta, 2*eb)
    e02 = np.where(delta >= 0,
                   np.divide(2*eb**2, discriminant+delta,
                             out=np.zeros_like(eb), where=discriminant+delta > 0),
                   (discriminant-delta)/2)
    drift = np.cross(e, b)/(b2+e02)[..., None]
    drift2 = np.sum(drift*drift, axis=-1)
    if np.any(drift2 > 1+64*np.finfo(float).eps):
        raise ValueError("Field invariants produced a superluminal drift")
    parallel = sign*np.sqrt(np.maximum(0., 1-drift2))
    return drift + parallel[..., None]*b/np.sqrt(b2)[..., None]


def select_branch(fields, position, time=0., flow=1):
    """Select the most outward (+1) or inward (-1) branch at injection.

    Keep the returned branch attached to that seed during transport. There
    need not be an inward branch outside the LC: a radial wind can drag both
    branches outward. An impossible radial choice raises, rather than
    relabeling an outward stream as inflow. A tie chooses the +B branch.
    """
    p = _vectors(position, "Position")
    radius = np.linalg.norm(p, axis=-1)
    if np.any(radius == 0):
        raise ValueError("Radial flow is undefined at the origin")
    flow = _signs(flow, p.shape[:-1], "flow")
    e, b = fields(p, time)
    plus = np.sum(p*prescribed_velocity(e, b, 1), axis=-1)/radius
    minus = np.sum(p*prescribed_velocity(e, b, -1), axis=-1)/radius
    sign = np.where(flow*plus >= flow*minus, 1., -1.)
    selected = np.where(sign == 1, plus, minus)
    if np.any(flow*selected <= 0):
        raise ValueError("Requested radial flow is unavailable in this field")
    return sign


def _velocity_at(fields, position, time, branch):
    return prescribed_velocity(*fields(position, time), branch=branch)


def _validate_field_stencil(fields, *positions):
    """Honor a provider's explicit discontinuity/validity gate, if supplied."""
    validator = getattr(fields, "validate_transport_stencil", None)
    if validator is not None:
        validator(np.stack(positions))


def _rk4(fields, position, time, branch, duration):
    h = np.asarray(duration, dtype=float)
    hv = h[..., None]
    k1 = _velocity_at(fields, position, time, branch)
    p2 = position+.5*hv*k1
    _validate_field_stencil(fields, position, p2)
    k2 = _velocity_at(fields, p2, time+.5*h, branch)
    p3 = position+.5*hv*k2
    _validate_field_stencil(fields, position, p3)
    k3 = _velocity_at(fields, p3, time+.5*h, branch)
    p4 = position+hv*k3
    _validate_field_stencil(fields, position, p4)
    k4 = _velocity_at(fields, p4, time+h, branch)
    result = position+hv*(k1+2*k2+2*k3+k4)/6
    _validate_field_stencil(fields, position, result)
    return result


def inertial_curvature(fields, position, time=0., branch=1, h=1e-4):
    """Laboratory trajectory radius in R_LC, including explicit time change.

    Differentiate beta(x,t) along (x +/- h beta, t +/- h), a second-order
    material derivative. This is not magnetic-line curvature. The small
    roundoff floor 32*epsilon/h suppresses only derivatives indistinguishable
    from floating-point cancellation (a straight trajectory has rho=inf).
    The provider must be valid over the derivative stencil; interpolate a
    resolved continuous field and vary h and orbit step in convergence tests.
    """
    if not np.isfinite(h) or h <= 0:
        raise ValueError("Positive finite curvature step required")
    p = _vectors(position, "Position")
    time = np.broadcast_to(np.asarray(time, dtype=float), p.shape[:-1])
    v = _velocity_at(fields, p, time, branch)
    _validate_field_stencil(fields, p-h*v, p, p+h*v)
    vp = _velocity_at(fields, p+h*v, time+h, branch)
    vm = _velocity_at(fields, p-h*v, time-h, branch)
    acceleration = (vp-vm)/(2*h)
    kappa = np.linalg.norm(np.cross(v, acceleration), axis=-1)
    resolved = kappa > 32*np.finfo(float).eps/h
    return np.divide(1., kappa, out=np.full_like(kappa, np.inf), where=resolved)


def arrival_phase(position, time, photon_direction):
    """Fixed-observer retarded phase (t - n dot x)/(2*pi), modulo one.

    photon_direction is a unit lab direction toward the observer. The
    omitted observer-distance light travel time is a common phase offset.
    This function does not decide whether a given particle beams to n.
    """
    p, n = np.broadcast_arrays(_vectors(position, "Position"),
                               _vectors(photon_direction, "Photon direction"))
    if not np.allclose(np.linalg.norm(n, axis=-1), 1., rtol=0., atol=1e-12):
        raise ValueError("Photon direction must be a unit vector")
    time = np.broadcast_to(np.asarray(time, dtype=float), p.shape[:-1])
    if not np.all(np.isfinite(time)):
        raise ValueError("Emission time must be finite")
    return ((time-np.sum(n*p, axis=-1))/(2*np.pi)) % 1


def sky_phase(position, time, direction):
    """Phase of a stationary pattern rotating at Omega, observer azimuth 0.

    For a trajectory sample at actual lab time t, rotate its entire emission
    event by -atan2(vy,vx), including the corresponding time translation.
    Hence phi=(t-azimuth(v)-v dot x)/(2*pi). This extra azimuth subtraction
    assumes rotating-pattern symmetry and is not valid for arbitrary
    transient fields. At the spin axis the azimuthal phase is undefined.
    """
    n = _vectors(direction, "Direction")
    if np.any(np.hypot(n[..., 0], n[..., 1]) < 1e-14):
        raise ValueError("Sky phase is undefined for an axial photon")
    return (arrival_phase(position, time, n)
            - np.arctan2(n[..., 1], n[..., 0])/(2*np.pi)) % 1


@dataclass(frozen=True)
class OuterTraceConfig:
    dt: float = .005
    max_time: float = 10.
    r_min: float = .05
    r_max: float = 2.5
    curvature_step: float = 1e-4

    def __post_init__(self):
        if (not all(np.isfinite(v) for v in asdict(self).values())
                or not (self.dt > 0 and self.max_time > 0
                        and 0 <= self.r_min < self.r_max
                        and self.curvature_step > 0)):
            raise ValueError("Invalid finite trace duration, bounds or derivative step")


def trace_outer(fields, seeds, branch=1, charge_sign=-1, *,
                config=OuterTraceConfig(), start_time=0., seed_weight=None):
    """RK4 lab trajectories with midpoint radiation segments and clipped exits.

    Returns padded arrays (segment,seed[,xyz]); dt is residence time in
    Omega^-1 and rho is in R_LC. e_work=charge_sign*E dot beta in statvolt/cm.
    Its sign is retained. No electric field is inferred from gamma or charge.
    'phase' assumes a rotating stationary pattern; 'fixed_phase' is the
    retarded time for each emitted ray along its own instantaneous direction.
    Axial rays have phase=NaN and sky_phase_valid=False, not an invented phase.

    The provider must cover each RK/curvature stencil, including a dt-wide
    buffer beyond termination surfaces. Discontinuous fields are unsupported;
    providers with validate_transport_stencil explicitly reject crossing RK
    or curvature stencils instead of producing unresolved interface radiation.
    r_max is spherical, not a new corotation limit. End statuses distinguish
    actual escape from 'horizon'; horizon does not prove physical escape.
    Branch stays fixed; no rule flips charge or branch to force acceleration.
    """
    p = _vectors(seeds, "Seeds").copy()
    if p.ndim == 1:
        p = p[None, :]
    if p.ndim != 2 or len(p) == 0:
        raise ValueError("Seeds must have shape (number_of_seeds,3)")
    n = len(p)
    signs = _signs(branch, (n,), "branch").copy()
    charge = _signs(charge_sign, (n,), "charge_sign").copy()
    start = np.broadcast_to(np.asarray(start_time, dtype=float), (n,)).copy()
    if not np.all(np.isfinite(start)):
        raise ValueError("Start times must be finite")
    radius = np.linalg.norm(p, axis=-1)
    if np.any(radius <= config.r_min) or np.any(radius >= config.r_max):
        raise ValueError("Seeds must be strictly inside the radial trace bounds")
    weights = (np.full(n, 1/n) if seed_weight is None
               else np.broadcast_to(np.asarray(seed_weight, dtype=float), (n,)).copy())
    if not np.all(np.isfinite(weights)) or np.any(weights < 0) or weights.sum() <= 0:
        raise ValueError("Seed weights must be finite, nonnegative, with positive sum")
    # Preserve supplied weights: their normalization belongs to injection.
    alive = np.ones(n, dtype=bool)
    status = np.full(n, "horizon", dtype="U7")
    elapsed = np.zeros(n)
    chunks = []

    def outside(x):
        r = np.linalg.norm(x, axis=-1)
        return (r <= config.r_min) | (r >= config.r_max)

    for iteration in range(int(np.ceil(config.max_time/config.dt))):
        if not np.any(alive):
            break
        ii = np.flatnonzero(alive)
        x, sign = p[ii], signs[ii]
        t = start[ii]+elapsed[ii]
        h = np.minimum(config.dt, config.max_time-elapsed[ii])
        end = _rk4(fields, x, t, sign, h)
        escaped = outside(end)
        if np.any(escaped):
            low = np.zeros(escaped.sum())
            high = h[escaped].copy()
            for _ in range(32):
                half = (low+high)/2
                out = outside(_rk4(fields, x[escaped], t[escaped], sign[escaped], half))
                high = np.where(out, half, high)
                low = np.where(out, low, half)
            h[escaped] = low
            end[escaped] = _rk4(fields, x[escaped], t[escaped], sign[escaped], low)
        mid = _rk4(fields, x, t, sign, h/2)
        tm = t+h/2
        e, b = fields(mid, tm)
        v = prescribed_velocity(e, b, sign)
        rho = inertial_curvature(fields, mid, tm, sign, config.curvature_step)
        work = charge[ii]*np.sum(e*v, axis=-1)
        # No inactive padded field queries, derivatives, or fabricated emission.
        sample = dict(position=p.copy(), direction=np.zeros_like(p),
                      rho=np.full(n, np.inf), dt=np.zeros(n),
                      time=start+elapsed, e_work=np.zeros(n),
                      phase=np.full(n, np.nan), fixed_phase=np.full(n, np.nan),
                      sky_phase_valid=np.zeros(n, dtype=bool))
        sample["position"][ii] = mid
        sample["direction"][ii] = v
        sample["rho"][ii] = rho
        sample["dt"][ii] = h
        sample["time"][ii] = tm
        sample["e_work"][ii] = work
        sample["fixed_phase"][ii] = arrival_phase(mid, tm, v)
        valid = np.hypot(v[:, 0], v[:, 1]) >= 1e-14
        sample["sky_phase_valid"][ii] = valid
        if np.any(valid):
            sample["phase"][ii[valid]] = sky_phase(mid[valid], tm[valid], v[valid])
        chunks.append(sample)
        p[ii] = end
        elapsed[ii] += h
        if np.any(escaped):
            r = np.linalg.norm(end[escaped], axis=-1)
            status[ii[escaped]] = np.where(abs(r-config.r_min) < abs(r-config.r_max),
                                           "r_min", "r_max")
            alive[ii[escaped]] = False
    result = {key: np.stack([sample[key] for sample in chunks]) for key in chunks[0]}
    result.update(seed_weight=weights, branch=signs, charge_sign=charge,
                  end_position=p, end_time=start+elapsed, status=status,
                  all_escaped=bool(np.all(status != "horizon")),
                  config=asdict(config), velocity_model="K14 Eq.12",
                  population_is_predicted=False)
    return result


def _primitive_low(y, reciprocal=False):
    """Stable integrals H(y)=int_0^y 1/(1+x^4), J(y)=int_0^y x^2/(1+x^4).

    Only 0<=y<=1 is evaluated. J(1/y) is the high-y tail of H, avoiding
    subtraction of nearly equal limiting primitives in strong cooling.
    """
    y = np.asarray(y, dtype=float)
    root2 = np.sqrt(2.)
    angle = np.arctan2(root2*y, 1-y*y)/(2*root2)
    logarithm = np.log1p(2*root2*y/(y*y-root2*y+1))/(4*root2)
    result = angle-logarithm if reciprocal else angle+logarithm
    # Direct difference for J loses precision at small argument.
    small = y < .1
    power = 3 if reciprocal else 1
    series = sum((-1.)**k*y**(power+4*k)/(power+4*k) for k in range(7))
    return np.where(small, series, result)


def _cooling_integral(lower, upper):
    """int_lower^upper dx/(1+x^4), for nonnegative lower<=upper."""
    lo, hi = np.broadcast_arrays(lower, upper)
    # Clipped arguments make all branches well-defined during vectorization.
    hlo = _primitive_low(np.minimum(lo, 1.))
    hhi = _primitive_low(np.minimum(hi, 1.))
    jlo = _primitive_low(1/np.maximum(lo, 1.), True)
    jhi = _primitive_low(1/np.maximum(hi, 1.), True)
    return np.where(hi <= 1, hhi-hlo,
                    np.where(lo >= 1, jlo-jhi,
                             _primitive_low(1.)-hlo+_primitive_low(1., True)-jhi))


def _advance_signed(gamma, acceleration, loss, ds, gamma_floor):
    """Frozen signed A-B gamma^4 flow, with an explicit validity boundary."""
    g, a, b, h = np.broadcast_arrays(gamma, acceleration, loss, ds)
    result = _advance_constant(g, np.maximum(a, 0.), b, h)
    negative = (a < 0) & (h > 0)
    linear = negative & (b == 0)
    result[linear] = g[linear]+a[linear]*h[linear]
    mixed = negative & (b > 0)
    if np.any(mixed):
        gm, am, bm, hm = g[mixed], -a[mixed], b[mixed], h[mixed]
        scale = (am/bm)**.25
        target = am*hm/scale
        upper = gm/scale
        minimum = gamma_floor/scale
        available = _cooling_integral(minimum, upper)
        if np.any(target > available):
            raise ValueError("Signed electric work/cooling crosses gamma_floor; prescribed-speed model invalid")
        lo, hi = minimum.copy(), upper.copy()
        for _ in range(56):
            mid = (lo+hi)/2
            time_to_mid = _cooling_integral(mid, upper)
            too_low = time_to_mid > target
            lo = np.where(too_low, mid, lo)
            hi = np.where(too_low, hi, mid)
        result[mixed] = scale*(lo+hi)/2
    if not np.all(np.isfinite(result)) or np.any(result < gamma_floor):
        raise ValueError("Signed electric work/cooling crosses gamma_floor; prescribed-speed model invalid")
    return result


def signed_transport_curvature(ds_cm, rho_cm, e_work_statvolt_cm, gamma0,
                               energy_bands_gev=None, spectral_nodes=1,
                               gamma_floor=10.):
    """Signed energy work and curvature losses along prescribed trajectories.

    d gamma/ds = e*E_work/(m_e*c^2)-2*e^2*gamma^4/(3*m_e*c^2*rho^2).
    E_work already includes the charge sign and the true velocity projection.
    A negative value decelerates. No absolute value or silent gamma clamping
    is allowed. Crossing gamma_floor raises: the caller must change the
    approximation or explicitly truncate the modeled stream. Default 10 is
    only an ultrarelativistic validity guard, not a physical injection model.

    Fields and curvature are frozen per segment. The positive branch reuses
    the existing exact solver; negative mixed work uses the monotone analytic
    cooling integral. Geometric and spectral quadrature require convergence.
    Returns the existing CurvatureTransport, with a conservative per-segment
    work/radiation/particle-energy ledger, in Gaussian CGS units.
    """
    if spectral_nodes not in (1, 3):
        raise ValueError("spectral_nodes must be 1 or 3")
    if not np.isfinite(gamma_floor) or gamma_floor < 1:
        raise ValueError("gamma_floor must be finite and at least 1")
    arrays = [np.asarray(v, dtype=float) for v in (ds_cm, rho_cm, e_work_statvolt_cm)]
    g0 = np.asarray(gamma0, dtype=float)
    ndim = max(1+g0.ndim, *(a.ndim for a in arrays), 1)
    arrays = [a.reshape(a.shape+(1,)*(ndim-1)) if a.ndim == 1 and ndim > 1
              else a for a in arrays]
    ds, rho, work, _ = np.broadcast_arrays(*arrays, g0[None, ...])
    if (np.any(~np.isfinite(ds)) or np.any(ds < 0)
            or np.any(np.isnan(rho)) or np.any(rho <= 0)
            or np.any(~np.isfinite(work)) or np.any(~np.isfinite(g0))
            or np.any(g0 < gamma_floor)):
        raise ValueError("Require ds>=0, rho>0, finite signed work and gamma0>=gamma_floor")
    a = E_STATC*work/MEC2_ERG
    b = 2*E_STATC**2/(3*MEC2_ERG*rho**2)
    gamma = np.empty((ds.shape[0]+1,)+ds.shape[1:])
    gamma[0] = g0
    gamma_mid = np.empty(ds.shape)
    work_erg = E_STATC*work*ds
    radiated = np.empty(ds.shape)
    ec = np.empty(ds.shape)
    photons = None if energy_bands_gev is None else np.empty(ds.shape+(len(energy_bands_gev),))
    fractions = .5*(1+np.array([-np.sqrt(3/5), 0., np.sqrt(3/5)]))
    weights = np.array([5/18, 4/9, 5/18])
    for i in range(ds.shape[0]):
        g, h = gamma[i], ds[i]
        gamma[i+1] = _advance_signed(g, a[i], b[i], h, gamma_floor)
        gm = _advance_signed(g, a[i], b[i], h/2, gamma_floor)
        gamma_mid[i] = gm
        radiated[i] = np.maximum(0., work_erg[i]-MEC2_ERG*(gamma[i+1]-g))
        # Straight segments have exactly zero radiation; suppress only their
        # subtraction roundoff, rather than fabricating tiny photon yields.
        radiated[i] = np.where(np.isinf(rho[i]), 0., radiated[i])
        ec[i] = critical_energy_gev(gm, rho[i])
        if photons is not None:
            if spectral_nodes == 1:
                photons[i] = curvature_photon_yield(radiated[i], ec[i], energy_bands_gev)
            else:
                node_g = np.stack([_advance_signed(g, a[i], b[i], h*f, gamma_floor)
                                   for f in fractions])
                node_ec = critical_energy_gev(node_g, rho[i])
                weighted = curvature_power_erg_s(node_g, rho[i])*weights.reshape((3,)+(1,)*g.ndim)
                normal = weighted.sum(axis=0)
                energy_fraction = np.divide(weighted, normal, out=np.zeros_like(weighted), where=normal > 0)
                photons[i] = curvature_photon_yield(
                    radiated[i]*energy_fraction, node_ec, energy_bands_gev).sum(axis=0)
    return CurvatureTransport(gamma, gamma_mid, work_erg, radiated, ec, photons)
