"""Independent physical identities and finite-difference verification."""
from dataclasses import replace
import numpy as np
from .model import (TAU, Parameters, fields, sheet_function, charge_current_4pi,
                    emitting_velocity, analytic_separation, peak_separation,
                    pulse_profile, observer_direction)


def differential_checks(points, t, chi, h):
    e, b = fields(points, t, chi)
    de = np.empty((len(points), 3, 3))
    db = np.empty_like(de)
    for axis in range(3):
        offset = np.eye(3)[axis] * h
        ep, bp = fields(points+offset, t, chi)
        em, bm = fields(points-offset, t, chi)
        de[:, axis, :] = (ep-em)/(2*h)
        db[:, axis, :] = (bp-bm)/(2*h)
    ep, bp = fields(points, t+h, chi)
    em, bm = fields(points, t-h, chi)
    dt_e, dt_b = (ep-em)/(2*h), (bp-bm)/(2*h)

    def curl(d):
        return np.stack((d[:, 1, 2]-d[:, 2, 1], d[:, 2, 0]-d[:, 0, 2],
                         d[:, 0, 1]-d[:, 1, 0]), axis=-1)

    rho, current = charge_current_4pi(points, t, chi)
    scale = np.linalg.norm(b, axis=-1) / np.linalg.norm(points, axis=-1)
    vecnorm = lambda v: float(np.max(np.linalg.norm(v, axis=-1)/scale))
    scalnorm = lambda v: float(np.max(np.abs(v)/scale))
    return {
        "div_B": scalnorm(np.trace(db, axis1=1, axis2=2)),
        "gauss_E": scalnorm(np.trace(de, axis1=1, axis2=2)-rho),
        "faraday": vecnorm(curl(de)+dt_b),
        "ampere": vecnorm(curl(db)-dt_e-current),
    }


def field_validation(p=Parameters()):
    rng = np.random.default_rng(20261012)
    er = rng.normal(size=(6000, 3)); er /= np.linalg.norm(er, axis=1)[:, None]
    xyz = er * rng.uniform(p.r_min, p.r_max, (6000, 1))
    t = 0.731
    xyz = xyz[np.abs(sheet_function(xyz, t, p.chi)) > 0.15][:512]
    e, b = fields(xyz, t, p.chi)
    b2 = np.sum(b*b, axis=1)
    e2 = np.sum(e*e, axis=1)
    r = np.linalg.norm(xyz, axis=1)
    rho, current = charge_current_4pi(xyz, t, p.chi)
    v = emitting_velocity(xyz, p.lorentz_factor)
    identities = {
        "E_dot_B_relative": float(np.max(np.abs(np.sum(e*b, axis=1))/b2)),
        "B2_minus_E2_identity_relative": float(np.max(np.abs((b2-e2)*r**4-1))),
        "force_free_relative": float(np.max(np.linalg.norm(rho[:, None]*e+np.cross(current,b), axis=1)/(b2/r))),
        "velocity_magnitude_error": float(np.max(np.abs(np.sum(v*v,axis=1)-(1-p.lorentz_factor**-2)))),
        "ideal_velocity_relative": float(np.max(np.linalg.norm(e+np.cross(v,b),axis=1)/np.sqrt(b2))),
    }
    resolutions = []
    for h in (0.008, 0.004, 0.002):
        resolutions.append({"step_RLC": h, **differential_checks(xyz,t,p.chi,h)})
    orders = {key: float(np.log2(resolutions[-2][key]/resolutions[-1][key]))
              for key in ("div_B", "gauss_E", "faraday", "ampere")}
    # Poynting luminosity in these units: integral sin^2(theta)/(4*pi) dOmega=2/3.
    mu, weight = np.polynomial.legendre.leggauss(48)
    ph = TAU*(np.arange(96)+0.5)/96
    er = np.stack(np.broadcast_arrays(np.sqrt(1-mu[:,None]**2)*np.cos(ph),
                   np.sqrt(1-mu[:,None]**2)*np.sin(ph), mu[:,None]), axis=-1)
    fluxes, net_b = [], []
    for rr in (1.5, 3., 5.):
        ee, bb = fields(rr*er, t, p.chi)
        angular_w = weight[:,None]*TAU/96
        fluxes.append(float(np.sum(np.sum(np.cross(ee,bb)*er,axis=-1)*angular_w)*rr**2/(4*np.pi)))
        net_b.append(float(np.sum(np.sum(bb*er,axis=-1)*angular_w)*rr**2))

    # Maxwell jump conditions, tested directly on the moving sheet.
    rr = np.linspace(1.5,5.,160)
    eta = np.linspace(0.1,TAU+0.1,160,endpoint=False)
    phase = t-rr
    tangent1 = np.stack((np.cos(p.chi)*np.cos(phase),np.cos(p.chi)*np.sin(phase),-np.sin(p.chi)*np.ones_like(rr)),axis=-1)
    tangent2 = np.stack((-np.sin(phase),np.cos(phase),np.zeros_like(rr)),axis=-1)
    surface = rr[:,None]*(np.cos(eta)[:,None]*tangent1+np.sin(eta)[:,None]*tangent2)
    h=1e-5
    grad = np.stack([(sheet_function(surface+np.eye(3)[k]*h,t,p.chi)-sheet_function(surface-np.eye(3)[k]*h,t,p.chi))/(2*h) for k in range(3)],axis=-1)
    ut=(sheet_function(surface,t+h,p.chi)-sheet_function(surface,t-h,p.chi))/(2*h)
    se,sb=fields(surface,t,p.chi,polarity=1.)
    jump_scale=np.linalg.norm(grad,axis=1)*np.linalg.norm(2*sb,axis=1)
    jump={
        "magnetic_normal_jump":float(np.max(np.abs(np.sum(grad*2*sb,axis=1))/jump_scale)),
        "moving_sheet_faraday_jump":float(np.max(np.linalg.norm(np.cross(grad,2*se)+ut[:,None]*2*sb,axis=1)/jump_scale)),
    }
    passed = (max(identities.values()) < 1e-11 and
              max(resolutions[-1][k] for k in orders) < 2e-5 and
              min(orders.values()) > 1.95 and
              max(abs(x-2/3) for x in fluxes) < 1e-12 and
              max(abs(x) for x in net_b) < 1e-12 and max(jump.values()) < 1e-8)
    return {"passed":bool(passed),"sample_points":len(xyz),"identities":identities,
            "finite_differences":resolutions,"observed_orders":orders,
            "poynting_luminosities":fluxes,"expected_luminosity":2/3,
            "net_magnetic_flux":net_b,"sheet_Maxwell_jumps":jump,
            "excludes":"Dissipative sheet force balance and inner dipolar magnetosphere"}


def timing_validation(p=Parameters()):
    n=observer_direction(p.zeta)
    rr=np.linspace(p.r_min,p.r_max,19)
    xyz=rr[:,None]*n
    phases=np.linspace(0,1,16384,endpoint=False)
    # Pencil-beam calculation using independently evaluated emission events.
    los=sheet_function(xyz[None,:,:], TAU*phases[:,None]+rr, p.chi)
    expected=np.cos(p.chi)*np.cos(p.zeta)+np.sin(p.chi)*np.sin(p.zeta)*np.cos(TAU*phases)
    error=float(np.max(np.abs(los-expected[:,None])))
    flux=np.exp(-0.5*(los[:,0]/0.002)**2)
    separation, peaks=peak_separation(phases,flux)
    predicted=analytic_separation(p.chi,p.zeta)
    no_delays=sheet_function(xyz[None,:,:],TAU*phases[:,None],p.chi)
    bad_error=float(np.max(np.abs(no_delays-expected[:,None])))
    # Repeatability and aligned-rotator invariance exercise the volume integrator.
    ph=np.arange(32)/32
    a=pulse_profile(p,ph,(6,20,40))
    b=pulse_profile(p,ph+1,(6,20,40))
    periodic=float(np.max(np.abs(a-b))/np.max(a))
    aligned=pulse_profile(replace(p,obliquity_deg=0),ph,(6,20,40))
    aligned_error=float(np.ptp(aligned)/np.mean(aligned))
    passed=error<1e-12 and abs(separation-predicted)<2/len(phases) and bad_error>0.1 and periodic<1e-12 and aligned_error<1e-12
    return {"passed":bool(passed),"arrival_phase_error":error,
            "analytic_pencil_beam_separation":predicted,"numerical_pencil_beam_separation":separation,
            "pencil_beam_peak_phases":peaks,"wrong_delay_negative_control_error":bad_error,
            "periodicity_relative_error":periodic,"aligned_rotator_variability":aligned_error}
