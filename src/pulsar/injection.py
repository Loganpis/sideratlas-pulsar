"""Declared primary-loading hypotheses for the legacy rotating-dipole model.

Leading-order Goldreich--Julian charge density is rho_GJ ~ -Omega.B/(2 pi c).
The particle-number flux through a seed cell is therefore proportional to
abs(B_z) abs(v.er) dA when the spin axis is z and |rho_GJ|/e is used as a
single-sign primary density. This is NOT a pair-multiplicity prescription,
an exact current solution, or a source term for the evolved Maxwell solver.

Source: Bai & Spitkovsky (2010), Eq. 3, https://arxiv.org/pdf/0910.5741 .
Its exact corotation-density expression also includes beta.curl(B)/(4 pi),
which is omitted in goldreich_julian mode. The leading expression is an inner-magnetosphere
loading hypothesis at the chosen r_seed=0.08 R_LC, not a general current-sheet
density formula. Absolute current and luminosity disappear on normalization.

The corotation_charge option retains the curl(B) term exactly for the assumed
corotation electric field, to finite-difference accuracy. This computes its
required net charge; it still does not solve particle extraction, pair supply,
actual current closure, or plasma feedback on the legacy magnetic geometry.
"""
import numpy as np

from .trajectories import Geometry, dipole_field, seed_positions, velocity


def _corotation_charge(position, chi, step):
    """Return 4 pi rho in c=Omega=1 units, using fourth-order curl(B)."""
    derivatives = []
    for axis in range(3):
        offset = np.zeros(3)
        offset[axis] = step
        derivative = (8*(dipole_field(position+offset,chi)
                         -dipole_field(position-offset,chi))
                      -(dipole_field(position+2*offset,chi)
                        -dipole_field(position-2*offset,chi)))/(12*step)
        derivatives.append(derivative)
    curl_b = np.column_stack((derivatives[1][:,2]-derivatives[2][:,1],
                              derivatives[2][:,0]-derivatives[0][:,2],
                              derivatives[0][:,1]-derivatives[1][:,0]))
    beta = np.column_stack((-position[:,1],position[:,0],np.zeros(len(position))))
    return -2*dipole_field(position,chi)[:,2]+np.sum(beta*curl_b,axis=-1)


def seed_injection_weights(bundle, geometry, mode="goldreich_julian"):
    """Return normalized seed weights without changing ``bundle``.

    ``uniform_flux`` returns the stored legacy unsigned-magnetic-flux weights.
    ``goldreich_julian`` reconstructs the actual injection-surface quadrature,
    never the first trajectory midpoint. It computes |Bz| |v.er| dA directly,
    equivalent to multiplying flux weights by |Bz| |v.er| / |Br|. The direct
    calculation avoids repeated reweighting if called twice on the same data.
    Magnetic tangencies are rejected instead of hiding a singular multiplier.

    ``corotation_charge`` uses |div[-(Omega cross r) cross B]| |v.er| dA,
    retaining the curl(B) contribution to the charge required for corotation.
    Its field curl uses fourth-order central differences at the true seeds.

    All modes are assumptions about primary loading, not fitted phase weights
    or self-consistent predictions of the electron/positron populations.
    """
    if mode not in ("uniform_flux", "goldreich_julian", "corotation_charge"):
        raise ValueError("Unknown injection mode")
    old = np.asarray(bundle["seed_weight"], dtype=float)
    if (old.ndim != 1 or not len(old) or not np.all(np.isfinite(old))
            or np.any(old < 0) or not np.isclose(old.sum(), 1., rtol=1e-12, atol=0)):
        raise ValueError("Stored seed weights must be nonnegative and normalized")
    if mode == "uniform_flux":
        return old.copy()
    if not isinstance(geometry, Geometry):
        raise TypeError("geometry must be a trajectory Geometry")
    stored = bundle.get("geometry", {})
    for name in ("inclination_deg", "seed_radius", "annulus_inner", "annulus_outer"):
        if name in stored and not np.isclose(stored[name], getattr(geometry,name), rtol=0, atol=1e-14):
            raise ValueError("Injection geometry differs from the traced bundle")
    azimuth = np.asarray(bundle["cap_azimuth"], dtype=float)
    rim = np.asarray(bundle["cap_rim"], dtype=float)
    rings = int(bundle["n_rings"])
    nphi = int(bundle["n_azimuth"])
    if (rings != bundle["n_rings"] or nphi != bundle["n_azimuth"]
            or rings < 1 or nphi < 1 or azimuth.shape != (nphi,)
            or rim.shape != (nphi,) or len(old) != 2*rings*nphi
            or not np.all(np.isfinite(azimuth)) or not np.all(np.isfinite(rim))
            or np.any(rim <= 0) or np.any(rim >= np.pi/2)):
        raise ValueError("Malformed seed quadrature")
    # The tracer uses a uniform periodic azimuthal quadrature. Arbitrary
    # nonuniform cells require their own cell widths, not this reconstruction.
    spacing = (np.roll(azimuth,-1)-azimuth) % (2*np.pi)
    if nphi > 1 and not np.allclose(spacing,2*np.pi/nphi,rtol=0,atol=1e-12):
        raise ValueError("Seed azimuth quadrature must be uniform and ordered")
    width = (geometry.annulus_outer-geometry.annulus_inner)/rings
    fractional_ring = geometry.annulus_inner+(np.arange(rings)+.5)*width
    theta = (fractional_ring[:,None]*rim[None,:]).ravel()
    phi = np.broadcast_to(azimuth,(rings,nphi)).ravel()
    north = seed_positions(theta,phi,geometry.chi,geometry.seed_radius)
    seeds = np.concatenate((north,-north))
    expected_pole = np.r_[np.ones(len(north)),-np.ones(len(north))]
    pole = np.asarray(bundle.get("pole",expected_pole),dtype=float)
    if pole.shape != expected_pole.shape or not np.array_equal(pole,expected_pole):
        raise ValueError("Seed hemisphere ordering differs from the tracer")
    field = dipole_field(seeds,geometry.chi)
    radial = seeds/geometry.seed_radius
    br = np.sum(field*radial,axis=-1)
    if np.any(np.abs(br) <= 64*np.finfo(float).eps*np.linalg.norm(field,axis=-1)):
        raise ValueError("Magnetic field tangent to the injection surface")
    direction = velocity(seeds,geometry.chi,pole)[0]
    vr = np.sum(direction*radial,axis=-1)
    if np.any(vr <= 0):
        raise ValueError("Injection quadrature includes inward or tangent trajectories")
    area_north = geometry.seed_radius**2*np.sin(theta)*np.tile(rim,rings)*width*(2*np.pi/nphi)
    area = np.r_[area_north,area_north]
    charge = field[:,2]
    if mode == "corotation_charge":
        # Small relative step stays away from the point-dipole singularity
        # even in narrow-cap limit tests. No trajectory midpoint is sampled.
        step = min(1e-5,1e-4*geometry.seed_radius)
        charge = _corotation_charge(seeds,geometry.chi,step)
    weight = np.abs(charge)*np.abs(vr)*area
    total = weight.sum()
    if not np.all(np.isfinite(weight)) or not np.isfinite(total) or total <= 0:
        raise ValueError("Charge loading has no finite nonzero normalization")
    return weight/total
