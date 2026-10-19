"""Explicitly normalized field snapshots and a rotating-pattern approximation.

Coordinates and time are x/R_LC and Omega*t. Field arrays must already be
collocated, have shape (3,Nx,Ny,Nz), and share a single E/B normalization
(Gaussian c=1 convention). No conversion of staggered solver data is hidden
here. A rotating snapshot is an approximation that needs a periodicity test.
Interpolation verification is not validation of the supplied magnetosphere.
"""
from dataclasses import dataclass
from pathlib import Path
import hashlib
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import map_coordinates,spline_filter


def rotate_z(points,angle):
    points=np.asarray(points,float)
    c,s=np.cos(angle),np.sin(angle)
    return np.stack((c*points[...,0]-s*points[...,1],
                     s*points[...,0]+c*points[...,1],points[...,2]),axis=-1)


@dataclass(frozen=True)
class SpinDownScale:
    """Model-derived dipole scale, not a direct measurement of Vela's field.

    Uses L=4 pi^2 I Pdot/P^3 and the force-free oblique-dipole torque
    L=mu^2 Omega^4/c^3 (1+sin^2 alpha), Spitkovsky (2006).
    I and the force-free torque are assumptions. No stellar radius required.
    """
    period_s:float
    period_derivative:float
    moment_of_inertia_g_cm2:float
    inclination_deg:float

    def __post_init__(self):
        values=(self.period_s,self.period_derivative,self.moment_of_inertia_g_cm2,self.inclination_deg)
        if not all(np.isfinite(values)) or min(values[:3])<=0 or not 0<=self.inclination_deg<=180:
            raise ValueError('Require finite positive timing/inertia and physical inclination')

    @property
    def luminosity_erg_s(self):
        return 4*np.pi**2*self.moment_of_inertia_g_cm2*self.period_derivative/self.period_s**3

    @property
    def rlc_cm(self):return 29979245800.*self.period_s/(2*np.pi)

    @property
    def dipole_scale_gauss(self):
        omega=2*np.pi/self.period_s
        return np.sqrt(self.luminosity_erg_s*omega**2/
                       (29979245800.**3*(1+np.sin(np.deg2rad(self.inclination_deg))**2)))


class RotatingSnapshot:
    """Interpolate a rigidly rotating field pattern with explicit physical scale.

    E/B are multiplied by amplitude_gauss exactly once on return. Cubic
    interpolation is a tensor-product B-spline; it does not preserve discrete
    divergence or magnetic dominance. Both remain measured diagnostics, not
    silently repaired. Linear interpolation is available for sensitivity.
    No extrapolation is permitted; cubic queries exclude a three-cell border.
    """
    def __init__(self,axes,e,b,*,amplitude_gauss,reference_time,
                 pattern_omega=1.,interpolation='cubic',verified=False,
                 provenance=None):
        self.axes=tuple(np.asarray(a,float).copy() for a in axes)
        if len(self.axes)!=3 or any(a.ndim!=1 or len(a)<8 or not np.all(np.isfinite(a))
            or np.any(np.diff(a)<=0) or not np.allclose(np.diff(a),np.diff(a)[0],rtol=1e-10,atol=1e-14)
            for a in self.axes):
            raise ValueError('Require three uniform increasing axes, at least eight cells')
        shape=(3,)+tuple(len(a) for a in self.axes)
        self.e=np.asarray(e,float);self.b=np.asarray(b,float)
        if self.e.shape!=shape or self.b.shape!=shape or not np.all(np.isfinite(self.e)) or not np.all(np.isfinite(self.b)):
            raise ValueError('Finite collocated E and B arrays must match axes')
        if not np.isfinite(amplitude_gauss) or amplitude_gauss<=0 or not np.isfinite(reference_time) or not np.isfinite(pattern_omega):
            raise ValueError('Explicit positive field scale and finite pattern time/rate required')
        if interpolation not in ('linear','cubic'):raise ValueError('Interpolation must be linear or cubic')
        self.amplitude_gauss=float(amplitude_gauss);self.reference_time=float(reference_time)
        self.pattern_omega=float(pattern_omega);self.interpolation=interpolation
        self.verified=bool(verified);self.provenance={} if provenance is None else dict(provenance)
        self.h=np.array([a[1]-a[0] for a in self.axes]);self.start=np.array([a[0] for a in self.axes])
        margin=3 if interpolation=='cubic' else 0
        self.minimum=np.array([a[margin] for a in self.axes]);self.maximum=np.array([a[-1-margin] for a in self.axes])
        if interpolation=='cubic':
            self._coefficients=tuple(tuple(spline_filter(v,order=3,mode='mirror') for v in f) for f in (self.e,self.b))
        else:
            self._linear=tuple(RegularGridInterpolator(self.axes,f.transpose(1,2,3,0),bounds_error=True) for f in (self.e,self.b))

    @classmethod
    def from_npz(cls,path,*,amplitude_gauss,interpolation='cubic'):
        import json
        path=Path(path)
        with np.load(path,allow_pickle=False) as d:
            axes=tuple(d[k] for k in ('x','y','z')) if 'x' in d else (d['axis'],)*3
            meta=json.loads(str(d['metadata'])) if 'metadata' in d else {}
            if 'metadata' in d and meta.get('coordinate_units')!='R_LC':
                raise ValueError('Canonical snapshot must explicitly use R_LC coordinates')
            # Existing local prototypes have normalized coordinates by design;
            # absence of verification metadata never grants release status.
            return cls(axes,d['e'],d['b'],amplitude_gauss=amplitude_gauss,
                       reference_time=float(d['time']),interpolation=interpolation,
                       verified=bool(meta.get('emission_release_allowed',False)),
                       provenance={**meta,'path':str(path.resolve()),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})

    def __call__(self,position,time):
        p=np.asarray(position,float)
        if p.shape[-1:]!=(3,) or not np.all(np.isfinite(p)) or not np.all(np.isfinite(time)):
            raise ValueError('Finite Cartesian positions and time required')
        try: t=np.broadcast_to(np.asarray(time,float),p.shape[:-1])
        except ValueError as exc:raise ValueError('Time must broadcast to the position batch') from exc
        angle=self.pattern_omega*(t-self.reference_time)
        q=rotate_z(p,-angle)
        if np.any(q<self.minimum) or np.any(q>self.maximum):
            raise ValueError('Field query leaves interpolation domain; extrapolation forbidden')
        flat=q.reshape(-1,3);result=[]
        for k in range(2):
            if self.interpolation=='cubic':
                index=((flat-self.start)/self.h).T
                value=np.stack([map_coordinates(c,index,order=3,mode='mirror',prefilter=False)
                                for c in self._coefficients[k]],axis=-1)
            else:value=self._linear[k](flat)
            result.append(rotate_z(value.reshape(p.shape),angle)*self.amplitude_gauss)
        return tuple(result)

    def require_emission_release(self):
        if not self.verified:
            raise RuntimeError('Snapshot has not passed field-release checks; only diagnostic transport is allowed')


def sample_quality(field,position,time,delta=1e-4):
    """Measure interpolated divB and E.B without projecting either away."""
    p=np.asarray(position,float);e,b=field(p,time)
    div=np.zeros(p.shape[:-1])
    for j in range(3):
        step=np.zeros(3);step[j]=delta
        bp=field(p+step,time)[1];bm=field(p-step,time)[1]
        div+=(bp[...,j]-bm[...,j])/(2*delta)
    eb=np.sum(e*b,axis=-1);b2=np.sum(b*b,axis=-1);e2=np.sum(e*e,axis=-1)
    return {'divB_per_RLC':div,'E_dot_B':eb,'B2_minus_E2':b2-e2}
