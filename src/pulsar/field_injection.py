"""Injection surfaces traced in the supplied numerical magnetic field.

No vacuum-dipole polar-cap rim is reused. This implements a star-shaped
open-cap approximation in magnetic coordinates on an explicit finite-radius
injection surface, NOT a pair-production model or a true-surface cascade.
Both poles are traced separately. The approximation and seed radius require
sensitivity checks before interpreting emission as a Vela prediction.
"""
from dataclasses import dataclass
import numpy as np
from .trajectories import seed_positions
from .snapshot_fields import rotate_z


@dataclass(frozen=True)
class CapConfig:
    inclination_deg:float=75.
    seed_radius:float=.35
    inner_radius:float=.34  # Retained transport/metadata radius; not the cap-return surface.
    field_time:float=0.
    magnetic_phase:float=0.
    step:float=.008
    max_length:float=15.
    outer_radius:float=2.5
    azimuths:int=64
    rings:int=3
    annulus_inner:float=.90
    annulus_outer:float=.96
    bisections:int=14
    bracket_samples:int=33
    rim_boundary:str='light_cylinder'

    def __post_init__(self):
        values=(self.inclination_deg,self.seed_radius,self.inner_radius,self.field_time,self.magnetic_phase,
                self.step,self.max_length,self.outer_radius,self.annulus_inner,self.annulus_outer)
        if not all(np.isfinite(values)) or not (0<=self.inclination_deg<=90 and
            0<self.inner_radius<self.seed_radius<1<self.outer_radius and
            self.step>0 and self.max_length>self.outer_radius and
            0<self.annulus_inner<self.annulus_outer<=1):
            raise ValueError('Invalid cap geometry')
        if any(isinstance(v,bool) or not isinstance(v,int) for v in (self.azimuths,self.rings,self.bisections,self.bracket_samples)) or self.azimuths<4 or self.rings<1 or self.bisections<4 or self.bracket_samples<9:
            raise ValueError('Invalid cap quadrature')
        if self.rim_boundary not in ('light_cylinder','outgoing_sphere'):
            raise ValueError('Unknown finite-domain cap boundary criterion')


def _classify(field,theta,azimuth,pole,cfg):
    """Trace magnetic lines until escape or re-entry through their launch sphere.

    A cap on r=seed_radius describes the outward-connected field there. A
    line that returns through that same sphere is closed for this mapping,
    even if it never reaches the separate inner transport cutoff. The small
    radial hysteresis avoids calling roundoff or a tangent orbit a re-entry.
    It does not assign an outcome to a line that remains unresolved.
    """
    q=rotate_z(pole*seed_positions(theta,azimuth,np.deg2rad(cfg.inclination_deg),cfg.seed_radius),cfg.magnetic_phase)
    alive=np.ones(len(q),bool);result=np.zeros(len(q),int)
    return_radius=cfg.seed_radius*(1.-1e-6)
    # 1: crossed the light cylinder, 2: escaped sphere without crossing it,
    # -1: re-entered launch sphere. Unresolved/null/domain exits raise.
    def direction(p):
        b=field(p,cfg.field_time)[1];bn=np.linalg.norm(b,axis=-1)
        if np.any(bn<=1e-30):raise ValueError('Magnetic null encountered in cap trace')
        return pole*b/bn[:,None]
    for _ in range(int(np.ceil(cfg.max_length/cfg.step))):
        if not alive.any():break
        ii=np.flatnonzero(alive);p=q[ii];h=cfg.step
        k1=direction(p);k2=direction(p+.5*h*k1);k3=direction(p+.5*h*k2);k4=direction(p+h*k3)
        q[ii]=p+h/6*(k1+2*k2+2*k3+k4)
        r=np.linalg.norm(q[ii],axis=-1);cyl=np.linalg.norm(q[ii,:2],axis=-1)
        cross=cyl>=1.;far=r>=cfg.outer_radius;returned=r<=return_radius
        if cfg.rim_boundary=='outgoing_sphere' and np.any(far&~cross):
            escaped=q[ii[far&~cross]]
            if np.any(np.sum(direction(escaped)*escaped,axis=-1)<=0):
                raise RuntimeError('Sphere escape is not outward along the traced branch')
        result[ii[returned]]=-1;result[ii[far]]=2;result[ii[cross]]=1
        alive[ii[cross|far|returned]]=False
    if alive.any():raise RuntimeError('Unclassified magnetic lines; increase trace horizon')
    return result


def trace_cap_rim(field,azimuth,pole,cfg):
    """Find a single connected cap; reject sampled nonmonotonic topology.

    By default the final rim requires an actual cylinder crossing. The
    explicit outgoing_sphere option maps a finite-domain boundary-open
    region instead: a traced branch may leave the outer sphere outward
    before crossing the cylinder. It does not establish opening at infinity.
    In that mode snapshot_injection requires every selected annular seed
    to reach the same outgoing finite-domain boundary separately.
    """
    az=np.asarray(azimuth,float)
    if pole not in (-1,1):raise ValueError('Pole must be +/-1')
    if az.ndim!=1 or not np.all(np.isfinite(az)):raise ValueError('Finite azimuth vector required')
    # A swept-back cap on a finite-radius surface can extend past 90 degrees
    # on some spokes. Find each spoke's actual transition rather than impose
    # a fixed maximum cap angle. Sample the full great-circle spoke first;
    # reject multiple open intervals instead of selecting a convenient one.
    angles=np.linspace(.005,np.pi-.005,cfg.bracket_samples)
    sampled=[]
    for angle in angles:
        sampled.append(_classify(field,np.full(az.shape,angle),az,pole,cfg)>0)
    sampled=np.array(sampled)
    if not sampled[0].all() or sampled[-1].any():raise RuntimeError('Unable to bracket open cap in supplied field')
    if np.any(np.diff(sampled.astype(int),axis=0)>0):raise RuntimeError('Open cap is not star-shaped on sampled spokes')
    first_closed=np.argmax(~sampled,axis=0)
    lo=angles[first_closed-1];hi=angles[first_closed]
    for _ in range(cfg.bisections):
        mid=(lo+hi)/2;opened=_classify(field,mid,az,pole,cfg)>0
        lo=np.where(opened,mid,lo);hi=np.where(opened,hi,mid)
    final=_classify(field,lo,az,pole,cfg)
    if cfg.rim_boundary=='light_cylinder' and np.any(final!=1):
        raise RuntimeError('Final rim not established by light-cylinder crossing')
    if cfg.rim_boundary=='outgoing_sphere' and np.any(final<=0):
        raise RuntimeError('Final rim not established by an outgoing finite-domain boundary')
    return (lo+hi)/2


def snapshot_injection(field,cfg,*,loading='magnetic_flux'):
    """Return finite-radius seeds, signed branches and geometric weights.

    magnetic_flux assigns equal injection per unsigned magnetic flux.
    surface_area assigns equal injection per surface area. Neither predicts
    actual net charge, primary current or pair multiplicity. Weights sum to
    unity, so absolute population normalization must be supplied separately.
    Annulus fractions are angular fractions of this surface's traced rim;
    they are not claimed identical to stellar-surface r_ovc in a paper.
    """
    if loading not in ('magnetic_flux','surface_area'):raise ValueError('Unknown injection loading')
    az=(np.arange(cfg.azimuths)+.5)*2*np.pi/cfg.azimuths
    du=(cfg.annulus_outer-cfg.annulus_inner)/cfg.rings
    u=cfg.annulus_inner+(np.arange(cfg.rings)+.5)*du
    seeds=[];weights=[];branches=[];rims=[]
    for pole in (1,-1):
        rim=trace_cap_rim(field,az,pole,cfg);rims.append(rim)
        theta=(u[:,None]*rim).ravel();phi=np.tile(az,cfg.rings)
        p=rotate_z(pole*seed_positions(theta,phi,np.deg2rad(cfg.inclination_deg),cfg.seed_radius),cfg.magnetic_phase)
        if cfg.rim_boundary=='outgoing_sphere' and np.any(_classify(field,theta,phi,pole,cfg)<=0):
            raise RuntimeError('Selected annular seeds must individually reach an outgoing finite-domain boundary')
        area=cfg.seed_radius**2*np.sin(theta)*np.tile(rim,cfg.rings)*du*(2*np.pi/cfg.azimuths)
        b=field(p,cfg.field_time)[1];br=np.sum(b*p,axis=-1)/cfg.seed_radius
        if np.any(pole*br<=0):raise RuntimeError('Inward/tangent B at selected injection seeds')
        weights.append(area*np.abs(br) if loading=='magnetic_flux' else area)
        seeds.append(p);branches.append(np.full(len(p),pole))
    weights=np.concatenate(weights)
    if not np.isfinite(weights).all() or weights.sum()<=0:raise ValueError('Invalid injection weights')
    return {'position':np.concatenate(seeds),'branch':np.concatenate(branches),
            'seed_weight':weights/weights.sum(),'cap_rim':np.array(rims),'cap_azimuth':az,
            'loading':loading,'population_is_predicted':False,
            'surface_mapping':'finite-radius angular fraction, not stellar-surface r_ovc'}
