"""Published macroscopic acceleration prescriptions with explicit assumptions.

NearForceFreeFIDO implements Kalapotharakos et al. (2014), Eq.11 and Sec.V:
force-free geometry, E_parallel from curl(B)/(4*pi*sigma), screened inside
the light cylinder. This is a postprocessed near-FFE approximation, not a
self-consistent finite-conductivity field evolution or pair cascade.
https://arxiv.org/html/1310.3545v2

Barnard2022Work records its alternative prescribed-primary acceleration.
It is not added to FIDO; the two prescriptions are mutually exclusive.
https://arxiv.org/html/2111.03405v1
"""
from dataclasses import dataclass
import numpy as np
from .radiation import MEC2_ERG,E_STATC


def field_curl(field,position,time,step):
    """Fourth-order curl of B in normalized x/RLC coordinates."""
    if not np.isfinite(step) or step<=0:raise ValueError('Positive derivative step required')
    p=np.asarray(position,float);d=[]
    for axis in range(3):
        offset=np.zeros(3);offset[axis]=step
        derivative=(8*(field(p+offset,time)[1]-field(p-offset,time)[1])-
                    (field(p+2*offset,time)[1]-field(p-2*offset,time)[1]))/(12*step)
        d.append(derivative)
    return np.stack((d[1][...,2]-d[2][...,1],d[2][...,0]-d[0][...,2],d[0][...,1]-d[1][...,0]),axis=-1)


class NearForceFreeFIDO:
    """Signed parallel E from a supplied force-free field and chosen sigma/Omega.

    Base E is projected perpendicular to B. The parallel component is the
    signed projection of curlB/(4*pi*sigma_over_omega), enabled where the
    cylindrical distance from the spin axis is >= RLC. No absolute electric
    work, manual bright spots or polarity-dependent luminosity is inserted.
    Fields are Gaussian compatible; coordinates use RLC so factors c/Omega
    cancel. The sigma=30 case is the published 2014 reference prescription.
    The grid must resolve current gradients; derivatives amplify grid error.
    The cylindrical screening switch is sharp. Transport must invoke
    validate_transport_stencil and reject stencils that meet this interface;
    otherwise its velocity jump would generate artificial curvature radiation.
    """
    def __init__(self,base,*,sigma_over_omega=30.,derivative_step):
        if not np.isfinite(sigma_over_omega) or sigma_over_omega<=0:
            raise ValueError('Finite positive sigma/Omega required')
        if not np.isfinite(derivative_step) or derivative_step<=0:
            raise ValueError('Finite positive derivative step required')
        self.base=base;self.sigma_over_omega=float(sigma_over_omega)
        self.derivative_step=float(derivative_step)
        self.provenance={'prescription':'K14 near-FFE FIDO Eq.11',
                         'sigma_over_omega':self.sigma_over_omega,
                         'screening_boundary':'cylindrical radius=RLC',
                         'kinetic_pair_model':False}

    def __call__(self,position,time):
        p=np.asarray(position,float);e,b=self.base(p,time)
        b2=np.sum(b*b,axis=-1)
        if np.any(b2<=np.finfo(float).tiny):
            raise ValueError('Parallel-field direction undefined at a magnetic null')
        bhat=b/np.sqrt(b2)[...,None]
        perpendicular=e-np.sum(e*bhat,axis=-1)[...,None]*bhat
        curl=field_curl(self.base,p,time,self.derivative_step)
        parallel=np.sum(curl*bhat,axis=-1)/(4*np.pi*self.sigma_over_omega)
        parallel=np.where(np.hypot(p[...,0],p[...,1])>=1,parallel,0.)
        return perpendicular+parallel[...,None]*bhat,b

    def require_emission_release(self):
        self.base.require_emission_release()

    def validate_transport_stencil(self,positions):
        """Reject a trajectory/derivative stencil at the sharp interface.

        The leading axis lists stencil points for each independent stream.
        Check radial sides and the straight segments from the first point
        to every other point. The latter also catches an outside-to-outside
        chord that passes through the cylinder. This is a conservative
        validity gate, not a transition model or an automatic smoothing rule.
        """
        p=np.asarray(positions,float)
        if p.ndim<2 or p.shape[0]<2 or p.shape[-1]!=3 or not np.all(np.isfinite(p)):
            raise ValueError('Finite transport stencil with at least two positions required')
        radius=np.linalg.norm(p[...,:2],axis=-1)
        a=p[0,...,:2];delta=p[1:,...,:2]-a
        distance2=np.sum(delta*delta,axis=-1)
        fraction=np.clip(np.divide(-np.sum(a*delta,axis=-1),distance2,
                                  out=np.zeros_like(distance2),where=distance2>0),0.,1.)
        closest=np.linalg.norm(a+fraction[...,None]*delta,axis=-1)
        minimum=np.minimum(np.min(radius,axis=0),np.min(closest,axis=0))
        maximum=np.max(radius,axis=0)
        tolerance=32*np.finfo(float).eps
        if np.any((minimum<=1.+tolerance)&(maximum>=1.-tolerance)):
            raise ValueError('FIDO transport stencil crosses or touches the sharp screening interface; '
                             'a resolved transition model is required for curvature transport')


@dataclass(frozen=True)
class Barnard2022Work:
    """Effective accelerating work per path length from Barnard Eq.(2).

    The selected scenario-2 rates enter dgamma/dl = R_acc - curvature losses
    directly, without an additional field-parallel velocity projection.
    The returned effective work field is R_acc*m_e*c^2/e, not a global vector
    electric field. This prescribes acceleration per primary and does not
    infer charge populations, current closure or Maxwell-consistent E.
    """
    rate_inside_cm_inverse:float=.04
    rate_outside_cm_inverse:float=.25
    transition_radius_rlc:float=1.

    def __post_init__(self):
        if not all(np.isfinite([self.rate_inside_cm_inverse,self.rate_outside_cm_inverse,self.transition_radius_rlc])) or min(self.rate_inside_cm_inverse,self.rate_outside_cm_inverse)<0 or self.transition_radius_rlc<=0:
            raise ValueError('Require nonnegative finite rates and a positive transition radius')

    def work_field(self,position,direction,b):
        p,v,b=np.broadcast_arrays(np.asarray(position,float),np.asarray(direction,float),np.asarray(b,float))
        norm=np.linalg.norm(b,axis=-1)
        if np.any(norm<=0) or not np.all(np.isfinite([p,v,b])):
            raise ValueError('Finite positions/directions and nonzero B required')
        if not np.allclose(np.linalg.norm(v,axis=-1),1,rtol=1e-10,atol=1e-10):
            raise ValueError('Direction must be the unit laboratory velocity')
        rate=np.where(np.linalg.norm(p,axis=-1)<self.transition_radius_rlc,
                      self.rate_inside_cm_inverse,self.rate_outside_cm_inverse)
        return rate*MEC2_ERG/E_STATC
