"""Numerical resistive Maxwell prototype; not a kinetic/PIC magnetosphere.

Units c=Omega=R_LC=1, rationalized EM units, J=div(E) v_d+kappa E_parallel.
Kappa is the damping rate 4*pi*sigma/Omega for the Gaussian conductivity
used in the source paper, NOT that paper's sigma/Omega itself.
The closure is Kalapotharakos et al. (2014), Eq.9, with finite conductivity
inside/outside LC. E_parallel damping is operator split and exact at fixed B.
B=curl(A) and compatible periodic central derivatives enforce div(B)=0.
The rotating conductor, outer sponge and numerical filter are explicit,
separately measured interventions. This is a verification prototype, not
the much finer published production calculation.
"""
from dataclasses import dataclass,asdict
from pathlib import Path
import hashlib
import numpy as np

# Capture the source when this module is loaded, not later at save time.
SOLVER_SOURCE_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


@dataclass(frozen=True)
class PlasmaConfig:
    n: int = 64
    half_extent: float = 4.
    star_radius: float = .35
    inclination_deg: float = 65.
    conductivity_inside: float = 80.
    conductivity_outside: float = 5.
    conductivity_transition: float = .12
    sponge_start: float = 2.8
    sponge_strength: float = 4.
    filter_strength: float = .08
    courant: float = .30
    spinup_time: float = 2.
    plasma: bool = True
    conductor: bool = True
    inside_magnetic_limit: bool = True

    def __post_init__(self):
        if not all(np.isfinite(value) for value in asdict(self).values()):
            raise ValueError('All numerical parameters must be finite')
        if not 0<=self.inclination_deg<=180:
            raise ValueError('Inclination must lie between zero and 180 degrees')
        if self.n<16 or self.n%2 or not 0<self.courant<=.35:
            raise ValueError('Use even N>=16 and a stable Courant factor <=.35')
        if not (0<self.star_radius<1<self.sponge_start<self.half_extent
                and self.conductivity_inside>=0 and self.conductivity_outside>=0
                and 0<=self.filter_strength<=.3 and self.sponge_strength>=0
                and self.spinup_time>0 and self.conductivity_transition>0):
            raise ValueError('Invalid domain/closure parameters')


def derivative(a,axis,h):
    return (np.roll(a,-1,axis=axis)-np.roll(a,1,axis=axis))/(2*h)


def curl(a,h):
    return np.stack((derivative(a[2],1,h)-derivative(a[1],2,h),
                     derivative(a[0],2,h)-derivative(a[2],0,h),
                     derivative(a[1],0,h)-derivative(a[0],1,h)))


def divergence(a,h):
    return sum(derivative(a[j],j,h) for j in range(3))


def cross(a,b):
    return np.stack((a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]))


def invariants(e,b):
    b2=np.sum(b*b,axis=0); e2=np.sum(e*e,axis=0); eb=np.sum(e*b,axis=0)
    delta=b2-e2; d=np.hypot(delta,2*eb)
    # Stable positive root when B^2>E^2 and E_parallel is tiny.
    e0sq=np.where(delta>=0,2*eb*eb/np.maximum(d+delta,1e-100),(d-delta)/2)
    return b2,e2,eb,e0sq


def drift_velocity(e,b):
    b2,_,_,e0sq=invariants(e,b)
    denominator=b2+e0sq
    return np.divide(cross(e,b),denominator[None],out=np.zeros_like(e),where=denominator[None]>1e-100)


def parallel_component(e,b):
    b2=np.sum(b*b,axis=0); eb=np.sum(e*b,axis=0)
    return b*np.divide(eb,b2,out=np.zeros_like(eb),where=b2>1e-100)[None]


def damp_parallel(e,b,kappa,dt):
    return e-parallel_component(e,b)*(-np.expm1(-np.asarray(kappa)*dt))[None]


def energy(e,b,h):
    return float(.5*np.sum(e*e+b*b)*h**3)


def filter4(a,strength):
    """Commuting linear KO-style fourth-difference filter; strength=0 disables."""
    if strength==0: return a
    for axis in (1,2,3):
        a=a-strength/16*(np.roll(a,2,axis)-4*np.roll(a,1,axis)+6*a
                         -4*np.roll(a,-1,axis)+np.roll(a,-2,axis))
    return a


class Magnetosphere:
    def __init__(self,config=PlasmaConfig(),dtype=np.float64):
        self.config=config; self.h=2*config.half_extent/config.n
        self.axis=(np.arange(config.n)+.5)*self.h-config.half_extent
        self.xyz=np.array(np.meshgrid(self.axis,self.axis,self.axis,indexing='ij'),dtype=dtype)
        self.r=np.sqrt(np.sum(self.xyz*self.xyz,axis=0))
        self.cyl=np.sqrt(self.xyz[0]**2+self.xyz[1]**2)
        self.beta=np.array([-self.xyz[1],self.xyz[0],np.zeros_like(self.r)])
        # Conductivity transitions at spherical R_LC as in outer-magnetosphere
        # models; the light cylinder drawn in diagnostics remains cylindrical.
        outer=.5*(1+np.tanh((self.r-1)/config.conductivity_transition))
        self.kappa=config.conductivity_inside*(1-outer)+config.conductivity_outside*outer
        s=np.clip((config.star_radius+self.h*.5-self.r)/self.h,0,1)
        self.star_mask=s*s*(3-2*s) if config.conductor else np.zeros_like(s)
        q=np.clip((self.r-config.sponge_start)/(config.half_extent-config.sponge_start),0,1)
        self.sponge=config.sponge_strength*q*q*q
        self.a=self.star_potential(0).astype(dtype)
        self.e=np.zeros_like(self.a)
        self.time=0.; self.step_number=0
        self.ledger={key:0. for key in ('conductor_exchange','ohmic_loss','inner_projection_loss','sponge_exchange','filter_loss','maxwell_drift_error')}
        self.initial_energy=energy(self.e,self.b,self.h)

    @property
    def b(self): return curl(self.a,self.h)

    def rotation(self,t):
        # Smooth angular velocity 0->1. Integrate analytically for phase.
        T=self.config.spinup_time; u=min(max(t/T,0),1)
        omega=u*u*(3-2*u)
        phase=T*(u**3-.5*u**4)+max(0,t-T)
        return phase,omega

    def star_potential(self,t):
        phase,_=self.rotation(t); chi=np.deg2rad(self.config.inclination_deg)
        moment=np.array([np.sin(chi)*np.cos(phase),np.sin(chi)*np.sin(phase),np.cos(chi)])[:,None,None,None]
        R=self.config.star_radius
        # C1 radial continuation avoids a singularity inside the conductor.
        f=np.where(self.r>=R,self.r**-3,(2.5-1.5*(self.r/R)**2)/R**3)
        return cross(np.broadcast_to(moment,self.xyz.shape),self.xyz)*f[None]

    def rhs(self,a,e):
        b=curl(a,self.h); de=curl(b,self.h)
        if self.config.plasma:
            de-=divergence(e,self.h)[None]*drift_velocity(e,b)
        return -e,de

    def maxwell_step(self,dt):
        a,e=self.a,self.e
        ka1,ke1=self.rhs(a,e)
        ka2,ke2=self.rhs(a+.5*dt*ka1,e+.5*dt*ke1)
        ka3,ke3=self.rhs(a+.5*dt*ka2,e+.5*dt*ke2)
        ka4,ke4=self.rhs(a+dt*ka3,e+dt*ke3)
        self.a=a+dt/6*(ka1+2*ka2+2*ka3+ka4)
        self.e=e+dt/6*(ke1+2*ke2+2*ke3+ke4)

    def step(self,dt=None):
        dt=self.config.courant*self.h if dt is None else dt
        if not 0<dt<=.35*self.h: raise ValueError('Timestep violates configured stability limit')
        b=self.b; old=energy(self.e,b,self.h)
        if self.config.plasma:
            self.e=damp_parallel(self.e,b,self.kappa,dt/2)
        now=energy(self.e,b,self.h); self.ledger['ohmic_loss']+=old-now; old=now
        self.maxwell_step(dt); b=self.b
        now=energy(self.e,b,self.h); self.ledger['maxwell_drift_error']+=now-old; old=now
        if self.config.plasma:
            self.e=damp_parallel(self.e,b,self.kappa,dt/2)
        now=energy(self.e,b,self.h); self.ledger['ohmic_loss']+=old-now; old=now
        if self.config.plasma and self.config.inside_magnetic_limit:
            e2=np.sum(self.e*self.e,axis=0); b2=np.sum(b*b,axis=0)
            factor=np.minimum(1,np.sqrt(b2/np.maximum(e2,1e-100)))
            factor=np.where((self.r<1)&(self.star_mask<.99),factor,1)
            self.e*=factor[None]
        now=energy(self.e,b,self.h); self.ledger['inner_projection_loss']+=old-now; old=now
        self.a=filter4(self.a,self.config.filter_strength)
        self.e=filter4(self.e,self.config.filter_strength); b=self.b
        now=energy(self.e,b,self.h); self.ledger['filter_loss']+=old-now; old=now
        damping=np.exp(-dt*self.sponge)
        self.a*=damping[None]; self.e*=damping[None]; b=self.b
        now=energy(self.e,b,self.h); self.ledger['sponge_exchange']+=now-old; old=now
        self.time+=dt; self.step_number+=1
        if self.config.conductor:
            self.a+=self.star_mask[None]*(self.star_potential(self.time)-self.a)
            b=self.b; _,omega=self.rotation(self.time)
            target=-cross(self.beta*omega,b)
            self.e+=self.star_mask[None]*(target-self.e)
        now=energy(self.e,b,self.h); self.ledger['conductor_exchange']+=now-old
        if not np.isfinite(now): raise FloatingPointError('Nonfinite field energy')
        return now

    def diagnostics(self):
        b=self.b; div=divergence(b,self.h)
        mask=(self.r>max(.65,self.config.star_radius+3*self.h))&(self.r<self.config.sponge_start-.3)
        ep=parallel_component(self.e,b)
        numer=np.linalg.norm(ep[:,mask]); denom=np.linalg.norm(self.e[:,mask])
        J=(divergence(self.e,self.h)[None]*drift_velocity(self.e,b)+self.kappa[None]*ep
           if self.config.plasma else np.zeros_like(self.e))
        final=energy(self.e,b,self.h)
        expected=(self.initial_energy+self.ledger['conductor_exchange']-self.ledger['ohmic_loss']
                  -self.ledger['inner_projection_loss']+self.ledger['sponge_exchange']
                  -self.ledger['filter_loss']+self.ledger['maxwell_drift_error'])
        return dict(time=self.time,steps=self.step_number,energy=final,
                    divB_max=float(np.max(np.abs(div))),
                    divB_relative=float(np.linalg.norm(div)*self.h/max(np.linalg.norm(b),1e-100)),
                    parallel_E_fraction=float(numer/max(denom,1e-100)),
                    max_drift_speed=float(np.max(np.linalg.norm(drift_velocity(self.e,b),axis=0))),
                    energy_ledger_closure=float((final-expected)/max(final,1)),
                    ohmic_power_active_domain=(float(np.sum(self.kappa[mask]*np.sum(ep[:,mask]**2,axis=0))*self.h**3)
                                              if self.config.plasma else 0.),
                    current_norm_active_domain=float(np.linalg.norm(J[:,mask])),ledger=self.ledger.copy())

    def save(self,path,compressed=True):
        import json
        writer=np.savez_compressed if compressed else np.savez
        writer(path,e=self.e,b=self.b,a=self.a,axis=self.axis,
                            time=self.time,config=json.dumps(asdict(self.config)),
                            initial_energy=self.initial_energy,
                            generation_source_sha256=SOLVER_SOURCE_SHA256,
                            diagnostics=json.dumps(self.diagnostics()))
