"""Reduced slot-gap geometry, not a self-consistent plasma solution.

Retarded point-dipole B is defined in the inertial frame. Following Bai &
Spitkovsky (2010), corotation supplies the perpendicular electric field and
v = Omega cross r + f Bhat is normalized to c. Curvature includes rotation.
Coordinates below use c=Omega=R_LC=1; q denotes rotating coordinates.
"""
from dataclasses import dataclass, asdict
import numpy as np


@dataclass(frozen=True)
class Geometry:
    period_s: float = .089
    inclination_deg: float = 75.
    observer_deg: float = 63.6
    seed_radius: float = .08
    annulus_inner: float = .95
    annulus_outer: float = .995
    cylindrical_cutoff: float = .95
    radial_cutoff: float = 1.2
    max_time: float = 8.
    viewing_sigma_deg: float = 1.0
    phase_sigma: float = .002

    def __post_init__(self):
        if not all(np.isfinite(x) for x in asdict(self).values()):
            raise ValueError('Geometry must be finite')
        if not (0 < self.period_s and 0 <= self.inclination_deg <= 90
                and 0 < self.observer_deg < 180
                and 0 < self.seed_radius < self.cylindrical_cutoff < 1
                and self.seed_radius < self.radial_cutoff
                and 0 < self.annulus_inner < self.annulus_outer <= 1
                and self.max_time > 0 and self.viewing_sigma_deg > 0
                and self.phase_sigma >= 0):
            raise ValueError('Geometry outside supported domain')

    @property
    def chi(self):
        return np.deg2rad(self.inclination_deg)

    @property
    def rlc_cm(self):
        return 2.99792458e10*self.period_s/(2*np.pi)


def rotate_z(x, angle):
    x=np.asarray(x); angle=np.asarray(angle)
    c,s=np.cos(angle),np.sin(angle)
    return np.stack((c*x[...,0]-s*x[...,1],s*x[...,0]+c*x[...,1],x[...,2]),axis=-1)


def dipole_field(x, chi, time=0.):
    """Exact retarded magnetic point dipole, unit moment; no stellar boundary."""
    x=np.asarray(x,dtype=float); r=np.linalg.norm(x,axis=-1)
    if np.any(r<=0):
        raise ValueError('Point-dipole origin is excluded')
    n=x/r[...,None]; p=time-r
    m=np.stack((np.sin(chi)*np.cos(p),np.sin(chi)*np.sin(p),np.full_like(p,np.cos(chi))),axis=-1)
    md=np.stack((-m[...,1],m[...,0],np.zeros_like(p)),axis=-1)
    mdd=np.stack((-m[...,0],-m[...,1],np.zeros_like(p)),axis=-1)
    def term(a,k):
        return k*n*np.sum(n*a,axis=-1)[...,None]-a
    return term(m,3)/r[...,None]**3+term(md,3)/r[...,None]**2+term(mdd,1)/r[...,None]


def velocity(q, chi, pole=1):
    """Inertial direction, rotating-coordinate velocity, |v dot Bhat|.

    The same root remains attached to each magnetic hemisphere. This is the
    ultrarelativistic guiding-centre approximation; acceleration does not
    feed back on these prescribed directions.
    """
    b=dipole_field(q,chi); b=b/np.linalg.norm(b,axis=-1)[...,None]
    beta=np.stack((-q[...,1],q[...,0],np.zeros_like(q[...,0])),axis=-1)
    dot=np.sum(beta*b,axis=-1)
    discriminant=dot*dot+1-np.sum(beta*beta,axis=-1)
    if np.any(discriminant < 0):
        raise ValueError('Corotation trajectory exceeded its real-velocity domain')
    projection=np.sqrt(discriminant)
    f=-dot+np.asarray(pole)*projection
    relative=f[...,None]*b
    return beta+relative, relative, projection


def curvature(q,chi,pole=1,h=1e-5):
    v,rel,projection=velocity(q,chi,pole)
    vp=velocity(q+h*rel,chi,pole)[0]
    vm=velocity(q-h*rel,chi,pole)[0]
    zcross=np.stack((-v[...,1],v[...,0],np.zeros_like(v[...,0])),axis=-1)
    dv=(vp-vm)/(2*h)+zcross
    # Remove the round-off component parallel to the unit direction.
    dv=dv-v*np.sum(dv*v,axis=-1)[...,None]
    norm=np.linalg.norm(dv,axis=-1)
    return np.divide(1.,norm,out=np.full_like(norm,np.inf),where=norm>0)


def seed_positions(theta,azimuth,chi,radius):
    theta,azimuth=np.broadcast_arrays(theta,azimuth)
    x=np.sin(theta)*np.cos(azimuth); y=np.sin(theta)*np.sin(azimuth); z=np.cos(theta)
    return radius*np.stack((np.cos(chi)*x+np.sin(chi)*z,y,-np.sin(chi)*x+np.cos(chi)*z),axis=-1)


def _open_lines(theta,azimuth,g,step=.01):
    """Geometric open means crosses cylindrical R_LC (vacuum lines close later).

    For the lower bisection bracket only, reaching spherical 4 R_LC is an
    open proxy. Rim convergence must verify that the final open bracket
    actually crosses the cylinder; no proxy is allowed for a returned rim.
    """
    q=seed_positions(theta,azimuth,g.chi,g.seed_radius)
    alive=np.ones(len(q),bool); opened=np.zeros(len(q),bool); proxy=np.zeros(len(q),bool)
    for _ in range(1800):
        if not alive.any(): break
        ii=np.flatnonzero(alive); x=q[ii]
        def direction(a):
            b=dipole_field(a,g.chi)
            return b/np.linalg.norm(b,axis=-1)[:,None]
        k1=direction(x); k2=direction(x+.5*step*k1)
        k3=direction(x+.5*step*k2); k4=direction(x+step*k3)
        q[ii]=x+step/6*(k1+2*k2+2*k3+k4)
        r=np.linalg.norm(q[ii],axis=-1); cyl=np.linalg.norm(q[ii,:2],axis=-1)
        cross=cyl>=1; far=r>=4
        opened[ii[cross|far]]=True; proxy[ii[far&~cross]]=True
        alive[ii[cross|far|(r<g.seed_radius*.99)]]=False
    if alive.any():
        raise RuntimeError('Unclassified cap traces; increase trace horizon')
    return opened,proxy


def polar_cap_rim(g,azimuth,trace_step=.01,iterations=18):
    azimuth=np.asarray(azimuth)
    lo=np.full_like(azimuth,.015); hi=np.full_like(azimuth,1.)
    if not np.all(_open_lines(lo,azimuth,g,trace_step)[0]):
        raise RuntimeError('Cap lower bracket not open')
    if np.any(_open_lines(hi,azimuth,g,trace_step)[0]):
        raise RuntimeError('Cap upper bracket not closed')
    for _ in range(iterations):
        mid=(lo+hi)/2; opened,_=_open_lines(mid,azimuth,g,trace_step)
        lo=np.where(opened,mid,lo); hi=np.where(opened,hi,mid)
    opened,proxy=_open_lines(lo,azimuth,g,trace_step)
    if not np.all(opened) or np.any(proxy):
        raise RuntimeError('Rim not established by actual cylinder crossings')
    return (lo+hi)/2


def trace_bundle(g,n_azimuth=128,n_rings=3,dt=.006,trace_step=.01):
    """Steady injection from both polar annuli, integrated in lab time.

    Each seed weight is proportional to unsigned magnetic flux at the outer
    injection surface. No fitted coefficient is attached to arrival phase.
    """
    if n_azimuth<8 or n_rings<1 or dt<=0: raise ValueError('Invalid quadrature')
    az=(np.arange(n_azimuth)+.5)*2*np.pi/n_azimuth
    rim=polar_cap_rim(g,az,trace_step)
    u=g.annulus_inner+(np.arange(n_rings)+.5)/n_rings*(g.annulus_outer-g.annulus_inner)
    th=(u[:,None]*rim[None,:]).ravel(); ps=np.broadcast_to(az,(n_rings,n_azimuth)).ravel()
    north=seed_positions(th,ps,g.chi,g.seed_radius)
    b=dipole_field(north,g.chi); er=north/g.seed_radius
    dtheta=np.tile(rim,n_rings)*(g.annulus_outer-g.annulus_inner)/n_rings
    weight=np.abs(np.sum(b*er,axis=-1))*g.seed_radius**2*np.sin(th)*dtheta*(2*np.pi/n_azimuth)
    weight=np.r_[weight,weight]; weight/=weight.sum()
    q=np.concatenate((north,-north)); pole=np.r_[np.ones(len(north)),-np.ones(len(north))]
    alive=np.ones(len(q),bool); chunks=[]
    for _ in range(int(np.ceil(g.max_time/dt))):
        if not alive.any(): break
        ii=np.flatnonzero(alive); x=q[ii]; p=pole[ii]
        def step(a,sign,h):
            h=np.asarray(h)[...,None]
            k1=velocity(a,g.chi,sign)[1]
            k2=velocity(a+.5*h*k1,g.chi,sign)[1]
            k3=velocity(a+.5*h*k2,g.chi,sign)[1]
            k4=velocity(a+h*k3,g.chi,sign)[1]
            return a+h/6*(k1+2*k2+2*k3+k4)
        def outside(a):
            r=np.linalg.norm(a,axis=-1)
            return ((np.linalg.norm(a[:,:2],axis=-1)>g.cylindrical_cutoff)
                    |(r>g.radial_cutoff)|(r<g.seed_radius*.99))
        # Clip the final residence time to the first escape surface, rather
        # than radiating a complete step beyond the declared emitting region.
        elapsed=np.full(len(x),dt); nxt=step(x,p,elapsed)
        escaped=outside(nxt)
        if escaped.any():
            low=np.zeros(escaped.sum()); high=np.full(escaped.sum(),dt)
            a=x[escaped]; sign=p[escaped]
            for _ in range(18):
                h=(low+high)/2
                out=outside(step(a,sign,h))
                high=np.where(out,h,high); low=np.where(out,low,h)
            elapsed[escaped]=low
            nxt[escaped]=step(a,sign,low)
        mid=step(x,p,.5*elapsed)
        v,_,proj=velocity(mid,g.chi,p); rho=curvature(mid,g.chi,p)
        # Rectangular padded arrays allow independent transport of every seed.
        pos=q.copy(); pos[ii]=mid
        direction=np.zeros_like(q); direction[ii]=v
        radcurv=np.full(len(q),np.inf); radcurv[ii]=rho
        duration=np.zeros(len(q)); duration[ii]=elapsed
        project=np.zeros(len(q)); project[ii]=proj
        chunks.append((pos,direction,radcurv,duration,project))
        q[ii]=nxt
        alive[ii[escaped]]=False
    if alive.any():
        raise RuntimeError('Particle trace horizon reached before an escape cutoff')
    pos,direction,rho,duration,projection=(np.stack([c[j] for c in chunks]) for j in range(5))
    arrival=(-np.arctan2(direction[...,1],direction[...,0])-np.sum(direction*pos,axis=-1))/(2*np.pi)%1
    return dict(position=pos,direction=direction,rho=rho,dt=duration,
                projection=projection,phase=arrival,seed_weight=weight,pole=pole,
                cap_azimuth=az,cap_rim=rim,geometry=asdict(g),step=dt,
                n_azimuth=n_azimuth,n_rings=n_rings,trace_step=trace_step)


def sky_profiles(bundle,photon_yield,g,nphase=512,view_sigma_deg=None):
    """Lab photon yield, CIC phase deposition and explicit sky-bin averaging.

    The Gaussian in mu is a numerical viewing bin, not a 1/gamma beam.
    Positive phase_sigma introduces a population spread assumption; zero
    leaves only phase deposition. The input may be raw photon yield or
    detector-area-weighted yield; this routine adds no detector response.
    """
    from scipy.ndimage import gaussian_filter1d
    sigdeg=g.viewing_sigma_deg if view_sigma_deg is None else view_sigma_deg
    zeta=np.deg2rad(g.observer_deg); sigma=np.sin(zeta)*np.deg2rad(sigdeg)
    mu=bundle['direction'][...,2]
    angular=np.exp(-.5*((mu-np.cos(zeta))/sigma)**2)/(np.sqrt(2*np.pi)*sigma)
    weight=angular*bundle['seed_weight'][None,:]*(bundle['dt']>0)
    p=bundle['phase']*nphase; left=np.floor(p).astype(int)%nphase; frac=p-np.floor(p)
    profiles=[]
    for j in range(photon_yield.shape[-1]):
        y=photon_yield[...,j]*weight
        a=np.bincount(left.ravel(),weights=(y*(1-frac)).ravel(),minlength=nphase)
        a+=np.bincount(((left+1)%nphase).ravel(),weights=(y*frac).ravel(),minlength=nphase)
        if g.phase_sigma>0:
            a=gaussian_filter1d(a,g.phase_sigma*nphase,mode='wrap')
        profiles.append(a/np.mean(a) if np.mean(a)>0 else a)
    return np.arange(nphase)/nphase,np.array(profiles)
