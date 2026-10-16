#!/usr/bin/env python3
"""Independent normalization/flux calibration and curl-conductor diagnostics.

No native fields are changed. The inner fit uses the existing actual D4-curl
source basis, calibrated here against a separately differentiated analytic
finite dipole. Surface diagnostics use angular quadrature plus interpolation;
the compatible cube flux separately checks the native D4 divergence theorem.
"""
from pathlib import Path
import argparse,hashlib,json
import numpy as np
from run_coffee_assessment import load,derivative,sample,rotation
from run_coffee_curl_audit import source_basis


def finite_radial_factor(r,radius):
    r=np.asarray(r,float);u=r*r/radius**2-1
    coeff=np.array([1.,-1.5,1.875,-2.1875,2.4609375])
    poly=np.polynomial.polynomial.polyval(u,coeff)/radius**3
    derivative_over_r=2*np.polynomial.polynomial.polyval(u,np.arange(1,5)*coeff[1:])/radius**5
    outside=r>=radius;safe=np.where(outside,r,radius)
    return np.where(outside,safe**-3,poly),np.where(outside,-3*safe**-5,derivative_over_r)


def analytic_finite_dipole(points,moment,radius=.5):
    """Closed-form curl of A=f(r) m cross r, independent of finite differences."""
    p=np.asarray(points,float);m=np.asarray(moment,float);r=np.linalg.norm(p,axis=-1)
    f,g=finite_radial_factor(r,radius)
    return (2*f+g*r*r)[...,None]*m-g[...,None]*np.sum(p*m,axis=-1)[...,None]*p


def surface_diagnostics(provider,radius,source_radius=.5,nu=64,nphi=128):
    """Signed flux, one-half unsigned flux, radial dipole harmonic and power."""
    u,w=np.polynomial.legendre.leggauss(nu);phi=np.arange(nphi)*2*np.pi/nphi
    u,phi=np.meshgrid(u,phi,indexing='ij');n=np.stack((np.sqrt(1-u*u)*np.cos(phi),np.sqrt(1-u*u)*np.sin(phi),u),axis=-1)
    weights=w[:,None]*(2*np.pi/nphi);e,b=provider(radius*n);br=np.sum(b*n,axis=-1)
    signed=float(radius**2*np.sum(br*weights));unsigned=float(radius**2*np.sum(abs(br)*weights))
    factor=float(finite_radial_factor(radius,source_radius)[0])
    moment=3/(8*np.pi*factor)*np.sum(br[...,None]*n*weights[...,None],axis=(0,1))
    power=float(radius**2*np.sum(np.sum(np.cross(e,b)*n,axis=-1)*weights)/(4*np.pi))
    return {'radius_rlc':radius,'signed_magnetic_flux':signed,'half_unsigned_magnetic_flux':unsigned/2,
            'net_over_unsigned_flux':signed/unsigned if unsigned else 0.,
            'radial_dipole_harmonic_moment':moment.tolist(),'luminosity_mu2_omega4_over_c3':power,
            'angular_quadrature':[nu,nphi]}


def compatible_cube_flux(b,axes,half_extent):
    """Exactly compatible closed-cube flux for the native centered D4 stencil.

    Face reconstruction F_(i+1/2)=(-B_(i+2)+7B_(i+1)+7B_i-B_(i-1))/12
    gives D4 B=(F_(i+1/2)-F_(i-1/2))/h. This is a discrete flux, distinct from
    interpolating B onto a physical spherical surface.
    """
    h=axes[0][1]-axes[0][0]
    limits=[(int(np.searchsorted(a,-half_extent)),int(np.searchsorted(a,half_extent,side='right')-1)) for a in axes]
    if any(lo<2 or hi>len(a)-3 or lo>hi for a,(lo,hi) in zip(axes,limits)):raise ValueError('Cube requires a two-node interior stencil buffer')
    signed=0.;unsigned=0.
    for j,(lo,hi) in enumerate(limits):
        def plane(i):
            index=[slice(l,r+1) for l,r in limits];index[j]=i
            return b[j][tuple(index)]
        high=(-plane(hi+2)+7*plane(hi+1)+7*plane(hi)-plane(hi-1))/12
        low=(-plane(lo+1)+7*plane(lo)+7*plane(lo-1)-plane(lo-2))/12
        signed+=float(np.sum(high-low)*h*h);unsigned+=float(np.sum(abs(high)+abs(low))*h*h)
    div=sum(derivative(b[j],h,j) for j in range(3));index=tuple(slice(lo,hi+1) for lo,hi in limits)
    volume=float(np.sum(div[index])*h**3)
    return {'requested_half_extent_rlc':half_extent,'node_index_limits':limits,
            'signed_discrete_flux':signed,'volume_integral_D4_divB':volume,
            'closure_error':signed-volume,'net_over_unsigned_discrete_flux':signed/unsigned if unsigned else 0.}


def fit_source(d,rmax=.35):
    xyz=np.stack(np.meshgrid(*d['axes'],indexing='ij'));r=np.linalg.norm(xyz,axis=0);mask=r<rmax
    h=d['axes'][0][1]-d['axes'][0][0];basis=source_basis(xyz[:,mask].T,h,d['metadata']['star_radius'])
    values=d['b'][:,mask].T.ravel();m=np.linalg.lstsq(basis,values,rcond=None)[0]
    cfg=d['metadata']['config'];initial=np.array([cfg['p1'],cfg['p2'],cfg['p3']],float);initial/=np.linalg.norm(initial)
    expected=rotation(initial,d['metadata'].get('magnetic_phase_rad',d['time']))
    phase=(np.angle((m[0]+1j*m[1])/(expected[0]+1j*expected[1])) if np.hypot(*expected[:2])>1e-12 and np.hypot(*m[:2])>1e-12 else None)
    beta=np.array([-xyz[1],xyz[0],np.zeros_like(r)]);target=-np.cross(beta,d['b'],axisa=0,axisb=0,axisc=0)
    den=np.linalg.norm(target[:,mask]);res=np.linalg.norm((d['e']-target)[:,mask])
    return {'fit_radius_rlc':rmax,'native_nodes':int(mask.sum()),'moment':m.tolist(),'magnitude':float(np.linalg.norm(m)),
            'relative_vector_error_to_prescribed':float(np.linalg.norm(m-expected)),
            'phase_error_degrees':float(np.degrees(phase)) if phase is not None else None,
            'inclination_degrees':float(np.degrees(np.arctan2(np.hypot(*m[:2]),m[2]))),
            'relative_fit_residual':float(np.linalg.norm(basis@m-values)/np.linalg.norm(values)),
            'basis_condition_number':float(np.linalg.cond(basis)),
            'conductor_E_relative_L2':float(res/den) if den else None}


def calibration():
    m=np.array([.3,-.4,np.sqrt(.75)]);rows=[]
    rng=np.random.default_rng(38231);p=rng.normal(size=(240,3));p*=rng.uniform(.03,.35,len(p))[:,None]/np.linalg.norm(p,axis=-1)[:,None]
    exact=analytic_finite_dipole(p,m)
    for n in (65,97,129,193,257):
        h=6/n;basis=source_basis(p,h,.5);actual=(basis@m).reshape(-1,3)
        fitted=np.linalg.lstsq(basis,exact.ravel(),rcond=None)[0]
        rows.append({'grid_n':n,'h':h,'D4_field_relative_error_to_analytic':float(np.linalg.norm(actual-exact)/np.linalg.norm(exact)),
                     'fit_analytic_field_moment_error':float(np.linalg.norm(fitted-m))})
    surfaces=[]
    for radius in (.25,.35,.5,.75,1.,2.):
        provider=lambda p:(np.zeros_like(p),analytic_finite_dipole(p,m))
        row=surface_diagnostics(provider,radius);row['moment_error']=float(np.linalg.norm(np.array(row['radial_dipole_harmonic_moment'])-m))
        if radius>=.5:row['half_unsigned_flux_relative_error']=row['half_unsigned_magnetic_flux']/(2*np.pi/radius)-1
        surfaces.append(row)
    return {'known_unit_moment':m.tolist(),'analytic_vs_D4':rows,'analytic_surface_checks':surfaces}


def audit_run(directory):
    rows=[]
    for path in sorted(directory.glob('snapshot-*.npz')):
        d=load(path);xyz=np.stack(np.meshgrid(*d['axes'],indexing='ij'));r=np.linalg.norm(xyz,axis=0);h=d['axes'][0][1]-d['axes'][0][0]
        div=sum(derivative(d['b'][j],h,j) for j in range(3));d2=sum(np.gradient(d['b'][j],h,axis=j,edge_order=2) for j in range(3));bm=np.linalg.norm(d['b'],axis=0)
        shells={}
        for lo,hi in ((0,.35),(.35,.75),(.75,1),(1,1.5),(1.5,2),(.75,2),(2,2.5)):
            mask=(r>lo)&(r<hi);den=np.sum(bm[mask]**2)
            shells[f'{lo:g}<r<{hi:g}']={'D4_div_rms_over_B_RLC':float(np.sqrt(np.sum(div[mask]**2)/den)),
                                     'D2_div_rms_over_B_RLC':float(np.sqrt(np.sum(d2[mask]**2)/den))}
        cfg=d['metadata']['config'];omega=cfg['omega'];lower=omega*np.array(cfg['lower']);size=omega*np.array(cfg['size']);width=np.array(cfg['pml'])*h
        points=xyz.transpose(1,2,3,0);sponge=np.any((points<lower+width)|(points>lower+size-width),axis=-1)
        interior=np.zeros(r.shape,bool);interior[2:-2,2:-2,2:-2]=True;exterior=interior&(r>.75);sponge&=exterior
        provider=lambda p:(sample(d,'e',p),sample(d,'b',p))
        rows.append({'snapshot':str(path.resolve()),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'time':d['time'],
            'source_fits':[fit_source(d,x) for x in (.25,.35,.45)],'divergence_shells':shells,
            'sponge':{'physical_widths_rlc':width.tolist(),'node_count':int(sponge.sum()),
                'D4_div_rms_over_B_RLC':float(np.sqrt(np.sum(div[sponge]**2)/np.sum(bm[sponge]**2))),
                'fraction_exterior_D4_div_squared_in_sponge':float(np.sum(div[sponge]**2)/np.sum(div[exterior]**2))},
            'spheres':[surface_diagnostics(provider,x,d['metadata']['star_radius']) for x in (.35,.5,.75,1.,1.5,2.)],
            'compatible_cube_fluxes':[compatible_cube_flux(d['b'],d['axes'],x) for x in (.35,.75,1.5,2.)]})
    if not rows:raise ValueError('No canonical snapshots available')
    first=rows[0]['source_fits'][1];last=rows[-1]['source_fits'][1]
    return {'run_dir':str(directory.resolve()),'metadata':d['metadata'],'rows':rows,
            'inner_moment_fractional_change':last['magnitude']/first['magnitude']-1,
            'final_inner_phase_error_degrees':last['phase_error_degrees'],
            'source_fit_scope':'Actual finite C4 source D4-curl basis; fitted moment is an inner-field diagnostic, not a global dipole fit',
            'surface_scope':'Linear interpolation and Gauss/azimuth quadrature; effective l=1 radial harmonic can include exterior currents',
            'D2_scope':'Different continuum derivative diagnostic only; D4 matches native solver and is the constraint test',
            'emission_release_allowed':False,'vela_validation':False}


def compare_final_fields(directory_a,directory_b):
    """Compare snapshots on first-run nodes; a directory selects its final file.

    Different times are permitted only at the same rotational phase, allowing
    a direct integer-period repeatability measurement without rotating fields.
    The second run is no exact solution.
    """
    paths=[p if p.is_file() else sorted(p.glob('snapshot-*.npz'))[-1] for p in (directory_a,directory_b)]
    a,b=[load(p) for p in paths]
    if abs(np.angle(np.exp(1j*(a['time']-b['time']))))>1e-12:raise ValueError('Rotational phases must match')
    for key in ('coordinate_units','field_units','star_radius','inclination_deg'):
        if a['metadata'][key]!=b['metadata'][key]:raise ValueError(f'Mismatched {key}')
    xyz=np.stack(np.meshgrid(*a['axes'],indexing='ij'));r=np.linalg.norm(xyz,axis=0);rows={}
    for lo,hi in ((0,.35),(.35,.75),(.75,1),(1,1.5),(1.5,2),(.75,2)):
        mask=(r>lo)&(r<hi);points=xyz[:,mask].T;ae=a['e'][:,mask].T;ab=a['b'][:,mask].T
        be=sample(b,'e',points);bb=sample(b,'b',points)
        angles=np.degrees(np.arccos(np.clip(np.sum(ab*bb,axis=-1)/(np.linalg.norm(ab,axis=-1)*np.linalg.norm(bb,axis=-1)),-1,1)))
        rows[f'{lo:g}<r<{hi:g}']={'cells':int(mask.sum()),
            'relative_L2_E':float(np.linalg.norm(ae-be)/np.linalg.norm(be)),
            'relative_L2_B':float(np.linalg.norm(ab-bb)/np.linalg.norm(bb)),
            'median_B_angle_degrees':float(np.median(angles)),'p90_B_angle_degrees':float(np.quantile(angles,.9))}
    return {'first_snapshot':str(paths[0].resolve()),'second_snapshot':str(paths[1].resolve()),
            'snapshot_sha256':[hashlib.sha256(p.read_bytes()).hexdigest() for p in paths],
            'times':[a['time'],b['time']],'grid_n':[len(d['axes'][0]) for d in (a,b)],
            'boundary_models':[d['metadata']['boundary_model'] for d in (a,b)],'shells':rows,
            'scope':'Second field linearly interpolated to first nodes; differences divided by second L2 norm; neither field is a validated reference.'}


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('directories',type=Path,nargs='*')
    ap.add_argument('--compare',type=Path,nargs=2,action='append',default=[],metavar=('FIRST','SECOND'))
    ap.add_argument('--output',type=Path,default=Path('outputs/curl-source-audit'));args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    result={'calibration':calibration(),'runs':[audit_run(p) for p in args.directories],
            'field_comparisons':[compare_final_fields(*pair) for pair in args.compare],
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'emission_release_allowed':False,'vela_validation':False}
    (args.output/'audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'output':str(args.output),'runs':[{'path':r['run_dir'],'moment_change':r['inner_moment_fractional_change'],'phase_error_degrees':r['final_inner_phase_error_degrees']} for r in result['runs']]},indent=2))


if __name__=='__main__':main()
