"""Finite spatial support for an already calculated directional-emission node.

Fits the eight r/u²/phi-ordered corners by a uniform affine parallelepiped and
integrates observer ray path length with deterministic 2x2 pixel quadrature.
This changes spatial quadrature, never the supplied angular radiation pattern.
Cell power is frozen at its node. Without source/arrival windows, spatial
retarded-time differences are not applied. With both windows, the source-time
overlap is integrated through each static subray. The caller must subdivide
moving geometry and audit unresolved/degenerate cubature separately.

Discrete normalization makes energy+croploss conservative, but is not evidence
of converged geometry. Arbitrary spatial subdivision need not commute exactly
with a fixed 2x2 raster: each cell receives its own quadrature normalization.
The pre-normalization volume/area errors expose that effect. The node's angular
power is not recalculated at interior points, and that error is not bounded here.
"""
from dataclasses import dataclass
from functools import lru_cache
import numpy as np

C_M_S = 299792458.0
SIGNS = np.array([[r,u,p] for r in (-1.,1.) for u in (-1.,1.) for p in (-1.,1.)])
_LOSS_NODES_1D,_LOSS_WEIGHTS_1D=np.polynomial.legendre.leggauss(4)
_LOSS_NODES=np.stack(np.meshgrid(*([_LOSS_NODES_1D]*3),indexing='ij'),axis=-1).reshape(-1,3)
_LOSS_WEIGHTS=np.prod(np.stack(np.meshgrid(*([_LOSS_WEIGHTS_1D/2]*3),indexing='ij'),axis=-1),axis=-1).ravel()
_GAUSS_FACTOR=1+SIGNS[:,None,:]*SIGNS[None,:,:]/np.sqrt(3)
_GAUSS_DERIVATIVE=np.empty((8,8,3))
for _j in range(3):
    _other=[k for k in range(3) if k!=_j]
    _GAUSS_DERIVATIVE[:,:,_j]=SIGNS[None,:,_j]*_GAUSS_FACTOR[:,:,_other[0]]*_GAUSS_FACTOR[:,:,_other[1]]/8


@dataclass
class CellPatch:
    y_slice: slice
    x_slice: slice
    patch_xyz: np.ndarray
    lost_xyz: np.ndarray
    metadata: dict
    window_xyz_energy: np.ndarray

    @property
    def image_slice(self):
        return self.y_slice, self.x_slice

    def deposit_into(self, image):
        if image.ndim != 3 or image.shape[-1] != 3:
            raise ValueError('Accumulation image must have shape (height,width,3)')
        image[self.image_slice] += self.patch_xyz
        return self.lost_xyz


@lru_cache(maxsize=32)
def _basis(observer_tuple, roll_deg):
    n=np.asarray(observer_tuple,float)
    norm=np.linalg.norm(n)
    if n.shape!=(3,) or not np.isfinite(n).all() or norm==0:
        raise ValueError('Observer must be one finite nonzero 3-vector')
    n=n/norm
    hint=np.array([0.,0.,1.])
    ey=hint-np.dot(hint,n)*n
    if np.linalg.norm(ey)<1e-10:
        hint=np.eye(3)[np.argmin(np.abs(n))]
        ey=hint-np.dot(hint,n)*n
    ey/=np.linalg.norm(ey)
    ex=np.cross(ey,n)
    angle=np.deg2rad(roll_deg)
    return np.cos(angle)*ex-np.sin(angle)*ey, np.sin(angle)*ex+np.cos(angle)*ey, n


def fit_affine_support(vertices_m):
    """Return centre, half-edge columns and explicit corner-fit diagnostics.

    The trilinear-hexahedron volume uses 2x2x2 Gauss integration. For an
    orientation-preserving trilinear cell its signed determinant is integrated
    exactly. A folded cell is flagged and must not be called an exact volume.
    """
    v=np.asarray(vertices_m,float)
    if v.shape!=(8,3) or not np.isfinite(v).all():
        raise ValueError('Expected eight finite r/u²/phi-ordered 3D corners')
    center=v.mean(axis=0)
    H=(SIGNS.T@(v-center)/8).T
    fitted=center+SIGNS@H.T
    residual=np.linalg.norm(v-fitted,axis=1)
    diameter=float(np.max(np.linalg.norm(v[:,None]-v[None,:],axis=-1)))
    det=np.linalg.det(H)
    jac_dets=np.linalg.det(np.einsum('vi,qvj->qij',v-center,_GAUSS_DERIVATIVE))
    trilinear_volume=float(np.sum(np.abs(jac_dets)))
    affine_volume=float(8*abs(det))
    folded=bool(np.any(jac_dets>0) and np.any(jac_dets<0))
    metadata={
      'support_model':'Uniform affine least-squares parallelepiped; corner order radial,u²,phi nested',
      'affine_corner_max_error_m':float(residual.max()),
      'affine_corner_rms_error_m':float(np.sqrt(np.mean(residual**2))),
      'affine_corner_max_error_over_diameter':float(residual.max()/diameter) if diameter else 0.,
      'affine_volume_m3':affine_volume,
      'trilinear_gauss_volume_m3':trilinear_volume,
      'trilinear_gauss_jacobian_sign_changes':folded,
      'interior_fold_certification':'Gauss-node signs checked only; absence of interior folds is not certified',
      'trilinear_volume_exact_if_orientation_preserving':not folded,
      'affine_volume_relative_error_vs_trilinear':float(affine_volume/trilinear_volume-1) if trilinear_volume else None,
    }
    return center,H,fitted,metadata


@dataclass
class AffineSupport:
    """Reference geometry fitted once; treat its arrays as immutable."""
    vertices_m: np.ndarray
    center_m: np.ndarray
    half_edges_m: np.ndarray
    fitted_vertices_m: np.ndarray
    metadata: dict
    svd_u: np.ndarray
    singular_values: np.ndarray
    svd_vh: np.ndarray
    inverse_half_edges: object


def prepare_affine_support(vertices_m):
    center,H,fitted,metadata=fit_affine_support(vertices_m)
    U,s,Vh=np.linalg.svd(H)
    inverse=np.linalg.inv(H) if s.min()>max(s.max(),1e-300)*1e-12 else None
    return AffineSupport(np.asarray(vertices_m,float).copy(),center,H,fitted,metadata,U,s,Vh,inverse)


def rasterize_prepared_cell(support,xyz_energy,observer_direction,*,rotation_z_rad=0.,**kwargs):
    """Reuse one reference fit through the caller's rigid spin-axis rotation."""
    return rasterize_affine_cell(None,xyz_energy,observer_direction,
                                prepared_support=support,rotation_z_rad=rotation_z_rad,**kwargs)


def _empty(energy,metadata):
    metadata.setdefault('window_xyz_energy',energy.tolist())
    metadata.setdefault('arrival_fraction',1.)
    return CellPatch(slice(0,0),slice(0,0),np.zeros((0,0,3)),energy.copy(),metadata,energy.copy())


def _from_weights(weights,x0,y0,energy,width,height,metadata,normalization_weights=None):
    # A time window only selects emitted energy. It must NEVER become the
    # normalization, which remains the full unwindowed geometric support.
    normal=weights if normalization_weights is None else normalization_weights
    total=float(normal.sum())
    window_total=float(weights.sum())
    if not np.isfinite(total) or total<=0 or not np.isfinite(window_total) or window_total<0:
        raise ValueError('Spatial support received invalid quadrature mass')
    arrival_fraction=window_total/total
    window_energy=energy*arrival_fraction
    metadata.update(arrival_fraction=arrival_fraction,window_xyz_energy=window_energy.tolist())
    h,w=weights.shape
    xa,xb=max(0,x0),min(width,x0+w)
    ya,yb=max(0,y0),min(height,y0+h)
    if xa>=xb or ya>=yb:
        metadata['retained_fraction']=0. if window_total else None
        return _empty(window_energy,metadata)
    inside=weights[ya-y0:yb-y0,xa-x0:xb-x0]
    fraction=inside/total
    patch=fraction[...,None]*energy
    lost=np.maximum(0.,window_energy-patch.sum(axis=(0,1)))
    metadata['retained_fraction']=float(inside.sum()/window_total) if window_total else None
    metadata['normalization']='Full UNWINDOWED support discrete quadrature sum, before cropping; window_xyz_energy = in-frame energy + lost_xyz'
    return CellPatch(slice(ya,yb),slice(xa,xb),patch,lost,metadata,window_energy)


def _time_window(source_interval_s,arrival_interval_s,metadata):
    if source_interval_s is None:
        if arrival_interval_s is not None:
            raise ValueError('Arrival window requires a source-time interval')
        return None
    source=np.asarray(source_interval_s,float)
    if source.shape!=(2,) or not np.isfinite(source).all() or source[1]<=source[0]:
        raise ValueError('Source interval must be two increasing finite times')
    metadata['source_interval_s']=source.tolist()
    if arrival_interval_s is None:return None
    arrival=np.asarray(arrival_interval_s,float)
    if arrival.shape!=(2,) or not np.isfinite(arrival).all() or arrival[1]<=arrival[0]:
        raise ValueError('Arrival interval must be two increasing finite times')
    metadata.update(arrival_interval_s=arrival.tolist(),spatial_retardation_applied=True,
                    temporal_geometry='Static affine support during this source subinterval; caller subdivides motion')
    return float(source[0]),float(source[1]-source[0]),float(arrival[0]-source[0]),float(arrival[1]-source[0])


def _windowed_ray_lengths(enter,leave,window,center_delay_s):
    """Exactly integrate source/arrival overlap through each static ray slab.

    The overlap is a piecewise-linear combination of four ReLUs. Integrating
    trapezoids on its four clipped breakpoints equals the ReLU-square primitive
    but avoids subtracting nearly equal large squares. Times are shifted by t0.
    """
    _,duration,a0,a1=window
    a0+=center_delay_s;a1+=center_delay_s
    valid=leave>enter
    lo=np.where(valid,enter,0.);hi=np.where(valid,leave,0.)
    breaks=C_M_S*np.sort(np.array([-a0,duration-a0,-a1,duration-a1]))
    inner=np.clip(breaks,lo[...,None],hi[...,None])
    points=np.concatenate((lo[...,None],inner,hi[...,None]),axis=-1)
    delay=points/C_M_S
    overlap=np.maximum(0.,np.minimum(duration,a1+delay)-np.maximum(0.,a0+delay))
    integral=np.sum(.5*(overlap[...,:-1]+overlap[...,1:])*np.diff(points,axis=-1),axis=-1)/duration
    return np.where(valid,np.minimum(np.maximum(integral,0.),hi-lo),0.)


def _unresolved_cic(H,center_pixel,X,Y,pixel,energy,width,height,metadata,max_points=1000000,window=None,center=None,n=None):
    """Disclosed subpixel/degenerate limit; finite support is still sampled.

    Midpoint tensor nodes have projected edge spacing <=0.5 pixel. CIC is a
    detector reconstruction footprint, not a physical angular-emission width.
    """
    projected=np.vstack((X@H,-Y@H))/pixel
    counts=np.maximum(1,np.ceil(4*np.linalg.norm(projected,axis=0)).astype(int))
    if window is not None:
        _,duration,a0,a1=window
        temporal_scale=min(duration,a1-a0)/8
        counts=np.maximum(counts,np.ceil(2*np.abs(n@H)/(C_M_S*temporal_scale)).astype(int))
        counts=np.maximum(1,counts)
    count=int(np.prod(counts,dtype=np.int64))
    if count>max_points:
        raise ValueError('Degenerate support needs too many midpoint nodes; refine its spatial representation')
    axes=[(np.arange(n)+.5)*2/n-1 for n in counts]
    node=np.stack(np.meshgrid(*axes,indexing='ij'),axis=-1).reshape(-1,3)
    uv=center_pixel+node@projected.T
    base=np.floor(uv).astype(int)
    fraction=uv-base
    x0,y0=base.min(axis=0);x1,y1=base.max(axis=0)+2
    weights=np.zeros((y1-y0,x1-x0));full_weights=np.zeros_like(weights)
    if window is None:
        time_fraction=np.ones(count)
    else:
        _,duration,a0,a1=window
        delay=(center+node@H.T)@n/C_M_S
        time_fraction=np.maximum(0.,np.minimum(duration,a1+delay)-np.maximum(0.,a0+delay))/duration
    for dx,dy in ((0,0),(1,0),(0,1),(1,1)):
        w=(fraction[:,0] if dx else 1-fraction[:,0])*(fraction[:,1] if dy else 1-fraction[:,1])
        ij=(base[:,1]+dy-y0,base[:,0]+dx-x0)
        np.add.at(full_weights,ij,w/count)
        np.add.at(weights,ij,w*time_fraction/count)
    metadata.update(method='Unresolved/degenerate affine support: midpoint tensor nodes plus CIC detector reconstruction',
                    midpoint_nodes=count,projected_node_spacing_bound_pixels=.5,
                    pre_normalization_volume_relative_error=None,
                    pre_normalization_projected_area_relative_error=None,
                    unresolved_support_fallback=True,
                    spatial_retardation_exact_for_static_affine_ray=False,
                    spatial_retardation_method='Midpoint cubature of unresolved/degenerate support; exact overlap at each node, not exact continuous volume',
                    fallback_temporal_delay_step_bound_s=float(np.sum(2*np.abs(n@H)/counts)/C_M_S) if window is not None else None)
    return _from_weights(weights,int(x0),int(y0),energy,width,height,metadata,full_weights)


def rasterize_affine_cell(vertices_m, xyz_energy, observer_direction, *,
                          meters_per_pixel, width, height, roll_deg=-55.,
                          center_fraction=(.5,.53), emission_node_m=None,
                          max_bbox_pixels=2000000, prepared_support=None, rotation_z_rad=0.,
                          source_interval_s=None,arrival_interval_s=None):
    """Return compact CellPatch; no full-frame allocation occurs per cell.

    Energy must already include the caller's physical number/source-time weight.
    Outputs retain that same energy unit per pixel, not radiance per square metre.
    The caller owns motion/time integration, final pixel-area/exposure division,
    observing-band emission and its angular variation inside the finite cell.
    Optional source_interval_s=(t0,t1), arrival_interval_s=(a0,a1) select energy
    using tau=t-n.r/c. xyz_energy must be the FULL source-subinterval energy;
    normalization always uses full unwindowed geometry. lost_xyz describes only
    cropped energy in this arrival window. Sum window_xyz_energy over a covering
    set of arrival bins/images, including periodic wraps, to recover the input.
    """
    energy=np.asarray(xyz_energy,float)
    if energy.shape!=(3,) or not np.isfinite(energy).all() or np.any(energy<0):
        raise ValueError('XYZ energy must be three finite nonnegative numbers')
    if not isinstance(width,(int,np.integer)) or not isinstance(height,(int,np.integer)) or min(width,height)<1:
        raise ValueError('Width and height must be positive integers')
    pixel=float(meters_per_pixel)
    if not np.isfinite(pixel) or pixel<=0 or not np.isfinite(roll_deg):
        raise ValueError('Positive pixel scale and finite roll are required')
    origin=np.asarray(center_fraction,float)*[width,height]
    if origin.shape!=(2,) or not np.isfinite(origin).all():
        raise ValueError('Expected two finite image-centre fractions')
    X,Y,n=_basis(tuple(np.asarray(observer_direction,float)),float(roll_deg))
    rotation=None
    if prepared_support is None:
        if rotation_z_rad!=0:
            raise ValueError('Use prepared geometry to request a reference z rotation')
        center,H,fitted,metadata=fit_affine_support(vertices_m)
        v=np.asarray(vertices_m,float)
        svd=None
    else:
        if not isinstance(prepared_support,AffineSupport) or not np.isfinite(rotation_z_rad):
            raise ValueError('Expected prepared AffineSupport and finite z rotation')
        q=prepared_support
        co,si=np.cos(rotation_z_rad),np.sin(rotation_z_rad)
        rotation=np.array([[co,-si,0.],[si,co,0.],[0.,0.,1.]])
        center=rotation@q.center_m;H=rotation@q.half_edges_m
        fitted=q.fitted_vertices_m@rotation.T;v=q.vertices_m@rotation.T
        metadata=dict(q.metadata)
        svd=(rotation@q.svd_u,q.singular_values,q.svd_vh)
    node=center if emission_node_m is None else np.asarray(emission_node_m,float)
    if node.shape!=(3,) or not np.isfinite(node).all():
        raise ValueError('Emission node must be a finite 3-vector')
    all_vertices=np.concatenate((v,fitted))
    delays=(all_vertices-node)@n/C_M_S
    metadata.update(
      supplied_xyz_energy=energy.tolist(),spatial_retardation_applied=False,
      spatial_delay_offset_min_s=float(delays.min()),spatial_delay_offset_max_s=float(delays.max()),
      spatial_delay_span_s=float(np.ptp(delays)),
      max_abs_spatial_delay_offset_s=float(np.max(np.abs(delays))),
      geometry_delay_error_bound_s=metadata['affine_corner_max_error_m']/C_M_S,
      power_frozen_at_node=True,angular_emissivity_variation_error_bound='Not estimated; no angular beam broadening applied',
      quadrature='Deterministic 2x2 pixel-centre offsets ±0.25 pixel',
    )
    window=_time_window(source_interval_s,arrival_interval_s,metadata)
    window_covers_support=False
    if window is not None:
        _,duration,a0,a1=window
        fitted_delays=fitted@n/C_M_S
        arrival_min=-float(fitted_delays.max())
        arrival_max=duration-float(fitted_delays.min())
        if a1<=arrival_min or a0>=arrival_max:
            metadata.update(method='Exact disjoint static-support arrival bound',arrival_fraction=0.,
                            window_xyz_energy=[0.,0.,0.],retained_fraction=None,
                            spatial_retardation_exact_for_static_affine_ray=True)
            return _empty(np.zeros(3),metadata)
        window_covers_support=a0<=arrival_min and a1>=arrival_max
        metadata['arrival_window_covers_entire_support']=bool(window_covers_support)
    center_pixel=origin+np.array([center@X,-center@Y])/pixel
    fitted_pixel=origin+np.column_stack((fitted@X,-fitted@Y))/pixel
    extent=np.ptp(fitted_pixel,axis=0)
    if np.linalg.norm(extent)<.5:
        return _unresolved_cic(H,center_pixel,X,Y,pixel,energy,width,height,metadata,window=window,center=center,n=n)
    lo=np.floor(fitted_pixel.min(axis=0)-.5).astype(int)
    hi=np.ceil(fitted_pixel.max(axis=0)+.5).astype(int)+1
    nx,ny=(hi-lo).tolist()
    if nx*ny>max_bbox_pixels:
        raise ValueError('Projected cell bounding box exceeds max_bbox_pixels; subdivide physical support')
    if hi[0]<=0 or lo[0]>=width or hi[1]<=0 or lo[1]>=height:
        fraction=1.
        if window is not None and not window_covers_support:
            _,duration,a0,a1=window
            delay=(center@n+_LOSS_NODES@(H.T@n))/C_M_S
            overlap=np.maximum(0.,np.minimum(duration,a1+delay)-np.maximum(0.,a0+delay))/duration
            fraction=float(np.clip(overlap@_LOSS_WEIGHTS,0,1))
        selected=energy*fraction
        metadata.update(method='Disjoint projected bounding box',retained_fraction=0.,
                        pre_normalization_volume_relative_error=None,
                        pre_normalization_projected_area_relative_error=None,
                        quadrature_skipped_because_fully_cropped=True,
                        window_xyz_energy=selected.tolist(),arrival_fraction=fraction,
                        spatial_retardation_exact_for_static_affine_ray=window_covers_support,
                        cropped_energy_timing_quadrature='4x4x4 volume Gauss nodes for invisible cropped energy only; covering arrival bins conserve source energy')
        return _empty(selected,metadata)
    U,s,Vh=np.linalg.svd(H) if svd is None else svd
    threshold=max(pixel*1e-8,float(s.max())*1e-12)
    regularized=np.maximum(s,threshold)
    Hray=(U*regularized)@Vh if np.any(s<threshold) else H
    metadata['degenerate_dimension_regularization_m']=float(max(0,threshold-s.min()))
    metadata['affine_rank_at_relative_1e12']=int(np.sum(s>max(s.max(),1e-300)*1e-12))
    if (prepared_support is not None and prepared_support.inverse_half_edges is not None
            and not np.any(s<threshold)):
        inverse=prepared_support.inverse_half_edges@rotation.T
    else:
        inverse=np.linalg.inv(Hray)
    aa=(inverse@X)*pixel;bb=-(inverse@Y)*pixel;direction=inverse@n
    offsets=np.array([[-.25,-.25],[.25,-.25],[-.25,.25],[.25,.25]])
    xx=(np.arange(lo[0],hi[0])[None,:,None]+offsets[:,0]-center_pixel[0])
    yy=(np.arange(lo[1],hi[1])[:,None,None]+offsets[:,1]-center_pixel[1])
    enter=np.full((ny,nx,4),-np.inf);leave=np.full((ny,nx,4),np.inf)
    parallel_tolerance=32*np.finfo(float).eps*max(np.linalg.norm(direction),1e-300)
    for j in range(3):
        coordinate=aa[j]*xx+bb[j]*yy
        if abs(direction[j])<=parallel_tolerance:
            outside=np.abs(coordinate)>1
            enter[outside]=np.inf;leave[outside]=-np.inf
        else:
            a=(-1-coordinate)/direction[j];b=(1-coordinate)/direction[j]
            np.maximum(enter,np.minimum(a,b),out=enter)
            np.minimum(leave,np.maximum(a,b),out=leave)
    lengths=np.maximum(0.,leave-enter)
    weights=lengths.mean(axis=-1)
    raster_volume=float(weights.sum()*pixel**2)
    if not np.isfinite(raster_volume) or raster_volume<=0:
        return _unresolved_cic(H,center_pixel,X,Y,pixel,energy,width,height,metadata,window=window,center=center,n=n)
    volume=float(8*abs(np.linalg.det(Hray)))
    px=X@Hray;py=Y@Hray
    area=4*(abs(px[0]*py[1]-px[1]*py[0])+abs(px[1]*py[2]-px[2]*py[1])+abs(px[2]*py[0]-px[0]*py[2]))
    raster_area=float(np.mean(lengths>0,axis=-1).sum()*pixel**2)
    metadata.update(method='Uniform affine volume: vectorized ray/slab path lengths',
                    unresolved_support_fallback=False,bounding_box_pixels=int(nx*ny),
                    ray_samples=int(nx*ny*4),quadrature_volume_m3=raster_volume,
                    ray_geometry_volume_m3=volume,
                    pre_normalization_volume_relative_error=float(raster_volume/volume-1),
                    projected_area_m2=float(area),quadrature_projected_area_m2=raster_area,
                    pre_normalization_projected_area_relative_error=float(raster_area/area-1) if area else None)
    window_weights=weights
    if window is not None and not window_covers_support:
        window_weights=_windowed_ray_lengths(enter,leave,window,float(center@n/C_M_S)).mean(axis=-1)
        metadata.update(spatial_retardation_exact_for_static_affine_ray=True,
                        spatial_retardation_method='Exact piecewise-linear source/arrival overlap integral along each static pixel subray; pixel quadrature remains 2x2')
    elif window is not None:
        metadata.update(spatial_retardation_exact_for_static_affine_ray=True,spatial_retardation_method='Arrival interval contains every source-time/depth combination')
    return _from_weights(window_weights,int(lo[0]),int(lo[1]),energy,width,height,metadata,weights)
