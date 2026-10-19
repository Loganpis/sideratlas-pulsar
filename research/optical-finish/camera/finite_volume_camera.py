"""Finite-support optical camera with explicitly piecewise-constant emissivity.

Directional powers remain those actually calculated at each quadrature node.
A uniform affine approximation to each injection-coordinate support supplies
spatial extent. Observer delays across this support are integrated, not added
as visual trails. Rigid motion is midpoint-quadrature within source intervals.
This does not make the physical emissivity or spatial sampling converged.
"""
import numpy as np
from scipy.ndimage import gaussian_filter

from render_optical import OpticalCamera
from event_camera import rotating_pattern_arrival_edges,periodic_interval_weights
from surface_camera import C_M_S
from volume_raster import prepare_affine_support,rasterize_prepared_cell,SIGNS
from compact_render_geometry import compact_render_events


class FiniteVolumeCamera(OpticalCamera):
    def __init__(self,*args,motion_arc_pixels=2.,omitted_xyz_fraction=1e-7,
                 center_fraction=(.5,.53),**kwargs):
        super().__init__(*args,**kwargs)
        self.center_fraction=tuple(center_fraction)
        if len(self.center_fraction)!=2 or not np.isfinite(self.center_fraction).all():
            raise ValueError('Two finite physical image-centre fractions required')
        if not np.isfinite(motion_arc_pixels) or motion_arc_pixels<=0:
            raise ValueError('Positive motion-quadrature arc required')
        if not 0<=omitted_xyz_fraction<1e-3:
            raise ValueError('Omission bound must lie in [0,1e-3)')
        self.motion_arc=motion_arc_pixels
        energy=self.power*self.number[:,None]*self.dt[:,None]
        self.all_source_energy=energy.sum(axis=0)
        score=np.sum(energy/np.maximum(self.all_source_energy,np.finfo(float).tiny),axis=1)
        order=np.argsort(score)[::-1]
        cumulative=np.cumsum(energy[order],axis=0)
        needed=max(np.searchsorted(cumulative[:,j],(1-omitted_xyz_fraction)*self.all_source_energy[j])+1
                   for j in range(3) if self.all_source_energy[j]>0)
        chosen=np.sort(order[:needed])
        self.omitted_energy=self.all_source_energy-energy[chosen].sum(axis=0)
        self.omitted_fraction=np.divide(self.omitted_energy,self.all_source_energy,
                                       out=np.zeros(3),where=self.all_source_energy>0)
        self.position=self.position[chosen];self.time=self.time[chosen];self.dt=self.dt[chosen]
        self.number=self.number[chosen];self.power=self.power[chosen];self.cell_index=self.cell_index[chosen]
        self.source_energy=energy[chosen].sum(axis=0)
        with np.load(self.scene,allow_pickle=False) as saved:
            reference=saved['reference_cell_positions_m']
        arrays,self.compaction_report=compact_render_events(
            reference_positions_m=reference,reference_vertices_m=self.reference_vertices,
            cell_index=self.cell_index,source_time_s=self.time,source_interval_s=self.dt,
            number_weights=self.number,xyz_lm_sr=self.power,grid_m=.25*self.pixel_m,
            period_s=self.period,positions_m=self.position)
        self.position=arrays['positions_m'];self.time=arrays['source_time_s']
        self.dt=arrays['source_interval_s'];self.number=arrays['number_weights']
        self.power=arrays['xyz_lm_sr'];self.cell_index=arrays['cell_index']
        self.reference_vertices=arrays['cell_footprint_vertices_reference_m']
        reference=arrays['reference_cell_positions_m']
        self.supports={}
        # Both the original and fitted affine corners are enclosed by this
        # rotation-invariant displacement bound. Candidate bins can therefore
        # include the full spatial-delay support, even across period wrap.
        v=self.reference_vertices
        center=v.mean(axis=1)
        H=np.einsum('vi,nvj->nji',SIGNS,v-center[:,None,:])/8
        fitted=center[:,None,:]+np.einsum('vi,nji->nvj',SIGNS,H)
        bound=np.maximum(np.linalg.norm(v-reference[:,None,:],axis=-1).max(axis=1),
                         np.linalg.norm(fitted-reference[:,None,:],axis=-1).max(axis=1))/C_M_S
        self.support_bound_s=bound
        self.support_radius=np.maximum(np.linalg.norm(v[:,:,:2],axis=-1).max(axis=1),
                                      np.linalg.norm(fitted[:,:,:2],axis=-1).max(axis=1))
        left,right=rotating_pattern_arrival_edges(self.time,self.dt,self.position,self.observer,self.omega)
        ids,bins,_,a,b=periodic_interval_weights(left-bound[self.cell_index],right+bound[self.cell_index],
                                               self.period,self.bins,return_edges=True)
        sort=np.argsort(bins,kind='stable')
        self.event_ids=ids[sort]
        width=self.period/self.bins
        self.arrival_left=(np.floor(((a+b)/2)/width)*width)[sort]
        self.offset=np.r_[0,np.cumsum(np.bincount(bins,minlength=self.bins))]
        self.volume_records={}
        self.exposure=None

    def support(self,index):
        index=int(index)
        if index not in self.supports:
            self.supports[index]=prepare_affine_support(self.reference_vertices[index])
        return self.supports[index]

    def radiance(self,bin_index):
        index=int(bin_index)%self.bins
        sl=slice(self.offset[index],self.offset[index+1])
        ids=self.event_ids[sl];arrival=self.arrival_left[sl]
        frame=np.zeros((self.height,self.width,3))
        lost=np.zeros(3);window=np.zeros(3)
        max_geometry=0.;max_volume_error=0.;max_raster_error=0.;pieces_count=0;culled=0
        for event,a0 in zip(ids,arrival):
            cell=int(self.cell_index[event])
            support=self.support(cell)
            count=max(1,int(np.ceil(self.omega*self.support_radius[cell]*self.dt[event]/(self.pixel_m*self.motion_arc))))
            if count>1024:raise ValueError('Refine source-time quadrature; excessive motion subintervals')
            edges=np.linspace(self.time[event]-self.dt[event]/2,self.time[event]+self.dt[event]/2,count+1)
            energy=self.power[event]*self.number[event]*(self.dt[event]/count)
            angle=self.omega*(edges[:-1]+edges[1:])/2
            co,si=np.cos(angle),np.sin(angle);n=self.observer
            body_n=np.stack((co*n[0]+si*n[1],-si*n[0]+co*n[1],np.full(count,n[2])),axis=1)
            projection=body_n@support.center_m
            extent=np.abs(body_n@support.half_edges_m).sum(axis=1)
            minimum=edges[:-1]-(projection+extent)/C_M_S
            maximum=edges[1:]-(projection-extent)/C_M_S
            relevant=(a0+self.period/self.bins>minimum)&(a0<maximum)
            culled+=int(np.count_nonzero(~relevant))
            for t0,t1 in zip(edges[:-1][relevant],edges[1:][relevant]):
                patch=rasterize_prepared_cell(support,energy,self.observer,
                        rotation_z_rad=self.omega*(t0+t1)/2,
                        meters_per_pixel=self.pixel_m,width=self.width,height=self.height,roll_deg=self.roll,
                        center_fraction=self.center_fraction,
                        source_interval_s=(t0,t1),arrival_interval_s=(a0,a0+self.period/self.bins))
                patch.deposit_into(frame)
                lost+=patch.lost_xyz
                window+=np.asarray(patch.metadata['window_xyz_energy'])
                md=patch.metadata
                max_geometry=max(max_geometry,md.get('affine_corner_max_error_m',0.))
                e=md.get('affine_volume_relative_error_vs_trilinear')
                if e is not None:max_volume_error=max(max_volume_error,abs(e))
                e=md.get('pre_normalization_volume_relative_error')
                if e is not None:max_raster_error=max(max_raster_error,abs(e))
                pieces_count+=1
        self.volume_records[index]={'arrival_xyz_energy':window.tolist(),'cropped_xyz_energy':lost.tolist(),
            'affine_corner_error_max_m':max_geometry,'affine_volume_relative_error_max':max_volume_error,
            'raw_raster_volume_relative_error_max':max_raster_error,'motion_pieces':pieces_count,
            'zero_arrival_overlap_pieces_analytically_culled':culled}
        self.crop_records[index]={'source_Y_energy':float(window[1]),'cropped_Y_energy':float(lost[1]),
            'retained_fraction':float(1-lost[1]/window[1]) if window[1]>0 else None}
        frame/=self.pixel_m**2*(self.period/self.bins)
        if self.surface_xyz is not None:
            frame+=self.surface_xyz((index+.5)*self.period/self.bins)
        if self.psf>0:
            before=frame.sum(axis=(0,1))
            frame=gaussian_filter(frame,(self.psf,self.psf,0),mode='constant',truncate=4)
            self.volume_records[index]['psf_edge_lost_xyz_radiance_sum']=(before-frame.sum(axis=(0,1))).tolist()
        return frame

    def volume_metadata(self):
        out={'spatial_model':'Uniform affine supports with node-calculated directional emissivity held constant within support',
             'source_xyz_energy':self.all_source_energy.tolist(),'retained_source_xyz_energy':self.source_energy.tolist(),
             'omitted_source_xyz_energy':self.omitted_energy.tolist(),'omitted_source_xyz_fraction':self.omitted_fraction.tolist(),
             'motion_quadrature_max_arc_pixels':self.motion_arc,
             'center_fraction':list(self.center_fraction),'subpixel_geometry_aggregation':self.compaction_report,
             'finite_support_retardation':'Analytic source/arrival-window overlap integrated along each pixel subray in a frozen motion-quadrature piece',
             'maximum_node_relative_spatial_delay_bound_s':float(self.support_bound_s[self.cell_index].max()),
             'piecewise_emissivity_error':'Node radiation and affine geometry are unresolved approximations; no calibrated or converged image claim',
             'rendered_detector_bins':len(self.volume_records),'prepared_supports':len(self.supports),
             'records':self.volume_records}
        if len(self.volume_records)==self.bins:
            arrival=np.sum([x['arrival_xyz_energy'] for x in self.volume_records.values()],axis=0)
            out['arrival_xyz_energy']=arrival.tolist()
            out['full_cycle_relative_energy_error']=float(np.linalg.norm(arrival-self.source_energy)/np.linalg.norm(self.source_energy))
        return out
