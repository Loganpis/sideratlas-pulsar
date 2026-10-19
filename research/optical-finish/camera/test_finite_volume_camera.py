import sys
from pathlib import Path
import numpy as np
sys.path[:0]=[str(Path(__file__).parent),str(Path(__file__).resolve().parents[3]/'src')]
from finite_volume_camera import FiniteVolumeCamera
from surface_camera import rotate_vectors
from volume_raster import SIGNS


def test_complete_period_energy_with_finite_support_delays(tmp_path):
    p=tmp_path/'events.npz'
    reference=np.array([[10000.,0,0]])
    phase=np.array([.25,.75])*2*np.pi
    position=rotate_vectors(np.repeat(reference,2,axis=0),[0,0,1],phase)
    np.savez(p,positions_m=position,source_time_s=np.array([.25,.75]),source_interval_s=np.full(2,.5),
             number_weights=np.ones(2),xyz_lm_sr=np.ones((2,3))*100,observer_direction=np.array([1.,0,0]),
             period_s=1.,cell_index=np.zeros(2,dtype=int),reference_cell_positions_m=reference,
             cell_footprint_vertices_reference_m=reference[:,None,:]+SIGNS[None,:,:]*1000)
    camera=FiniteVolumeCamera(p,width=100,arrival_bins=8,meters_per_pixel=1000,psf_sigma_pixels=0,
                              motion_arc_pixels=3,omitted_xyz_fraction=0)
    for i in range(8):camera.radiance(i)
    md=camera.volume_metadata()
    assert md['full_cycle_relative_energy_error']<1e-8
    assert np.max(np.abs(np.asarray(md['arrival_xyz_energy'])-100))<1e-6
    assert np.max(md['omitted_source_xyz_fraction'])<1e-14
