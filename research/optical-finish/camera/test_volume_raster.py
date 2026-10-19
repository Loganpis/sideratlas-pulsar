import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).parent))
from volume_raster import rasterize_affine_cell,fit_affine_support,SIGNS


def box(half=(1.5,1.5,1.),center=(0,0,0)):
    return np.asarray(center)+SIGNS*np.asarray(half)


def draw(vertices,energy=(2.,3.,5.),**kwargs):
    return rasterize_affine_cell(vertices,energy,[0,0,1],meters_per_pixel=1.,width=8,height=8,roll_deg=0,center_fraction=(.5,.5),**kwargs)


def canvas(p):
    x=np.zeros((8,8,3));p.deposit_into(x);return x


def test_analytic_axis_box_silhouette_and_absolute_integral():
    p=draw(box())
    image=canvas(p)
    expected=np.zeros((8,8,3));expected[3:6,3:6]=np.array([2,3,5])/9
    np.testing.assert_allclose(image,expected,atol=2e-15)
    np.testing.assert_allclose(p.lost_xyz,0,atol=2e-15)
    assert abs(p.metadata['pre_normalization_volume_relative_error'])<1e-14
    assert abs(p.metadata['pre_normalization_projected_area_relative_error'])<1e-14


def test_crop_ledger_and_half_box():
    p=rasterize_affine_cell(box(),[2,3,5],[0,0,1],meters_per_pixel=1.,width=8,height=8,roll_deg=0,center_fraction=(-.5/8,.5))
    np.testing.assert_allclose(p.patch_xyz.sum(axis=(0,1)),np.array([2,3,5])*.5,atol=2e-14)
    np.testing.assert_allclose(p.patch_xyz.sum(axis=(0,1))+p.lost_xyz,[2,3,5],atol=2e-14)
    outside=draw(box(center=(100,100,0)))
    assert outside.patch_xyz.size==0
    np.testing.assert_array_equal(outside.lost_xyz,[2,3,5])


def test_split_aligned_volume_is_additive():
    whole=canvas(draw(box()))
    # Split along physical x; volume-proportional energies. The cut crosses
    # the middle detector pixel, where each half receives two of four rays.
    a=canvas(draw(box(half=(.75,1.5,1),center=(-.75,0,0)),energy=(1,1.5,2.5)))
    b=canvas(draw(box(half=(.75,1.5,1),center=(.75,0,0)),energy=(1,1.5,2.5)))
    np.testing.assert_allclose(a+b,whole,atol=2e-15)


def test_energy_linear_in_particle_split_for_oblique_cell():
    H=np.array([[1.7,.4,.2],[.2,1.3,.3],[.4,.3,.9]])
    vertices=SIGNS@H.T+[.2,.3,.1]
    p=rasterize_affine_cell(vertices,[7,2,3],[1,.5,.8],meters_per_pixel=.3,width=48,height=64,roll_deg=32)
    q=rasterize_affine_cell(vertices,[.7,.2,.3],[1,.5,.8],meters_per_pixel=.3,width=48,height=64,roll_deg=32)
    np.testing.assert_allclose(p.patch_xyz,10*q.patch_xyz,atol=2e-14)
    np.testing.assert_allclose(p.patch_xyz.sum(axis=(0,1))+p.lost_xyz,[7,2,3],atol=3e-14)
    assert p.metadata['angular_emissivity_variation_error_bound'].startswith('Not estimated')


def test_point_limit_and_exact_flat_sheet_are_finite():
    point=box(half=(0,0,0),center=(.2,.3,1))
    tiny=box(half=(1e-12,1e-12,1e-12),center=(.2,.3,1))
    np.testing.assert_allclose(canvas(draw(point)),canvas(draw(tiny)),atol=2e-14)
    flat=draw(box(half=(1.5,1.5,0)))
    np.testing.assert_allclose(canvas(flat),canvas(draw(box())),rtol=1e-8,atol=1e-13)
    assert flat.metadata['degenerate_dimension_regularization_m']>0
    np.testing.assert_allclose(flat.patch_xyz.sum(axis=(0,1))+flat.lost_xyz,[2,3,5],atol=2e-14)


def test_affine_fit_reports_nonaffine_geometry_and_retarded_bound():
    v=box();v[7]+=[.2,.1,.3]
    _,_,_,m=fit_affine_support(v)
    assert m['affine_corner_max_error_m']>0
    p=draw(v,emission_node_m=[0,0,0])
    assert not p.metadata['spatial_retardation_applied']
    assert p.metadata['max_abs_spatial_delay_offset_s']>=1.3/299792458.
    assert p.metadata['geometry_delay_error_bound_s']>0


def test_invalid_shapes_and_huge_support_rejected():
    with pytest.raises(ValueError):draw(np.zeros((7,3)))
    with pytest.raises(ValueError):draw(box(),energy=(-1,2,3))
    with pytest.raises(ValueError,match='bounding box'):
        draw(box(half=(1000,1000,1)),max_bbox_pixels=100)


def test_prepared_geometry_rotates_like_fresh_corners():
    from volume_raster import prepare_affine_support,rasterize_prepared_cell
    H=np.array([[13.,3.,2.],[2.,7.,1.],[1.,2.,4.]])
    v=SIGNS@H.T+[20,10,4]
    support=prepare_affine_support(v)
    phase=.47;co,si=np.cos(phase),np.sin(phase)
    rot=np.array([[co,-si,0],[si,co,0],[0,0,1.]])
    kwargs=dict(meters_per_pixel=1.,width=1080,height=1920,roll_deg=-55.)
    a=rasterize_affine_cell(v@rot.T,[1,2,3],[1,.3,.7],**kwargs)
    b=rasterize_prepared_cell(support,[1,2,3],[1,.3,.7],rotation_z_rad=phase,**kwargs)
    assert a.image_slice==b.image_slice
    np.testing.assert_allclose(a.patch_xyz,b.patch_xyz,atol=1e-14)
    np.testing.assert_allclose(a.lost_xyz,b.lost_xyz,atol=1e-14)
    assert abs(a.metadata['spatial_delay_span_s']-b.metadata['spatial_delay_span_s'])<1e-20


def test_edge_on_flat_support_uses_disclosed_line_reconstruction():
    p=rasterize_affine_cell(box(half=(1.5,1.5,0)),[2,3,5],[1,0,0],
                            meters_per_pixel=1.,width=8,height=8,roll_deg=0,
                            center_fraction=(.5,.5))
    assert p.metadata['unresolved_support_fallback']
    assert p.metadata['midpoint_nodes']>1
    np.testing.assert_allclose(p.patch_xyz.sum(axis=(0,1))+p.lost_xyz,[2,3,5],atol=2e-14)


def test_static_slab_time_window_has_exact_analytic_fraction():
    from volume_raster import C_M_S
    p=draw(box(half=(1.5,1.5,.1*C_M_S)),source_interval_s=(0.,.02),arrival_interval_s=(-.005,.025))
    # Arrival PDF is constant 1/.2 s in this central interval; width=.03 s.
    np.testing.assert_allclose(p.window_xyz_energy,np.array([2,3,5])*.15,atol=2e-14)
    np.testing.assert_allclose(p.patch_xyz.sum(axis=(0,1))+p.lost_xyz,p.window_xyz_energy,atol=2e-14)
    assert p.metadata['spatial_retardation_exact_for_static_affine_ray']


def test_periodic_arrival_bins_cover_all_source_energy_without_renormalizing_window():
    from volume_raster import C_M_S,prepare_affine_support,rasterize_prepared_cell
    support=prepare_affine_support(box(half=(1.5,1.5,.1*C_M_S)))
    image=np.zeros((8,8,3));lost=np.zeros(3);window=np.zeros(3)
    for cycle in (-1,0,1):
        for i in range(10):
            p=rasterize_prepared_cell(support,[2,3,5],[0,0,1],meters_per_pixel=1.,width=8,height=8,
                                     roll_deg=0,center_fraction=(.5,.5),source_interval_s=(0.,1.),
                                     arrival_interval_s=(cycle+i/10,cycle+(i+1)/10))
            p.deposit_into(image);lost+=p.lost_xyz;window+=p.window_xyz_energy
    np.testing.assert_allclose(window,[2,3,5],atol=2e-14)
    np.testing.assert_allclose(image.sum(axis=(0,1))+lost,[2,3,5],atol=2e-14)


def test_nearer_support_arrives_earlier_and_zero_overlap_stays_zero():
    from volume_raster import C_M_S
    near=draw(box(half=(1.5,1.5,.005*C_M_S),center=(0,0,.2*C_M_S)),
              source_interval_s=(0.,.02),arrival_interval_s=(-.21,-.17))
    far=draw(box(half=(1.5,1.5,.005*C_M_S),center=(0,0,-.2*C_M_S)),
             source_interval_s=(0.,.02),arrival_interval_s=(-.21,-.17))
    np.testing.assert_allclose(near.window_xyz_energy,[2,3,5],atol=2e-14)
    np.testing.assert_array_equal(far.window_xyz_energy,[0,0,0])
    np.testing.assert_array_equal(far.patch_xyz,0)
    np.testing.assert_array_equal(far.lost_xyz,[0,0,0])


def test_time_window_crop_loss_is_only_window_energy():
    from volume_raster import C_M_S
    p=rasterize_affine_cell(box(half=(1.5,1.5,.1*C_M_S)),[2,3,5],[0,0,1],
                            meters_per_pixel=1.,width=8,height=8,roll_deg=0,center_fraction=(-.5/8,.5),
                            source_interval_s=(0.,.02),arrival_interval_s=(-.005,.025))
    expected=np.array([2,3,5])*.15
    np.testing.assert_allclose(p.window_xyz_energy,expected,atol=2e-14)
    np.testing.assert_allclose(p.lost_xyz,expected*.5,atol=2e-14)
    np.testing.assert_allclose(p.patch_xyz.sum(axis=(0,1)),expected*.5,atol=2e-14)


def test_window_time_origin_shift_and_point_limit():
    from volume_raster import C_M_S
    v=box(half=(1.5,1.5,.1*C_M_S))
    a=draw(v,source_interval_s=(0.,.02),arrival_interval_s=(-.005,.025))
    b=draw(v,source_interval_s=(1000.,1000.02),arrival_interval_s=(999.995,1000.025))
    np.testing.assert_allclose(a.window_xyz_energy,b.window_xyz_energy,atol=3e-12)
    point=draw(box(half=(0,0,0),center=(0,0,.2*C_M_S)),
               source_interval_s=(0.,.02),arrival_interval_s=(-.2,-.19))
    np.testing.assert_allclose(point.window_xyz_energy,np.array([2,3,5])*.5,atol=2e-14)
    with pytest.raises(ValueError,match='source-time'):
        draw(v,arrival_interval_s=(0,1))


def test_fully_cropped_window_ledger_conserves_across_arrival_bins():
    from volume_raster import C_M_S
    v=box(half=(1.5,1.5,.1*C_M_S),center=(100.,0,0))
    edges=np.linspace(-.2,.2,17)
    total=np.zeros(3)
    for a,b in zip(edges[:-1],edges[1:]):
        p=draw(v,source_interval_s=(0.,.02),arrival_interval_s=(a,b))
        np.testing.assert_array_equal(p.patch_xyz,0)
        np.testing.assert_allclose(p.window_xyz_energy,p.lost_xyz)
        total+=p.lost_xyz
    np.testing.assert_allclose(total,[2,3,5],atol=2e-14)
