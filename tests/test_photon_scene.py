import numpy as np
from pulsar.photon_scene import (TAU,source_frame,emission_event,propagate,
    observer_event,wavefront_crossing,received_flux,systematic_sample)


def test_rotating_source_recovers_original_events():
    rng=np.random.default_rng(120)
    r=rng.normal(size=(500,3));n=rng.normal(size=(500,3));n/=np.linalg.norm(n,axis=1)[:,None]
    t=rng.uniform(-20,20,500)
    r0,n0=source_frame(r,n,t)
    restored=emission_event(r0,n0,t)
    np.testing.assert_allclose(restored[0],r,atol=3e-15)
    np.testing.assert_allclose(restored[1],n,atol=3e-15)


def test_photons_do_not_corotate_after_emission():
    r,n=emission_event(np.array([[1.,0,0]]),np.array([[0.,1,0]]),np.array([.7]))
    points=np.array([propagate(r,n,age)[0] for age in [0.,1.,2.,3.]])
    np.testing.assert_allclose(np.diff(points,axis=0),np.broadcast_to(n,(3,3)),atol=1e-15)
    np.testing.assert_allclose(np.linalg.norm(np.diff(points,axis=0),axis=1),1.)
    # Rotation moves a newly emitted ray, but cannot bend the old one.
    _,later_direction=emission_event(np.array([[1.,0,0]]),np.array([[0.,1,0]]),np.array([1.7]))
    assert np.linalg.norm(later_direction-n)>.9


def test_far_observer_phase_equals_original_sky_equation():
    rng=np.random.default_rng(87)
    r=rng.normal(size=(500,3));n=rng.normal(size=(500,3));n/=np.linalg.norm(n,axis=1)[:,None]
    t=rng.uniform(-10,10,500)
    r0,n0=source_frame(r,n,t)
    te,re,ne,phase=observer_event(r0,n0)
    original=(t-np.arctan2(n[:,1],n[:,0])-np.sum(r*n,axis=1))/TAU%1
    np.testing.assert_allclose((phase-original+.5)%1-.5,0,atol=1e-15)
    np.testing.assert_allclose(ne[:,1],0,atol=1e-15)
    assert np.all(ne[:,0]>0)
    crossing,arrival=wavefront_crossing(re,ne,te,20.)
    np.testing.assert_allclose(np.sum(crossing*ne,axis=1),20.,atol=2e-14)
    np.testing.assert_allclose(((arrival-20.)/TAU-phase+.5)%1-.5,0,atol=1e-15)


def test_receiver_has_flight_delay_and_integrates_frame_shutter():
    edges=np.arange(101)/100;profile=np.zeros(100);profile[20]=10
    # A delay of 0.3 turns places arrival at lab phase 0.505.
    assert np.isclose(received_flux(edges,profile,.505,.3*TAU,.01)[0],10)
    assert np.isclose(received_flux(edges,profile,.505,.3*TAU,.02)[0],5)
    assert received_flux(edges,profile,.205,.3*TAU,.01)[0]==0
    phases=np.arange(100)/100+.005
    flux=received_flux(edges,profile,phases,.3*TAU,.01)
    assert np.argmax(flux)==50
    np.testing.assert_allclose(flux.mean(),profile.mean(),atol=1e-14)


def test_packet_sampling_preserves_total_number_without_changing_directions():
    ids,weight=systematic_sample([0.,1.,2.,0.,7.],100)
    assert set(ids)=={1,2,4}
    np.testing.assert_array_equal(np.bincount(ids,minlength=5),[0,10,20,0,70])
    assert weight*len(ids)==10
