"""Conservative transfer checks independent of an equilibrium solution."""
import numpy as np
from conservative_magnetic_column import ConservativeColumn
from magnetic_column import Parameters,SIGMA
from export_conservative_surface import emergent


def test_conservative_energy_ledger_arbitrary_profile():
    col=ConservativeColumn(Parameters(depth_points=16,energy_points=20,angles_per_hemisphere=3,column_max=1e3),10)
    t=np.geomspace(2e5,5e6,len(col.mass));out=col.transfer(t,True)
    jump=np.diff(out['flux']);heat=4*np.pi*out['heating']*col.width
    target=SIGMA*col.p.teff_k**4
    np.testing.assert_allclose(jump/target,heat/target,atol=1e-7,rtol=3e-9)


def test_explicit_angle_export_matches_quadrature_rays():
    col=ConservativeColumn(Parameters(depth_points=12,energy_points=12,angles_per_hemisphere=3,column_max=1e3),8)
    t=np.geomspace(2e5,5e6,len(col.mass));out=col.transfer(t,True)
    explicit=emergent(col,t,out['u'],col.mu)
    quadrature=out['surface_i_nu_cgs'].reshape((-1,2,len(col.mu))).sum(axis=1)
    np.testing.assert_allclose(explicit,quadrature,atol=0,rtol=2e-13)


def test_grazing_limit_is_positive_and_continuous():
    col=ConservativeColumn(Parameters(depth_points=12,energy_points=12,angles_per_hemisphere=3,column_max=1e3),8)
    t=np.geomspace(2e5,5e6,len(col.mass));out=col.transfer(t,True)
    # The tenuous outermost X-mode cell may have optical depth far below
    # 1e-12; use an actual asymptotic mu to test this mathematical limit.
    intens=emergent(col,t,out['u'],np.array([0,1e-50,.1,1]))
    assert np.all(np.isfinite(intens)) and np.all(intens>=0)
    np.testing.assert_allclose(intens[:,0],intens[:,1],rtol=1e-5,atol=1e-100)
