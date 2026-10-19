import sys
from pathlib import Path
import numpy as np
import pytest
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'src'))
from atmosphere_camera import TabulatedAtmosphere,SurfaceComponent


def table(tmp_path,**changes):
    p=tmp_path/'surface.npz'
    values=dict(wavelength_m=np.linspace(300,650,100)*1e-9,mu=np.linspace(0,1,31),
                i_lambda_si=np.full((100,31),2e10),units='W m^-2 sr^-1 m^-1',
                numerical_equilibrium_pass=True,scientific_validation_pass=False)
    values.update(changes)
    np.savez(p,**values)
    return p


def test_units_and_explicit_approximation(tmp_path):
    p=table(tmp_path)
    with pytest.raises(ValueError,match='opt-in'):TabulatedAtmosphere(p)
    a=TabulatedAtmosphere(p,allow_unvalidated_model=True)
    assert np.allclose(a(wavelength_emitted_nm=np.array([[400,500]]),cos_emission_angle=[.5]),20)
    with pytest.raises(ValueError):a(wavelength_emitted_nm=np.array([[700,800]]),cos_emission_angle=[.5])


def test_reject_failed_equilibrium_or_missing_angles(tmp_path):
    with pytest.raises(ValueError,match='equilibrium'):
        TabulatedAtmosphere(table(tmp_path,numerical_equilibrium_pass=False),allow_unvalidated_model=True)
    with pytest.raises(ValueError,match='angular support'):
        TabulatedAtmosphere(table(tmp_path,mu=np.linspace(.001,1,31)),allow_unvalidated_model=True)


def test_resolved_area_and_surface_rotation_constant(tmp_path):
    a=TabulatedAtmosphere(table(tmp_path),allow_unvalidated_model=True)
    c=SurfaceComponent(a,width=100,height=180,meters_per_pixel=1000,
                        observer_direction=[.9,0,.4358898944],period_s=.089371,subpixels=5)
    first=c(0)
    assert first.shape==(180,100,3)
    assert np.all(first>=0)
    assert first.sum()>0
    assert np.array_equal(first,c(.04))
    assert np.count_nonzero(first[:,:,1])<1000
    assert c.metadata()['scientific_validation_pass'] is False
