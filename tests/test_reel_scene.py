"""Energy conservation in the movie's display-grid transfer."""
from pathlib import Path
import sys,unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from prepare_reel_scene import deposit_voxels,continuous_display_volume


class DisplayDepositionTests(unittest.TestCase):
    def test_conserves_energy_and_signed_color_moment(self):
        rng=np.random.default_rng(23);p=rng.uniform(-1.8,1.8,(103,3))
        energy=rng.uniform(0,4,103);values=np.column_stack((energy,energy*rng.uniform(-2,2,103)))
        grid=np.zeros((25,25,25,2))
        deposit_voxels(grid,p,values,2.)
        np.testing.assert_allclose(grid.sum(axis=(0,1,2)),values.sum(axis=0),rtol=1e-13,atol=1e-12)
        self.assertTrue(np.all(grid[...,0]>=0))

    def test_rejects_unmodeled_out_of_volume_emission(self):
        with self.assertRaises(ValueError):
            deposit_voxels(np.zeros((25,25,25,2)),np.array([[3.,0.,0.]]),np.ones((1,2)),2.)

    def test_angular_display_keeps_energy_and_constant_energy_color(self):
        v=np.zeros((25,25,25,2));v[18,12,12]=[2,1];v[12,18,12]=[3,1.5]
        out,meta=continuous_display_volume(v,32)
        self.assertAlmostEqual(out[...,0].sum(),5)
        np.testing.assert_allclose(out[...,1],.5*out[...,0],rtol=1e-12,atol=1e-12)
        self.assertTrue(np.all(out[...,0]>=0))
        self.assertGreater(np.count_nonzero(out[...,0]),np.count_nonzero(v[...,0]))
        self.assertAlmostEqual(meta['half_width_deg'],11.25)
