import hashlib
import json
from pathlib import Path
import unittest
import numpy as np
from astropy.io import fits
from pulsar.observations import load_vela,compare_to_vela

ROOT=Path(__file__).resolve().parents[1]


class ObservationChecks(unittest.TestCase):
    def test_original_data_hash_and_rebin_conservation(self):
        path=ROOT/'data/J0835-4510_3PC_data.fits'
        provenance=json.loads((ROOT/'data/provenance.json').read_text())
        self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),provenance['files'][path.name])
        data=load_vela(path)
        self.assertTrue(data['duplicate_display_cycle_removed'])
        self.assertEqual(len(data['counts']),200)
        with fits.open(path) as f:
            d=f['GAMMA_LC'].data
            raw=np.array(d['GT100_WtCnt'][:800],float)
            errors=np.array(d['Unc_GT100_WtCnt'][:800],float)
        self.assertAlmostEqual(data['counts'].sum(),raw.sum(),places=6)
        self.assertAlmostEqual(np.sum(data['error']**2),np.sum(errors**2),places=6)

    def test_phase_scale_and_constant_recovered_from_independent_bin_integrals(self):
        phase=np.arange(4096)/4096
        profile=1+.5*np.cos(2*np.pi*phase)+.2*np.sin(4*np.pi*phase)
        edges=np.linspace(0,1,65); lo=edges[:-1]; hi=edges[1:]; shift=.125
        # Exact antiderivatives of the synthetic source, independent of the fitter.
        primitive=lambda x:x+.5*np.sin(2*np.pi*(x-shift))/(2*np.pi)-.2*np.cos(4*np.pi*(x-shift))/(4*np.pi)
        averaged=(primitive(hi)-primitive(lo))/(hi-lo)
        data={'phase_lo':lo,'phase_hi':hi,'counts':2*averaged+5,'error':np.ones(64)*.1}
        fit=compare_to_vela(phase,profile,data,n_shifts=256)
        self.assertEqual(fit['phase_shift_cycles'],shift)
        self.assertAlmostEqual(fit['scale'],2,places=4)
        self.assertAlmostEqual(fit['constant_component'],5,places=4)
        self.assertLess(fit['chi_squared'],1e-6)


if __name__=='__main__':
    unittest.main()
