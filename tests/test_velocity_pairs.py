"""Physical sign, degenerate-gauge and provenance checks of PAW charge pairs."""
import copy
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tools'))
from kubo_pairs import pair_hall_spectrum, read_pairs, write_pairs
from velocity_pairs import velocity_pair_data


def example(double=False):
    # One or two identical flat occupied/empty sectors with D_x=sigma_x,
    # D_y=sigma_y: each occupied curvature is -1/2 A^2.
    count = 2 if double else 1
    nb = 2*count
    e = np.tile([-1.]*count+[1.]*count, (4,1))
    d = np.zeros((4,3,nb,nb), dtype=np.complex128)
    for i in range(count):
        d[:,0,i,i+count] = d[:,0,i+count,i] = 1
        d[:,1,i,i+count] = -1j
        d[:,1,i+count,i] = 1j
    data = dict(energies_eV=e, velocity_eVA=d, lattice_A=2*np.pi*np.eye(3),
        kpoints_fractional=np.array([[0.,0.,0.],[0.,.5,0.],[.5,0.,0.],[.5,.5,0.]]),
        weights=np.full(4,.25), band_indices=np.arange(1,nb+1))
    metadata = dict(operator_accuracy='paw_full_velocity', physical_paw_validated=True,
        producer_status='PASS',full_source_band_coverage=True,diagonal_available=True,
        degenerate_blocks_available=True,normalization_applied=False,source_nbands=nb,nkpoints=4,
        producer='Analytic operator fixture, not material data',operator_scope='unit test',
        source_wavecar_sha256='a'*64,source_sha256={'WAVECAR':'a'*64})
    return data, metadata


def pairs(data, meta, **kw):
    return velocity_pair_data(data,meta,sampling=kw.get('sampling',
        {'kind':'uniform_full_2d','mesh':[2,2],'plane_axes':[0,1]}),energy_reference='test zero')


def sigma(data):
    rows,_ = pair_hall_spectrum(data,[0.],[0.],mu_reference=0.)
    return next(r['sigma_e2_over_h'] for r in rows if r['region']=='total')


class VelocityPairTests(unittest.TestCase):
    def test_analytic_sign_units_and_no_spin_factor(self):
        data,meta=example()
        result=pairs(data,meta)
        self.assertAlmostEqual(sigma(result),1/(4*np.pi),places=14)
        self.assertEqual(result.metadata['source_operator']['kind'],'paw_full_velocity')
        self.assertEqual(result.metadata['spin_multiplicity'],1)
        self.assertEqual(result.metadata['source_wavecar_sha256'],'a'*64)

    def test_complete_degenerate_subspaces_are_unitary_invariant(self):
        data,meta=example(True)
        expected=1/(2*np.pi)
        rng=np.random.default_rng(317)
        for k in range(4):
            u=np.zeros((4,4),dtype=complex)
            for sl in (slice(0,2),slice(2,4)):
                x=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
                u[sl,sl]=np.linalg.qr(x)[0]
            data['velocity_eVA'][k]=u.conj().T@data['velocity_eVA'][k]@u
        self.assertAlmostEqual(sigma(pairs(data,meta)),expected,places=14)

    def test_time_reversal_conjugation_reverses_charge_curvature(self):
        data,meta=example()
        original=sigma(pairs(data,meta))
        data['velocity_eVA']=-data['velocity_eVA'].conj()
        self.assertAlmostEqual(sigma(pairs(data,meta)),-original,places=14)

    def test_roundtrip_preserves_operator_and_original_energies(self):
        data,meta=example(True)
        result=pairs(data,meta)
        with tempfile.TemporaryDirectory() as t:
            out=Path(t)/'pairs'
            write_pairs(out,result)
            loaded=read_pairs(out)
            np.testing.assert_array_equal(loaded.energies_eV,data['energies_eV'])
            self.assertEqual(loaded.metadata['source_operator'],result.metadata['source_operator'])
            self.assertEqual(sigma(loaded),sigma(result))

    def test_incomplete_or_nonphysical_operator_is_rejected(self):
        data,meta=example()
        for key,value in [('diagonal_available',False),('degenerate_blocks_available',False),
                          ('normalization_applied',True),('operator_accuracy','canonical_momentum')]:
            bad=copy.deepcopy(meta);bad[key]=value
            with self.assertRaisesRegex(ValueError,'complete physical'):
                pairs(data,bad)
        data['velocity_eVA'][0,0,0,1]+=.1j
        with self.assertRaisesRegex(ValueError,'Hermitian'):
            pairs(data,meta)

    def test_wrong_mesh_and_split_source_are_rejected(self):
        data,meta=example()
        with self.assertRaises(ValueError):
            pairs(data,meta,sampling={'kind':'uniform_full_2d','mesh':[3,3],'plane_axes':[0,1]})
        data['band_indices']=np.array([2,3])
        with self.assertRaisesRegex(ValueError,'all source bands'):
            pairs(data,meta)


if __name__=='__main__':
    unittest.main()
