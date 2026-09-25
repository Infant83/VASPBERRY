"""Independent Pauli, physical metric, provenance and cache checks."""
from dataclasses import replace
import json
from pathlib import Path
import struct
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
import spin_operators as so

PAULI = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]])
BASIS = 'Cartesian Pauli basis: up/down along +z; VASP SAXIS=(0,0,1)'


def wavecar_fixture(path, coefficients):
    """Real RTAG45200 binary, cubic1A cell/1eV cutoff: only Gamma G=0."""
    c = np.asarray(coefficients, dtype='<c8')
    n = len(c); stride = 256
    data = bytearray((3+n)*stride)
    struct.pack_into('<3d', data, 0, stride, 1., 45200.)
    struct.pack_into('<12d', data, stride, 1., n, 1., *np.eye(3).ravel())
    headers = [2., 0., 0., 0.]
    for i in range(n): headers.extend([float(i), 0., 1.])
    struct.pack_into('<'+str(len(headers))+'d', data, 2*stride, *headers)
    for i, row in enumerate(c):data[(3+i)*stride:(3+i)*stride+row.nbytes]=row.tobytes()
    path.write_bytes(data)
    return c.astype(np.complex128)


def aug_metadata(source):
    return dict(schema=so.AUGMENTATION_SCHEMA, version=1, complete=True, producer_status='PASS',
        source_wavecar_sha256=source.metadata['source_wavecar_sha256'], spin_basis=BASIS,
        gauge=so.GAUGE, units=so.UNITS.copy(), producer='Independent exact enlarged-Hilbert-space fixture',
        source_sha256={'fixture': 'b'*64})


class SpinOperatorTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.path=Path(self.temp.name); self.wave=self.path/'WAVECAR'
        self.coeff=wavecar_fixture(self.wave, [[.4+.1j, .2+.3j], [.15-.2j, .6+.05j]])

    def load(self, **kwargs):
        return so.extract_wavecar_spin(self.wave, spin_basis=BASIS, **kwargs)

    def augment(self, source, ds, do, **changes):
        kw=dict(metadata=aug_metadata(source), kpoints_fractional=source.kpoints_fractional,
            energies_eV=source.energies_eV, lattice_A=source.lattice_A, band_indices=source.band_indices)
        kw.update(changes)
        return so.add_paw_augmentation(source, ds, do, **kw)

    def exact_fixture(self):
        wavecar_fixture(self.wave, [[np.sqrt(.7), 0], [0, np.sqrt(.8)]])
        source=self.load()
        c=so.Wavecar(self.wave, spinor_components=2).coefficients(0,[1,2]).astype(complex)
        p=np.zeros((1,2,2,1),complex)
        p[0,0,0,0]=np.sqrt(1-source.overlap[0,0,0].real)
        p[0,1,1,0]=np.sqrt(1-source.overlap[0,1,1].real)
        ds,do=so.projector_spin_correction(p,np.ones((1,1)))
        # Actual orthonormal physical vectors in a larger Hilbert space:
        # spin-major ordering, G orbital followed by augmentation orbital.
        physical=np.concatenate((c,p[0]),axis=2).reshape(2,4)
        oracle=np.array([physical.conj()@np.kron(s,np.eye(2))@physical.T for s in PAULI])
        return source,ds,do,oracle

    def test_complex_pauli_direct_kronecker_oracle_without_normalization(self):
        rng=np.random.default_rng(942)
        c=rng.normal(size=(4,2,7))+1j*rng.normal(size=(4,2,7))
        spin,overlap=so.pauli_from_coefficients(c)
        flat=c.reshape(4,14)
        oracle=np.array([flat.conj()@np.kron(s,np.eye(7))@flat.T for s in PAULI])
        np.testing.assert_allclose(spin,oracle,atol=1e-13)
        np.testing.assert_allclose(overlap,flat.conj()@flat.T,atol=1e-13)
        self.assertGreater(abs(overlap[0,0]-1),1)
        self.assertGreater(abs(spin[1,0,1].imag),.01)

    def test_actual_wavecar_record_order_phase_and_band_selection(self):
        data=self.load(bands=[2,1]); c=self.coeff[[1,0]]
        for a,s in enumerate(PAULI):
            np.testing.assert_allclose(data.spin_pauli[0,a],c.conj()@s@c.T,atol=1e-15)
        np.testing.assert_allclose(data.overlap[0],c.conj()@c.T,atol=1e-15)
        np.testing.assert_array_equal(data.energies_eV,[[1.,0.]])
        np.testing.assert_array_equal(data.band_indices,[2,1])
        self.assertFalse(data.metadata['physical_paw_validated'])
        self.assertFalse(data.metadata['normalization_applied'])

    def test_band_gauge_covariance(self):
        rng=np.random.default_rng(821)
        c=rng.normal(size=(3,2,5))+1j*rng.normal(size=(3,2,5))
        u,_=np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))
        before,o=so.pauli_from_coefficients(c)
        after,on=so.pauli_from_coefficients(np.einsum('mn,msg->nsg',u,c))
        np.testing.assert_allclose(after,np.array([u.conj().T@s@u for s in before]),atol=1e-13)
        np.testing.assert_allclose(on,u.conj().T@o@u,atol=1e-13)

    def test_projector_general_complex_metric_independent_kron_oracle(self):
        rng=np.random.default_rng(6)
        p=rng.normal(size=(2,3,2,4))+1j*rng.normal(size=(2,3,2,4))
        z=rng.normal(size=(2,4,4))+1j*rng.normal(size=(2,4,4));q=z+z.conj().swapaxes(-1,-2)
        ds,do=so.projector_spin_correction(p,q)
        for k in range(2):
            flat=p[k].reshape(3,8)
            for a,s in enumerate(PAULI):np.testing.assert_allclose(ds[k,a],flat.conj()@np.kron(s,q[k])@flat.T,atol=2e-13)
            np.testing.assert_allclose(do[k],flat.conj()@np.kron(np.eye(2),q[k])@flat.T,atol=2e-13)

    def test_exact_physical_metric_restored_with_no_renormalization(self):
        source,ds,do,oracle=self.exact_fixture()
        full=self.augment(source,ds,do)
        np.testing.assert_allclose(full.overlap[0],np.eye(2),atol=1e-15)
        np.testing.assert_allclose(full.spin_pauli[0],oracle,atol=1e-15)
        np.testing.assert_array_equal(full.pseudo_overlap,source.overlap)
        self.assertTrue(full.metadata['physical_paw_validated'])
        self.assertLess(full.spin_pauli[0,0,0,1].real,1.)
        with self.assertRaisesRegex(ValueError,'once'):self.augment(full,ds,do)

    def test_physical_metric_defect_is_rejected_not_repaired(self):
        source,ds,do,_=self.exact_fixture()
        do[0,0,0]+=.01
        with self.assertRaisesRegex(ValueError,'not identity'):self.augment(source,ds,do)

    def test_hermitian_but_unphysical_augmented_spin_is_rejected(self):
        source,ds,do,_=self.exact_fixture()
        ds[0,0,0,0]+=2
        with self.assertRaisesRegex(ValueError,'spin bound'):self.augment(source,ds,do)

    def test_augmentation_association_and_completion_guards(self):
        source,ds,do,_=self.exact_fixture()
        for key,value in [('source_wavecar_sha256','c'*64),('spin_basis','other'),('gauge','different'),
                          ('complete',False),('producer_status','TIMEOUT'),('source_sha256',{})]:
            with self.subTest(key=key):
                meta=aug_metadata(source);meta[key]=value
                with self.assertRaises(ValueError):self.augment(source,ds,do,metadata=meta)
        for key,value in [('kpoints_fractional',source.kpoints_fractional+.01),
                          ('energies_eV',source.energies_eV+.01),('lattice_A',source.lattice_A*1.1),
                          ('band_indices',source.band_indices[::-1])]:
            with self.subTest(key=key),self.assertRaisesRegex(ValueError,'mismatch'):
                self.augment(source,ds,do,**{key:value})

    def test_nonhermitian_or_incomplete_inputs_rejected(self):
        source,ds,do,_=self.exact_fixture()
        broken=ds.copy();broken[0,1,0,1]+=.02j
        with self.assertRaisesRegex(ValueError,'Hermitian'):self.augment(source,broken,do)
        with self.assertRaises(ValueError):self.augment(source,ds[:,:,:1,:],do)
        with self.assertRaises(ValueError):so.projector_spin_correction(np.ones((1,2,2,3)),np.ones((1,2,3,3)))
        bad=ds.copy();bad[0,0,0,0]=np.nan
        with self.assertRaises(ValueError):self.augment(source,bad,do)

    def test_extraction_preflight_invalid_spin_basis_and_bands(self):
        with self.assertRaisesRegex(ValueError,'spin_basis'):so.extract_wavecar_spin(self.wave,spin_basis='')
        for bands in [[],[1,1],[0],[3],[True],[1.5]]:
            with self.subTest(bands=bands),self.assertRaises(ValueError):self.load(bands=bands)

    def test_cache_preserves_raw_and_full_matrices_and_guards_hash(self):
        source,ds,do,_=self.exact_fixture();full=self.augment(source,ds,do)
        directory=self.path/'cache';so.write_spin(directory,full);restored=so.read_spin(directory)
        for name in so.ARRAYS:np.testing.assert_array_equal(getattr(full,name),getattr(restored,name))
        with self.assertRaisesRegex(ValueError,'exists'):so.write_spin(directory,full)
        p=directory/'spin.npz';p.write_bytes(p.read_bytes()+b'tamper')
        with self.assertRaisesRegex(ValueError,'checksum'):so.read_spin(directory)

    def test_delta_file_import_preserves_source_and_rejects_extra_arrays(self):
        source,ds,do,_=self.exact_fixture();path=self.path/'delta.npz';mp=self.path/'delta.json'
        arrays=dict(delta_spin_pauli=ds,delta_overlap=do,
            **{name:getattr(source,name) for name in ['kpoints_fractional','energies_eV','lattice_A','band_indices']})
        np.savez(path,**arrays);meta=aug_metadata(source);meta['data_npz_sha256']=so.sha256(path)
        mp.write_text(json.dumps(meta));full=so.import_paw_augmentation(source,path,mp)
        np.testing.assert_allclose(full.overlap[0],np.eye(2),atol=1e-15)
        self.assertEqual(full.metadata['augmentation']['imported_files']['npz'],'delta.npz')
        np.savez(path,**arrays,extra=np.zeros(1));meta['data_npz_sha256']=so.sha256(path);mp.write_text(json.dumps(meta))
        with self.assertRaisesRegex(ValueError,'archive'):so.import_paw_augmentation(source,path,mp)

    def test_raw_cache_cannot_claim_physical_paw_or_hide_renormalization(self):
        source=self.load()
        for key,value in [('physical_paw_validated',True),('normalization_applied',True)]:
            with self.subTest(key=key),self.assertRaises(ValueError):
                so.validate_spin(replace(source,metadata=dict(source.metadata,**{key:value})))
        changed=source.spin_pauli.copy();changed*=2
        with self.assertRaises(ValueError):so.validate_spin(replace(source,spin_pauli=changed))


if __name__=='__main__':unittest.main()
