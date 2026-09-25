"""Independent occupation-pair physics, degeneracy, geometry and output checks."""
import copy
import csv
from dataclasses import replace
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
import kubo_pairs as kp
from berry_data import CONDUCTANCE_QUANTUM_S, write_hall
from vaspberry_transport import fermi_dirac


def pair_data(energies, numerators, *, mesh=2):
    energies=np.array(energies,dtype=np.float64)
    if energies.ndim==1:
        energies=np.tile(energies,(mesh*mesh,1))
    nk,nb=energies.shape
    x,y=np.meshgrid(np.arange(mesh)/mesh,np.arange(mesh)/mesh,indexing='ij')
    q=np.column_stack((x.ravel(),y.ravel(),np.zeros(nk)))
    n,m=np.triu_indices(nb,1)
    values=np.array(numerators,dtype=np.float64)
    if values.ndim==2:
        values=np.tile(values,(nk,1,1))
    meta=dict(schema=kp.SCHEMA,version=1,complete=True,source_nbands=nb,
              normalization='numerator=-2Im(D_a_nm*D_b_mn)',components=['yz','zx','xy'],
              units={'energy':'eV','numerator':'eV^2 Angstrom^2','lattice':'Angstrom','reciprocal':'1/Angstrom'},
              sampling={'kind':'uniform_full_2d','mesh':[mesh,mesh],'plane_axes':[0,1]},
              spin_multiplicity=1,source_nspin=1,spinor_components=1,spin_channel_1based=1,
              energy_reference='analytic test zero',source_operator={'kind':'analytic_test','accuracy_status':'validated'})
    return kp.PairData(meta,np.arange(1,nk+1),np.arange(1,nb+1),n,m,q,np.full(nk,1/nk),
                       energies,values,np.eye(3),2*np.pi*np.eye(3))


def from_vertices(energies,vertices):
    n,m=np.triu_indices(vertices.shape[-1],1)
    values=np.stack([-2*np.imag(vertices[:,a,n,m]*vertices[:,b,m,n])
                     for a,b in ((1,2),(2,0),(0,1))],axis=-1)
    return pair_data(energies,values)


def random_vertices(nb=3):
    rng=np.random.default_rng(85293)
    x=rng.normal(size=(4,3,nb,nb))+1j*rng.normal(size=(4,3,nb,nb))
    return (x+x.conj().swapaxes(-1,-2))/2


def qwz(mesh=16):
    x,y=np.meshgrid(np.arange(mesh)/mesh*2*np.pi,np.arange(mesh)/mesh*2*np.pi,indexing='ij')
    x,y=x.ravel(),y.ravel()
    d=np.column_stack((np.sin(x),np.sin(y),-1+np.cos(x)+np.cos(y)))
    dx=np.column_stack((np.cos(x),np.zeros_like(x),-np.sin(x)))
    dy=np.column_stack((np.zeros_like(y),np.cos(y),-np.sin(y)))
    norm=np.linalg.norm(d,axis=1)
    omega=np.einsum('ki,ki->k',d,np.cross(dx,dy))/(2*norm**3)
    numerator=np.zeros((mesh*mesh,1,3));numerator[:,0,2]=omega*(2*norm)**2
    return pair_data(np.column_stack((-norm,norm)),numerator,mesh=mesh),omega


def total(rows):
    return [r for r in rows if r['region']=='total']


class PairPhysicsTests(unittest.TestCase):
    def test_pair_band_limit_preserves_source_and_separates_virtual_window(self):
        vertices=random_vertices(4); energies=np.tile([-.5,.2,1.5,3.],(4,1))
        data=from_vertices(energies,vertices)
        before=data.energies_eV.copy()
        rows,meta=kp.pair_hall_spectrum(data,[-.1,0.],[0.,300.],mu_reference=0.,pair_band_max=3)
        for row in total(rows):
            f=fermi_dirac(energies[:,:3],row['mu_eV'],row['temperature_K'])
            expected=np.zeros(4)
            for n in range(3):
                for m in range(3):
                    if n!=m:
                        expected += f[:,n]*(-2*np.imag(vertices[:,0,n,m]*vertices[:,1,m,n]))/(energies[:,n]-energies[:,m])**2
            self.assertAlmostEqual(row['sigma_e2_over_h'],-2*np.pi*expected.mean(),places=12)
        self.assertEqual(meta['source_nbands'],4)
        self.assertEqual(meta['pair_band_window'],[1,3])
        self.assertEqual(meta['scope'],'truncated_pair_space')
        np.testing.assert_array_equal(data.energies_eV,before)
        full,_=kp.pair_hall_spectrum(data,[0.],[0.],mu_reference=0.)
        explicit,_=kp.pair_hall_spectrum(data,[0.],[0.],mu_reference=0.,pair_band_max=4)
        self.assertEqual(full,explicit)
        for bad in [1,5,2.5,True]:
            with self.assertRaises(ValueError):
                kp.pair_hall_spectrum(data,[0.],[0.],mu_reference=0.,pair_band_max=bad)
        with self.assertRaisesRegex(ValueError,'highest selected'):
            kp.pair_hall_spectrum(data,[2.],[0.],mu_reference=0.,pair_band_max=3)
        degenerate=replace(data,energies_eV=np.tile([-.5,.2,1.5,1.5],(4,1)))
        with self.assertRaisesRegex(ValueError,'cutoff splits'):
            kp.pair_hall_spectrum(degenerate,[0.],[0.],mu_reference=0.,pair_band_max=3)
        # Averaging a bounded near-degenerate group must not disguise a
        # raw unresolved separation across the virtual-state cutoff.
        clustered=replace(data,energies_eV=np.tile([-.5,1.,1.+.75e-7,1.+1.5e-7],(4,1)))
        with self.assertRaisesRegex(ValueError,'cutoff splits'):
            kp.pair_hall_spectrum(clustered,[0.],[0.],mu_reference=0.,pair_band_max=3,
                                  degeneracy_policy='coalesce',degeneracy_threshold_eV=1e-7)

    def test_two_band_qwz_chern_sign_factor_and_finite_temperature(self):
        data,omega=qwz(32)
        rows,_=kp.pair_hall_spectrum(data,[.2,-.2,0.],[0.,300.],mu_reference=0.,mu_chunk=2)
        for row in total(rows):
            f=fermi_dirac(data.energies_eV,row['mu_eV'],row['temperature_K'])
            expected=-2*np.pi*np.mean((f[:,0]-f[:,1])*omega)
            self.assertAlmostEqual(row['sigma_e2_over_h'],expected,places=12)
            self.assertAlmostEqual(row['sigma_S'],expected*CONDUCTANCE_QUANTUM_S,places=16)
        self.assertAlmostEqual(total(rows)[0]['sigma_e2_over_h'],-1.,delta=1e-8)
        self.assertEqual([r['mu_eV'] for r in total(rows)[:3]],[.2,-.2,0.])

    def test_finite_T_pair_sum_equals_independent_resolved_band_sum(self):
        vertices=random_vertices();energies=np.tile([-.3,.2,2.],(4,1))
        energies[:,:2]+=np.arange(4)[:,None]*.01
        data=from_vertices(energies,vertices)
        omega=np.zeros((4,3))
        for k in range(4):
            for n in range(3):
                for m in range(3):
                    if n!=m:
                        omega[k,n]+=-2*np.imag(vertices[k,0,n,m]*vertices[k,1,m,n])/(energies[k,n]-energies[k,m])**2
        rows,_=kp.pair_hall_spectrum(data,[.4,-.4,0.],[0.,150.,400.],mu_reference=.05,mu_chunk=1)
        for row in total(rows):
            occ=fermi_dirac(energies,row['mu_eV'],row['temperature_K'])
            ref=fermi_dirac(energies,.05,row['temperature_K'])
            self.assertAlmostEqual(row['sigma_e2_over_h'],-2*np.pi*np.mean((occ*omega).sum(axis=1)),places=12)
            self.assertAlmostEqual(row['delta_sigma_e2_over_h'],-2*np.pi*np.mean(((occ-ref)*omega).sum(axis=1)),places=12)
            self.assertAlmostEqual(row['electrons_per_cell'],np.mean(occ.sum(axis=1)),places=13)
        larger,_=kp.pair_hall_spectrum(data,[.4,-.4,0.],[0.,150.,400.],mu_reference=.05,mu_chunk=10)
        np.testing.assert_allclose([r['sigma_e2_over_h'] for r in rows],[r['sigma_e2_over_h'] for r in larger],atol=1e-12)

    def test_exact_degenerate_subspace_rotation_leaves_response_invariant(self):
        vertices=random_vertices();data=from_vertices([-1.,-1.,2.],vertices)
        rng=np.random.default_rng(525)
        rotation=np.eye(3,dtype=complex)
        rotation[:2,:2]=np.linalg.qr(rng.normal(size=(2,2))+1j*rng.normal(size=(2,2)))[0]
        rotated=rotation.conj().T@vertices@rotation
        second=from_vertices([-1.,-1.,2.],rotated)
        first_rows,_=kp.pair_hall_spectrum(data,[0.,-.8],[0.,300.],mu_reference=0.)
        second_rows,_=kp.pair_hall_spectrum(second,[0.,-.8],[0.,300.],mu_reference=0.)
        self.assertTrue(np.any(abs(data.numerator_eV2_A2[:,0,2])>1e-3))
        np.testing.assert_allclose([r['sigma_e2_over_h'] for r in first_rows],
                                   [r['sigma_e2_over_h'] for r in second_rows],atol=1e-12)

    def test_unresolved_unequal_occupations_reject_and_coalescing_is_explicit(self):
        data=from_vertices([-1e-8,1e-8,2.],random_vertices());original=data.energies_eV.copy()
        for temperature in (0.,300.):
            with self.subTest(temperature=temperature),self.assertRaisesRegex(ValueError,'unequal occupation'):
                kp.pair_hall_spectrum(data,[0.],[temperature],mu_reference=.01,degeneracy_threshold_eV=1e-7)
        with self.assertRaisesRegex(ValueError,'cuts unresolved'):
            kp.pair_hall_spectrum(data,[0.],[0.],mu_reference=.01,degeneracy_policy='coalesce')
        rows,meta=kp.pair_hall_spectrum(data,[0.],[300.],mu_reference=.01,degeneracy_policy='coalesce')
        self.assertEqual(meta['degeneracy']['policy'],'coalesce')
        self.assertEqual(meta['degeneracy']['coalesced_groups'],4)
        self.assertAlmostEqual(meta['degeneracy']['max_energy_shift_eV'],1e-8)
        self.assertGreater(meta['degeneracy']['max_occupation_shift'],0.)
        np.testing.assert_array_equal(data.energies_eV,original)
        exact=replace(data,energies_eV=np.tile([0.,0.,2.],(4,1)))
        expected,_=kp.pair_hall_spectrum(exact,[0.],[300.],mu_reference=.01)
        self.assertAlmostEqual(total(rows)[0]['sigma_e2_over_h'],total(expected)[0]['sigma_e2_over_h'],places=13)
        effective,groups=kp.coalesce_energies(np.array([0.,.75e-7,1.5e-7]),1e-7)
        self.assertEqual(len(groups),1)
        self.assertLessEqual(groups[0][3]-groups[0][2],1e-7)
        self.assertEqual(effective[-1],1.5e-7)

    def test_regions_orientation_multiplicity_and_delta_survive_filled_baseline(self):
        data=pair_data([-10.,.25,2.],[[0,0,0],[0,0,1e20],[0,0,1.]])
        spec={'regions':[{'name':'left','k_ids':[1]},{'name':'right','k_ids':[2]}]}
        rows,_=kp.pair_hall_spectrum(data,[.5],[0.,300.],mu_reference=0.,region_spec=spec,
                                   differences=[('contrast','left','right')])
        for t in (0.,300.):
            selected={r['region']:r for r in rows if r['temperature_K']==t}
            self.assertEqual(selected['contrast']['sigma_e2_over_h'],
                             selected['left']['sigma_e2_over_h']-selected['right']['sigma_e2_over_h'])
            self.assertEqual(selected['total']['sigma_e2_over_h'],
                             sum(selected[r]['sigma_e2_over_h'] for r in ('left','right','rest')))
            delta_f=(fermi_dirac(np.array([.25,2.]),.5,t)-fermi_dirac(np.array([.25,2.]),0.,t))
            expected=-2*np.pi*(delta_f[0]-delta_f[1])/1.75**2
            # At finite T the huge lower-to-upper pair has a tiny, physically
            # real upper-band occupation change; isolate the exact T=0 case.
            if t==0:
                self.assertAlmostEqual(selected['total']['delta_sigma_e2_over_h'],expected,places=13)
        regular,_=qwz(8)
        before,_=kp.pair_hall_spectrum(regular,[0.],[0.],mu_reference=0.)
        rot=np.array([[0.,1,0],[0,0,1],[1,0,0]])
        rotated=replace(regular,lattice_A=regular.lattice_A@rot,reciprocal_inv_A=regular.reciprocal_inv_A@rot,
                        numerator_eV2_A2=regular.numerator_eV2_A2@rot)
        after,_=kp.pair_hall_spectrum(rotated,[0.],[0.],mu_reference=0.)
        self.assertAlmostEqual(total(before)[0]['sigma_e2_over_h'],total(after)[0]['sigma_e2_over_h'],places=13)
        meta=copy.deepcopy(rotated.metadata);meta['sampling']['plane_axes']=[1,0]
        reverse,_=kp.pair_hall_spectrum(replace(rotated,metadata=meta),[0.],[0.],mu_reference=0.)
        self.assertAlmostEqual(total(before)[0]['sigma_e2_over_h'],-total(reverse)[0]['sigma_e2_over_h'],places=13)
        doubled=replace(regular,metadata=dict(regular.metadata,spin_multiplicity=2))
        result,_=kp.pair_hall_spectrum(doubled,[0.],[0.],mu_reference=0.)
        self.assertAlmostEqual(total(result)[0]['sigma_e2_over_h'],2*total(before)[0]['sigma_e2_over_h'],places=13)

    def test_invalid_sampling_occupation_window_and_spin_count_reject(self):
        data,_=qwz(4)
        with self.assertRaisesRegex(ValueError,'highest source band occupied'):
            kp.pair_hall_spectrum(data,[4.],[0.],mu_reference=0.)
        _,meta=kp.pair_hall_spectrum(data,[4.],[0.],mu_reference=0.,allow_partial_bands=True)
        self.assertEqual(meta['scope'],'partial_band_contribution')
        bad=replace(data,metadata=dict(data.metadata,sampling={'kind':'points'}))
        with self.assertRaisesRegex(ValueError,'full uniform'):
            kp.pair_hall_spectrum(bad,[0.],[0.],mu_reference=0.)
        q=data.kpoints_fractional.copy();q[0]=q[1]
        with self.assertRaisesRegex(ValueError,'duplicate'):
            kp.validate_pairs(replace(data,kpoints_fractional=q))
        with self.assertRaisesRegex(ValueError,'multiplicity one'):
            kp.validate_pairs(replace(data,metadata=dict(data.metadata,spinor_components=2,spin_multiplicity=2)))
        for mus,ts in (([0.,0.],[0.]),([0.],[-1.]),([np.nan],[0.])):
            with self.subTest(mus=mus,ts=ts),self.assertRaises(ValueError):
                kp.pair_hall_spectrum(data,mus,ts,mu_reference=0.)

    def test_pair_cache_roundtrip_and_all_hall_formats_have_identical_rows(self):
        data,_=qwz(4)
        with tempfile.TemporaryDirectory() as directory:
            base=Path(directory)
            kp.write_pairs(base/'pairs',data)
            loaded=kp.read_pairs(base/'pairs')
            for key in kp.ARRAYS:
                np.testing.assert_array_equal(getattr(data,key),getattr(loaded,key))
            rows,meta=kp.pair_hall_spectrum(loaded,[.1,0.],[0.,300.],mu_reference=0.)
            written=write_hall(base/'hall',rows,meta,formats=('csv','dat','npz'))
            self.assertEqual(written['output_formats'],['csv','dat','npz'])
            with (base/'hall/conductivity.csv').open() as f:
                csv_rows=list(csv.DictReader(f))
            dat=(base/'hall/conductivity.dat').read_text()
            self.assertTrue(dat.startswith('# '))
            dat_rows=list(csv.DictReader(io.StringIO(dat[2:]),delimiter='\t'))
            self.assertEqual(csv_rows,dat_rows)
            with np.load(base/'hall/conductivity.npz',allow_pickle=False) as arrays:
                self.assertEqual(set(arrays),set(rows[0]))
                for key in rows[0]:
                    np.testing.assert_array_equal(arrays[key],[r[key] for r in rows])
            sidecar=json.loads((base/'hall/conductivity.json').read_text())
            self.assertEqual(set(sidecar['output_sha256']),{'conductivity.csv','conductivity.dat','conductivity.npz'})
            write_hall(base/'dat-only',rows,meta,formats=('dat',))
            self.assertEqual({p.name for p in (base/'dat-only').iterdir()},{'conductivity.dat','conductivity.json'})
            with self.assertRaisesRegex(ValueError,'exists'):
                kp.write_pairs(base/'pairs',data)
            with (base/'pairs/pairs.npz').open('ab') as f:f.write(b'damaged')
            with self.assertRaisesRegex(ValueError,'checksum'):
                kp.read_pairs(base/'pairs')

class NativePairContractTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.data=from_vertices([-1.,-1.,2.],random_vertices())
        self.wave=self.root/'WAVECAR';self.wave.write_bytes(b'synthetic parser-contract fixture')
        self.fake=SimpleNamespace(energies=self.data.energies_eV,kpoints=self.data.kpoints_fractional,
            header=SimpleNamespace(ispin=1,lattice=self.data.lattice_A,reciprocal=self.data.reciprocal_inv_A),
            coefficients=lambda k,b:np.ones((1,1,1)))
        self.options=dict(spin=1,spinor_components=1,spin_multiplicity=1,
                          sampling=self.data.metadata['sampling'],energy_reference='fixture zero')

    def native(self,name='pairs.csv',overrides=None,row_transform=None,footer=True):
        d=self.data;nk,nb=d.energies_eV.shape
        metadata=dict(schema=kp.NATIVE_SCHEMA,result_kind='UNORDERED_INTERBAND_NUMERATORS',
            normalization='STANDARD_MINUS_TWO_IM',operator=kp.OPERATOR,berry_connection=kp.CONNECTION,
            occupation_weighting='NONE',denominator_weighting='NONE',pair_order='n_lt_m',
            components='yz,zx,xy',numerator_units='eV^2*Angstrom^2',gap_definition='ABS_EN_MINUS_EM',index_base='1',
            source_nbands=nb,source_nkpoints=nk,source_nspin=1,spinor_components=1,
            pairs_per_k=len(d.pair_n),expected_rows=nk*len(d.pair_n),reciprocal_convention='2pi')
        for i in range(3):
            metadata[f'lattice_A_{i+1}']=','.join(map(str,d.lattice_A[i]))
            metadata[f'reciprocal_inv_A_{i+1}']=','.join(map(str,d.reciprocal_inv_A[i]))
        metadata.update(overrides or {})
        rows=[]
        for k in range(nk):
            for p,(n,m) in enumerate(zip(d.pair_n,d.pair_m)):
                rows.append(dict(spin=1,k_index=k+1,n_band=n+1,m_band=m+1,
                    **dict(zip(('kx_frac','ky_frac','kz_frac'),d.kpoints_fractional[k])),
                    energy_n_eV=d.energies_eV[k,n],energy_m_eV=d.energies_eV[k,m],
                    **dict(zip((f'numerator_{a}_eV2_A2' for a in ('yz','zx','xy')),d.numerator_eV2_A2[k,p])),
                    gap_eV=abs(d.energies_eV[k,n]-d.energies_eV[k,m])))
        if row_transform:rows=row_transform(rows)
        path=self.root/name
        with path.open('w',newline='') as f:
            for key,value in metadata.items():f.write(f'# {key}={value}\n')
            writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
            if footer:f.write('# result_status=PASS\n')
        return path

    def load(self,path):
        with patch.object(kp,'Wavecar',return_value=self.fake):
            return kp.import_native_pairs(path,self.wave,**self.options)

    def test_native_import_requires_matching_schema_completion_and_coordinates(self):
        loaded=self.load(self.native())
        np.testing.assert_array_equal(loaded.numerator_eV2_A2,self.data.numerator_eV2_A2)
        np.testing.assert_array_equal(loaded.energies_eV,self.data.energies_eV)
        with self.assertRaisesRegex(ValueError,'complete native'):
            self.load(self.native('no-footer.csv',footer=False))
        for key,value in [('components','xy,yz,zx'),('numerator_units','eV*Angstrom'),
                          ('result_kind','BAND_CURVATURE'),('gap_definition','SIGNED'),('expected_rows',13)]:
            with self.subTest(key=key),self.assertRaises(ValueError):
                self.load(self.native('bad-meta.csv',{key:value}))
        for transform in (lambda rows:rows+rows[:1],lambda rows:rows[:-1],
                          lambda rows:[dict(rows[0],energy_n_eV=-.9)]+rows[1:],
                          lambda rows:[dict(rows[0],kx_frac=.1)]+rows[1:],
                          lambda rows:[dict(rows[0],gap_eV=1.)]+rows[1:]):
            with self.subTest(transform=transform),self.assertRaises(ValueError):
                self.load(self.native('bad-row.csv',row_transform=transform))

    def test_integer_indices_and_pair_coverage_are_required(self):
        for field in ('k_ids','band_ids','pair_n','pair_m'):
            with self.subTest(field=field),self.assertRaisesRegex(ValueError,'integer'):
                kp.validate_pairs(replace(self.data,**{field:getattr(self.data,field).astype(float)}))
        with self.assertRaisesRegex(ValueError,'complete unique'):
            kp.validate_pairs(replace(self.data,pair_m=self.data.pair_m[::-1].copy()))

    def test_chunked_import_rejects_cross_chunk_duplicate_and_truncation(self):
        # 4 k points * C(129,2) = 33,024 rows spans the actual 32,768-row
        # parsing boundary. The second chunk must check prior-chunk coverage.
        self.data=pair_data(np.arange(129,dtype=float),np.zeros((129*128//2,3)))
        self.fake.energies=self.data.energies_eV
        loaded=self.load(self.native('large.csv'))
        np.testing.assert_array_equal(loaded.energies_eV,self.data.energies_eV)
        np.testing.assert_array_equal(loaded.numerator_eV2_A2,self.data.numerator_eV2_A2)
        with self.assertRaisesRegex(ValueError,'duplicate'):
            self.load(self.native('duplicate-across-chunks.csv',
                row_transform=lambda rows:rows[:32768]+rows[:1]+rows[32769:]))
        with self.assertRaisesRegex(ValueError,'total row count'):
            self.load(self.native('truncated-last-chunk.csv',row_transform=lambda rows:rows[:-1]))
        path=self.native('wrong-columns.csv')
        path.write_text(path.read_text().replace('numerator_yz_eV2_A2','numerator_yz_wrong'))
        with self.assertRaisesRegex(ValueError,'column names/order'):
            self.load(path)

    def test_fixed_bundle_is_only_a_constant_T0_global_gap_result(self):
        path=self.root/'bundle.csv'
        metadata=dict(schema='VASPBERRY_BARE_MOMENTUM_KUBO_BUNDLE_V1',normalization='STANDARD_MINUS_TWO_IM',
            operator=kp.OPERATOR,berry_connection=kp.CONNECTION,result_status='PASS',band_min=1,band_max=2,
            band_rank=2,source_nbands=3,intermediate_bands='EXTERNAL_TO_SELECTED_BUNDLE_WITHIN_SOURCE_NBANDS')
        omega=np.array([.2,.4,.6,.8])
        with path.open('w',newline='') as f:
            for k,v in metadata.items():f.write(f'# {k}={v}\n')
            fields=['spin','k_index','kx_frac','ky_frac','kz_frac','omega_z_A2','min_external_gap_eV']
            writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
            for k in range(4):
                writer.writerow(dict(spin=1,k_index=k+1,**dict(zip(fields[2:5],self.data.kpoints_fractional[k])),
                                     omega_z_A2=omega[k],min_external_gap_eV=3.))
        with patch.object(kp,'Wavecar',return_value=self.fake):
            rows,meta=kp.bundle_hall_spectrum(path,self.wave,[1.,-.5,0.],occupied=2,mu_reference=.25,**self.options)
            self.assertEqual([r['mu_eV'] for r in total(rows)],[1.,-.5,0.])
            for row in total(rows):
                self.assertAlmostEqual(row['sigma_e2_over_h'],-2*np.pi*np.mean(omega))
                self.assertEqual(row['delta_sigma_e2_over_h'],0.)
                self.assertEqual(row['electrons_per_cell'],2.)
            self.assertEqual(meta['global_gap_eV'],3.)
            for mu in (-1.,2.,3.):
                with self.subTest(mu=mu),self.assertRaisesRegex(ValueError,'global insulating gap'):
                    kp.bundle_hall_spectrum(path,self.wave,[mu],occupied=2,mu_reference=.25,**self.options)
