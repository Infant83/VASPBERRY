"""Full CLI/table checks against an analytic spin-conserving QSH fixture.

The fixture exercises a producer contract; it is not a VASP/PAW material
dataset and its explicitly labelled provenance must never be published as one.
"""
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import csv
import json
import os
import subprocess
import struct
import sys
import tempfile
import unittest

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
import spin_hall_workflow as workflow
from exported_matrix_kubo import sha256
from spin_hall import CONDUCTANCE_QUANTUM_S

X=np.array([[0,1],[1,0]],complex)
Y=np.array([[0,-1j],[1j,0]],complex)
Z=np.diag([1.,-1.]).astype(complex)


def fixture(directory):
    q=np.array([(i/4+.037,j/4+.021,0.) for i in range(4) for j in range(4)])
    e=[];v=[];s=[];oracle=[]
    for point in q:
        kx,ky=2*np.pi*point[:2]
        vector=np.array([np.sin(kx),np.sin(ky),-1+np.cos(kx)+np.cos(ky)])
        upper=vector[0]*X+vector[1]*Y+vector[2]*Z
        lower=-vector[0]*X+vector[1]*Y+vector[2]*Z
        h=np.zeros((4,4),complex);h[:2,:2]=upper;h[2:,2:]=lower
        d=np.zeros((3,4,4),complex)
        d[0,:2,:2]=np.cos(kx)*X-np.sin(kx)*Z
        d[1,:2,:2]=np.cos(ky)*Y-np.sin(ky)*Z
        d[0,2:,2:]=-np.cos(kx)*X-np.sin(kx)*Z
        d[1,2:,2:]=np.cos(ky)*Y-np.sin(ky)*Z
        energies,u=np.linalg.eigh(h)
        e.append(energies);v.append(u.conj().T@d@u)
        s.append(np.array([u.conj().T@np.kron(pauli,np.eye(2))@u for pauli in [X,Y,Z]]))
        # Analytic occupied spin curvature from two opposite Chern blocks.
        oracle.append((np.cos(kx)+np.cos(ky)-np.cos(kx)*np.cos(ky))/np.linalg.norm(vector)**3)
    arrays=dict(energies_eV=np.array(e),kpoints_fractional=q,weights=np.full(16,1/16),
        lattice_A=np.eye(3),spin_pauli=np.array(s),overlap=np.broadcast_to(np.eye(4,dtype=complex),(16,4,4)).copy(),
        velocity_eVA=np.array(v),band_indices=np.arange(1,5,dtype=np.int64))
    metadata=dict(schema=workflow.MATRIX_SCHEMA,version=1,complete=True,producer_status='PASS',
        source_wavecar_sha256='0'*64,operator_accuracy='paw_full_velocity',physical_paw_validated=True,
        spin_basis_id='cartesian',velocity_basis='cartesian',gauge='same_eigenvectors_as_source_wavecar',
        normalization_applied=False,full_source_band_coverage=True,diagonal_available=True,degenerate_blocks_available=True,
        units={'spin':'dimensionless Pauli','velocity':'eV angstrom = hbar*v','energies':'eV','lattice':'angstrom','overlap':'dimensionless'},
        spinor_components=2,spin_multiplicity=1,
        producer='Analytic QSH regression fixture; no VASP/PAW material result',
        operator_scope='Exact finite test space emulating the producer interface; no physical PAW claim',
        source_sha256={'WAVECAR':'0'*64,'analytic_fixture':'1'*64},nkpoints=16,source_nbands=4)
    matrix=directory/'matrices.npz';meta=directory/'matrices.json'
    save(matrix,meta,arrays,metadata)
    return arrays,metadata,np.array(oracle),matrix,meta


def save(path,meta_path,arrays,metadata):
    np.savez_compressed(path,**arrays)
    metadata=dict(metadata,data_npz_sha256=sha256(path))
    meta_path.write_text(json.dumps(metadata,indent=2)+'\n')


class SpinHallWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.path=Path(self.tmp.name)
        self.arrays,self.meta,self.oracle,self.matrix,self.metadata=fixture(self.path)

    def args(self,name='result',**changes):
        values=dict(matrices=self.matrix,metadata=self.metadata,output_dir=self.path/name,
            mesh=[4,4],plane_axes=[0,1],occupied=2,source_band_limit=None,mu_eV=None,
            gap_threshold_eV=1e-8,k_chunk=3,memory_limit_mib=128.,time_limit=30.,formats=['csv','dat','npz'])
        values.update(changes)
        return SimpleNamespace(**values)

    def run_cli(self,out='cli-result',extra=()):
        env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
        return subprocess.run([sys.executable,str(ROOT/'tools/vaspberry_kubo.py'),'spin-hall',
            '--matrices',str(self.matrix),'--metadata',str(self.metadata),'--mesh','4','4','--occupied','2',
            '--k-chunk','3','--formats','csv','dat','npz','--output-dir',str(self.path/out),*extra],
            env=env,capture_output=True,text=True,timeout=30)

    def test_cli_analytic_qsh_full_formats_units_and_no_coarse_rounding(self):
        p=self.run_cli();self.assertEqual(p.returncode,0,p.stderr+p.stdout)
        out=self.path/'cli-result';meta=json.loads((out/'spin_hall.json').read_text())
        run=json.loads((out/'run.json').read_text());self.assertEqual(run['status'],'PASS')
        self.assertEqual(run['completed_points'],16)
        with np.load(out/'spin_hall.npz',allow_pickle=False) as z:
            np.testing.assert_allclose(z['omega_spin_A2'][:,2,0,1],self.oracle,atol=3e-14)
            np.testing.assert_allclose(z['omega_spin_A2'][:,2,1,0],-self.oracle,atol=3e-14)
            expected=np.pi*np.mean(self.oracle)
            self.assertAlmostEqual(z['sigma_hbar_over_e_e2_over_h'][2,0,1],expected,places=13)
            self.assertGreater(abs(expected-1),1e-3) # four-by-four is deliberately coarse
            np.testing.assert_allclose(z['charge_sigma_e2_over_h'],0,atol=2e-14)
            a=list(csv.DictReader((out/'spin_hall.csv').open()))
            b=list(csv.DictReader((out/'spin_hall.dat').open(),delimiter=' '))
            self.assertEqual(a,b);self.assertEqual(len(a),27)
            for row in a:
                index=tuple('xyz'.index(row[key]) for key in ['spin','current','electric_field'])
                value=float(row['sigma_hbar_over_e_e2_over_h'])
                self.assertEqual(value,z['sigma_hbar_over_e_e2_over_h'][index])
                self.assertAlmostEqual(float(row['sigma_hbar_over_e_S']),value*CONDUCTANCE_QUANTUM_S,places=19)
            curvature=list(csv.DictReader((out/'spin_curvature.csv').open()))
            self.assertEqual(len(curvature),16*27)
            for row in curvature:
                index=(int(row['k_index'])-1,)+tuple('xyz'.index(row[key]) for key in ['spin','current','electric_field'])
                self.assertEqual(float(row['omega_spin_A2']),z['omega_spin_A2'][index])
        self.assertEqual(meta['spin_current_method'],'finite_band_projected_product')
        self.assertFalse(meta['integer_rounding_applied'])
        self.assertIn('analytic',meta['source_metadata']['operator_scope'].lower()+meta['source_metadata']['producer'].lower())
        for name,digest in meta['output_sha256'].items():self.assertEqual(sha256(out/name),digest)

    def test_k_chunk_invariance_and_repeat_output_guard(self):
        workflow.spin_hall_command(self.args('one',k_chunk=1))
        workflow.spin_hall_command(self.args('seven',k_chunk=7))
        with np.load(self.path/'one/spin_hall.npz') as a,np.load(self.path/'seven/spin_hall.npz') as b:
            for key in a.files:np.testing.assert_allclose(a[key],b[key],rtol=0,atol=3e-14)
        with self.assertRaisesRegex(ValueError,'output exists'):workflow.spin_hall_command(self.args('one'))

    def test_explicit_band_limit_keeps_original_source_count(self):
        a=dict(self.arrays)
        a['energies_eV']=np.concatenate((a['energies_eV'],np.tile([10.,11.],(16,1))),axis=1)
        for name in ['spin_pauli','velocity_eVA']:
            x=np.zeros((16,3,6,6),complex);x[:,:,:4,:4]=a[name];a[name]=x
        a['spin_pauli'][:,:,4:,4:]=np.array([X,Y,Z])[None]
        a['overlap']=np.broadcast_to(np.eye(6,dtype=complex),(16,6,6)).copy()
        a['band_indices']=np.arange(1,7,dtype=np.int64)
        save(self.matrix,self.metadata,a,dict(self.meta,source_nbands=6))
        result=workflow.spin_hall_command(self.args('capped',source_band_limit=4))
        self.assertEqual(result['source_band_count'],6)
        self.assertEqual(result['retained_band_count'],4)
        self.assertAlmostEqual(result['sigma_hbar_over_e_e2_over_h'][2][0][1],np.pi*np.mean(self.oracle),places=13)
        self.assertIn('finite-band',result['limitation'])

    def test_spin_matrix_cli_raw_subset_and_unpaired_augmentation_rejection(self):
        wave=self.path/'WAVECAR';stride=256;data=bytearray(5*stride)
        struct.pack_into('<3d',data,0,stride,1.,45200.)
        struct.pack_into('<12d',data,stride,1.,2.,1.,*np.eye(3).ravel())
        struct.pack_into('<10d',data,2*stride,2.,0.,0.,0.,-1.,0.,1.,2.,0.,0.)
        for band,row in enumerate([[.7,.1j],[.2+.1j,.8]]):
            c=np.array(row,dtype='<c8');data[(3+band)*stride:(3+band)*stride+c.nbytes]=c.tobytes()
        wave.write_bytes(data)
        argv=[sys.executable,str(ROOT/'tools/vaspberry_kubo.py'),'spin-matrix','--wavecar',str(wave),
              '--spin-basis','Analytic binary test spinor axes','--bands','2:2','--output-dir',str(self.path/'raw')]
        p=subprocess.run(argv,capture_output=True,text=True,timeout=30)
        self.assertEqual(p.returncode,0,p.stderr)
        meta=json.loads((self.path/'raw/spin.json').read_text())
        self.assertEqual(meta['operator_scope'],'raw_pseudo');self.assertFalse(meta['physical_paw_validated'])
        self.assertEqual(meta['source_nbands'],2)
        with np.load(self.path/'raw/spin.npz') as z:
            np.testing.assert_array_equal(z['band_indices'],[2]);self.assertEqual(z['energies_eV'][0,0],2.)
            self.assertNotEqual(z['overlap'][0,0,0],1.)
        argv[-1]=str(self.path/'unpaired');argv+=['--augmentation',str(self.path/'missing.npz')]
        p=subprocess.run(argv,capture_output=True,text=True,timeout=30)
        self.assertNotEqual(p.returncode,0);self.assertIn('supplied together',p.stderr)
        self.assertFalse((self.path/'unpaired').exists())

    def test_invalid_source_coverage_gauge_units_and_producer_rejected(self):
        variants=[('complete',False),('producer_status','TIMEOUT'),('physical_paw_validated',False),
            ('operator_accuracy','canonical_pseudo'),('gauge','unassociated'),('spin_basis_id','unknown'),
            ('normalization_applied',True),('full_source_band_coverage',False),('diagonal_available',False),
            ('degenerate_blocks_available',False),('spinor_components',1),('spin_multiplicity',2),
            ('source_sha256',{'WAVECAR':'f'*64}),('units',dict(self.meta['units'],spin='hbar/2'))]
        for i,(key,value) in enumerate(variants):
            with self.subTest(key=key):
                save(self.matrix,self.metadata,self.arrays,dict(self.meta,**{key:value}))
                with self.assertRaises(ValueError):workflow.spin_hall_command(self.args(f'bad-{i}'))
                self.assertFalse((self.path/f'bad-{i}').exists())
                self.assertFalse((self.path/f'bad-{i}.partial').exists())

    def test_matrix_geometry_sampling_and_metric_rejections(self):
        changes=[]
        weights=self.arrays['weights'].copy();weights[:2]=[.08,.045];changes.append(('weights',weights))
        points=self.arrays['kpoints_fractional'].copy();points[1]=points[0];changes.append(('kpoints_fractional',points))
        bands=self.arrays['band_indices'][::-1];changes.append(('band_indices',bands))
        overlap=self.arrays['overlap'].copy();overlap[0,0,0]=.9;changes.append(('overlap',overlap))
        velocity=self.arrays['velocity_eVA'].copy();velocity[0,0,0,1]+=.1j;changes.append(('velocity_eVA',velocity))
        spin=self.arrays['spin_pauli'].copy();spin[0,0]=2*np.eye(4);changes.append(('spin_pauli',spin))
        skew=self.arrays['lattice_A'].copy();skew[2,0]=.4;changes.append(('lattice_A',skew))
        for i,(key,value) in enumerate(changes):
            with self.subTest(key=key):
                save(self.matrix,self.metadata,dict(self.arrays,**{key:value}),self.meta)
                with self.assertRaises(ValueError):workflow.spin_hall_command(self.args(f'matrix-bad-{i}'))
                self.assertFalse((self.path/f'matrix-bad-{i}.partial').exists())

    def test_memory_mu_and_source_group_cutoff_preflight(self):
        for options in [dict(memory_limit_mib=.001),dict(mu_eV=100.),dict(source_band_limit=3),
                        dict(occupied=4),dict(k_chunk=0),dict(formats=['csv','csv']),
                        dict(gap_threshold_eV=-1),dict(gap_threshold_eV=float('nan'))]:
            with self.subTest(options=options),self.assertRaises(ValueError):workflow.spin_hall_command(self.args(**options))
            self.assertFalse((self.path/'result').exists());self.assertFalse((self.path/'result.partial').exists())

    def test_midrun_failure_and_timeout_leave_failed_partial_only(self):
        original=workflow.occupied_spin_curvature;calls=[]
        def fail_second(*args,**kwargs):
            calls.append(1)
            if len(calls)==2:raise RuntimeError('independent injected second-chunk failure')
            return original(*args,**kwargs)
        with patch.object(workflow,'occupied_spin_curvature',side_effect=fail_second):
            with self.assertRaisesRegex(RuntimeError,'second-chunk'):workflow.spin_hall_command(self.args('broken',k_chunk=1))
        record=json.loads((self.path/'broken.partial/run.json').read_text())
        self.assertEqual(record['status'],'FAILED');self.assertEqual(record['completed_points'],1)
        self.assertFalse((self.path/'broken').exists());self.assertFalse((self.path/'broken.partial/spin_hall.json').exists())
        with self.assertRaisesRegex(ValueError,'time limit'):workflow.spin_hall_command(self.args('timeout',time_limit=1e-12))
        self.assertEqual(json.loads((self.path/'timeout.partial/run.json').read_text())['status'],'FAILED')
        self.assertFalse((self.path/'timeout').exists())

    def test_input_mutation_is_not_published_as_success(self):
        original=workflow.occupied_spin_curvature
        def mutate(*args,**kwargs):
            result=original(*args,**kwargs)
            self.metadata.write_text(self.metadata.read_text()+' ')
            return result
        with patch.object(workflow,'occupied_spin_curvature',side_effect=mutate):
            with self.assertRaisesRegex(ValueError,'changed during'):workflow.spin_hall_command(self.args('mutated'))
        self.assertFalse((self.path/'mutated').exists())
        self.assertEqual(json.loads((self.path/'mutated.partial/run.json').read_text())['status'],'FAILED')

    def test_hash_and_declared_dimensions_checked_before_array_use(self):
        self.matrix.write_bytes(self.matrix.read_bytes()+b'changed')
        with self.assertRaisesRegex(ValueError,'checksum'):workflow.read_matrices(self.matrix,self.metadata)
        save(self.matrix,self.metadata,self.arrays,dict(self.meta,source_nbands=100000))
        with self.assertRaisesRegex(ValueError,'memory'):workflow.read_matrices(self.matrix,self.metadata)
        save(self.matrix,self.metadata,dict(self.arrays,unexpected=np.zeros(1)),self.meta)
        with self.assertRaisesRegex(ValueError,'archive'):workflow.read_matrices(self.matrix,self.metadata)


if __name__=='__main__':unittest.main()
