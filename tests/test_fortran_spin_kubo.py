"""End-to-end native spin-sector Kubo proxy tests; synthetic WAVECARs only."""
import csv
import io
import os
from pathlib import Path
import shutil
import struct
import subprocess
import tempfile
import unittest

import numpy as np
from test_fortran_spin_chern import fixture, table
from test_fortran_spinor_auto import write_wavecar

ROOT=Path(__file__).resolve().parents[1]


def response_fixture(case):
    """Random finite-band canonical response with nonconserved projected spin.

    The energy-degenerate two-state parent is isolated; the four retained
    plane-wave eigenvectors are orthonormal. No material/model Chern claim.
    """
    fixture(case,n=4)
    path=case/'WAVECAR';data=bytearray(path.read_bytes())
    recl=int(struct.unpack_from('d',data,0)[0]);rng=np.random.default_rng(716)
    angle=.6;rotation=np.array([[np.cos(angle/2),-np.sin(angle/2)],
                              [np.sin(angle/2),np.cos(angle/2)]])
    for k in range(16):
        rec=2+k*5;ng=int(struct.unpack_from('d',data,rec*recl)[0])//2
        c=np.zeros((2*ng,4),complex)
        a=rng.normal(size=ng)+1j*rng.normal(size=ng);a/=np.linalg.norm(a)
        b=rng.normal(size=ng)+1j*rng.normal(size=ng);b/=np.linalg.norm(b)
        c[:ng,0]=a;c[ng:,1]=b
        for n in (2,3):
            v=rng.normal(size=2*ng)+1j*rng.normal(size=2*ng)
            v-=c[:,:n]@(c[:,:n].conj().T@v);c[:,n]=v/np.linalg.norm(v)
        c=np.einsum('st,tgn->sgn',rotation,c.reshape(2,ng,4)).reshape(2*ng,4)
        for n in range(4):
            arr=c[:,n].astype(np.complex64);start=(rec+n+1)*recl
            data[start:start+arr.nbytes]=arr.tobytes()
    path.write_bytes(data)


@unittest.skipUnless(shutil.which('gfortran'),'gfortran required')
class NativeSpinKuboTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp=tempfile.TemporaryDirectory(prefix='vaspberry-spin-kubo-');cls.work=Path(cls.tmp.name)
        flags=['-cpp','-O1','-g','-fcheck=all','-ffixed-line-length-none','-fallow-argument-mismatch']
        cls.serial=cls.work/'vaspberry';cls.parallel=cls.work/'vaspberry-mpi'
        cls.invoke(['gfortran',*flags,str(ROOT/'vaspberry.f'),'-llapack','-lblas','-o',str(cls.serial)],cls.work)
        cls.has_mpi=bool(shutil.which('mpifort') and shutil.which('mpiexec'))
        if cls.has_mpi:cls.invoke(['mpifort',*flags,'-DMPI_USE',str(ROOT/'vaspberry.f'),'-llapack','-lblas','-o',str(cls.parallel)],cls.work)

    @classmethod
    def tearDownClass(cls):cls.tmp.cleanup()

    @staticmethod
    def invoke(cmd,cwd,success=True):
        env=dict(os.environ,OMPI_ALLOW_RUN_AS_ROOT='1',OMPI_ALLOW_RUN_AS_ROOT_CONFIRM='1')
        p=subprocess.run(cmd,cwd=cwd,capture_output=True,text=True,env=env,timeout=120)
        if success and p.returncode:raise AssertionError(f'{cmd}\n{p.stdout}\n{p.stderr}')
        return p

    def execute(self,name,flags=(),mpi=False,mode='response',success=True):
        case=self.work/(self._testMethodName+'-'+name)
        if mode=='response':response_fixture(case)
        else:fixture(case,n=4,mode=mode)
        cmd=[str(self.parallel if mpi else self.serial),'--task','spin-kubo', '--kubo-source', 'wavecar','--bands','1:2',*flags]
        if mpi:cmd=['mpiexec','-n','2',*cmd]
        p=self.invoke(cmd,case,success)
        if not success:
            self.assertNotEqual(p.returncode,0);self.assertFalse((case/'SPIN_KUBO.csv').exists())
        return case,p

    def values(self,case,key):return np.array([float(row[key]) for row in table(case/'SPIN_KUBO.csv')[1]])

    def test_point_mode_three_components_and_sector_derivative(self):
        case,_=self.execute('point')
        self.assertFalse((case/'SPIN_KUBO_INTEGRAL.csv').exists())
        meta,rows=table(case/'SPIN_KUBO.csv')
        self.assertEqual(meta['integration'],'NONE_ARBITRARY_K_LIST')
        self.assertEqual(meta['derivative'],'CANONICAL_MOMENTUM_DERIVATIVE_PROXY')
        self.assertEqual(meta['physical_paw_validated'],'false')
        self.assertEqual(meta['result_status'],'PASS')
        self.assertEqual(len(rows),16);self.assertEqual(len(table(case/'SPIN_KUBO_SPECTRUM.csv')[1]),32)
        nonzero=[]
        for component in ('yz','zx','xy'):
            plus=self.values(case,f'omega_plus_{component}_A2');minus=self.values(case,f'omega_minus_{component}_A2')
            parent=self.values(case,f'omega_parent_{component}_A2')
            np.testing.assert_allclose(plus+minus,parent,rtol=1e-11,atol=1e-11)
            mix=self.values(case,f'spin_mixing_plus_{component}_A2')
            np.testing.assert_allclose(mix,-self.values(case,f'spin_mixing_minus_{component}_A2'),atol=1e-11)
            np.testing.assert_allclose(plus,self.values(case,f'parent_projected_plus_{component}_A2')+mix,atol=1e-11)
            nonzero.append(np.max(abs(mix)))
        self.assertGreater(max(nonzero),1e-6)

    def test_axis_reversal_swaps_sector_curvature(self):
        normal,_=self.execute('normal');reverse,_=self.execute('reverse',flags=['--spin-axis','0,0,-1'])
        for component in ('yz','zx','xy'):
            np.testing.assert_allclose(self.values(normal,f'omega_plus_{component}_A2'),self.values(reverse,f'omega_minus_{component}_A2'),rtol=1e-10,atol=1e-10)

    def test_sum_bands_default_and_invalid_boundaries(self):
        normal,_=self.execute('default')
        explicit,_=self.execute('all',flags=['--sum-bands','4'])
        self.assertEqual((normal/'SPIN_KUBO.csv').read_bytes(),
                         (explicit/'SPIN_KUBO.csv').read_bytes())
        meta,_=table(explicit/'SPIN_KUBO.csv')
        self.assertEqual(meta['source_nbands'],'4')
        self.assertEqual(meta['sum_band_max'],'4')
        self.assertEqual(meta['summed_external_bands'],'2')
        for cutoff in ('0','1','2','5','-1','3.5'):
            self.execute('invalid-'+cutoff,flags=['--sum-bands',cutoff],success=False)
        _,p=self.execute('split',flags=['--sum-bands','3'],success=False)
        self.assertIn('cuts an unresolved external degeneracy',p.stderr)
        case=self.work/(self._testMethodName+'-wrong-task');response_fixture(case)
        p=self.invoke([str(self.serial),'--task','spin-chern','--bands','1:2',
                       '--mesh','4,4','--sum-bands','4'],case,False)
        self.assertNotEqual(p.returncode,0)
        self.assertIn('only for --task spin-kubo',p.stderr)
        self.assertFalse((case/'SPIN_CHERN.csv').exists())

    def test_partial_sum_matches_zero_omitted_couplings(self):
        # One nondegenerate omitted external state: removing its coupling must
        # equal the truncated sum without changing the selected spin projector.
        cases={}
        for name in ('full','partial','zero-coupling'):
            case=self.work/(self._testMethodName+'-'+name);response_fixture(case)
            path=case/'WAVECAR';data=bytearray(path.read_bytes())
            recl=int(struct.unpack_from('d',data,0)[0])
            for k in range(16):
                rec=2+k*5
                struct.pack_into('d',data,rec*recl+8*(4+3*3),2.)
                if name=='zero-coupling':
                    nc=int(struct.unpack_from('d',data,rec*recl)[0])
                    start=(rec+4)*recl;data[start:start+8*nc]=bytes(8*nc)
            path.write_bytes(data)
            args=['--sum-bands','3'] if name=='partial' else []
            self.invoke([str(self.serial),'--task','spin-kubo', '--kubo-source', 'wavecar','--bands','1:2',*args],case)
            cases[name]=case
        change=0.
        for sector in ('plus','minus','parent'):
            for component in ('yz','zx','xy'):
                key=f'omega_{sector}_{component}_A2'
                a=self.values(cases['partial'],key)
                np.testing.assert_allclose(a,self.values(cases['zero-coupling'],key),
                                           rtol=2e-12,atol=2e-12)
                change=max(change,float(abs(a-self.values(cases['full'],key)).max()))
        self.assertGreater(change,1e-6)
        meta,_=table(cases['partial']/'SPIN_KUBO.csv')
        self.assertEqual(meta['source_nbands'],'4')
        self.assertEqual(meta['sum_band_max'],'3')
        self.assertEqual(meta['summed_external_bands'],'1')
        self.assertAlmostEqual(float(meta['sum_boundary_gap_eV']),1.)

    def test_explicit_mesh_enables_unrounded_vector_integral(self):
        case,_=self.execute('mesh',flags=['--mesh','4,4'])
        meta,rows=table(case/'SPIN_KUBO_INTEGRAL.csv');row=rows[0]
        self.assertEqual(meta['geometric_integer_certified'],'false')
        # The fixture's reciprocal cell is identity, so area vector=(0,0,1/16).
        for sector in ('plus','minus','parent'):
            expected=self.values(case,f'omega_{sector}_xy_A2').sum()/(16*2*np.pi)
            self.assertAlmostEqual(float(row[f'C_est_{sector}']),expected,places=12)
        _,p=self.execute('bad-mesh',flags=['--mesh','3,3'],success=False)
        self.assertIn('complete NX by NY',p.stderr)

    def test_tilted_cell_uses_all_curvature_components_in_integral(self):
        normal,_=self.execute('normal',flags=['--mesh','4,4'])
        tilted=self.work/(self._testMethodName+'-tilted');response_fixture(tilted)
        angle=.43;rot=np.array([[1,0,0],[0,np.cos(angle),-np.sin(angle)],
                              [0,np.sin(angle),np.cos(angle)]])
        lattice=2*np.pi*rot.T
        path=tilted/'WAVECAR';data=bytearray(path.read_bytes());recl=int(struct.unpack_from('d',data,0)[0])
        struct.pack_into('9d',data,recl+24,*lattice.ravel());path.write_bytes(data)
        outcar=(tilted/'OUTCAR').read_text().splitlines()
        for i,row in enumerate(lattice):outcar[2+i]=' '.join(f'{x:.12f}' for x in row)+' 0 0 0'
        (tilted/'OUTCAR').write_text('\n'.join(outcar)+'\n')
        self.invoke([str(self.serial),'--task','spin-kubo', '--kubo-source', 'wavecar','--bands','1:2','--mesh','4,4'],tilted)
        for sector in ('plus','minus','parent'):
            expected=np.array([self.values(normal,f'omega_{sector}_{x}_A2') for x in ('yz','zx','xy')]).T@rot.T
            observed=np.array([self.values(tilted,f'omega_{sector}_{x}_A2') for x in ('yz','zx','xy')]).T
            np.testing.assert_allclose(observed,expected,rtol=2e-10,atol=2e-10)
            x=float(table(normal/'SPIN_KUBO_INTEGRAL.csv')[1][0][f'C_est_{sector}'])
            y=float(table(tilted/'SPIN_KUBO_INTEGRAL.csv')[1][0][f'C_est_{sector}'])
            self.assertAlmostEqual(x,y,places=11)

    def test_duplicate_k_points_allowed_for_path_rejected_for_mesh(self):
        case,_=self.execute('path',mode='duplicate')
        self.assertEqual(len(table(case/'SPIN_KUBO.csv')[1]),16)
        _,p=self.execute('mesh',mode='duplicate',flags=['--mesh','4,4'],success=False)
        self.assertIn('duplicate mesh',p.stderr)

    def test_gap_and_scalar_errors_before_output(self):
        for mode,message in [('closed-energy','energy isolated'),('closed-spin','zero gap'),('singular-gram','Gram metric')]:
            _,p=self.execute(mode,mode=mode,success=False);self.assertIn(message,p.stderr)
        case=self.work/(self._testMethodName+'-scalar');case.mkdir();write_wavecar(case/'WAVECAR',components=1)
        p=self.invoke([str(self.serial),'--task','spin-kubo', '--kubo-source', 'wavecar','--bands','1'],case,False)
        self.assertNotEqual(p.returncode,0);self.assertIn('two-component spinors',p.stderr)

    def test_exact_degenerate_gauge_rotation(self):
        case,_=self.execute('base');rot=self.work/(self._testMethodName+'-rot');shutil.copytree(case,rot)
        for p in rot.glob('SPIN*'):p.unlink()
        path=rot/'WAVECAR';data=bytearray(path.read_bytes());recl=int(struct.unpack_from('d',data,0)[0]);rng=np.random.default_rng(981)
        for k in range(16):
            rec=2+5*k;nc=int(struct.unpack_from('d',data,rec*recl)[0])
            c=np.array([np.frombuffer(data,dtype=np.complex64,count=nc,offset=(rec+n+1)*recl).copy() for n in range(4)]).T
            for sl in (slice(0,2),slice(2,4)):
                u,_=np.linalg.qr(rng.normal(size=(2,2))+1j*rng.normal(size=(2,2)));c[:,sl]=c[:,sl]@u
            for n in range(4):
                a=c[:,n].astype(np.complex64);start=(rec+n+1)*recl;data[start:start+a.nbytes]=a.tobytes()
        path.write_bytes(data)
        self.invoke([str(self.serial),'--task','spin-kubo', '--kubo-source', 'wavecar','--bands','1:2'],rot)
        for sector in ('plus','minus','parent'):
            for component in ('yz','zx','xy'):
                np.testing.assert_allclose(self.values(case,f'omega_{sector}_{component}_A2'),self.values(rot,f'omega_{sector}_{component}_A2'),rtol=2e-5,atol=1e-7)

    def test_mpi_equality_and_collective_failure(self):
        if not self.has_mpi:self.skipTest('MPI required')
        serial,_=self.execute('serial');mpi,_=self.execute('mpi',mpi=True)
        for name in ('SPIN_KUBO.csv','SPIN_KUBO_SPECTRUM.csv'):
            self.assertEqual((serial/name).read_bytes(),(mpi/name).read_bytes())
        _,p=self.execute('bad-mpi',mpi=True,mode='closed-spin',success=False)
        self.assertIn('zero gap',p.stderr)


if __name__=='__main__':unittest.main()
