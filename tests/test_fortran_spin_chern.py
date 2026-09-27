"""Native executable tests with explicitly synthetic plane-wave Chern fixtures.

The fixture embeds two opposite QWZ Chern blocks into smooth plane-wave
orbital envelopes. This is a model WAVECAR format test, not a VASP result.
"""
import csv
import io
import itertools
import json
import os
from pathlib import Path
import shutil
import struct
import subprocess
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def fixture(folder, n=8, mode='normal', energy_scale=1., frame=None):
    folder.mkdir(parents=True, exist_ok=True)
    lat = 2*np.pi*np.eye(3)
    cutoff, nb, recl = 100., 4, 16384
    points = np.array([(i/n, j/n, 0.) for i in range(n) for j in range(n)])
    if mode == 'duplicate':
        points[-1] = points[0]
    buf = bytearray(recl*(2+len(points)*(nb+1)))
    struct.pack_into('3d', buf, 0, recl, 1, 45200)
    struct.pack_into('12d', buf, recl, len(points), nb, cutoff, *lat.ravel())
    seq = list(range(9))+list(range(-8, 0))
    xyz = np.array([(x,y,z) for z in seq for y in seq for x in seq])
    rng = np.random.default_rng(412)
    for k, point in enumerate(points):
        g = xyz[np.sum((xyz+point)**2, axis=1)/.262465831 < cutoff]
        ng = len(g)
        orbital = np.zeros((2, ng), complex)
        for o in range(2):
            orbital[o] = (g[:,2] == o)*np.exp(-np.sum((g[:,:2]+point[:2])**2, axis=1)/.7)
            orbital[o] /= np.linalg.norm(orbital[o])
        x,y = 2*np.pi*point[:2]
        h = np.array([[ -1+np.cos(x)+np.cos(y), np.sin(x)-1j*np.sin(y)],
                      [np.sin(x)+1j*np.sin(y), 1-np.cos(x)-np.cos(y)]])
        _, u = np.linalg.eigh(h)
        c = np.zeros((4,2,ng), complex)
        c[0,0] = u[:,0]@orbital
        c[1,1] = u[:,0].conj()@orbital
        c[2,0] = u[:,1]@orbital
        c[3,1] = u[:,1].conj()@orbital
        if mode == 'spin-mixed':
            angle = .5
            spinrotation = np.array([[np.cos(angle/2),-np.sin(angle/2)],
                                     [np.sin(angle/2),np.cos(angle/2)]])
            c = np.einsum('st,ntg->nsg',spinrotation,c)
        if mode in ('rotated','metric'):
            w,_ = np.linalg.qr(rng.normal(size=(2,2))+1j*rng.normal(size=(2,2)))
            if mode == 'metric':
                w = w@np.array([[1.4,.2j],[0.,.6]])
            c[:2] = np.einsum('mn,msg->nsg',w,c[:2])
        if mode == 'closed-spin':
            c[:2] = orbital[:,None,:]/np.sqrt(2)*np.ones((1,2,1))
        if mode == 'rank-change' and k == len(points)-1:
            c[1,0] = u[:,1]@orbital
            c[1,1] = 0
        if mode == 'singular-gram':
            c[1] = c[0]
        if mode == 'nonfinite':
            c[0,0,0] = np.nan
        e = energy_scale*np.array([-1.,-1.,1.,1.])
        if mode == 'closed-energy':
            e[2] = e[1]
        rec = 2+k*(nb+1)
        triples = []
        for band in range(nb):
            triples.extend((e[band], 0., (.9 if mode == 'smeared' else 1.) if band < 2 else (.1 if mode == 'smeared' else 0.)))
        struct.pack_into('16d',buf,rec*recl,2*ng,*point,*triples)
        for band in range(nb):
            arr = c[band].astype(np.complex64).reshape(-1)
            start=(rec+band+1)*recl
            buf[start:start+arr.nbytes]=arr.tobytes()
    (folder/'WAVECAR').write_bytes(buf)
    rot = np.eye(3) if frame is None else np.asarray(frame)
    text = f' NKPTS = {len(points)} k-points in BZ NBANDS= {nb}\n'
    text += ' direct lattice vectors                 reciprocal lattice vectors\n'
    text += ''.join(' '.join(f'{v:.12f}' for v in row)+'  0 0 0\n' for row in lat)
    text += ' transformation matrix from SAXIS to cartesian coordinates\n -----------------\n'
    text += ''.join(f' {row[0]:.12f} m_x {row[1]:.12f} m_y {row[2]:.12f} m_z\n' for row in rot)
    (folder/'OUTCAR').write_text(text)
    return n


def table(path):
    text = path.read_text()
    meta = dict(line[2:].split('=',1) for line in text.splitlines() if line.startswith('# ') and '=' in line)
    rows = list(csv.DictReader(io.StringIO('\n'.join(line for line in text.splitlines() if not line.startswith('#')))))
    return meta, rows


@unittest.skipUnless(shutil.which('gfortran'), 'gfortran required')
class NativeSpinChernTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix='vaspberry-spin-chern-')
        cls.work = Path(cls.temp.name)
        cls.serial = cls.work/'vaspberry'
        cls.parallel = cls.work/'vaspberry-mpi'
        flags=['-cpp','-O1','-g','-fcheck=all','-ffixed-line-length-none','-fallow-argument-mismatch']
        cls.invoke(['gfortran',*flags,str(ROOT/'vaspberry.f'),'-llapack','-lblas','-o',str(cls.serial)],cls.work)
        cls.has_mpi=bool(shutil.which('mpifort') and shutil.which('mpiexec'))
        if cls.has_mpi:
            cls.invoke(['mpifort',*flags,'-DMPI_USE',str(ROOT/'vaspberry.f'),'-llapack','-lblas','-o',str(cls.parallel)],cls.work)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    @staticmethod
    def invoke(command,cwd,success=True):
        env=dict(os.environ,OMPI_ALLOW_RUN_AS_ROOT='1',OMPI_ALLOW_RUN_AS_ROOT_CONFIRM='1')
        result=subprocess.run(command,cwd=cwd,capture_output=True,text=True,env=env,timeout=120)
        if success and result.returncode:
            raise AssertionError(f'{command}\n{result.stdout}\n{result.stderr}')
        return result

    def execute(self,name,mode='normal',flags=(),mpi=False,success=True,**kw):
        case=self.work/(self._testMethodName+'-'+name)
        n=fixture(case,mode=mode,**kw)
        cmd=[str(self.parallel if mpi else self.serial),'--task','spin-chern','--bands','1:2','--mesh',f'{n},{n}',*flags]
        if mpi:cmd=['mpiexec','-n','2',*cmd]
        result=self.invoke(cmd,case,success)
        if not success:
            self.assertNotEqual(result.returncode,0)
            self.assertFalse((case/'SPIN_CHERN.csv').exists())
        return case,result

    def test_opposite_chern_sectors_outputs_and_default_input(self):
        case,result=self.execute('normal')
        meta,rows=table(case/'SPIN_CHERN.csv'); row=rows[0]
        self.assertEqual(meta['result_status'],'PASS')
        self.assertEqual(meta['representation'],'WAVECAR_PSEUDO_GRAM')
        self.assertEqual(meta['physical_paw_validated'],'false')
        self.assertAlmostEqual(abs(float(row['C_plus'])),1.,places=10)
        self.assertAlmostEqual(float(row['C_plus']),-float(row['C_minus']),places=10)
        self.assertAlmostEqual(float(row['C_charge']),0.,places=10)
        self.assertEqual((row['rank_plus'],row['rank_minus']),('1','1'))
        self.assertEqual(len(table(case/'SPIN_BERRY.csv')[1]),64)
        self.assertEqual(len(table(case/'SPIN_SPECTRUM.csv')[1]),128)
        self.assertIn('result_status=PASS',result.stdout)

    def test_arbitrary_local_gauge_and_nonorthogonal_metric(self):
        baseline,_=self.execute('normal')
        target=table(baseline/'SPIN_CHERN.csv')[1][0]
        for mode in ('rotated','metric'):
            case,_=self.execute(mode,mode=mode)
            row=table(case/'SPIN_CHERN.csv')[1][0]
            for name in ('C_plus','C_minus','C_parent'):
                self.assertAlmostEqual(float(row[name]),float(target[name]),places=10)

    def test_axis_reversal_and_rotated_source_frame(self):
        case,_=self.execute('reversed',flags=['--spin-axis','0,0,-1'])
        baseline,_=self.execute('normal')
        a=table(case/'SPIN_CHERN.csv')[1][0];b=table(baseline/'SPIN_CHERN.csv')[1][0]
        self.assertAlmostEqual(float(a['C_plus']),float(b['C_minus']),places=10)
        # Local z is Cartesian x under this proper rotation.
        rotated,_=self.execute('frame',frame=[[0,0,1],[0,1,0],[-1,0,0]],flags=['--spin-axis','x'])
        c=table(rotated/'SPIN_CHERN.csv')[1][0]
        self.assertAlmostEqual(float(c['C_plus']),float(b['C_plus']),places=10)

    def test_spin_mixed_gapped_sectors_and_transverse_closure(self):
        baseline,_=self.execute('normal')
        mixed,_=self.execute('mixed',mode='spin-mixed')
        a=table(baseline/'SPIN_CHERN.csv')[1][0]
        meta,b=table(mixed/'SPIN_CHERN.csv')
        self.assertAlmostEqual(float(a['C_plus']),float(b[0]['C_plus']),places=10)
        self.assertAlmostEqual(float(a['C_minus']),float(b[0]['C_minus']),places=10)
        ev=np.array([float(row['projected_pauli_eigenvalue']) for row in table(mixed/'SPIN_SPECTRUM.csv')[1]])
        self.assertGreater(np.min(abs(ev)),.8)
        self.assertLess(np.min(abs(ev)),.95)
        _,failed=self.execute('transverse',flags=['--spin-axis','x'],success=False)
        self.assertIn('zero gap',failed.stderr)

    def test_parser_rejects_malformed_axis_and_tolerances(self):
        for i,flags in enumerate((['--spin-axis','1,,0'],['--spin-axis','1,0,0,0'],
                                  ['--spin-axis','1,1,0'],['--spin-gap-tol','NaN'],
                                  ['--energy-gap-tol','1e-8,2'])):
            _,failed=self.execute(str(i),flags=flags,success=False)
            self.assertIn('invalid value',failed.stderr)

    def test_failures_before_outputs(self):
        for mode,message in [('closed-spin','zero gap'),('closed-energy','energy isolated'),
                             ('singular-gram','Gram metric'),('rank-change','rank changes'),
                             ('nonfinite','nonfinite coefficients'),('duplicate','duplicate mesh')]:
            with self.subTest(mode=mode):
                _,result=self.execute(mode,mode=mode,success=False)
                self.assertIn(message,result.stderr)

    def test_small_energy_gap_and_fractional_occupations(self):
        case,_=self.execute('micro-eV',energy_scale=2e-6,mode='smeared')
        meta,_=table(case/'SPIN_CHERN.csv')
        self.assertAlmostEqual(float(meta['min_external_energy_gap_eV']),4e-6)
        self.assertEqual(meta['source_occupations_match_integer_projector'],'F')
        _,result=self.execute('threshold',energy_scale=2e-6,flags=['--energy-gap-tol','1e-5'],success=False)
        self.assertIn('energy isolated',result.stderr)
        _,result=self.execute('cut-parent',flags=['--bands','1'],success=False)
        self.assertIn('energy isolated',result.stderr)

    def test_frame_provenance_and_no_clobber(self):
        for mode in ('missing','wrong-dim','unknown-frame'):
            case=self.work/(self._testMethodName+mode);fixture(case)
            if mode=='missing':(case/'OUTCAR').unlink()
            elif mode=='wrong-dim':(case/'OUTCAR').write_text((case/'OUTCAR').read_text().replace('NBANDS= 4','NBANDS= 5'))
            else:(case/'OUTCAR').write_text('SAXIS = 0 0 1\n')
            result=self.invoke([str(self.serial),'--task','spin-chern','--bands','1:2','--mesh','8,8'],case,False)
            self.assertNotEqual(result.returncode,0);self.assertFalse((case/'SPIN_CHERN.csv').exists())
        case,_=self.execute('valid');before=(case/'SPIN_CHERN.csv').read_bytes()
        result=self.invoke([str(self.serial),'--task','spin-chern','--bands','1:2','--mesh','8,8'],case,False)
        self.assertNotEqual(result.returncode,0);self.assertEqual(before,(case/'SPIN_CHERN.csv').read_bytes())

    def test_mpi_agrees_and_collective_failure(self):
        if not self.has_mpi:self.skipTest('MPI required')
        serial,_=self.execute('serial',mode='rotated')
        parallel,_=self.execute('mpi',mode='rotated',mpi=True)
        for name in ('CHERN','BERRY','SPECTRUM'):
            self.assertEqual((serial/f'SPIN_{name}.csv').read_bytes(),(parallel/f'SPIN_{name}.csv').read_bytes())
        _,failed=self.execute('bad-mpi',mode='closed-spin',mpi=True,success=False)
        self.assertIn('zero gap',failed.stderr)


if __name__=='__main__':unittest.main()
