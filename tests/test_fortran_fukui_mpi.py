"""Full-program MPI Fukui output regression using public synthetic PW records.

A smooth positive PW envelope times a two-component gapped lattice-model
spinor gives nonzero plaquette curvatures with no proprietary input files.
Both extrema must be found on the complete gathered native mesh.
"""
import hashlib
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def make_wavecar(path):
    n = 4
    cutoff = 31.713
    c = .262465831
    sequence = [0, 1, 2, 3, 4, -4, -3, -2, -1]
    basis = np.array([(x,y,z) for z in sequence for y in sequence for x in sequence])
    axis = np.array([0., .25, -.5, -.25])
    kpoints = np.array([(x,y,0.) for x in axis for y in axis])
    states = []
    for k in kpoints:
        g = basis[np.sum((basis+k)**2, axis=1)/c < cutoff]
        qx, qy = 2*np.pi*k[:2]
        d = [np.sin(qx)+.13*np.cos(qy), np.sin(qy)+.07*np.sin(qx),
             -1+np.cos(qx)+np.cos(qy)+.09*np.sin(qx+2*qy)]
        h = np.array([[d[2], d[0]-1j*d[1]], [d[0]+1j*d[1], -d[2]]])
        energies, u = np.linalg.eigh(h)
        envelope = np.exp(-.3*np.sum((g+k)**2, axis=1))
        envelope /= np.linalg.norm(envelope)
        states.append((k, g, energies, u, envelope))
    stride = max(256, max(2*len(g)*8 for _,g,_,_,_ in states))
    stride = ((stride+255)//256)*256
    records = bytearray((2+3*len(states))*stride)
    def put(record, values, dtype):
        payload = np.asarray(values, dtype=dtype).tobytes()
        assert len(payload) <= stride
        records[record*stride:record*stride+len(payload)] = payload
    put(0, [stride,1,45200], '<f8')
    put(1, [len(states),2,cutoff,*(2*np.pi*np.eye(3)).ravel()], '<f8')
    for ik,(k,g,energies,u,envelope) in enumerate(states):
        put(2+3*ik, [2*len(g),*k,energies[0],0,1,energies[1],0,0], '<f8')
        for band in range(2):
            put(3+3*ik+band, (u[:,band,None]*envelope).ravel(), '<c8')
    path.write_bytes(records)


def read_result(path):
    text = path.read_text()
    rows = np.loadtxt(path)
    headers = {}
    for kind in ('MAXVAL', 'MINVAL'):
        line = next(line for line in text.splitlines() if kind in line)
        headers[kind] = np.fromstring(line.split('=')[-1], sep=' ')
    chern = float(re.search(r'^# Chern Number =\s*(.*)$',text,re.M).group(1))
    return rows, headers, chern


@unittest.skipUnless(all(shutil.which(tool) for tool in ('gfortran','mpifort','mpiexec')),
                     'GNU Fortran and MPI required for the full-program MPI regression')
class NativeFukuiMpiExtremaTests(unittest.TestCase):
    def test_two_rank_curvatures_chern_and_header_extrema_match_serial(self):
        with tempfile.TemporaryDirectory(prefix='vaspberry-fukui-mpi-') as temporary:
            directory = Path(temporary)
            source = ROOT/'vaspberry.f'
            env = {**os.environ, 'OMP_NUM_THREADS':'1', 'OPENBLAS_NUM_THREADS':'1',
                   'OMPI_ALLOW_RUN_AS_ROOT':'1', 'OMPI_ALLOW_RUN_AS_ROOT_CONFIRM':'1'}
            wavecar = directory/'WAVECAR'
            make_wavecar(wavecar)
            digest = hashlib.sha256(wavecar.read_bytes()).hexdigest()
            outputs = []
            for mode,compiler,defines,prefix in [
                ('serial','gfortran',[],[]),
                ('mpi','mpifort',['-DMPI_USE'],['mpiexec','-n','2']),
            ]:
                case = directory/mode
                case.mkdir()
                binary = case/'vaspberry'
                build = subprocess.run([compiler,'-cpp','-O2',
                    '-ffixed-line-length-none','-fallow-argument-mismatch',*defines,
                    str(source),'-o',str(binary),'-llapack','-lblas'],
                    capture_output=True,text=True,env=env,timeout=120)
                self.assertEqual(build.returncode,0,build.stderr)
                run = subprocess.run([*prefix,str(binary),'--task','chern',
                    '--wavecar',str(wavecar),'--mesh','4,4','--bands','1','--spinor','2','-kp','1'],
                    cwd=case,capture_output=True,text=True,env=env,timeout=60)
                self.assertEqual(run.returncode,0,run.stdout[-1000:]+run.stderr)
                rows, headers, chern = read_result(case/'BERRYCURV.dat')
                self.assertEqual(rows.shape,(16,7))
                self.assertTrue(np.isfinite(rows).all())
                self.assertGreater(float(np.max(abs(rows[:,3]))), .1)
                for kind,index in [('MAXVAL',int(rows[:,3].argmax())),
                                   ('MINVAL',int(rows[:,3].argmin()))]:
                    self.assertAlmostEqual(headers[kind][3],rows[index,3],places=4)
                    np.testing.assert_allclose(headers[kind][:3],rows[index,4:7],atol=5.1e-5,rtol=0)
                    # This fixture deliberately avoids tied extrema, so the
                    # serial scan and MPI maxloc/minloc choose the same row.
                    self.assertEqual(np.count_nonzero(rows[:,3]==rows[index,3]),1)
                    # Native x indices 2 and 3 belong to MPI rank 1;
                    # their plaquette centers are negative. Both extrema
                    # therefore require the gathered non-root data.
                    self.assertLess(headers[kind][0], 0.)
                outputs.append((rows,headers,chern))
            np.testing.assert_array_equal(outputs[0][0],outputs[1][0])
            self.assertEqual(outputs[0][2],outputs[1][2])
            for kind in ('MAXVAL','MINVAL'):
                np.testing.assert_array_equal(outputs[0][1][kind],outputs[1][1][kind])
            self.assertEqual(hashlib.sha256(wavecar.read_bytes()).hexdigest(),digest)


if __name__ == '__main__':
    unittest.main()
