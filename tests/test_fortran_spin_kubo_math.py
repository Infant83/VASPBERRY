"""Independent analytic Hamiltonian vs native sector-tangent kernel.

These are four-state model tests, not VASP/material calculations. Exact model
Hamiltonian derivatives validate the geometry separately from the production
canonical-momentum derivative approximation.
"""
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.diag([1., -1.])
IDENTITY = np.eye(2)
GENERATOR = np.kron(SX, SY)
SPIN = np.kron(SZ, IDENTITY)


def model(k, mixing=.22):
    """TR-related QWZ blocks with a smooth, periodic orbital/spin rotation."""
    x, y, z = k
    h = np.sin(x)*SX + np.sin(y)*SY + (-1+np.cos(x)+np.cos(y))*SZ
    down = -np.sin(x)*SX + np.sin(y)*SY + (-1+np.cos(x)+np.cos(y))*SZ
    zero = np.zeros((2, 2), complex)
    h0 = np.block([[h, zero], [zero, down]])
    derivatives = [
        np.block([[np.cos(x)*SX-np.sin(x)*SZ, zero],
                  [zero, -np.cos(x)*SX-np.sin(x)*SZ]]),
        np.block([[np.cos(y)*SY-np.sin(y)*SZ, zero],
                  [zero, np.cos(y)*SY-np.sin(y)*SZ]]),
        np.zeros((4, 4), complex),
    ]
    f = mixing*(np.sin(x)+np.sin(y)+.3*np.sin(z))
    df = mixing*np.array([np.cos(x), np.cos(y), .3*np.cos(z)])
    rotation = np.cos(f)*np.eye(4) + 1j*np.sin(f)*GENERATOR
    hamiltonian = rotation @ h0 @ rotation.conj().T
    derivative = np.array([rotation @ (dh+1j*df[a]*(GENERATOR@h0-h0@GENERATOR))
                           @ rotation.conj().T for a, dh in enumerate(derivatives)])
    energy, states = np.linalg.eigh(hamiltonian)
    c, outside = states[:, :2], states[:, 2:]
    dc = np.array([outside @ ((outside.conj().T @ dh @ c)
                             / (energy[:2][None, :]-energy[2:][:, None]))
                   for dh in derivative])
    return c, dc


def sectors(k, mixing=.22):
    c, _ = model(k, mixing)
    val, vec = np.linalg.eigh(c.conj().T @ SPIN @ c)
    if abs(val).min() < .1:
        raise AssertionError('model projected spin gap unexpectedly small')
    return [c@vec[:, val > 0], c@vec[:, val < 0], c]


def geometric_curvature(k, step=2e-4, mixing=.22):
    """Centered tiny-loop determinant phases, without derivative formulae."""
    result = np.empty((3, 3))
    for component, (a, b) in enumerate([(1, 2), (2, 0), (0, 1)]):
        corners = []
        for da, db in [(-.5, -.5), (.5, -.5), (.5, .5), (-.5, .5)]:
            q = np.array(k, float)
            q[a] += da*step
            q[b] += db*step
            corners.append(sectors(q, mixing))
        for sector in range(3):
            product = 1.+0j
            for j in range(4):
                product *= np.linalg.slogdet(corners[j][sector].conj().T
                                            @ corners[(j+1)%4][sector])[0]
            result[sector, component] = -np.angle(product)/step**2
    return result


def extract_subroutine(path, name):
    source = path.read_text()
    match = re.search(r'(?im)^\s*subroutine '+name+r'\b', source)
    if not match:
        raise AssertionError(f'missing {name} in {path}')
    end = re.search(r'(?im)^\s*end subroutine '+name+r'\b', source[match.start():])
    if not end:
        raise AssertionError(f'missing end of {name}')
    return source[match.start():match.start()+end.end()]+'\n'


@unittest.skipUnless(shutil.which('gfortran'), 'gfortran required')
class SpinKuboAnalyticMathTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix='spin-kubo-analytic-')
        cls.folder = Path(cls.temp.name)
        main = '''program probe
implicit none
integer ncase,i,j,l,t
complex*16 c(4,2),dc(4,2,3)
real*8 axis(3),omega(3,3),base(2,3),mix(2,3),eval(2),gmin,gcond
read(*,*) ncase
do t=1,ncase
 read(*,*) axis
 read(*,*) c
 read(*,*) dc
 call spinkubo_curvature(c,dc,2,2,axis,1d-6,omega,base,mix,eval,gmin,gcond,t)
 write(*,'(25ES25.16)') omega,base,mix,eval,gmin,gcond
enddo
end program
'''
        (cls.folder/'probe.f90').write_text(main)
        kernel = extract_subroutine(ROOT/'vaspberry_spin_kubo.inc', 'spinkubo_curvature')
        kernel += extract_subroutine(ROOT/'vaspberry_spin_kubo.inc', 'spinkubo_sigma')
        kernel += extract_subroutine(ROOT/'vaspberry_spin_kubo.inc', 'spinkubo_error')
        kernel += extract_subroutine(ROOT/'vaspberry_spin_chern.inc', 'spinchern_split')
        kernel += '''      subroutine spinchern_error(message,k)
      implicit none
      character*(*) message
      integer k
      write(0,*)message,k
      stop 1
      end subroutine spinchern_error
      subroutine vaspberry_fail
      stop 1
      end subroutine vaspberry_fail
'''
        (cls.folder/'kernel.f').write_text(kernel)
        cls.binary = cls.folder/'probe'
        result = subprocess.run(['gfortran','-O1','-fcheck=all','-ffixed-line-length-none',
                                 str(cls.folder/'probe.f90'),str(cls.folder/'kernel.f'),
                                 '-llapack','-lblas','-o',str(cls.binary)],capture_output=True,text=True)
        if result.returncode:
            raise AssertionError(result.stdout+result.stderr)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def run_kernel(self, cases):
        lines = [str(len(cases))]
        for c, dc, axis in cases:
            lines.append(' '.join(f'{v:.17e}' for v in axis))
            for array in [c, np.moveaxis(dc, 0, -1)]:
                lines.append(' '.join(f'({v.real:.17e},{v.imag:.17e})' for v in array.flatten(order='F')))
        input_path=self.folder/'input.txt'
        output_path=self.folder/'output.txt'
        error_path=self.folder/'error.txt'
        input_path.write_text('\n'.join(lines)+'\n')
        with input_path.open() as inp, output_path.open('w') as out, error_path.open('w') as err:
            result = subprocess.run([str(self.binary)],stdin=inp,stdout=out,stderr=err,
                                    env=dict(os.environ,OPENBLAS_NUM_THREADS='1'),timeout=30)
        stdout,stderr=output_path.read_text(),error_path.read_text()
        self.assertEqual(result.returncode, 0, stdout+stderr)
        rows = np.array([[float(x) for x in line.split()] for line in stdout.splitlines()])
        self.assertEqual(rows.shape, (len(cases), 25))
        return [{'omega':row[:9].reshape((3,3),order='F'),
                 'base':row[9:15].reshape((2,3),order='F'),
                 'mix':row[15:21].reshape((2,3),order='F'),
                 'eval':row[21:23], 'gmin':row[23], 'gcond':row[24]} for row in rows]

    def test_mixed_spin_curvature_matches_independent_loops_all_components(self):
        points = [(.4,.7,.2), (1.1,-.8,.3), (2.0,1.3,-.2), (-2.1,-.4,.5)]
        results = self.run_kernel([(*model(k),[0,0,1]) for k in points])
        correction = 0.
        for k, result in zip(points, results):
            np.testing.assert_allclose(result['omega'],geometric_curvature(k),atol=3e-7,rtol=2e-5)
            np.testing.assert_allclose(result['omega'][:2],result['base']+result['mix'],atol=2e-13)
            np.testing.assert_allclose(result['omega'][:2].sum(axis=0),result['omega'][2],atol=2e-13)
            correction = max(correction,abs(result['mix']).max())
        self.assertGreater(correction,1e-3, 'fixture must expose the missing-sector-derivative bug')

    def test_nonorthogonal_frame_vertical_gauge_and_axis_reversal(self):
        k = (.4,.7,.2)
        c, dc = model(k)
        matrix = np.array([[1.4,.2+.3j],[-.1j,.65]])
        rng = np.random.default_rng(937)
        vertical = rng.normal(size=(3,2,2))+1j*rng.normal(size=(3,2,2))
        altered = np.array([d@matrix+c@g for d,g in zip(dc,vertical)])
        ordinary, gauged, reversed_axis = self.run_kernel([
            (c,dc,[0,0,1]),(c@matrix,altered,[0,0,1]),(c,dc,[0,0,-1])])
        np.testing.assert_allclose(gauged['omega'],ordinary['omega'],atol=3e-13)
        np.testing.assert_allclose(reversed_axis['omega'],ordinary['omega'][[1,0,2]],atol=3e-13)
        self.assertGreater(gauged['gcond'],2.)

    def test_conserved_spin_limit_matches_analytic_two_level_curvature(self):
        points = [(.4,.7,0), (1.1,-.8,0), (2.0,1.3,0)]
        results = self.run_kernel([(*model(k,mixing=0),[0,0,1]) for k in points])
        for k,result in zip(points,results):
            x,y,_=k
            d=np.array([np.sin(x),np.sin(y),-1+np.cos(x)+np.cos(y)])
            dx=np.array([np.cos(x),0,-np.sin(x)])
            dy=np.array([0,np.cos(y),-np.sin(y)])
            exact=np.dot(d,np.cross(dx,dy))/(2*np.linalg.norm(d)**3)
            np.testing.assert_allclose(result['omega'][:,2],[exact,-exact,0],atol=3e-13)
            np.testing.assert_allclose(result['mix'],0,atol=3e-13)

    def test_exact_model_full_bz_integral_converges_to_fukui_sector_invariant(self):
        estimates = []
        for n in [8,24]:
            points = [(2*np.pi*i/n,2*np.pi*j/n,0.) for i in range(n) for j in range(n)]
            results = self.run_kernel([(*model(k),[0,0,1]) for k in points])
            estimates.append(2*np.pi*np.mean([v['omega'][:,2] for v in results],axis=0))
        n=12
        frames=[[sectors((2*np.pi*i/n,2*np.pi*j/n,0.)) for j in range(n)] for i in range(n)]
        phases=np.zeros((3,n,n))
        for i in range(n):
            for j in range(n):
                corners=[frames[i][j],frames[(i+1)%n][j],frames[(i+1)%n][(j+1)%n],frames[i][(j+1)%n]]
                for sector in range(3):
                    product=1.+0j
                    for t in range(4):
                        product*=np.linalg.slogdet(corners[t][sector].conj().T@corners[(t+1)%4][sector])[0]
                    phases[sector,i,j]=-np.angle(product)
        chern=phases.sum(axis=(1,2))/(2*np.pi)
        np.testing.assert_allclose(chern,[1,-1,0],atol=1e-12)
        np.testing.assert_allclose(estimates[-1],chern,atol=5e-6)
        self.assertLess(np.max(abs(estimates[-1]-chern)),np.max(abs(estimates[0]-chern))/50)


if __name__ == '__main__':
    unittest.main()
