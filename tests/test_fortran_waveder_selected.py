"""Native selected WAVEDER geometry with an independent full-matrix oracle.

Fixtures are synthetic Hermitian optical connections. The tests contract all
source bands directly, independently of the production unordered-pair kernel.
"""
import math
import os
from pathlib import Path
import re
import shutil
import struct
import subprocess
import tempfile
import unittest

import numpy as np
from test_fortran_waveder_default import ROOT, table, write_sources


def write_selected_sources(path, *, nd=3, energies=(-3., -1., 1., 3.),
                           components=2, channels=1, q=((0., 0., 0.), (.5, 0., 0.)),
                           full=None):
    """Create a same-run synthetic quartet; return C[bra,ket,k,spin,xyz]."""
    write_sources(path, nd=nd, components=components, channels=channels, q=q)
    nk, nb, recl = len(q), 4, 2048
    energies = np.broadcast_to(np.asarray(energies, float), (channels, nk, nb))
    blob = bytearray((path/'WAVECAR').read_bytes())
    for s in range(channels):
        for k in range(nk):
            rec = 2+(s*nk+k)*(nb+1)
            for band in range(nb):
                struct.pack_into('<d', blob, rec*recl+32+band*24, energies[s, k, band])
    (path/'WAVECAR').write_bytes(blob)
    lines = (path/'OUTCAR').read_text().splitlines()
    s, k = 0, -1
    for index, line in enumerate(lines):
        if line.startswith(' spin component'):
            s = int(line.split()[-1])-1;k = -1
        if 'band No.' in line:
            k += 1
            for band in range(nb):
                parts = lines[index+1+band].split()
                parts[1] = f'{energies[s, k, band]:.4f}'
                lines[index+1+band] = ' '.join(parts)
    (path/'OUTCAR').write_text('\n'.join(lines)+'\n')
    if full is None:
        full = np.zeros((nb, nb, nk, channels, 3), np.complex64)
        for i in range(nb):
            for j in range(i+1, nb):
                for k in range(nk):
                    for s in range(channels):
                        full[j, i, k, s] = [1+.125j*(i+1), .25*(j+1)+1j*(k+s+j+1), -.5j]
                        full[i, j, k, s] = full[j, i, k, s].conj()
    full = np.asarray(full, dtype=np.complex64)
    def record(payload):
        return struct.pack('<i', len(payload))+payload+struct.pack('<i', len(payload))
    blob = record(struct.pack('<4i', nb, nd, nk, channels))
    blob += record(struct.pack('<d', 0.))+record(np.zeros((3, 3), '<f8').tobytes())
    blob += record(full[:, :nd].tobytes(order='F'))
    (path/'WAVEDER').write_bytes(blob)
    return full


def oracle(full, bands):
    c = full.astype(np.complex128)
    # All M=NBANDS, separately for every target n. No energy denominator.
    per = -2*np.imag((c[..., 0].conj()*c[..., 1]).sum(axis=0))
    return per[np.asarray(bands)-1].transpose(2, 1, 0)


@unittest.skipUnless(shutil.which('gfortran'), 'gfortran required')
class NativeSelectedWavederTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix='vaspberry-selected-')
        cls.work = Path(cls.temp.name)
        cls.binary = cls.work/'vaspberry'
        cls.mpi = cls.work/'vaspberry-mpi'
        flags = ['-cpp', '-O0', '-g', '-fcheck=all', '-ffixed-line-length-none', '-fallow-argument-mismatch']
        cls.invoke(['gfortran', *flags, str(ROOT/'vaspberry.f'), '-llapack', '-lblas', '-o', str(cls.binary)], cls.work)
        cls.has_mpi = bool(shutil.which('mpifort') and shutil.which('mpiexec'))
        if cls.has_mpi:
            cls.invoke(['mpifort', *flags, '-DMPI_USE', str(ROOT/'vaspberry.f'), '-llapack', '-lblas', '-o', str(cls.mpi)], cls.work)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    @staticmethod
    def invoke(command, cwd):
        env = dict(os.environ, OMPI_ALLOW_RUN_AS_ROOT='1', OMPI_ALLOW_RUN_AS_ROOT_CONFIRM='1',
                   OMPI_MCA_rmaps_base_oversubscribe='1')
        result = subprocess.run(command, cwd=cwd, capture_output=True, text=True, env=env, timeout=90)
        if result.returncode:
            raise AssertionError(f'{command}\n{result.stdout}\n{result.stderr}')
        return result

    def case(self, name, **kw):
        path = self.work/(self._testMethodName+'-'+name)
        full = write_selected_sources(path, **kw)
        return path, full

    def run_case(self, path, bands, *args, reject=None, mpi=False):
        command = [str(self.mpi if mpi else self.binary), '--task', 'kubo', '--bands', bands, *args]
        if mpi:
            command = ['mpiexec', '-n', '2', *command]
        if reject is not None:
            result = subprocess.run(command, cwd=path, capture_output=True, text=True, timeout=90)
            self.assertNotEqual(result.returncode, 0, result.stdout)
            self.assertIn(reject, result.stderr)
            self.assertFalse((path/'KUBO_WAVEDER.csv').exists())
            self.assertFalse((path/'KUBO_WAVEDER.csv.partial').exists())
            return result
        self.invoke(command, path)
        return table(path/'KUBO_WAVEDER.csv')

    def test_noncontiguous_trace_and_per_band_use_all_virtual_states(self):
        path, full = self.case('trace')
        meta, rows = self.run_case(path, '3,1')
        self.assertEqual(meta['band_ids'], '1,3')
        self.assertEqual(meta['band_rank'], '2')
        self.assertEqual(meta['schema'], 'VASPBERRY_WAVEDER_KUBO_BUNDLE_V1')
        np.testing.assert_allclose([float(r['omega_z_A2']) for r in rows], oracle(full, [1, 3]).sum(axis=2).ravel(), atol=1e-14)
        # Band 4 is outside S but makes a large, essential contribution.
        truncated = full.copy();truncated[3] = 0;truncated[:, 3] = 0
        self.assertGreater(np.max(np.abs(oracle(full, [1, 3])-oracle(truncated, [1, 3]))), 1.)
        path, full = self.case('bands')
        meta, rows = self.run_case(path, '1,3:4', '--per-band', '1')
        self.assertEqual(meta['schema'], 'VASPBERRY_WAVEDER_KUBO_BAND_V1')
        self.assertEqual(meta['band_ids'], '1,3,4')
        self.assertEqual([r['band'] for r in rows], ['1','3','4']*2)
        np.testing.assert_allclose([float(r['omega_z_A2']) for r in rows], oracle(full, [1,3,4]).ravel(), atol=1e-14)

    def test_reverse_orientation_and_missing_high_high_coverage(self):
        path, full = self.case('reverse', nd=3)
        _, rows = self.run_case(path, '4', '--per-band', '1')
        np.testing.assert_allclose([float(r['omega_z_A2']) for r in rows], oracle(full,[4]).ravel(), atol=1e-14)
        for perband in (False, True):
            path, _ = self.case('missing'+str(perband), nd=2)
            args = ('--per-band','1') if perband else ()
            self.run_case(path,'3',*args,reject='high-high block')
        path, full = self.case('all-high',nd=2)
        _, rows = self.run_case(path,'3:4')
        np.testing.assert_allclose([float(r['omega_z_A2']) for r in rows],oracle(full,[3,4]).sum(axis=2).ravel(),atol=1e-14)

    def test_complete_cluster_trace_and_full_source_transitive_bridge(self):
        energies=(-2.,-1.9985,-1.997,3.)
        for spec, perband in [('1,3',False),('3',False),('1:3',True)]:
            path, _ = self.case(spec.replace(':','_')+str(perband), energies=energies)
            self.run_case(path,spec,*(['--per-band','1'] if perband else []),reject='producer cluster')
        path, full = self.case('complete',energies=energies)
        _, rows = self.run_case(path,'1:3')
        np.testing.assert_allclose([float(r['omega_z_A2']) for r in rows],oracle(full,[1,2,3]).sum(axis=2).ravel(),atol=1e-14)

    def test_internal_degenerate_basis_rotation_leaves_trace(self):
        path, full = self.case('original',energies=(-2.,-2.,1.,3.))
        _, rows = self.run_case(path,'1:2')
        target=np.array([float(r['omega_z_A2']) for r in rows])
        u=np.eye(4,dtype=complex);u[:2,:2]=np.array([[1,1j],[1j,1]])/np.sqrt(2)
        rotated=np.einsum('ai,abksd,bj->ijksd',u.conj(),full.astype(complex),u)
        path, _=self.case('rotated',energies=(-2.,-2.,1.,3.),full=rotated)
        _, rows=self.run_case(path,'1:2')
        np.testing.assert_allclose([float(r['omega_z_A2']) for r in rows],target,rtol=2e-7,atol=2e-7)

    def test_required_hermiticity_rejects_only_needed_pairs(self):
        path, full=self.case('internal')
        full[0,1,:,0,0]+=100
        write_selected_sources(path,full=full)
        self.run_case(path,'1:2') # bad pair is internal and cancels
        path,_=self.case('required',full=full)
        self.run_case(path,'1',reject='Hermiticity')

    def test_full_basis_zero_has_explicit_truncated_scope(self):
        path,_=self.case('full',nd=1)
        meta,rows=self.run_case(path,'1:4')
        self.assertEqual(meta['no_external_states'],'true')
        self.assertEqual(meta['zero_trace_scope'],'TRUNCATED_WAVEDER_BASIS')
        self.assertTrue(all(float(r['omega_z_A2'])==0 and r['min_external_gap_eV']=='NA' for r in rows))

    def test_selected_integral_is_unit_weight_geometry_not_total_hall(self):
        q=((0.,0.,0.),(.5,0.,0.),(0.,.5,0.),(.5,.5,0.))
        for channels in (1,2):
            path,full=self.case(str(channels),components=1,channels=channels,q=q)
            meta,_=self.run_case(path,'3:4','--mesh','2,2')
            expected=oracle(full,[3,4]).sum()/(4*2*math.pi)
            self.assertAlmostEqual(float(meta['selected_chern_sum']),expected,delta=1e-13)
            self.assertNotIn('total_chern',meta);self.assertNotIn('sheet_hall_e2_over_h',meta)
            self.assertEqual(meta['occupation_weighting'],'NONE')

    def test_duplicates_bounds_and_nonstandard_list_routing_fail(self):
        for value,text in [('1,1','duplicate'),('1:3,2','duplicate'),('0','invalid'),('1:5','outside'),('1,,3','invalid'),('3:1','invalid')]:
            path,_=self.case(value.replace(':','_').replace(',','x'))
            self.run_case(path,value,reject=text)
        path,_=self.case('wavecar')
        self.run_case(path,'1,3','--kubo-source','wavecar',reject='requires')

    def test_incar_system_semicolon_text_preserves_following_assignments(self):
        path,_=self.case('system-title')
        incar=path/'INCAR'
        incar.write_text('SYSTEM = selected optical probe; no exporter patch; LPEAD=F\n'+incar.read_text())
        self.run_case(path,'3')
        path,_=self.case('bad-fragment')
        incar=path/'INCAR';incar.write_text(incar.read_text()+'LPEAD=F; invalid fragment\n')
        self.run_case(path,'3',reject='assignment')
        path,_=self.case('system-branch')
        incar=path/'INCAR';incar.write_text('SYSTEM=descriptive; title; LPEAD=T\nLOPTICS=T; LNABLA=F\n')
        self.run_case(path,'3',reject='LPEAD')

    def test_band_selector_overrides_preserve_historical_order(self):
        for name,flags in [('modern',['--bands','3']),('legacy',['-is','3'])]:
            path,full=self.case(name)
            meta,rows=self.run_case(path,'1:2',*flags)
            self.assertEqual(meta['band_ids'],'3')
            np.testing.assert_allclose([float(r['omega_z_A2']) for r in rows],oracle(full,[3]).ravel(),atol=1e-14)
        path,_=self.case('ambiguous')
        self.run_case(path,'1,3','-if','4',reject='comma-list')

    def test_serial_and_mpi_match_noncontiguous_selection(self):
        if not self.has_mpi:self.skipTest('MPI unavailable')
        a,_=self.case('serial');b,_=self.case('mpi')
        ma,ra=self.run_case(a,'1,3')
        mb,rb=self.run_case(b,'1,3',mpi=True)
        self.assertEqual(ra,rb)
        self.assertEqual(ma['band_ids'],mb['band_ids'])


if __name__=='__main__':unittest.main()
