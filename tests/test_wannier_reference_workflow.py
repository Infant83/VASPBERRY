"""Bounded workflow tests; synthetic process fixtures are not material results."""
from pathlib import Path
import argparse
import importlib.util
import json
import os
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
path = ROOT/'examples/materials/mnbi2te4-qah/wannier_reference.py'
spec = importlib.util.spec_from_file_location('tested_wannier_reference', path)
wr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(wr)


@unittest.skipUnless(os.name == "posix", "executable fixtures require a POSIX host")
class Workflow(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='vaspberry-wannier-workflow-test-')
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.here = self.root/'example'
        self.package = self.here/'inputs/wannier/operators'
        self.package.mkdir(parents=True)
        toolchain = self.here/'inputs/wannier/toolchain'
        toolchain.mkdir()
        (toolchain/'effective-reader-modern31.patch').write_bytes((path.parent/'inputs/wannier/toolchain/effective-reader-modern31.patch').read_bytes())
        self.ops = self.root/'restored'
        self.ops.mkdir()
        lattice = [[4.336, 0, 0], [-2.168, 3.755086150809, 0], [0, 0, 52.891132533333]]
        meta = dict(status='VALIDATED_FULL_OPERATOR_INPUT', lattice_A=lattice,
                    source_use_ws_distance=True, source_transl_inv=False, num_wann=138)
        for key in ['HH', 'AA']:
            f = self.ops/f'wannier90_{key}_R.dat'
            f.write_text('SYNTHETIC WORKFLOW TEST INPUT, NO PHYSICAL OPERATOR\n'+key+'\n')
            meta[key] = dict(sha256=wr.sha(f))
        wr.save(self.package/'export.json', meta)
        direct = self.here/'reference/direct-dft'
        direct.mkdir(parents=True)
        wr.save(direct/'metadata.json', dict(lattice_A=lattice, vbm_eV=3.2, cbm_eV=3.3,
                                           energy_reference='synthetic test zero'))
        self.patcher1 = patch.object(wr, 'HERE', self.here)
        self.patcher2 = patch.object(wr, 'PACKAGE', self.package)
        self.patcher1.start(); self.patcher2.start()
        self.addCleanup(self.patcher1.stop); self.addCleanup(self.patcher2.stop)
        self.out = self.root/'calculation'
        wr.prepare(argparse.Namespace(output_dir=self.out, chunks=2, base=4, refine=3,
                                      radius=.18, operators_dir=self.ops))

    def fake_binary(self, *, exit_code=0, footer=True, version='3.1.0', components=True):
        program = self.root/'fake_postw90.py'
        program.write_text(f'''#!{sys.executable}
import os,re,sys
from pathlib import Path
for key in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']:
 assert os.environ[key]=='1'
t=Path('wannier90.win').read_text()
get=lambda key:float(re.search(r'^'+key+r'\\s*=\\s*(.*)$',t,re.M).group(1))
lo,hi=get('fermi_energy_min'),get('fermi_energy_max')
mus=[lo,(lo+hi)/2,hi]
rows=[' '.join(f'{{v:.6f}}' for v in [mu,.1,.2,.3]) for mu in mus]
Path('wannier90-ahc-fermiscan.dat').write_text('\\n'.join(rows)+'\\n')
text='Release: {version}\\n'
if {components!r}:
 for mu in mus:
  text+='J0 term : 0.1000 0.2000 0.3000\\nJ1 term : 0.0000 0.0000 0.0000\\nJ2 term : 0.0000 0.0000 0.0000\\n'
if {footer!r}:text+='All done: postw90 exiting\\n'
Path('wannier90.wpout').write_text(text)
sys.exit({exit_code})
''')
        program.chmod(0o755)
        return program

    def run_case(self, *, workers=2, **kwargs):
        return wr.run(argparse.Namespace(output_dir=self.out, postw90=self.fake_binary(**kwargs),
                                         workers=workers, timeout_seconds=15))

    def test_default_gap_range_and_partition_weights(self):
        plan = wr.read_plan(self.out)
        np.testing.assert_allclose(plan['requested_mu_eV'], [3.205, 3.25, 3.295], atol=1e-14)
        self.assertAlmostEqual(sum(p['weight_sum'] for p in plan['parts']), 1)
        self.assertEqual(plan['temperature_K'], 0)
        self.assertEqual(plan['reader_patch']['sha256'], wr.sha(self.here/plan['reader_patch']['source_relative_path']))
        self.assertTrue(all((self.out/p['directory']/'wannier90_AA_R.dat').is_file() for p in plan['parts']))

    def test_success_and_format_parity(self):
        self.assertEqual(self.run_case()['status'], 'FINISHED')
        wr.collect(argparse.Namespace(output_dir=self.out, result_dir=None))
        result = self.out/'results'
        csv_data = np.loadtxt(result/'conductivity.csv', skiprows=1, delimiter=',')
        dat_data = np.loadtxt(result/'conductivity.dat')
        np.testing.assert_array_equal(csv_data, dat_data)
        with np.load(result/'conductivity.npz', allow_pickle=False) as data:
            np.testing.assert_array_equal(csv_data, np.column_stack([data[k] for k in wr.COLUMNS]))
            self.assertIn('full J0+J1+J2', str(data['operator']))
        self.assertAlmostEqual(csv_data[0, -1], .6*52.891132533333*1e-8/(1.602176634e-19**2/6.62607015e-34))
        with self.assertRaises(ValueError):wr.collect(argparse.Namespace(output_dir=self.out, result_dir=None))

    def test_fortran_stop_zero_without_footer_is_failure(self):
        with self.assertRaises(ValueError):self.run_case(footer=False, workers=1)
        r=json.loads((self.out/'part00/run.json').read_text())
        self.assertEqual(r['status'], 'FAILED'); self.assertEqual(r['returncode'], 0)
        with self.assertRaises(ValueError):wr.collect(argparse.Namespace(output_dir=self.out, result_dir=None))
        self.assertFalse((self.out/'results').exists())

    def test_nonzero_with_footer_is_failure(self):
        with self.assertRaises(ValueError):self.run_case(exit_code=9)
        self.assertEqual(json.loads((self.out/'execution.json').read_text())['status'], 'FAILED')

    def test_missing_connection_terms_is_failure(self):
        with self.assertRaises(ValueError):self.run_case(components=False)

    def test_wrong_version_is_failure(self):
        with self.assertRaises(ValueError):self.run_case(version='2.1.0')

    def test_stale_output_rejected_before_execution(self):
        (self.out/'part00/wannier90.wpout').write_text('stale')
        with self.assertRaises(ValueError):self.run_case()
        self.assertFalse((self.out/'execution.json').exists())

    def test_modified_operator_rejected_before_execution(self):
        (self.out/'part00/wannier90_AA_R.dat').write_text('wrong operator')
        with self.assertRaises(ValueError):self.run_case()
        self.assertFalse((self.out/'execution.json').exists())

    def test_modified_output_rejected_before_results(self):
        self.run_case()
        with (self.out/'part00/wannier90-ahc-fermiscan.dat').open('a') as f:f.write('0 0 0 0\n')
        with self.assertRaises(ValueError):wr.collect(argparse.Namespace(output_dir=self.out, result_dir=None))
        self.assertFalse((self.out/'results').exists())

    def test_invalid_workers_do_not_launch(self):
        with self.assertRaises(ValueError):wr.run(argparse.Namespace(output_dir=self.out, postw90=self.fake_binary(), workers=9, timeout_seconds=15))
        self.assertFalse((self.out/'execution.json').exists())


if __name__ == '__main__':
    unittest.main(verbosity=2)
