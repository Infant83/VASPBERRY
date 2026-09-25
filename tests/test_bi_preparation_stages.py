"""Ordinary Bi WAVECAR preparation must preserve the physical spin-run inputs."""
import contextlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'examples/materials/bi-spin-hall/prepare_vasp.py'
spec = importlib.util.spec_from_file_location('bi_prepare', SCRIPT)
prepare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prepare)


def settings(path):
    return {key.strip(): value.strip() for line in path.read_text().splitlines()
            if '=' in line for key, value in [line.split('=', 1)]}


class BiPreparationStages(unittest.TestCase):
    def test_wavecar_preserves_hamiltonian_and_complete_mesh(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            potential = folder / 'POTCAR'
            potential.write_text('test-only potential identity fixture\n')
            charge = folder / 'CHGCAR'
            charge.write_bytes((prepare.HERE / 'inputs/POSCAR').read_bytes())
            identity = json.loads((prepare.HERE / 'inputs/provenance.json').read_text())['potcar_sha256']
            actual_sha = prepare.sha
            def fixture_sha(path):
                return identity if Path(path) == potential else actual_sha(path)
            for stage in ('spin', 'wavecar'):
                args = [str(SCRIPT), '--stage', stage, '--potcar', str(potential),
                        '--charge', str(charge), '--mesh', '6', '8', '--nbands', '48',
                        '--output-dir', str(folder/stage)]
                with mock.patch.object(sys, 'argv', args), \
                        mock.patch.object(prepare, 'sha', side_effect=fixture_sha), \
                        contextlib.redirect_stdout(io.StringIO()):
                    prepare.main()
            spin = settings(folder/'spin/INCAR')
            ordinary = settings(folder/'wavecar/INCAR')
            optional = {'LOPTICS', 'LPEAD', 'LNABLA', 'LBERRY_EXPORT', 'LSPIN_EXPORT'}
            self.assertTrue(optional <= spin.keys())
            self.assertFalse(optional & ordinary.keys())
            self.assertEqual({k:v for k,v in spin.items() if k not in optional|{'SYSTEM'}},
                             {k:v for k,v in ordinary.items() if k != 'SYSTEM'})
            for name in ('POSCAR', 'POTCAR', 'CHGCAR'):
                self.assertEqual((folder/'spin'/name).read_bytes(),
                                 (folder/'wavecar'/name).read_bytes())
            points = (folder/'wavecar/KPOINTS').read_text().splitlines()[1:]
            self.assertEqual(points, (folder/'spin/KPOINTS').read_text().splitlines()[1:])
            self.assertEqual(points[0], '48')
            self.assertEqual(len(points[2:]), 48)
            self.assertEqual(len(set(points[2:])), 48)
            record=json.loads((folder/'wavecar/input_manifest.json').read_text())
            self.assertEqual(record['stage'], 'wavecar')
            self.assertIn('CHGCAR', record['input_sha256'])


if __name__ == '__main__':
    unittest.main()
