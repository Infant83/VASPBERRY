"""Portable process-boundary checks; these fixtures do not simulate physics."""
import json
from pathlib import Path
import tempfile
import unittest

from run_vasp_spin_producer import run_case
import vasp544_spin_instrument as instrument


class ProducerProcessTests(unittest.TestCase):
    def prepare(self, root, script, executable=True):
        case = root/'case'; case.mkdir()
        for name in ('INCAR', 'POSCAR', 'KPOINTS', 'POTCAR', 'CHGCAR'):
            (case/name).write_text('synthetic input\n')
        binary = root/'fake-producer'; binary.write_text('#!/bin/sh\n'+script+'\n')
        binary.chmod(0o700 if executable else 0o600)
        manifest = root/'producer.json'
        manifest.write_text(json.dumps({'schema': 'vaspberry.vasp544-spin-producer', 'version': 1}))
        return case, binary, manifest

    def test_launch_failure_leaves_terminal_record(self):
        with tempfile.TemporaryDirectory() as directory:
            case, binary, manifest = self.prepare(Path(directory), 'exit 0', executable=False)
            with self.assertRaises(OSError): run_case(case, binary, manifest)
            record = json.loads((case/'run.json').read_text())
            self.assertEqual(record['status'], 'FAILED_LAUNCH')
            self.assertIn('finished_utc', record)
            self.assertIn('error', record)

    def test_zero_exit_without_outputs_is_failure_and_cannot_clobber(self):
        with tempfile.TemporaryDirectory() as directory:
            case, binary, manifest = self.prepare(Path(directory), 'exit 0')
            record = run_case(case, binary, manifest)
            self.assertEqual(record['returncode'], 0)
            self.assertEqual(record['status'], 'FAILED_OUTPUT')
            before = (case/'run.json').read_bytes()
            with self.assertRaisesRegex(ValueError, 'already contains output'):
                run_case(case, binary, manifest)
            self.assertEqual((case/'run.json').read_bytes(), before)

    def test_timeout_terminates_process_group_and_preserves_failure(self):
        with tempfile.TemporaryDirectory() as directory:
            case, binary, manifest = self.prepare(Path(directory), 'sleep 10')
            record = run_case(case, binary, manifest, timeout=.01)
            self.assertEqual(record['status'], 'TIMEOUT')
            self.assertIsNotNone(record['returncode'])
            self.assertFalse(record['completed_export'])
            self.assertIn('finished_utc', record)

    def test_instrumenter_rejects_wrong_source_without_modifying_files(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); (root/'src').mkdir()
            for name in instrument.BASE_HASHES:
                (root/'src'/name).write_text('unsupported synthetic source\n')
            before = {p.name: p.read_bytes() for p in (root/'src').iterdir()}
            with self.assertRaises(ValueError): instrument.instrument(root)
            self.assertEqual(before, {p.name: p.read_bytes() for p in (root/'src').iterdir()})
            self.assertFalse((root/'vaspberry-spin-producer.json').exists())


if __name__ == '__main__':
    unittest.main()
