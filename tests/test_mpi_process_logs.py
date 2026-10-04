"""Direct MPI process logs preserve strict rejection checks and literal argv."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

from mpi_validation_launcher import mpi_payload_command, mpi_process_logs, mpi_validation_command
import run_intel_mpi_validation as intel

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('process_logs_spin_runner',
    ROOT / '.github/scripts/validate_native_spin.py')
spin = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(spin)


class MPIProcessLogTests(unittest.TestCase):
    def invoke(self, runner, command, case, diagnostic='required diagnostic'):
        if runner == 'main':
            intel.run_recorded(command, case, expected_success=False, diagnostic=diagnostic)
        else:
            case.mkdir()
            spin.invoke(command, case, expected_error=diagnostic)

    def command(self, program):
        return mpi_payload_command([sys.executable, '-c', program], rank_logs=True)

    def test_literal_argv_default_unchanged_and_real_exit_status(self):
        values = ['path with spaces', 'a"b\'c', '$(printf injected)', '`printf injected`', '$HOME', '; exit 99']
        for code in (0, 1, 7):
            with self.subTest(code=code), tempfile.TemporaryDirectory() as folder:
                command = [sys.executable, '-c',
                    f'import json,sys; print(json.dumps(sys.argv[1:])); print("stderr marker",file=sys.stderr); sys.exit({code})',
                    *values]
                saved = command.copy()
                result = subprocess.run(mpi_payload_command(command, rank_logs=True), cwd=folder,
                                        capture_output=True, text=True)
                self.assertEqual(result.returncode, code)
                self.assertEqual((result.stdout, result.stderr), ('', ''))
                self.assertEqual(json.loads(mpi_process_logs(Path(folder), stream='stdout')[0].read_text()), values)
                self.assertIn('stderr marker', mpi_process_logs(Path(folder), stream='stderr')[0].read_text())
                self.assertEqual(command, saved)
                self.assertEqual(mpi_payload_command(command), command)

    def test_both_drivers_require_real_exit_one_and_actual_diagnostic(self):
        for runner in ('main', 'spin'):
            for code, diagnostic in ((1, 'required diagnostic'), (7, 'required diagnostic'), (1, 'wrong text')):
                with self.subTest(runner=runner, code=code, diagnostic=diagnostic), tempfile.TemporaryDirectory() as folder:
                    case = Path(folder) / 'case'
                    command = self.command(f'import sys; print({diagnostic!r},file=sys.stderr); sys.exit({code})')
                    if code == 1 and diagnostic == 'required diagnostic':
                        self.invoke(runner, command, case)
                        receipt = json.loads((case / 'command.json').read_text())
                        self.assertEqual(receipt['status'], 'EXPECTED_REJECTION')
                        self.assertEqual(receipt['command'], command)
                        self.assertEqual(len(receipt['mpi_process_logs']), 2)
                        for name, record in receipt['mpi_process_logs'].items():
                            data = (case / name).read_bytes()
                            self.assertEqual(record, {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()})
                        self.assertEqual((case / 'stderr.log').read_text(), '')
                    else:
                        with self.assertRaises((ValueError, AssertionError)):
                            self.invoke(runner, command, case)
                        self.assertEqual(json.loads((case / 'command.json').read_text())['status'], 'FAIL')

    def test_completed_results_still_rejected_by_both_drivers(self):
        for runner in ('main', 'spin'):
            with self.subTest(runner=runner), tempfile.TemporaryDirectory() as folder:
                case = Path(folder) / 'case'
                command = self.command('from pathlib import Path; import sys; '
                    'Path("SPIN_BAD.csv").write_text("# result_status=PASS"); '
                    'print("required diagnostic",file=sys.stderr); sys.exit(1)')
                with self.assertRaises((ValueError, AssertionError)):
                    self.invoke(runner, command, case)

    def test_spin_still_requires_stderr_not_stdout_for_diagnostic(self):
        with tempfile.TemporaryDirectory() as folder:
            command = self.command('import sys; print("required diagnostic"); sys.exit(1)')
            with self.assertRaises(AssertionError):
                self.invoke('spin', command, Path(folder) / 'case')

    def test_spin_wraps_only_mpi_and_records_output_policy(self):
        for enabled in (False, True):
            with self.subTest(enabled=enabled), tempfile.TemporaryDirectory() as folder:
                args = argparse.Namespace(serial=sys.executable, mpi=sys.executable,
                    mpiexec=sys.executable, mpi_arg=['--launcher-argument'], timeout=120,
                    mpi_ignore_sigpipe=False, mpi_rank_logs=enabled)
                calls = []
                with mock.patch.object(spin, 'invoke', side_effect=lambda command, *a: calls.append(command)):
                    with self.assertRaises(FileNotFoundError):
                        spin.validate(args, Path(folder))
                binary = str(Path(sys.executable).resolve())
                self.assertEqual(calls[0][0], binary)
                self.assertEqual(calls[1], mpi_validation_command([binary, '--launcher-argument', '-n', '2',
                    *mpi_payload_command(calls[0], rank_logs=enabled)]))
                environment = json.loads((Path(folder) / 'environment.json').read_text())
                self.assertEqual(environment['mpi_output_policy'],
                                 'direct_process_files' if enabled else 'launcher_streams')


if __name__ == '__main__':
    unittest.main()
