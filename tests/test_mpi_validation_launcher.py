"""Signal policy must preserve argv, exit codes and strict rejection guards."""
import argparse
import importlib.util
import json
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

from mpi_validation_launcher import mpi_validation_command
import run_intel_mpi_validation as intel

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('signal_policy_spin_runner',
    ROOT/'.github/scripts/validate_native_spin.py')
spin = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(spin)


class MPILaunchPolicyTests(unittest.TestCase):
    def test_opt_in_changes_only_sigpipe_and_keeps_real_status(self):
        command = ['/bin/sh', '-c', 'kill -PIPE $$; exit 17']
        self.assertEqual(subprocess.run(mpi_validation_command(command)).returncode,
                         -signal.SIGPIPE)
        self.assertEqual(subprocess.run(mpi_validation_command(command, ignore_sigpipe=True)).returncode,17)
        for code in (0,1,7):
            command = ['/bin/sh','-c',f'exit {code}']
            self.assertEqual(subprocess.run(mpi_validation_command(command,ignore_sigpipe=True)).returncode,code)
        term = ['/bin/sh','-c','kill -TERM $$']
        self.assertEqual(subprocess.run(mpi_validation_command(term,ignore_sigpipe=True)).returncode,
                         -signal.SIGTERM)

    def test_arguments_are_not_reparsed_as_shell_code(self):
        values=['path with spaces', 'a"b\'c', '$(printf injected)', '`printf injected`', '$HOME', '; exit 99']
        command=[sys.executable,'-c','import json,sys; print(json.dumps(sys.argv[1:]))',*values]
        saved=command.copy()
        result=subprocess.run(mpi_validation_command(command,ignore_sigpipe=True),
                              capture_output=True,text=True,check=True)
        self.assertEqual(json.loads(result.stdout),values)
        self.assertEqual(command,saved)
        self.assertEqual(mpi_validation_command(command),command)

    def test_intel_guard_keeps_exit_diagnostic_and_completion_checks(self):
        for code in (1,7):
            with tempfile.TemporaryDirectory() as directory:
                case=Path(directory)/'case'
                command=mpi_validation_command([sys.executable,'-c',
                    f'import sys; sys.stderr.write("expected diagnostic"); sys.exit({code})'],ignore_sigpipe=True)
                if code==1:
                    intel.run_recorded(command,case,expected_success=False,diagnostic='expected diagnostic')
                    self.assertEqual(json.loads((case/'command.json').read_text())['status'],'EXPECTED_REJECTION')
                else:
                    with self.assertRaisesRegex(ValueError,'expected controlled exit 1'):
                        intel.run_recorded(command,case,expected_success=False,diagnostic='expected diagnostic')
                    self.assertEqual(json.loads((case/'command.json').read_text())['status'],'FAIL')
                self.assertEqual(json.loads((case/'command.json').read_text())['command'],command)
        with tempfile.TemporaryDirectory() as directory:
            case=Path(directory)/'case'
            with self.assertRaisesRegex(ValueError,'expected diagnostic'):
                intel.run_recorded(mpi_validation_command(['/bin/sh','-c','exit 1'],ignore_sigpipe=True),
                                   case,expected_success=False,diagnostic='expected diagnostic')
        with tempfile.TemporaryDirectory() as directory:
            case=Path(directory)/'case'
            command=mpi_validation_command([sys.executable,'-c',
                'from pathlib import Path; import sys; Path("BAD.csv").write_text("# result_status=PASS"); sys.exit(1)'],
                ignore_sigpipe=True)
            with self.assertRaisesRegex(ValueError,'completed CSV'):
                intel.run_recorded(command,case,expected_success=False)

    def test_spin_rejection_requires_exact_one_even_with_matching_diagnostic(self):
        for code in (0,1,7):
            with tempfile.TemporaryDirectory() as directory:
                case=Path(directory)
                command=mpi_validation_command([sys.executable,'-c',
                    f'import sys; sys.stderr.write("bad gap"); sys.exit({code})'],ignore_sigpipe=True)
                if code==1:
                    spin.invoke(command,case,expected_error='bad gap')
                else:
                    with self.assertRaisesRegex(AssertionError,'expected rejection'):
                        spin.invoke(command,case,expected_error='bad gap')
                    self.assertEqual(json.loads((case/'command.json').read_text())['status'],'FAIL')
        with tempfile.TemporaryDirectory() as directory:
            case=Path(directory)
            command=mpi_validation_command(['/bin/sh','-c','echo "bad gap" >&2; kill -TERM $$'],ignore_sigpipe=True)
            with self.assertRaisesRegex(AssertionError,'expected rejection'):
                spin.invoke(command,case,expected_error='bad gap')
            self.assertEqual(json.loads((case/'command.json').read_text())['exit_code'],-signal.SIGTERM)

    def test_spin_wraps_only_mpi_and_records_explicit_policy(self):
        for enabled in (False,True):
            with tempfile.TemporaryDirectory() as directory:
                args=argparse.Namespace(serial=sys.executable,mpi=sys.executable,
                    mpiexec=sys.executable,mpi_arg=['--launcher-argument'],timeout=120,
                    mpi_ignore_sigpipe=enabled)
                calls=[]
                with mock.patch.object(spin,'invoke',side_effect=lambda command,*a:calls.append(command)):
                    with self.assertRaises(FileNotFoundError):
                        spin.validate(args,Path(directory))
                binary=str(Path(sys.executable).resolve())
                self.assertEqual(calls[0][0],binary)
                raw=[binary,'--launcher-argument','-n','2',*calls[0]]
                self.assertEqual(calls[1],mpi_validation_command(raw,ignore_sigpipe=enabled))
                environment=json.loads((Path(directory)/'environment.json').read_text())
                self.assertEqual(environment['mpi_signal_policy'],
                                 'ignore_sigpipe' if enabled else 'default')


if __name__=='__main__':
    unittest.main()
