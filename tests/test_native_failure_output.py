"""Check production fatal I/O before an abrupt MPI-style process exit.

The MPI_ABORT test double uses POSIX _exit so language-runtime finalization
cannot hide a missing FLUSH. Actual MPI exit/diagnostic contracts remain covered
by the full-executable validation drivers.
"""
import os
import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


@unittest.skipUnless(os.name == 'posix' and shutil.which('gfortran'),
                     'POSIX and gfortran required for abrupt-exit I/O regression')
class NativeFailureOutputTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix='vaspberry-fatal-output-')
        cls.work = Path(cls.temp.name)
        cls.addClassCleanup(cls.temp.cleanup)
        source = (ROOT / 'vaspberry.f').read_text()
        match = re.search(r'(?ims)^ {6}subroutine vaspberry_fail\b.*?'
                          r'^ {6}end subroutine vaspberry_fail', source)
        if match is None:
            raise AssertionError('production fatal routine not found')
        (cls.work / 'fail.f').write_text(match.group(0) + '\n')
        (cls.work / 'mpif.h').write_text('      integer MPI_COMM_WORLD\n'
                                        '      parameter (MPI_COMM_WORLD=0)\n')
        (cls.work / 'driver.f90').write_text('''program buffered_fatal
 implicit none
 character(20) mode
 call get_command_argument(1,mode)
 if(trim(mode)=='closed')then
   close(0)
   close(6)
 else
   write(0,'(A)') 'REQUIRED STDERR BEFORE ABORT'
   write(6,'(A)') 'REQUIRED STDOUT BEFORE ABORT'
 endif
 call vaspberry_fail
end program
subroutine MPI_ABORT(comm, code, ierr)
 use iso_c_binding, only: c_int
 implicit none
 integer comm, code, ierr
 interface
   subroutine immediate_exit(status) bind(C,name="_exit")
     import c_int
     integer(c_int),value :: status
   end subroutine
 end interface
 call immediate_exit(int(code,c_int))
end subroutine
''')
        cls.binary = cls.work / 'fatal-output'
        result = subprocess.run(['gfortran', '-cpp', '-DMPI_USE',
                                 '-ffixed-line-length-none', '-O0', '-I.',
                                 'fail.f', 'driver.f90', '-o', str(cls.binary)],
                                cwd=cls.work, capture_output=True, text=True)
        if result.returncode:
            raise AssertionError(result.stdout + result.stderr)

    def run_fatal(self, mode):
        folder = self.work / mode
        folder.mkdir()
        env = os.environ.copy()
        env['GFORTRAN_UNBUFFERED_PRECONNECTED'] = 'n'
        env['GFORTRAN_UNBUFFERED_ALL'] = 'n'
        with (folder / 'stdout.log').open('wb') as stdout, \
             (folder / 'stderr.log').open('wb') as stderr:
            result = subprocess.run([str(self.binary), mode], cwd=folder,
                                    stdout=stdout, stderr=stderr, env=env, timeout=10)
        return result.returncode, (folder / 'stdout.log').read_text(), \
            (folder / 'stderr.log').read_text()

    def test_buffered_diagnostics_survive_abrupt_exit_with_file_redirection(self):
        code, stdout, stderr = self.run_fatal('buffered')
        self.assertEqual(code, 1)
        self.assertIn('REQUIRED STDOUT BEFORE ABORT', stdout)
        self.assertIn('REQUIRED STDERR BEFORE ABORT', stderr)

    def test_unavailable_output_units_do_not_prevent_abort_exit_one(self):
        code, _, _ = self.run_fatal('closed')
        self.assertEqual(code, 1)


if __name__ == '__main__':
    unittest.main()
