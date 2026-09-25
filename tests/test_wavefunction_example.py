"""Exercise complete Fortran wavefunction paths with actual generated records."""
import importlib.util
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "wavefunction_example", ROOT / "validation/models/wavefunction/run.py")
EXAMPLE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(EXAMPLE)


@unittest.skipUnless(shutil.which("gfortran"), "gfortran required for the production bounds regression")
class WavefunctionBoundsRegression(unittest.TestCase):
    def test_complete_current_and_historical_scalar_gamma_paths(self):
        with tempfile.TemporaryDirectory(prefix="vaspberry-wavefunction-") as temporary:
            directory = Path(temporary)
            for source in ("vaspberry.f", "vaspberry_gfortran_serial.f"):
                with self.subTest(source=source):
                    case = directory / Path(source).stem
                    case.mkdir()
                    binary = case / "vaspberry"
                    # A default integer value of 1 makes the old uninitialized
                    # ifinish skip its atom-count loop; bounds checks detect the
                    # old npmax-versus-ncnt vector assignment independently.
                    built = subprocess.run([
                        "gfortran", "-cpp", "-O0", "-fcheck=all", "-finit-integer=1",
                        "-ffixed-line-length-none", "-fallow-argument-mismatch",
                        str(ROOT / source), "-o", str(binary), "-llapack", "-lblas",
                    ], capture_output=True, text=True)
                    self.assertEqual(built.returncode, 0, built.stderr)
                    before = EXAMPLE.fixture(case)
                    result = subprocess.run([
                        str(binary), "-f", "WAVECAR.synthetic", "-s", "1",
                        "-kx", "1", "-ky", "1", "-wf", "1", "-k", "1",
                        "-ng", "8,8,10", "-im", "1",
                    ], cwd=case, capture_output=True, text=True, timeout=30)
                    self.assertEqual(result.returncode, 0, result.stderr)
                    real_file = case / "PARCHG-W-K001-E001-SPIN1"
                    imag_file = case / "PARCHG-W-K001-E001-IM-SPIN1"
                    for path in (real_file, imag_file):
                        lines = path.read_text().splitlines()
                        self.assertEqual(lines[5].strip(), "H")
                        self.assertEqual(lines[6].strip(), "1")
                        self.assertEqual(lines[7].strip(), "Direct")
                    actual = (EXAMPLE.read_grid(real_file, (8, 8, 10)) +
                              1j*EXAMPLE.read_grid(imag_file, (8, 8, 10))) / np.sqrt((2*np.pi)**3)
                    axis = np.arange(8)/8
                    expected = (1+np.exp(2j*np.pi*axis)[None, :]+
                                np.exp(2j*np.pi*axis)[:, None])/np.sqrt(3)
                    np.testing.assert_allclose(actual, np.broadcast_to(expected, actual.shape),
                                               atol=2e-6, rtol=0)
                    for name, digest in before.items():
                        self.assertEqual(EXAMPLE.sha(case/name), digest)


if __name__ == "__main__":
    unittest.main()
