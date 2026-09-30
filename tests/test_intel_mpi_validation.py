"""Failure guards for the real-input Intel MPI numerical release check."""
import copy
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import struct
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("intel_mpi_check", ROOT / "tests/run_intel_mpi_validation.py")
check = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(check)

# Identical stderr captured from hosted ifx/ifort runs on the public MoS2 path.
INTEL_WRAPPED_GAP_STDERR = (
    " *** error - individual Kubo band is not \n"
    " isolated: spin,k_index,band,other_band,min_gap_eV =            1          23\n"
    "          17          18  4.666991826107747E-006\n"
    " *** requires all gaps > 1e-5 eV; explicitly \n"
    " choose a separated --bands FIRST:LAST subspace \n"
    " without --per-band 1; selection was not enlarged\n"
    "1\n"
)


class IntelMPIValidationTests(unittest.TestCase):
    def diagnostic_command(self, stderr, exit_code=1):
        return [sys.executable, "-c", f"import sys; sys.stderr.write({stderr!r}); sys.exit({exit_code})"]

    def test_recorded_rejection_accepts_intel_wrapped_diagnostic(self):
        with tempfile.TemporaryDirectory() as directory:
            for compiler in ("ifx", "ifort"):
                output=Path(directory)/compiler
                check.run_recorded(self.diagnostic_command(INTEL_WRAPPED_GAP_STDERR),output,
                                   expected_success=False,diagnostic="not isolated")
                record=json.loads((output/"command.json").read_text())
                self.assertEqual(record["exit_code"],1)
                self.assertEqual(record["status"],"EXPECTED_REJECTION")
                self.assertEqual((output/"stderr.log").read_text(),INTEL_WRAPPED_GAP_STDERR)

    def test_recorded_rejection_still_requires_the_requested_diagnostic(self):
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory)/"wrong-diagnostic"
            with self.assertRaisesRegex(ValueError,"expected diagnostic"):
                check.run_recorded(self.diagnostic_command("*** error - cannot open WAVECAR\n"),output,
                                   expected_success=False,diagnostic="not isolated")
            self.assertEqual(json.loads((output/"command.json").read_text())["status"],"FAIL")

    def test_wrapped_diagnostic_does_not_override_wrong_exit_status(self):
        with tempfile.TemporaryDirectory() as directory:
            for exit_code, error in ((0,RuntimeError),(2,ValueError)):
                output=Path(directory)/f"exit-{exit_code}"
                with self.subTest(exit_code=exit_code),self.assertRaises(error):
                    check.run_recorded(self.diagnostic_command(INTEL_WRAPPED_GAP_STDERR,exit_code),output,
                                       expected_success=False,diagnostic="not isolated")
                record=json.loads((output/"command.json").read_text())
                self.assertEqual(record["exit_code"],exit_code)
                self.assertEqual(record["status"],"FAIL")

    @classmethod
    def setUpClass(cls):
        cls.pairs, cls.bundle = [], []
        for k in range(1, 49):
            terms = []
            for n in range(1, 33):
                for m in range(n + 1, 33):
                    numerator = (n - m) / 16.0
                    cls.pairs.append(dict(spin=1., k_index=float(k), n_band=float(n), m_band=float(m),
                                          kx_frac=k / 48., ky_frac=0., kz_frac=0.,
                                          energy_n_eV=float(n), energy_m_eV=float(m), gap_eV=float(m - n),
                                          numerator_xy_eV2_A2=numerator))
                    if n <= 18 < m:
                        terms.append(numerator / (n - m)**2)
            cls.bundle.append(dict(spin=1., k_index=float(k), kx_frac=k / 48., ky_frac=0., kz_frac=0.,
                                   omega_z_A2=math.fsum(terms), min_external_gap_eV=1.))

    def test_complete_pair_sum_recovers_occupied_bundle(self):
        self.assertEqual(check.pair_bundle_check(self.bundle, self.pairs)["max_absolute_difference"], 0.)

    def test_missing_pair_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "incomplete"):
            check.pair_bundle_check(self.bundle, self.pairs[:-1])

    def test_duplicate_pair_cannot_replace_missing_pair(self):
        with self.assertRaisesRegex(ValueError, "duplicate"):
            check.pair_bundle_check(self.bundle, self.pairs[:-1] + self.pairs[:1])

    def test_duplicate_bundle_k_point_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "distinct"):
            check.pair_bundle_check(self.bundle[:-1] + self.bundle[:1], self.pairs)

    def test_wrong_pair_coordinate_is_rejected(self):
        pairs = [dict(self.pairs[0], kx_frac=.123), *self.pairs[1:]]
        with self.assertRaisesRegex(ValueError, "coordinates"):
            check.pair_bundle_check(self.bundle, pairs)

    def test_wrong_gap_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "energy/gap"):
            check.pair_bundle_check(self.bundle, [dict(self.pairs[0], gap_eV=999.), *self.pairs[1:]])

    def test_wrong_normalization_is_rejected(self):
        bundle = [dict(row, omega_z_A2=2 * row["omega_z_A2"]) for row in self.bundle]
        with self.assertRaisesRegex(ValueError, "omega_z_A2 differs"):
            check.pair_bundle_check(bundle, self.pairs)

    def test_comparison_rejects_nan_and_large_finite_error(self):
        for invalid in (float("nan"), 1.):
            with self.subTest(value=invalid), self.assertRaises(ValueError):
                check.compare_rows([{"value": 0.}], [{"value": invalid}], label="guard")

    def test_read_export_requires_final_pass_and_version(self):
        text = ("# schema=VASPBERRY_BARE_MOMENTUM_KUBO_PAIRS_V1\n"
                "# normalization=STANDARD_MINUS_TWO_IM\n# vaspberry_version=1.2.3\n"
                "# source_nkpoints=1\n# source_nspin=1\n# source_nbands=2\n# expected_rows=1\n"
                "spin,k_index,n_band,m_band\n1,1,1,2\n# result_status=PASS\n")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "pairs.csv"
            path.write_text(text)
            self.assertEqual(check.read_export(path, "PAIRS", "1.2.3"),
                             [dict(spin=1.,k_index=1.,n_band=1.,m_band=2.)])
            with self.assertRaisesRegex(ValueError, "vaspberry_version"):
                check.read_export(path, "PAIRS", "9.9.9")
            path.write_text(text.replace("# result_status=PASS\n", ""))
            with self.assertRaisesRegex(ValueError, "result_status"):
                check.read_export(path, "PAIRS", "1.2.3")
            path.write_text(text + "# data appended after completion\n")
            with self.assertRaisesRegex(ValueError, "terminal"):
                check.read_export(path, "PAIRS", "1.2.3")

    def test_synthetic_layout_fixtures_differ_only_in_logical_stride(self):
        with tempfile.TemporaryDirectory() as directory:
            paths = [Path(directory)/name for name in ("byte", "word4")]
            for path, word4 in zip(paths, (False, True)):
                check.write_synthetic_wavecar(path, word4=word4)
            byte, word4 = [path.read_bytes() for path in paths]
            self.assertEqual(len(byte), 18432)
            self.assertEqual(struct.unpack_from("<d", byte)[0],1024.)
            self.assertEqual(struct.unpack_from("<d", word4)[0],256.)
            self.assertEqual(byte[8:],word4[8:])
            for corruption in ("truncated", "nan", "inf", "dimensions"):
                path=Path(directory)/corruption
                check.write_synthetic_wavecar(path,corruption=corruption)
                self.assertNotEqual(path.read_bytes(),byte)
                if corruption=="truncated":
                    self.assertEqual(path.read_bytes(),byte[:-1024])

    def test_harness_help_needs_only_the_python_standard_library(self):
        import sys
        result=subprocess.run([sys.executable,"-S",str(ROOT/"tests/run_intel_mpi_validation.py"),"--help"],
                              capture_output=True,text=True,check=True)
        self.assertIn("--serial",result.stdout)

    def test_make_checks_use_matching_wrapper_and_intel_launcher(self):
        for compiler in ("ifx", "ifort"):
            with self.subTest(compiler=compiler):
                result = subprocess.run(["make", "-n", f"check-{compiler}-mpi"], cwd=ROOT,
                                        text=True, capture_output=True, check=True)
                self.assertIn(f"mpi{compiler} -fpp -O2 -extend-source -assume byterecl -DMPI_USE", result.stdout)
                self.assertIn(f"mpi{compiler} -O2 -o build/test-{compiler}-mpi-runtime", result.stdout)
                self.assertIn(f"mpiexec.hydra  -n 2 build/test-{compiler}-mpi-runtime", result.stdout)
                self.assertIn(f"mpiexec.hydra  -n 2 build/vaspberry-{compiler}-mpi --help", result.stdout)

    def test_make_intel_launcher_override_is_honored(self):
        result = subprocess.run(["make", "-n", "check-ifx-mpi", "INTEL_MPIEXEC=site-mpiexec",
                                 "INTEL_MPIEXEC_FLAGS=-verbose"], cwd=ROOT,
                                text=True, capture_output=True, check=True)
        self.assertIn("site-mpiexec -verbose -n 2 build/test-ifx-mpi-runtime", result.stdout)


if __name__ == "__main__":
    unittest.main()
