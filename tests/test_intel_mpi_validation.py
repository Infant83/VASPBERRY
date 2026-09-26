"""Failure guards for the real-input Intel MPI numerical release check."""
import copy
import importlib.util
import math
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("intel_mpi_check", ROOT / "tests/run_intel_mpi_validation.py")
check = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(check)


class IntelMPIValidationTests(unittest.TestCase):
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
                "x,y\n1,2\n# result_status=PASS\n")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "pairs.csv"
            path.write_text(text)
            self.assertEqual(check.read_export(path, "PAIRS", "1.2.3"), [{"x": 1., "y": 2.}])
            with self.assertRaisesRegex(ValueError, "vaspberry_version"):
                check.read_export(path, "PAIRS", "9.9.9")
            path.write_text(text.replace("# result_status=PASS\n", ""))
            with self.assertRaisesRegex(ValueError, "result_status"):
                check.read_export(path, "PAIRS", "1.2.3")

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
