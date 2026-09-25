"""Native aliases and full-executable output contracts.

Parser tests reject malformed commands before file I/O. Independent tiny
WAVECAR fixtures check occupation-free pairs and SI velocity expectations;
actual material equivalence runs are recorded separately in runs/.
"""
import re
import os
import time
import shutil
import struct
import cmath
import math
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


@unittest.skipUnless(shutil.which("gfortran"), "gfortran required for native parser")
class NativeCliTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix="vaspberry-native-cli-")
        cls.work = Path(cls.temp.name)
        source = (ROOT / "vaspberry.f").read_text()
        match = re.search(r"(?ims)^ {6}subroutine parse\b.*?^ {6}end subroutine parse\s*$", source)
        if not match:
            raise AssertionError("production parser not found")
        (cls.work / "parse.f").write_text(match.group(0) + "\n")
        driver = (ROOT / "tests/fortran/test_kubo_parser.f90").read_text()
        driver = driver.replace("  write(*,'(A)')trim(foname)", """
  write(*,'(A)')trim(filename),trim(foname),trim(kubo_csv),trim(kubo_pairs)
  write(*,*)nkx,nky,ispinor,icd,ixt,ivel,ikubo,iz,ihf,nini,nmax,nn,kperiod,it,iskp,ine
  write(*,*)iwf,ikwf,ng,imag,nediv,ikubo_bundle
  write(*,*)theta,phi,init_e,fina_e,sigma,rs
""")
        driver = driver.replace("  stop 1", "  print *, 'HELP'\n  stop 0")
        (cls.work / "driver.f90").write_text(driver)
        cls.binary = cls.work / "parser"
        built = subprocess.run(["gfortran", "-cpp", "-O0", "-fcheck=all",
                                "-ffixed-line-length-none", "parse.f", "driver.f90",
                                "-o", str(cls.binary)], cwd=cls.work, capture_output=True, text=True)
        if built.returncode:
            raise AssertionError(built.stdout + built.stderr)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def run_parser(self, arguments, success=True):
        result = subprocess.run([str(self.binary), *arguments], cwd=self.work,
                                capture_output=True, text=True)
        if success:
            self.assertEqual(result.returncode, 0, result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def equivalent(self, new, legacy):
        self.assertEqual(self.run_parser(new).stdout, self.run_parser(legacy).stdout)

    def test_each_task_dispatch_and_argument_order(self):
        for task, legacy, needed in [
            ("chern", [], []), ("z2", ["-z2", "1"], []),
            ("kubo", ["-kubo", "2"], []), ("kubo-line", ["-kubo", "2"], []),
            ("kubo-integral", ["-kubo", "1"], []),
            ("kubo-pairs", ["-kubo", "2"], ["-kubo_pairs", "pairs.csv"]),
            ("optical", ["-cd", "1"], []), ("spectrum", ["-cd", "2"], []),
            ("velocity", ["-vel", "1"], []),
            ("wavefunction", [], ["-wf", "18"]),
        ]:
            with self.subTest(task=task):
                self.equivalent(["--task", task, *needed], [*legacy, *needed])
                self.equivalent([*needed, "--task", task], [*legacy, *needed])
                self.equivalent(["--task", task, *legacy, *needed], [*legacy, *needed])

    def test_mesh_bands_paths_spinor_and_output_aliases(self):
        self.equivalent(["--task", "chern", "--wavecar", "/tmp/a b/WAVECAR",
                         "--mesh", "12,16", "--bands", "1:18", "--spinor", "1",
                         "--output", "sample"],
                        ["-f", "/tmp/a b/WAVECAR", "-kx", "12", "-ky", "16",
                         "-ii", "1", "-if", "18", "-s", "1", "-o", "sample"])
        self.equivalent(["--bands", "18"], ["-is", "18"])

    def test_curvature_bundle_and_pairs_aliases(self):
        self.equivalent(["--task", "kubo", "--bands", "1:18", "--bundle", "1",
                         "--curvature-csv", "/tmp/a b/curvature.csv"],
                        ["-kubo", "2", "-ii", "1", "-if", "18", "-kubo_bundle", "1",
                         "-kubo_csv", "/tmp/a b/curvature.csv"])
        self.equivalent(["--pairs-csv", "pairs.csv", "--task", "kubo-pairs"],
                        ["-kubo", "2", "-kubo_pairs", "pairs.csv"])

    def test_wavefunction_and_angle_aliases(self):
        self.equivalent(["--task", "wavefunction", "--wavefunction-band", "18",
                         "--kpoint", "2", "--real-grid", "8,10,12", "--imaginary", "1"],
                        ["-wf", "18", "-k", "2", "-ng", "8,10,12", "-im", "1"])
        self.equivalent(["--task", "optical", "--theta", "-45.5", "--phi", "1e2"],
                        ["-cd", "1", "-theta", "-45.5", "-phi", "1e2"])

    def test_task_conflicts_rejected_in_both_orders(self):
        for task, flags in [("chern", ["-kubo", "2"]), ("kubo", ["-kubo", "1"]),
                            ("z2", ["-cd", "1"]), ("optical", ["-vel", "1"]),
                            ("wavefunction", ["-wf", "18", "-z2", "1"]),
                            ("kubo", ["-wf", "18"]), ("chern", ["-ixt", "1"]),
                            ("chern", ["-t", "1"])]:
            for args in [["--task", task, *flags], [*flags, "--task", task]]:
                with self.subTest(args=args):
                    result = self.run_parser(args, False)
                    self.assertIn("conflicts with legacy task", result.stderr)

    def test_missing_requirements_and_existing_csv_guards(self):
        for args in [["--task", "wavefunction"], ["--task", "kubo-pairs"],
                     ["--task", "kubo", "--bundle", "1"],
                     ["--task", "chern", "--curvature-csv", "c.csv"],
                     ["--task", "kubo-pairs", "--pairs-csv", "p.csv", "--bands", "1:2"],
                     ["--task", "kubo", "--wavecar", "same", "--curvature-csv", "same"],
                     ["--task", "kubo-pairs", "--pairs-csv", "p.csv", "--curvature-csv", "c.csv"]]:
            with self.subTest(args=args):
                self.run_parser(args, False)

    def test_bad_task_and_malformed_values_rejected(self):
        for args in [["--task", "unknown"], ["--task", "chern", "--task", "chern"],
                     ["--task"], ["--mesh", "--task", "chern"], ["--unknown", "1"],
                     ["--task", "chern", "-typo", "1"],
                     ["-typo", "1", "--task", "chern"]]:
            self.run_parser(args, False)
        for option, bad_values in {
            "--mesh": ["0,2", "-1,2", "2", "2,3,4", "2,,3", "2,", ",2", "2/3", "2,3x", "2, 3"],
            "--bands": ["0", "3:2", "1:", ":2", "1::2", "1:2:3", "1,2", "1/2", "1x"],
            "--real-grid": ["1,2", "1,2,3,4", "1,,2", "1,0,2"],
            "--spinor": ["0", "3", "1,2", "1/", "abc"],
            "--wavefunction-band": ["0", "-1", "1.2"],
            "--kpoint": ["0", "-1"], "--bundle": ["2", "0,1"],
            "--imaginary": ["2"], "--theta": ["abc", "1,2", "1/", "NaN", "Inf", "1e999"],
        }.items():
            for value in bad_values:
                with self.subTest(option=option, value=value):
                    self.run_parser([option, value], False)

    def test_serial_build_replaces_newer_legacy_binary_with_alias(self):
        project = self.work / "build-contract"
        project.mkdir()
        shutil.copy2(ROOT / "Makefile", project / "Makefile")
        (project / "vaspberry.f").write_text("      program fixture\n      end\n")
        command = ["make", "serial", "FC=gfortran", "GNU_FLAGS=-O0", "GNU_LIBS="]
        for attempt in range(2):
            result = subprocess.run(command, cwd=project, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            canonical = project / "build/vaspberry"
            compatibility = project / "build/vaspberry-gfortran"
            self.assertTrue(canonical.is_file())
            self.assertTrue(compatibility.is_symlink())
            self.assertEqual(compatibility.resolve(), canonical.resolve())
            if attempt == 0:
                compatibility.unlink()
                compatibility.write_text("older build under the former canonical name")
                future = time.time() + 600
                os.utime(compatibility, (future, future))

    def test_standalone_help_is_success(self):
        self.assertIn("HELP", self.run_parser(["--help"]).stdout)
        self.assertIn("HELP", self.run_parser(["-h"]).stdout)


@unittest.skipUnless(shutil.which("gfortran"), "gfortran required for native execution")
class NativeExecutionTests(unittest.TestCase):
    """Run the complete executable on independent, tiny WAVECAR records.

    Pair vertices do not use occupations, unlike automatic occupied-subspace
    selection. These are format/dispatch fixtures, not material calculations.
    """

    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix="vaspberry-native-pair-occupations-")
        cls.work = Path(cls.temp.name)
        cls.binary = cls.work / "vaspberry"
        built = subprocess.run([
            "gfortran", "-cpp", "-O0", "-fcheck=all", "-ffixed-line-length-none",
            "-fallow-argument-mismatch", str(ROOT / "vaspberry.f"),
            "-llapack", "-lblas", "-o", str(cls.binary),
        ], capture_output=True, text=True)
        if built.returncode:
            raise AssertionError(built.stdout + built.stderr)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def write_fixture(self, path, variable_occupations, spinor_components=2):
        stride, nk, nb = 512, 2, 2
        data = bytearray(stride * (2 + nk * (nb + 1)))
        struct.pack_into("<3d", data, 0, stride, 1, 45200)
        struct.pack_into("<12d", data, stride, nk, nb, 120.,
                         1., 0., 0., 0., 1., 0., 0., 0., 1.)
        for ik, point in enumerate([(0., 0., 0.), (.25, .25, 0.)]):
            # Independently enumerate the scalar plane waves; double for SOC.
            count = sum(sum((2 * math.pi * (point[j] + g[j])) ** 2 for j in range(3))
                        / .262465831 < 120.
                        for g in [(x, y, z) for x in range(-2, 3)
                                  for y in range(-2, 3) for z in range(-2, 3)])
            npw = spinor_components * count
            record = 2 + ik * (nb + 1)
            occ = .1 if ik == 1 and variable_occupations else 1.
            struct.pack_into("<10d", data, record * stride,
                             npw, *point, -1., 0., occ, 1., 0., 0.)
            for band in range(nb):
                # Two orthonormal Fourier states over the retained coefficients.
                for ipw in range(npw):
                    coeff = cmath.exp(2j * math.pi * band * ipw / npw) / math.sqrt(npw)
                    struct.pack_into("<2f", data, (record + 1 + band) * stride + ipw * 8,
                                     coeff.real, coeff.imag)
        path.write_bytes(data)

    def test_pairs_ignore_source_occupations_but_subspace_guard_remains(self):
        outputs = []
        for variable in (False, True):
            case = self.work / ("variable" if variable else "uniform")
            case.mkdir()
            self.write_fixture(case / "WAVECAR", variable)
            result = subprocess.run([
                str(self.binary), "--task", "kubo-pairs", "--wavecar", "WAVECAR",
                "--spinor", "2", "--pairs-csv", "PAIRS.csv",
            ], cwd=case, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn("Source occupations: not used for pair export", result.stdout)
            outputs.append((case / "PAIRS.csv").read_bytes())
            self.assertTrue(outputs[-1].endswith(b"# result_status=PASS\n"))
        self.assertEqual(outputs[0], outputs[1])
        # Pair independence must not loosen the occupied-subspace validation.
        result = subprocess.run([
            str(self.binary), "--task", "chern", "--mesh", "2,1", "--bands", "1",
        ], cwd=self.work / "variable", capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("inconsistent occupations", result.stderr)

    def test_velocity_matches_independent_si_plane_wave_expectation(self):
        for spinor in (1, 2):
            with self.subTest(spinor=spinor):
                self.check_velocity(spinor)

    def check_velocity(self, spinor):
        case = self.work / f"velocity-{spinor}"
        case.mkdir()
        self.write_fixture(case / "WAVECAR", False, spinor_components=spinor)
        result = subprocess.run([
            str(self.binary), "--task", "velocity", "--bands", "1",
            "--spinor", str(spinor), "--mesh", "2,1", "-kp", "1",
        ], cwd=case, capture_output=True, text=True)
        # Full executable is compiled with -fcheck=all: this also guards the
        # historical extrema calculation reading ik=nk+1 after the k loop.
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        output = (case / "VEL_EXPT.dat").read_text()
        rows = [list(map(float, line.split())) for line in output.splitlines()
                if line.strip() and not line.startswith("#")]
        self.assertEqual(len(rows), 2)
        # At k=(1/4,1/4,0), the equally weighted G are 0,-x,-y. Thus
        # <kx>=<ky>=-pi/6 A^-1. Use SI constants, independent of the
        # production conversion through electron rest energy in eV.
        hbar_si, mass_si = 1.054571817e-34, 9.1093837015e-31
        amplitude = struct.unpack("<f", struct.pack("<f", 1 / math.sqrt(3 * spinor)))[0]
        expected = hbar_si / mass_si * 1.e10 * (-math.pi / 6) * (3 * spinor * amplitude**2)
        for row in rows:
            target = expected if abs(row[5] - .25) < 1.e-6 else 0.
            self.assertAlmostEqual(row[3], target, delta=5.e-2)
            self.assertAlmostEqual(row[4], target, delta=5.e-2)
        self.assertGreater(abs(expected), 1.e4)
        for axis in ("x", "y"):
            min_line = next(line for line in output.splitlines()
                            if f"MINVAL of VEL_EXPT <v_{axis}>" in line)
            max_line = next(line for line in output.splitlines()
                            if f"MAXVAL of VEL_EXPT <v_{axis}>" in line)
            self.assertAlmostEqual(float(min_line.split()[-1]), expected, delta=5.e-2)
            self.assertAlmostEqual(float(max_line.split()[-1]), 0., delta=1.e-6)


if __name__ == "__main__":
    unittest.main()
