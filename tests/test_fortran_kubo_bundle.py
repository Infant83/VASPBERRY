"""The compiled native bundle trace must agree with analytic plane-wave states.

Includes internal degeneracy, a basis rotation in that degenerate subspace,
external-gap rejection, scalar/spinor states, CSV metadata, and MPI reduction.
"""
import csv
import io
import os
import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def routine(source, name, kind="subroutine"):
    prefix = r"real\*8\s+" if kind == "function" else ""
    match = re.search(
        rf"(?ims)^ {{6}}{prefix}{kind}\s+{name}\b.*?"
        rf"^ {{6}}end\s*{kind}(?:[ \t]+{name})?[ \t]*$", source
    )
    if not match:
        raise AssertionError(f"missing production routine {name}")
    return match.group(0) + "\n"


def csv_rows(path):
    text = path.read_text()
    metadata = dict(line[2:].split("=", 1) for line in text.splitlines()
                    if line.startswith("# ") and "=" in line)
    rows = list(csv.DictReader(io.StringIO("\n".join(
        line for line in text.splitlines() if not line.startswith("#")))))
    return metadata, rows


@unittest.skipUnless(shutil.which("gfortran"), "gfortran required for native bundle tests")
class CompiledBundleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix="vaspberry-bundle-")
        cls.work = Path(cls.temp.name)
        source = (ROOT / "vaspberry.f").read_text()
        helpers = routine(source, "kubo_interband_term", "function")
        helpers += "\n".join(routine(source, name) for name in (
            "kubo_bundle_gap", "check_kubo_bundle_gaps", "kubo_bundle_curvature",
            "write_kubo_bundle_csv", "kubo_berry_curvature", "mpi_job_distribution_chain",
        ))
        (cls.work / "production.f").write_text(helpers)
        flags = ["-cpp", "-O2", "-fcheck=all", "-ffixed-line-length-none"]
        inputs = [str(cls.work / "production.f"), str(ROOT / "tests/fortran/test_kubo_bundle.F90")]
        cls.run_command(["gfortran", *flags, *inputs, "-o", str(cls.work / "bundle")], cls.work)
        cls.has_mpi = bool(shutil.which("mpifort") and shutil.which("mpiexec"))
        if cls.has_mpi:
            cls.run_command(["mpifort", *flags, "-DMPI_USE", *inputs,
                             "-o", str(cls.work / "bundle-mpi")], cls.work)
        (cls.work / "parse.f").write_text(routine(source, "parse"))
        cls.run_command(["gfortran", "-cpp", "-O2", "-ffixed-line-length-none",
                         str(cls.work / "parse.f"), str(ROOT / "tests/fortran/test_kubo_parser.f90"),
                         "-o", str(cls.work / "parse")], cls.work)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    @classmethod
    def run_command(cls, command, cwd, expect_success=True):
        env = dict(os.environ, OMPI_ALLOW_RUN_AS_ROOT="1", OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1")
        result = subprocess.run(command, cwd=cwd, capture_output=True, text=True, env=env, timeout=60)
        if expect_success and result.returncode:
            raise AssertionError(f"{command}\n{result.stdout}\n{result.stderr}")
        return result

    def execute(self, mode, mpi=False):
        case = self.work / (self._testMethodName + "-" + mode + ("-mpi" if mpi else ""))
        case.mkdir()
        command = [str(self.work / ("bundle-mpi" if mpi else "bundle")), mode]
        if mpi:
            command = ["mpiexec", "-n", "2", *command]
        result = self.run_command(command, case)
        self.assertIn("KUBO_BUNDLE_PASS", result.stdout)
        return case

    def test_bundle_matches_band_sum_internal_pairs_cancel(self):
        self.execute("normal")

    def test_internal_degeneracy_and_rotation_are_finite_and_invariant(self):
        direct = self.execute("degenerate")
        rotated = self.execute("rotated")
        _, rows = csv_rows(direct / "bundle.csv")
        _, rotated_rows = csv_rows(rotated / "bundle.csv")
        np.testing.assert_allclose([float(r["omega_z_A2"]) for r in rows],
                                   [float(r["omega_z_A2"]) for r in rotated_rows], rtol=3e-6)

    def test_bundle_csv_range_external_gaps_and_empty_external_space(self):
        case = self.execute("degenerate")
        meta, rows = csv_rows(case / "bundle.csv")
        expected = {
            "schema": "VASPBERRY_BARE_MOMENTUM_KUBO_BUNDLE_V1",
            "normalization": "STANDARD_MINUS_TWO_IM", "band_min": "1", "band_max": "2",
            "band_rank": "2", "source_nbands": "3", "gap_threshold_eV": "1e-5",
            "result_kind": "ISOLATED_BUNDLE_TRACE", "result_status": "PASS",
            "no_external_states": "false", "internal_transitions": "EXCLUDED_ANALYTICALLY",
        }
        for key, value in expected.items():
            self.assertEqual(meta[key], value)
        self.assertEqual([(int(r["spin"]), int(r["k_index"])) for r in rows],
                         [(s, k) for s in (1, 2) for k in (1, 2, 3)])
        self.assertEqual([float(r["min_external_gap_eV"]) for r in rows], [2, 3, 4]*2)
        all_meta, all_rows = csv_rows(case / "all.csv")
        self.assertEqual(all_meta["no_external_states"], "true")
        self.assertEqual(all_meta["zero_trace_scope"], "TRUNCATED_WAVECAR_BASIS")
        self.assertEqual([float(r["omega_z_A2"]) for r in all_rows], [0, 0, 0])
        self.assertEqual([r["min_external_gap_eV"] for r in all_rows], ["NA"]*3)
        # Opt-in outputs must not overwrite existing results.
        original = (case / "bundle.csv").read_bytes()
        failed = self.run_command([str(self.work / "bundle"), "degenerate"], case, False)
        self.assertNotEqual(failed.returncode, 0)
        self.assertIn("cannot open new Kubo bundle CSV", failed.stderr)
        self.assertEqual((case / "bundle.csv").read_bytes(), original)

    def test_external_gap_rejection_happens_before_csv_for_all_spins(self):
        for mode in ("external-gap", "near-gap", "later-spin", "nonfinite"):
            with self.subTest(mode=mode):
                case = self.work / (self._testMethodName + mode)
                case.mkdir()
                failed = self.run_command([str(self.work / "bundle"), mode], case, False)
                self.assertNotEqual(failed.returncode, 0)
                self.assertIn("Kubo bundle is not isolated", failed.stderr)
                self.assertFalse((case / "bundle.csv").exists())

    def test_mpi_reduces_same_three_k_points_as_serial(self):
        if not self.has_mpi:
            self.skipTest("mpifort and mpiexec required")
        serial = self.execute("degenerate")
        parallel = self.execute("degenerate", mpi=True)
        self.assertEqual((serial / "bundle.csv").read_bytes(), (parallel / "bundle.csv").read_bytes())
        self.assertEqual((serial / "all.csv").read_bytes(), (parallel / "all.csv").read_bytes())

    def test_parser_requires_kubo_only_and_csv_without_changing_default(self):
        invalid = [(["-kubo_bundle", "1"], "needs -kubo_csv"),
                   (["-kubo_bundle", "2"], "must be 0 or 1"),
                   (["-kubo_bundle", "1", "-kubo_csv", "x.csv"], "requires Kubo-only"),
                   (["-kubo_bundle", "1", "-kubo", "2", "-kubo_csv", "x.csv", "-z2", "1"],
                    "requires Kubo-only")]
        for args, message in invalid:
            with self.subTest(args=args):
                result = self.run_command([str(self.work / "parse"), *args], self.work, False)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(message, result.stderr)
        self.run_command([str(self.work / "parse"), "-kubo", "2", "-kubo_bundle", "1",
                          "-ii", "1", "-if", "2", "-kubo_csv", "x.csv"], self.work)


if __name__ == "__main__":
    unittest.main()
