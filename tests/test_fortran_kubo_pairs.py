"""Native all-pair export: raw finite vertices, complete schema and MPI parity."""
import csv
import io
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def routine(source, name, kind="subroutine"):
    prefix = r"real\*8\s+" if kind == "function" else ""
    match = re.search(rf"(?ims)^ {{6}}{prefix}{kind}\s+{name}\b.*?"
                      rf"^ {{6}}end\s*{kind}(?:[ \t]+{name})?[ \t]*$", source)
    if not match:
        raise AssertionError(f"missing production routine {name}")
    return match.group(0) + "\n"


def read_csv(path):
    text = path.read_text()
    metadata = dict(line[2:].split("=", 1) for line in text.splitlines()
                    if line.startswith("# ") and "=" in line)
    rows = list(csv.DictReader(io.StringIO("\n".join(
        line for line in text.splitlines() if not line.startswith("#")))))
    return metadata, rows


@unittest.skipUnless(shutil.which("gfortran"), "gfortran required for native pair tests")
class CompiledPairTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix="vaspberry-pairs-")
        cls.work = Path(cls.temp.name)
        source = (ROOT / "vaspberry.f").read_text()
        production = routine(source, "kubo_interband_term", "function")
        production += "\n".join(routine(source, name) for name in (
            "kubo_pair_numerators", "write_kubo_pairs_csv", "kubo_bundle_gap",
            "kubo_bundle_curvature", "mpi_job_distribution_chain",
        ))
        (cls.work / "production.f").write_text(production)
        flags = ["-cpp", "-O2", "-fcheck=all", "-ffixed-line-length-none"]
        inputs = [str(cls.work / "production.f"), str(ROOT / "tests/fortran/test_kubo_pairs.F90")]
        cls.run_command(["gfortran", *flags, *inputs, "-o", str(cls.work / "pairs")], cls.work)
        cls.has_mpi = bool(shutil.which("mpifort") and shutil.which("mpiexec"))
        if cls.has_mpi:
            cls.run_command(["mpifort", *flags, "-DMPI_USE", *inputs, "-o", str(cls.work / "pairs-mpi")], cls.work)
        (cls.work / "parse.f").write_text(routine(source, "parse"))
        cls.run_command(["gfortran", "-cpp", "-O2", "-ffixed-line-length-none", str(cls.work / "parse.f"),
                         str(ROOT / "tests/fortran/test_kubo_parser.f90"), "-o", str(cls.work / "parse")], cls.work)

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
        command = [str(self.work / ("pairs-mpi" if mpi else "pairs")), mode]
        if mpi:
            command = ["mpiexec", "-n", "2", *command]
        result = self.run_command(command, case)
        self.assertIn("KUBO_PAIRS_PASS", result.stdout)
        return case

    def test_all_three_numerators_match_plane_wave_matrix_oracle(self):
        case = self.execute("normal")
        # Independently form the momentum matrices of three Fourier states
        # on the 0,Gx,Gy,Gz plane-wave basis, including single-precision input.
        basis = np.array([[1, 1, 1], [1, 1j, -1], [1, -1, 1], [1, -1j, -1]], dtype=np.complex64)
        unit = (float(np.float32(6.582119569e-16))**4 / float(np.float32(.510998950e6))**2
                * float(np.float32(1e10))**4 * float(np.float32(2.99792458e8))**4)
        for spinor, name, channels in [(1, "pairs-scalar.csv", 2), (2, "pairs-spinor.csv", 1)]:
            metadata, rows = read_csv(case / name)
            u = (basis/np.sqrt(np.float32(4*spinor))).astype(np.complex128)
            momentum = np.array([spinor*u.conj().T @ np.diag(np.eye(4)[axis+1]) @ u for axis in range(3)])
            self.assertEqual(len(rows), 3*3*channels)
            self.assertEqual([(int(r["spin"]), int(r["k_index"]), int(r["n_band"]), int(r["m_band"])) for r in rows],
                             [(s,k,n,m) for s in range(1,channels+1) for k in (1,2,3) for n,m in ((1,2),(1,3),(2,3))])
            for row in rows:
                n,m = int(row["n_band"])-1,int(row["m_band"])-1
                for component,a,b in (("yz",1,2),("zx",2,0),("xy",0,1)):
                    expected = -2*np.imag(momentum[a,n,m]*momentum[b,m,n])*unit
                    self.assertAlmostEqual(float(row[f"numerator_{component}_eV2_A2"]), expected, delta=1e-12)
                self.assertAlmostEqual(float(row["gap_eV"]), abs(float(row["energy_n_eV"])-float(row["energy_m_eV"])))
            self.assertEqual(metadata["source_nspin"], str(channels))
            self.assertEqual(metadata["spinor_components"], str(spinor))

    def test_degenerate_pair_is_finite_nonzero_and_cache_fallback_agrees(self):
        case = self.execute("degenerate")
        meta, rows = read_csv(case / "pairs-scalar.csv")
        zero_gap = [row for row in rows if float(row["gap_eV"]) == 0]
        self.assertEqual(len(zero_gap), 6)
        self.assertTrue(all(abs(float(r["numerator_xy_eV2_A2"])) > 0 for r in zero_gap))
        self.assertEqual(meta["denominator_weighting"], "NONE")
        # The compiled fixture compares cache/fallback bitwise and reconstructs
        # the occupied-bundle trace from the external pairs at all three k.

    def test_complete_metadata_and_no_clobber(self):
        case = self.execute("degenerate")
        p = case / "pairs-scalar.csv"
        meta, rows = read_csv(p)
        expected = {"schema":"VASPBERRY_BARE_MOMENTUM_KUBO_PAIRS_V1",
                    "result_kind":"UNORDERED_INTERBAND_NUMERATORS", "result_status":"PASS",
                    "normalization":"STANDARD_MINUS_TWO_IM", "source_nbands":"3",
                    "source_nkpoints":"3", "pairs_per_k":"3", "expected_rows":"18",
                    "components":"yz,zx,xy", "pair_order":"n_lt_m", "occupation_weighting":"NONE",
                    "gap_definition":"ABS_EN_MINUS_EM", "reciprocal_convention":"2pi"}
        for key, value in expected.items():
            self.assertEqual(meta[key], value)
        self.assertTrue(p.read_text().endswith("# result_status=PASS\n"))
        lattice=np.array([[float(x) for x in meta[f"lattice_A_{i}"].split(",")] for i in (1,2,3)])
        reciprocal=np.array([[float(x) for x in meta[f"reciprocal_inv_A_{i}"].split(",")] for i in (1,2,3)])
        np.testing.assert_allclose(lattice @ reciprocal.T,2*np.pi*np.eye(3),atol=1e-15)
        before=p.read_bytes()
        failed=self.run_command([str(self.work/"pairs"),"degenerate"],case,False)
        self.assertNotEqual(failed.returncode,0)
        self.assertIn("cannot open new Kubo pairs CSV",failed.stderr)
        self.assertEqual(p.read_bytes(),before)

    def test_failed_exports_do_not_have_pass_footer_and_count_overflow_is_guarded(self):
        for mode in ("nonfinite-energy","nonfinite-coeff","overflow"):
            with self.subTest(mode=mode):
                case=self.work/(self._testMethodName+mode);case.mkdir()
                result=self.run_command([str(self.work/"pairs"),mode],case,False)
                self.assertNotEqual(result.returncode,0)
                self.assertIn("integer limit" if mode=="overflow" else "nonfinite",result.stderr)
                for path in case.glob("*.csv"):
                    self.assertNotIn("result_status=PASS",path.read_text())

    def test_mpi_parallel_batches_match_serial_in_source_order(self):
        if not self.has_mpi:
            self.skipTest("mpifort and mpiexec required")
        serial=self.execute("degenerate")
        parallel=self.execute("degenerate",mpi=True)
        for name in ("pairs-scalar.csv","pairs-spinor.csv"):
            self.assertEqual((serial/name).read_bytes(),(parallel/name).read_bytes())

    def test_parser_pair_export_is_explicit_and_exclusive(self):
        valid=["-kubo","2","-kubo_pairs","pairs.csv"]
        self.run_command([str(self.work/"parse"),*valid],self.work)
        invalid=[(["-kubo_pairs","pairs.csv"],"requires Kubo-only mode 2"),
                 (["-kubo","1","-kubo_pairs","pairs.csv"],"requires Kubo-only mode 2"),
                 ([*valid,"-kubo_csv","band.csv"],"exclusive"),
                 ([*valid,"-kubo_bundle","1","-kubo_csv","bundle.csv"],"exclusive"),
                 ([*valid,"-ii","1"],"omit -ii/-if/-is/-nn"),
                 ([*valid,"-if","3"],"omit -ii/-if/-is/-nn"),
                 ([*valid,"-is","1"],"omit -ii/-if/-is/-nn"),
                 ([*valid,"-nn","1"],"omit -ii/-if/-is/-nn"),
                 ([*valid,"-z2","1"],"requires Kubo-only mode 2"),
                 ([*valid,"-f","pairs.csv"],"cannot overwrite WAVECAR")]
        for args,message in invalid:
            with self.subTest(args=args):
                result=self.run_command([str(self.work/"parse"),*args],self.work,False)
                self.assertNotEqual(result.returncode,0)
                self.assertIn(message,result.stderr)


if __name__ == "__main__":
    unittest.main()
