import os
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_VERSION = "1.2.0"
ARCHIVED_V1_DOI = "10.5281/zenodo.1402593"


V110_HISTORY = """! version 1.1.0 safety fixes for the legacy Z2 candidate
!               : 2026. Aug. 31. H.-J. Kim"""

V111_HISTORY = """! version 1.1.1 correct TRIM reciprocal-G mapping and
!               separate physical flux from gauge-dependent n-field
!               : 2026. Aug. 31. H.-J. Kim"""

V120_HISTORY = """! version 1.2.0 formal Fukui-Hatsugai n-field Z2 interface
!               and GNU/OpenMPI portability checks
!               : 2026. Sep. 04."""


class VersionMetadataTests(unittest.TestCase):
    def test_version_file_is_authoritative(self):
        self.assertEqual(
            (ROOT / "VERSION").read_text(encoding="utf-8").strip(),
            EXPECTED_VERSION,
        )

    def test_current_citation_has_version_and_no_v1_doi(self):
        citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
        self.assertRegex(citation, rf"(?m)^version:\s*{re.escape(EXPECTED_VERSION)}\s*$")
        self.assertRegex(citation, r"(?m)^date-released:\s*2026-09-04\s*$")
        self.assertNotIn(ARCHIVED_V1_DOI, citation)
        self.assertNotRegex(citation, r"(?m)^doi:\s*")

    def test_changelog_records_scientific_release_boundaries(self):
        changelog = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
        self.assertIn("## [1.2.0] - 2026-09-04", changelog)
        self.assertIn("FUKUI_HATSUGAI_NFIELD_Z2", changelog)
        self.assertIn("## [1.1.1] - 2026-08-31", changelog)
        self.assertIn("reciprocal-G", changelog)
        self.assertIn("WAVECAR_PSEUDO_NO_PAW_AUGMENTATION", changelog)
        self.assertIn("## [1.1.0] - 2026-08-31", changelog)
        self.assertIn("legacy Z2 candidate", changelog)
        self.assertIn("unconverged mesh", changelog)
        self.assertIn(ARCHIVED_V1_DOI, changelog)
        self.assertIn("does not create", changelog)
        self.assertIn("GitHub Release", changelog)
        self.assertIn("transport_t0.csv", changelog)
        self.assertIn("python tools/wavecar_z2.py WAVECAR", changelog)
        self.assertIn("z2_wilson_wcc.csv", changelog)
        self.assertIn("`z2: null`", changelog)

    def test_ci_covers_python_serial_and_mpi_paths(self):
        workflow = (ROOT / ".github/workflows/ci.yml").read_text(encoding="utf-8")
        makefile = (ROOT / "Makefile").read_text(encoding="utf-8")
        self.assertIn("python -m unittest discover -s tests -v", workflow)
        self.assertIn("gfortran-11", workflow)
        self.assertIn("gfortran-13", workflow)
        self.assertIn("make gnu", workflow)
        self.assertIn("make check-gnu", workflow)
        self.assertIn("MPIEXEC_FLAGS=--oversubscribe", workflow)
        self.assertIn(
            "fortran-lang/setup-fortran@"
            "be037f0a45b1160f139d4ccd0b96f9e9bfd6a682 # v1",
            workflow,
        )
        self.assertIn('compiler: intel\n', workflow)
        self.assertIn('compiler: intel-classic\n', workflow)
        self.assertIn('version: "2025.0"', workflow)
        self.assertIn('version: "2021.10"', workflow)
        self.assertIn("IFX_MKL_FLAGS='-llapack -lblas'", workflow)
        self.assertIn("IFORT_MKL_FLAGS='-llapack -lblas'", workflow)
        self.assertIn("vaspberry_gfortran_serial.f", makefile)
        self.assertIn("-fallow-argument-mismatch", makefile)
        self.assertIn("-DMPI_USE", makefile)
        self.assertIn("tests/fortran/test_mpi_runtime.f90", makefile)
        self.assertIn("tests/check_fortran_help.sh", makefile)
        self.assertIn("-assume byterecl", makefile)
        self.assertIn("-qmkl=sequential", makefile)
        self.assertIn("ifx-mpi:", makefile)
        self.assertIn("ifort-mpi:", makefile)
        self.assertIn("check-build-dir:", makefile)
        self.assertIn("SAFE_BUILD_DIR", makefile)
        self.assertNotIn("-qmkl-ilp64", makefile)
        self.assertNotIn("-fdefault-integer-8", makefile)
        self.assertIn("mpifort -cpp -DMPI_USE", workflow)
        self.assertIn("vaspberry.f", workflow)
        self.assertIn("test_z2_helpers.f90", workflow)
        self.assertIn("objcopy --redefine-sym", workflow)

    def test_make_clean_rejects_a_build_directory_outside_the_repository(self):
        with tempfile.TemporaryDirectory(dir=ROOT.parent) as temp_dir:
            outside = Path(temp_dir) / "outside-build"
            outside.mkdir()
            sentinel = outside / "must-survive.txt"
            sentinel.write_text("keep\n", encoding="utf-8")
            completed = subprocess.run(
                ["make", "-s", "clean", f"BUILD_DIR={outside}"],
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(completed.returncode, 0)
            self.assertTrue(sentinel.is_file())
            self.assertIn(
                "BUILD_DIR must be build or a build-* directory",
                completed.stdout + completed.stderr,
            )

        protected = ROOT / "examples" / "README.md"
        completed = subprocess.run(
            ["make", "-s", "clean", "BUILD_DIR=examples"],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertTrue(protected.is_file())

    def test_bi_workflow_compares_serial_and_two_rank_mpi_results(self):
        workflow = (
            ROOT / ".github/workflows/bi-z2-validation.yml"
        ).read_text(encoding="utf-8")
        runner = (
            ROOT / "examples/Bi_Z2/scripts/run_z2.sh"
        ).read_text(encoding="utf-8")
        for required in (
            "make gnu",
            "results-z2-serial",
            "results-z2-mpi",
            "VASPBERRY_MPI_NPROCS=2",
            "tests/compare_z2_fields.py",
        ):
            self.assertIn(required, workflow)
        self.assertIn('run_command=("$launcher" -n "$mpi_nprocs"', runner)
        self.assertNotIn("eval ", runner)

    def test_bi_runner_rejects_unsafe_mpi_process_counts(self):
        runner = ROOT / "examples/Bi_Z2/scripts/run_z2.sh"
        for value in ("0", "145", "001", "abc", "2;echo unsafe"):
            with self.subTest(value=value):
                completed = subprocess.run(
                    ["bash", str(runner)],
                    cwd=ROOT,
                    env={**os.environ, "VASPBERRY_MPI_NPROCS": value},
                    check=False,
                    capture_output=True,
                    text=True,
                )
                self.assertEqual(completed.returncode, 2)
                self.assertIn(
                    "VASPBERRY_MPI_NPROCS must be an integer from 1 to 144",
                    completed.stderr,
                )

    def test_fortran_sources_report_current_version(self):
        for name in ("vaspberry.f", "vaspberry_gfortran_serial.f"):
            with self.subTest(source=name):
                source = (ROOT / name).read_text(encoding="utf-8", errors="replace")
                self.assertIn(f"Version {EXPECTED_VERSION}", source)
                self.assertIn(f"Ver {EXPECTED_VERSION}", source)


    def test_fortran_version_history_is_preserved_exactly(self):
        for name in ("vaspberry.f", "vaspberry_gfortran_serial.f"):
            source = (ROOT / name).read_text(
                encoding="utf-8", errors="strict"
            )
            with self.subTest(source=name):
                self.assertEqual(source.count(V110_HISTORY), 1)
                self.assertEqual(source.count(V111_HISTORY), 1)
                self.assertEqual(source.count(V120_HISTORY), 1)
                self.assertLess(
                    source.index(V110_HISTORY), source.index(V111_HISTORY)
                )
                self.assertLess(
                    source.index(V111_HISTORY), source.index(V120_HISTORY)
                )

    def test_python_tools_report_current_version(self):
        for relative in (
            "tools/vaspberry_transport.py",
            "tools/wavecar_fukui.py",
        ):
            with self.subTest(tool=relative):
                completed = subprocess.run(
                    [sys.executable, str(ROOT / relative), "--version"],
                    cwd=ROOT,
                    check=False,
                    capture_output=True,
                    text=True,
                )
                self.assertEqual(completed.returncode, 0, completed.stderr)
                self.assertRegex(
                    completed.stdout.strip(),
                    rf"(?<![0-9.]){re.escape(EXPECTED_VERSION)}(?![0-9.])",
                )


if __name__ == "__main__":
    unittest.main()
