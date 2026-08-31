import re
import subprocess
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_VERSION = "1.1.0"
ARCHIVED_V1_DOI = "10.5281/zenodo.1402593"


class VersionMetadataTests(unittest.TestCase):
    def test_version_file_is_authoritative(self):
        self.assertEqual(
            (ROOT / "VERSION").read_text(encoding="utf-8").strip(),
            EXPECTED_VERSION,
        )

    def test_current_citation_has_version_and_no_v1_doi(self):
        citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
        self.assertRegex(citation, rf"(?m)^version:\s*{re.escape(EXPECTED_VERSION)}\s*$")
        self.assertRegex(citation, r"(?m)^date-released:\s*2026-08-31\s*$")
        self.assertNotIn(ARCHIVED_V1_DOI, citation)
        self.assertNotRegex(citation, r"(?m)^doi:\s*")

    def test_changelog_records_scientific_release_boundaries(self):
        changelog = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
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
        self.assertIn("python -m unittest discover -s tests -v", workflow)
        self.assertIn("vaspberry_gfortran_serial.f", workflow)
        self.assertEqual(workflow.count("-fallow-argument-mismatch"), 2)
        self.assertIn('grep -F "Ver 1.1.0"', workflow)
        self.assertIn("mpifort -cpp -DMPI_USE", workflow)
        self.assertIn("vaspberry.f", workflow)

    def test_fortran_sources_report_current_version(self):
        for name in ("vaspberry.f", "vaspberry_gfortran_serial.f"):
            with self.subTest(source=name):
                source = (ROOT / name).read_text(encoding="utf-8", errors="replace")
                self.assertIn("Version 1.1.0", source)
                self.assertIn("Ver 1.1.0", source)

    def test_python_tools_report_current_version(self):
        for relative in (
            "tools/vaspberry_transport.py",
            "tools/wavecar_fukui.py",
            "tools/wavecar_z2.py",
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
