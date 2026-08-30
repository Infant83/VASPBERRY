import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_PATH = ROOT / "vaspberry.f"
GFORTRAN_SOURCE_PATH = ROOT / "vaspberry_gfortran_serial.f"


def source_text() -> str:
    return SOURCE_PATH.read_text(encoding="utf-8", errors="strict")


def parse_subroutine(source: str) -> str:
    match = re.search(
        r"(?ims)^\s*subroutine\s+parse\b.*?^\s*end\s+subroutine\s+parse\s*$",
        source,
    )
    if match is None:
        raise AssertionError("could not locate the Fortran parse subroutine")
    return match.group(0)


class FortranSourceRegressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.source = source_text()
        cls.parse = parse_subroutine(cls.source)

    def test_kubo_remains_opt_in_and_default_fukui_branch_is_present(self):
        self.assertRegex(self.parse, r"(?i)\bikubo\s*=\s*0\b")
        self.assertIn(
            "if(icd+ivel+iz+iwf+ikubo.eq.0)then",
            self.source.replace(" ", "").lower(),
        )

    def test_fortran_command_line_interface_is_unchanged(self):
        expected_options = {
            "-atlist",
            "-cd",
            "-f",
            "-fbz",
            "-fen",
            "-h",
            "-hf",
            "-ien",
            "-if",
            "-ii",
            "-im",
            "-is",
            "-ishift",
            "-ixt",
            "-k",
            "-klist",
            "-kp",
            "-kubo",
            "-kx",
            "-ky",
            "-ne",
            "-nediv",
            "-ng",
            "-nn",
            "-o",
            "-phi",
            "-s",
            "-sigma",
            "-skp",
            "-sw_file",
            "-t",
            "-theta",
            "-vel",
            "-wf",
            "-z2",
        }
        # Limit the extraction to option comparisons.  The parse routine also
        # contains output-name fragments such as "-E00", which are not CLI
        # switches.
        actual_options = set(
            re.findall(
                r'(?i)(?:else\s+)?if\s*\(\s*option\s*==\s*"(-[A-Za-z0-9_]+)"',
                self.parse,
            )
        )
        self.assertEqual(actual_options, expected_options)
        self.assertNotRegex(self.parse, r"(?i)transport|hall|ahc")

    def test_kubo_total_is_initialized_immediately_after_allocation(self):
        allocation_and_initialization = re.compile(
            r"(?im)^\s*if\s*\(\s*ikubo\.ge\.1\s*\)\s*"
            r"allocate\s*\(\s*berrycurv_kubo_tot\s*\(\s*nk\s*\)\s*\)\s*$"
            r"\n^\s*if\s*\(\s*ikubo\.ge\.1\s*\)\s*"
            r"berrycurv_kubo_tot\s*=\s*0d0\s*$"
        )
        self.assertRegex(self.source, allocation_and_initialization)

    def test_kubo_total_is_reset_for_each_spin_before_band_loop(self):
        spin_loop_match = re.search(
            r"(?ims)^\s*do\s+isp\s*=\s*1\s*,\s*ispin\s*!\s*ispin start\s*$"
            r"(?P<body>.*?)"
            r"^\s*enddo\s*!\s*ispin loop over\s*$",
            self.source,
        )
        self.assertIsNotNone(spin_loop_match)
        spin_loop = spin_loop_match.group("body")

        kubo_branch_match = re.search(
            r"(?ims)^\s*elseif\s*\(\s*ikubo\.ge\.1\s*\)\s*then\s*$"
            r"(?P<body>.*?)"
            r"^\s*endif\s*!\s*icd over\s*$",
            spin_loop,
        )
        self.assertIsNotNone(kubo_branch_match)
        kubo_branch = kubo_branch_match.group("body")

        reset_before_band_loop = re.compile(
            r"(?ims)^\s*berrycurv_kubo_tot\s*=\s*0d0\s*$"
            r".*?"
            r"^\s*do\s+ie\s*=\s*nini\s*,\s*nmax\s*$"
        )
        self.assertRegex(kubo_branch, reset_before_band_loop)

    def test_legacy_fukui_formula_and_chern_accumulation_are_unchanged(self):
        compact = re.sub(r"\s+", "", self.source).lower()
        self.assertIn(
            "berrycurv(ik)=-1.*aimag(log(detloop))/dskxky",
            compact,
        )
        self.assertIn(
            "chernnumber=chernnumber+berrycurv(ik)*dskxky/(2.*pi)",
            compact,
        )

    def test_legacy_result_headers_and_numeric_format_are_preserved(self):
        required_fragments = (
            '"# NKPOINT          : "',
            '"#  K-GRID          : "',
            '"# NBANDS           : "',
            '"# Chern Number =   "',
            '"# Berry Curvature (A^2) :',
            '"# (cart) kx        ky        kz(A^-1)',
            "'(3F11.6,A,F20.4,A,3F11.6)'",
            "'.UP.dat'",
            "'.DN.dat'",
            "'.EIG-'",
        )
        for fragment in required_fragments:
            with self.subTest(fragment=fragment):
                self.assertIn(fragment, self.source)

    def test_both_supported_help_paths_point_to_guarded_transport_tools(self):
        required = (
            "*VALLEY TRANSPORT POSTPROCESSING:",
            "python3 tools/vaspberry_transport.py -h",
            "python3 tools/wavecar_fukui.py -h",
            "A full uniform 2D BZ is required for transport.",
            "A CBM-only result is an incremental contribution.",
            "Total sigma_xy also needs the valence manifold.",
        )
        for source_path in (SOURCE_PATH, GFORTRAN_SOURCE_PATH):
            source = source_path.read_text(encoding="utf-8", errors="strict")
            for fragment in required:
                with self.subTest(source=source_path.name, fragment=fragment):
                    self.assertIn(fragment, source)


if __name__ == "__main__":
    unittest.main()
