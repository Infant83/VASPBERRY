import re
import unittest
from pathlib import Path

import numpy as np


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


def fortran_subroutine(source: str, name: str) -> str:
    match = re.search(
        rf"(?ims)^\s*subroutine\s+{re.escape(name)}\b"
        rf".*?^\s*end\s+subroutine\s+{re.escape(name)}\s*$",
        source,
    )
    if match is None:
        raise AssertionError(f"could not locate Fortran subroutine {name}")
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

    def test_z2_kramers_partner_branches_are_reachable(self):
        for source_path in (SOURCE_PATH, GFORTRAN_SOURCE_PATH):
            source = source_path.read_text(encoding="utf-8", errors="strict")
            compact = re.sub(r"\s+", "", source).lower()
            with self.subTest(source=source_path.name):
                self.assertIn("if(itrim(iilp).ge.1.and.mod(nni,2).eq.0)then", compact)
                self.assertIn("source_band=nni-1", compact)
                self.assertIn("itheta=1", compact)
                self.assertIn("callz2_apply_theta(", compact)
                theta = re.sub(
                    r"\s+", "", fortran_subroutine(source, "z2_apply_theta")
                ).lower()
                self.assertIn("coutup=-conjg(cindn)", theta)
                self.assertIn("coutdn=conjg(cinup)", theta)

    def test_z2_accumulators_are_reset_inside_each_spin_loop(self):
        reset_sequence = re.compile(
            r"(?ims)^\s*do\s+isp\s*=\s*1\s*,\s*ispin\s*!\s*ispin start\s*$"
            r".*?^\s*recivec\s*=\s*0d0\s*$"
            r"\n^\s*recilat\s*=\s*0d0\s*$"
            r".*?^\s*rnnfield\s*=\s*0d0\s*$"
            r"\n^\s*rnnfield_bottom\s*=\s*0d0\s*$"
            r".*?^\s*rnfield\s*=\s*0d0\s*$"
        )
        for source_path in (SOURCE_PATH, GFORTRAN_SOURCE_PATH):
            source = source_path.read_text(encoding="utf-8", errors="strict")
            with self.subTest(source=source_path.name):
                self.assertRegex(source, reset_sequence)

    def test_z2_guards_reject_unsupported_inputs(self):
        required = (
            "Z2 option must be 0 or 1",
            "legacy Z2 needs a full k mesh",
            "legacy Z2 needs ISPIN=1 spinors",
            "legacy Z2 grid must be at least 4x4",
            "legacy Z2 grid must be even",
            "legacy Z2 needs the full 2D mesh",
            "legacy Z2 needs one of each TRIM",
            "legacy Z2 needs 2 <= NE < NBANDS",
            "legacy Z2 needs bands 1 through NE",
            "legacy Z2 needs an even band rank",
        )
        for source_path in (SOURCE_PATH, GFORTRAN_SOURCE_PATH):
            source = source_path.read_text(encoding="utf-8", errors="strict")
            for fragment in required:
                with self.subTest(source=source_path.name, fragment=fragment):
                    self.assertIn(fragment, source)

    def test_z2_index_and_lapack_failures_are_guarded(self):
        required = (
            "incomplete five-point k loop",
            "selected k outside range",
            "ZGETRF failed, INFO=",
        )
        for source_path in (SOURCE_PATH, GFORTRAN_SOURCE_PATH):
            source = source_path.read_text(encoding="utf-8", errors="strict")
            compact = re.sub(r"\s+", "", source).lower()
            for fragment in required:
                with self.subTest(source=source_path.name, fragment=fragment):
                    self.assertIn(fragment, source)
            self.assertNotIn("kpoint.gt.nk", compact)
            self.assertRegex(
                source,
                r"(?is)call\s+zgetrf\s*\(.*?\)\s*"
                r"if\s*\(\s*info\s*\.ne\.\s*0\s*\)\s*then",
            )

    def test_z2_parity_uses_nonnegative_modulo(self):
        for source_path in (SOURCE_PATH, GFORTRAN_SOURCE_PATH):
            source = source_path.read_text(encoding="utf-8", errors="strict")
            candidate_lines = [
                line
                for line in source.splitlines()
                if "modulo(nint(rnnfield" in line
                or "modulo(nint(rvari" in line
            ]
            with self.subTest(source=source_path.name):
                self.assertEqual(len(candidate_lines), 4)

    def test_legacy_z2_output_is_explicitly_a_candidate(self):
        required = (
            "Legacy Z2 candidate",
            "Do not report it as an invariant.",
            "tools/wavecar_z2.py",
            "Z2 only when diagnostics PASS.",
        )
        for source_path in (SOURCE_PATH, GFORTRAN_SOURCE_PATH):
            source = source_path.read_text(encoding="utf-8", errors="strict")
            for fragment in required:
                with self.subTest(source=source_path.name, fragment=fragment):
                    self.assertIn(fragment, source)

    def test_new_z2_failures_return_nonzero_status(self):
        for source_path in (SOURCE_PATH, GFORTRAN_SOURCE_PATH):
            source = source_path.read_text(encoding="utf-8", errors="strict")
            lines = source.splitlines()
            error_indices = [
                index
                for index, line in enumerate(lines)
                if "*** error - Z2" in line
                or "*** error - legacy Z2" in line
                or "*** error - incomplete five-point" in line
            ]
            with self.subTest(source=source_path.name):
                self.assertTrue(error_indices)
                for index in error_indices:
                    following_lines = (
                        line.strip().lower()
                        for line in lines[index + 1 :]
                        if line.strip()
                        and not line.lstrip().startswith("&")
                    )
                    self.assertEqual(next(following_lines), "call vaspberry_fail")

                self.assertRegex(
                    source,
                    r"(?im)^\s*subroutine\s+vaspberry_fail\s*$",
                )
                self.assertRegex(source, r"(?im)^\s*stop\s+1\s*$")

        mpi_source = SOURCE_PATH.read_text(encoding="utf-8", errors="strict")
        self.assertIn("call MPI_ABORT(MPI_COMM_WORLD,1,ierr)", mpi_source)

    def test_z2_uses_exact_folded_reciprocal_mapping_and_flux_schema(self):
        required = (
            "subroutine z2_map_gvector",
            "delta=-wks(j)-wkt(j)",
            "gt(j)=-gs(j)+ishift",
            "source_band=nni-1",
            "subroutine z2_apply_theta",
            "coutup=-conjg(cindn)",
            "coutdn=conjg(cinup)",
            "rflux=-wilson_phase",
            "pi=4d0*datan(1d0)",
            "WAVECAR_PSEUDO_NO_PAW_AUGMENTATION",
            "physical_tr_rule=",
            "nfield_note=gauge_and_log_branch_dependent",
            "max_flux_tr_odd_residual_rad",
        )
        for source_path in (SOURCE_PATH, GFORTRAN_SOURCE_PATH):
            source = source_path.read_text(encoding="utf-8", errors="strict")
            compact = re.sub(r"\s+", "", source).lower()
            with self.subTest(source=source_path.name):
                for fragment in required:
                    self.assertIn(
                        re.sub(r"\s+", "", fragment).lower(),
                        compact,
                    )

    def test_z2_field_writer_checks_partner_map_without_gating_nfield_pairs(self):
        for source_path in (SOURCE_PATH, GFORTRAN_SOURCE_PATH):
            source = source_path.read_text(encoding="utf-8", errors="strict")
            writer = re.search(
                r"(?ims)^\s*subroutine\s+write_z2_field_csv\b"
                r".*?^\s*end\s+subroutine\s+write_z2_field_csv\s*$",
                source,
            )
            self.assertIsNotNone(writer)
            body = re.sub(r"\s+", "", writer.group(0)).lower()
            self.assertIn("recilat(m,i)+recilat(m,j)", body)
            self.assertIn("partner(partner(i)).ne.i", body)
            self.assertIn("maxfluxres.gt.fluxtol", body)
            self.assertIn("maxnres.gt.1d-10", body)
            self.assertIn("minsv.lt.minsvtol", body)
            self.assertNotIn("maxnpairres.gt.", body)
            self.assertIn("i=1,nnk", body)


    def test_complex8_norm_promotes_components_before_squaring(self):
        value = np.complex64(0.1 + 0.2j)
        legacy = float(np.abs(value)) ** 2
        promoted = float(value.real) ** 2 + float(value.imag) ** 2
        mapped_complex16 = abs(complex(value)) ** 2
        self.assertGreater(abs(legacy - promoted), 1.0e-10)
        self.assertAlmostEqual(promoted, mapped_complex16, places=15)

        for source_path in (SOURCE_PATH, GFORTRAN_SOURCE_PATH):
            source = source_path.read_text(encoding="utf-8", errors="strict")
            helper = re.sub(
                r"\s+", "", fortran_subroutine(source, "z2_raw_norm")
            ).lower()
            with self.subTest(source=source_path.name):
                self.assertIn(
                    "norm=norm+dble(real(raw(i)))**2"
                    "+dble(aimag(raw(i)))**2",
                    helper,
                )
                self.assertNotIn("dble(abs(raw(i)))**2", source.lower())
                self.assertIn("call z2_raw_norm(raw,nraw,norm_in)", source)

    def test_determinant_magnitude_is_not_a_conditioning_proxy(self):
        hidden_singularity = np.diag([1.0, 1.0e-8]).astype(np.complex128)
        high_rank = 0.1 * np.eye(18, dtype=np.complex128)
        self.assertGreater(abs(np.linalg.det(hidden_singularity)), 1.0e-14)
        self.assertLess(
            np.linalg.svd(hidden_singularity, compute_uv=False)[-1],
            1.0e-6,
        )
        self.assertLess(abs(np.linalg.det(high_rank)), 1.0e-14)
        self.assertGreater(
            np.linalg.svd(high_rank, compute_uv=False)[-1],
            1.0e-6,
        )

    def test_z2_link_quality_uses_svd_and_robust_lu_phase(self):
        for source_path in (SOURCE_PATH, GFORTRAN_SOURCE_PATH):
            source = source_path.read_text(encoding="utf-8", errors="strict")
            helper = re.sub(
                r"\s+",
                "",
                fortran_subroutine(source, "z2_link_svd_phase"),
            ).lower()
            nfield = re.sub(
                r"\s+", "", fortran_subroutine(source, "get_nfield")
            ).lower()
            writer = re.sub(
                r"\s+",
                "",
                fortran_subroutine(source, "write_z2_field_csv"),
            ).lower()
            with self.subTest(source=source_path.name):
                self.assertIn("callzgesvd(", helper)
                self.assertIn("asvd=sijt", helper)
                self.assertIn("alu=sijt", helper)
                self.assertIn("callzgetrf(", helper)
                self.assertIn("if(ipiv(i).ne.i)phase=phase+pi", helper)
                self.assertIn("callz2_link_svd_phase(", nfield)
                self.assertNotIn("detmag.le.1d-14", nfield)
                self.assertIn("threshold_min_link_singular", writer)
                self.assertIn("minimum_link_singular_value", writer)
                self.assertIn("minsv.lt.minsvtol", writer)

    def test_z2_output_lifecycle_prevents_stale_pass_files(self):
        for source_path in (SOURCE_PATH, GFORTRAN_SOURCE_PATH):
            source = source_path.read_text(encoding="utf-8", errors="strict")
            main_end = re.search(r"(?im)^\s*end\s+program\b", source)
            self.assertIsNotNone(main_end)
            main = re.sub(r"\s+", "", source[: main_end.start()]).lower()
            begin = main.index("callbegin_z2_field_outputs(filename)")
            input_read = main.index("callinforead(")
            self.assertLess(begin, input_read)

            prepare = re.sub(
                r"\s+",
                "",
                fortran_subroutine(source, "begin_z2_field_outputs"),
            ).lower()
            writer = re.sub(
                r"\s+",
                "",
                fortran_subroutine(source, "write_z2_field_csv"),
            ).lower()
            with self.subTest(source=source_path.name):
                self.assertIn(
                    "trim(filename).eq.'z2_field.csv'", prepare
                )
                self.assertIn("call delete_z2_file('z2_field.csv')", prepare)
                self.assertIn("#result_status=incomplete", prepare)
                self.assertIn("file='z2_field.tmp'", writer)
                self.assertIn("callz2_atomic_replace(", writer)
                self.assertIn("#result_status=pass", writer)
                self.assertIn("#result_status=invalid", writer)
                self.assertIn("#reportable_invariant=0", writer)
                self.assertIn(
                    "#input_trs_independently_verified=0", writer
                )

    def test_mpi_z2_field_status_is_broadcast_before_common_gate(self):
        source = SOURCE_PATH.read_text(encoding="utf-8", errors="strict")
        main_end = re.search(r"(?im)^\s*end\s+program\b", source)
        self.assertIsNotNone(main_end)
        main = re.sub(r"\s+", "", source[: main_end.start()]).lower()
        writer = main.index("callwrite_z2_field_csv(")
        broadcast_token = (
            "callmpi_bcast(ifieldok,1,mpi_integer,0,"
            "mpi_comm_world,ierr)"
        )
        broadcast = main.index(broadcast_token)
        gate = main.index("if(ifieldok.ne.1)then")
        self.assertLess(writer, broadcast)
        self.assertLess(broadcast, gate)
        self.assertEqual(main.count(broadcast_token), 1)
        self.assertIn(
            "callmpi_reduce(rminsv,rminsv_tot,iz2*nk,mpi_real8,"
            "mpi_sum,0,mpi_comm_world,ierr)",
            main,
        )


    def test_z2_helper_implementations_do_not_drift_between_sources(self):
        mpi_source = SOURCE_PATH.read_text(encoding="utf-8", errors="strict")
        serial_source = GFORTRAN_SOURCE_PATH.read_text(
            encoding="utf-8", errors="strict"
        )
        for name in (
            "begin_z2_field_outputs",
            "delete_z2_file",
            "z2_atomic_replace",
            "write_z2_field_csv",
            "get_nfield",
            "z2_link_svd_phase",
            "z2_map_gvector",
            "z2_raw_norm",
            "z2_apply_theta",
            "get_z2_state",
        ):
            with self.subTest(subroutine=name):
                self.assertEqual(
                    fortran_subroutine(mpi_source, name),
                    fortran_subroutine(serial_source, name),
                )

    def test_fortran_version_banner_is_1_1_1(self):
        for source_path in (SOURCE_PATH, GFORTRAN_SOURCE_PATH):
            source = source_path.read_text(encoding="utf-8", errors="strict")
            with self.subTest(source=source_path.name):
                self.assertIn("PROGRAM VASPBERRY Version 1.1.1", source)
                self.assertIn("# VASPBERRY (Ver 1.1.1)", source)


if __name__ == "__main__":
    unittest.main()
