"""Regression checks for the public example layout and Bi Z2 fixture."""

from __future__ import annotations

import importlib.util
import os
import re
import stat
import struct
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"
BI = EXAMPLES / "Bi_Z2"
MOS2 = EXAMPLES / "1H-MoS2"
PLOT = BI / "scripts" / "plot_nfield.py"
RUNNER = BI / "scripts" / "run_z2.sh"
CURRENT_FIELD = BI / "reference-v1.2.0-12x12" / "Z2_FIELD.csv"
FIGURE_PDF = BI / "Z2_nfield_12x12.pdf"
FIGURE_PNG = BI / "Z2_nfield_12x12.png"


def _read_incar(path: Path) -> dict[str, str]:
    """Read the simple KEY=VALUE tags used by the reviewed templates."""
    tags: dict[str, str] = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.split("!", 1)[0].split("#", 1)[0]
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        tags[key.strip().upper()] = value.strip().upper()
    return tags


def _read_kpoints(path: Path) -> tuple[str, tuple[int, ...], tuple[int, ...]]:
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if len(lines) < 5:
        raise AssertionError(f"incomplete KPOINTS template: {path}")
    return lines[2].lower(), tuple(map(int, lines[3].split())), tuple(
        map(int, lines[4].split())
    )


def _load_plot_module():
    spec = importlib.util.spec_from_file_location("example_plot_nfield", PLOT)
    if spec is None or spec.loader is None:
        raise AssertionError(f"could not import {PLOT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _read_band_10_11_gaps(path: Path) -> tuple[int, float, float]:
    valence: list[float] = []
    conduction: list[float] = []
    direct: list[float] = []
    band_10: float | None = None
    for line in path.read_text(encoding="utf-8").splitlines():
        columns = line.split()
        if len(columns) != 3:
            continue
        try:
            band = int(columns[0])
            energy = float(columns[1])
            occupation = float(columns[2])
        except ValueError:
            continue
        if not 0.0 <= occupation <= 1.0:
            continue
        if band == 10:
            band_10 = energy
            valence.append(energy)
        elif band == 11 and band_10 is not None:
            conduction.append(energy)
            direct.append(energy - band_10)
    if not direct or len(valence) != len(conduction):
        raise AssertionError(f"could not read band-10/band-11 gaps from {path}")
    return len(direct), min(direct), min(conduction) - max(valence)


class ExamplesLayoutTests(unittest.TestCase):
    def test_examples_index_and_current_directories_exist(self) -> None:
        self.assertTrue((EXAMPLES / "README.md").is_file())
        self.assertTrue(BI.is_dir())
        self.assertTrue(MOS2.is_dir())

    def test_old_top_level_example_paths_are_absent(self) -> None:
        for old_path in ("Bi_Z2", "1H-MoS2"):
            with self.subTest(path=old_path):
                self.assertFalse((ROOT / old_path).exists())

        tracked = subprocess.check_output(
            ["git", "ls-files", "-z"], cwd=ROOT
        ).decode("utf-8").split("\0")
        offenders = sorted(
            path
            for path in tracked
            if path.startswith(("Bi_Z2/", "1H-MoS2/"))
        )
        self.assertEqual(offenders, [])

    def test_bi_scripts_have_stable_paths_and_executable_modes(self) -> None:
        for path in (RUNNER, PLOT):
            with self.subTest(path=path.relative_to(ROOT)):
                self.assertTrue(path.is_file())
                self.assertTrue(
                    path.stat().st_mode & stat.S_IXUSR,
                    f"{path.relative_to(ROOT)} must be executable",
                )
                self.assertTrue(os.access(path, os.X_OK))

        readme = (BI / "README.md").read_text(encoding="utf-8")
        self.assertIn("examples/Bi_Z2/scripts/run_z2.sh", readme)
        self.assertIn("examples/Bi_Z2/scripts/plot_nfield.py", readme)

    def test_bi_runner_has_valid_bash_syntax(self) -> None:
        subprocess.run(["bash", "-n", str(RUNNER)], check=True, cwd=ROOT)
        self.assertNotIn("--legacy", RUNNER.read_text(encoding="utf-8"))

    def test_plot_script_has_valid_python_syntax(self) -> None:
        source = PLOT.read_text(encoding="utf-8")
        compile(source, str(PLOT), "exec")

        help_result = subprocess.run(
            [sys.executable, str(PLOT), "-h"],
            check=True,
            cwd=ROOT,
            capture_output=True,
            text=True,
        )
        self.assertIn("--output", help_result.stdout)
        self.assertNotIn("legacy", help_result.stdout.lower())

    def test_reviewed_inputs_use_the_same_full_gamma_mesh(self) -> None:
        for stage in ("01_scf", "02_z2_nscf"):
            with self.subTest(stage=stage):
                mode, mesh, shift = _read_kpoints(
                    BI / "inputs" / stage / "KPOINTS"
                )
                self.assertEqual(mode, "gamma")
                self.assertEqual(mesh, (12, 12, 1))
                self.assertEqual(shift, (0, 0, 0))

    def test_reviewed_scf_and_z2_nscf_tags(self) -> None:
        scf = _read_incar(BI / "inputs" / "01_scf" / "INCAR")
        nscf = _read_incar(BI / "inputs" / "02_z2_nscf" / "INCAR")

        for name, tags in (("SCF", scf), ("Z2 NSCF", nscf)):
            with self.subTest(stage=name):
                self.assertEqual(tags["ISPIN"], "1")
                self.assertEqual(tags["LSORBIT"], ".TRUE.")
                self.assertEqual(tags["MAGMOM"], "6*0.0")
                self.assertEqual(tags["ISYM"], "-1")
                self.assertEqual(tags["LREAL"], ".FALSE.")
                self.assertEqual(tags["LASPH"], ".TRUE.")

        self.assertEqual(scf["ICHARG"], "2")
        self.assertEqual(scf["LORBIT"], "11")
        self.assertEqual(scf["LWAVE"], ".FALSE.")
        self.assertEqual(scf["LCHARG"], ".TRUE.")
        self.assertEqual(nscf["ICHARG"], "11")
        self.assertEqual(nscf["LWAVE"], ".TRUE.")
        self.assertEqual(nscf["LCHARG"], ".FALSE.")
        self.assertEqual(nscf["NBANDS"], "18")

    def test_archived_outcar_is_not_assigned_unverified_incar_provenance(self) -> None:
        archive = BI / "archive-2016-run"
        self.assertTrue((archive / "OUTCAR").is_file())
        self.assertTrue((archive / "INCAR.unverified").is_file())
        self.assertFalse((archive / "INCAR").exists())

        archive_note = (archive / "README.md").read_text(encoding="utf-8")
        self.assertIn("INCAR.unverified", archive_note)
        self.assertRegex(
            archive_note.lower(), r"not\s+asserted[\s\S]{0,80}provenance"
        )

        outcar = (archive / "OUTCAR").read_text(
            encoding="utf-8", errors="replace"
        )
        unverified = _read_incar(archive / "INCAR.unverified")
        self.assertRegex(outcar, r"(?m)^\s*ICHARG\s*=\s*11\b")
        self.assertEqual(unverified["ICHARG"], "2")

    def test_archived_eigenval_sampled_gaps_match_documentation(self) -> None:
        count, direct_gap, global_gap = _read_band_10_11_gaps(
            BI / "archive-2016-run" / "EIGENVAL"
        )
        self.assertEqual(count, 144)
        self.assertAlmostEqual(direct_gap, 0.592449, places=6)
        self.assertAlmostEqual(global_gap, 0.510045, places=6)

    def test_plot_parser_validates_current_result_contract(self) -> None:
        plot = _load_plot_module()
        field, metadata = plot.read_field(CURRENT_FIELD)

        self.assertEqual(len(field), 144)
        self.assertEqual(metadata["nkx"], "12")
        self.assertEqual(metadata["nky"], "12")
        self.assertEqual(metadata["schema_version"], "2")
        self.assertEqual(metadata["result_status"], "PASS")
        self.assertEqual(metadata["reportable_invariant"], "1")
        self.assertEqual(plot.reported_half_sums(metadata), (-3, 3))
        self.assertEqual(plot.half_sums(field), (-3, 3))
        self.assertEqual(tuple(map(plot.parity, (-3, 3))), (1, 1))
        self.assertEqual(
            plot.validate_result(field, metadata, CURRENT_FIELD), (12, 12, 1)
        )

    def test_plot_rejects_noncurrent_or_inconsistent_results(self) -> None:
        plot = _load_plot_module()
        field, original_metadata = plot.read_field(CURRENT_FIELD)
        cases = (
            ("old schema", "schema_version", "1", "version 2 or newer"),
            ("invalid status", "result_status", "INVALID", "must be PASS"),
            ("nonreportable", "reportable_invariant", "0", "not reportable"),
            ("wrong kind", "result_kind", "OTHER", "unsupported Z2 result kind"),
            ("small mesh", "nkx", "2", "greater than or equal to 4"),
            ("wrong sum", "half_top_nfield_sum", "-2", "do not reproduce"),
        )
        for label, key, value, message in cases:
            with self.subTest(case=label):
                changed_metadata = dict(original_metadata)
                changed_metadata[key] = value
                with self.assertRaisesRegex(ValueError, message):
                    plot.validate_result(field, changed_metadata, CURRENT_FIELD)

        even_metadata = dict(original_metadata)
        even_metadata.update(
            z2_invariant="0",
            half_top_z2_parity="0",
            half_bottom_z2_parity="0",
        )
        with self.assertRaisesRegex(ValueError, "do not reproduce Z2 parities"):
            plot.validate_result(field, even_metadata, CURRENT_FIELD)

    def test_plot_color_scale_preserves_larger_integer_fields(self) -> None:
        plot = _load_plot_module()
        field = {
            (-0.25, -0.25): -2,
            (0.25, -0.25): -1,
            (-0.25, 0.25): 1,
            (0.25, 0.25): 2,
        }
        cmap, norm, legend = plot.integer_style(field)

        mapped = [int(norm(value)) for value in (-2, -1, 0, 1, 2)]
        self.assertEqual(cmap.N, 5)
        self.assertEqual(len(set(mapped)), 5)
        self.assertGreaterEqual(min(mapped), 0)
        self.assertLess(max(mapped), cmap.N)
        self.assertEqual(len(legend), 5)

    def test_plot_cli_writes_single_panel_pdf_and_png(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            output = Path(temporary_directory) / "nfield.pdf"
            subprocess.run(
                [
                    sys.executable,
                    str(PLOT),
                    str(CURRENT_FIELD),
                    "--output",
                    str(output),
                ],
                check=True,
                cwd=ROOT,
                capture_output=True,
                text=True,
            )
            png_output = output.with_suffix(".png")
            self.assertTrue(output.is_file())
            self.assertTrue(png_output.is_file())

            png_header = png_output.read_bytes()[:24]
            self.assertEqual(png_header[:8], b"\x89PNG\r\n\x1a\n")
            width, height = struct.unpack(">II", png_header[16:24])
            self.assertGreater(width, 1000)
            self.assertGreater(height, 1000)
            self.assertGreater(width / height, 0.9)
            self.assertLess(width / height, 1.2)

    def test_checked_in_figure_is_a_single_page_pdf(self) -> None:
        payload = FIGURE_PDF.read_bytes()
        self.assertTrue(payload.startswith(b"%PDF-"))
        self.assertTrue(payload.rstrip().endswith(b"%%EOF"))
        self.assertEqual(len(re.findall(rb"/Type\s*/Page\b", payload)), 1)

        png_header = FIGURE_PNG.read_bytes()[:24]
        self.assertEqual(png_header[:8], b"\x89PNG\r\n\x1a\n")
        width, height = struct.unpack(">II", png_header[16:24])
        self.assertGreater(width / height, 0.9)
        self.assertLess(width / height, 1.2)

    def test_new_example_markdown_has_no_broken_local_links(self) -> None:
        markdown_files = (
            EXAMPLES / "README.md",
            BI / "README.md",
            BI / "inputs" / "README.md",
            BI / "archive-2016-run" / "README.md",
        )
        link_pattern = re.compile(r"(?<!!)\[[^]]+\]\(([^)]+)\)")
        broken: list[str] = []
        for markdown in markdown_files:
            text = markdown.read_text(encoding="utf-8")
            for target in link_pattern.findall(text):
                if target.startswith(("#", "http://", "https://", "mailto:")):
                    continue
                local_target = target.split("#", 1)[0]
                resolved = (markdown.parent / local_target).resolve()
                if local_target and not resolved.exists():
                    broken.append(f"{markdown.relative_to(ROOT)} -> {target}")
        self.assertEqual(broken, [])


if __name__ == "__main__":
    unittest.main()
