"""Check both actual Bi runners publish PNG/PDF from stored native outputs.

The native execution boundary is replaced with archived public results; the
runner's validation, plotting and result-manifest paths execute unchanged.
"""
import contextlib
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import shutil
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "examples/features/fukui-chern"))


class BiFigureOutputs(unittest.TestCase):
    def check_runner(self, feature):
        directory = ROOT / "examples/features" / feature
        spec = importlib.util.spec_from_file_location("bi_figure_" + feature.replace("-", "_"), directory / "run.py")
        runner = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(runner)
        reference = directory / "reference"
        saved = json.loads((reference / "result.json").read_text())
        lines = (ROOT / "examples/Bi_Z2/archive-2016-run/EIGENVAL").read_text().splitlines()
        points = np.array([[float(x) for x in lines[7 + ik * 20].split()[:3]] for ik in range(144)])
        poscar = (ROOT / "examples/Bi_Z2/inputs/POSCAR").read_text().splitlines()
        lattice = float(poscar[1]) * np.array([[float(x) for x in line.split()] for line in poscar[2:5]])
        reciprocal = 2 * np.pi * np.linalg.inv(lattice).T
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp)
            for source in reference.iterdir():
                if source.is_file() and source.suffix in {".dat", ".csv", ".log"}:
                    shutil.copyfile(source, output / source.name)
            shutil.copyfile(reference / "command.json", output / "command.json")
            gaps = {key: saved["summary"][key] for key in ("minimum_direct_gap_eV", "sampled_global_gap_eV")}
            args = SimpleNamespace(wavecar=Path("unused-native-execution-is-replaced"))
            with mock.patch.object(runner, "arguments", return_value=args), \
                 mock.patch.object(runner, "execute", return_value=(output, {}, None, gaps)), \
                 mock.patch.object(runner, "wavecar_geometry", return_value=(points, reciprocal)), \
                 contextlib.redirect_stdout(io.StringIO()):
                runner.main()
            result = json.loads((output / "result.json").read_text())
            self.assertEqual(result["status"], "PASS")
            for name, signature in (("figure.png", b"\x89PNG"), ("figure.pdf", b"%PDF")):
                payload = (output / name).read_bytes()
                self.assertTrue(payload.startswith(signature), name)
                self.assertEqual(result["output_sha256"][name], hashlib.sha256(payload).hexdigest())

    def test_fukui_runner_publishes_png_and_pdf(self):
        self.check_runner("fukui-chern")

    def test_z2_runner_publishes_png_and_pdf(self):
        self.check_runner("z2")


if __name__ == "__main__":
    unittest.main()
