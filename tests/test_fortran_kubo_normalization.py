"""Compile the production Kubo term and WAVECAR-reading loop, then test physics.

No VASP source or material fixture is used. The QWZ oracle is analytic;
the separate plane-wave fixture exercises the actual Fortran loop and CSV.
"""
import csv
import ctypes
import io
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


class Complex16(ctypes.Structure):
    _fields_ = [("real", ctypes.c_double), ("imag", ctypes.c_double)]


@unittest.skipUnless(shutil.which("gfortran"), "gfortran required for compiled oracle")
class CompiledKuboTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix="vaspberry-kubo-")
        cls.work = Path(cls.temp.name)
        source = (ROOT / "vaspberry.f").read_text()
        kernel = routine(source, "kubo_interband_term", "function")
        (cls.work / "kernel.f").write_text(kernel)
        cls.run_command([
            "gfortran", "-O2", "-ffixed-line-length-none", "-fPIC", "-shared",
            "kernel.f", "-o", "libkubo.so",
        ])
        cls.library = ctypes.CDLL(str(cls.work / "libkubo.so"))
        cls.term = cls.library.kubo_interband_term_
        cls.term.restype = ctypes.c_double
        cls.term.argtypes = [ctypes.POINTER(Complex16), ctypes.POINTER(Complex16),
                             ctypes.POINTER(ctypes.c_double)]
        helpers = kernel + "\n".join(routine(source, name) for name in (
            "kubo_berry_curvature", "mpi_job_distribution_chain",
            "write_kubo_metadata", "write_kubo_band_csv",
        ))
        (cls.work / "production.f").write_text(helpers)
        cls.run_command([
            "gfortran", "-cpp", "-O2", "-fcheck=bounds", "-ffixed-line-length-none",
            "production.f", str(ROOT / "tests/fortran/test_kubo_loop.f90"),
            "-o", "test-kubo-loop",
        ])
        (cls.work / "parse.f").write_text(routine(source, "parse"))
        cls.run_command([
            "gfortran", "-cpp", "-O2", "-ffixed-line-length-none", "parse.f",
            str(ROOT / "tests/fortran/test_kubo_parser.f90"), "-o", "test-parser",
        ])

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    @classmethod
    def run_command(cls, command):
        result = subprocess.run(command, cwd=cls.work, capture_output=True, text=True)
        if result.returncode:
            raise AssertionError(f"{command}\n{result.stdout}\n{result.stderr}")
        return result.stdout

    @classmethod
    def omega(cls, px, py, gap):
        x, y = Complex16(px.real, px.imag), Complex16(py.real, py.imag)
        e = ctypes.c_double(gap)
        return cls.term(ctypes.byref(x), ctypes.byref(y), ctypes.byref(e))

    def test_qwz_pointwise_sign_factor_phase_and_energy_scale(self):
        sigma = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]])
        for kx, ky, mass in [(0.37, -1.21, -1), (1.8, .72, 1), (-.44, 2.4, 3)]:
            d = np.array([np.sin(kx), np.sin(ky), mass+np.cos(kx)+np.cos(ky)])
            dx = np.array([np.cos(kx), 0, -np.sin(kx)])
            dy = np.array([0, np.cos(ky), -np.sin(ky)])
            energies, u = np.linalg.eigh(np.einsum("a,aij->ij", d, sigma))
            vx = u.conj().T @ np.einsum("a,aij->ij", dx, sigma) @ u
            vy = u.conj().T @ np.einsum("a,aij->ij", dy, sigma) @ u
            expected = np.dot(d, np.cross(dx, dy))/(2*np.linalg.norm(d)**3)
            gap = energies[0]-energies[1]
            actual = self.omega(vx[0, 1], vy[0, 1], gap)
            self.assertAlmostEqual(actual, expected, delta=2e-13)
            self.assertAlmostEqual(self.omega(vx[1, 0], vy[1, 0], -gap), -expected, delta=2e-13)
            self.assertAlmostEqual(self.omega(vy[0, 1], vx[0, 1], gap), -expected, delta=2e-13)
            phase = np.exp(1.237j)
            self.assertAlmostEqual(self.omega(phase*vx[0, 1], phase*vy[0, 1], gap), expected, delta=2e-13)
            self.assertAlmostEqual(self.omega(7*vx[0, 1], 7*vy[0, 1], 7*gap), expected, delta=2e-13)

    def test_qwz_integral_is_one_not_two_and_trivial_phase_zero(self):
        sigma = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]])
        n = 64
        axis = (np.arange(n)+.5)*2*np.pi/n-np.pi
        kx, ky = np.meshgrid(axis, axis, indexing="ij")
        for mass, expected in [(-1, 1), (1, -1), (3, 0)]:
            d = np.stack((np.sin(kx), np.sin(ky), mass+np.cos(kx)+np.cos(ky)), axis=-1)
            dx = np.stack((np.cos(kx), np.zeros_like(kx), -np.sin(kx)), axis=-1)
            dy = np.stack((np.zeros_like(ky), np.cos(ky), -np.sin(ky)), axis=-1)
            energies, u = np.linalg.eigh(np.einsum("...a,aij->...ij", d, sigma))
            bra = u.conj().swapaxes(-1, -2)
            vx = bra @ np.einsum("...a,aij->...ij", dx, sigma) @ u
            vy = bra @ np.einsum("...a,aij->...ij", dy, sigma) @ u
            values = [self.omega(x, y, g) for x, y, g in zip(
                vx[..., 0, 1].ravel(), vy[..., 0, 1].ravel(), (energies[..., 0]-energies[..., 1]).ravel())]
            chern = sum(values)*(2*np.pi/n)**2/(2*np.pi)
            self.assertAlmostEqual(chern, expected, delta=1e-9)

    def test_actual_loop_scalar_spinor_and_normalized_csv(self):
        stdout = self.run_command([str(self.work / "test-kubo-loop")])
        self.assertIn("KUBO_LOOP_PASS", stdout)
        text = (self.work / "fixture.csv").read_text()
        self.assertIn("normalization=STANDARD_MINUS_TWO_IM", text)
        self.assertIn("WAVECAR_BARE_MOMENTUM_NO_PAW_NONLOCAL_VELOCITY", text)
        rows = list(csv.DictReader(io.StringIO("\n".join(
            line for line in text.splitlines() if not line.startswith("#")))))
        self.assertEqual([(int(r["spin"]), int(r["k_index"]), int(r["band"])) for r in rows],
                         [(1, 1, 1), (1, 1, 2)])
        self.assertEqual([float(r["energy_eV"]) for r in rows], [0, 2])
        self.assertEqual([float(r["min_gap_eV"]) for r in rows], [2, 2])
        self.assertAlmostEqual(float(rows[0]["omega_z_A2"]), -float(rows[1]["omega_z_A2"]), delta=1e-12)
        # An explicit new output path must never clobber an earlier run (or an
        # aliased input). The failed second run must preserve the first CSV.
        repeated = subprocess.run([str(self.work / "test-kubo-loop")],
                                  cwd=self.work, capture_output=True, text=True)
        self.assertNotEqual(repeated.returncode, 0)
        self.assertIn("cannot open new Kubo CSV", repeated.stderr)
        self.assertEqual((self.work / "fixture.csv").read_text(), text)

    def test_compiled_parser_preserves_output_basename_and_csv_contract(self):
        for arguments, expected in [(["-o", "sample"], "BERRYCURV.sample"),
                                    (["-kubo", "1"], "BERRYCURV_KUBO"),
                                    (["-o", "sample", "-z2", "1"], "sample"),
                                    (["-o", "sample", "-vel", "1"], "VEL_EXPT.sample")]:
            actual = self.run_command([str(self.work / "test-parser"), *arguments])
            self.assertEqual(actual.strip(), expected)
        invalid = subprocess.run([str(self.work / "test-parser"), "-kubo_csv", "out.csv"],
                                 cwd=self.work, capture_output=True, text=True)
        self.assertNotEqual(invalid.returncode, 0)
        self.assertIn("requires Kubo-only mode", invalid.stderr)


if __name__ == "__main__":
    unittest.main()
