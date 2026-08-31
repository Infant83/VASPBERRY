"""Time-reversal regression for the public 1H-MoS2 Berry map.

This checks the published occupied-subspace example, not a rerun of the
corrected Fortran Z2 path.  The file contains a 3x3 periodic replication of
the fundamental 12x12 map, so the test folds its plaquette centers modulo one.
"""

from __future__ import annotations

import math
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "1H-MoS2" / "BERRYCURV.dat"


def _mod1(value: float) -> float:
    folded = value - math.floor(value)
    return 0.0 if math.isclose(folded, 1.0, abs_tol=5.0e-6) else folded


def _key(q1: float, q2: float) -> tuple[float, float]:
    return (round(_mod1(q1), 6), round(_mod1(q2), 6))


class MoS2ExampleTimeReversalTests(unittest.TestCase):
    def test_occupied_subspace_curvature_is_tr_odd(self) -> None:
        lines = EXAMPLE.read_text(encoding="utf-8").splitlines()
        rows = []
        for line in lines:
            fields = line.split()
            if len(fields) < 7:
                continue
            try:
                values = [float(value) for value in fields[:7]]
            except ValueError:
                continue
            rows.append(values)

        self.assertEqual(len(rows), 1296)
        fundamental: dict[tuple[float, float], float] = {}
        for values in rows:
            fundamental.setdefault(_key(values[4], values[5]), values[3])
        self.assertEqual(len(fundamental), 144)

        residuals = []
        for (q1, q2), curvature in fundamental.items():
            partner = _key(-q1, -q2)
            self.assertIn(partner, fundamental)
            residuals.append(abs(curvature + fundamental[partner]))
        self.assertLessEqual(max(residuals), 5.0e-6)

        chern_lines = [
            line for line in lines if line.startswith("# Chern Number =")
        ]
        self.assertEqual(len(chern_lines), 1)
        self.assertAlmostEqual(float(chern_lines[0].split("=")[1]), 0.0)


if __name__ == "__main__":
    unittest.main()
