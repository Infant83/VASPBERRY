"""Time-reversal regression for the public 1H-MoS2 Berry map.

This checks the published occupied-subspace example, not a rerun of the
corrected Fortran Z2 path.  The file contains a 3x3 periodic replication of
the fundamental 12x12 map, so the test folds its plaquette centers modulo one.
"""

from __future__ import annotations

import math
import unittest
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples" / "1H-MoS2" / "BERRYCURV.dat"


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
        groups: dict[
            tuple[float, float], list[tuple[int, int, float]]
        ] = defaultdict(list)
        for values in rows:
            q1, q2 = values[4], values[5]
            folded_q1, folded_q2 = _mod1(q1), _mod1(q2)
            shift_q1 = int(round(q1 - folded_q1))
            shift_q2 = int(round(q2 - folded_q2))
            self.assertAlmostEqual(q1, folded_q1 + shift_q1, delta=1.0e-12)
            self.assertAlmostEqual(q2, folded_q2 + shift_q2, delta=1.0e-12)
            groups[_key(q1, q2)].append(
                (shift_q1, shift_q2, values[3])
            )
        self.assertEqual(len(groups), 144)

        fundamental: dict[tuple[float, float], float] = {}
        for key, copies in groups.items():
            with self.subTest(cell=key):
                self.assertEqual(len(copies), 9)
                shifts = {(sx, sy) for sx, sy, _ in copies}
                x_images = {sx for sx, _, _ in copies}
                y_images = {sy for _, sy, _ in copies}
                self.assertEqual(len(x_images), 3)
                self.assertEqual(len(y_images), 3)
                self.assertEqual(
                    shifts,
                    {
                        (sx, sy)
                        for sx in x_images
                        for sy in y_images
                    },
                )
                curvatures = {curvature for _, _, curvature in copies}
                self.assertEqual(len(curvatures), 1)
                fundamental[key] = copies[0][2]

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
