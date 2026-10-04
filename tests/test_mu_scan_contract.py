"""Public Hall entry points agree with the INI single-point/range contract."""
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))

import kubo_hall_workflow as pairs
import vaspberry_kubo as curvature
import waveder_hall as optical
import wannier_workflow as wannier


class ReachedCalculation(Exception):
    """Stop before independent source-loading or integration work."""


class MuScanContractTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.output = Path(temporary.name)/'new-output'

    def arguments(self, lower, upper, count):
        return SimpleNamespace(mu_min=lower, mu_max=upper, mu_num=count,
            mu_reference=0., regions=None, difference=[], temperatures=[0.],
            output_dir=self.output, curvature=Path('unused-curvature'),
            operators=Path('unused-operators'), workers=1, batch_size=1,
            time_limit=1., gap_threshold=1e-7, formats=['csv'],
            allow_partial_bands=False, band_resolved=False, mu_chunk=32)

    def invoke(self, route, args):
        if route == 'pairs':
            pairs.scan_options(args)
            raise ReachedCalculation
        if route == 'waveder':
            with patch.object(optical, 'resolve_cli_inputs', side_effect=ReachedCalculation):
                optical.command(args)
        elif route == 'curvature':
            # This route loads its saved curvature before checking the scan.
            with patch.object(curvature, 'read_curvature', return_value=object()), \
                    patch.object(curvature, 'hall_spectrum', side_effect=ReachedCalculation):
                curvature.hall_command(args)
        elif route == 'wannier':
            with patch.object(wannier, 'sha256', side_effect=ReachedCalculation):
                wannier.hall_command(args)

    def test_equal_endpoints_need_exactly_one_point_on_every_route(self):
        for route in ('pairs', 'waveder', 'curvature', 'wannier'):
            with self.subTest(route=route):
                with self.assertRaisesRegex(ValueError, 'MIN = MAX with N = 1'):
                    self.invoke(route, self.arguments(0., 0., 3))
                self.assertFalse(self.output.exists())

    def test_single_point_and_increasing_range_reach_calculation(self):
        for route in ('pairs', 'waveder', 'curvature', 'wannier'):
            for scan in ((0., 0., 1), (-.5, .5, 3)):
                with self.subTest(route=route, scan=scan), self.assertRaises(ReachedCalculation):
                    self.invoke(route, self.arguments(*scan))
        single, _ = pairs.scan_options(self.arguments(0., 0., 1))
        increasing, _ = pairs.scan_options(self.arguments(-.5, .5, 3))
        self.assertEqual(single.tolist(), [0.])
        self.assertEqual(increasing.tolist(), [-.5, 0., .5])

    def test_reverse_nonfinite_and_unequal_single_point_reject_on_every_route(self):
        for route in ('pairs', 'waveder', 'curvature', 'wannier'):
            for scan in ((1., 0., 2), (float('nan'), 1., 2),
                         (0., float('inf'), 2), (0., 1., 1), (0., 1., 0)):
                with self.subTest(route=route, scan=scan), self.assertRaises(ValueError):
                    self.invoke(route, self.arguments(*scan))
                self.assertFalse(self.output.exists())


if __name__ == '__main__':
    unittest.main()
