"""Real trivial/nontrivial fields and structural plot-only input guards."""
import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "examples/features/z2/run.py"
SPEC = importlib.util.spec_from_file_location("z2_comparison_validator", SOURCE)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)
MOS2 = ROOT / "examples/features/z2/mos2/reference/Z2_FIELD.csv"
BI = ROOT / "examples/materials/bi-spin-hall/reference/z2/Z2_FIELD.csv"


class Z2ComparisonTests(unittest.TestCase):
    def test_actual_mos2_and_fresh_bi_parities_from_native_rows(self):
        for path, rank, parity, sums in [(MOS2, 18, 0, (0, 0)), (BI, 10, 1, (-3, 3))]:
            with self.subTest(material=path):
                _, metadata, summary = MODULE.validate_rows(path)
                with path.open() as stream:
                    rows = list(csv.DictReader(line for line in stream if not line.startswith('#')))
                top = sum(int(row['nfield_int']) for row in rows if float(row['q2']) % 1 < .5)
                bottom = sum(int(row['nfield_int']) for row in rows if float(row['q2']) % 1 >= .5)
                self.assertEqual((top, bottom), sums)
                self.assertEqual(top % 2, parity)
                self.assertEqual(bottom % 2, parity)
                self.assertEqual(summary['z2'], parity)
                self.assertEqual((summary['mesh_nx'], summary['mesh_ny'], len(rows)), (12, 12, 144))
                self.assertEqual(int(metadata['band_rank']), rank)

    def test_incompatible_occupied_bundle_metadata_is_rejected(self):
        original = MOS2.read_text()
        cases = {
            'odd_rank': {'band_max': '17', 'band_rank': '17'},
            'inconsistent_rank': {'band_rank': '20'},
            'nonoccupied_window': {'band_min': '3', 'band_rank': '16'},
            'zero_rank': {'band_max': '0', 'band_rank': '0'},
            'scalar_states': {'spinor_components': '1'},
            'invalid_integer': {'band_rank': '18.5'},
        }
        with tempfile.TemporaryDirectory(prefix='vaspberry-z2-contract-') as temporary:
            for name, changes in cases.items():
                with self.subTest(case=name):
                    lines = []
                    for line in original.splitlines():
                        if line.startswith('# ') and '=' in line:
                            key = line[2:].split('=', 1)[0]
                            if key in changes:
                                line = '# ' + key + '=' + changes[key]
                        lines.append(line)
                    path = Path(temporary) / (name + '.csv')
                    path.write_text('\n'.join(lines) + '\n')
                    with self.assertRaisesRegex(ValueError, 'bundle metadata|positive even rank'):
                        MODULE.validate_rows(path)


if __name__ == '__main__':
    unittest.main()
