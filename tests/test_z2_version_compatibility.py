"""Keep schema-2 golden data immutable while accepting the new producer."""
import tempfile
import unittest
from pathlib import Path

from compare_z2_fields import Z2ComparisonError, compare_z2_fields, load_z2_field

ROOT = Path(__file__).resolve().parents[1]
GOLDEN = ROOT / "examples/Bi_Z2/reference-v1.2.0-12x12/Z2_FIELD.csv"


class Z2VersionCompatibilityTests(unittest.TestCase):
    def test_current_producer_matches_unchanged_historical_golden(self):
        original = GOLDEN.read_text()
        with tempfile.TemporaryDirectory() as work:
            path = Path(work) / "current.csv"
            for version in ("1.3.0", "1.4.0", "1.4.1", "1.4.2"):
                with self.subTest(version=version):
                    path.write_text(original.replace("# vaspberry_version=1.2.0", f"# vaspberry_version={version}"))
                    compare_z2_fields(load_z2_field(GOLDEN), load_z2_field(path), rtol=1e-11, atol=1e-12)
        self.assertEqual(GOLDEN.read_text(), original)

    def test_unreviewed_producer_version_is_rejected(self):
        with tempfile.TemporaryDirectory() as work:
            path = Path(work) / "future.csv"
            path.write_text(GOLDEN.read_text().replace("# vaspberry_version=1.2.0", "# vaspberry_version=9.9.9"))
            with self.assertRaises(Z2ComparisonError):
                load_z2_field(path)

    def test_other_metadata_comparisons_are_preserved(self):
        with tempfile.TemporaryDirectory() as work:
            path = Path(work) / "changed.csv"
            path.write_text(GOLDEN.read_text().replace("# nfield_note=", "# nfield_note=CHANGED "))
            with self.assertRaises(Z2ComparisonError):
                compare_z2_fields(load_z2_field(GOLDEN), load_z2_field(path), rtol=1e-11, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
