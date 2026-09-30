"""Corruption and completion checks independent of the native CSV writer."""
import csv
import io
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import native_kubo_csv as native


def export_text(kind="trace", *, spins=2, legacy=False):
    trace = kind == "trace"
    schema = (native.LEGACY_TRACE_SCHEMA if legacy else native.TRACE_SCHEMA) if trace else (
        native.LEGACY_BAND_SCHEMA if legacy else native.BAND_SCHEMA)
    meta = dict(schema=schema, normalization="STANDARD_MINUS_TWO_IM",
                band_min=1, band_max=2, band_rank=2, source_nbands=3)
    if not legacy:
        meta.update(source_nkpoints=2, source_nspin=spins,
                    expected_rows=2*spins*(1 if trace else 2), result_status="INCOMPLETE")
    elif trace:
        meta["result_status"] = "PASS"
    out = io.StringIO()
    for key, value in meta.items():
        out.write(f"# {key}={value}\n")
    writer = csv.writer(out, lineterminator="\n")
    writer.writerow(native.TRACE_COLUMNS if trace else native.BAND_COLUMNS)
    for spin in range(1, spins+1):
        for k in (1, 2):
            if trace:
                writer.writerow([spin, k, .25*k, 0., 0., .125*spin*k, 2.])
            else:
                for band in (1, 2):
                    writer.writerow([spin, k, band, .25*k, 0., 0., -3.+band, .125*spin*k, 1.])
    if not legacy:
        out.write("# result_status=PASS\n")
    return out.getvalue()


class NativeCurvatureCSVTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / "KUBO.csv"

    def read(self, text, kind="trace"):
        self.path.write_text(text)
        return native.read_curvature_csv(self.path, kind=kind)

    def test_current_trace_and_band_cover_both_spins(self):
        for kind, expected in (("trace", 4), ("band", 8)):
            with self.subTest(kind=kind):
                meta, rows = self.read(export_text(kind), kind)
                self.assertEqual(meta["completion_validation"], "terminal_pass_and_full_coverage")
                self.assertEqual(len(rows), expected)
                self.assertEqual({row["spin"] for row in rows}, {1, 2})

    def test_legacy_compatibility_is_explicitly_unverified(self):
        for kind in ("trace", "band"):
            meta, rows = self.read(export_text(kind, legacy=True), kind)
            self.assertEqual(meta["completion_validation"], "legacy_unverified")
            self.assertTrue(rows)

    def test_missing_terminal_footer_or_appended_data_is_rejected(self):
        text = export_text()
        for bad in (text.rsplit("# result_status=PASS", 1)[0],
                    text + "# interrupted copy\n", text + "1,1,0,0,0,0,2\n",
                    text.replace("# result_status=INCOMPLETE", "# result_status=PASS")):
            with self.subTest(text=bad), self.assertRaisesRegex(ValueError, "result_status"):
                self.read(bad)

    def test_copied_footer_does_not_make_truncated_or_duplicate_rows_complete(self):
        for kind in ("trace", "band"):
            lines = export_text(kind).splitlines()
            data = [i for i, line in enumerate(lines) if line[0].isdigit()]
            cases = [lines[:data[-1]] + lines[data[-1]+1:],
                     [line for line in lines if not line.startswith("2,")],
                     lines[:data[-1]] + [lines[data[0]]] + lines[data[-1]+1:]]
            for bad in cases:
                with self.subTest(kind=kind), self.assertRaisesRegex(ValueError, "count|coverage"):
                    self.read("\n".join(bad)+"\n", kind)

    def test_nonfinite_values_in_unselected_spin_or_band_are_rejected(self):
        for kind in ("trace", "band"):
            for value in ("NaN", "Infinity", "-Infinity"):
                lines = export_text(kind).splitlines()
                fields = lines[-2].split(",")
                for column in range(2 if kind == "trace" else 3, len(fields)):
                    changed = fields.copy(); changed[column] = value
                    bad = lines[:-2] + [",".join(changed), lines[-1]]
                    with self.subTest(kind=kind, column=column, value=value), self.assertRaisesRegex(ValueError, "nonfinite"):
                        self.read("\n".join(bad)+"\n", kind)

    def test_bad_counts_ranges_indices_and_partial_rows_are_rejected(self):
        text = export_text()
        changes = [("# source_nspin=2", "# source_nspin=3"),
                   ("# source_nkpoints=2", "# source_nkpoints=nan"),
                   ("# expected_rows=4", "# expected_rows=999999999999999999999999"),
                   ("# source_nkpoints=2\n", ""),
                   ("# band_rank=2", "# band_rank=3"),
                   ("# band_max=2", "# band_max=4"),
                   ("1,1,0.25", "1,0,0.25"),
                   ("1,1,0.25", "1,1.5,0.25"),
                   ("2,2,0.5,0.0,0.0,0.5,2.0", "2,2,0.5,0.0,0.0,0.5")]
        for old, new in changes:
            self.assertIn(old, text)
            with self.subTest(change=(old,new)), self.assertRaises(ValueError):
                self.read(text.replace(old, new))

    def test_conflicting_metadata_cannot_override_source_counts(self):
        with self.assertRaisesRegex(ValueError, "conflicting"):
            self.read(export_text().replace("# source_nspin=2", "# source_nspin=1\n# source_nspin=2"))

    def test_full_retained_trace_has_undefined_gap_not_nonfinite_curvature(self):
        text=export_text().replace('# band_max=2','# band_max=3').replace('# band_rank=2','# band_rank=3')
        text='# no_external_states=true\n# zero_trace_scope=TRUNCATED_WAVECAR_BASIS\n'+text
        lines=text.splitlines()
        for i,line in enumerate(lines):
            if line[0].isdigit():
                fields=line.split(',');fields[-2:]=['0.0','NA'];lines[i]=','.join(fields)
        text='\n'.join(lines)+'\n'
        meta,rows=self.read(text)
        self.assertEqual(meta['completion_validation'],'terminal_pass_and_full_coverage')
        self.assertTrue(all(row['omega_z_A2']==0. and row['min_external_gap_eV'] is None for row in rows))
        for bad in (text.replace('no_external_states=true','no_external_states=false'),
                    text.replace('# band_max=3','# band_max=2'),
                    text.replace('0.0,NA','0.1,NA'),text.replace('0.0,NA','nan,NA'),
                    text.replace('0.0,NA','0.0,nan')):
            with self.subTest(text=bad),self.assertRaises(ValueError):self.read(bad)

    def test_historical_public_trace_and_band_files_remain_readable(self):
        base = ROOT / "examples/features/kubo-curvature/reference"
        for path, kind, count in ((base/"bundle/KUBO_mesh.csv", "trace",144),
                                  (base/"KUBO.csv", "band",96)):
            meta, rows = native.read_curvature_csv(path, kind=kind)
            self.assertEqual(meta["completion_validation"], "legacy_unverified")
            self.assertEqual(len(rows), count)


if __name__ == "__main__":
    unittest.main()
