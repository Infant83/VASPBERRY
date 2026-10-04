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


def optical_export_text():
    text = export_text().replace(native.TRACE_SCHEMA, native.PAW_TRACE_SCHEMA)
    metadata = dict(kubo_source='WAVEDER', source_operator=native.WAVEDER_OPERATOR,
        result_kind='T0_INSULATING_OCCUPIED_BUNDLE_TRACE', complete_occupied_bundle='true',
        source_ndbands=2, spinor_components=1, physical_spin_multiplicity=1,
        occupied_bands_spin_1=2, occupied_bands_spin_2=1, producer_cluster_threshold_eV='0.002')
    for key in ('band_min', 'band_max', 'band_rank'):
        text = '\n'.join(line for line in text.split('\n') if not line.startswith('# '+key+'='))
    return ''.join(f'# {key}={value}\n' for key, value in metadata.items()) + text


def selected_optical_text(kind="trace", *, ids=(1, 3), nb=4, nd=3):
    """Independent schema fixture; gaps/values do not come from the producer."""
    trace = kind == "trace"
    full = trace and len(ids) == nb
    meta = dict(schema=native.PAW_BUNDLE_SCHEMA if trace else native.PAW_BAND_SCHEMA,
        normalization="STANDARD_MINUS_TWO_IM", kubo_source="WAVEDER",
        source_operator=native.WAVEDER_OPERATOR, source_nbands=nb, source_ndbands=nd,
        source_nkpoints=2, source_nspin=1, spinor_components=2, physical_spin_multiplicity=1,
        expected_rows=2*(1 if trace else len(ids)), band_ids=','.join(map(str, ids)),
        band_min=min(ids), band_max=max(ids), band_rank=len(ids),
        result_kind="ISOLATED_SELECTED_BUNDLE_TRACE" if trace else "ISOLATED_SELECTED_BAND_CURVATURE",
        occupation_weighting="NONE", required_pair_coverage="PASS",
        intermediate_bands="EXTERNAL_TO_SELECTED_BUNDLE_WITHIN_SOURCE_NBANDS" if trace else "ALL_SOURCE_BANDS_EXCEPT_SELF",
        producer_cluster_threshold_eV="0.002", no_external_states=str(full).lower(),
        result_status="INCOMPLETE")
    if full:
        meta['zero_trace_scope'] = 'TRUNCATED_WAVEDER_BASIS'
    out = io.StringIO()
    for key, value in meta.items():
        out.write(f'# {key}={value}\n')
    writer = csv.writer(out, lineterminator='\n')
    writer.writerow(native.TRACE_COLUMNS if trace else native.BAND_COLUMNS)
    for k in (1, 2):
        if trace:
            writer.writerow([1, k, k*.25, 0., 0., 0. if full else .25*k, 'NA' if full else 1.])
        else:
            for band in ids:
                writer.writerow([1, k, band, k*.25, 0., 0., float(band), .25*k, 1.])
    return out.getvalue() + '# result_status=PASS\n'


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

    def test_paw_trace_preserves_operator_and_different_collinear_occupations(self):
        meta, rows = self.read(optical_export_text())
        self.assertEqual(meta['schema'], native.PAW_TRACE_SCHEMA)
        self.assertEqual(meta['source_operator'], native.WAVEDER_OPERATOR)
        self.assertEqual(meta['completion_validation'], 'terminal_pass_and_full_coverage')
        self.assertEqual((meta['occupied_bands_spin_1'], meta['occupied_bands_spin_2']), ('2', '1'))
        self.assertEqual(len(rows), 4)
        self.assertNotIn('operator', meta)  # No canonical-operator label is invented.

    def test_paw_trace_rejects_missing_coverage_or_changed_operator_and_threshold(self):
        text = optical_export_text()
        replacements = [('source_ndbands=2', 'source_ndbands=1'),
                        ('occupied_bands_spin_2=1', 'occupied_bands_spin_2=3'),
                        ('source_operator='+native.WAVEDER_OPERATOR, 'source_operator=canonical'),
                        ('complete_occupied_bundle=true', 'complete_occupied_bundle=false'),
                        ('producer_cluster_threshold_eV=0.002', 'producer_cluster_threshold_eV=0.0001'),
                        ('physical_spin_multiplicity=1', 'physical_spin_multiplicity=2'),
                        ('# result_status=PASS\n', ''), ('0.125,2.0', '0.125,0.001')]
        for old, new in replacements:
            self.assertIn(old, text)
            with self.subTest(field=old), self.assertRaises(ValueError):
                self.read(text.replace(old, new))

    def test_paw_trace_cannot_be_imported_as_canonical_bundle_hall(self):
        from kubo_pairs import bundle_hall_spectrum
        self.path.write_text(optical_export_text())
        with self.assertRaisesRegex(ValueError, 'kubo-hall --run-dir'):
            bundle_hall_spectrum(self.path, 'not-opened-WAVECAR', [0.], occupied=2,
                sampling={'kind': 'uniform_full_2d', 'mesh': [2, 1]},
                energy_reference='test', mu_reference=0.)

    def test_selected_optical_trace_and_disjoint_band_rows_keep_exact_selection(self):
        for kind in ('trace', 'band'):
            meta, rows = self.read(selected_optical_text(kind), kind)
            self.assertEqual(native.selected_band_ids(meta), [1, 3])
            self.assertEqual(len(rows), 2 if kind == 'trace' else 4)
            self.assertEqual(meta['completion_validation'], 'terminal_pass_and_full_coverage')
            self.assertNotIn('sheet_hall_e2_over_h', meta)
        text = selected_optical_text('band')
        with self.assertRaisesRegex(ValueError, 'coverage'):
            self.read(text.replace('1,1,3,', '1,1,2,'), 'band')

    def test_selected_optical_reverse_coverage_and_empty_complement(self):
        self.read(selected_optical_text(ids=(3, 4), nd=2))
        self.read(selected_optical_text('band', ids=(4,), nd=3), 'band')
        for kind in ('trace', 'band'):
            with self.assertRaisesRegex(ValueError, 'missing high-high'):
                self.read(selected_optical_text(kind, ids=(3,), nd=2), kind)
        _, rows = self.read(selected_optical_text(ids=(1, 2, 3, 4), nd=2))
        self.assertTrue(all(r['omega_z_A2'] == 0. and r['min_external_gap_eV'] is None for r in rows))

    def test_selected_optical_metadata_cannot_claim_occupation_or_missing_pairs(self):
        text = selected_optical_text()
        for old, new in [('band_ids=1,3', 'band_ids=1,2'),
                         ('band_ids=1,3', 'band_ids=3,1'),
                         ('occupation_weighting=NONE', 'occupation_weighting=FERMI'),
                         ('required_pair_coverage=PASS', 'required_pair_coverage=PARTIAL'),
                         ('intermediate_bands=EXTERNAL_TO_SELECTED_BUNDLE_WITHIN_SOURCE_NBANDS', 'intermediate_bands=1:3'),
                         ('no_external_states=false', 'no_external_states=true')]:
            with self.subTest(field=old), self.assertRaises(ValueError):
                self.read(text.replace(old, new))

    def test_selected_optical_trace_is_not_canonical_bundle_hall(self):
        from kubo_pairs import bundle_hall_spectrum
        self.path.write_text(selected_optical_text())
        with self.assertRaisesRegex(ValueError, 'kubo-hall --run-dir'):
            bundle_hall_spectrum(self.path, 'not-opened-WAVECAR', [0.], occupied=2,
                sampling={'kind': 'uniform_full_2d', 'mesh': [2, 1]},
                energy_reference='test', mu_reference=0.)


if __name__ == "__main__":
    unittest.main()
