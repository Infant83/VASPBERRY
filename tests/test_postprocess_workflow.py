"""End-to-end study runner: validation, failure isolation, reuse and plot-only use.

The native executable is replaced only by an assigned analytic CSV fixture.
Import, integration, PROCAR association and plot commands run the real tools.
The fixture is not a material calculation or a native Fortran physics test.
"""
import contextlib
import importlib.util
import io
import json
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
import vaspberry_post as post
from postprocess_config import load_settings
from kubo_pairs import import_native_pairs, pair_hall_spectrum

spec = importlib.util.spec_from_file_location('postprocess_analytic_fixture',
    ROOT/'examples/features/procar-character/make_fixture.py')
fixture = importlib.util.module_from_spec(spec); spec.loader.exec_module(fixture)
REAL_RUN = subprocess.run


class PostprocessWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='vaspberry workflow ')
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name).resolve()
        self.inputs = self.root/'input files'; fixture.create(self.inputs)
        self.binary = self.root/'native executable'
        self.binary.write_text('#!/bin/sh\nexit 93\n') # Must never execute directly in these tests.
        self.binary.chmod(0o755)
        self.config = self.root/'study.ini'
        self.native_calls = []

    def settings(self, name='result', *, projection=False):
        text = '''[run]
wavecar = input files/WAVECAR
binary = native executable
output = NAME
mesh = 2 2
spin_mode = soc
energy_reference = synthetic analytic fixture zero
[hall]
mu = -1.2 0.2 8
reference = 0
temperatures = 0 300
[region first]
k_ids = 1
[differences]
contrast = first rest
[plot]
hall_regions = total first contrast
hall_temperatures = 300
hall_quantity = delta-sigma
'''.replace('output = NAME', 'output = '+name)
        if projection:
            text += '''
[projection]
bands = 1
axis = 1 0 0
[group lower]
ions = 1
[group upper]
ions = 2
'''
        self.config.write_text(text)
        return load_settings(self.config)

    def native_fixture(self, *, returncode=0, footer=True, launch_error=False):
        def execute(command, **kwargs):
            if command[0] != str(self.binary):
                return REAL_RUN(command, **kwargs)
            self.native_calls.append(command)
            if launch_error:
                raise OSError('synthetic native launch failure')
            kwargs['stdout'].write('analytic assigned pair export; no material calculation\n')
            kwargs['stderr'].write('native test diagnostic\n')
            text = (self.inputs/'PAIRS.csv').read_text()
            if not footer: text = text.replace('# result_status=PASS\n', '')
            (Path(kwargs['cwd'])/'PAIRS.csv').write_text(text)
            return SimpleNamespace(returncode=returncode)
        return execute

    def calculate(self, settings, **options):
        with contextlib.redirect_stdout(io.StringIO()), \
             patch.object(post.subprocess, 'run', side_effect=self.native_fixture(**options)):
            return post.run_calculation(settings)

    def assert_preflight_failure(self, settings, error=None, reuse=None):
        with patch.object(post.subprocess, 'run') as execute, contextlib.redirect_stdout(io.StringIO()):
            if error is None:
                with self.assertRaises((ValueError, OSError)):
                    post.run_calculation(settings, reuse)
            else:
                with self.assertRaisesRegex((ValueError, OSError), error):
                    post.run_calculation(settings, reuse)
            execute.assert_not_called()
        self.assertFalse(Path(settings['run']['output']).exists())

    def test_full_workflow_matches_direct_numerical_kernel_and_projection_analytics(self):
        settings = self.settings(projection=True)
        record = self.calculate(settings)
        out = Path(settings['run']['output'])
        self.assertEqual((record['status'], record['complete']), ('PASS', True))
        self.assertEqual(len(self.native_calls), 1)
        self.assertEqual([s['name'] for s in record['stages']],
            ['native-pairs', 'import-pairs', 'charge-hall', 'procar-character', 'projected-hall'])
        self.assertTrue(all(s['status'] == 'PASS' and s['exit_code'] == 0 for s in record['stages']))
        pairs = import_native_pairs(self.inputs/'PAIRS.csv', self.inputs/'WAVECAR', spin=1,
            spinor_components=2, spin_multiplicity=1,
            sampling={'kind': 'uniform_full_2d', 'mesh': [2, 2], 'plane_axes': [0, 1]},
            energy_reference='synthetic analytic fixture zero')
        rows, _ = pair_hall_spectrum(pairs, np.linspace(-1.2, .2, 8), [0., 300.], mu_reference=0.,
            region_spec={'regions': [{'name': 'first', 'k_ids': [1]}]}, differences=[('contrast', 'first', 'rest')])
        with np.load(out/'hall/conductivity.npz', allow_pickle=False) as saved:
            self.assertEqual(set(saved.files), set(rows[0]))
            for key in saved.files:
                np.testing.assert_array_equal(saved[key], np.array([row[key] for row in rows]))
        with np.load(out/'character-hall/character_hall.npz', allow_pickle=False) as saved:
            common = ((saved['region'] == 'total') & (saved['temperature_K'] == 0)
                      & np.isclose(saved['mu_eV'], 0, rtol=0, atol=1e-14))
            for group, component, expected in [('lower', 'charge', -.4*np.pi),
                    ('lower', 'plus', -.4*np.pi), ('lower', 'minus', 0.),
                    ('upper', 'minus', -1.2*np.pi), ('$unweighted', 'charge', -2*np.pi),
                    ('$unprojected_residual', 'charge', -.4*np.pi)]:
                mask = common & (saved['group'] == group) & (saved['component'] == component)
                self.assertEqual(mask.sum(), 1)
                self.assertAlmostEqual(float(saved['attribution_e2_over_h'][mask][0]), expected, places=12)
        self.assertEqual((out/'settings.ini').read_bytes(), self.config.read_bytes())
        self.assertEqual(json.loads((out/'groups.json').read_text()), settings['groups'])
        for name, checksum in record['outputs'].items():
            self.assertEqual(post.sha256(out/name), checksum)

    def test_existing_output_and_sentinel_are_preserved_without_launch(self):
        settings = self.settings()
        out = Path(settings['run']['output']); out.mkdir()
        sentinel = out/'prior work'; sentinel.write_bytes(b'untouched\x00data')
        with patch.object(post.subprocess, 'run') as execute, self.assertRaisesRegex(ValueError, 'exists'):
            post.run_calculation(settings)
        execute.assert_not_called()
        self.assertEqual(list(out.iterdir()), [sentinel])
        self.assertEqual(sentinel.read_bytes(), b'untouched\x00data')

    def test_missing_invalid_or_wrong_spin_source_fails_before_output_or_native(self):
        settings = self.settings()
        original = (self.inputs/'WAVECAR').read_bytes()
        for content in (b'', b'not a WAVECAR'):
            with self.subTest(content=content):
                (self.inputs/'WAVECAR').write_bytes(content)
                self.assert_preflight_failure(settings)
        (self.inputs/'WAVECAR').write_bytes(original)
        settings['run']['expected_nspin'] = 2
        self.assert_preflight_failure(settings, 'spin_mode')
        settings['run']['expected_nspin'] = 1
        (self.inputs/'WAVECAR').unlink()
        self.assert_preflight_failure(settings, 'WAVECAR does not exist')

    def test_mesh_count_and_geometry_fail_before_native(self):
        settings = self.settings()
        settings['run']['mesh'] = [2, 3]
        self.assert_preflight_failure(settings, 'mesh')
        settings['run']['mesh'] = [2, 2]
        original = (self.inputs/'WAVECAR').read_bytes()
        # Change the last k point while preserving a parseable source and NKPTS.
        for q in ((0., 0., 0.), (.25, .5, 0.), (.5, .5, .1)):
            with self.subTest(kpoint=q):
                changed = bytearray(original)
                struct.pack_into('<3d', changed, (2+3*3)*512+8, *q)
                (self.inputs/'WAVECAR').write_bytes(changed)
                self.assert_preflight_failure(settings)
        (self.inputs/'WAVECAR').write_bytes(original)
        settings['run']['plane_axes'] = [0, 2]
        self.assert_preflight_failure(settings)

    def test_missing_executable_launcher_or_projection_input_fails_before_native(self):
        settings = self.settings()
        self.binary.chmod(0o644)
        self.assert_preflight_failure(settings, 'not executable')
        self.binary.chmod(0o755)
        settings['run']['mpi_procs'] = 2
        settings['run']['mpi_launcher'] = str(self.root/'missing MPI launcher')
        self.assert_preflight_failure(settings, 'MPI launcher')
        settings = self.settings(projection=True)
        (self.inputs/'PROCAR').unlink()
        self.assert_preflight_failure(settings, 'PROCAR does not exist')

    def test_source_band_windows_and_map_band_fail_before_native(self):
        settings = self.settings()
        settings['hall']['pair_band_max'] = 3
        self.assert_preflight_failure(settings, 'pair_band_max')
        settings = self.settings(projection=True)
        settings['projection']['bands'] = [3]
        self.assert_preflight_failure(settings, 'bands')
        settings['projection']['bands'] = [1]
        settings['plot']['map_band'] = 3
        self.assert_preflight_failure(settings, 'map_band')

    def test_native_exit_and_launch_failure_remain_recorded_without_hall(self):
        for index, options in enumerate(({'returncode': 7}, {'launch_error': True})):
            with self.subTest(options=options):
                settings = self.settings('failure'+str(index))
                with self.assertRaises((ValueError, OSError)):
                    self.calculate(settings, **options)
                out = Path(settings['run']['output'])
                record = json.loads((out/'run.json').read_text())
                self.assertEqual((record['status'], record['complete']), ('FAILED', False))
                self.assertEqual(len(record['stages']), 1)
                self.assertEqual(record['stages'][0]['status'], 'FAILED')
                self.assertIn('error', record)
                self.assertFalse((out/'hall').exists())
                self.assertFalse((out/'pairs').exists())
                self.assertTrue((out/'native/stdout.log').is_file())
                self.assertTrue((out/'native/stderr.log').is_file())
                if 'returncode' in options:
                    self.assertEqual(record['stages'][0]['exit_code'], 7)
                    self.assertIn('native test diagnostic', (out/'native/stderr.log').read_text())

    def test_exit_zero_with_incomplete_native_export_is_a_failed_import(self):
        settings = self.settings()
        with self.assertRaisesRegex(ValueError, 'import-pairs failed'):
            self.calculate(settings, footer=False)
        out = Path(settings['run']['output'])
        record = json.loads((out/'run.json').read_text())
        self.assertEqual(record['status'], 'FAILED')
        self.assertEqual([s['status'] for s in record['stages']], ['PASS', 'FAILED'])
        self.assertFalse((out/'hall').exists())
        self.assertFalse((out/'pairs').exists())
        self.assertIn('complete native physical output required', (out/'logs/import-pairs.stderr.log').read_text())

    def test_reuse_skips_native_without_binary_and_retains_exact_numerical_arrays(self):
        first = self.settings('first'); self.calculate(first)
        first_dir = Path(first['run']['output'])
        second = self.settings('second')
        self.binary.unlink()
        second['run']['mpi_procs'] = 8
        second['run']['mpi_launcher'] = '/missing/launcher'
        with contextlib.redirect_stdout(io.StringIO()), patch.object(post.subprocess, 'run', side_effect=REAL_RUN) as execute:
            record = post.run_calculation(second, first_dir)
        out = Path(second['run']['output'])
        self.assertEqual(len(self.native_calls), 1)
        self.assertEqual([s['name'] for s in record['stages']], ['charge-hall'])
        self.assertEqual(execute.call_count, 1)
        self.assertEqual(record['reused_pair_cache'], str(first_dir/'pairs'))
        self.assertFalse((out/'native').exists())
        for name in ('pairs/pairs.npz', 'pairs/pairs.json', 'hall/conductivity.csv', 'hall/conductivity.dat'):
            self.assertEqual((first_dir/name).read_bytes(), (out/name).read_bytes())
        with np.load(first_dir/'hall/conductivity.npz') as first_data, np.load(out/'hall/conductivity.npz') as second_data:
            for key in first_data.files: np.testing.assert_array_equal(first_data[key], second_data[key])

    def test_reuse_rejects_changed_source_settings_bytes_or_incomplete_run_without_launch(self):
        first = self.settings('first'); self.calculate(first)
        first_dir = Path(first['run']['output'])
        second = self.settings('second')
        second['run']['energy_reference'] = 'different zero'
        self.assert_preflight_failure(second, 'cannot change source setting', first_dir)
        second = self.settings('second')
        with (self.inputs/'WAVECAR').open('ab') as handle: handle.write(b'new source bytes')
        self.assert_preflight_failure(second, 'same WAVECAR', first_dir)
        (self.inputs/'WAVECAR').write_bytes((self.inputs/'WAVECAR').read_bytes()[:-16])
        record_path = first_dir/'run.json'
        record = json.loads(record_path.read_text()); record.update(status='FAILED', complete=False)
        record_path.write_text(json.dumps(record))
        self.assert_preflight_failure(second, 'did not complete', first_dir)

    def test_reuse_rejects_changed_pair_cache_or_metadata_before_output(self):
        first = self.settings('first'); self.calculate(first)
        first_dir = Path(first['run']['output']); second = self.settings('second')
        for relative in ('pairs/pairs.npz', 'pairs/pairs.json'):
            path = first_dir/relative; original = path.read_bytes()
            path.write_bytes(original+b' ')
            self.assert_preflight_failure(second, 'saved output changed', first_dir)
            path.write_bytes(original)

    def test_plot_only_needs_saved_results_and_preserves_independent_band_scope(self):
        settings = self.settings(projection=True); self.calculate(settings)
        out = Path(settings['run']['output'])
        shutil.rmtree(self.inputs); self.config.unlink(); self.binary.unlink()
        with contextlib.redirect_stdout(io.StringIO()), patch.object(post.subprocess, 'run', side_effect=REAL_RUN) as execute:
            record = post.plot_results(out, temperature=300., group='upper', band=2, quantity='delta-sigma')
        self.assertEqual((record['status'], record['complete']), ('PASS', True))
        self.assertEqual(execute.call_count, 2)
        self.assertEqual([s['name'] for s in record['stages']], ['charge-hall-figures', 'character-figures'])
        self.assertEqual(record['preferences']['map_band'], 2)
        metadata = json.loads((out/'character-hall/character_hall.json').read_text())
        self.assertEqual(metadata['selected_bands'], [1])
        for relative in ('charge-hall/hall', 'character/character', 'character/character_hall'):
            for suffix in ('png', 'pdf', 'svg'):
                self.assertGreater((out/'figures'/(relative+'.'+suffix)).stat().st_size, 100)
        self.assertEqual(len(self.native_calls), 1)

    def test_plot_rejects_tampered_hall_table_before_creating_figures(self):
        settings = self.settings(); self.calculate(settings)
        out = Path(settings['run']['output'])
        with (out/'hall/conductivity.csv').open('a') as handle: handle.write('\n')
        with patch.object(post.subprocess, 'run') as execute, self.assertRaisesRegex(ValueError, 'saved output changed'):
            post.plot_results(out)
        execute.assert_not_called()
        self.assertFalse((out/'figures').exists())

    def test_plot_rejects_map_band_above_saved_band_count_before_output(self):
        settings = self.settings(projection=True); self.calculate(settings)
        out = Path(settings['run']['output'])
        shutil.rmtree(self.inputs); self.config.unlink(); self.binary.unlink()
        with patch.object(post.subprocess, 'run') as execute, self.assertRaisesRegex(ValueError, 'band'):
            post.plot_results(out, band=3)
        execute.assert_not_called()
        self.assertFalse((out/'figures').exists())

    def test_existing_figure_directory_is_preserved_without_plot_launch(self):
        settings = self.settings(); self.calculate(settings)
        out = Path(settings['run']['output']); figures = out/'figures'; figures.mkdir()
        marker = figures/'existing.png'; marker.write_bytes(b'prior figure')
        with patch.object(post.subprocess, 'run') as execute, self.assertRaisesRegex(ValueError, 'figure directory exists'):
            post.plot_results(out)
        execute.assert_not_called()
        self.assertEqual(list(figures.iterdir()), [marker])
        self.assertEqual(marker.read_bytes(), b'prior figure')


if __name__ == '__main__':
    unittest.main()
