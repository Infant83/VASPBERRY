"""WAVEDER is the execution default; source failures must not select another operator.

These tests parse actual Fortran-record WAVEDER fixtures and evaluate their
optical matrices. Only WAVECAR record reading is mocked, as in its existing
same-run adapter tests; this is a source-routing regression, not a material test.
"""
import contextlib
import io
import json
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
import test_waveder_hall as optical_fixture
import kubo_hall_workflow as workflow
from postprocess_config import load_settings
import vaspberry_kubo as cli
import vaspberry_post as post


class KuboSourceDefaultsTests(unittest.TestCase):
    def setUp(self):
        self.fixture = optical_fixture.WaveDerHallTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.root = self.fixture.root
        # A rectangular 3 x 2 optical matrix is sufficient for the two-band
        # occupied bundle; absent empty-empty columns are never fabricated.
        optical_fixture.write_waveder(self.fixture.run/'WAVEDER', self.fixture.c[..., :2])

    def cli_args(self, output='cli', extra=()):
        return cli.parser().parse_args(['kubo-hall', *self.fixture.cli_arguments(output), *extra])

    def settings(self, output='study'):
        path = self.root/'study.ini'
        path.write_text(f'''[run]
input_dir = run
wavecar = run/WAVECAR
output = {output}
mesh = 2 2
spin_mode = soc
energy_reference = unchanged fixture VASP zero
[hall]
occupied = 2
mu = -0.2 0.2 3
reference = 0.25
temperatures = 0
''')
        return load_settings(path)

    def test_general_cli_defaults_to_real_rectangular_optical_matrices(self):
        args = self.cli_args()
        self.assertEqual(args.kubo_source, 'waveder')
        with patch.object(workflow.subprocess, 'run') as native:
            meta = workflow.kubo_hall_command(args)
        native.assert_not_called()
        self.assertEqual(meta['method'], 'standard_waveder_occupied_bundle_T0')
        self.assertEqual(meta['source_metadata']['source_ndbands'], 2)
        with np.load(args.output_dir/'conductivity.npz') as actual:
            self.assertGreater(np.max(np.abs(actual['sigma_e2_over_h'])), 0.)
        optical_fixture.write_waveder(self.fixture.run/'WAVEDER', np.zeros_like(self.fixture.c[..., :2]))
        workflow.kubo_hall_command(self.cli_args('zero'))
        with np.load(self.root/'zero/conductivity.npz') as actual:
            np.testing.assert_array_equal(actual['sigma_e2_over_h'], 0.)

    def test_default_settings_use_optical_route_without_native_binary(self):
        settings = self.settings()
        self.assertIsNone(settings['run']['binary'])
        self.assertEqual(settings['run']['kubo_source'], 'waveder')
        with patch.object(post.subprocess, 'run') as native, contextlib.redirect_stdout(io.StringIO()):
            record = post.run_calculation(settings)
        native.assert_not_called()
        self.assertEqual((record['status'], record['complete']), ('PASS', True))
        self.assertEqual([s['name'] for s in record['stages']], ['waveder-charge-hall'])
        self.assertEqual(set(record['inputs']), {'waveder', 'wavecar', 'incar', 'outcar'})
        out = Path(settings['run']['output'])
        self.assertFalse((out/'pairs').exists())
        meta = json.loads((out/'hall/conductivity.json').read_text())
        self.assertEqual(meta['method'], 'standard_waveder_occupied_bundle_T0')
        self.assertEqual(meta['provenance']['kubo_source'], 'waveder')

    def assert_rejected_without_fallback(self, args, settings, error):
        with patch.object(workflow.subprocess, 'run') as native:
            with self.assertRaisesRegex((ValueError, OSError), error):
                workflow.kubo_hall_command(args)
            with self.assertRaisesRegex((ValueError, OSError), error):
                post.run_calculation(settings)
        native.assert_not_called()
        self.assertFalse(args.output_dir.exists())
        self.assertFalse(Path(settings['run']['output']).exists())

    def test_missing_waveder_errors_before_output_without_fallback(self):
        (self.fixture.run/'WAVEDER').unlink()
        self.assert_rejected_without_fallback(self.cli_args(), self.settings(), 'WAVEDER input is missing')

    def test_invalid_existing_waveder_errors_without_fallback(self):
        (self.fixture.run/'WAVEDER').write_bytes(b'invalid optical matrix')
        self.assert_rejected_without_fallback(self.cli_args(), self.settings(), '.')

    def test_finite_temperature_is_not_rerouted_to_wavecar(self):
        args = self.cli_args(extra=['--temperatures', '0', '300'])
        settings = self.settings(); settings['hall']['temperatures'] = [0., 300.]
        self.assert_rejected_without_fallback(args, settings, 'only insulating T=0')

    def test_scan_into_a_band_is_not_rerouted_to_wavecar(self):
        args = self.cli_args(extra=['--mu-max', '3'])
        settings = self.settings(); settings['hall']['mu_max'] = 3.
        self.assert_rejected_without_fallback(args, settings, 'global insulating gap')

    def test_general_cli_dispatch_uses_default_handler(self):
        with contextlib.redirect_stdout(io.StringIO()) as stdout:
            cli.main(['kubo-hall', *self.fixture.cli_arguments()])
        self.assertEqual(json.loads(stdout.getvalue())['schema'], 'vaspberry.hall-spectrum')

    def test_source_selection_rejects_misspelling(self):
        self.settings()
        path = self.root/'study.ini'
        path.write_text(path.read_text().replace('[run]', '[run]\nkubo_source = auto'))
        with self.assertRaisesRegex(ValueError, 'kubo_source'):
            load_settings(path)

    def test_waveder_projection_or_pair_cutoff_is_not_silently_ignored(self):
        settings = self.settings()
        settings['hall']['pair_band_max'] = 2
        with self.assertRaisesRegex(ValueError, 'pair_band_max'):
            post.run_calculation(settings)
        settings['hall']['pair_band_max'] = None
        settings['projection'] = {'bands': [1]}
        with self.assertRaisesRegex(ValueError, r'\[projection\]'):
            post.run_calculation(settings)
        self.assertFalse(Path(settings['run']['output']).exists())

    def test_wavecar_settings_do_not_silently_ignore_optical_occupied_count(self):
        settings = self.settings()
        settings['run']['kubo_source'] = 'wavecar'
        with patch.object(post.subprocess, 'run') as native, \
             self.assertRaisesRegex(ValueError, 'occupied/bands is only supported'):
            post.run_calculation(settings)
        native.assert_not_called()
        self.assertFalse(Path(settings['run']['output']).exists())

    def test_waveder_cli_does_not_ignore_native_execution_options(self):
        for flags in (['--binary', 'not-used'], ['--mpi-procs', '2'],
                      ['--mpi-launcher', 'custom-mpiexec']):
            args = self.cli_args(extra=flags)
            with self.subTest(flags=flags), patch.object(workflow.subprocess, 'run') as native, \
                 self.assertRaisesRegex(ValueError, 'native MPI options apply only'):
                workflow.kubo_hall_command(args)
            native.assert_not_called()
            self.assertFalse(args.output_dir.exists())

    def test_waveder_settings_do_not_ignore_native_execution_controls(self):
        for key, value in (('binary', '/not-used'), ('mpi_procs', 2),
                           ('mpi_launcher', 'custom-mpiexec')):
            settings = self.settings()
            settings['run'][key] = value
            with self.subTest(setting=key), patch.object(post.subprocess, 'run') as native, \
                 self.assertRaisesRegex(ValueError, 'native MPI settings apply only'):
                post.run_calculation(settings)
            native.assert_not_called()
            self.assertFalse(Path(settings['run']['output']).exists())

    def test_wavecar_settings_reject_explicit_optical_directory(self):
        self.settings()
        path = self.root/'study.ini'
        path.write_text(path.read_text().replace('[run]',
            '[run]\nkubo_source = wavecar\noptical_run_dir = unused-optical-run'))
        with self.assertRaisesRegex(ValueError, 'optical_run_dir applies only'):
            load_settings(path)
