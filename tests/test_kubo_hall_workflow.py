"""CLI failure isolation, reusable scans and independently usable Hall formats."""
import contextlib
import csv
import io
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
import kubo_hall_workflow as workflow
import kubo_pairs as pairs
import plot_hall
import vaspberry_kubo as cli
from berry_data import write_hall


class HallWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.binary = self.root/'vaspberry'; self.binary.write_bytes(b'mocked native binary')
        self.wavecar = self.root/'WAVECAR'; self.wavecar.write_bytes(b'mocked WAVECAR')
        self.q = np.array([[0., 0., 0.], [0., .5, 0.], [.5, 0., 0.], [.5, .5, 0.]])
        self.energy = np.tile([-1., 2.], (4, 1))
        self.fake_wavecar = SimpleNamespace(energies=self.energy, kpoints=self.q,
            header=SimpleNamespace(ispin=1, lattice=np.eye(3), reciprocal=2*np.pi*np.eye(3)),
            coefficients=lambda k, b: np.ones((1, 1, 1)))
        self.addCleanup(patch.stopall)
        patch.object(pairs, 'Wavecar', return_value=self.fake_wavecar).start()

    def native_csv(self, footer=True):
        meta = dict(schema=pairs.NATIVE_SCHEMA, result_kind='UNORDERED_INTERBAND_NUMERATORS',
            normalization='STANDARD_MINUS_TWO_IM', operator=pairs.OPERATOR,
            berry_connection=pairs.CONNECTION, occupation_weighting='NONE',
            denominator_weighting='NONE', pair_order='n_lt_m', gap_definition='ABS_EN_MINUS_EM',
            numerator_units='eV^2*Angstrom^2', components='yz,zx,xy', reciprocal_convention='2pi',
            index_base=1, source_nbands=2, source_nkpoints=4, source_nspin=1,
            spinor_components=1, pairs_per_k=1, expected_rows=4)
        for i in range(3):
            meta[f'lattice_A_{i+1}'] = ','.join(map(str, np.eye(3)[i]))
            meta[f'reciprocal_inv_A_{i+1}'] = ','.join(map(str, (2*np.pi*np.eye(3))[i]))
        out = io.StringIO()
        for key, value in meta.items():
            out.write(f'# {key}={value}\n')
        writer = csv.writer(out)
        writer.writerow(['spin', 'k_index', 'n_band', 'm_band', 'kx_frac', 'ky_frac', 'kz_frac',
                         'energy_n_eV', 'energy_m_eV', 'numerator_yz_eV2_A2',
                         'numerator_zx_eV2_A2', 'numerator_xy_eV2_A2', 'gap_eV'])
        for k, q in enumerate(self.q):
            writer.writerow([1, k+1, 1, 2, *q, -1., 2., 0., 0., .1*(k+1), 3.])
        if footer:
            out.write('# result_status=PASS\n')
        return out.getvalue()

    def source_arguments(self):
        return ['--wavecar', str(self.wavecar), '--spinor-components', '1',
                '--spin-multiplicity', '1', '--mesh', '2', '2', '--energy-reference', 'fixture zero']

    def scan_arguments(self):
        # A distant reference must not be appended to the requested output grid.
        return ['--mu-min', '-.2', '--mu-max', '.2', '--mu-num', '3',
                '--mu-reference', '.75', '--temperatures', '0', '300']

    def wave_args(self, name='run', extra=()):
        return cli.parser().parse_args(['wavecar-hall', '--binary', str(self.binary),
            *self.source_arguments(), *self.scan_arguments(), '--output-dir', str(self.root/name), *extra])

    def native_process(self, *, returncode=0, footer=True, omit_csv=False):
        def run(argv, *, cwd, stdout, stderr, check):
            self.assertFalse(check)
            self.assertEqual(argv[-4:], ['-kubo', '2', '-kubo_pairs', 'PAIRS.csv'])
            stdout.write('native diagnostic retained\n')
            stderr.write('STOP 0\n' if not footer else '')
            if not omit_csv:
                (Path(cwd)/'PAIRS.csv').write_text(self.native_csv(footer))
            return SimpleNamespace(returncode=returncode)
        return run

    def run_success(self, name='run', extra=()):
        args = self.wave_args(name, extra)
        with patch.object(workflow.subprocess, 'run', side_effect=self.native_process()) as native:
            meta = workflow.wavecar_hall_command(args)
        self.assertEqual(native.call_count, 1)
        return args, meta

    def test_all_in_one_npz_only_and_reusable_cache_rescan_match(self):
        args, meta = self.run_success(extra=['--formats', 'npz'])
        out = args.output_dir
        self.assertEqual(set(p.name for p in (out/'hall').iterdir()), {'conductivity.npz', 'conductivity.json'})
        self.assertEqual(meta['output_formats'], ['npz'])
        manifest = json.loads((out/'workflow.json').read_text())
        self.assertEqual((manifest['status'], manifest['complete'], manifest['native_exit_code']), ('PASS', True, 0))
        self.assertTrue((out/'pairs/pairs.npz').is_file())
        self.assertIn('native diagnostic', (out/'native/stdout.log').read_text())
        rescan = cli.parser().parse_args(['pair-hall', '--pairs-dir', str(out/'pairs'),
            *self.scan_arguments(), '--mu-chunk', '1', '--formats', 'npz',
            '--output-dir', str(self.root/'rescan')])
        with patch.object(workflow.subprocess, 'run', side_effect=AssertionError('cache rescan invoked native')):
            workflow.pair_hall_command(rescan)
        with np.load(out/'hall/conductivity.npz', allow_pickle=False) as first, \
             np.load(self.root/'rescan/conductivity.npz', allow_pickle=False) as second:
            self.assertEqual(first.files, second.files)
            for key in first.files:
                np.testing.assert_array_equal(first[key], second[key])
            np.testing.assert_array_equal(np.unique(first['mu_eV']), [-.2, 0., .2])
            self.assertNotIn(.75, first['mu_eV'])

    def test_composed_import_and_scan_equal_all_in_one(self):
        args, _ = self.run_success()
        imported = cli.parser().parse_args(['import-pairs', '--csv', str(args.output_dir/'native/PAIRS.csv'),
            *self.source_arguments(), '--output-dir', str(self.root/'imported')])
        workflow.import_pairs_command(imported)
        scan = cli.parser().parse_args(['pair-hall', '--pairs-dir', str(imported.output_dir),
            *self.scan_arguments(), '--output-dir', str(self.root/'composed')])
        workflow.pair_hall_command(scan)
        for suffix in ('csv', 'dat'):
            self.assertEqual((args.output_dir/f'hall/conductivity.{suffix}').read_bytes(),
                             (scan.output_dir/f'conductivity.{suffix}').read_bytes())

    def test_existing_output_is_untouched_and_never_runs_native(self):
        args = self.wave_args(); args.output_dir.mkdir()
        sentinel = args.output_dir/'old-result'; sentinel.write_text('retain')
        with patch.object(workflow.subprocess, 'run') as native, self.assertRaisesRegex(ValueError, 'exists'):
            workflow.wavecar_hall_command(args)
        native.assert_not_called()
        self.assertEqual(list(args.output_dir.iterdir()), [sentinel])
        self.assertEqual(sentinel.read_text(), 'retain')

    def test_native_failures_keep_failed_manifest_and_no_hall(self):
        cases = [('exit', dict(returncode=7), ValueError),
                 ('stop-zero', dict(footer=False), ValueError),
                 ('no-csv', dict(omit_csv=True), OSError)]
        for name, options, error in cases:
            with self.subTest(name=name):
                args = self.wave_args(name)
                with patch.object(workflow.subprocess, 'run', side_effect=self.native_process(**options)), \
                     self.assertRaises(error):
                    workflow.wavecar_hall_command(args)
                manifest = json.loads((args.output_dir/'workflow.json').read_text())
                self.assertEqual(manifest['status'], 'FAILED')
                self.assertIs(manifest['complete'], False)
                self.assertEqual(manifest['native_exit_code'], options.get('returncode', 0))
                self.assertTrue(manifest['error'])
                self.assertFalse((args.output_dir/'hall').exists())
                self.assertFalse((args.output_dir/'pairs').exists())

    def test_process_launch_error_and_late_kernel_error_keep_manifest(self):
        args = self.wave_args('exec-error')
        with patch.object(workflow.subprocess, 'run', side_effect=PermissionError('not executable')), \
             self.assertRaises(PermissionError):
            workflow.wavecar_hall_command(args)
        manifest = json.loads((args.output_dir/'workflow.json').read_text())
        self.assertEqual(manifest['status'], 'FAILED')
        self.assertNotIn('native_exit_code', manifest)
        args = self.wave_args('occupied-top', ['--mu-min', '3', '--mu-max', '4'])
        with patch.object(workflow.subprocess, 'run', side_effect=self.native_process()), \
             self.assertRaisesRegex(ValueError, 'highest source'):
            workflow.wavecar_hall_command(args)
        manifest = json.loads((args.output_dir/'workflow.json').read_text())
        self.assertEqual((manifest['status'], manifest['complete']), ('FAILED', False))
        self.assertTrue((args.output_dir/'pairs/pairs.json').exists())
        self.assertFalse((args.output_dir/'hall').exists())

    def test_mpi_setup_preflight_and_launcher_argv(self):
        args = self.wave_args('no-launcher', ['--mpi-procs', '2'])
        with patch.object(workflow.shutil, 'which', return_value=None), \
             patch.object(workflow.subprocess, 'run') as native, \
             self.assertRaisesRegex(ValueError, 'MPI launcher'):
            workflow.wavecar_hall_command(args)
        native.assert_not_called()
        self.assertFalse(args.output_dir.exists())
        args = self.wave_args('mpi', ['--mpi-procs', '2', '--mpi-launcher', 'chosen-mpiexec'])
        with patch.object(workflow.shutil, 'which', return_value='/fake/mpiexec') as which, \
             patch.object(workflow.subprocess, 'run', side_effect=self.native_process()) as native:
            workflow.wavecar_hall_command(args)
        which.assert_called_once_with('chosen-mpiexec')
        self.assertEqual(native.call_args.args[0][:3], ['/fake/mpiexec', '-np', '2'])

    def test_cli_dispatch_and_legacy_parser_defaults(self):
        parser = cli.parser()
        hall = parser.parse_args(['hall', '--curvature', 'cache', *self.scan_arguments(), '--output-dir', 'out'])
        self.assertEqual(hall.formats, ['csv', 'npz'])
        self.assertFalse(hall.band_resolved)
        self.assertEqual(self.wave_args().formats, ['csv', 'dat', 'npz'])
        argv = ['pair-hall', '--pairs-dir', 'cache', *self.scan_arguments(), '--output-dir', str(self.root/'dispatch')]
        with patch.dict(cli.PAIR_COMMANDS, {'pair-hall': lambda args: {'schema': 'test-dispatch', 'version': 1}}), \
             contextlib.redirect_stdout(io.StringIO()) as stdout:
            cli.main(argv)
        self.assertEqual(json.loads(stdout.getvalue())['schema'], 'test-dispatch')
        bundle = parser.parse_args(['bundle-hall', '--csv', 'bundle.csv', '--occupied', '1',
            *self.source_arguments(), *self.scan_arguments(), '--output-dir', str(self.root/'bundle')])
        with self.assertRaisesRegex(ValueError, 'only T=0'):
            workflow.bundle_hall_command(bundle)


class HallPlotTests(unittest.TestCase):
    def setUp(self):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        self.plt = plt
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.addCleanup(plt.close, 'all')
        self.root = Path(self.tmp.name)
        rows = []
        for temperature in (0., 300.):
            for region, sign in [('total', 1.), ('valley', -1.)]:
                for mu in (.2, -.2, 0.):
                    # Include other band rows to ensure only total-subspace rows are plotted.
                    for band in (0, 1):
                        sigma = sign*(mu + temperature/1000.) + 100*band
                        rows.append(dict(mu_eV=mu, mu_minus_reference_eV=mu-.75,
                            temperature_K=temperature, region=region, band_id=band,
                            sigma_e2_over_h=sigma, delta_sigma_e2_over_h=sign*mu+100*band))
        self.rows = rows
        write_hall(self.root/'tables', rows, {'schema': 'vaspberry.hall-spectrum', 'version': 1},
                   formats=['csv', 'dat', 'npz'])

    def test_csv_dat_npz_numeric_tables_are_equivalent(self):
        tables = [plot_hall.read_table(self.root/f'tables/conductivity.{fmt}') for fmt in ('csv', 'dat', 'npz')]
        for other in tables[1:]:
            self.assertEqual(tables[0].keys(), other.keys())
            for key in tables[0]:
                np.testing.assert_array_equal(tables[0][key], other[key])

    def test_all_input_formats_plot_same_curves_and_export_all_figure_formats(self):
        from matplotlib.axes import Axes
        original = Axes.plot
        curves = []
        def capture(ax, x, y, *args, **kwargs):
            curves.append((np.array(x), np.array(y), kwargs['label']))
            return original(ax, x, y, *args, **kwargs)
        for fmt in ('csv', 'dat', 'npz'):
            out = self.root/f'plot-{fmt}'
            with patch.object(Axes, 'plot', capture):
                plot_hall.main([str(self.root/f'tables/conductivity.{fmt}'), '--output-dir', str(out),
                    '--regions', 'total', 'valley', '--quantity', 'delta-sigma',
                    '--energy-origin-eV', '-1', '--energy-label', 'mu - VBM (eV)'])
            self.assertEqual((out/'hall.png').read_bytes()[:8], b'\x89PNG\r\n\x1a\n')
            self.assertTrue((out/'hall.pdf').read_bytes().startswith(b'%PDF-'))
            self.assertIn('<svg', (out/'hall.svg').read_text())
            manifest = json.loads((out/'plot.json').read_text())
            self.assertEqual(set(manifest['output_sha256']), {'hall.png', 'hall.pdf', 'hall.svg'})
        self.assertEqual(len(curves), 12)
        for i in range(4):
            x, y, label = curves[i]
            np.testing.assert_allclose(x, [.8, 1., 1.2])
            np.testing.assert_allclose(y, [-.2, 0., .2] if i < 2 else [.2, 0., -.2])
            for offset in (4, 8):
                np.testing.assert_array_equal(x, curves[i+offset][0])
                np.testing.assert_array_equal(y, curves[i+offset][1])
                self.assertEqual(label, curves[i+offset][2])

    def test_invalid_plot_requests_do_not_create_output_or_overwrite(self):
        source = self.root/'tables/conductivity.npz'
        figures_before = self.plt.get_fignums()
        for name, extra in [('region', ['--regions', 'missing']),
                            ('duplicate-format', ['--formats', 'png', 'png']),
                            ('temperature', ['--temperatures', '10'])]:
            with self.subTest(name=name), contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as exc:
                plot_hall.main([str(source), '--output-dir', str(self.root/name), *extra])
            self.assertEqual(exc.exception.code, 2)
            self.assertFalse((self.root/name).exists())
            self.assertEqual(self.plt.get_fignums(), figures_before)
        out = self.root/'existing'; out.mkdir(); (out/'sentinel').write_text('retain')
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            plot_hall.main([str(source), '--output-dir', str(out)])
        self.assertEqual([x.name for x in out.iterdir()], ['sentinel'])

    def test_nonfinite_and_duplicate_curve_data_reject(self):
        bad = dict(self.rows[0], sigma_e2_over_h=float('nan'))
        write_hall(self.root/'nonfinite', [bad], {'schema': 'test'}, formats=['csv'])
        with self.assertRaisesRegex(ValueError, 'nonfinite'):
            plot_hall.read_table(self.root/'nonfinite/conductivity.csv')
        write_hall(self.root/'duplicates', [self.rows[0], self.rows[0]], {'schema': 'test'}, formats=['npz'])
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as exc:
            plot_hall.main([str(self.root/'duplicates/conductivity.npz'), '--output-dir', str(self.root/'bad-plot')])
        self.assertEqual(exc.exception.code, 2)
        self.assertFalse((self.root/'bad-plot').exists())
