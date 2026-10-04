"""Independent directory and file selection preserve the optical source contract."""
import contextlib
import io
import json
import os
from pathlib import Path
import shutil
import sys
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
import test_waveder_hall as fixtures
import kubo_hall_workflow as workflow
from postprocess_config import load_settings
import vaspberry_kubo as cli
import vaspberry_post as post
import waveder_hall as wh


@contextlib.contextmanager
def changed_directory(path):
    previous = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


class KuboInputDirectoryTests(unittest.TestCase):
    def setUp(self):
        self.fixture = fixtures.WaveDerHallTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.root, self.run = self.fixture.root.resolve(), self.fixture.run.resolve()
        self.cwd = self.root/'invocation'; self.cwd.mkdir()

    def arguments(self, output='answer', extra=()):
        return ['--occupied', '2', '--mesh', '2', '2', '--energy-reference', 'unchanged fixture zero',
                '--mu-min', '-.2', '--mu-max', '.2', '--mu-num', '3', '--mu-reference', '.25',
                '--output-dir', str(self.root/output), *map(str, extra)]

    def calculate(self, output='answer', extra=(), command='kubo-hall'):
        args = cli.parser().parse_args([command, *self.arguments(output, extra)])
        return workflow.kubo_hall_command(args) if command == 'kubo-hall' else wh.command(args)

    def settings(self, extra='', output='ini-result'):
        path = self.root/'study.ini'
        path.write_text(f'''[run]
output = {output}
mesh = 2 2
spin_mode = soc
energy_reference = unchanged fixture zero
{extra}
[hall]
occupied = 2
mu = -0.2 0.2 3
reference = 0.25
temperatures = 0
''')
        return path

    def test_default_cli_uses_cwd_and_wavecar_override_does_not_redirect_auxiliaries(self):
        external = self.root/'wave-only'; external.mkdir()
        override = external/'renamed-wave'; shutil.copyfile(self.run/'WAVECAR', override)
        for name in ('WAVEDER', 'INCAR', 'OUTCAR'):
            (external/name).write_bytes(b'wrong auxiliary trap')
        with changed_directory(self.run):
            meta = self.calculate(extra=['--wavecar', '../wave-only/renamed-wave'])
            self.assertEqual(Path.cwd(), self.run)
        paths = meta['source_metadata']['source_paths']
        self.assertEqual(paths['WAVECAR'], str(override))
        for name in ('WAVEDER', 'INCAR', 'OUTCAR'):
            self.assertEqual(paths[name], str(self.run/name))

    def test_input_dir_works_from_another_cwd_and_relative_overrides_use_cwd(self):
        override = self.cwd/'custom-OUTCAR'; shutil.copyfile(self.run/'OUTCAR', override)
        with changed_directory(self.cwd):
            meta = self.calculate(extra=['--input-dir', '../run', '--outcar', 'custom-OUTCAR'])
            self.assertEqual(Path.cwd(), self.cwd)
        paths = meta['source_metadata']['source_paths']
        self.assertEqual(paths['OUTCAR'], str(override))
        self.assertEqual(paths['WAVEDER'], str(self.run/'WAVEDER'))

    def test_standalone_defaults_and_all_four_explicit_overrides(self):
        with changed_directory(self.run), contextlib.redirect_stdout(io.StringIO()):
            wh.main(self.arguments('standalone-default'))
        metadata = json.loads((self.root/'standalone-default/conductivity.json').read_text())
        self.assertEqual(metadata['source_metadata']['source_paths']['WAVEDER'], str(self.run/'WAVEDER'))
        overrides = []
        for name in wh.SOURCE_NAMES:
            copied = self.cwd/(name.lower()+'.custom'); shutil.copyfile(self.run/name, copied)
            overrides += ['--'+name.lower(), copied.name]
        with changed_directory(self.cwd):
            meta = self.calculate('all-overrides', ['--input-dir', '.', *overrides], 'waveder-hall')
        self.assertEqual(meta['source_metadata']['source_paths'],
            {name: str(self.cwd/(name.lower()+'.custom')) for name in wh.SOURCE_NAMES})

    def test_missing_cwd_auxiliary_is_reported_despite_wavecar_in_valid_run(self):
        with changed_directory(self.cwd), patch.object(workflow.subprocess, 'run') as native:
            with self.assertRaisesRegex(ValueError, str(self.cwd/'WAVEDER')):
                self.calculate(extra=['--wavecar', self.run/'WAVECAR'])
        native.assert_not_called()
        self.assertFalse((self.root/'answer').exists())

    def test_missing_input_dir_is_rejected_even_with_all_file_overrides(self):
        flags = ['--input-dir', self.root/'absent']
        for name in wh.SOURCE_NAMES:
            flags += ['--'+name.lower(), self.run/name]
        with self.assertRaisesRegex(ValueError, 'input directory does not exist'):
            self.calculate(extra=flags)
        self.assertFalse((self.root/'answer').exists())

    def test_explicit_file_override_retains_same_run_consistency_checks(self):
        bad = self.cwd/'OUTCAR'; bad.write_text((self.run/'OUTCAR').read_text())
        # Change the effective electron count, independently of fixture smearing.
        bad.write_text(bad.read_text().replace('NELECT = 2', 'NELECT = 3'))
        with self.assertRaisesRegex(ValueError, 'NELECT'):
            self.calculate(extra=['--input-dir', self.run, '--outcar', bad])
        self.assertFalse((self.root/'answer').exists())

    def test_legacy_single_run_accepts_file_override_but_multiple_runs_reject_it(self):
        override = self.cwd/'der'; shutil.copyfile(self.run/'WAVEDER', override)
        meta = self.calculate(extra=['--run-dir', self.run, '--waveder', override])
        self.assertEqual(meta['source_metadata']['source_paths']['WAVEDER'], str(override))
        with self.assertRaisesRegex(ValueError, 'multiple --run-dir'):
            self.calculate('rejected', ['--run-dir', self.run, self.cwd, '--waveder', override])
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            cli.parser().parse_args(['kubo-hall', *self.arguments('rejected',
                ['--input-dir', self.run, '--run-dir', self.run])])
        self.assertFalse((self.root/'rejected').exists())

    def test_ini_defaults_to_invocation_cwd_without_required_wavecar(self):
        path = self.settings()
        with changed_directory(self.run), contextlib.redirect_stdout(io.StringIO()):
            settings = load_settings(path)
            record = post.run_calculation(settings)
            self.assertEqual(Path.cwd(), self.run)
        self.assertEqual(settings['run']['input_dir'], str(self.run))
        for name in wh.SOURCE_NAMES:
            self.assertEqual(record['inputs'][name.lower()]['path'], str(self.run/name))

    def test_ini_explicit_paths_are_config_relative_and_override_only_one_file(self):
        override = self.root/'renamed-der'; shutil.copyfile(self.run/'WAVEDER', override)
        path = self.settings('input_dir = run\nwaveder = renamed-der')
        with changed_directory(self.cwd), contextlib.redirect_stdout(io.StringIO()):
            settings = load_settings(path); record = post.run_calculation(settings)
        self.assertEqual(record['inputs']['waveder']['path'], str(override))
        self.assertEqual(record['inputs']['wavecar']['path'], str(self.run/'WAVECAR'))

    def test_ini_wavecar_override_does_not_select_auxiliary_directory(self):
        path = self.settings('wavecar = run/WAVECAR')
        with changed_directory(self.cwd):
            settings = load_settings(path)
            self.assertEqual(settings['run']['waveder'], str(self.cwd/'WAVEDER'))
            with self.assertRaisesRegex(ValueError, str(self.cwd/'WAVEDER')):
                post.run_calculation(settings)

    def test_ini_alias_and_exclusive_directories(self):
        path = self.settings('optical_run_dir = run')
        settings = load_settings(path)
        self.assertEqual(settings['run']['input_dir'], str(self.run))
        path = self.settings('input_dir = run\noptical_run_dir = run')
        with self.assertRaisesRegex(ValueError, 'mutually exclusive'):
            load_settings(path)

    def test_ini_projection_defaults_stay_at_input_dir_when_wavecar_is_elsewhere(self):
        path = self.settings('input_dir = run\nwavecar = foreign/WAVECAR\nkubo_source = wavecar')
        path.write_text(path.read_text().replace('occupied = 2\n', '')+'''
[projection]
bands = 1
axis = 0 0 1
[group one]
ions = 1
''')
        settings = load_settings(path)
        self.assertEqual(settings['run']['wavecar'], str(self.root/'foreign/WAVECAR'))
        self.assertEqual(settings['projection']['procar'], str(self.run/'PROCAR'))
        self.assertEqual(settings['projection']['outcar'], str(self.run/'OUTCAR'))

    def test_ini_requires_input_directory_even_with_all_source_overrides(self):
        explicit = 'input_dir = missing\n'+'\n'.join(name.lower()+' = run/'+name for name in wh.SOURCE_NAMES)
        settings = load_settings(self.settings(explicit))
        with self.assertRaisesRegex(ValueError, 'input directory does not exist'):
            post.run_calculation(settings)
        self.assertFalse(Path(settings['run']['output']).exists())

    def test_wavecar_opt_in_uses_input_directory_default_and_records_resolved_path(self):
        # Stop immediately after subprocess launch; source selection happens
        # before any native data import. This tests command construction only.
        binary = self.cwd/'binary'; binary.write_bytes(b'not executed')
        args = cli.parser().parse_args(['kubo-hall', '--kubo-source', 'wavecar',
            '--input-dir', str(self.run), '--binary', str(binary), '--mesh', '2', '2',
            '--energy-reference', 'fixture', '--mu-min', '0', '--mu-max', '0', '--mu-num', '1',
            '--mu-reference', '0', '--output-dir', str(self.root/'wavecar-output')])
        with patch.object(workflow.subprocess, 'run', side_effect=RuntimeError('stop after selection')) as native, \
             contextlib.redirect_stderr(io.StringIO()), self.assertRaisesRegex(RuntimeError, 'stop after selection'):
            workflow.kubo_hall_command(args)
        argv = native.call_args.args[0]
        self.assertEqual(argv[argv.index('--wavecar')+1], str(self.run/'WAVECAR'))
        manifest = json.loads((args.output_dir/'workflow.json').read_text())
        self.assertEqual(manifest['source_paths'], {'WAVECAR': str(self.run/'WAVECAR')})
