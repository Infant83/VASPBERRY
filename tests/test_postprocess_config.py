"""Human study inputs preserve scientific choices and fail before execution."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from postprocess_config import load_settings


BASE = '''# Paths and scientific choices, once per study.
[run]
wavecar = input files/WAVECAR
binary = build/vaspberry
output = results/first
mesh = 14 14
spin_mode = soc
energy_reference = unchanged VASP eigenvalue zero

[hall]
mu = -3.15 -1.05 101
reference = -2.5319
temperatures = 0 100 300
'''
PROJECTION = '''
[projection]
bands = 31 32
axis = 0 0 1

[group Cr_lower_Cl_side]
ions = 1
orbitals = dz2 x2-y2 dxy

[group Cr_upper_S_side]
ions = 2
'''
REGIONS = '''
[region K]
center = 1/3 1/3 0 ; fractional reciprocal coordinates
radius = 0.35 # inverse angstrom

[region Kprime]
center = -1/3 -1/3 0
radius = 0.35

[differences]
valley = K Kprime
'''


class PostprocessConfigTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='vaspberry settings ')
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name).resolve()
        self.path = self.root/'study.ini'

    def load(self, text=BASE):
        self.path.write_text(text)
        return load_settings(self.path)

    def test_minimal_settings_are_explicit_serializable_and_do_not_touch_data(self):
        settings = self.load()
        self.assertEqual(settings['schema'], 'vaspberry.postprocess-settings')
        self.assertEqual(settings['version'], 1)
        self.assertEqual(settings['config_path'], str(self.path))
        run = settings['run']
        self.assertEqual(run['wavecar'], str(self.root/'input files/WAVECAR'))
        self.assertEqual(run['output'], str(self.root/'results/first'))
        self.assertEqual((run['spin'], run['spinor_components'], run['spin_multiplicity'], run['expected_nspin']), (1, 2, 1, 1))
        self.assertEqual((run['mpi_procs'], run['mpi_launcher'], run['plane_axes']), (1, 'mpiexec', [0, 1]))
        self.assertEqual(settings['hall'], dict(mu_min=-3.15, mu_max=-1.05, mu_num=101,
            mu_reference=-2.5319, temperatures=[0., 100., 300.], pair_band_max=None))
        self.assertIsNone(settings['projection'])
        self.assertEqual(settings['groups'], {'groups': []})
        self.assertEqual(settings['regions'], {'regions': []})
        self.assertEqual(settings['differences'], [])
        self.assertEqual(settings['plot'], dict(hall_regions=['total'],
            hall_temperatures=[0., 100., 300.], hall_quantity='sigma'))
        self.assertEqual(json.loads(json.dumps(settings, allow_nan=False)), settings)
        self.assertEqual(list(self.root.iterdir()), [self.path])

    def test_paths_are_config_relative_with_spaces_and_without_environment_interpolation(self):
        config = BASE.replace('input files/WAVECAR', '$SOURCE/input files/WAVECAR')
        config = config.replace('build/vaspberry', '~/bin/my vaspberry')
        config = config.replace('mesh = 14 14', 'mesh = 14 14\nmpi_launcher = ./MPI tools/mpiexec\nmpi_procs = 4')
        with patch.dict('os.environ', {'SOURCE': '/do/not/expand'}):
            settings = self.load(config)
        self.assertEqual(settings['run']['wavecar'], str(self.root/'$SOURCE/input files/WAVECAR'))
        self.assertEqual(settings['run']['binary'], str(Path('~/bin/my vaspberry').expanduser().resolve()))
        self.assertEqual(settings['run']['mpi_launcher'], str(self.root/'MPI tools/mpiexec'))
        self.assertEqual(settings['run']['mpi_procs'], 4)

    def test_all_spin_modes_have_distinct_physical_mappings(self):
        mappings = {'soc': (1, 2, 1, 1), 'scalar-degenerate': (1, 1, 2, 1),
                    'collinear-up': (1, 1, 1, 2), 'collinear-down': (2, 1, 1, 2)}
        for mode, expected in mappings.items():
            with self.subTest(mode=mode):
                run = self.load(BASE.replace('spin_mode = soc', 'spin_mode = '+mode))['run']
                self.assertEqual(tuple(run[k] for k in ('spin', 'spinor_components', 'spin_multiplicity', 'expected_nspin')), expected)

    def test_region_fraction_group_case_and_projection_inheritance(self):
        settings = self.load(BASE+PROJECTION+REGIONS)
        self.assertEqual(settings['regions']['regions'][0], dict(name='K',
            center_fractional=[1/3, 1/3, 0.], radius_inv_A=.35))
        self.assertEqual(settings['differences'], [['valley', 'K', 'Kprime']])
        self.assertEqual(settings['groups']['groups'][0], dict(name='Cr_lower_Cl_side',
            ions=[1], orbitals=['dz2', 'x2-y2', 'dxy']))
        proj = settings['projection']
        self.assertEqual(proj['procar'], str(self.root/'input files/PROCAR'))
        self.assertEqual(proj['outcar'], str(self.root/'input files/OUTCAR'))
        self.assertEqual(proj['bands'], [31, 32])
        self.assertEqual(proj['axis'], [0., 0., 1.])
        for key in ('mu_min', 'mu_max', 'mu_num', 'mu_reference', 'temperatures'):
            self.assertEqual(proj[key], settings['hall'][key])
        self.assertIsNot(proj['temperatures'], settings['hall']['temperatures'])
        self.assertEqual(settings['plot']['character_group'], 'Cr_lower_Cl_side')
        self.assertEqual(settings['plot']['map_band'], 31)
        self.assertEqual(settings['plot']['character_temperature'], 0.)
        self.assertFalse(settings['plot']['character_delta'])

    def test_projection_scan_and_plot_can_differ_without_inheriting_pair_cap(self):
        config = BASE.replace('[hall]', '[hall]\npair_band_max = 40')+PROJECTION.replace(
            'axis = 0 0 1', 'axis = 0 0 1\nmu = -2.1 -1.05 101\nreference = -2.5\n'
            'temperatures = 100 200\nprocar = projected/PROCAR\noutcar = projected/OUTCAR')+REGIONS+'''
[plot]
hall_regions = total K Kprime valley
hall_temperatures = 100
hall_quantity = delta-sigma
character_group = Cr_upper_S_side
map_band = 30
character_temperature = 200
character_region = K
character_delta = true
'''
        settings = self.load(config)
        self.assertEqual(settings['hall']['pair_band_max'], 40)
        self.assertNotIn('pair_band_max', settings['projection'])
        self.assertEqual(settings['projection']['mu_min'], -2.1)
        self.assertEqual(settings['projection']['mu_reference'], -2.5)
        self.assertEqual(settings['projection']['temperatures'], [100., 200.])
        self.assertEqual(settings['projection']['procar'], str(self.root/'projected/PROCAR'))
        self.assertEqual(settings['plot']['map_band'], 30) # Map is independent of Hall-selected bands.
        self.assertEqual(settings['projection']['bands'], [31, 32])
        self.assertEqual(settings['plot']['hall_regions'], ['total', 'K', 'Kprime', 'valley'])
        self.assertEqual(settings['plot']['hall_temperatures'], [100.])
        self.assertEqual(settings['plot']['character_temperature'], 200.)
        self.assertTrue(settings['plot']['character_delta'])

    def test_k_ids_and_ordered_difference_channels_use_existing_schema(self):
        settings = self.load(BASE+'''\n[region Pocket]\nk_ids = 1 3 7\n
[differences]
contrast = Pocket rest
second = contrast total
[plot]
hall_regions = second contrast
''')
        self.assertEqual(settings['regions'], {'regions': [dict(name='Pocket', k_ids=[1, 3, 7])]})
        self.assertEqual(settings['differences'], [['contrast', 'Pocket', 'rest'], ['second', 'contrast', 'total']])

    def test_one_point_total_scan_is_supported_but_projection_requires_increasing_grid(self):
        config = BASE.replace('mu = -3.15 -1.05 101', 'mu = -1 -1 1')
        self.assertEqual(self.load(config)['hall']['mu_num'], 1)
        with self.assertRaisesRegex(ValueError, 'projection.*mu'):
            self.load(config+PROJECTION)

    def test_unknown_missing_duplicate_and_empty_settings_fail(self):
        cases = [
            (BASE+'\n[unknown]\nx=1\n', 'unknown section'),
            (BASE.replace('mesh = 14 14', 'MESH = 14 14'), 'unknown key'),
            (BASE.replace('mesh = 14 14', ''), 'missing key'),
            (BASE.replace('mesh = 14 14', 'mesh = 14 14\nmesh = 18 18'), 'invalid study input'),
            (BASE+'\n[hall]\nmu = -1 1 3\n', 'invalid study input'),
            (BASE+'\n[DEFAULT]\n', 'DEFAULT'),
            (BASE.replace('spin_mode = soc', 'spin_mode ='), 'must not be empty'),
            (BASE.replace('[run]', '[other]'), 'unknown section'),
            (BASE.split('[hall]')[0], 'run.*hall.*required'),
            ('mesh = 14 14\n'+BASE, 'invalid study input'),
        ]
        for config, error in cases:
            with self.subTest(error=error, text=config):
                with self.assertRaisesRegex(ValueError, error): self.load(config)
                self.assertFalse((self.root/'results').exists())

    def test_invalid_scientific_values_fail_in_parser(self):
        replacements = [
            ('mesh = 14 14', 'mesh = 1 14'), ('mesh = 14 14', 'mesh = 2.0 14'),
            ('mesh = 14 14', 'mesh = 2 2 2'), ('spin_mode = soc', 'spin_mode = auto'),
            ('mu = -3.15 -1.05 101', 'mu = -1 -3 101'),
            ('mu = -3.15 -1.05 101', 'mu = -1 -1 101'),
            ('mu = -3.15 -1.05 101', 'mu = -1 1 1'),
            ('mu = -3.15 -1.05 101', 'mu = -1 1 3.0'),
            ('reference = -2.5319', 'reference = nan'),
            ('reference = -2.5319', 'reference = inf'),
            ('temperatures = 0 100 300', 'temperatures = -1 100'),
            ('temperatures = 0 100 300', 'temperatures = 100 100.0'),
            ('temperatures = 0 100 300', 'temperatures = 0 NaN'),
        ]
        for old, new in replacements:
            with self.subTest(new=new), self.assertRaises(ValueError): self.load(BASE.replace(old, new))
        for extra in ('plane_axes=0 0', 'plane_axes=0 3', 'mpi_procs=0', 'mpi_launcher=mpiexec -n 4'):
            with self.subTest(extra=extra), self.assertRaises(ValueError):
                self.load(BASE.replace('[run]', '[run]\n'+extra))
        with self.assertRaises(ValueError): self.load(BASE.replace('[hall]', '[hall]\npair_band_max=1'))

    def test_group_projection_scope_errors_fail_without_silently_changing_physics(self):
        cases = [
            BASE+'\n[group lower]\nions=1\n',
            BASE+'\n[projection]\nbands=1\naxis=0 0 1\n',
            BASE.replace('spin_mode = soc', 'spin_mode = scalar-degenerate')+PROJECTION,
            BASE+PROJECTION.replace('ions = 1', 'ions = 0'),
            BASE+PROJECTION.replace('ions = 1', 'ions = 1 1'),
            BASE+PROJECTION.replace('bands = 31 32', 'bands = 31 31'),
            BASE+PROJECTION.replace('axis = 0 0 1', 'axis = 0 0 2'),
            BASE+PROJECTION.replace('axis = 0 0 1', 'axis = 0 nan 1'),
            BASE+PROJECTION.replace('orbitals = dz2 x2-y2 dxy', 'orbitals = dz2 dz2'),
            BASE+PROJECTION.replace('Cr_lower_Cl_side', '$unweighted'),
            BASE+PROJECTION.replace('Cr_lower_Cl_side', 'lower plane'),
            BASE+PROJECTION.replace('[projection]', '[projection]\npair_band_max=40'),
        ]
        for config in cases:
            with self.subTest(text=config), self.assertRaises(ValueError): self.load(config)

    def test_invalid_regions_and_differences_fail(self):
        suffixes = [
            '[region total]\nk_ids=1', '[region rest]\nk_ids=1',
            '[region K]\ncenter=0 0 0', '[region K]\ncenter=0 0 0\nradius=0',
            '[region K]\ncenter=1/0 0 0\nradius=0.3',
            '[region K]\ncenter=nan 0 0\nradius=0.3',
            '[region K]\ncenter=0 0\nradius=0.3',
            '[region K]\ncenter=0 0 0\nradius=.3\nk_ids=1',
            '[region K]\nk_ids=1 1', '[region K]\nk_ids=0',
            '[differences]\nvalley=K Kprime', '[differences]\nrest=total rest',
            '[differences]\ncontrast=total total',
        ]
        for suffix in suffixes:
            with self.subTest(text=suffix), self.assertRaises(ValueError): self.load(BASE+'\n'+suffix+'\n')

    def test_plot_choices_must_exist_in_saved_scan(self):
        cases = [
            'hall_regions=missing', 'hall_regions=total total', 'hall_temperatures=200',
            'hall_quantity=conductivity', 'character_group=missing', 'map_band=0',
            'character_temperature=200', 'character_region=valley', 'character_delta=yes',
        ]
        for field in cases:
            with self.subTest(field=field), self.assertRaises(ValueError):
                self.load(BASE+PROJECTION+REGIONS+'\n[plot]\n'+field+'\n')
        with self.assertRaisesRegex(ValueError, 'require.*projection'):
            self.load(BASE+'\n[plot]\nmap_band=1\n')

    def test_parser_import_and_load_need_only_the_standard_library(self):
        self.path.write_text(BASE)
        code = ('import sys; sys.path.insert(0, sys.argv[1]); '
                'from postprocess_config import load_settings; '
                's=load_settings(sys.argv[2]); '
                'assert s["run"]["mesh"]==[14,14]; '
                'assert "numpy" not in sys.modules; assert "matplotlib" not in sys.modules')
        result = subprocess.run([sys.executable, '-I', '-S', '-c', code, str(ROOT/'tools'), str(self.path)],
                                cwd='/', text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == '__main__':
    unittest.main()
