"""Public example INIs remain valid when the CI harness relocates them."""
from pathlib import Path
import sys
import tempfile
import unittest

from run_simple_postprocess_example import adapt_example_config

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
from postprocess_config import load_settings


class SimplePostprocessExampleTests(unittest.TestCase):
    def test_copied_examples_resolve_actual_input_independently_of_ini_location(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            inputs = root / 'actual input'
            inputs.mkdir()
            wavecar = inputs / 'WAVECAR'
            wavecar.write_bytes(b'path-resolution fixture; no calculation')
            binary = root / 'native'
            binary.write_bytes(b'path-resolution fixture; no execution')
            relocated = root / 'unrelated' / 'results' / 'copied-ini'
            relocated.mkdir(parents=True)
            for name in ('bi', 'bi-rescan'):
                with self.subTest(example=name):
                    config = adapt_example_config(
                        ROOT / f'examples/features/simple-postprocess/{name}.ini',
                        wavecar=wavecar, binary=binary, output=relocated/name)
                    ini = relocated / (name + '.ini')
                    with ini.open('w') as stream:
                        config.write(stream)
                    run = load_settings(ini)['run']
                    self.assertEqual(Path(run['input_dir']), inputs.resolve())
                    self.assertEqual(Path(run['wavecar']), wavecar.resolve())
                    self.assertEqual(Path(run['binary']), binary.resolve())
                    self.assertEqual(Path(run['output']), (relocated/name).resolve())
                    self.assertEqual(run['kubo_source'], 'wavecar')
                    self.assertEqual(run['mpi_procs'], 1)
                    self.assertTrue(Path(run['input_dir']).is_dir())


if __name__ == '__main__':
    unittest.main()
