"""Check that help presentation retains parsing and exposes scientific contracts."""
import argparse
from pathlib import Path
import sys
import unittest

TOOLS = Path(__file__).resolve().parents[1] / 'tools'
sys.path.insert(0, str(TOOLS))

from cli_help import configure_help
import vaspberry_kubo
import vaspberry_post
import procar_character
import vaspberry_transport


class CliHelpTests(unittest.TestCase):
    def test_grouping_preserves_required_flags_mutex_defaults_and_types(self):
        parser = argparse.ArgumentParser()
        sub = parser.add_subparsers(dest='command', required=True)
        child = sub.add_parser('charge', help='charge analysis')
        choice = child.add_mutually_exclusive_group(required=True)
        choice.add_argument('--occupied', type=int)
        choice.add_argument('--bands')
        child.add_argument('--output-dir', type=Path, required=True)
        child.add_argument('--mesh', type=int, nargs=2, required=True)
        child.add_argument('--mu-min', type=float, required=True)
        child.add_argument('--mu-max', type=float, required=True)
        child.add_argument('--temperatures', type=float, nargs='+', default=[0.])
        child.add_argument('--spin', type=int, default=1)
        child.add_argument('--formats', nargs='+', default=['csv', 'npz'])
        argv = ['charge', '--occupied', '8', '--output-dir', 'new-result',
                '--mesh', '12', '12', '--mu-min', '-1', '--mu-max', '1']
        expected = vars(parser.parse_args(argv)).copy()
        actions = tuple(child._actions)
        configure_help(parser)
        self.assertEqual(tuple(child._actions), actions)
        self.assertEqual(vars(parser.parse_args(argv)), expected)
        for bad in (argv + ['--bands', '1:8'], argv[:3], argv[0:1] + argv[3:]):
            with self.assertRaises(SystemExit):
                parser.parse_args(bad)
        text = child.format_help()
        for value in ('Inputs and selection:', 'Response scan and regions:',
                      'Output:', 'default: 1', 'new result directory', '(required)'):
            self.assertIn(value, text)

    def test_direct_subcommand_help_keeps_source_and_output_contracts(self):
        parsers = [(vaspberry_kubo.parser(), 'kubo-hall',
                    ('exactly one of --occupied or --bands', 'computes in Python',
                     'conductivity.csv', 'new result directory')),
                   (vaspberry_kubo.parser(), 'spin-hall',
                    ('gapped 2D', 'full PAW', 'instrumented VASP')),
                   (procar_character.parser(), 'hall',
                    ('attribution', 'not a spin, layer or orbital-current response',
                     'character.npz', 'pairs.npz')),
                   (vaspberry_post.parser(), 'plot',
                    ('run.json', 'RESULT/figures', 'workflow.json')),
                   (vaspberry_transport.build_parser(), 'sigma',
                    ('whitespace DAT', 'removed before input validation', 'KUBO.csv'))]
        for parser, name, expected in parsers:
            action = next(a for a in parser._actions
                          if isinstance(a, argparse._SubParsersAction))
            child = action.choices[name]
            text = ' '.join(child.format_help().split())
            with self.subTest(command=name):
                for value in expected:
                    self.assertIn(value, text)

    def test_all_dispatcher_commands_have_standalone_descriptions(self):
        for factory in (vaspberry_kubo.parser, vaspberry_post.parser,
                        vaspberry_transport.build_parser, procar_character.parser):
            parser = factory()
            action = next(a for a in parser._actions
                          if isinstance(a, argparse._SubParsersAction))
            for name, child in action.choices.items():
                with self.subTest(parser=parser.prog, command=name):
                    self.assertTrue(child.description)
                    self.assertNotIn('\x1b', child.format_help())


if __name__ == '__main__':
    unittest.main()
