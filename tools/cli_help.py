"""Shared plain-text CLI presentation; parser values and execution are unchanged."""
from __future__ import annotations

import argparse
import textwrap


class ReadableHelpFormatter(argparse.RawDescriptionHelpFormatter):
    """Use predictable wrapping and show useful defaults next to their options."""

    def __init__(self, prog):
        super().__init__(prog, width=80, max_help_position=25)

    def _get_help_string(self, action):
        help_text = action.help or ''
        if action.required and action.option_strings:
            help_text += ' (required)'
        value = action.default
        if (not action.required and value is not None and value is not argparse.SUPPRESS
                and not isinstance(value, bool) and value != []
                and 'default' not in help_text.lower() and '%(default)' not in help_text):
            if isinstance(value, (list, tuple)):
                value = ' '.join(str(item) for item in value)
            help_text += ' (default: '+str(value)+')'
        return help_text

    def _format_usage(self, usage, actions, groups, prefix):
        formatted = super()._format_usage(usage, actions, groups, prefix)
        if usage is None and len(formatted.splitlines()) > 5:
            # Long workflows are easier to enter through their grouped options
            # than through a screen of nested optional brackets. Required
            # actions are explicitly marked by _get_help_string above.
            return (prefix or 'usage: ')+self._prog+' [OPTIONS]\n\n'
        return formatted

    def _fill_text(self, text, width, indent):
        # Keep paragraph/command boundaries in descriptions without wide prose.
        output = []
        for line in text.splitlines():
            if not line.strip():
                output.append('')
            elif line.startswith('  '):
                output.append(indent+line)
            else:
                output.append(textwrap.fill(line, width, initial_indent=indent,
                                            subsequent_indent=indent,
                                            break_long_words=False, break_on_hyphens=False))
        return '\n'.join(output)


OPTION_HELP = {
    '--output-dir': 'new result directory; refuses an existing directory',
    '--run-dir': 'completed same-run VASP optical input directory',
    '--mpi-launcher': 'MPI launcher from the same runtime as the native executable',
    '--occupied': 'leading occupied-band count, starting at source band 1',
    '--energy-reference': 'description of the unchanged source energy zero',
    '--band': 'one-based selected band index',
    '--matrices': 'NPZ operator matrix bundle; pair with --metadata',
    '--metadata': 'matching JSON operator metadata and provenance',
    '--n-bands': 'inclusive one-based FIRST:LAST target-band range',
    '--m-bands': 'inclusive one-based FIRST:LAST intermediate-band range',
    '--pairs-dir': 'cache with pairs.npz and pairs.json, from import-pairs or velocity-pairs',
    '--character-dir': 'character.npz/character.json cache from procar_character.py project',
    '--hall-dir': 'saved attributed-Hall tables from procar_character.py hall',
    '--operators': 'operator cache from wannier-import',
    '--hh': 'Wannier Hamiltonian HH_R file',
    '--aa': 'matching Wannier position/connection AA_R file',
    '--mesh': 'full periodic 2D mesh, NX NY; match the source sampling',
    '--nx': 'number of stored mesh points along the first reciprocal axis',
    '--ny': 'number of stored mesh points along the second reciprocal axis',
    '--spin': 'one-based source spin channel; ISPIN=2 channels are separate',
    '--spinor-components': 'declared source component count: 1 scalar, 2 spinor',
    '--spin-multiplicity': '1 for resolved spin/channel; 2 for scalar spin-degenerate states',
    '--mu-min': 'lower chemical-potential endpoint in eV, in the source energy zero',
    '--mu-max': 'upper chemical-potential endpoint in eV, in the source energy zero',
    '--mu-num': 'number of chemical-potential points; equal endpoints require 1',
    '--mu-points': 'number of chemical-potential scan points',
    '--mu-reference': 'reference chemical potential in eV for subtracted response',
    '--temperatures': 'nonnegative response temperatures in K',
    '--temperature': 'response temperature in K',
    '--regions': 'named region JSON; region labels are user-defined',
    '--difference': 'signed region difference NAME:LEFT:RIGHT; no implicit factor 1/2',
    '--formats': 'saved output formats; choose one or more',
    '--memory-limit-mib': 'working-memory limit in MiB',
    '--time-limit': 'execution time limit in seconds',
    '--batch-size': 'number of k points processed in one batch',
    '--k-chunk': 'number of k points processed in one batch',
    '--mu-chunk': 'chemical-potential batch size',
    '--gap-threshold': 'minimum accepted energy gap in eV',
    '--gap-threshold-eV': 'minimum accepted energy gap in eV',
    '--degeneracy-threshold-eV': 'energy separation threshold in eV',
    '--degeneracy-policy': 'handling of unresolved energy degeneracies',
    '--band-resolved': 'also save individual-band contributions',
    '--plane-axes': 'ordered zero-based reciprocal axes defining the oriented 2D plane',
    '--photon-min': 'lower photon-energy endpoint in eV',
    '--photon-max': 'upper photon-energy endpoint in eV',
    '--photon-num': 'number of photon-energy points',
    '--initial': 'inclusive one-based initial-band endpoints, FIRST LAST',
    '--beam-vector': 'Cartesian propagation direction, X Y Z',
    '--polarization-axis': 'Cartesian direction defining the transverse polarization frame',
    '--relative-intensity-floor': 'dimensionless floor for unresolved optical selectivity',
    '--points-per-segment': 'number of sampled k points per path segment',
    '--vertices': 'fractional path vertices, K1 K2 K3 for each vertex',
    '--labels': 'labels of successive path vertices',
    '--workers': 'number of worker processes',
    '--refine': 'integer local mesh refinement factor',
    '--refine-radius': 'refinement disk radius in inverse Angstrom',
    '--refine-center': 'fractional in-plane refinement center; repeat for multiple centers',
    '--length': 'number of cells in the finite ribbon direction',
    '--direction': 'lattice direction of the ribbon',
    '--edge-width': 'number of boundary cells included in the edge weight',
    '--edge-cells': 'number of boundary cells included in the edge weight',
    '--periodic-axis': 'zero-based periodic lattice axis (0, 1 or 2)',
    '--open-axis': 'zero-based open lattice axis (0, 1 or 2)',
    '--q-min': 'lower fractional reciprocal-coordinate endpoint of the ribbon path',
    '--q-max': 'upper fractional reciprocal-coordinate endpoint of the ribbon path',
    '--edge-threshold': 'minimum edge weight used to label an edge state',
    '--kpoints': 'number of k points along the periodic ribbon direction',
    '--run-dirs': 'same-density instrumented VASP k-chunk run directories',
    '--augmentation': 'audited PAW spin-augmentation NPZ file',
    '--augmentation-metadata': 'matching PAW augmentation JSON provenance',
    '--poscar': 'matching POSCAR defining the lattice',
    '--eigenval': 'matching EIGENVAL supplying band energies',
    '--k-center': 'fractional in-plane coordinates K1 K2 of the first valley center',
    '--kp-center': 'fractional in-plane coordinates K1 K2 of the second valley center',
    '--plot': 'output figure filename; existing file policy is stated below',
    '--groups': 'named atom/orbital-group JSON; group labels are user-defined',
    '--procar': 'same-run PROCAR supplying atom/orbital/spin character',
    '--wavecar': 'same-run WAVECAR wavefunction and energy source',
    '--outcar': 'matching OUTCAR source metadata and spin frame',
    '--title': 'figure title',
    '--dpi': 'raster resolution in dots per inch',
    '--samples': 'number of plotted interpolation samples',
    '--energy-band': 'one-based energy-reference band; example-oriented default 19',
    '--valley-k': 'fractional coordinates of one named valley center',
    '--valley-kp': 'fractional coordinates of the other valley center',
    '--min-link-sv': 'minimum accepted link singular value (dimensionless)',
    '--min-neighbor-gap-ev': 'minimum neighbor-band separation in eV',
    '--max-abs-phi': 'maximum accepted absolute plaquette flux in radians',
    '--mass': 'analytic QWZ model mass parameter',
}


GROUPS = (
    ('Inputs and selection', {'--input', '--input-dir', '--run-dir', '--run-dirs',
      '--wavecar', '--waveder', '--incar', '--outcar', '--csv', '--pairs-dir',
      '--curvature', '--character-dir', '--hall-dir', '--matrices', '--metadata',
      '--operators', '--hh', '--aa', '--poscar', '--procar', '--groups', '--bands',
      '--occupied', '--n-bands', '--m-bands', '--kubo-source', '--mesh', '--nx',
      '--ny', '--plane-axes', '--spin', '--spinor-components', '--spin-multiplicity',
      '--energy-reference', '--normalization', '--augmentation', '--augmentation-metadata',
      '--axis', '--spin-basis'}),
    ('Response scan and regions', {'--mu-min', '--mu-max', '--mu-num', '--mu-points',
      '--mu-reference', '--mu-eV', '--temperatures', '--temperature', '--regions',
      '--difference', '--band-resolved'}),
    ('WAVECAR execution', {'--binary', '--mpi-procs', '--mpi-launcher'}),
    ('Output', {'--output-dir', '--output', '--plot', '--summary', '--formats'}),
)


def configure_help(parser, *, descriptions=None, epilogs=None, option_help=None):
    """Annotate and group existing actions without modifying their parse contract."""
    descriptions = descriptions or {}
    epilogs = epilogs or {}
    option_help = option_help or {}

    def visit(current, command=None, command_help=None):
        current.formatter_class = ReadableHelpFormatter
        if command in descriptions:
            current.description = descriptions[command]
        elif not current.description and command_help:
            current.description = command_help
        if command in epilogs:
            current.epilog = epilogs[command]
        for action in current._actions:
            if isinstance(action, argparse._SubParsersAction):
                help_by_name = {choice.dest: choice.help for choice in action._choices_actions}
                action.metavar = 'COMMAND'
                for name, child in action.choices.items():
                    visit(child, name, help_by_name.get(name))
                continue
            for spelling in action.option_strings:
                override = option_help.get((command, spelling))
                if override is not None:
                    action.help = override
                elif not action.help and spelling in OPTION_HELP:
                    action.help = OPTION_HELP[spelling]
        if command is None or len(current._actions) < 10:
            return
        # Help groups contain the same action objects; parser._actions,
        # required flags, defaults and mutual-exclusion groups remain intact.
        for title, flags in GROUPS:
            selected = [action for action in current._actions
                        if flags.intersection(action.option_strings)]
            if not selected:
                continue
            for group in current._action_groups:
                group._group_actions[:] = [a for a in group._group_actions if a not in selected]
            current.add_argument_group(title)._group_actions.extend(selected)
        remaining = [action for action in current._optionals._group_actions
                     if action.dest not in ('help', 'version')]
        if remaining:
            current._optionals._group_actions[:] = [a for a in current._optionals._group_actions
                                                    if a not in remaining]
            current.add_argument_group('Other controls')._group_actions.extend(remaining)

    visit(parser)
    return parser
