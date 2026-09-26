"""Strict, standard-library reader for commented postprocessing study inputs.

Lists are whitespace-separated. Paths are relative to the input file; only
``~`` is expanded. Full-line and whitespace-prefixed ``#``/``;`` comments are
accepted. Names and keys are case-sensitive. No shell or expression is run.
"""
from __future__ import annotations

import configparser
from fractions import Fraction
import math
from pathlib import Path
import re


_NAME = re.compile(r'[A-Za-z][A-Za-z0-9_.-]*\Z')
_SPINS = {
    'soc': (1, 2, 1, 1),
    'scalar-degenerate': (1, 1, 2, 1),
    'collinear-up': (1, 1, 1, 2),
    'collinear-down': (2, 1, 1, 2),
}
_CHARACTER_KEYS = {'character_group', 'map_band', 'character_temperature',
                   'character_region', 'character_delta'}


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _keys(section, allowed, required=()):
    unknown = set(section) - set(allowed)
    _require(not unknown, f'[{section.name}] unknown key(s): '+', '.join(sorted(unknown)))
    missing = set(required) - set(section)
    _require(not missing, f'[{section.name}] missing key(s): '+', '.join(sorted(missing)))
    for key, value in section.items():
        _require(bool(value.strip()), f'[{section.name}] {key} must not be empty')


def _name(value, context):
    _require(bool(_NAME.fullmatch(value)),
             f'{context}: use a name beginning with a letter, then letters, digits, _, . or -')
    return value


def _float(value, context, *, fraction=False):
    try:
        result = float(Fraction(value)) if fraction else float(value)
    except (ValueError, ZeroDivisionError, OverflowError):
        raise ValueError(f'{context}: finite number required') from None
    _require(math.isfinite(result), f'{context}: finite number required')
    return result


def _int(value, context, minimum=1):
    _require(bool(re.fullmatch(r'[+-]?\d+', value)), f'{context}: integer required')
    result = int(value)
    _require(result >= minimum, f'{context}: integer >= {minimum} required')
    return result


def _list(value, convert, context, *, length=None, unique=False):
    result = [convert(token, context) for token in value.split()]
    _require(bool(result), f'{context}: nonempty list required')
    if length is not None:
        _require(len(result) == length, f'{context}: exactly {length} values required')
    if unique:
        _require(len(set(result)) == len(result), f'{context}: duplicate values')
    return result


def _temperatures(value, context):
    values = _list(value, _float, context, unique=True)
    _require(all(t >= 0 for t in values), f'{context}: temperatures must be >= 0 K')
    return values


def _scan(value, context, *, projection=False):
    tokens = value.split()
    _require(len(tokens) == 3, f'{context}: use MIN MAX N (energies in eV)')
    lo, hi = (_float(v, context) for v in tokens[:2])
    count = _int(tokens[2], context, minimum=2 if projection else 1)
    _require((count == 1 and lo == hi) or (count > 1 and lo < hi),
             f'{context}: use increasing MIN MAX, or equal energies with N=1 for total Hall')
    return dict(mu_min=lo, mu_max=hi, mu_num=count)


def _path(value, base):
    path = Path(value).expanduser()
    return str((base/path).resolve())


def _boolean(value, context):
    _require(value.lower() in ('true', 'false'), f'{context}: use true or false')
    return value.lower() == 'true'


def load_settings(path):
    """Return validated, JSON-serializable settings without probing data files.

    ``groups`` and ``regions`` already use the existing JSON schemas.
    ``differences`` is a list of [name, left, right] triples. ``projection`` is
    None when omitted. All data/output paths are absolute; a bare MPI launcher
    remains a command name for PATH resolution. Existence and source-data checks
    belong to the caller, so cache-only and plot-only use remain possible.
    """
    path = Path(path).expanduser().resolve()
    text = path.read_text(encoding='utf-8')
    _require(not re.search(r'^\s*\[DEFAULT\]', text, re.MULTILINE),
             '[DEFAULT] is not supported; put each setting in its named section')
    parser = configparser.ConfigParser(interpolation=None, strict=True,
        delimiters=('=',), inline_comment_prefixes=('#', ';'), empty_lines_in_values=False)
    parser.optionxform = str
    try:
        parser.read_string(text, source=str(path))
    except configparser.Error as exc:
        raise ValueError(f'invalid study input: {exc}') from None
    known = {'run', 'hall', 'projection', 'differences', 'plot'}
    for section in parser.sections():
        _require(section in known or section.startswith(('group ', 'region ')),
                 f'unknown section [{section}]')
    _require('run' in parser and 'hall' in parser, '[run] and [hall] are required')
    base = path.parent
    src = parser['run']
    _keys(src, {'wavecar', 'binary', 'output', 'mesh', 'spin_mode', 'energy_reference',
                'mpi_procs', 'mpi_launcher', 'plane_axes'},
          {'wavecar', 'binary', 'output', 'mesh', 'spin_mode', 'energy_reference'})
    run = {key: _path(src[key], base) for key in ('wavecar', 'binary', 'output')}
    run['mesh'] = _list(src['mesh'], lambda v, c: _int(v, c, 2), '[run] mesh', length=2)
    axes = _list(src.get('plane_axes', '0 1'), lambda v, c: _int(v, c, 0),
                 '[run] plane_axes', length=2, unique=True)
    _require(all(i <= 2 for i in axes), '[run] plane_axes: use distinct axes from 0 1 2')
    run['plane_axes'] = axes
    run['spin_mode'] = mode = src['spin_mode']
    _require(mode in _SPINS, '[run] spin_mode: choose '+', '.join(_SPINS))
    run.update(zip(('spin', 'spinor_components', 'spin_multiplicity', 'expected_nspin'), _SPINS[mode]))
    run['energy_reference'] = src['energy_reference']
    run['mpi_procs'] = _int(src.get('mpi_procs', '1'), '[run] mpi_procs')
    launcher = src.get('mpi_launcher', 'mpiexec')
    _require('/' in launcher or launcher.startswith('~') or len(launcher.split()) == 1,
             '[run] mpi_launcher: use one executable name or path, without launcher options')
    run['mpi_launcher'] = _path(launcher, base) if '/' in launcher or launcher.startswith('~') else launcher

    src = parser['hall']
    _keys(src, {'mu', 'reference', 'temperatures', 'pair_band_max'}, {'mu', 'reference', 'temperatures'})
    hall = _scan(src['mu'], '[hall] mu')
    hall.update(mu_reference=_float(src['reference'], '[hall] reference'),
                temperatures=_temperatures(src['temperatures'], '[hall] temperatures'),
                pair_band_max=_int(src['pair_band_max'], '[hall] pair_band_max', 2)
                if 'pair_band_max' in src else None)

    groups = []
    regions = []
    for section in parser.sections():
        src = parser[section]
        if section.startswith('group '):
            name = _name(section[6:], f'[{section}]')
            _keys(src, {'ions', 'orbitals'}, {'ions'})
            group = dict(name=name, ions=_list(src['ions'], _int, f'[{section}] ions', unique=True))
            if 'orbitals' in src:
                group['orbitals'] = _list(src['orbitals'], lambda v, c: v,
                                          f'[{section}] orbitals', unique=True)
            groups.append(group)
        elif section.startswith('region '):
            name = _name(section[7:], f'[{section}]')
            _require(name not in ('total', 'rest'), f'[{section}] reserved region name')
            _keys(src, {'center', 'radius', 'k_ids'})
            if 'k_ids' in src:
                _require(set(src) == {'k_ids'}, f'[{section}] use k_ids OR center and radius')
                regions.append(dict(name=name, k_ids=_list(src['k_ids'], _int,
                                                          f'[{section}] k_ids', unique=True)))
            else:
                _require(set(src) == {'center', 'radius'}, f'[{section}] center and radius required')
                center = _list(src['center'], lambda v, c: _float(v, c, fraction=True),
                               f'[{section}] center', length=3)
                radius = _float(src['radius'], f'[{section}] radius')
                _require(radius > 0, f'[{section}] radius must be > 0 inverse angstrom')
                regions.append(dict(name=name, center_fractional=center, radius_inv_A=radius))

    projection = None
    if 'projection' in parser:
        src = parser['projection']
        _keys(src, {'bands', 'axis', 'procar', 'outcar', 'mu', 'reference', 'temperatures'}, {'bands', 'axis'})
        _require(mode == 'soc', '[projection] requires [run] spin_mode = soc')
        _require(bool(groups), '[projection] requires at least one [group NAME]')
        projection = dict(bands=_list(src['bands'], _int, '[projection] bands', unique=True),
                          axis=_list(src['axis'], _float, '[projection] axis', length=3))
        _require(abs(math.hypot(*projection['axis'])-1) <= 1e-10,
                 '[projection] axis must be a Cartesian unit vector')
        for key in ('procar', 'outcar'):
            projection[key] = _path(src[key], base) if key in src else str(Path(run['wavecar']).parent/key.upper())
        projection.update(_scan(src['mu'], '[projection] mu', projection=True) if 'mu' in src else
                          {key: hall[key] for key in ('mu_min', 'mu_max', 'mu_num')})
        _require(projection['mu_num'] >= 2 and projection['mu_min'] < projection['mu_max'],
                 '[projection] mu requires increasing energies and at least two samples')
        projection['mu_reference'] = _float(src['reference'], '[projection] reference') if 'reference' in src else hall['mu_reference']
        projection['temperatures'] = _temperatures(src['temperatures'], '[projection] temperatures') if 'temperatures' in src else list(hall['temperatures'])
    else:
        _require(not groups, '[group NAME] requires a [projection] section')

    region_names = {'total', 'rest'} | {r['name'] for r in regions}
    differences = []
    if 'differences' in parser:
        src = parser['differences']
        for name, value in src.items():
            _name(name, '[differences]')
            _require(name not in region_names, f'[differences] {name}: duplicate/reserved region name')
            names = value.split()
            _require(len(names) == 2 and len(set(names)) == 2 and set(names) <= region_names,
                     f'[differences] {name}: use two distinct existing region names')
            differences.append([name, *names])
            region_names.add(name)

    src = parser['plot'] if 'plot' in parser else None
    values = dict(src) if src is not None else {}
    if src is not None:
        _keys(src, {'hall_regions', 'hall_temperatures', 'hall_quantity'} | _CHARACTER_KEYS)
    plot = dict(hall_regions=_list(values.get('hall_regions', 'total'), lambda v, c: v,
                                  '[plot] hall_regions', unique=True),
                hall_temperatures=_temperatures(values['hall_temperatures'], '[plot] hall_temperatures')
                if 'hall_temperatures' in values else list(hall['temperatures']),
                hall_quantity=values.get('hall_quantity', 'sigma'))
    _require(set(plot['hall_regions']) <= region_names, '[plot] hall_regions contains an unknown region')
    _require(set(plot['hall_temperatures']) <= set(hall['temperatures']),
             '[plot] hall_temperatures must be saved [hall] temperatures')
    _require(plot['hall_quantity'] in ('sigma', 'delta-sigma'), '[plot] hall_quantity: use sigma or delta-sigma')
    if projection is None:
        _require(not set(values) & _CHARACTER_KEYS, '[plot] character settings require [projection]')
    else:
        plot.update(character_group=values.get('character_group', groups[0]['name']),
                    map_band=_int(values['map_band'], '[plot] map_band') if 'map_band' in values else projection['bands'][0],
                    character_temperature=_float(values['character_temperature'], '[plot] character_temperature')
                    if 'character_temperature' in values else projection['temperatures'][0],
                    character_region=values.get('character_region', 'total'),
                    character_delta=_boolean(values.get('character_delta', 'false'), '[plot] character_delta'))
        _require(plot['character_group'] in {g['name'] for g in groups}, '[plot] character_group is unknown')
        _require(plot['character_temperature'] in projection['temperatures'],
                 '[plot] character_temperature must be a saved [projection] temperature')
        _require(plot['character_region'] in {'total', 'rest'} | {r['name'] for r in regions},
                 '[plot] character_region must be a region; difference channels are total-Hall only')
    return dict(schema='vaspberry.postprocess-settings', version=1, config_path=str(path),
                run=run, hall=hall, projection=projection, groups={'groups': groups},
                regions={'regions': regions}, differences=differences, plot=plot)
