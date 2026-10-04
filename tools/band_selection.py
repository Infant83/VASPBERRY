"""Shared standard-library parser for explicit selected source bands."""
import re

def _require(condition, message):
    if not condition:
        raise ValueError(message)


def parse_bands(value):
    """One-based IDs, inclusive ranges and comma lists; duplicates are errors."""
    _require(isinstance(value, str) and bool(value.strip()), 'bands require IDs/ranges, e.g. 31,33:34')
    result = []
    for token in value.split(','):
        token = token.strip()
        _require(re.fullmatch(r'[1-9][0-9]*(?::[1-9][0-9]*)?', token) is not None,
                'bands require positive IDs or inclusive ranges, e.g. 31,33:34')
        ends = [int(v) for v in token.split(':')]
        lo, hi = ends[0], ends[-1]
        _require(hi >= lo and hi-lo < 1000000, 'bands range must be increasing and bounded')
        result.extend(range(lo, hi+1))
    _require(len(set(result)) == len(result), 'duplicate selected bands')
    return sorted(result)

