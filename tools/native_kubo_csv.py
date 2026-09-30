"""Validate native trace/band CSVs before selecting or plotting any rows.

New schemas certify whole-file completion. Older saved schemas have no source
row counts or terminal completion contract; their visible rows are validated
but their metadata explicitly reports legacy_unverified completeness.
"""
from __future__ import annotations

import csv
import math
from pathlib import Path

TRACE_SCHEMA = "VASPBERRY_BARE_MOMENTUM_KUBO_BUNDLE_V2"
LEGACY_TRACE_SCHEMA = "VASPBERRY_BARE_MOMENTUM_KUBO_BUNDLE_V1"
BAND_SCHEMA = "VASPBERRY_BARE_MOMENTUM_KUBO_V3"
LEGACY_BAND_SCHEMA = "VASPBERRY_BARE_MOMENTUM_KUBO_V2"
TRACE_COLUMNS = ["spin", "k_index", "kx_frac", "ky_frac", "kz_frac",
                 "omega_z_A2", "min_external_gap_eV"]
BAND_COLUMNS = ["spin", "k_index", "band", "kx_frac", "ky_frac", "kz_frac",
                "energy_eV", "omega_z_A2", "min_gap_eV"]


def metadata(lines):
    result = {}
    for line in lines:
        if line.startswith("#") and "=" in line:
            key, value = (part.strip() for part in line[1:].split("=", 1))
            if key in result and key != "result_status" and result[key] != value:
                raise ValueError("conflicting native metadata: " + key)
            result[key] = value
    return result


def require_terminal_pass(lines):
    nonempty = [line.strip() for line in lines if line.strip()]
    if not nonempty or nonempty[-1] != "# result_status=PASS":
        raise ValueError("native CSV requires terminal result_status=PASS")
    statuses = [line.split("=", 1)[1].strip() for line in nonempty
                if line.startswith("#") and "=" in line
                and line[1:].split("=", 1)[0].strip() == "result_status"]
    if statuses not in (["PASS"], ["INCOMPLETE", "PASS"]):
        raise ValueError("native result_status must certify completion only at the end")


def positive_integer(value, name):
    try:
        result = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError("native " + name + " must be a positive integer") from exc
    if result < 1:
        raise ValueError("native " + name + " must be a positive integer")
    return result


def read_curvature_csv(path, *, kind):
    """Return metadata and all finite rows, with exact current-schema coverage."""
    lines = Path(path).read_text().splitlines()
    meta = metadata(lines)
    trace = kind == "trace"
    if kind not in ("trace", "band"):
        raise ValueError("native curvature kind must be trace or band")
    current = TRACE_SCHEMA if trace else BAND_SCHEMA
    legacy = LEGACY_TRACE_SCHEMA if trace else LEGACY_BAND_SCHEMA
    if meta.get("schema") not in (current, legacy):
        raise ValueError("unsupported native " + kind + " CSV schema")
    strict = meta["schema"] == current
    if strict:
        require_terminal_pass(lines)
    elif trace and meta.get("result_status") != "PASS":
        raise ValueError("legacy trace requires result_status=PASS")
    if meta.get("normalization") != "STANDARD_MINUS_TWO_IM":
        raise ValueError("native CSV must use the current physical normalization")
    columns = TRACE_COLUMNS if trace else BAND_COLUMNS
    gap_name = "min_external_gap_eV" if trace else "min_gap_eV"
    # An empty external space has no gap to report. NA is a defined absence,
    # not a NaN escape hatch, and the finite-basis trace must be exactly zero.
    if trace:
        no_external = (meta.get("no_external_states") == "true"
                       and meta.get("zero_trace_scope") == "TRUNCATED_WAVECAR_BASIS"
                       and meta.get("band_min") == "1"
                       and meta.get("band_max") == meta.get("band_rank") == meta.get("source_nbands")
                       and meta.get("source_nbands") is not None)
    else:
        no_external = False
    reader = csv.DictReader(line for line in lines if line.strip() and not line.startswith("#"))
    if reader.fieldnames != columns:
        raise ValueError("native curvature column names/order disagree")
    rows = []
    for raw in reader:
        if set(raw) != set(columns) or any(value is None for value in raw.values()):
            raise ValueError("incomplete native curvature row")
        try:
            row = {name: (None if name == gap_name and raw[name].strip() == "NA" and no_external
                          else float(raw[name])) for name in columns}
        except (TypeError, ValueError) as exc:
            raise ValueError("native curvature rows must be numeric") from exc
        if not all(value is None or math.isfinite(value) for value in row.values()):
            raise ValueError("nonfinite native curvature data")
        if no_external and (row[gap_name] is not None or row["omega_z_A2"] != 0.):
            raise ValueError("empty external space requires gap NA and zero finite-basis curvature")
        for name in ("spin", "k_index") if trace else ("spin", "k_index", "band"):
            row[name] = positive_integer(raw[name], name)
        rows.append(row)
    if not rows:
        raise ValueError("empty native curvature CSV")
    if strict:
        nk, ns, count = [positive_integer(meta.get(name), name) for name in
                         ("source_nkpoints", "source_nspin", "expected_rows")]
        first, last, rank, nb = [positive_integer(meta.get(name), name) for name in
                                ("band_min", "band_max", "band_rank", "source_nbands")]
        if ns not in (1, 2) or last < first or rank != last-first+1 or last > nb:
            raise ValueError("invalid native spin count or selected band range")
        expected = nk * ns * (1 if trace else rank)
        if count != expected or len(rows) != expected:
            raise ValueError("native expected_rows or complete row count disagrees")
    else:
        nk = max(row["k_index"] for row in rows)
        ns = max(row["spin"] for row in rows)
        first = min(row["band"] for row in rows) if not trace else 1
        last = max(row["band"] for row in rows) if not trace else 1
        expected = nk * ns * (last-first+1)
        if ns not in (1, 2) or len(rows) != expected:
            raise ValueError("legacy curvature needs complete visible spin/k/band coverage")
    seen = set()
    for row in rows:
        key = (row["spin"], row["k_index"]) + (() if trace else (row["band"],))
        if (not 1 <= key[0] <= ns or not 1 <= key[1] <= nk
                or (not trace and not first <= key[2] <= last) or key in seen):
            raise ValueError("native curvature needs unique complete spin/k/band coverage")
        seen.add(key)
    # Count, bounds and uniqueness imply every declared index is covered.
    meta["completion_validation"] = "terminal_pass_and_full_coverage" if strict else "legacy_unverified"
    return meta, rows
