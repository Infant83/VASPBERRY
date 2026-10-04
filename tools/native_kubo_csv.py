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
PAW_TRACE_SCHEMA = "VASPBERRY_WAVEDER_KUBO_OCCUPIED_V1"
PAW_BUNDLE_SCHEMA = "VASPBERRY_WAVEDER_KUBO_BUNDLE_V1"
PAW_BAND_SCHEMA = "VASPBERRY_WAVEDER_KUBO_BAND_V1"
PAW_SCHEMAS = (PAW_TRACE_SCHEMA, PAW_BUNDLE_SCHEMA, PAW_BAND_SCHEMA)
WAVEDER_OPERATOR = "VASP_5.4.4_LONGITUDINAL_PAW_OPTICAL"
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


def selected_band_ids(meta):
    """Read the explicit, sorted native selection without filling holes."""
    value = meta.get("band_ids", "")
    ids = [positive_integer(item, "band_ids") for item in value.split(",")]
    nb = positive_integer(meta.get("source_nbands"), "source_nbands")
    first, last, rank = [positive_integer(meta.get(name), name) for name in
                         ("band_min", "band_max", "band_rank")]
    if ids != sorted(set(ids)) or ids[0] != first or ids[-1] != last or len(ids) != rank or last > nb:
        raise ValueError("native band_ids must match the exact sorted selected band set")
    return ids


def read_curvature_csv(path, *, kind):
    """Return metadata and all finite rows, with exact current-schema coverage."""
    lines = Path(path).read_text().splitlines()
    meta = metadata(lines)
    trace = kind == "trace"
    if kind not in ("trace", "band"):
        raise ValueError("native curvature kind must be trace or band")
    current = TRACE_SCHEMA if trace else BAND_SCHEMA
    legacy = LEGACY_TRACE_SCHEMA if trace else LEGACY_BAND_SCHEMA
    selected_optical = meta.get("schema") == (PAW_BUNDLE_SCHEMA if trace else PAW_BAND_SCHEMA)
    optical = (trace and meta.get("schema") == PAW_TRACE_SCHEMA) or selected_optical
    if meta.get("schema") not in (current, legacy) and not optical:
        raise ValueError("unsupported native " + kind + " CSV schema")
    strict = meta["schema"] == current or optical
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
    ids = selected_band_ids(meta) if selected_optical else None
    if trace and selected_optical:
        no_external = (meta.get("no_external_states") == "true"
                       and meta.get("zero_trace_scope") == "TRUNCATED_WAVEDER_BASIS"
                       and len(ids) == positive_integer(meta.get("source_nbands"), "source_nbands"))
        if (len(ids) == int(meta["source_nbands"])) != no_external or (
                not no_external and meta.get("no_external_states") != "false"):
            raise ValueError("WAVEDER empty-complement metadata disagrees with the selected bands")
    elif trace and not optical:
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
        if optical:
            nb, nd, components, mult = [positive_integer(meta.get(name), name) for name in
                ("source_nbands", "source_ndbands", "spinor_components", "physical_spin_multiplicity")]
            if (ns not in (1, 2) or components not in (1, 2) or (ns == 2 and components != 1)
                    or mult != (2 if ns == 1 and components == 1 else 1) or nd > nb):
                raise ValueError("invalid WAVEDER source dimensions or physical spin multiplicity")
            if (meta.get("kubo_source") != "WAVEDER" or meta.get("source_operator") != WAVEDER_OPERATOR
                    or meta.get("producer_cluster_threshold_eV") != "0.002"):
                raise ValueError("unsupported WAVEDER operator metadata")
            if not no_external and any(row[gap_name] is None or row[gap_name] <= .002 for row in rows):
                raise ValueError("WAVEDER required gap must exceed the producer 0.002 eV threshold")
            if selected_optical:
                expected_kind = "ISOLATED_SELECTED_BUNDLE_TRACE" if trace else "ISOLATED_SELECTED_BAND_CURVATURE"
                intermediate = ("EXTERNAL_TO_SELECTED_BUNDLE_WITHIN_SOURCE_NBANDS" if trace
                                else "ALL_SOURCE_BANDS_EXCEPT_SELF")
                if (meta.get("result_kind") != expected_kind or meta.get("occupation_weighting") != "NONE"
                        or meta.get("required_pair_coverage") != "PASS"
                        or meta.get("intermediate_bands") != intermediate):
                    raise ValueError("unsupported WAVEDER selected-band metadata")
                high = set(range(nd+1, nb+1))
                selected_high = high.intersection(ids)
                if (trace and selected_high and selected_high != high) or (
                        not trace and selected_high and len(high) > 1):
                    raise ValueError("WAVEDER selection requires missing high-high matrix pairs")
                first, last, rank = ids[0], ids[-1], len(ids)
            else:
                occupations = [positive_integer(meta.get("occupied_bands_spin_"+str(s)), "occupied bands")
                               for s in range(1, ns+1)]
                if any(value >= nb or value > nd for value in occupations):
                    raise ValueError("WAVEDER occupied-ket coverage or stored empty bands are incomplete")
                if (meta.get("result_kind") != "T0_INSULATING_OCCUPIED_BUNDLE_TRACE"
                        or meta.get("complete_occupied_bundle") != "true"):
                    raise ValueError("unsupported WAVEDER occupied-bundle operator metadata")
                # Different collinear channels may have different occupied ranks.
                first, last, rank = 1, max(occupations), max(occupations)
        else:
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
    allowed_bands = set(ids) if ids is not None else range(first, last+1)
    for row in rows:
        key = (row["spin"], row["k_index"]) + (() if trace else (row["band"],))
        if (not 1 <= key[0] <= ns or not 1 <= key[1] <= nk
                or (not trace and key[2] not in allowed_bands) or key in seen):
            raise ValueError("native curvature needs unique complete spin/k/band coverage")
        seen.add(key)
    # Count, bounds and uniqueness imply every declared index is covered.
    meta["completion_validation"] = "terminal_pass_and_full_coverage" if strict else "legacy_unverified"
    return meta, rows
