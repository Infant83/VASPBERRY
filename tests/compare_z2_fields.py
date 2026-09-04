#!/usr/bin/env python3
"""Compare serial and MPI VASPBERRY v1.2 Z2 field CSV files."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


EXPECTED_COLUMNS = (
    "cell_id",
    "q1",
    "q2",
    "q3",
    "kx_A-1",
    "ky_A-1",
    "kz_A-1",
    "half_bz",
    "tr_partner",
    "berry_flux_rad",
    "berry_curvature_A2",
    "min_link_singular_value",
    "nfield_raw",
    "nfield_int",
    "nfield_integer_residual",
    "nfield_tr_pair_residual",
    "flux_tr_odd_residual_rad",
    "pair_nfield_int_sum",
    "plaquette_checks_pass",
)

DISCRETE_COLUMNS = (
    "cell_id",
    "half_bz",
    "tr_partner",
    "nfield_int",
    "pair_nfield_int_sum",
    "plaquette_checks_pass",
)

CONTINUOUS_COLUMNS = tuple(
    column for column in EXPECTED_COLUMNS if column not in DISCRETE_COLUMNS
)

CONTINUOUS_METADATA = (
    "threshold_partner_fractional",
    "threshold_min_link_singular",
    "threshold_flux_tr_odd_rad",
    "threshold_nfield_integer",
    "threshold_total_chern",
    "threshold_phase_margin_fraction_pi",
    "threshold_max_abs_flux_rad",
    "total_chern",
    "minimum_link_singular_value",
    "max_flux_tr_odd_residual_rad",
    "max_nfield_integer_residual",
    "max_nfield_pair_residual",
    "max_abs_flux_rad",
)

REQUIRED_EXACT_METADATA = {
    "schema": "VASPBERRY_Z2_FIELD",
    "schema_version": "2",
    "vaspberry_version": "1.2.0",
    "result_status": "PASS",
    "result_kind": "FUKUI_HATSUGAI_NFIELD_Z2",
    "reportable_invariant": "1",
    "z2_invariant": "1",
    "half_top_z2_parity": "1",
    "half_bottom_z2_parity": "1",
    "half_bz_parity_consistent": "1",
    "nkx": "12",
    "nky": "12",
    "band_min": "1",
    "band_max": "10",
    "band_rank": "10",
    "spinor_components": "2",
    "index_base": "1",
    "input_trs_assumed": "1",
    "input_trs_independently_verified": "0",
    "numerical_self_consistency_checks_pass": "1",
}

REQUIRED_METADATA = frozenset(
    {
        *REQUIRED_EXACT_METADATA,
        *CONTINUOUS_METADATA,
        "plaquette_orientation",
        "nfield_definition",
        "overlap_backend",
        "check_scope",
        "check_scope_excludes",
        "physical_tr_rule",
        "nfield_note",
        "threshold_policy",
        "half_top_nfield_sum",
        "half_bottom_nfield_sum",
    }
)

EXPECTED_CELL_IDS = frozenset(range(1, 12 * 12 + 1))
MAX_REPORTED_MISMATCHES = 20


class Z2ComparisonError(ValueError):
    """A malformed input or semantic mismatch between inputs."""


@dataclass(frozen=True)
class Z2Field:
    """Parsed, type-checked v1.2 Z2 field output."""

    path: Path
    metadata: dict[str, str]
    diagnostics: dict[str, float]
    rows: dict[int, dict[str, int | float]]


def _parse_integer(value: str, *, location: str) -> int:
    try:
        return int(value.strip())
    except ValueError as error:
        raise Z2ComparisonError(
            f"{location}: expected an integer, found {value!r}"
        ) from error


def _parse_finite_float(value: str, *, location: str) -> float:
    try:
        result = float(value.strip())
    except ValueError as error:
        raise Z2ComparisonError(
            f"{location}: expected a floating-point number, found {value!r}"
        ) from error
    if not math.isfinite(result):
        raise Z2ComparisonError(
            f"{location}: value must be finite, found {value!r}"
        )
    return result


def _read_records(path: Path) -> tuple[dict[str, str], list[str], list[list[str]]]:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as error:
        raise Z2ComparisonError(f"cannot read {path}: {error}") from error

    metadata: dict[str, str] = {}
    header: list[str] | None = None
    records: list[list[str]] = []

    for line_number, line in enumerate(lines, start=1):
        stripped = line.strip()
        if not stripped:
            continue
        if header is None and stripped.startswith("#"):
            body = stripped[1:].strip()
            if "=" not in body:
                raise Z2ComparisonError(
                    f"{path}:{line_number}: malformed metadata line"
                )
            key, value = (part.strip() for part in body.split("=", 1))
            if not key:
                raise Z2ComparisonError(
                    f"{path}:{line_number}: empty metadata key"
                )
            if key in metadata:
                raise Z2ComparisonError(
                    f"{path}:{line_number}: duplicate metadata key {key!r}"
                )
            metadata[key] = value
            continue

        try:
            fields = next(csv.reader([line], strict=True))
        except csv.Error as error:
            raise Z2ComparisonError(
                f"{path}:{line_number}: invalid CSV: {error}"
            ) from error
        fields = [field.strip() for field in fields]
        if header is None:
            header = fields
        else:
            records.append(fields)

    if header is None:
        raise Z2ComparisonError(f"{path}: missing CSV header")
    return metadata, header, records


def load_z2_field(path: Path) -> Z2Field:
    """Load and validate one reportable 12 x 12 v1.2 Z2 field."""
    metadata, header, records = _read_records(path)

    if len(header) != len(set(header)):
        duplicates = sorted({column for column in header if header.count(column) > 1})
        raise Z2ComparisonError(
            f"{path}: duplicate CSV column(s): {', '.join(duplicates)}"
        )
    missing_columns = sorted(set(EXPECTED_COLUMNS) - set(header))
    unexpected_columns = sorted(set(header) - set(EXPECTED_COLUMNS))
    if missing_columns or unexpected_columns:
        details = []
        if missing_columns:
            details.append(f"missing {', '.join(missing_columns)}")
        if unexpected_columns:
            details.append(f"unexpected {', '.join(unexpected_columns)}")
        raise Z2ComparisonError(
            f"{path}: schema-v2 CSV columns are invalid ({'; '.join(details)})"
        )

    for key, expected in REQUIRED_EXACT_METADATA.items():
        if key not in metadata:
            raise Z2ComparisonError(
                f"{path}: missing required metadata key {key!r}"
            )
        actual = metadata[key]
        if actual != expected:
            raise Z2ComparisonError(
                f"{path}: metadata {key!r} must be {expected!r}, found {actual!r}"
            )

    missing_metadata = sorted(REQUIRED_METADATA - set(metadata))
    if missing_metadata:
        raise Z2ComparisonError(
            f"{path}: missing schema-v2 metadata key(s): "
            f"{', '.join(missing_metadata)}"
        )

    diagnostics = {
        key: _parse_finite_float(
            metadata[key], location=f"{path}: metadata {key!r}"
        )
        for key in CONTINUOUS_METADATA
    }

    top_sum = _parse_integer(
        metadata["half_top_nfield_sum"],
        location=f"{path}: metadata 'half_top_nfield_sum'",
    )
    bottom_sum = _parse_integer(
        metadata["half_bottom_nfield_sum"],
        location=f"{path}: metadata 'half_bottom_nfield_sum'",
    )
    if abs(top_sum) % 2 != 1 or abs(bottom_sum) % 2 != 1:
        raise Z2ComparisonError(
            f"{path}: half-zone sums must both have odd parity for z2_invariant=1"
        )

    rows: dict[int, dict[str, int | float]] = {}
    for record_number, fields in enumerate(records, start=1):
        if len(fields) != len(header):
            raise Z2ComparisonError(
                f"{path}: data row {record_number} has {len(fields)} fields; "
                f"expected {len(header)}"
            )
        raw_row = dict(zip(header, fields))
        cell_id = _parse_integer(
            raw_row["cell_id"],
            location=f"{path}: data row {record_number} column 'cell_id'",
        )
        if cell_id in rows:
            raise Z2ComparisonError(f"{path}: duplicate cell_id {cell_id}")

        row: dict[str, int | float] = {}
        for column in DISCRETE_COLUMNS:
            row[column] = _parse_integer(
                raw_row[column],
                location=f"{path}: cell_id {cell_id} column {column!r}",
            )
        for column in CONTINUOUS_COLUMNS:
            row[column] = _parse_finite_float(
                raw_row[column],
                location=f"{path}: cell_id {cell_id} column {column!r}",
            )
        rows[cell_id] = row

    if set(rows) != EXPECTED_CELL_IDS:
        missing = sorted(EXPECTED_CELL_IDS - set(rows))
        unexpected = sorted(set(rows) - EXPECTED_CELL_IDS)
        details = []
        if missing:
            details.append(f"missing IDs {missing}")
        if unexpected:
            details.append(f"unexpected IDs {unexpected}")
        raise Z2ComparisonError(
            f"{path}: expected exactly 144 unique cell_id rows 1..144 "
            f"({'; '.join(details)})"
        )

    for cell_id, row in rows.items():
        if row["half_bz"] not in (-1, 1):
            raise Z2ComparisonError(
                f"{path}: cell_id {cell_id} has invalid half_bz {row['half_bz']}"
            )
        if row["tr_partner"] not in EXPECTED_CELL_IDS:
            raise Z2ComparisonError(
                f"{path}: cell_id {cell_id} has invalid tr_partner "
                f"{row['tr_partner']}"
            )
        if row["plaquette_checks_pass"] != 1:
            raise Z2ComparisonError(
                f"{path}: cell_id {cell_id} has plaquette_checks_pass="
                f"{row['plaquette_checks_pass']}, expected 1"
            )

    calculated_top_sum = sum(
        int(row["nfield_int"])
        for row in rows.values()
        if row["half_bz"] == 1
    )
    calculated_bottom_sum = sum(
        int(row["nfield_int"])
        for row in rows.values()
        if row["half_bz"] == -1
    )
    if calculated_top_sum != top_sum:
        raise Z2ComparisonError(
            f"{path}: half_top_nfield_sum={top_sum}, but row data sum to "
            f"{calculated_top_sum}"
        )
    if calculated_bottom_sum != bottom_sum:
        raise Z2ComparisonError(
            f"{path}: half_bottom_nfield_sum={bottom_sum}, but row data sum to "
            f"{calculated_bottom_sum}"
        )

    for cell_id, row in rows.items():
        partner_id = int(row["tr_partner"])
        partner_row = rows[partner_id]
        if partner_row["tr_partner"] != cell_id:
            raise Z2ComparisonError(
                f"{path}: TR partner map is not involutive at cell_id {cell_id} "
                f"(partner {partner_id} maps to {partner_row['tr_partner']})"
            )
        expected_pair_sum = int(row["nfield_int"]) + int(
            partner_row["nfield_int"]
        )
        if row["pair_nfield_int_sum"] != expected_pair_sum:
            raise Z2ComparisonError(
                f"{path}: cell_id {cell_id} has pair_nfield_int_sum="
                f"{row['pair_nfield_int_sum']}, expected {expected_pair_sum}"
            )

    return Z2Field(path=path, metadata=metadata, diagnostics=diagnostics, rows=rows)


def compare_z2_fields(
    serial: Z2Field,
    mpi: Z2Field,
    *,
    rtol: float,
    atol: float,
) -> None:
    """Raise Z2ComparisonError unless both parsed outputs agree semantically."""
    serial_keys = set(serial.metadata)
    mpi_keys = set(mpi.metadata)
    if serial_keys != mpi_keys:
        only_serial = sorted(serial_keys - mpi_keys)
        only_mpi = sorted(mpi_keys - serial_keys)
        details = []
        if only_serial:
            details.append(f"serial-only keys {only_serial}")
        if only_mpi:
            details.append(f"MPI-only keys {only_mpi}")
        raise Z2ComparisonError(
            "metadata key sets differ (" + "; ".join(details) + ")"
        )

    mismatches: list[str] = []

    for key in sorted(serial_keys - set(CONTINUOUS_METADATA)):
        serial_value = serial.metadata[key]
        mpi_value = mpi.metadata[key]
        if serial_value != mpi_value:
            mismatches.append(
                f"metadata {key!r}: serial={serial_value!r}, MPI={mpi_value!r}"
            )

    for key in CONTINUOUS_METADATA:
        serial_value = serial.diagnostics[key]
        mpi_value = mpi.diagnostics[key]
        if not math.isclose(serial_value, mpi_value, rel_tol=rtol, abs_tol=atol):
            mismatches.append(
                f"metadata {key!r}: serial={serial_value:.17g}, "
                f"MPI={mpi_value:.17g}"
            )

    for cell_id in sorted(EXPECTED_CELL_IDS):
        serial_row = serial.rows[cell_id]
        mpi_row = mpi.rows[cell_id]
        for column in DISCRETE_COLUMNS:
            serial_value = serial_row[column]
            mpi_value = mpi_row[column]
            if serial_value != mpi_value:
                mismatches.append(
                    f"cell_id {cell_id} column {column!r}: "
                    f"serial={serial_value}, MPI={mpi_value}"
                )
        for column in CONTINUOUS_COLUMNS:
            serial_value = float(serial_row[column])
            mpi_value = float(mpi_row[column])
            if not math.isclose(
                serial_value, mpi_value, rel_tol=rtol, abs_tol=atol
            ):
                mismatches.append(
                    f"cell_id {cell_id} column {column!r}: "
                    f"serial={serial_value:.17g}, MPI={mpi_value:.17g}"
                )

    if mismatches:
        shown = mismatches[:MAX_REPORTED_MISMATCHES]
        suffix = ""
        if len(mismatches) > len(shown):
            suffix = (
                f"\n... {len(mismatches) - len(shown)} additional mismatch(es)"
            )
        raise Z2ComparisonError(
            f"serial and MPI fields differ at {len(mismatches)} value(s):\n"
            + "\n".join(f"  - {message}" for message in shown)
            + suffix
        )


def _nonnegative_finite(value: str) -> float:
    try:
        result = float(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("must be a number") from error
    if not math.isfinite(result) or result < 0.0:
        raise argparse.ArgumentTypeError("must be finite and nonnegative")
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Semantically compare serial and MPI VASPBERRY v1.2 "
            "Z2_FIELD.csv outputs for the 12 x 12 Bi regression."
        )
    )
    parser.add_argument("serial_csv", type=Path, help="serial Z2_FIELD.csv")
    parser.add_argument("mpi_csv", type=Path, help="MPI Z2_FIELD.csv")
    parser.add_argument(
        "--rtol",
        type=_nonnegative_finite,
        default=1.0e-11,
        help="relative tolerance for finite continuous values (default: 1e-11)",
    )
    parser.add_argument(
        "--atol",
        type=_nonnegative_finite,
        default=1.0e-12,
        help="absolute tolerance for finite continuous values (default: 1e-12)",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        serial = load_z2_field(args.serial_csv)
        mpi = load_z2_field(args.mpi_csv)
        compare_z2_fields(serial, mpi, rtol=args.rtol, atol=args.atol)
    except Z2ComparisonError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1

    print(
        "PASS: serial and MPI Z2 fields are semantically equivalent "
        f"(144 cells, rtol={args.rtol:g}, atol={args.atol:g})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
