#!/usr/bin/env python3
"""Guarded Wilson-loop Z2 diagnostics for a spinor VASP WAVECAR.

The legacy Fortran ``-z2`` path writes a gauge-dependent lattice-field
diagnostic.  This program instead builds unitary Wilson loops for the occupied
subspace and reports a Z2 invariant only after independent gap, time-reversal,
link-quality, topology, and fixed-grid convergence checks have passed.

An INVALID calculation is still useful as a diagnostic: its largest-gap
crossing parity is emitted as ``candidate_z2``, while the reportable ``z2``
field is JSON null.  Hard input or I/O errors leave no planned output files.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable, Iterable, Sequence

import numpy as np

try:
    from wavecar_fukui import (
        Grid,
        Wavecar,
        infer_uniform_grid,
        json_safe,
        link_matrix,
        wrap_angle,
    )
except ImportError:  # pragma: no cover - permits package-style imports in tests
    from tools.wavecar_fukui import (  # type: ignore
        Grid,
        Wavecar,
        infer_uniform_grid,
        json_safe,
        link_matrix,
        wrap_angle,
    )


__version__ = "1.1.0"

DEFAULT_CSV = "z2_wilson_wcc.csv"
DEFAULT_JSON = "z2_diagnostics.json"
DEFAULT_PLOT = "z2_wilson_wcc.png"


@dataclass(frozen=True)
class Thresholds:
    """Conservative defaults for deciding whether a parity is reportable."""

    min_direct_gap_ev: float = 1.0e-5
    min_indirect_gap_ev: float = 1.0e-5
    max_energy_tr_error_ev: float = 1.0e-5
    max_trim_kramers_energy_split_ev: float = 1.0e-5
    max_projector_tr_residual: float = 1.0e-6
    min_tr_g_basis_coverage: float = 1.0 - 1.0e-12
    min_tr_raw_norm_coverage: float = 0.999999
    max_chern_residual: float = 1.0e-4
    max_flux_odd_residual_rad: float = 1.0e-5
    min_link_singular_value: float = 1.0e-6
    min_plane_wave_coverage: float = 0.90
    max_unitarity_residual: float = 1.0e-8
    max_wcc_tr_set_distance: float = 1.0e-5
    max_endpoint_kramers_split: float = 1.0e-5
    max_wcc_position_change: float = 0.01
    min_combined_gap_ratio: float = 0.30
    max_neighbor_move_ratio: float = 0.30


@dataclass(frozen=True)
class LinkField:
    """Polar-unitary occupied-subspace links on a full pump grid."""

    axis: int
    stride: int
    unitary: np.ndarray
    min_singular: np.ndarray
    max_singular: np.ndarray
    coverage_left: np.ndarray
    coverage_right: np.ndarray
    unitarity_residual: np.ndarray


@dataclass(frozen=True)
class WilsonSpectrum:
    loop_axis: int
    pump_axis: int
    pump: np.ndarray
    wcc: np.ndarray


def polar_unitary(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return the SVD polar factor and singular values of an overlap matrix."""

    matrix = np.asarray(matrix, dtype=np.complex128)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("occupied-subspace overlap must be a square matrix")
    if not np.all(np.isfinite(matrix)):
        raise ValueError("occupied-subspace overlap contains nonfinite values")
    left, singular, right_h = np.linalg.svd(matrix, full_matrices=False)
    return left @ right_h, singular


def circular_delta(left: np.ndarray | float, right: np.ndarray | float):
    """Signed shortest displacement ``left-right`` in cycles."""

    return (np.asarray(left) - np.asarray(right) + 0.5) % 1.0 - 0.5


def circular_set_distance(left: Sequence[float], right: Sequence[float]) -> float:
    """Bottleneck distance between two unordered point sets on a circle.

    For equal-cardinality one-dimensional circular sets, an optimal
    orientation-preserving matching is one of the cyclic shifts of the sorted
    sets.  Testing those shifts avoids attaching physical meaning to sorted
    Wilson-loop branch labels.
    """

    a = np.sort(np.mod(np.asarray(left, dtype=float), 1.0))
    b = np.sort(np.mod(np.asarray(right, dtype=float), 1.0))
    if a.ndim != 1 or b.ndim != 1 or a.size == 0 or a.size != b.size:
        raise ValueError("circular sets must be nonempty and have equal size")
    return min(
        float(np.max(np.abs(circular_delta(a, np.roll(b, shift)))))
        for shift in range(a.size)
    )


def largest_circular_gap(wcc: Sequence[float]) -> tuple[float, float, int]:
    """Return (center, size, lower-index) of a deterministic largest gap."""

    values = np.sort(np.mod(np.asarray(wcc, dtype=float), 1.0))
    if values.ndim != 1 or values.size == 0 or not np.all(np.isfinite(values)):
        raise ValueError("WCC values must be a nonempty finite one-dimensional set")
    extended = np.r_[values, values[0] + 1.0]
    gaps = np.diff(extended)
    index = int(np.argmax(gaps))
    return float((values[index] + 0.5 * gaps[index]) % 1.0), float(gaps[index]), index


def z2pack_largest_gap_parity(wcc_half: np.ndarray) -> tuple[int, np.ndarray, np.ndarray]:
    """Evaluate the largest-gap crossing parity without tracking WCC branches.

    The reference point is the center of the largest WCC gap on every pump
    line. Between neighboring lines, the Z2Pack 2.2.1 convention counts the
    right-line WCCs in the linear interval between the two normalized gap
    positions. This set-based construction remains well-defined when sorted
    branches exchange labels.
    """

    values = np.asarray(wcc_half, dtype=float)
    if values.ndim != 2 or values.shape[0] < 2 or values.shape[1] < 2:
        raise ValueError("half-BZ WCC data require at least two lines and two bands")
    centers = np.empty(values.shape[0])
    sizes = np.empty(values.shape[0])
    for index, row in enumerate(values):
        centers[index], sizes[index], _ = largest_circular_gap(row)

    parity = 0
    for index in range(values.shape[0] - 1):
        left = centers[index]
        right = centers[index + 1]
        lower, upper = sorted((left, right))
        if upper - lower <= np.finfo(float).eps:
            continue
        # Z2Pack's discrete invariant counts the WCCs on the right line
        # inside the linear interval between the two normalized gap positions.
        # It does not replace that interval by the shortest periodic arc.
        row = values[index + 1]
        crossings = int(np.count_nonzero((row > lower) & (row < upper)))
        parity ^= crossings & 1
    return int(parity), np.mod(centers, 1.0), sizes


def endpoint_kramers_split(wcc: Sequence[float]) -> float:
    """Return the best adjacent-pair Kramers splitting on a circular spectrum."""

    values = np.sort(np.mod(np.asarray(wcc, dtype=float), 1.0))
    if values.size == 0 or values.size % 2:
        raise ValueError("a Kramers spectrum must contain a positive even number of WCCs")
    best = math.inf
    for offset in (0, 1):
        rolled = np.roll(values, -offset)
        pair_distance = np.abs(circular_delta(rolled[0::2], rolled[1::2]))
        best = min(best, float(np.max(pair_distance)))
    return best


def _pump_partner_indices(pump: np.ndarray, tolerance: float = 1.0e-7) -> np.ndarray:
    result = np.empty(pump.size, dtype=int)
    for i, value in enumerate(pump):
        distances = np.abs(circular_delta(pump, -value))
        j = int(np.argmin(distances))
        if distances[j] > tolerance:
            raise ValueError("pump grid is not closed under time reversal")
        result[i] = j
    return result


def wcc_time_reversal_diagnostics(spectrum: WilsonSpectrum) -> tuple[float, float]:
    """Return unordered q/-q set distance and endpoint Kramers splitting."""

    partners = _pump_partner_indices(spectrum.pump)
    set_error = max(
        circular_set_distance(spectrum.wcc[i], spectrum.wcc[j])
        for i, j in enumerate(partners)
    )
    endpoint_errors: list[float] = []
    for target in (0.0, 0.5):
        distances = np.abs(circular_delta(spectrum.pump, target))
        index = int(np.argmin(distances))
        if distances[index] > 1.0e-7:
            raise ValueError("pump mesh must contain the q=0 and q=1/2 TR lines")
        endpoint_errors.append(endpoint_kramers_split(spectrum.wcc[index]))
    return float(set_error), float(max(endpoint_errors))


def half_bz_spectrum(spectrum: WilsonSpectrum) -> tuple[np.ndarray, np.ndarray]:
    """Select the closed pump interval q=0...1/2 in increasing order."""

    q = np.mod(spectrum.pump, 1.0)
    mask = (q >= -1.0e-10) & (q <= 0.5 + 1.0e-10)
    indices = np.flatnonzero(mask)
    order = indices[np.argsort(q[indices])]
    if order.size < 3 or abs(q[order[0]]) > 1.0e-7 or abs(q[order[-1]] - 0.5) > 1.0e-7:
        raise ValueError("pump mesh must resolve the closed half Brillouin zone")
    return q[order], spectrum.wcc[order]


def fixed_grid_convergence(
    full: WilsonSpectrum,
    stride_two: WilsonSpectrum,
) -> dict[str, float]:
    """Evaluate fixed-grid analogues of Z2Pack POS, GAP, and MOVE controls.

    MOVE uses Z2Pack's union-largest-gap cut before sorting neighboring WCC
    sets, divided by the smaller of their largest-gap sizes.  GAP checks both
    directions: neighbor WCCs must remain away from the other line's largest-
    gap center by a fraction of that gap size.  POS compares full loops with
    loops built from every second k point.  Unlike adaptive Z2Pack, this
    function does not insert pump lines, so these are fixed-grid guards rather
    than a convergence proof or an error bar.
    """

    if full.loop_axis != stride_two.loop_axis:
        raise ValueError("full and stride-two spectra use different loop axes")
    if full.pump.shape != stride_two.pump.shape or np.max(
        np.abs(circular_delta(full.pump, stride_two.pump))
    ) > 1.0e-10:
        raise ValueError("full and stride-two spectra use different pump lines")
    position = max(
        z2pack_max_move(left, right)
        for left, right in zip(full.wcc, stride_two.wcc)
    )
    _, half = half_bz_spectrum(full)
    move_ratios: list[float] = []
    gap_ratios: list[float] = []
    for left, right in zip(half[:-1], half[1:]):
        move_ratios.append(z2pack_move_ratio(left, right))
        gap_ratios.append(z2pack_gap_ratio(left, right))
    return {
        "full_vs_stride2_wcc_distance_cycles": float(position),
        "neighbor_move_ratio": float(max(move_ratios)),
        "neighbor_to_largest_gap_ratio": float(min(gap_ratios)),
    }


def z2pack_max_move(left: Sequence[float], right: Sequence[float]) -> float:
    """Z2Pack neighbor matcher: sort after cutting the union's largest gap."""

    left_values = np.mod(np.asarray(left, dtype=float), 1.0)
    right_values = np.mod(np.asarray(right, dtype=float), 1.0)
    if left_values.shape != right_values.shape or left_values.ndim != 1:
        raise ValueError("neighbor WCC sets must be one-dimensional and equal-sized")
    cut, _, _ = largest_circular_gap(np.r_[left_values, right_values])
    left_sorted = np.sort(np.mod(left_values - cut, 1.0))
    right_sorted = np.sort(np.mod(right_values - cut, 1.0))
    displacement = np.abs(left_sorted - right_sorted)
    return float(np.max(np.minimum(displacement, 1.0 - displacement)))


def z2pack_move_ratio(left: Sequence[float], right: Sequence[float]) -> float:
    """Maximum matched move divided by the smaller neighboring gap size."""

    left_gap = largest_circular_gap(left)[1]
    right_gap = largest_circular_gap(right)[1]
    denominator = min(left_gap, right_gap)
    return z2pack_max_move(left, right) / denominator if denominator > 0.0 else math.inf


def z2pack_gap_ratio(left: Sequence[float], right: Sequence[float]) -> float:
    """Bidirectional neighbor-to-largest-gap clearance ratio."""

    left_values = np.mod(np.asarray(left, dtype=float), 1.0)
    right_values = np.mod(np.asarray(right, dtype=float), 1.0)
    left_center, left_gap, _ = largest_circular_gap(left_values)
    right_center, right_gap, _ = largest_circular_gap(right_values)
    if left_gap <= 0.0 or right_gap <= 0.0:
        return 0.0
    # GAP(l1,l2) and GAP(l2,l1), as in the adaptive Z2Pack control.
    right_away_from_left = float(np.min(np.abs(right_values - left_center)) / left_gap)
    left_away_from_right = float(np.min(np.abs(left_values - right_center)) / right_gap)
    return min(right_away_from_left, left_away_from_right)


def _orthonormal_basis(frame: np.ndarray, rank: int) -> np.ndarray:
    left, singular, _ = np.linalg.svd(
        np.asarray(frame, dtype=np.complex128), full_matrices=False
    )
    if singular.size < rank or singular[rank - 1] <= 1.0e-12:
        raise ValueError("restricted occupied frame loses rank")
    return left[:, :rank]


def projector_distance(left: np.ndarray, right: np.ndarray) -> float:
    """Spectral distance between equal-rank projectors without forming D x D.

    If ``s_min`` is the smallest singular value of ``Q_left^H Q_right``, the
    projector distance is ``sqrt(1-s_min**2)``.  A residual evaluation is used
    close to coincident subspaces to avoid cancellation in that square root;
    both paths require only O(D*r) memory for basis dimension D and rank r.
    """

    if left.shape[1] != right.shape[1]:
        raise ValueError("projector frames have different ranks")
    qleft = _orthonormal_basis(left, left.shape[1])
    qright = _orthonormal_basis(right, right.shape[1])
    overlap = qleft.conj().T @ qright
    min_singular = float(np.linalg.svd(overlap, compute_uv=False)[-1])
    min_singular = min(1.0, max(0.0, min_singular))
    distance_squared = max(0.0, (1.0 - min_singular) * (1.0 + min_singular))
    if distance_squared > 1.0e-12:
        return float(math.sqrt(distance_squared))
    # Stable near P_left == P_right, where 1-s_min loses relative precision.
    residual = qright - qleft @ overlap
    return float(np.linalg.norm(residual, 2))


def frame_link_field(frames: np.ndarray, axis: int, stride: int = 1) -> LinkField:
    """Build link polar factors from analytic orthonormal occupied frames."""

    values = np.asarray(frames, dtype=np.complex128)
    if values.ndim != 4:
        raise ValueError("frames must have shape [nx, ny, basis, occupied]")
    nx, ny, _, nocc = values.shape
    length = nx if axis == 0 else ny
    if axis not in (0, 1) or stride < 1 or length % stride:
        raise ValueError("loop stride must divide the selected mesh dimension")
    shape = (nx, ny)
    unitary = np.empty(shape + (nocc, nocc), dtype=np.complex128)
    min_sv = np.empty(shape)
    max_sv = np.empty(shape)
    residual = np.empty(shape)
    for ix in range(nx):
        for iy in range(ny):
            jx = (ix + stride) % nx if axis == 0 else ix
            jy = (iy + stride) % ny if axis == 1 else iy
            overlap = values[ix, iy].conj().T @ values[jx, jy]
            polar, singular = polar_unitary(overlap)
            unitary[ix, iy] = polar
            min_sv[ix, iy] = singular[-1]
            max_sv[ix, iy] = singular[0]
            residual[ix, iy] = np.linalg.norm(
                polar.conj().T @ polar - np.eye(nocc), ord=2
            )
    ones = np.ones(shape)
    return LinkField(axis, stride, unitary, min_sv, max_sv, ones, ones, residual)


def wavecar_link_field(
    wavecar: Wavecar,
    grid: Grid,
    bands: Sequence[int],
    axis: int,
    stride: int = 1,
) -> LinkField:
    """Build SVD polar links using the shared WAVECAR/G-vector reader."""

    length = grid.nx if axis == 0 else grid.ny
    if axis not in (0, 1) or stride < 1 or length % stride:
        raise ValueError("loop stride must divide the selected mesh dimension")
    nocc = len(bands)
    shape = (grid.nx, grid.ny)
    unitary = np.empty(shape + (nocc, nocc), dtype=np.complex128)
    min_sv = np.empty(shape)
    max_sv = np.empty(shape)
    cov_left = np.empty(shape)
    cov_right = np.empty(shape)
    residual = np.empty(shape)
    step = np.zeros(3)
    step[axis] = stride / length
    for ix in range(grid.nx):
        for iy in range(grid.ny):
            jx = (ix + stride) % grid.nx if axis == 0 else ix
            jy = (iy + stride) % grid.ny if axis == 1 else iy
            ik = int(grid.index[ix, iy])
            jk = int(grid.index[jx, jy])
            overlap, left_coverage, right_coverage = link_matrix(
                wavecar, ik, jk, bands, wavecar.kpoints[ik] + step
            )
            polar, singular = polar_unitary(overlap)
            unitary[ix, iy] = polar
            min_sv[ix, iy] = singular[-1]
            max_sv[ix, iy] = singular[0]
            cov_left[ix, iy] = left_coverage
            cov_right[ix, iy] = right_coverage
            residual[ix, iy] = np.linalg.norm(
                polar.conj().T @ polar - np.eye(nocc), ord=2
            )
    return LinkField(
        axis, stride, unitary, min_sv, max_sv, cov_left, cov_right, residual
    )


def spectrum_from_links(field: LinkField, pump_offset: float = 0.0) -> WilsonSpectrum:
    """Multiply a closed string of polar links and return eigenphases in cycles."""

    nx, ny, nocc, _ = field.unitary.shape
    loop_length = nx if field.axis == 0 else ny
    pump_length = ny if field.axis == 0 else nx
    indices = range(0, loop_length, field.stride)
    wcc = np.empty((pump_length, nocc))
    for pump in range(pump_length):
        wilson = np.eye(nocc, dtype=np.complex128)
        for loop in indices:
            matrix = field.unitary[loop, pump] if field.axis == 0 else field.unitary[pump, loop]
            wilson = wilson @ matrix
        wcc[pump] = np.sort(np.mod(np.angle(np.linalg.eigvals(wilson)) / (2.0 * np.pi), 1.0))
    pump_values = np.mod(pump_offset + np.arange(pump_length) / pump_length, 1.0)
    return WilsonSpectrum(field.axis, 1 - field.axis, pump_values, wcc)


def link_flux_and_chern(xlinks: LinkField, ylinks: LinkField) -> tuple[np.ndarray, float]:
    if xlinks.stride != 1 or ylinks.stride != 1:
        raise ValueError("Chern flux requires nearest-neighbor links")
    nx, ny = xlinks.unitary.shape[:2]
    phase_x = np.angle(np.linalg.det(xlinks.unitary))
    phase_y = np.angle(np.linalg.det(ylinks.unitary))
    flux = np.empty((nx, ny))
    for ix in range(nx):
        for iy in range(ny):
            loop = (
                phase_x[ix, iy]
                + phase_y[(ix + 1) % nx, iy]
                - phase_x[ix, (iy + 1) % ny]
                - phase_y[ix, iy]
            )
            flux[ix, iy] = -float(wrap_angle(loop))
    return flux, float(np.sum(flux) / (2.0 * np.pi))


def flux_odd_residual(flux: np.ndarray) -> float:
    nx, ny = flux.shape
    return float(
        max(
            abs(float(wrap_angle(flux[ix, iy] + flux[(-ix - 1) % nx, (-iy - 1) % ny])))
            for ix in range(nx)
            for iy in range(ny)
        )
    )


def frame_projector_tr_residual(frames: np.ndarray, theta: np.ndarray) -> float:
    """Projector TR residual for analytic frames and a unitary Theta matrix."""

    values = np.asarray(frames, dtype=np.complex128)
    theta = np.asarray(theta, dtype=np.complex128)
    nx, ny, basis, _ = values.shape
    if theta.shape != (basis, basis):
        raise ValueError("Theta matrix shape does not match the frame basis")
    return float(
        max(
            projector_distance(theta @ values[ix, iy].conj(), values[-ix % nx, -iy % ny])
            for ix in range(nx)
            for iy in range(ny)
        )
    )


def _grid_tr_partner(grid: Grid, wavecar: Wavecar, ix: int, iy: int) -> tuple[int, int]:
    """Return the in-plane TR partner after kz=0 modulo G is hard-gated."""

    target = np.mod(-wavecar.kpoints[int(grid.index[ix, iy]), :2] - grid.offset[:2], 1.0)
    x_float, y_float = target * np.asarray((grid.nx, grid.ny))
    jx, jy = int(round(x_float)) % grid.nx, int(round(y_float)) % grid.ny
    if max(abs(x_float - round(x_float)), abs(y_float - round(y_float))) > 1.0e-6:
        raise ValueError("uniform mesh is not closed under k -> -k")
    return jx, jy


def wavecar_time_reversal_diagnostics(
    wavecar: Wavecar,
    grid: Grid,
    bands: Sequence[int],
) -> dict[str, float]:
    """Evaluate energy, Kramers-energy, and occupied-projector TR residuals."""

    if wavecar.spinor_components != 2 or wavecar.header.ispin != 1:
        raise ValueError("Z2 validation requires an ISPIN=1 two-component spinor WAVECAR")
    nocc = len(bands)
    max_energy_error = 0.0
    max_kramers_split = 0.0
    max_projector_error = 0.0
    min_tr_coverage = 1.0
    min_tr_norm_coverage = 1.0
    checked_pairs: set[tuple[int, int]] = set()
    sentinel = min(max(bands) + 1, wavecar.header.nbands)
    energy_columns = np.arange(0, sentinel)

    for ix in range(grid.nx):
        for iy in range(grid.ny):
            jx, jy = _grid_tr_partner(grid, wavecar, ix, iy)
            ik, jk = int(grid.index[ix, iy]), int(grid.index[jx, jy])
            max_energy_error = max(
                max_energy_error,
                float(np.max(np.abs(wavecar.energies[ik, energy_columns] - wavecar.energies[jk, energy_columns]))),
            )
            if ik == jk:
                occupied_energy = wavecar.energies[ik, np.asarray(bands) - 1]
                max_kramers_split = max(
                    max_kramers_split,
                    float(np.max(np.abs(occupied_energy[0::2] - occupied_energy[1::2]))),
                )
            key = tuple(sorted((ik, jk)))
            if key in checked_pairs:
                continue
            checked_pairs.add(key)

            source_g = wavecar.g_vectors(ik).astype(np.int64)
            target_g = wavecar.g_vectors(jk).astype(np.int64)
            shift = np.rint(-wavecar.kpoints[ik] - wavecar.kpoints[jk]).astype(np.int64)
            lookup = {tuple(int(x) for x in g): index for index, g in enumerate(target_g)}
            source_indices: list[int] = []
            target_indices: list[int] = []
            for source_index, gvector in enumerate(source_g):
                target_index = lookup.get(tuple(int(x) for x in shift - gvector))
                if target_index is not None:
                    source_indices.append(source_index)
                    target_indices.append(target_index)
            if not source_indices:
                raise ValueError("time-reversal projector comparison has no common plane waves")
            min_tr_coverage = min(
                min_tr_coverage,
                len(source_indices) / source_g.shape[0],
                len(target_indices) / target_g.shape[0],
            )
            source_full = wavecar.coefficients(ik, bands)
            target_full = wavecar.coefficients(jk, bands)
            source = source_full[:, :, source_indices]
            target_frame = target_full[:, :, target_indices]
            for restricted, complete in (
                (source, source_full),
                (target_frame, target_full),
            ):
                complete_norm = np.sum(np.abs(complete) ** 2, axis=(1, 2))
                restricted_norm = np.sum(np.abs(restricted) ** 2, axis=(1, 2))
                if np.any(complete_norm <= 0.0):
                    raise ValueError("occupied spinor coefficient has zero raw norm")
                min_tr_norm_coverage = min(
                    min_tr_norm_coverage,
                    float(np.min(restricted_norm / complete_norm)),
                )
            theta_source = np.stack((-source[:, 1].conj(), source[:, 0].conj()), axis=1)
            theta_columns = theta_source.transpose(1, 2, 0).reshape(-1, nocc)
            target_columns = target_frame.transpose(1, 2, 0).reshape(-1, nocc)
            max_projector_error = max(
                max_projector_error, projector_distance(theta_columns, target_columns)
            )
    return {
        "max_energy_tr_error_ev": float(max_energy_error),
        "max_trim_kramers_energy_split_ev": float(max_kramers_split),
        "max_projector_tr_residual": float(max_projector_error),
        "min_tr_g_basis_coverage": float(min_tr_coverage),
        "min_tr_raw_norm_coverage": float(min_tr_norm_coverage),
    }


def _check(name: str, value: float | bool, threshold: float | None, relation: str) -> dict[str, object]:
    if isinstance(value, (bool, np.bool_)):
        passed = bool(value)
    elif relation == ">=":
        passed = bool(float(value) >= float(threshold))
    elif relation == "<=":
        passed = bool(float(value) <= float(threshold))
    elif relation == ">":
        passed = bool(float(value) > float(threshold))
    elif relation == "<":
        passed = bool(float(value) < float(threshold))
    elif relation == "==":
        passed = bool(float(value) == float(threshold))
    else:  # pragma: no cover - internal programming guard
        raise ValueError(f"unknown check relation {relation}")
    return {
        "name": name,
        "pass": passed,
        "value": value,
        "relation": relation,
        "threshold": threshold,
    }


def _axis_name(axis: int) -> str:
    return "x" if axis == 0 else "y"


def analyze_link_fields(
    xlinks: LinkField,
    ylinks: LinkField,
    stride2: dict[int, LinkField],
    pump_offsets: Sequence[float],
    thresholds: Thresholds,
    requested_axes: Sequence[int],
    extra_metrics: dict[str, float],
) -> dict[str, object]:
    """Combine Wilson-loop results and all strict validation gates."""

    flux, chern = link_flux_and_chern(xlinks, ylinks)
    spectra: dict[int, WilsonSpectrum] = {}
    axis_reports: dict[str, object] = {}
    checks: list[dict[str, object]] = []

    min_link = min(float(np.min(xlinks.min_singular)), float(np.min(ylinks.min_singular)))
    min_coverage = min(
        float(np.min(xlinks.coverage_left)),
        float(np.min(xlinks.coverage_right)),
        float(np.min(ylinks.coverage_left)),
        float(np.min(ylinks.coverage_right)),
    )
    max_unitarity = max(
        float(np.max(xlinks.unitarity_residual)), float(np.max(ylinks.unitarity_residual))
    )
    checks.extend(
        (
            _check("occupied_subspace_chern_zero", abs(chern), thresholds.max_chern_residual, "<="),
            _check("berry_flux_time_reversal_odd", flux_odd_residual(flux), thresholds.max_flux_odd_residual_rad, "<="),
            _check("link_subspace_nonsingular", min_link, thresholds.min_link_singular_value, ">="),
            _check("plane_wave_link_coverage", min_coverage, thresholds.min_plane_wave_coverage, ">="),
            _check("polar_link_unitarity", max_unitarity, thresholds.max_unitarity_residual, "<="),
        )
    )

    candidates: list[int] = []
    for axis in requested_axes:
        full_field = xlinks if axis == 0 else ylinks
        full = spectrum_from_links(full_field, pump_offsets[1 - axis])
        coarse = spectrum_from_links(stride2[axis], pump_offsets[1 - axis])
        spectra[axis] = full
        half_pump, half_wcc = half_bz_spectrum(full)
        candidate, gap_center, gap_size = z2pack_largest_gap_parity(half_wcc)
        candidates.append(candidate)
        tr_set, kramers = wcc_time_reversal_diagnostics(full)
        convergence = fixed_grid_convergence(full, coarse)
        prefix = f"{_axis_name(axis)}_loop"
        axis_checks = [
            _check(f"{prefix}_wcc_time_reversal_set", tr_set, thresholds.max_wcc_tr_set_distance, "<="),
            _check(f"{prefix}_endpoint_kramers_pairing", kramers, thresholds.max_endpoint_kramers_split, "<="),
            _check(
                f"{prefix}_wcc_position_converged",
                convergence["full_vs_stride2_wcc_distance_cycles"],
                thresholds.max_wcc_position_change,
                "<",
            ),
            _check(
                f"{prefix}_wcc_gap_neighbor_converged",
                convergence["neighbor_to_largest_gap_ratio"],
                thresholds.min_combined_gap_ratio,
                ">",
            ),
            _check(
                f"{prefix}_wcc_move_neighbor_converged",
                convergence["neighbor_move_ratio"],
                thresholds.max_neighbor_move_ratio,
                "<",
            ),
        ]
        checks.extend(axis_checks)
        axis_reports[_axis_name(axis)] = {
            "loop_axis": _axis_name(axis),
            "pump_axis": _axis_name(1 - axis),
            "candidate_z2": int(candidate),
            "max_wcc_tr_set_distance_cycles": tr_set,
            "max_endpoint_kramers_split_cycles": kramers,
            **convergence,
            "half_bz_pump": half_pump,
            "largest_gap_center": gap_center,
            "largest_gap_size": gap_size,
        }

    if len(candidates) > 1:
        checks.append(_check("loop_axis_z2_agreement", len(set(candidates)) == 1, None, "=="))
    checks.append(
        _check("both_loop_axes_evaluated", set(requested_axes) == {0, 1}, None, "==")
    )

    metric_checks = (
        ("direct_gap_open", "direct_gap_ev", thresholds.min_direct_gap_ev, ">="),
        ("indirect_gap_open", "indirect_gap_ev", thresholds.min_indirect_gap_ev, ">="),
        ("band_energy_time_reversal", "max_energy_tr_error_ev", thresholds.max_energy_tr_error_ev, "<="),
        (
            "trim_kramers_energy_pairing",
            "max_trim_kramers_energy_split_ev",
            thresholds.max_trim_kramers_energy_split_ev,
            "<=",
        ),
        ("occupied_projector_time_reversal", "max_projector_tr_residual", thresholds.max_projector_tr_residual, "<="),
        (
            "time_reversal_g_basis_complete",
            "min_tr_g_basis_coverage",
            thresholds.min_tr_g_basis_coverage,
            ">=",
        ),
        (
            "time_reversal_raw_norm_coverage",
            "min_tr_raw_norm_coverage",
            thresholds.min_tr_raw_norm_coverage,
            ">=",
        ),
    )
    for check_name, metric_name, threshold, relation in metric_checks:
        if metric_name in extra_metrics:
            checks.append(_check(check_name, extra_metrics[metric_name], threshold, relation))

    valid = all(bool(item["pass"]) for item in checks)
    candidate = candidates[0] if candidates and len(set(candidates)) == 1 else None
    failed = [str(item["name"]) for item in checks if not item["pass"]]
    return {
        "valid": valid,
        "z2": int(candidate) if valid and candidate is not None else None,
        "candidate_z2": int(candidate) if candidate is not None else None,
        "failed_guards": failed,
        "checks": checks,
        "metrics": {
            **extra_metrics,
            "chern": chern,
            "max_flux_odd_residual_rad": flux_odd_residual(flux),
            "min_link_singular_value": min_link,
            "min_plane_wave_link_coverage": min_coverage,
            "max_polar_unitarity_residual": max_unitarity,
        },
        "axes": axis_reports,
        "_spectra": spectra,
    }


def analyze_analytic_frames(
    frames: np.ndarray,
    energies: np.ndarray,
    theta: np.ndarray,
    axes: Sequence[int] = (0, 1),
    thresholds: Thresholds | None = None,
) -> dict[str, object]:
    """Public model oracle used by tests and method-validation notebooks."""

    thresholds = thresholds or Thresholds()
    nx, ny, _, nocc = frames.shape
    if nx < 4 or ny < 4 or nx % 2 or ny % 2 or nocc % 2:
        raise ValueError("Z2 model mesh must be even, at least 4x4, with even occupied rank")
    if energies.shape[:2] != (nx, ny) or energies.shape[2] <= nocc:
        raise ValueError("model energies require at least one unoccupied sentinel band")
    xlinks = frame_link_field(frames, 0)
    ylinks = frame_link_field(frames, 1)
    stride2 = {0: frame_link_field(frames, 0, 2), 1: frame_link_field(frames, 1, 2)}
    direct = float(np.min(energies[:, :, nocc] - energies[:, :, nocc - 1]))
    indirect = float(np.min(energies[:, :, nocc]) - np.max(energies[:, :, nocc - 1]))
    energy_tr = max(
        float(np.max(np.abs(energies[ix, iy] - energies[-ix % nx, -iy % ny])))
        for ix in range(nx)
        for iy in range(ny)
    )
    kramers_energy = max(
        float(
            np.max(
                np.abs(
                    energies[ix, iy, :nocc:2]
                    - energies[ix, iy, 1:nocc:2]
                )
            )
        )
        for ix in (0, nx // 2)
        for iy in (0, ny // 2)
    )
    metrics = {
        "direct_gap_ev": direct,
        "indirect_gap_ev": indirect,
        "max_energy_tr_error_ev": energy_tr,
        "max_trim_kramers_energy_split_ev": kramers_energy,
        "max_projector_tr_residual": frame_projector_tr_residual(frames, theta),
        "min_tr_g_basis_coverage": 1.0,
        "min_tr_raw_norm_coverage": 1.0,
    }
    return analyze_link_fields(
        xlinks, ylinks, stride2, (0.0, 0.0), thresholds, axes, metrics
    )


def _validate_wavecar_inputs(wavecar: Wavecar, grid: Grid, occupied: int) -> None:
    if wavecar.header.ispin != 1 or wavecar.spinor_components != 2:
        raise ValueError("Z2 validation requires ISPIN=1 and --spinor-components 2")
    if occupied < 2 or occupied % 2:
        raise ValueError("--occupied-bands must be a positive even rank of at least 2")
    if occupied >= wavecar.header.nbands:
        raise ValueError("Z2 validation requires band occupied+1 as an unoccupied sentinel")
    if grid.nx < 4 or grid.ny < 4 or grid.nx % 2 or grid.ny % 2:
        raise ValueError("--nx and --ny must be even and at least 4 for Z2 validation")
    kz_offset = float(np.mod(grid.offset[2], 1.0))
    if min(abs(kz_offset), abs(kz_offset - 1.0)) > 1.0e-7:
        raise ValueError(
            "strict 2D Z2 validation requires the Gamma-centered kz=0 plane"
        )
    for axis, length in enumerate((grid.nx, grid.ny)):
        offset = float(np.mod(grid.offset[axis], 1.0))
        if min(abs(offset), abs(offset - 1.0)) > 1.0e-7:
            raise ValueError(
                "strict Wilson-loop validation requires a Gamma-centered mesh "
                "containing the q=0 and q=1/2 time-reversal lines"
            )
        if length % 2:
            raise ValueError("time-reversal endpoint validation requires even grid sizes")


def analyze_wavecar(
    wavecar: Wavecar,
    nx: int,
    ny: int,
    occupied: int,
    axes: Sequence[int],
    thresholds: Thresholds,
) -> dict[str, object]:
    grid = infer_uniform_grid(wavecar.kpoints, nx, ny)
    _validate_wavecar_inputs(wavecar, grid, occupied)
    bands = list(range(1, occupied + 1))
    xlinks = wavecar_link_field(wavecar, grid, bands, 0)
    ylinks = wavecar_link_field(wavecar, grid, bands, 1)
    stride2 = {
        0: wavecar_link_field(wavecar, grid, bands, 0, 2),
        1: wavecar_link_field(wavecar, grid, bands, 1, 2),
    }
    occupied_energy = wavecar.energies[:, occupied - 1]
    sentinel_energy = wavecar.energies[:, occupied]
    metrics = {
        "direct_gap_ev": float(np.min(sentinel_energy - occupied_energy)),
        "indirect_gap_ev": float(np.min(sentinel_energy) - np.max(occupied_energy)),
        **wavecar_time_reversal_diagnostics(wavecar, grid, bands),
    }
    result = analyze_link_fields(
        xlinks,
        ylinks,
        stride2,
        (float(grid.offset[0]), float(grid.offset[1])),
        thresholds,
        axes,
        metrics,
    )
    result["grid"] = {
        "nx": nx,
        "ny": ny,
        "offset": grid.offset,
        "occupied_bands": occupied,
    }
    return result


def _atomic_path(path: Path) -> Path:
    return path.with_name(f".{path.name}.{os.getpid()}.tmp")


def _atomic_write_json(path: Path, payload: object) -> None:
    temporary = _atomic_path(path)
    try:
        with temporary.open("w", encoding="utf-8") as handle:
            json.dump(json_safe(payload), handle, indent=2, allow_nan=False)
            handle.write("\n")
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def write_wcc_csv(path: Path, spectra: dict[int, WilsonSpectrum]) -> None:
    temporary = _atomic_path(path)
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(
                ("loop_axis", "pump_axis", "pump_fractional", "pump_signed", "wcc_index", "wcc_cycles")
            )
            for axis in sorted(spectra):
                spectrum = spectra[axis]
                for pump, row in zip(spectrum.pump, spectrum.wcc):
                    signed = float(circular_delta(pump, 0.0))
                    for branch, wcc in enumerate(row, start=1):
                        writer.writerow(
                            (
                                _axis_name(axis),
                                _axis_name(1 - axis),
                                f"{float(pump):.16g}",
                                f"{signed:.16g}",
                                branch,
                                f"{float(wcc):.16g}",
                            )
                        )
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def _matplotlib_pyplot():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def plot_wcc(path: Path, result: dict[str, object], spectra: dict[int, WilsonSpectrum]) -> None:
    """Plot unordered WCC sets; sorted points are intentionally not connected."""

    plt = _matplotlib_pyplot()
    axes_order = sorted(spectra)
    figure, plot_axes = plt.subplots(
        len(axes_order), 2, figsize=(11.0, 4.2 * len(axes_order)), squeeze=False
    )
    failed = list(result["failed_guards"])
    for row_index, axis in enumerate(axes_order):
        spectrum = spectra[axis]
        axis_report = result["axes"][_axis_name(axis)]
        signed = np.asarray(circular_delta(spectrum.pump, 0.0), dtype=float)
        order = np.argsort(signed)
        left = plot_axes[row_index, 0]
        for branch in range(spectrum.wcc.shape[1]):
            left.scatter(signed[order], spectrum.wcc[order, branch], s=13, color="#2457A7")
        left.set_xlim(-0.52, 0.52)
        left.set_ylim(-0.03, 1.03)
        left.set_xlabel(f"pump q{_axis_name(1-axis)}")
        left.set_ylabel("hybrid WCC (cycles)")
        left.set_title(
            f"{_axis_name(axis)} Wilson loop: q and -q WCC sets\n"
            f"max set distance = {axis_report['max_wcc_tr_set_distance_cycles']:.3e} cycles"
        )
        left.grid(alpha=0.2)

        half_q, half_wcc = half_bz_spectrum(spectrum)
        right = plot_axes[row_index, 1]
        for branch in range(half_wcc.shape[1]):
            right.scatter(half_q, half_wcc[:, branch], s=15, color="#2457A7")
        centers = np.asarray(axis_report["largest_gap_center"])
        right.scatter(half_q, centers, s=18, marker="x", color="#C43B3B", label="largest-gap center")
        right.set_xlim(-0.01, 0.51)
        right.set_ylim(-0.03, 1.03)
        right.set_xlabel(f"half-BZ pump q{_axis_name(1-axis)}")
        right.set_ylabel("hybrid WCC (cycles)")
        status = (
            f"PASS: Z2 = {result['z2']}"
            if result["valid"]
            else f"INVALID (candidate Z2 = {axis_report['candidate_z2']})"
        )
        axis_failures = [item for item in failed if item.startswith(f"{_axis_name(axis)}_loop")]
        suffix = "" if not axis_failures else "\nfailed: " + ", ".join(
            item.removeprefix(f"{_axis_name(axis)}_loop_") for item in axis_failures
        )
        right.set_title(status + suffix)
        right.grid(alpha=0.2)
        right.legend(loc="best", fontsize=8)
    figure.suptitle("Wilson-loop time-reversal and fixed-grid diagnostics")
    figure.tight_layout()
    temporary = path.with_name(f".{path.stem}.{os.getpid()}.tmp{path.suffix}")
    try:
        figure.savefig(temporary, dpi=180)
        temporary.replace(path)
    finally:
        plt.close(figure)
        temporary.unlink(missing_ok=True)


def _resolved(path: Path) -> Path:
    return path.expanduser().absolute().resolve(strict=False)


def validate_output_targets(input_path: Path, outputs: Iterable[Path]) -> None:
    resolved_input = _resolved(input_path)
    resolved_outputs = [_resolved(path) for path in outputs]
    if any(path == resolved_input for path in resolved_outputs):
        raise ValueError("an output path would overwrite the input WAVECAR")
    if len(set(resolved_outputs)) != len(resolved_outputs):
        raise ValueError("Z2 output paths must be distinct")


def _parse_axes(value: str) -> tuple[int, ...]:
    return {"x": (0,), "y": (1,), "both": (0, 1)}[value]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Guarded Wilson-loop Z2 validation from a full-mesh spinor WAVECAR."
    )
    parser.add_argument("wavecar", type=Path)
    parser.add_argument("--nx", type=int, required=True)
    parser.add_argument("--ny", type=int, required=True)
    parser.add_argument("--occupied-bands", type=int, required=True, metavar="N")
    parser.add_argument("--axis", choices=("x", "y", "both"), default="both")
    parser.add_argument("--spinor-components", type=int, choices=(1, 2), default=2)
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument("--csv-name", default=DEFAULT_CSV)
    parser.add_argument("--json-name", default=DEFAULT_JSON)
    parser.add_argument("--plot-name", default=DEFAULT_PLOT)
    parser.add_argument("--plot", action="store_true")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    output_dir = args.output_dir
    csv_path = output_dir / args.csv_name
    json_path = output_dir / args.json_name
    plot_path = output_dir / args.plot_name
    planned = [csv_path, json_path] + ([plot_path] if args.plot else [])
    cleanup_authorized = False
    try:
        validate_output_targets(args.wavecar, planned)
        cleanup_authorized = True
        if not args.wavecar.is_file():
            raise ValueError(f"WAVECAR does not exist or is not a file: {args.wavecar}")
        output_dir.mkdir(parents=True, exist_ok=True)
        for path in planned:
            path.unlink(missing_ok=True)
        wavecar = Wavecar(args.wavecar, spinor_components=args.spinor_components)
        result = analyze_wavecar(
            wavecar,
            args.nx,
            args.ny,
            args.occupied_bands,
            _parse_axes(args.axis),
            Thresholds(),
        )
        spectra = result.pop("_spectra")
        public_result = {
            "schema_version": 1,
            "tool": "wavecar_z2",
            "tool_version": __version__,
            "input": args.wavecar.name,
            "method": "SVD polar-unitary Wilson loops with largest-gap parity",
            "interpretation": (
                "z2 is reportable only when valid is true; candidate_z2 is diagnostic only"
            ),
            "thresholds": asdict(Thresholds()),
            **result,
        }
        write_wcc_csv(csv_path, spectra)
        _atomic_write_json(json_path, public_result)
        if args.plot:
            plot_wcc(plot_path, public_result, spectra)
        print(
            f"{'PASS' if result['valid'] else 'INVALID'}: "
            f"z2={result['z2']}, candidate_z2={result['candidate_z2']}"
        )
        if result["failed_guards"]:
            print("failed guards: " + ", ".join(result["failed_guards"]))
        return 0 if result["valid"] else 2
    except Exception as error:
        if cleanup_authorized:
            for path in planned:
                path.unlink(missing_ok=True)
                _atomic_path(path).unlink(missing_ok=True)
        print(f"error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
