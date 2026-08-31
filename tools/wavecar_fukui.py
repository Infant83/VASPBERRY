#!/usr/bin/env python3
"""Read a VASP WAVECAR and evaluate 2D Fukui plaquette phases.

This is a deliberately small, read-only validation program.  Its conventions
match VASPBERRY 1.0:

* VASP's g3/g2/g1 plane-wave ordering is reproduced exactly.
* Spinor records are interpreted as all up coefficients followed by all down.
* The loop is k -> k+dk1 -> k+dk1+dk2 -> k+dk2 -> k.
* The reported flux is -arg(product of link determinants), and
  Omega_z = flux / (|b1 x b2| / (nkx*nky)).

The program handles both byte-RECL and legacy four-byte-word-RECL WAVECARs by
validating record 2 at both possible physical strides.  It never writes to the
WAVECAR.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np


__version__ = "1.1.1"
KINETIC_C = 0.262465831  # 2 m_e / hbar^2, 1/(eV Angstrom^2)
VALLEY_DISTANCE_TIE_ATOL = 1.0e-12  # inverse Angstrom
TRANSPORT_OUTPUT_NAMES = (
    "transport_t0.csv",
    "transport_t0_diagnostics.json",
    "wavecar_fukui_sigma_mu.png",
)
SHORTEST_K_KP_LINE_TITLE = (
    "Shortest periodic-image K$\\rightarrow$K' cut "
    "(not the K$\\rightarrow\\Gamma\\rightarrow$K' path); "
    "V+1 is diagnostic if link quality fails"
)
FULL_TRANSPORT_OUTPUT_NAMES = (
    "transport_full_t0.csv",
    "transport_full_t0_diagnostics.json",
    "wavecar_fukui_sigma_full_mu.png",
)


def json_safe(value: object) -> object:
    """Recursively replace nonfinite numerics so emitted JSON is strict."""

    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return json_safe(value.tolist())
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        number = float(value)
        return number if math.isfinite(number) else None
    return value


def write_strict_json(path: Path, payload: object) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(json_safe(payload), handle, indent=2, allow_nan=False)
        handle.write("\n")


def fortran_implicit_integer(value: float) -> int:
    """Match assignment to VASPBERRY's implicit INTEGER n* variables."""

    return int(value)


def signed_vasp_sequence(limit: int) -> list[int]:
    """Return 0,1,...,N,-N,...,-1, VASP/VASPBERRY ordering."""

    return list(range(limit + 1)) + list(range(-limit, 0))


@dataclass(frozen=True)
class WavecarHeader:
    logical_recl: int
    stride_bytes: int
    ispin: int
    rtag: int
    nkpoints: int
    nbands: int
    encut_ev: float
    lattice: np.ndarray
    reciprocal: np.ndarray


class Wavecar:
    def __init__(self, path: Path, spinor_components: int = 2, spin: int = 1):
        self.path = path
        self.spinor_components = spinor_components
        self.spin = spin
        self.header = self._read_header()
        if self.header.rtag != 45200:
            raise ValueError(
                f"RTAG={self.header.rtag} is unsupported; this validator expects "
                "single-precision complex coefficients (RTAG=45200)."
            )
        if self.header.ispin == 2 and spinor_components != 1:
            raise ValueError("ISPIN=2 WAVECAR requires --spinor-components 1")
        if not (1 <= spin <= self.header.ispin):
            raise ValueError(f"spin channel {spin} is invalid for ISPIN={self.header.ispin}")
        self.kpoints, self.nplane, self.energies, self.occupations = (
            self._read_k_headers()
        )
        if np.any(self.nplane <= 0):
            bad = int(np.flatnonzero(self.nplane <= 0)[0] + 1)
            raise ValueError(f"NPLANE must be positive at k index {bad}")
        for label, values in (
            ("k-point coordinates", self.kpoints),
            ("band energies", self.energies),
            ("band occupations", self.occupations),
        ):
            if not np.all(np.isfinite(values)):
                raise ValueError(f"WAVECAR {label} must all be finite")
        self.nbmax = self._reciproperty_bounds()
        self._g_cache: dict[int, np.ndarray] = {}

    def _read_header(self) -> WavecarHeader:
        with self.path.open("rb") as handle:
            rec1 = np.fromfile(handle, dtype="<f8", count=3)
        if rec1.size != 3:
            raise ValueError("WAVECAR is too short for record 1")
        logical_recl, ispin, rtag = (int(round(x)) for x in rec1)
        if ispin not in (1, 2):
            raise ValueError(f"invalid ISPIN={ispin}; expected 1 or 2")
        candidates: list[tuple[int, np.ndarray]] = []
        file_size = self.path.stat().st_size
        for stride in dict.fromkeys((logical_recl, 4 * logical_recl)):
            if stride < 96 or stride >= file_size:
                continue
            with self.path.open("rb") as handle:
                handle.seek(stride)
                rec2 = np.fromfile(handle, dtype="<f8", count=12)
            if rec2.size != 12:
                continue
            nk, nb = int(round(rec2[0])), int(round(rec2[1]))
            encut = float(rec2[2])
            lattice = rec2[3:12].reshape(3, 3)
            volume = float(np.dot(lattice[0], np.cross(lattice[1], lattice[2])))
            expected_min = stride * (
                2 + ispin * max(nk, 0) * (max(nb, 0) + 1)
            )
            plausible = (
                0 < nk < 10_000_000
                and 0 < nb < 1_000_000
                and 0.0 < encut < 100_000.0
                and abs(volume) > 1.0e-12
                and expected_min <= file_size
            )
            if plausible:
                candidates.append((stride, rec2))
        if len(candidates) != 1:
            strides = [item[0] for item in candidates]
            raise ValueError(f"could not identify a unique physical RECL stride: {strides}")
        stride, rec2 = candidates[0]
        lattice = rec2[3:12].reshape(3, 3).copy()
        volume = float(np.dot(lattice[0], np.cross(lattice[1], lattice[2])))
        reciprocal = np.stack(
            (
                2.0 * np.pi * np.cross(lattice[1], lattice[2]) / volume,
                2.0 * np.pi * np.cross(lattice[2], lattice[0]) / volume,
                2.0 * np.pi * np.cross(lattice[0], lattice[1]) / volume,
            )
        )
        return WavecarHeader(
            logical_recl=logical_recl,
            stride_bytes=stride,
            ispin=ispin,
            rtag=rtag,
            nkpoints=int(round(rec2[0])),
            nbands=int(round(rec2[1])),
            encut_ev=float(rec2[2]),
            lattice=lattice,
            reciprocal=reciprocal,
        )

    def record_number(self, ik: int, band: int | None = None, spin: int = 1) -> int:
        """Return VASPBERRY's one-based direct-access record number."""

        h = self.header
        if not (0 <= ik < h.nkpoints):
            raise IndexError(ik)
        if not (1 <= spin <= h.ispin):
            raise IndexError(spin)
        rec = 3 + ik * (h.nbands + 1) + h.nkpoints * (h.nbands + 1) * (spin - 1)
        if band is not None:
            if not (1 <= band <= h.nbands):
                raise IndexError(band)
            rec += band
        return rec

    def _read_k_headers(self) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        h = self.header
        kpoints = np.empty((h.nkpoints, 3), dtype=np.float64)
        nplane = np.empty(h.nkpoints, dtype=np.int64)
        energies = np.empty((h.nkpoints, h.nbands), dtype=np.float64)
        occupations = np.empty_like(energies)
        with self.path.open("rb") as handle:
            for ik in range(h.nkpoints):
                offset = (self.record_number(ik, spin=self.spin) - 1) * h.stride_bytes
                handle.seek(offset)
                raw = np.fromfile(handle, dtype="<f8", count=4 + 3 * h.nbands)
                if raw.size != 4 + 3 * h.nbands:
                    raise ValueError(f"truncated k header at k index {ik + 1}")
                nplane[ik] = int(round(raw[0]))
                kpoints[ik] = raw[1:4]
                energies[ik] = raw[4::3]
                occupations[ik] = raw[6::3]
        return kpoints, nplane, energies, occupations

    def _reciproperty_bounds(self) -> tuple[int, int, int]:
        b1, b2, b3 = self.header.reciprocal
        bmag = np.linalg.norm(self.header.reciprocal, axis=1)
        root = math.sqrt(self.header.encut_ev * KINETIC_C)

        phi12 = math.acos(float(np.dot(b1, b2) / (bmag[0] * bmag[1])))
        vtmp = np.cross(b1, b2)
        sin123 = float(np.dot(b3, vtmp) / (np.linalg.norm(vtmp) * bmag[2]))
        bound_a = (
            root / (bmag[0] * abs(math.sin(phi12))) + 1.0,
            root / (bmag[1] * abs(math.sin(phi12))) + 1.0,
            root / (bmag[2] * abs(sin123)) + 1.0,
        )

        phi13 = math.acos(float(np.dot(b1, b3) / (bmag[0] * bmag[2])))
        vtmp = np.cross(b1, b3)
        sin123 = float(np.dot(b2, vtmp) / (np.linalg.norm(vtmp) * bmag[1]))
        bound_b = (
            root / (bmag[0] * abs(math.sin(phi13))) + 1.0,
            root / (bmag[1] * abs(sin123)) + 1.0,
            root / (bmag[2] * abs(math.sin(phi13))) + 1.0,
        )

        phi23 = math.acos(float(np.dot(b2, b3) / (bmag[1] * bmag[2])))
        vtmp = np.cross(b2, b3)
        sin123 = float(np.dot(b1, vtmp) / (np.linalg.norm(vtmp) * bmag[0]))
        bound_c = (
            root / (bmag[0] * abs(sin123)) + 1.0,
            root / (bmag[1] * abs(math.sin(phi23))) + 1.0,
            root / (bmag[2] * abs(math.sin(phi23))) + 1.0,
        )

        return tuple(
            max(
                fortran_implicit_integer(bound_a[i]),
                fortran_implicit_integer(bound_b[i]),
                fortran_implicit_integer(bound_c[i]),
            )
            for i in range(3)
        )

    def g_vectors(self, ik: int) -> np.ndarray:
        if ik in self._g_cache:
            return self._g_cache[ik]
        n1, n2, n3 = self.nbmax
        kpoint = self.kpoints[ik]
        reciprocal = self.header.reciprocal
        accepted: list[tuple[int, int, int]] = []
        for g3 in signed_vasp_sequence(n3):
            for g2 in signed_vasp_sequence(n2):
                for g1 in signed_vasp_sequence(n1):
                    cart = (kpoint + np.array((g1, g2, g3))) @ reciprocal
                    if float(np.dot(cart, cart)) / KINETIC_C < self.header.encut_ev:
                        accepted.append((g1, g2, g3))
        result = np.asarray(accepted, dtype=np.int16)
        expected = int(self.nplane[ik]) // self.spinor_components
        if expected * self.spinor_components != int(self.nplane[ik]):
            raise ValueError(f"NPLANE at k={ik + 1} is incompatible with spinor count")
        if result.shape[0] != expected:
            raise ValueError(
                f"G-vector count mismatch at k={ik + 1}: generated {result.shape[0]}, "
                f"expected {expected}"
            )
        self._g_cache[ik] = result
        return result

    def coefficients(self, ik: int, bands: Sequence[int]) -> np.ndarray:
        """Return [band, spinor, G] coefficients in WAVECAR ordering."""

        ng = self.g_vectors(ik).shape[0]
        result = np.empty((len(bands), self.spinor_components, ng), dtype=np.complex64)
        with self.path.open("rb") as handle:
            for row, band in enumerate(bands):
                offset = (
                    self.record_number(ik, band, self.spin) - 1
                ) * self.header.stride_bytes
                handle.seek(offset)
                raw = np.fromfile(handle, dtype="<c8", count=int(self.nplane[ik]))
                if raw.size != int(self.nplane[ik]):
                    raise ValueError(f"truncated coefficient record k={ik + 1}, band={band}")
                if not np.all(np.isfinite(raw)):
                    raise ValueError(
                        f"nonfinite coefficient record k={ik + 1}, band={band}"
                    )
                result[row] = raw.reshape(self.spinor_components, ng)
        return result


@dataclass(frozen=True)
class Grid:
    nx: int
    ny: int
    index: np.ndarray
    offset: np.ndarray


def infer_uniform_grid(kpoints: np.ndarray, nx: int, ny: int, tolerance: float = 1e-7) -> Grid:
    if nx < 2 or ny < 2:
        raise ValueError("--nx and --ny must both be at least 2 for plaquettes")
    if not np.all(np.isfinite(kpoints)):
        raise ValueError("k-point coordinates must all be finite")
    if kpoints.shape[0] != nx * ny:
        raise ValueError(f"NKPTS={kpoints.shape[0]} is not nx*ny={nx * ny}")
    if np.max(np.abs(kpoints[:, 2] - kpoints[0, 2])) > tolerance:
        raise ValueError("kz is not constant; this is not a 2D mesh")
    steps = np.asarray((1.0 / nx, 1.0 / ny))
    modulo = np.mod(kpoints[0, :2], steps)
    modulo[np.isclose(modulo, steps, atol=tolerance)] = 0.0
    offset = np.asarray((modulo[0], modulo[1], kpoints[0, 2]))
    index = np.full((nx, ny), -1, dtype=np.int64)
    for ik, kpoint in enumerate(kpoints):
        ix_float = ((kpoint[0] - offset[0]) % 1.0) * nx
        iy_float = ((kpoint[1] - offset[1]) % 1.0) * ny
        ix = int(round(ix_float)) % nx
        iy = int(round(iy_float)) % ny
        if abs(ix_float - round(ix_float)) > tolerance * nx:
            raise ValueError(f"kx at index {ik + 1} is off the requested uniform grid")
        if abs(iy_float - round(iy_float)) > tolerance * ny:
            raise ValueError(f"ky at index {ik + 1} is off the requested uniform grid")
        if index[ix, iy] >= 0:
            raise ValueError("duplicate uniform-grid point")
        index[ix, iy] = ik
    if np.any(index < 0):
        raise ValueError("uniform grid is incomplete")
    return Grid(nx=nx, ny=ny, index=index, offset=offset)


def common_g_indices(g_left: np.ndarray, g_right: np.ndarray, right_shift: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Align left G with right extended-zone G = stored G - shift."""

    right_lookup = {
        tuple(int(x) for x in (g - right_shift)): i for i, g in enumerate(g_right)
    }
    left_indices: list[int] = []
    right_indices: list[int] = []
    for i, g in enumerate(g_left):
        j = right_lookup.get(tuple(int(x) for x in g))
        if j is not None:
            left_indices.append(i)
            right_indices.append(j)
    return np.asarray(left_indices), np.asarray(right_indices)


def link_matrix(
    wavecar: Wavecar,
    ik_left: int,
    ik_right: int,
    bands: Sequence[int],
    desired_right: np.ndarray,
) -> tuple[np.ndarray, float, float]:
    shift = np.rint(desired_right - wavecar.kpoints[ik_right]).astype(np.int64)
    residual = desired_right - wavecar.kpoints[ik_right] - shift
    if np.max(np.abs(residual)) > 1e-6:
        raise ValueError("neighbor is not related by a reciprocal-lattice shift")
    g_left = wavecar.g_vectors(ik_left)
    g_right = wavecar.g_vectors(ik_right)
    il, ir = common_g_indices(g_left, g_right, shift)
    left = wavecar.coefficients(ik_left, bands)
    right = wavecar.coefficients(ik_right, bands)
    overlap = np.zeros((len(bands), len(bands)), dtype=np.complex128)
    for component in range(wavecar.spinor_components):
        # WAVECAR coefficients are complex64, but link overlaps can contain
        # severe cancellation across plane waves.  Cast both operands before
        # matrix multiplication so products and their reduction are evaluated
        # in complex128 rather than merely accumulated into a complex128 array.
        left_component = left[:, component, il].astype(np.complex128, copy=False)
        right_component = right[:, component, ir].astype(np.complex128, copy=False)
        overlap += left_component.conj() @ right_component.T
    coverage_left = il.size / g_left.shape[0]
    coverage_right = ir.size / g_right.shape[0]
    return overlap, coverage_left, coverage_right


@dataclass
class LinkSet:
    phase: np.ndarray
    logabs: np.ndarray
    min_singular: np.ndarray
    max_singular: np.ndarray
    coverage_left: np.ndarray
    coverage_right: np.ndarray


def compute_links(wavecar: Wavecar, grid: Grid, bands: Sequence[int], axis: int) -> LinkSet:
    shape = (grid.nx, grid.ny)
    phase = np.empty(shape)
    logabs = np.empty(shape)
    min_singular = np.empty(shape)
    max_singular = np.empty(shape)
    coverage_left = np.empty(shape)
    coverage_right = np.empty(shape)
    step = np.zeros(3)
    step[axis] = 1.0 / (grid.nx if axis == 0 else grid.ny)
    for ix in range(grid.nx):
        for iy in range(grid.ny):
            ik_left = int(grid.index[ix, iy])
            jx = (ix + 1) % grid.nx if axis == 0 else ix
            jy = (iy + 1) % grid.ny if axis == 1 else iy
            ik_right = int(grid.index[jx, jy])
            desired_right = wavecar.kpoints[ik_left] + step
            matrix, cov_l, cov_r = link_matrix(
                wavecar, ik_left, ik_right, bands, desired_right
            )
            sign, lad = np.linalg.slogdet(matrix)
            if sign == 0:
                phase[ix, iy] = np.nan
                logabs[ix, iy] = -np.inf
            else:
                phase[ix, iy] = float(np.angle(sign))
                logabs[ix, iy] = float(lad)
            singular = np.linalg.svd(matrix, compute_uv=False)
            min_singular[ix, iy] = float(singular[-1])
            max_singular[ix, iy] = float(singular[0])
            coverage_left[ix, iy] = cov_l
            coverage_right[ix, iy] = cov_r
    return LinkSet(
        phase=phase,
        logabs=logabs,
        min_singular=min_singular,
        max_singular=max_singular,
        coverage_left=coverage_left,
        coverage_right=coverage_right,
    )


def wrap_angle(value: np.ndarray | float) -> np.ndarray | float:
    return np.angle(np.exp(1j * value))


def plaquette_flux(xlinks: LinkSet, ylinks: LinkSet, grid: Grid) -> np.ndarray:
    flux = np.empty((grid.nx, grid.ny))
    for ix in range(grid.nx):
        for iy in range(grid.ny):
            ix1, iy1 = (ix + 1) % grid.nx, (iy + 1) % grid.ny
            loop_arg = (
                xlinks.phase[ix, iy]
                + ylinks.phase[ix1, iy]
                - xlinks.phase[ix, iy1]
                - ylinks.phase[ix, iy]
            )
            flux[ix, iy] = -float(wrap_angle(loop_arg))
    return flux


def link_loop_diagnostics(links_x: LinkSet, links_y: LinkSet, ix: int, iy: int, grid: Grid) -> dict[str, float]:
    ix1, iy1 = (ix + 1) % grid.nx, (iy + 1) % grid.ny
    logabs = (
        links_x.logabs[ix, iy],
        links_y.logabs[ix1, iy],
        links_x.logabs[ix, iy1],
        links_y.logabs[ix, iy],
    )
    min_sv = (
        links_x.min_singular[ix, iy],
        links_y.min_singular[ix1, iy],
        links_x.min_singular[ix, iy1],
        links_y.min_singular[ix, iy],
    )
    max_sv = (
        links_x.max_singular[ix, iy],
        links_y.max_singular[ix1, iy],
        links_x.max_singular[ix, iy1],
        links_y.max_singular[ix, iy],
    )
    coverage = (
        links_x.coverage_left[ix, iy], links_x.coverage_right[ix, iy],
        links_y.coverage_left[ix1, iy], links_y.coverage_right[ix1, iy],
        links_x.coverage_left[ix, iy1], links_x.coverage_right[ix, iy1],
        links_y.coverage_left[ix, iy], links_y.coverage_right[ix, iy],
    )
    return {
        **{f"link{j + 1}_logabs_det": float(logabs[j]) for j in range(4)},
        **{
            f"link{j + 1}_abs_det": (
                float(math.exp(logabs[j])) if logabs[j] > math.log(np.finfo(float).tiny) else 0.0
            )
            for j in range(4)
        },
        **{f"link{j + 1}_min_sv": float(min_sv[j]) for j in range(4)},
        **{f"link{j + 1}_max_sv": float(max_sv[j]) for j in range(4)},
        "min_link_sv": float(min(min_sv)),
        "min_pw_coverage": float(min(coverage)),
    }


def compute_band_map(wavecar: Wavecar, grid: Grid, bands: Sequence[int]) -> tuple[np.ndarray, LinkSet, LinkSet]:
    xlinks = compute_links(wavecar, grid, bands, axis=0)
    ylinks = compute_links(wavecar, grid, bands, axis=1)
    return plaquette_flux(xlinks, ylinks, grid), xlinks, ylinks


def _empty_link_set(shape: tuple[int, int]) -> LinkSet:
    """Allocate one link container used by the cumulative-band batch path."""

    return LinkSet(
        phase=np.empty(shape),
        logabs=np.empty(shape),
        min_singular=np.empty(shape),
        max_singular=np.empty(shape),
        coverage_left=np.empty(shape),
        coverage_right=np.empty(shape),
    )


def compute_cumulative_links(
    wavecar: Wavecar,
    grid: Grid,
    max_band: int,
    axis: int,
) -> dict[int, LinkSet]:
    """Compute links for every leading subspace 1:n in one WAVECAR pass.

    Reading the overlap matrix for bands 1 through ``max_band`` once per link
    avoids rereading the same coefficient records for every cumulative
    subspace.  The leading principal block is the overlap matrix for 1:n.
    Internal degeneracies are therefore kept inside the determinant subspace;
    isolation from band n+1 is checked later wherever that subspace is used.
    """

    if axis not in (0, 1):
        raise ValueError("cumulative link axis must be 0 or 1")
    if not (1 <= max_band <= wavecar.header.nbands):
        raise ValueError(
            f"cumulative max band {max_band} is outside 1:{wavecar.header.nbands}"
        )
    shape = (grid.nx, grid.ny)
    result = {band: _empty_link_set(shape) for band in range(1, max_band + 1)}
    bands = list(range(1, max_band + 1))
    step = np.zeros(3)
    step[axis] = 1.0 / (grid.nx if axis == 0 else grid.ny)
    for ix in range(grid.nx):
        for iy in range(grid.ny):
            ik_left = int(grid.index[ix, iy])
            jx = (ix + 1) % grid.nx if axis == 0 else ix
            jy = (iy + 1) % grid.ny if axis == 1 else iy
            ik_right = int(grid.index[jx, jy])
            desired_right = wavecar.kpoints[ik_left] + step
            matrix, coverage_left, coverage_right = link_matrix(
                wavecar, ik_left, ik_right, bands, desired_right
            )
            for band, links in result.items():
                block = matrix[:band, :band]
                sign, logabs = np.linalg.slogdet(block)
                if sign == 0:
                    links.phase[ix, iy] = np.nan
                    links.logabs[ix, iy] = -np.inf
                else:
                    links.phase[ix, iy] = float(np.angle(sign))
                    links.logabs[ix, iy] = float(logabs)
                singular = np.linalg.svd(block, compute_uv=False)
                links.min_singular[ix, iy] = float(singular[-1])
                links.max_singular[ix, iy] = float(singular[0])
                links.coverage_left[ix, iy] = coverage_left
                links.coverage_right[ix, iy] = coverage_right
    return result


def compute_cumulative_band_maps(
    wavecar: Wavecar,
    grid: Grid,
    max_band: int,
) -> dict[int, tuple[np.ndarray, LinkSet, LinkSet]]:
    """Return Fukui maps for every gauge-covariant subspace 1:n, n<=N."""

    xlinks = compute_cumulative_links(wavecar, grid, max_band, axis=0)
    ylinks = compute_cumulative_links(wavecar, grid, max_band, axis=1)
    return {
        band: (plaquette_flux(xlinks[band], ylinks[band], grid),
               xlinks[band], ylinks[band])
        for band in range(1, max_band + 1)
    }


def vertex_indices(grid: Grid, ix: int, iy: int) -> tuple[int, int, int, int]:
    ix1, iy1 = (ix + 1) % grid.nx, (iy + 1) % grid.ny
    return (
        int(grid.index[ix, iy]),
        int(grid.index[ix1, iy]),
        int(grid.index[ix1, iy1]),
        int(grid.index[ix, iy1]),
    )


def write_table(
    output: Path,
    wavecar: Wavecar,
    grid: Grid,
    label: str,
    bands: Sequence[int],
    flux: np.ndarray,
    xlinks: LinkSet,
    ylinks: LinkSet,
    energy_band: int,
) -> None:
    area = float(np.linalg.norm(np.cross(*wavecar.header.reciprocal[:2]))) / (grid.nx * grid.ny)
    rows: list[dict[str, int | float | str]] = []
    for iy in range(grid.ny):
        for ix in range(grid.nx):
            vertices = vertex_indices(grid, ix, iy)
            base = wavecar.kpoints[vertices[0]]
            center = base + np.array((0.5 / grid.nx, 0.5 / grid.ny, 0.0))
            energies = [float(wavecar.energies[ik, energy_band - 1]) for ik in vertices]
            gap_below = [
                float(wavecar.energies[ik, energy_band - 1] - wavecar.energies[ik, energy_band - 2])
                if energy_band > 1 else math.nan
                for ik in vertices
            ]
            gap_above = [
                float(wavecar.energies[ik, energy_band] - wavecar.energies[ik, energy_band - 1])
                if energy_band < wavecar.header.nbands else math.nan
                for ik in vertices
            ]
            energy_above = [
                float(wavecar.energies[ik, energy_band])
                if energy_band < wavecar.header.nbands else math.nan
                for ik in vertices
            ]
            subspace_edge = bands[-1]
            subspace_gap_above = [
                float(
                    wavecar.energies[ik, subspace_edge]
                    - wavecar.energies[ik, subspace_edge - 1]
                )
                if subspace_edge < wavecar.header.nbands else math.nan
                for ik in vertices
            ]
            row: dict[str, int | float | str] = {
                "map": label,
                "band_first": int(bands[0]),
                "band_last": int(bands[-1]),
                "band": energy_band,
                "cell_id": iy * grid.nx + ix + 1,
                "ix": ix,
                "iy": iy,
                "k0_index": vertices[0] + 1,
                "k1_index": vertices[1] + 1,
                "k2_index": vertices[2] + 1,
                "k3_index": vertices[3] + 1,
                "k0": vertices[0] + 1,
                "k1": vertices[1] + 1,
                "k2": vertices[2] + 1,
                "k3": vertices[3] + 1,
                "k0_frac_x": float(base[0]),
                "k0_frac_y": float(base[1]),
                "k0_frac_z": float(base[2]),
                "center_frac_x": float(center[0]),
                "center_frac_y": float(center[1]),
                "center_frac_z": float(center[2]),
                "q_center_frac_x": float(center[0]),
                "q_center_frac_y": float(center[1]),
                "q_center_frac_z": float(center[2]),
                "q1": float(center[0]),
                "q2": float(center[1]),
                "q3": float(center[2]),
                "energy1_ev": energies[0],
                "energy2_ev": energies[1],
                "energy3_ev": energies[2],
                "energy4_ev": energies[3],
                "energy_above1_ev": energy_above[0],
                "energy_above2_ev": energy_above[1],
                "energy_above3_ev": energy_above[2],
                "energy_above4_ev": energy_above[3],
                "e0_eV": energies[0],
                "e1_eV": energies[1],
                "e2_eV": energies[2],
                "e3_eV": energies[3],
                "min_gap_below_ev": min(gap_below),
                "min_gap_above_ev": min(gap_above),
                "gap_below1_ev": gap_below[0],
                "gap_below2_ev": gap_below[1],
                "gap_below3_ev": gap_below[2],
                "gap_below4_ev": gap_below[3],
                "gap_above1_ev": gap_above[0],
                "gap_above2_ev": gap_above[1],
                "gap_above3_ev": gap_above[2],
                "gap_above4_ev": gap_above[3],
                "subspace_min_gap_above_ev": min(subspace_gap_above),
                "fukui_flux_rad": float(flux[ix, iy]),
                "phi_rad": float(flux[ix, iy]),
                "plaquette_area_inv_ang2": area,
                "dS_A-2": area,
                "berry_curvature_ang2": float(flux[ix, iy] / area),
                "omega_ang2": float(flux[ix, iy] / area),
                "omega_A2": float(flux[ix, iy] / area),
                "neighbor_gap_below_ev": min(gap_below),
                "neighbor_gap_above_ev": min(gap_above),
            }
            row.update(link_loop_diagnostics(xlinks, ylinks, ix, iy, grid))
            row["link_quality_min_sv"] = row["min_link_sv"]
            rows.append(row)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def map_summary(
    label: str,
    bands: Sequence[int],
    flux: np.ndarray,
    xlinks: LinkSet,
    ylinks: LinkSet,
    min_pw_coverage_required: float,
) -> dict[str, object]:
    min_coverage = float(
        min(
            np.min(xlinks.coverage_left),
            np.min(xlinks.coverage_right),
            np.min(ylinks.coverage_left),
            np.min(ylinks.coverage_right),
        )
    )
    return {
        "label": label,
        "bands": [int(bands[0]), int(bands[-1])],
        "chern": float(np.sum(flux) / (2.0 * np.pi)),
        "flux_min_rad": float(np.min(flux)),
        "flux_max_rad": float(np.max(flux)),
        "max_abs_flux_rad": float(np.max(np.abs(flux))),
        "min_link_singular_value": float(
            min(np.min(xlinks.min_singular), np.min(ylinks.min_singular))
        ),
        "max_link_condition_number": float(
            max(
                np.max(xlinks.max_singular / xlinks.min_singular),
                np.max(ylinks.max_singular / ylinks.min_singular),
            )
        ),
        "min_plane_wave_coverage": min_coverage,
        "plane_wave_coverage_required": min_pw_coverage_required,
        "plane_wave_coverage_pass": min_coverage >= min_pw_coverage_required,
    }


def cell_min_link_singular_value(xlinks: LinkSet, ylinks: LinkSet, grid: Grid) -> np.ndarray:
    result = np.empty((grid.nx, grid.ny))
    for ix in range(grid.nx):
        for iy in range(grid.ny):
            result[ix, iy] = link_loop_diagnostics(
                xlinks, ylinks, ix, iy, grid
            )["min_link_sv"]
    return result


def periodic_reciprocal_distance(
    qpoint: np.ndarray, center: np.ndarray, reciprocal: np.ndarray
) -> float:
    """Shortest in-plane distance over every reciprocal-lattice image.

    A fixed image stencil is not sufficient for a skew or non-reduced basis:
    the closest image can have integer coefficients well outside [-1, 1].  We
    Gauss-reduce the two lattice vectors and solve the resulting two-dimensional
    closest-vector problem exactly with the four-candidate reduced-basis rule.
    """

    base_delta = np.asarray(qpoint, dtype=float) - np.asarray(center, dtype=float)
    base_delta[:2] -= np.rint(base_delta[:2])
    lattice_basis = np.asarray(reciprocal[:2], dtype=float)
    reduced_basis = gauss_reduce_2d_lattice_basis(lattice_basis)
    target = -(base_delta @ reciprocal)
    return float(np.linalg.norm(closest_2d_lattice_residual(target, reduced_basis)))


def gauss_reduce_2d_lattice_basis(basis: np.ndarray) -> np.ndarray:
    """Return a Gauss-reduced basis for the same rank-two lattice."""

    reduced = np.asarray(basis, dtype=float).copy()
    if reduced.shape[0] != 2 or not np.all(np.isfinite(reduced)):
        raise ValueError("in-plane reciprocal basis must be finite and rank two")
    while True:
        norm0 = float(np.dot(reduced[0], reduced[0]))
        norm1 = float(np.dot(reduced[1], reduced[1]))
        if (
            not math.isfinite(norm0)
            or not math.isfinite(norm1)
            or norm0 <= 0.0
            or norm1 <= 0.0
        ):
            raise ValueError("in-plane reciprocal basis must be finite and rank two")
        if norm1 < norm0:
            reduced[[0, 1]] = reduced[[1, 0]]
            norm0, norm1 = norm1, norm0
        coefficient = int(np.rint(np.dot(reduced[0], reduced[1]) / norm0))
        if coefficient == 0:
            singular = np.linalg.svd(reduced, compute_uv=False)
            rank_tolerance = (
                np.finfo(float).eps * max(reduced.shape) * singular[0]
            )
            if singular[-1] <= rank_tolerance:
                raise ValueError(
                    "in-plane reciprocal basis must be finite and rank two"
                )
            return reduced
        updated = reduced[1] - coefficient * reduced[0]
        updated_norm = float(np.dot(updated, updated))
        if (
            not np.all(np.isfinite(updated))
            or not math.isfinite(updated_norm)
            or updated_norm >= norm1
        ):
            raise RuntimeError(
                "Gauss reduction lost numerical progress for reciprocal basis"
            )
        reduced[1] = updated


def closest_2d_lattice_residual(target: np.ndarray, basis: np.ndarray) -> np.ndarray:
    """Exact nearest-lattice residual for a Gauss-reduced rank-two basis.

    For a reduced basis u,v, completing the quadratic form shows that only
    floor/ceil of the continuous v coordinate, followed by floor/ceil of the
    conditional u coordinate, can contain the integer optimum.  Thus at most
    four full Cartesian residuals are evaluated, independent of skew or scale.
    """

    basis = np.asarray(basis, dtype=float)
    matrix = basis.T
    target = np.asarray(target, dtype=float)
    # The reduced basis may be related to the input by a large integer shear.
    # Shift the target into its nearby parallelogram before QR to avoid losing
    # precision through cancellation of large lattice coordinates.
    target_coordinates = np.linalg.lstsq(matrix, target, rcond=None)[0]
    target = target - np.rint(target_coordinates) @ basis
    vector_u, vector_v = basis
    norm_u = float(np.dot(vector_u, vector_u))
    norm_v = float(np.dot(vector_v, vector_v))
    cross_term = float(np.dot(vector_u, vector_v))
    target_u = float(np.dot(vector_u, target))
    target_v = float(np.dot(vector_v, target))
    determinant = norm_u * norm_v - cross_term * cross_term
    if not math.isfinite(determinant) or determinant <= 0.0:
        raise ValueError("in-plane reciprocal basis must be rank two")
    second_real = (
        norm_u * target_v - cross_term * target_u
    ) / determinant
    second_candidates = {math.floor(second_real), math.ceil(second_real)}
    best_squared = math.inf
    best_residual: np.ndarray | None = None
    for second in second_candidates:
        first_real = (target_u - cross_term * second) / norm_u
        for first in {math.floor(first_real), math.ceil(first_real)}:
            residual = first * vector_u + second * vector_v - target
            squared = float(np.dot(residual, residual))
            if squared < best_squared:
                best_squared = squared
                best_residual = residual
    if best_residual is None:
        raise RuntimeError("closest-lattice candidate enumeration produced no result")
    return best_residual


def shortest_periodic_fractional_displacement(
    qpoint: np.ndarray, center: np.ndarray, reciprocal: np.ndarray
) -> np.ndarray:
    """Return q-center in the globally shortest in-plane periodic image."""

    base_delta = np.asarray(qpoint, dtype=float) - np.asarray(center, dtype=float)
    base_delta[:2] -= np.rint(base_delta[:2])
    reduced_basis = gauss_reduce_2d_lattice_basis(
        np.asarray(reciprocal[:2], dtype=float)
    )
    residual = closest_2d_lattice_residual(
        -(base_delta @ reciprocal), reduced_basis
    )
    try:
        return np.linalg.solve(np.asarray(reciprocal, dtype=float).T, residual)
    except np.linalg.LinAlgError as error:
        raise ValueError("reciprocal basis must be rank three") from error


def validate_valley_definition(
    k_center: np.ndarray,
    kp_center: np.ndarray,
    reciprocal: np.ndarray,
    radius_inv_angstrom: float | None,
) -> float:
    """Validate valley inputs and return their shortest periodic separation."""

    if (
        np.asarray(k_center).shape != (3,)
        or np.asarray(kp_center).shape != (3,)
        or not np.all(np.isfinite(k_center))
        or not np.all(np.isfinite(kp_center))
    ):
        raise ValueError("K and K' centers must be finite three-component vectors")
    if radius_inv_angstrom is not None:
        if not math.isfinite(radius_inv_angstrom):
            raise ValueError("--valley-radius must be finite")
        if radius_inv_angstrom <= 0.0:
            raise ValueError("--valley-radius must be positive")
    valley_separation = periodic_reciprocal_distance(
        k_center, kp_center, reciprocal
    )
    if valley_separation <= VALLEY_DISTANCE_TIE_ATOL:
        raise ValueError("K and K' centers must be periodically distinct")
    if (
        radius_inv_angstrom is not None
        and 2.0 * radius_inv_angstrom >= valley_separation
    ):
        raise ValueError(
            "K and K' radius masks overlap analytically: 2*radius must be "
            f"smaller than their periodic separation ({valley_separation:.9g} A^-1)"
        )
    return valley_separation


def make_valley_partition(
    wavecar: Wavecar,
    grid: Grid,
    k_center: np.ndarray,
    kp_center: np.ndarray,
    radius_inv_angstrom: float | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    validate_valley_definition(
        k_center, kp_center, wavecar.header.reciprocal, radius_inv_angstrom
    )
    mask_k = np.zeros((grid.nx, grid.ny), dtype=bool)
    mask_kp = np.zeros_like(mask_k)
    mask_outside = np.zeros_like(mask_k)
    distance_k = np.empty_like(mask_k, dtype=float)
    distance_kp = np.empty_like(mask_k, dtype=float)
    for ix in range(grid.nx):
        for iy in range(grid.ny):
            ik = int(grid.index[ix, iy])
            center = wavecar.kpoints[ik] + np.asarray(
                (0.5 / grid.nx, 0.5 / grid.ny, 0.0)
            )
            dk = periodic_reciprocal_distance(
                center, k_center, wavecar.header.reciprocal
            )
            dkp = periodic_reciprocal_distance(
                center, kp_center, wavecar.header.reciprocal
            )
            distance_k[ix, iy], distance_kp[ix, iy] = dk, dkp
            if radius_inv_angstrom is None:
                if math.isclose(
                    dk,
                    dkp,
                    rel_tol=1.0e-12,
                    abs_tol=VALLEY_DISTANCE_TIE_ATOL,
                ):
                    mask_outside[ix, iy] = True
                elif dk < dkp:
                    mask_k[ix, iy] = True
                else:
                    mask_kp[ix, iy] = True
            else:
                in_k, in_kp = dk <= radius_inv_angstrom, dkp <= radius_inv_angstrom
                if in_k and in_kp:
                    raise ValueError(
                        "K and K' radius masks overlap; reduce --valley-radius"
                    )
                mask_k[ix, iy], mask_kp[ix, iy] = in_k, in_kp
                mask_outside[ix, iy] = not (in_k or in_kp)
    if not np.any(mask_k) or not np.any(mask_kp):
        mode = (
            "--valley-radius"
            if radius_inv_angstrom is not None
            else "nearest-center Voronoi partition"
        )
        raise ValueError(
            f"{mode} selects no plaquette centers for K or K'; "
            "adjust the valley definition or use a denser mesh"
        )
    return mask_k, mask_kp, mask_outside, distance_k, distance_kp


def parse_fractional_vector(text: str) -> np.ndarray:
    try:
        values = np.asarray([float(item) for item in text.split(",")], dtype=float)
    except ValueError as error:
        raise argparse.ArgumentTypeError("expected q1,q2,q3") from error
    if values.shape != (3,):
        raise argparse.ArgumentTypeError("expected q1,q2,q3")
    if not np.all(np.isfinite(values)):
        raise argparse.ArgumentTypeError("q1,q2,q3 must all be finite")
    return values


def cumulative_t0_effective_phi(
    phi_v: np.ndarray,
    phi_v1: np.ndarray,
    phi_v2: np.ndarray,
    energy_v1_vertices: np.ndarray,
    energy_v2_vertices: np.ndarray,
    mu_ev: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return Sawahata-style phase and the 0/1/2 cumulative-state counts."""

    if not math.isfinite(float(mu_ev)):
        raise ValueError("chemical potential must be finite")
    if (
        energy_v1_vertices.shape != energy_v2_vertices.shape
        or energy_v1_vertices.ndim != 3
        or energy_v1_vertices.shape[2] != 4
        or energy_v1_vertices.shape[:2] != phi_v.shape
        or phi_v1.shape != phi_v.shape
        or phi_v2.shape != phi_v.shape
    ):
        raise ValueError("cumulative transport requires four vertices per plaquette")
    if not all(
        np.all(np.isfinite(values))
        for values in (
            phi_v,
            phi_v1,
            phi_v2,
            energy_v1_vertices,
            energy_v2_vertices,
        )
    ):
        raise ValueError("cumulative transport phase and energy arrays must be finite")

    state0 = mu_ev < energy_v1_vertices
    state1 = (mu_ev >= energy_v1_vertices) & (mu_ev < energy_v2_vertices)
    state2 = mu_ev >= energy_v2_vertices
    count0 = np.sum(state0, axis=2)
    count1 = np.sum(state1, axis=2)
    count2 = np.sum(state2, axis=2)
    if not np.all(count0 + count1 + count2 == 4):
        raise RuntimeError(
            "cumulative transport state-count invariant failed: every cell must "
            "assign exactly four vertices"
        )
    effective_phi = (
        count0 * phi_v + count1 * phi_v1 + count2 * phi_v2
    ) / 4.0
    return effective_phi, count0, count1, count2


def full_cumulative_t0_effective_phi(
    phi_by_occupied_count: np.ndarray,
    band_energy_vertices: np.ndarray,
    mu_ev: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Average cumulative-subspace phases selected at all four vertices.

    ``phi_by_occupied_count[n]`` is the plaquette phase of bands 1:n and
    index zero is the exact empty-subspace phase.  ``band_energy_vertices``
    has shape ``(max_band, nx, ny, 4)``.  Counting all energies ``E <= mu``
    makes an exactly degenerate group enter together and never forms an
    arbitrary band-by-band Berry-curvature sum inside that group.
    """

    if not math.isfinite(float(mu_ev)):
        raise ValueError("chemical potential must be finite")
    if (
        phi_by_occupied_count.ndim != 3
        or band_energy_vertices.ndim != 4
        or band_energy_vertices.shape[0] + 1 != phi_by_occupied_count.shape[0]
        or band_energy_vertices.shape[1:3] != phi_by_occupied_count.shape[1:]
        or band_energy_vertices.shape[3] != 4
    ):
        raise ValueError(
            "full cumulative transport requires phases 0:N and four energy "
            "vertices for bands 1:N"
        )
    if not np.all(np.isfinite(band_energy_vertices)):
        raise ValueError("full cumulative transport band energies must be finite")
    if not np.all(phi_by_occupied_count[0] == 0.0):
        raise ValueError("empty occupied subspace must have exactly zero phase")
    if np.any(np.diff(band_energy_vertices, axis=0) < 0.0):
        raise ValueError("band energies must be nondecreasing at every vertex")

    occupied_count = np.sum(
        band_energy_vertices <= float(mu_ev), axis=0, dtype=np.int64
    )
    nx, ny = phi_by_occupied_count.shape[1:]
    ix = np.arange(nx)[:, None]
    iy = np.arange(ny)[None, :]
    effective_phi = np.zeros((nx, ny), dtype=float)
    for vertex in range(4):
        # A nonfinite map is permitted to exist outside the requested
        # occupation window. If selected, the writer's active-cell quality
        # gate records and refuses it instead of treating NaN as a passing
        # comparison.
        effective_phi += phi_by_occupied_count[
            occupied_count[:, :, vertex], ix, iy
        ]
    effective_phi /= 4.0
    return effective_phi, occupied_count


def remove_stale_transport_outputs(output_dir: Path) -> None:
    """Remove only transport artifacts that could be mistaken for this run."""

    remove_planned_outputs(output_dir, TRANSPORT_OUTPUT_NAMES)


def remove_planned_outputs(output_dir: Path, output_names: Iterable[str]) -> None:
    """Remove the exact fixed outputs planned for one invocation, and no others."""

    output_dir.mkdir(parents=True, exist_ok=True)
    for name in output_names:
        (output_dir / name).unlink(missing_ok=True)


def planned_output_names(
    maps: Sequence[tuple[str, Sequence[int]]],
    plot: bool,
    transport_t0: bool,
    transport_full_t0: bool = False,
) -> set[str]:
    """Return every fixed artifact name that this invocation can overwrite."""

    names = {"diagnostics.json"}
    names.update(f"fukui_{label}.csv" for label, _ in maps)
    if plot:
        names.update(
            ("wavecar_fukui_kresolved.png", "wavecar_fukui_line_K_Kp.png")
        )
    if transport_t0:
        # The sigma plot is removed as stale even when this run omits --plot.
        names.update(TRANSPORT_OUTPUT_NAMES)
    if transport_full_t0:
        # Likewise, never leave an old full-window curve after a refused run.
        names.update(FULL_TRANSPORT_OUTPUT_NAMES)
    return names


def validate_input_output_collision(
    input_path: Path, output_dir: Path, output_names: Iterable[str]
) -> None:
    """Refuse any run whose resolved WAVECAR path is a planned output."""

    resolved_input = input_path.resolve(strict=False)
    for name in output_names:
        candidate = (output_dir / name).resolve(strict=False)
        if candidate == resolved_input:
            raise ValueError(
                f"WAVECAR input resolves to planned output {candidate}; refusing "
                "to overwrite or remove the input"
            )


def write_t0_transport(
    output_dir: Path,
    wavecar: Wavecar,
    grid: Grid,
    map_results: dict[str, tuple[list[int], np.ndarray, LinkSet, LinkSet]],
    energy_band: int,
    mu_min: float,
    mu_max: float,
    mu_num: int,
    valley_k: np.ndarray,
    valley_kp: np.ndarray,
    valley_radius: float | None,
    min_link_sv: float,
    min_neighbor_gap_ev: float,
    max_abs_phi: float,
    allow_invalid: bool,
) -> None:
    wavecar_path = getattr(wavecar, "path", None)
    if wavecar_path is not None:
        validate_input_output_collision(
            Path(wavecar_path), output_dir, TRANSPORT_OUTPUT_NAMES
        )
    remove_stale_transport_outputs(output_dir)
    if energy_band < 2:
        raise ValueError("target energy band must be at least 2")
    if not (math.isfinite(mu_min) and math.isfinite(mu_max)):
        raise ValueError("mu_min and mu_max must be finite")
    if mu_num < 2 or mu_max <= mu_min:
        raise ValueError("transport requires mu_num>=2 and mu_max>mu_min")
    for option, value in (
        ("min_link_sv", min_link_sv),
        ("min_neighbor_gap_ev", min_neighbor_gap_ev),
        ("max_abs_phi", max_abs_phi),
    ):
        if not math.isfinite(value):
            raise ValueError(f"{option} must be finite")
    if min_link_sv < 0.0 or min_neighbor_gap_ev < 0.0:
        raise ValueError("link and neighbor-gap thresholds must be nonnegative")
    if not (0.0 < max_abs_phi < np.pi):
        raise ValueError("max_abs_phi must be in (0, pi)")
    if energy_band + 2 > wavecar.header.nbands:
        raise ValueError("T=0 cumulative transport needs target and target+1 bands")
    if not np.all(
        np.isfinite(
            wavecar.energies[:, energy_band - 2 : energy_band + 2]
        )
    ):
        raise ValueError("transport energy window contains nonfinite band energies")
    valence_band_max = float(np.max(wavecar.energies[:, energy_band - 2]))
    if mu_min <= valence_band_max:
        raise ValueError(
            f"mu_min={mu_min:.9g} eV is not above the band {energy_band - 1} "
            f"maximum ({valence_band_max:.9g} eV); the fixed valence baseline "
            "would contain holes"
        )
    next_unrepresented_min = float(np.min(wavecar.energies[:, energy_band + 1]))
    if mu_max >= next_unrepresented_min:
        raise ValueError(
            f"mu_max={mu_max:.9g} eV reaches band {energy_band + 2} "
            f"(minimum {next_unrepresented_min:.9g} eV); add another cumulative "
            "subspace before extending the transport window"
        )

    by_range: dict[tuple[int, int], tuple[str, np.ndarray, LinkSet, LinkSet]] = {}
    for label, (bands, flux, xlinks, ylinks) in map_results.items():
        by_range[(bands[0], bands[-1])] = (label, flux, xlinks, ylinks)
    required = ((1, energy_band - 1), (1, energy_band), (1, energy_band + 1))
    missing = [item for item in required if item not in by_range]
    if missing:
        raise ValueError(f"transport requires cumulative maps {required}; missing {missing}")

    label_v, phi_v, x_v, y_v = by_range[required[0]]
    label_v1, phi_v1, x_v1, y_v1 = by_range[required[1]]
    label_v2, phi_v2, x_v2, y_v2 = by_range[required[2]]
    link_v1 = cell_min_link_singular_value(x_v1, y_v1, grid)
    link_v2 = cell_min_link_singular_value(x_v2, y_v2, grid)
    link_v = cell_min_link_singular_value(x_v, y_v, grid)
    baseline_min_link = float(np.min(link_v))
    baseline_min_gap = float(
        np.min(
            wavecar.energies[:, energy_band - 1]
            - wavecar.energies[:, energy_band - 2]
        )
    )
    baseline_max_abs_phi = float(np.max(np.abs(phi_v)))
    baseline_quality_pass = bool(
        baseline_min_link >= min_link_sv
        and baseline_min_gap >= min_neighbor_gap_ev
        and baseline_max_abs_phi <= max_abs_phi
    )

    e1 = np.empty((grid.nx, grid.ny, 4))
    e2 = np.empty_like(e1)
    e3 = np.empty_like(e1)
    for ix in range(grid.nx):
        for iy in range(grid.ny):
            vertices = vertex_indices(grid, ix, iy)
            e1[ix, iy] = [wavecar.energies[ik, energy_band - 1] for ik in vertices]
            e2[ix, iy] = [wavecar.energies[ik, energy_band] for ik in vertices]
            if energy_band + 1 < wavecar.header.nbands:
                e3[ix, iy] = [wavecar.energies[ik, energy_band + 1] for ik in vertices]
            else:
                e3[ix, iy] = math.nan

    mask_k, mask_kp, mask_outside, distance_k, distance_kp = make_valley_partition(
        wavecar, grid, valley_k, valley_kp, valley_radius
    )
    tie_cell_count = (
        int(np.count_nonzero(mask_outside)) if valley_radius is None else 0
    )
    mask_cell_count_sum = int(
        np.count_nonzero(mask_k)
        + np.count_nonzero(mask_kp)
        + np.count_nonzero(mask_outside)
    )
    if (
        mask_cell_count_sum != grid.nx * grid.ny
        or np.any(mask_k & mask_kp)
        or np.any(mask_k & mask_outside)
        or np.any(mask_kp & mask_outside)
    ):
        raise RuntimeError("valley masks do not form a complete disjoint partition")
    mask_partition_complete = bool(
        np.all(
            mask_k.astype(np.int8)
            + mask_kp.astype(np.int8)
            + mask_outside.astype(np.int8)
            == 1
        )
    )
    if not mask_partition_complete:
        raise RuntimeError("valley masks are not exclusive and exhaustive per cell")
    two_pi = 2.0 * np.pi
    baseline_chern = float(np.sum(phi_v) / two_pi)
    rows: list[dict[str, float | int | str]] = []
    invalid_rows: list[dict[str, object]] = []
    for mu in np.linspace(mu_min, mu_max, mu_num):
        effective_phi, count0, count1, count2 = cumulative_t0_effective_phi(
            phi_v, phi_v1, phi_v2, e1, e2, float(mu)
        )
        delta_phi = effective_phi - phi_v
        active = (count1 + count2) > 0
        used_v1 = count1 > 0
        used_v2 = count2 > 0

        quality_link = np.full((grid.nx, grid.ny), math.inf)
        quality_link[used_v1] = np.minimum(quality_link[used_v1], link_v1[used_v1])
        quality_link[used_v2] = np.minimum(quality_link[used_v2], link_v2[used_v2])
        cell_gap_v1 = np.min(e2 - e1, axis=2)
        cell_gap_v2 = np.min(e3 - e2, axis=2)
        used_cell_gaps: list[float] = []
        if np.any(used_v1):
            used_cell_gaps.extend(cell_gap_v1[used_v1].tolist())
        if np.any(used_v2):
            used_cell_gaps.extend(cell_gap_v2[used_v2].tolist())
        worst_gap = min(used_cell_gaps) if used_cell_gaps else math.inf
        worst_link = float(np.min(quality_link[active])) if np.any(active) else math.inf
        used_phi_values: list[np.ndarray] = []
        if np.any(used_v1):
            used_phi_values.append(np.abs(phi_v1[used_v1]))
        if np.any(used_v2):
            used_phi_values.append(np.abs(phi_v2[used_v2]))
        worst_phi = (
            float(max(np.max(item) for item in used_phi_values))
            if used_phi_values else 0.0
        )
        quality_pass = bool(
            baseline_quality_pass
            and worst_link >= min_link_sv
            and worst_gap >= min_neighbor_gap_ev
            and worst_phi <= max_abs_phi
        )
        invalid_cells = np.argwhere(
            active
            & (
                (quality_link < min_link_sv)
                | (used_v1 & (cell_gap_v1 < min_neighbor_gap_ev))
                | (used_v2 & (cell_gap_v2 < min_neighbor_gap_ev))
                | ((used_v1) & (np.abs(phi_v1) > max_abs_phi))
                | ((used_v2) & (np.abs(phi_v2) > max_abs_phi))
            )
        )
        if not quality_pass:
            invalid_rows.append(
                {
                    "mu_ev": float(mu),
                    "worst_link_sv": worst_link,
                    "worst_neighbor_gap_ev": worst_gap,
                    "max_abs_phi_rad": worst_phi,
                    "baseline_quality_pass": baseline_quality_pass,
                    "invalid_cell_ids": [
                        int(iy * grid.nx + ix + 1) for ix, iy in invalid_cells[:20]
                    ],
                    "invalid_cell_count": int(len(invalid_cells)),
                }
            )

        delta_c = float(np.sum(delta_phi) / two_pi)
        delta_c_k = float(np.sum(delta_phi[mask_k]) / two_pi)
        delta_c_kp = float(np.sum(delta_phi[mask_kp]) / two_pi)
        delta_c_out = float(np.sum(delta_phi[mask_outside]) / two_pi)
        total_c = baseline_chern + delta_c
        rows.append(
            {
                "mu_eV": float(mu),
                "quality": "PASS" if quality_pass else "INVALID",
                "active_cell_count": int(np.count_nonzero(active)),
                "worst_active_link_sv": worst_link,
                "worst_active_neighbor_gap_eV": worst_gap,
                "max_active_abs_phi_rad": worst_phi,
                "chern_total": total_c,
                "sigma_total_e2_over_h": -total_c,
                "delta_chern": delta_c,
                "delta_sigma_e2_over_h": -delta_c,
                "delta_chern_K": delta_c_k,
                "delta_chern_Kp": delta_c_kp,
                "delta_chern_outside": delta_c_out,
                "delta_sigma_K_e2_over_h": -delta_c_k,
                "delta_sigma_Kp_e2_over_h": -delta_c_kp,
                "delta_sigma_outside_e2_over_h": -delta_c_out,
                "delta_sigma_valley_contrast_e2_over_h": -delta_c_k + delta_c_kp,
            }
        )

    diagnostics = {
        "method": "Sawahata-style T=0 cumulative-subspace plaquette average",
        "maps": {"V": label_v, "V_plus_1": label_v1, "V_plus_2": label_v2},
        "energy_band": energy_band,
        "mu_range_ev": [mu_min, mu_max, mu_num],
        "energy_window": {
            "valence_band_index": energy_band - 1,
            "valence_band_max_ev": valence_band_max,
            "mu_min_ev": mu_min,
            "mu_min_above_valence_max": mu_min > valence_band_max,
            "next_unrepresented_band_index": energy_band + 2,
            "next_unrepresented_band_min_ev": next_unrepresented_min,
            "mu_max_ev": mu_max,
            "mu_max_below_next_unrepresented_band_min": (
                mu_max < next_unrepresented_min
            ),
        },
        "baseline_chern": baseline_chern,
        "baseline_quality": {
            "pass": baseline_quality_pass,
            "minimum_link_singular_value": baseline_min_link,
            "minimum_gap_to_target_band_ev": baseline_min_gap,
            "maximum_abs_phi_rad": baseline_max_abs_phi,
            "valence_band_max_ev": valence_band_max,
            "mu_min_ev": mu_min,
            "mu_min_above_valence_max": mu_min > valence_band_max,
        },
        "valley_partition": {
            "K_fractional": valley_k.tolist(),
            "Kp_fractional": valley_kp.tolist(),
            "radius_inv_angstrom": valley_radius,
            "mode": "radius" if valley_radius is not None else "nearest-center Voronoi",
            "K_cell_count": int(np.count_nonzero(mask_k)),
            "Kp_cell_count": int(np.count_nonzero(mask_kp)),
            "outside_cell_count": int(np.count_nonzero(mask_outside)),
            "tie_cell_count": tie_cell_count,
            "mask_cell_count_sum": mask_cell_count_sum,
            "expected_cell_count": grid.nx * grid.ny,
            "mask_partition_complete": mask_partition_complete,
            "tie_tolerance_inv_angstrom": VALLEY_DISTANCE_TIE_ATOL,
            "min_center_distance_K_inv_angstrom": float(np.min(distance_k)),
            "min_center_distance_Kp_inv_angstrom": float(np.min(distance_kp)),
        },
        "quality_thresholds": {
            "min_link_singular_value": min_link_sv,
            "min_neighbor_gap_ev": min_neighbor_gap_ev,
            "max_abs_phi_rad": max_abs_phi,
            "active_gap_scope": (
                "all four vertices of each plaquette whose effective flux uses "
                "the V+1 or V+2 cumulative manifold"
            ),
        },
        "validated": baseline_quality_pass and not invalid_rows,
        "invalid_mu_count": len(invalid_rows),
        "invalid_mu": invalid_rows,
    }
    diagnostics_path = output_dir / "transport_t0_diagnostics.json"
    write_strict_json(diagnostics_path, diagnostics)
    if not baseline_quality_pass:
        raise RuntimeError(
            "fixed valence baseline failed its global link/gap/phase quality "
            f"gate (min link={baseline_min_link:.6g}, min gap="
            f"{baseline_min_gap:.6g} eV, max |phi|={baseline_max_abs_phi:.6g} rad); "
            f"{diagnostics_path.name} was written, but transport output is refused"
        )
    if invalid_rows and not allow_invalid:
        raise RuntimeError(
            "transport validation failed on active plaquettes; raw maps and "
            f"{diagnostics_path.name} were written, but transport_t0.csv was refused. "
            "Use --allow-invalid-transport only for explicitly unvalidated diagnostics."
        )
    with (output_dir / "transport_t0.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_full_t0_transport(
    output_dir: Path,
    wavecar: Wavecar,
    grid: Grid,
    max_band: int,
    mu_min: float,
    mu_max: float,
    mu_num: int,
    valley_k: np.ndarray,
    valley_kp: np.ndarray,
    valley_radius: float | None,
    min_link_sv: float,
    min_neighbor_gap_ev: float,
    min_pw_coverage: float,
    max_abs_phi: float,
    allow_invalid: bool,
) -> None:
    """Write absolute T=0 sigma(mu) from all cumulative subspaces 1:n.

    Band ``max_band + 1`` is a mandatory unoccupied sentinel.  The method can
    therefore span any represented valence bands, an insulating gap, and
    selected conduction bands without assigning Berry curvature to individual
    members of a degenerate group.
    """

    wavecar_path = getattr(wavecar, "path", None)
    if wavecar_path is not None:
        validate_input_output_collision(
            Path(wavecar_path), output_dir, FULL_TRANSPORT_OUTPUT_NAMES
        )
    remove_planned_outputs(output_dir, FULL_TRANSPORT_OUTPUT_NAMES)
    if not (1 <= max_band < wavecar.header.nbands):
        raise ValueError(
            "--transport-full-t0 MAX_BAND requires 1 <= MAX_BAND < NBANDS "
            "because MAX_BAND+1 is the upper occupation sentinel"
        )
    if not (math.isfinite(mu_min) and math.isfinite(mu_max)):
        raise ValueError("mu_min and mu_max must be finite")
    if mu_num < 2 or mu_max <= mu_min:
        raise ValueError("full transport requires mu_num>=2 and mu_max>mu_min")
    for option, value in (
        ("min_link_sv", min_link_sv),
        ("min_neighbor_gap_ev", min_neighbor_gap_ev),
        ("min_pw_coverage", min_pw_coverage),
        ("max_abs_phi", max_abs_phi),
    ):
        if not math.isfinite(value):
            raise ValueError(f"{option} must be finite")
    if min_link_sv < 0.0 or min_neighbor_gap_ev < 0.0:
        raise ValueError("link and neighbor-gap thresholds must be nonnegative")
    if not (0.0 <= min_pw_coverage <= 1.0):
        raise ValueError("min_pw_coverage must be in [0, 1]")
    if not (0.0 < max_abs_phi < np.pi):
        raise ValueError("max_abs_phi must be in (0, pi)")

    represented_energies = np.asarray(
        wavecar.energies[:, : max_band + 1], dtype=float
    )
    if not np.all(np.isfinite(represented_energies)):
        raise ValueError("full transport energy window contains nonfinite band energies")
    if np.any(np.diff(represented_energies, axis=1) < 0.0):
        raise ValueError("WAVECAR band energies are not nondecreasing")
    sentinel_min = float(np.min(represented_energies[:, max_band]))
    max_bundle_energy_max = float(
        np.max(represented_energies[:, max_band - 1])
    )
    indirect_gap_above_max_bundle = sentinel_min - max_bundle_energy_max
    if mu_max >= sentinel_min:
        raise ValueError(
            f"mu_max={mu_max:.9g} eV reaches sentinel band {max_band + 1} "
            f"(minimum {sentinel_min:.9g} eV); increase MAX_BAND before "
            "extending the transport window"
        )
    cumulative_maps = compute_cumulative_band_maps(wavecar, grid, max_band)
    phi_by_count = np.zeros((max_band + 1, grid.nx, grid.ny), dtype=float)
    link_by_count = np.full_like(phi_by_count, math.inf)
    gap_by_count = np.full_like(phi_by_count, math.inf)
    phase_abs_by_count = np.zeros_like(phi_by_count)
    finite_by_count = np.ones_like(phi_by_count, dtype=bool)
    map_summaries: list[dict[str, object]] = []
    cumulative_map_numeric_issues: list[dict[str, object]] = []
    minimum_coverage = math.inf
    nonfinite_coverage = False
    for band in range(1, max_band + 1):
        flux, xlinks, ylinks = cumulative_maps[band]
        phi_by_count[band] = flux
        link_by_count[band] = cell_min_link_singular_value(xlinks, ylinks, grid)
        phase_abs_by_count[band] = np.abs(flux)
        for ix in range(grid.nx):
            for iy in range(grid.ny):
                vertices = vertex_indices(grid, ix, iy)
                gap_by_count[band, ix, iy] = min(
                    float(
                        wavecar.energies[ik, band]
                        - wavecar.energies[ik, band - 1]
                    )
                    for ik in vertices
                )
        finite_by_count[band] = (
            np.isfinite(flux)
            & np.isfinite(link_by_count[band])
            & (link_by_count[band] > 0.0)
        )
        bad_cells = np.argwhere(~finite_by_count[band])
        if len(bad_cells):
            cumulative_map_numeric_issues.append(
                {
                    "occupied_band_count": band,
                    "invalid_cell_count": int(len(bad_cells)),
                    "invalid_cell_ids": [
                        int(iy * grid.nx + ix + 1)
                        for ix, iy in bad_cells[:20]
                    ],
                }
            )
        coverage_values = np.concatenate(
            (
                xlinks.coverage_left.ravel(), xlinks.coverage_right.ravel(),
                ylinks.coverage_left.ravel(), ylinks.coverage_right.ravel(),
            )
        )
        nonfinite_coverage = nonfinite_coverage or not bool(
            np.all(np.isfinite(coverage_values))
        )
        finite_coverage = coverage_values[np.isfinite(coverage_values)]
        if finite_coverage.size:
            minimum_coverage = min(
                minimum_coverage, float(np.min(finite_coverage))
            )
        with np.errstate(divide="ignore", invalid="ignore"):
            summary = map_summary(
                f"occ_{band}", list(range(1, band + 1)), flux,
                xlinks, ylinks, min_pw_coverage
            )
        map_summaries.append(summary)

    max_bundle_phase_finite = bool(
        np.all(np.isfinite(phi_by_count[max_band]))
    )
    max_bundle_link_finite_nonzero = bool(
        np.all(np.isfinite(link_by_count[max_band]))
        and np.all(link_by_count[max_band] > 0.0)
    )
    max_bundle_min_link = float(np.min(link_by_count[max_band]))
    max_bundle_min_gap = float(np.min(gap_by_count[max_band]))
    max_bundle_max_abs_phi = float(np.max(phase_abs_by_count[max_band]))
    max_bundle_quality = {
        "phase_finite": max_bundle_phase_finite,
        "links_finite_and_nonsingular": max_bundle_link_finite_nonzero,
        "min_link_singular_value": max_bundle_min_link,
        "min_link_pass": bool(
            max_bundle_link_finite_nonzero
            and max_bundle_min_link >= min_link_sv
        ),
        "min_gap_to_sentinel_ev": max_bundle_min_gap,
        "min_gap_to_sentinel_pass": bool(
            math.isfinite(max_bundle_min_gap)
            and max_bundle_min_gap >= min_neighbor_gap_ev
        ),
        "max_abs_phi_rad": max_bundle_max_abs_phi,
        "principal_branch_margin_pass": bool(
            math.isfinite(max_bundle_max_abs_phi)
            and max_bundle_max_abs_phi <= max_abs_phi
        ),
    }
    max_bundle_quality_pass = bool(
        max_bundle_quality["phase_finite"]
        and max_bundle_quality["min_link_pass"]
        and max_bundle_quality["min_gap_to_sentinel_pass"]
        and max_bundle_quality["principal_branch_margin_pass"]
    )
    max_bundle_quality["quality_pass"] = max_bundle_quality_pass

    diagnostics_path = output_dir / "transport_full_t0_diagnostics.json"
    if nonfinite_coverage:
        write_strict_json(
            diagnostics_path,
            {
                "method": "all-band Sawahata-style T=0 cumulative-subspace average",
                "max_cumulative_band": max_band,
                "validated": False,
                "failure": "nonfinite common plane-wave coverage",
                "cumulative_map_numeric_issues": cumulative_map_numeric_issues,
                "max_bundle_quality": max_bundle_quality,
                "cumulative_maps": map_summaries,
            },
        )
        raise RuntimeError(
            "full cumulative transport encountered nonfinite common plane-wave "
            f"coverage; {diagnostics_path.name} was written but output was refused"
        )
    if minimum_coverage < min_pw_coverage:
        write_strict_json(
            diagnostics_path,
            {
                "method": "all-band Sawahata-style T=0 cumulative-subspace average",
                "max_cumulative_band": max_band,
                "minimum_plane_wave_coverage": minimum_coverage,
                "required_plane_wave_coverage": min_pw_coverage,
                "validated": False,
                "failure": "plane-wave coverage below threshold",
                "cumulative_map_numeric_issues": cumulative_map_numeric_issues,
                "max_bundle_quality": max_bundle_quality,
                "cumulative_maps": map_summaries,
            },
        )
        raise RuntimeError(
            "full cumulative transport failed the common plane-wave coverage "
            f"gate ({minimum_coverage:.6g} < {min_pw_coverage:.6g}); "
            f"{diagnostics_path.name} was written but transport output was refused"
        )
    if not max_bundle_phase_finite:
        write_strict_json(
            diagnostics_path,
            {
                "method": (
                    "all-band Sawahata-style T=0 cumulative-subspace average"
                ),
                "max_cumulative_band": max_band,
                "validated": False,
                "failure": (
                    "full MAX_BAND reference bundle has undefined/nonfinite "
                    "plaquette phase"
                ),
                "cumulative_map_numeric_issues": cumulative_map_numeric_issues,
                "max_bundle_quality": max_bundle_quality,
                "cumulative_maps": map_summaries,
            },
        )
        raise RuntimeError(
            "full cumulative transport has an undefined/nonfinite MAX_BAND "
            "reference phase; the relative conductivity baseline cannot be "
            f"formed, so {diagnostics_path.name} was written and output was refused"
        )

    energy_vertices = np.empty((max_band, grid.nx, grid.ny, 4), dtype=float)
    for band in range(max_band):
        for ix in range(grid.nx):
            for iy in range(grid.ny):
                vertices = vertex_indices(grid, ix, iy)
                energy_vertices[band, ix, iy] = [
                    wavecar.energies[ik, band] for ik in vertices
                ]

    mask_k, mask_kp, mask_outside, distance_k, distance_kp = make_valley_partition(
        wavecar, grid, valley_k, valley_kp, valley_radius
    )
    if not np.all(
        mask_k.astype(np.int8)
        + mask_kp.astype(np.int8)
        + mask_outside.astype(np.int8)
        == 1
    ):
        raise RuntimeError("full-transport valley masks are not exclusive and exhaustive")

    # A uniform mu grid can step over a very narrow avoided crossing.  Audit
    # every occupation-event interval in the requested continuous range, in
    # addition to validating the explicitly sampled chemical potentials.
    continuous_invalid_events: list[dict[str, object]] = []
    invalid_interval_fragments: list[tuple[float, float]] = []
    for band in range(1, max_band + 1):
        finite_failure_for_band = ~finite_by_count[band]
        link_failure_for_band = link_by_count[band] < min_link_sv
        gap_failure_for_band = gap_by_count[band] < min_neighbor_gap_ev
        phase_failure_for_band = phase_abs_by_count[band] > max_abs_phi
        failed_cells = np.argwhere(
            finite_failure_for_band
            | link_failure_for_band
            | gap_failure_for_band
            | phase_failure_for_band
        )
        for ix, iy in failed_cells:
            lower = energy_vertices[band - 1, ix, iy]
            if band < max_band:
                upper = energy_vertices[band, ix, iy]
            else:
                vertices = vertex_indices(grid, int(ix), int(iy))
                upper = np.asarray(
                    [wavecar.energies[ik, max_band] for ik in vertices],
                    dtype=float,
                )
            for vertex in range(4):
                # Occupation count n is selected on [E_n, E_(n+1)).  The
                # requested scan interval is closed at both endpoints.
                if (
                    lower[vertex] < upper[vertex]
                    and mu_max >= lower[vertex]
                    and mu_min < upper[vertex]
                ):
                    lo = max(mu_min, float(lower[vertex]))
                    hi = min(mu_max, float(upper[vertex]))
                    if lo <= hi:
                        invalid_interval_fragments.append((lo, hi))
                        continuous_invalid_events.append(
                            {
                                "occupied_band_count": band,
                                "cell_id": int(iy * grid.nx + ix + 1),
                                "vertex": vertex,
                                "mu_interval_lower_ev": lo,
                                "mu_interval_upper_ev": hi,
                                "lower_bound_inclusive": True,
                                "upper_bound_inclusive": bool(
                                    hi == mu_max and mu_max < upper[vertex]
                                ),
                                "occupation_interval_convention": "[E_n,E_(n+1))",
                                "nonfinite_or_singular_failure": bool(
                                    finite_failure_for_band[ix, iy]
                                ),
                                "link_failure": bool(link_failure_for_band[ix, iy]),
                                "gap_failure": bool(gap_failure_for_band[ix, iy]),
                                "phase_failure": bool(phase_failure_for_band[ix, iy]),
                                "min_link_sv": float(link_by_count[band, ix, iy]),
                                "min_neighbor_gap_ev": float(gap_by_count[band, ix, iy]),
                                "abs_phi_rad": float(phase_abs_by_count[band, ix, iy]),
                            }
                        )

    merged_invalid_intervals: list[list[float]] = []
    for lo, hi in sorted(invalid_interval_fragments):
        if not merged_invalid_intervals or lo > merged_invalid_intervals[-1][1]:
            merged_invalid_intervals.append([lo, hi])
        else:
            merged_invalid_intervals[-1][1] = max(
                merged_invalid_intervals[-1][1], hi
            )
    continuous_range_validated = not continuous_invalid_events

    rows: list[dict[str, float | int | str]] = []
    invalid_rows: list[dict[str, object]] = []
    two_pi = 2.0 * np.pi
    max_bundle_chern = float(np.sum(phi_by_count[max_band]) / two_pi)
    max_bundle_chern_k = float(np.sum(phi_by_count[max_band][mask_k]) / two_pi)
    max_bundle_chern_kp = float(np.sum(phi_by_count[max_band][mask_kp]) / two_pi)
    max_bundle_chern_outside = float(
        np.sum(phi_by_count[max_band][mask_outside]) / two_pi
    )
    max_bundle_chern_nearest_integer = int(round(max_bundle_chern))
    max_bundle_chern_nearest_integer_residual = abs(
        max_bundle_chern - max_bundle_chern_nearest_integer
    )
    max_bundle_sigma = -max_bundle_chern
    max_bundle_sigma_k = -max_bundle_chern_k
    max_bundle_sigma_kp = -max_bundle_chern_kp
    max_bundle_sigma_outside = -max_bundle_chern_outside
    max_bundle_sigma_valley_contrast = max_bundle_sigma_k - max_bundle_sigma_kp
    for mu in np.linspace(mu_min, mu_max, mu_num):
        effective_phi, occupied_count = full_cumulative_t0_effective_phi(
            phi_by_count, energy_vertices, float(mu)
        )
        used = np.zeros_like(phi_by_count, dtype=bool)
        for band in range(max_band + 1):
            used[band] = np.any(occupied_count == band, axis=2)
        used_nonempty = used[1:]
        finite_failure = np.any(
            used_nonempty & ~finite_by_count[1:], axis=0
        )
        link_failure = np.any(
            used_nonempty & (link_by_count[1:] < min_link_sv), axis=0
        )
        gap_failure = np.any(
            used_nonempty & (gap_by_count[1:] < min_neighbor_gap_ev), axis=0
        )
        phase_failure = np.any(
            used_nonempty & (phase_abs_by_count[1:] > max_abs_phi), axis=0
        )
        invalid_cells = finite_failure | link_failure | gap_failure | phase_failure
        quality_pass = not bool(np.any(invalid_cells))

        selected_links = np.where(used_nonempty, link_by_count[1:], math.inf)
        selected_gaps = np.where(used_nonempty, gap_by_count[1:], math.inf)
        selected_phases = np.where(used_nonempty, phase_abs_by_count[1:], 0.0)
        worst_link = float(np.min(selected_links))
        worst_gap = float(np.min(selected_gaps))
        worst_phi = float(np.max(selected_phases))
        if not quality_pass:
            ids = np.argwhere(invalid_cells)
            invalid_rows.append(
                {
                    "mu_ev": float(mu),
                    "worst_link_sv": worst_link,
                    "worst_neighbor_gap_ev": worst_gap,
                    "max_abs_phi_rad": worst_phi,
                    "link_failure_cell_count": int(np.count_nonzero(link_failure)),
                    "gap_failure_cell_count": int(np.count_nonzero(gap_failure)),
                    "phase_failure_cell_count": int(np.count_nonzero(phase_failure)),
                    "nonfinite_or_singular_failure_cell_count": int(
                        np.count_nonzero(finite_failure)
                    ),
                    "invalid_cell_count": int(len(ids)),
                    "invalid_cell_ids": [
                        int(iy * grid.nx + ix + 1) for ix, iy in ids[:20]
                    ],
                }
            )

        transport_finite = bool(np.all(np.isfinite(effective_phi)))
        if transport_finite:
            chern_total = float(np.sum(effective_phi) / two_pi)
            chern_k = float(np.sum(effective_phi[mask_k]) / two_pi)
            chern_kp = float(np.sum(effective_phi[mask_kp]) / two_pi)
            chern_outside = float(np.sum(effective_phi[mask_outside]) / two_pi)
        else:
            chern_total = chern_k = chern_kp = chern_outside = None

        def finite_difference(value: float | None, baseline: float) -> float | str:
            if value is None or not math.isfinite(baseline):
                return ""
            return value - baseline

        rows.append(
            {
                "mu_eV": float(mu),
                "quality": "PASS" if quality_pass else "INVALID",
                "continuous_scan_quality": (
                    "PASS" if continuous_range_validated
                    else "CONTAINS_INVALID_SUBINTERVALS"
                ),
                "occupied_band_count_min": int(np.min(occupied_count)),
                "occupied_band_count_max": int(np.max(occupied_count)),
                "active_cell_count": int(np.count_nonzero(np.any(used_nonempty, axis=0))),
                "nonfinite_or_singular_active_cell_count": int(
                    np.count_nonzero(finite_failure)
                ),
                "link_failure_active_cell_count": int(
                    np.count_nonzero(link_failure)
                ),
                "gap_failure_active_cell_count": int(
                    np.count_nonzero(gap_failure)
                ),
                "phase_failure_active_cell_count": int(
                    np.count_nonzero(phase_failure)
                ),
                "worst_active_link_sv": (
                    worst_link if math.isfinite(worst_link) else ""
                ),
                "worst_active_neighbor_gap_eV": (
                    worst_gap if math.isfinite(worst_gap) else ""
                ),
                "max_active_abs_phi_rad": (
                    worst_phi if math.isfinite(worst_phi) else ""
                ),
                "chern_total": chern_total if chern_total is not None else "",
                "sigma_total_e2_over_h": (
                    -chern_total if chern_total is not None else ""
                ),
                "chern_K": chern_k if chern_k is not None else "",
                "chern_Kp": chern_kp if chern_kp is not None else "",
                "chern_outside": chern_outside if chern_outside is not None else "",
                "sigma_K_e2_over_h": -chern_k if chern_k is not None else "",
                "sigma_Kp_e2_over_h": -chern_kp if chern_kp is not None else "",
                "sigma_outside_e2_over_h": (
                    -chern_outside if chern_outside is not None else ""
                ),
                "sigma_valley_contrast_e2_over_h": (
                    -chern_k + chern_kp if chern_k is not None else ""
                ),
                "sigma_relative_to_full_max_bundle_e2_over_h": (
                    finite_difference(
                        -chern_total if chern_total is not None else None,
                        max_bundle_sigma,
                    )
                ),
                "sigma_K_relative_to_full_max_bundle_e2_over_h": (
                    finite_difference(
                        -chern_k if chern_k is not None else None,
                        max_bundle_sigma_k,
                    )
                ),
                "sigma_Kp_relative_to_full_max_bundle_e2_over_h": (
                    finite_difference(
                        -chern_kp if chern_kp is not None else None,
                        max_bundle_sigma_kp,
                    )
                ),
                "sigma_outside_relative_to_full_max_bundle_e2_over_h": (
                    finite_difference(
                        -chern_outside if chern_outside is not None else None,
                        max_bundle_sigma_outside,
                    )
                ),
                "sigma_valley_contrast_relative_to_full_max_bundle_e2_over_h": (
                    finite_difference(
                        (-chern_k + chern_kp) if chern_k is not None else None,
                        max_bundle_sigma_valley_contrast,
                    )
                ),
            }
        )

    plateau_rows = [
        row for row in rows
        if row["occupied_band_count_min"] == max_band
        and row["occupied_band_count_max"] == max_band
    ]
    plateau_validated_rows = [
        row for row in plateau_rows
        if row["quality"] == "PASS" and row["chern_total"] != ""
    ]
    plateau_max_chern_error = (
        max(
            abs(float(row["chern_total"]) - max_bundle_chern)
            for row in plateau_validated_rows
        )
        if plateau_validated_rows else None
    )
    plateau_quality_pass = (
        bool(plateau_rows)
        and len(plateau_validated_rows) == len(plateau_rows)
        and max_bundle_quality_pass
    )
    validation_failure_reasons: list[str] = []
    if not max_bundle_quality_pass:
        validation_failure_reasons.append("MAX_BAND reference bundle quality gate")
    if invalid_rows:
        validation_failure_reasons.append("sampled active-cell quality gate")
    if continuous_invalid_events:
        validation_failure_reasons.append(
            "continuous occupation-event interval quality gate"
        )
    diagnostics = {
        "method": "all-band Sawahata-style T=0 cumulative-subspace plaquette average",
        "occupation_convention": "E <= mu is occupied",
        "gauge_treatment": (
            "determinant links of cumulative subspaces 1:n; no isolated-band "
            "Berry-curvature weighting inside degeneracies"
        ),
        "max_cumulative_band": max_band,
        "sentinel_band": max_band + 1,
        "sentinel_min_ev": sentinel_min,
        "max_bundle_energy_max_ev": max_bundle_energy_max,
        "indirect_gap_above_max_bundle_ev": indirect_gap_above_max_bundle,
        "positive_indirect_gap_above_max_bundle": (
            indirect_gap_above_max_bundle > 0.0
        ),
        "max_bundle_plateau_in_requested_range": (
            indirect_gap_above_max_bundle > 0.0
            and mu_max >= max_bundle_energy_max
        ),
        "max_bundle_plateau_sample_count": len(plateau_rows),
        "max_bundle_plateau_validated_sample_count": len(
            plateau_validated_rows
        ),
        "max_bundle_plateau_quality_pass": plateau_quality_pass,
        "max_bundle_chern": max_bundle_chern,
        "max_bundle_chern_nearest_integer": max_bundle_chern_nearest_integer,
        "max_bundle_chern_nearest_integer_residual": (
            max_bundle_chern_nearest_integer_residual
        ),
        "max_bundle_sigma_e2_over_h": max_bundle_sigma,
        "max_bundle_sigma_valley_contrast_e2_over_h": (
            max_bundle_sigma_valley_contrast
        ),
        "max_bundle_plateau_max_chern_error": plateau_max_chern_error,
        "mu_range_ev": [mu_min, mu_max, mu_num],
        "mu_max_below_sentinel_min": mu_max < sentinel_min,
        "quality_thresholds": {
            "min_link_singular_value": min_link_sv,
            "min_neighbor_gap_ev": min_neighbor_gap_ev,
            "min_plane_wave_coverage": min_pw_coverage,
            "max_abs_phi_rad": max_abs_phi,
        },
        "minimum_plane_wave_coverage": minimum_coverage,
        "max_bundle_quality": max_bundle_quality,
        "cumulative_map_numeric_issues": cumulative_map_numeric_issues,
        "cumulative_maps": map_summaries,
        "valley_partition": {
            "K_fractional": valley_k.tolist(),
            "Kp_fractional": valley_kp.tolist(),
            "radius_inv_angstrom": valley_radius,
            "mode": "radius" if valley_radius is not None else "nearest-center Voronoi",
            "K_cell_count": int(np.count_nonzero(mask_k)),
            "Kp_cell_count": int(np.count_nonzero(mask_kp)),
            "outside_cell_count": int(np.count_nonzero(mask_outside)),
            "min_center_distance_K_inv_angstrom": float(np.min(distance_k)),
            "min_center_distance_Kp_inv_angstrom": float(np.min(distance_kp)),
        },
        "validation_scope": (
            "all requested discrete mu points plus every occupation-event "
            "subinterval in the closed requested mu range"
        ),
        "validated": (
            max_bundle_quality_pass
            and continuous_range_validated
            and not invalid_rows
        ),
        "validation_failure_reasons": validation_failure_reasons,
        "max_bundle_reference_validated": max_bundle_quality_pass,
        "requested_mu_points_validated": not invalid_rows,
        "continuous_mu_range_validated": continuous_range_validated,
        "continuous_invalid_event_count": len(continuous_invalid_events),
        "continuous_invalid_events": continuous_invalid_events,
        "merged_invalid_mu_intervals_ev": merged_invalid_intervals,
        "invalid_mu_count": len(invalid_rows),
        "invalid_mu": invalid_rows,
    }
    write_strict_json(diagnostics_path, diagnostics)
    if (
        not max_bundle_quality_pass
        or invalid_rows
        or continuous_invalid_events
    ) and not allow_invalid:
        raise RuntimeError(
            "full cumulative transport validation failed for the MAX_BAND "
            "reference, at sampled mu points, or within an occupation-event "
            "subinterval; "
            f"{diagnostics_path.name} was written but transport_full_t0.csv was "
            "refused. Use --allow-invalid-transport only for explicitly "
            "unvalidated diagnostics."
        )
    with (output_dir / "transport_full_t0.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def periodic_bilinear_sample(
    values: np.ndarray,
    qpoints: np.ndarray,
    grid: Grid,
    centered: bool,
) -> np.ndarray:
    """Periodic bilinear interpolation on a vertex or plaquette-center grid."""

    result = np.empty(len(qpoints), dtype=float)
    half = 0.5 if centered else 0.0
    for iq, qpoint in enumerate(qpoints):
        x = ((qpoint[0] - grid.offset[0]) % 1.0) * grid.nx - half
        y = ((qpoint[1] - grid.offset[1]) % 1.0) * grid.ny - half
        ix0, iy0 = math.floor(x), math.floor(y)
        tx, ty = x - ix0, y - iy0
        ix1, iy1 = (ix0 + 1) % grid.nx, (iy0 + 1) % grid.ny
        ix0, iy0 = ix0 % grid.nx, iy0 % grid.ny
        result[iq] = (
            (1.0 - tx) * (1.0 - ty) * values[ix0, iy0]
            + tx * (1.0 - ty) * values[ix1, iy0]
            + (1.0 - tx) * ty * values[ix0, iy1]
            + tx * ty * values[ix1, iy1]
        )
    return result


def closest_periodic_image(
    start: np.ndarray, endpoint: np.ndarray, reciprocal: np.ndarray
) -> np.ndarray:
    return np.asarray(start, dtype=float) + shortest_periodic_fractional_displacement(
        endpoint, start, reciprocal
    )


def _matplotlib_pyplot():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def periodic_image_in_plot_domain(
    fractional_point: np.ndarray, grid_offset: np.ndarray
) -> np.ndarray:
    """Map a periodic q point into [offset, offset+1) for shifted-mesh plots."""

    result = np.asarray(fractional_point, dtype=float).copy()
    result[:2] = np.mod(result[:2] - grid_offset[:2], 1.0) + grid_offset[:2]
    return result


def fractional_plot_extent(grid_offset: np.ndarray) -> tuple[float, float, float, float]:
    """Return the one-period plotting domain anchored at the mesh offset."""

    return (
        float(grid_offset[0]),
        float(grid_offset[0] + 1.0),
        float(grid_offset[1]),
        float(grid_offset[1] + 1.0),
    )


def reciprocal_plane_frame(reciprocal: np.ndarray) -> np.ndarray:
    """Return an oriented orthonormal frame spanning reciprocal b1 and b2.

    Rows of the returned 2x3 array are the Cartesian axes used for physical
    first-Brillouin-zone plots.  This keeps the plotting geometry valid even
    when the two-dimensional reciprocal plane is not the Cartesian xy plane.
    """

    basis = np.asarray(reciprocal, dtype=float)
    if basis.shape != (3, 3) or not np.all(np.isfinite(basis)):
        raise ValueError("reciprocal basis must be a finite 3x3 array")
    vector_x = basis[0]
    norm_x = float(np.linalg.norm(vector_x))
    normal = np.cross(basis[0], basis[1])
    norm_normal = float(np.linalg.norm(normal))
    if norm_x <= 0.0 or norm_normal <= np.finfo(float).eps * norm_x * max(
        float(np.linalg.norm(basis[1])), 1.0
    ):
        raise ValueError("in-plane reciprocal basis must be rank two")
    vector_x = vector_x / norm_x
    normal = normal / norm_normal
    vector_y = np.cross(normal, vector_x)
    # The cross-product construction should already orient b2 along +y.  Keep
    # that convention explicit so K/K' plots are deterministic under roundoff.
    if float(np.dot(basis[1], vector_y)) < 0.0:
        vector_y = -vector_y
    return np.stack((vector_x, vector_y))


def polygon_signed_area(vertices: np.ndarray) -> float:
    """Signed area of a two-dimensional polygon."""

    polygon = np.asarray(vertices, dtype=float)
    if polygon.ndim != 2 or polygon.shape[1] != 2 or polygon.shape[0] < 3:
        raise ValueError("polygon must contain at least three 2D vertices")
    return 0.5 * float(
        np.sum(
            polygon[:, 0] * np.roll(polygon[:, 1], -1)
            - polygon[:, 1] * np.roll(polygon[:, 0], -1)
        )
    )


def clip_convex_polygon_halfplane(
    vertices: np.ndarray,
    normal: np.ndarray,
    limit: float,
    tolerance: float,
) -> np.ndarray:
    """Clip a convex polygon to points satisfying dot(point, normal)<=limit."""

    polygon = np.asarray(vertices, dtype=float)
    clipped: list[np.ndarray] = []
    for start, stop in zip(polygon, np.roll(polygon, -1, axis=0)):
        start_value = float(np.dot(start, normal) - limit)
        stop_value = float(np.dot(stop, normal) - limit)
        start_inside = start_value <= tolerance
        stop_inside = stop_value <= tolerance
        if start_inside:
            clipped.append(start)
        if start_inside != stop_inside:
            denominator = start_value - stop_value
            if abs(denominator) <= np.finfo(float).tiny:
                continue
            fraction = start_value / denominator
            clipped.append(start + fraction * (stop - start))
    if not clipped:
        return np.empty((0, 2), dtype=float)
    result = np.asarray(clipped, dtype=float)
    scale = float(np.max(np.linalg.norm(result, axis=1)))
    if not math.isfinite(scale) or scale <= 0.0:
        raise RuntimeError("half-plane clipping produced a degenerate polygon")
    keep = np.ones(result.shape[0], dtype=bool)
    for index in range(result.shape[0]):
        duplicate_tolerance = 64.0 * np.finfo(float).eps * scale
        if np.linalg.norm(result[index] - result[index - 1]) <= duplicate_tolerance:
            keep[index] = False
    return result[keep]


def wigner_seitz_polygon_2d(reciprocal: np.ndarray) -> np.ndarray:
    """Construct the first 2D Brillouin zone as a Wigner--Seitz polygon.

    The reciprocal basis is first Gauss-reduced.  In two dimensions the
    Voronoi-relevant vectors are contained among short combinations of a
    reduced basis; the [-2,2]^2 enumeration also supplies redundant guards and
    is independent of how skew the input b1,b2 happen to be.
    """

    reciprocal = np.asarray(reciprocal, dtype=float)
    frame = reciprocal_plane_frame(reciprocal)
    reduced_cartesian = gauss_reduce_2d_lattice_basis(reciprocal[:2])
    reduced = reduced_cartesian @ frame.T
    scale = float(np.max(np.linalg.norm(reduced, axis=1)))
    if not math.isfinite(scale) or scale <= 0.0:
        raise ValueError("in-plane reciprocal basis must be finite and rank two")
    radius = 4.0 * scale
    polygon = np.asarray(
        ((-radius, -radius), (radius, -radius), (radius, radius), (-radius, radius)),
        dtype=float,
    )
    candidates: list[np.ndarray] = []
    for first in range(-2, 3):
        for second in range(-2, 3):
            if first == 0 and second == 0:
                continue
            vector = first * reduced[0] + second * reduced[1]
            candidates.append(vector)
    candidates.sort(key=lambda vector: float(np.dot(vector, vector)))
    tolerance = 256.0 * np.finfo(float).eps * scale * scale
    for vector in candidates:
        limit = 0.5 * float(np.dot(vector, vector))
        polygon = clip_convex_polygon_halfplane(
            polygon, vector, limit, tolerance
        )
        if polygon.shape[0] < 3:
            raise RuntimeError("Wigner--Seitz half-plane intersection became empty")
    area = polygon_signed_area(polygon)
    if area < 0.0:
        polygon = polygon[::-1].copy()
        area = -area
    expected_area = float(np.linalg.norm(np.cross(reciprocal[0], reciprocal[1])))
    area_tolerance = 256.0 * np.finfo(float).eps * expected_area
    if not math.isclose(
        area, expected_area, rel_tol=2.0e-10, abs_tol=area_tolerance
    ):
        raise RuntimeError(
            "Wigner--Seitz polygon area does not match the reciprocal primitive "
            f"cell ({area:.12g} versus {expected_area:.12g} A^-2)"
        )
    return polygon


def fold_fractional_to_first_bz(
    fractional_point: np.ndarray, reciprocal: np.ndarray
) -> np.ndarray:
    """Fold a fractional reciprocal point into the first 2D Brillouin zone."""

    point = np.asarray(fractional_point, dtype=float)
    if point.shape != (3,) or not np.all(np.isfinite(point)):
        raise ValueError(
            "fractional reciprocal point must have three finite components"
        )
    reciprocal = np.asarray(reciprocal, dtype=float)
    frame = reciprocal_plane_frame(reciprocal)
    reduced = gauss_reduce_2d_lattice_basis(reciprocal[:2])
    in_plane_cartesian = point[:2] @ reciprocal[:2]
    folded_cartesian = closest_2d_lattice_residual(
        -in_plane_cartesian, reduced
    )
    return folded_cartesian @ frame.T


def periodic_images_in_first_bz(
    fractional_point: np.ndarray,
    reciprocal: np.ndarray,
    polygon: np.ndarray | None = None,
) -> np.ndarray:
    """Return every periodically equivalent image lying in the closed first BZ.

    An interior point has one image.  Points on Wigner--Seitz edges or vertices
    have multiple closed-cell representatives; for a hexagonal reciprocal
    lattice a K class therefore appears at three alternating BZ corners.
    """

    reciprocal = np.asarray(reciprocal, dtype=float)
    if polygon is None:
        polygon = wigner_seitz_polygon_2d(reciprocal)
    else:
        polygon = np.asarray(polygon, dtype=float)
    frame = reciprocal_plane_frame(reciprocal)
    reduced = gauss_reduce_2d_lattice_basis(reciprocal[:2]) @ frame.T
    representative = fold_fractional_to_first_bz(
        fractional_point, reciprocal
    )
    lattice_scale = float(np.max(np.linalg.norm(reduced, axis=1)))
    if not math.isfinite(lattice_scale) or lattice_scale <= 0.0:
        raise ValueError("in-plane reciprocal basis must be finite and rank two")
    length_tolerance = 2.0e-10 * lattice_scale
    # point_in_convex_polygon_2d compares a 2D cross product, so its tolerance
    # has reciprocal-area dimensions rather than reciprocal-length dimensions.
    area_tolerance = 2.0e-10 * lattice_scale * lattice_scale
    images: list[np.ndarray] = []
    for first in range(-2, 3):
        for second in range(-2, 3):
            candidate = representative + first * reduced[0] + second * reduced[1]
            if not bool(
                point_in_convex_polygon_2d(
                    candidate, polygon, tolerance=area_tolerance
                )
            ):
                continue
            if any(
                np.linalg.norm(candidate - existing) <= length_tolerance
                for existing in images
            ):
                continue
            images.append(candidate)
    if not images:
        raise RuntimeError("no periodic representative was found in the first BZ")
    result = np.asarray(images)
    angles = np.arctan2(result[:, 1], result[:, 0])
    return result[np.argsort(angles)]


def point_in_convex_polygon_2d(
    points: np.ndarray,
    polygon: np.ndarray,
    tolerance: float | None = None,
) -> np.ndarray:
    """Return a mask for points inside or on a counter-clockwise polygon."""

    values = np.asarray(points, dtype=float)
    vertices = np.asarray(polygon, dtype=float)
    edges = np.roll(vertices, -1, axis=0) - vertices
    if tolerance is None:
        edge_scale = float(np.max(np.linalg.norm(edges, axis=1)))
        tolerance = 64.0 * np.finfo(float).eps * edge_scale * edge_scale
    if not math.isfinite(tolerance) or tolerance < 0.0:
        raise ValueError("polygon cross-product tolerance must be finite and nonnegative")
    relative = values[..., None, :] - vertices
    cross = edges[:, 0] * relative[..., :, 1] - edges[:, 1] * relative[..., :, 0]
    return np.all(cross >= -tolerance, axis=-1)


def first_bz_display_grid(
    data: np.ndarray,
    grid: Grid,
    reciprocal: np.ndarray,
    polygon: np.ndarray,
    maximum_pixels: int = 420,
) -> tuple[np.ndarray, tuple[float, float, float, float]]:
    """Rasterize plaquette data periodically in Cartesian first-BZ coordinates.

    Each display pixel receives the value of the plaquette containing its
    inverse-mapped fractional coordinate.  No interpolation is used, so this
    changes only the display domain and not the underlying nx-by-ny data.
    """

    values = np.asarray(data, dtype=float)
    if values.shape != (grid.nx, grid.ny):
        raise ValueError("first-BZ display data must match the uniform grid")
    frame = reciprocal_plane_frame(reciprocal)
    basis_2d = np.asarray(reciprocal[:2], dtype=float) @ frame.T
    determinant = float(np.linalg.det(basis_2d))
    if abs(determinant) <= np.finfo(float).eps:
        raise ValueError("in-plane reciprocal basis must be rank two")
    x_min, y_min = np.min(polygon, axis=0)
    x_max, y_max = np.max(polygon, axis=0)
    width, height = float(x_max - x_min), float(y_max - y_min)
    longest = max(width, height)
    pixels_x = max(2, int(round(maximum_pixels * width / longest)))
    pixels_y = max(2, int(round(maximum_pixels * height / longest)))
    x = np.linspace(x_min, x_max, pixels_x, endpoint=False) + 0.5 * width / pixels_x
    y = np.linspace(y_min, y_max, pixels_y, endpoint=False) + 0.5 * height / pixels_y
    xx, yy = np.meshgrid(x, y, indexing="xy")
    plane_points = np.stack((xx, yy), axis=-1)
    fractional = plane_points @ np.linalg.inv(basis_2d)
    ix = np.floor(
        np.mod(fractional[..., 0] - grid.offset[0], 1.0) * grid.nx
    ).astype(int)
    iy = np.floor(
        np.mod(fractional[..., 1] - grid.offset[1], 1.0) * grid.ny
    ).astype(int)
    raster = values[ix, iy]
    raster[~point_in_convex_polygon_2d(plane_points, polygon)] = np.nan
    return raster, (float(x_min), float(x_max), float(y_min), float(y_max))


def padded_cartesian_plot_limits(
    extent: tuple[float, float, float, float], padding_fraction: float = 0.05
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Add display-only padding around a Cartesian reciprocal-space extent."""

    x_min, x_max, y_min, y_max = (float(value) for value in extent)
    if not all(math.isfinite(value) for value in (x_min, x_max, y_min, y_max)):
        raise ValueError("Cartesian plot extent must be finite")
    if x_max <= x_min or y_max <= y_min:
        raise ValueError("Cartesian plot extent must have positive width and height")
    if not math.isfinite(padding_fraction) or padding_fraction < 0.0:
        raise ValueError("Cartesian plot padding must be finite and nonnegative")
    pad_x = padding_fraction * (x_max - x_min)
    pad_y = padding_fraction * (y_max - y_min)
    return ((x_min - pad_x, x_max + pad_x), (y_min - pad_y, y_max + pad_y))


def symmetric_finite_color_limits(values: np.ndarray) -> tuple[float, float]:
    """Return zero-centered limits, with a safe fallback for all-zero data."""

    finite = np.abs(np.asarray(values, dtype=float)[np.isfinite(values)])
    maximum = float(np.max(finite)) if finite.size else 0.0
    if maximum == 0.0:
        maximum = 1.0
    return -maximum, maximum


def plot_k_resolved_maps(
    output_dir: Path,
    wavecar: Wavecar,
    grid: Grid,
    map_results: dict[str, tuple[list[int], np.ndarray, LinkSet, LinkSet]],
    energy_band: int,
    valley_k: np.ndarray,
    valley_kp: np.ndarray,
    plot_domain: str = "fractional",
) -> Path:
    plt = _matplotlib_pyplot()
    by_range = {
        (bands[0], bands[-1]): (label, flux, xlinks, ylinks)
        for label, (bands, flux, xlinks, ylinks) in map_results.items()
    }
    required = ((1, energy_band - 1), (1, energy_band), (1, energy_band + 1))
    if any(item not in by_range for item in required):
        raise ValueError("k-resolved plot requires V, V+1, and V+2 cumulative maps")
    _, phi_v, _, _ = by_range[required[0]]
    _, phi_v1, x_v1, y_v1 = by_range[required[1]]
    _, phi_v2, _, _ = by_range[required[2]]
    area = float(np.linalg.norm(np.cross(*wavecar.header.reciprocal[:2]))) / (
        grid.nx * grid.ny
    )
    mean_energy = np.empty((grid.nx, grid.ny))
    min_gap = np.empty_like(mean_energy)
    for ix in range(grid.nx):
        for iy in range(grid.ny):
            vertices = vertex_indices(grid, ix, iy)
            e1 = np.asarray(
                [wavecar.energies[ik, energy_band - 1] for ik in vertices]
            )
            e2 = np.asarray([wavecar.energies[ik, energy_band] for ik in vertices])
            mean_energy[ix, iy] = float(np.mean(e1))
            min_gap[ix, iy] = float(np.min(e2 - e1))
    link_quality = cell_min_link_singular_value(x_v1, y_v1, grid)
    panels = (
        (mean_energy, f"Band {energy_band} mean vertex energy (eV)", "viridis"),
        (phi_v / area, f"Cumulative 1:{energy_band - 1} $\\Omega$ ($\\AA^2$)", "RdBu_r"),
        (phi_v1 / area, f"Cumulative 1:{energy_band} raw $\\Omega$ ($\\AA^2$)", "RdBu_r"),
        (phi_v2 / area, f"Cumulative 1:{energy_band + 1} $\\Omega$ ($\\AA^2$)", "RdBu_r"),
        (np.log10(np.maximum(link_quality, 1e-300)), "V+1 log10(min link singular value)", "magma"),
        (np.log10(np.maximum(min_gap, 1e-300)), f"log10 min E{energy_band + 1}-E{energy_band} (eV)", "magma"),
    )
    if plot_domain not in ("fractional", "first-bz"):
        raise ValueError("plot_domain must be 'fractional' or 'first-bz'")
    if plot_domain == "fractional":
        plot_extent = fractional_plot_extent(grid.offset)
        plot_k_images = periodic_image_in_plot_domain(
            valley_k, grid.offset
        )[None, :2]
        plot_kp_images = periodic_image_in_plot_domain(
            valley_kp, grid.offset
        )[None, :2]
        bz_polygon = None
    else:
        bz_polygon = wigner_seitz_polygon_2d(wavecar.header.reciprocal)
        plot_extent = None
        plot_k_images = periodic_images_in_first_bz(
            valley_k, wavecar.header.reciprocal, bz_polygon
        )
        plot_kp_images = periodic_images_in_first_bz(
            valley_kp, wavecar.header.reciprocal, bz_polygon
        )
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
    for axis, (data, title, cmap) in zip(axes.flat, panels):
        color_limits = (
            symmetric_finite_color_limits(data) if cmap == "RdBu_r" else (None, None)
        )
        if plot_domain == "fractional":
            display_data = data.T
            display_extent = plot_extent
        else:
            assert bz_polygon is not None
            display_data, display_extent = first_bz_display_grid(
                data, grid, wavecar.header.reciprocal, bz_polygon
            )
        image = axis.imshow(
            display_data,
            origin="lower",
            interpolation="nearest",
            aspect="equal",
            cmap=cmap,
            extent=display_extent,
            vmin=color_limits[0],
            vmax=color_limits[1],
        )
        axis.set_title(title)
        if plot_domain == "fractional":
            axis.set_xlabel("fractional $q_1$")
            axis.set_ylabel("fractional $q_2$")
        else:
            assert bz_polygon is not None
            boundary = np.vstack((bz_polygon, bz_polygon[0]))
            axis.plot(boundary[:, 0], boundary[:, 1], color="black", linewidth=1.0)
            axis.set_xlabel(r"$k_x$ ($\AA^{-1}$)")
            axis.set_ylabel(r"$k_y$ ($\AA^{-1}$)")
            assert display_extent is not None
            x_limits, y_limits = padded_cartesian_plot_limits(display_extent)
            axis.set_xlim(*x_limits)
            axis.set_ylim(*y_limits)
            axis.set_aspect("equal", adjustable="box")
        axis.scatter(
            plot_k_images[:, 0],
            plot_k_images[:, 1],
            marker="*",
            s=70,
            color="gold",
            edgecolor="black",
            linewidth=0.5,
        )
        axis.scatter(
            plot_kp_images[:, 0],
            plot_kp_images[:, 1],
            marker="x",
            s=50,
            color="black",
            linewidth=1.2,
        )
        for label, markers in (("K", plot_k_images), ("K'", plot_kp_images)):
            for marker in markers:
                if plot_domain == "fractional":
                    annotation_offset = (4.0, 4.0)
                else:
                    # Boundary labels point inward toward Gamma, avoiding axes
                    # clipping and making the alternating K/K' corners clear.
                    marker_norm = float(np.linalg.norm(marker))
                    annotation_offset = (
                        tuple(-10.0 * marker / marker_norm)
                        if marker_norm > 0.0
                        else (4.0, 4.0)
                    )
                axis.annotate(
                    label,
                    (marker[0], marker[1]),
                    xytext=annotation_offset,
                    textcoords="offset points",
                    ha="center",
                    va="center",
                    zorder=6,
                )
        fig.colorbar(image, ax=axis, shrink=0.82)
    domain_title = (
        "fractional reciprocal primitive cell"
        if plot_domain == "fractional"
        else "Cartesian first Brillouin zone"
    )
    fig.suptitle(
        "WAVECAR-direct Fukui maps in " + domain_title
        + " (raw maps include quality diagnostics)"
    )
    output = output_dir / "wavecar_fukui_kresolved.png"
    fig.savefig(output, dpi=180)
    plt.close(fig)
    return output


def plot_k_to_kp_line(
    output_dir: Path,
    wavecar: Wavecar,
    grid: Grid,
    map_results: dict[str, tuple[list[int], np.ndarray, LinkSet, LinkSet]],
    energy_band: int,
    valley_k: np.ndarray,
    valley_kp: np.ndarray,
) -> Path:
    plt = _matplotlib_pyplot()
    endpoint = closest_periodic_image(
        valley_k, valley_kp, wavecar.header.reciprocal
    )
    parameter = np.linspace(0.0, 1.0, 401)
    qpoints = valley_k[None, :] + parameter[:, None] * (endpoint - valley_k)[None, :]
    segment_cart = (endpoint - valley_k) @ wavecar.header.reciprocal
    distance = parameter * float(np.linalg.norm(segment_cart))
    energy1_grid = np.empty((grid.nx, grid.ny))
    energy2_grid = np.empty_like(energy1_grid)
    for ix in range(grid.nx):
        for iy in range(grid.ny):
            ik = int(grid.index[ix, iy])
            energy1_grid[ix, iy] = wavecar.energies[ik, energy_band - 1]
            energy2_grid[ix, iy] = wavecar.energies[ik, energy_band]
    line_e1 = periodic_bilinear_sample(energy1_grid, qpoints, grid, centered=False)
    line_e2 = periodic_bilinear_sample(energy2_grid, qpoints, grid, centered=False)

    by_range = {
        (bands[0], bands[-1]): (label, flux)
        for label, (bands, flux, _xlinks, _ylinks) in map_results.items()
    }
    area = float(np.linalg.norm(np.cross(*wavecar.header.reciprocal[:2]))) / (
        grid.nx * grid.ny
    )
    line_fields: list[tuple[str, np.ndarray, str]] = []
    for band_range, color in (
        ((1, energy_band - 1), "black"),
        ((1, energy_band), "tab:red"),
        ((1, energy_band + 1), "tab:blue"),
    ):
        if band_range in by_range:
            label, flux = by_range[band_range]
            line_fields.append(
                (
                    label,
                    periodic_bilinear_sample(flux / area, qpoints, grid, centered=True),
                    color,
                )
            )

    fig, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=True, constrained_layout=True)
    axes[0].plot(distance, line_e1, label=f"band {energy_band}")
    axes[0].plot(distance, line_e2, label=f"band {energy_band + 1}")
    axes[0].set_ylabel("energy (eV)")
    axes[0].legend()
    axes[0].grid(alpha=0.25)
    for label, values, color in line_fields:
        axes[1].plot(distance, values, label=label, color=color)
    axes[1].axhline(0.0, color="0.5", linewidth=0.8)
    axes[1].set_ylabel(r"raw $\Omega$ ($\AA^2$)")
    axes[1].set_xlabel(r"periodic shortest-path distance K$\rightarrow$K' ($\AA^{-1}$)")
    axes[1].legend()
    axes[1].grid(alpha=0.25)
    fig.suptitle(SHORTEST_K_KP_LINE_TITLE)
    output = output_dir / "wavecar_fukui_line_K_Kp.png"
    fig.savefig(output, dpi=180)
    plt.close(fig)
    return output


def plot_transport_csv(csv_path: Path, output_path: Path) -> Path:
    plt = _matplotlib_pyplot()
    with csv_path.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    mu = np.asarray([float(row["mu_eV"]) for row in rows])
    total = np.asarray([float(row["delta_sigma_e2_over_h"]) for row in rows])
    kpart = np.asarray([float(row["delta_sigma_K_e2_over_h"]) for row in rows])
    kppart = np.asarray([float(row["delta_sigma_Kp_e2_over_h"]) for row in rows])
    contrast = np.asarray(
        [float(row["delta_sigma_valley_contrast_e2_over_h"]) for row in rows]
    )
    invalid = np.asarray([row["quality"] != "PASS" for row in rows])
    fig, axis = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    series = (
        (total, "black", "net increment", "-"),
        (kpart, "tab:blue", "K increment", "-"),
        (kppart, "tab:orange", "K' increment", "-"),
        (contrast, "tab:green", "K - K' contrast", "--"),
    )
    for values, color, label, style in series:
        axis.plot(
            mu,
            np.where(~invalid, values, np.nan),
            color=color,
            linestyle=style,
            label=label,
        )
        if np.any(invalid):
            axis.plot(
                mu,
                np.where(invalid, values, np.nan),
                color=color,
                linestyle=":",
                alpha=0.45,
            )
    if np.any(invalid):
        ymin, ymax = axis.get_ylim()
        axis.scatter(
            mu[invalid],
            np.full(np.count_nonzero(invalid), ymin + 0.03 * (ymax - ymin)),
            marker="x",
            color="red",
            s=24,
            label="INVALID active-cell quality",
            zorder=5,
        )
        axis.set_title("UNVALIDATED DIAGNOSTIC: red x marks rejected chemical potentials")
    else:
        axis.set_title("Validated T=0 cumulative-subspace transport")
    axis.set_xlabel("chemical potential (eV)")
    axis.set_ylabel(r"$\Delta\sigma_{xy}$ ($e^2/h$)")
    axis.grid(alpha=0.25)
    axis.legend(ncol=2)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    return output_path


def plot_full_transport_csv(
    csv_path: Path,
    output_path: Path,
    diagnostics_path: Path | None = None,
) -> Path:
    """Plot absolute and full-bundle-relative sigma(mu) with quality marks."""

    plt = _matplotlib_pyplot()
    with csv_path.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError("full transport CSV contains no rows")

    def optional_float(row: dict[str, str], field: str) -> float:
        value = row[field].strip()
        return float(value) if value else math.nan

    mu = np.asarray([float(row["mu_eV"]) for row in rows])
    total = np.asarray(
        [optional_float(row, "sigma_total_e2_over_h") for row in rows]
    )
    kpart = np.asarray(
        [optional_float(row, "sigma_K_e2_over_h") for row in rows]
    )
    kppart = np.asarray(
        [optional_float(row, "sigma_Kp_e2_over_h") for row in rows]
    )
    contrast = np.asarray(
        [
            optional_float(row, "sigma_valley_contrast_e2_over_h")
            for row in rows
        ]
    )
    relative_total = np.asarray(
        [
            optional_float(
                row, "sigma_relative_to_full_max_bundle_e2_over_h"
            )
            for row in rows
        ]
    )
    relative_k = np.asarray(
        [
            optional_float(
                row, "sigma_K_relative_to_full_max_bundle_e2_over_h"
            )
            for row in rows
        ]
    )
    relative_kp = np.asarray(
        [
            optional_float(
                row, "sigma_Kp_relative_to_full_max_bundle_e2_over_h"
            )
            for row in rows
        ]
    )
    relative_contrast = np.asarray(
        [
            optional_float(
                row,
                "sigma_valley_contrast_relative_to_full_max_bundle_e2_over_h",
            )
            for row in rows
        ]
    )
    invalid = np.asarray([row["quality"] != "PASS" for row in rows])
    continuous_bad = rows[0]["continuous_scan_quality"] != "PASS"
    diagnostics_bad = False
    invalid_intervals: list[list[float]] = []
    if diagnostics_path is not None and diagnostics_path.exists():
        with diagnostics_path.open(encoding="utf-8") as handle:
            payload = json.load(handle)
        diagnostics_bad = not bool(payload.get("validated", False))
        invalid_intervals = payload.get("merged_invalid_mu_intervals_ev", [])

    fig, axes = plt.subplots(
        2, 1, figsize=(10, 9), sharex=True, constrained_layout=True
    )
    panels = (
        (
            axes[0],
            (total, kpart, kppart, contrast),
            r"absolute $\sigma_{xy}$ ($e^2/h$)",
            "Absolute conductivity of represented occupied bands",
        ),
        (
            axes[1],
            (relative_total, relative_k, relative_kp, relative_contrast),
            r"$\Delta\sigma_{xy}$ from full MAX_BAND bundle ($e^2/h$)",
            "Hole/electron change relative to the fully occupied MAX_BAND bundle",
        ),
    )
    styles = (
        ("black", "net", "-"),
        ("tab:blue", "K", "-"),
        ("tab:orange", "K'", "-"),
        ("tab:green", "K - K' contrast", "--"),
    )
    for axis, values_set, ylabel, title in panels:
        for values, (color, label, style) in zip(values_set, styles):
            axis.plot(
                mu, np.where(~invalid, values, np.nan), color=color,
                linestyle=style, label=label,
            )
            if np.any(invalid):
                axis.plot(
                    mu, np.where(invalid, values, np.nan), color=color,
                    linestyle=":", alpha=0.45,
                )
        for index, interval in enumerate(invalid_intervals):
            lo, hi = (float(interval[0]), float(interval[1]))
            label = "INVALID occupation-event interval" if index == 0 else None
            if hi > lo:
                axis.axvspan(lo, hi, color="red", alpha=0.10, label=label)
            else:
                axis.axvline(
                    lo, color="red", alpha=0.25, linewidth=1.0, label=label
                )
        if np.any(invalid):
            ymin, ymax = axis.get_ylim()
            axis.scatter(
                mu[invalid],
                np.full(
                    np.count_nonzero(invalid), ymin + 0.03 * (ymax - ymin)
                ),
                marker="x", color="red", s=18,
                label="INVALID sampled mu", zorder=5,
            )
        axis.set_title(title)
        axis.set_ylabel(ylabel)
        axis.grid(alpha=0.25)
        axis.legend(ncol=2)
    axes[1].set_xlabel("chemical potential (eV)")
    if diagnostics_bad or continuous_bad or np.any(invalid):
        fig.suptitle(
            "UNVALIDATED DIAGNOSTIC: rejected points/intervals shown in red"
        )
    else:
        fig.suptitle("Validated all-band T=0 cumulative-subspace transport")
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    return output_path


def parse_map(specification: str) -> tuple[str, list[int]]:
    try:
        label, range_text = specification.split("=", 1)
        if ":" in range_text:
            start_text, stop_text = range_text.split(":", 1)
            start, stop = int(start_text), int(stop_text)
        else:
            start = stop = int(range_text)
    except ValueError as error:
        raise argparse.ArgumentTypeError("map must have LABEL=N or LABEL=N:M") from error
    if (
        re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.+-]*", label) is None
        or start < 1
        or stop < start
    ):
        raise argparse.ArgumentTypeError("map must have LABEL=N or LABEL=N:M")
    return label, list(range(start, stop + 1))


def default_map_specifications(energy_band: int) -> list[tuple[str, list[int]]]:
    if energy_band < 2:
        raise ValueError("target energy band must be at least 2")
    return [
        ("V", list(range(1, energy_band))),
        ("VC", list(range(1, energy_band + 1))),
        ("VCC", list(range(1, energy_band + 2))),
        ("C", [energy_band]),
        ("Cpair", [energy_band, energy_band + 1]),
    ]


def validate_map_specifications(
    maps: Sequence[tuple[str, Sequence[int]]],
    energy_band: int,
    nbands: int,
    require_cumulative: bool,
    require_transport_headroom: bool = False,
) -> None:
    """Validate explicit/default map bounds before any expensive overlap work."""

    if energy_band < 2:
        raise ValueError("target energy band must be at least 2")
    if energy_band + 1 > nbands:
        raise ValueError(
            "--energy-band must leave target+1 for raw-map energy/gap export"
        )
    if require_transport_headroom and energy_band + 2 > nbands:
        raise ValueError(
            "--transport-t0 requires target+2 for cumulative-subspace "
            "isolation and energy-window checks"
        )
    labels: set[str] = set()
    ranges: set[tuple[int, int]] = set()
    for label, bands in maps:
        if label in labels:
            raise ValueError(f"duplicate map label {label!r}")
        labels.add(label)
        if not bands or bands[0] < 1 or any(
            right != left + 1 for left, right in zip(bands, bands[1:])
        ):
            raise ValueError(f"map {label!r} must contain a contiguous positive band range")
        if bands[-1] > nbands:
            raise ValueError(f"map {label} exceeds WAVECAR NBANDS={nbands}")
        band_range = (bands[0], bands[-1])
        if band_range in ranges:
            raise ValueError(f"duplicate map band range {band_range}")
        ranges.add(band_range)
    if require_cumulative:
        required = {
            (1, energy_band - 1),
            (1, energy_band),
            (1, energy_band + 1),
        }
        missing = sorted(required - ranges)
        if missing:
            raise ValueError(
                "plot/transport requires cumulative maps "
                f"1:{energy_band - 1}, 1:{energy_band}, and 1:{energy_band + 1}; "
                f"missing {missing}"
            )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("wavecar", type=Path)
    parser.add_argument("--nx", type=int, required=True)
    parser.add_argument("--ny", type=int, required=True)
    parser.add_argument("--spinor-components", type=int, choices=(1, 2), default=2)
    parser.add_argument(
        "--spin",
        type=int,
        choices=(1, 2),
        default=1,
        help="ISPIN channel (SOC ISPIN=1 uses spin=1 with two spinor components)",
    )
    parser.add_argument("--energy-band", type=int, default=19)
    parser.add_argument(
        "--map",
        dest="maps",
        action="append",
        type=parse_map,
        default=None,
        help="LABEL=N or LABEL=N:M; repeat for multiple maps",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--transport-t0",
        action="store_true",
        help="evaluate guarded T=0 cumulative-subspace transport",
    )
    parser.add_argument(
        "--transport-full-t0",
        type=int,
        metavar="MAX_BAND",
        help=(
            "evaluate absolute T=0 sigma(mu) using every cumulative subspace "
            "1:n through MAX_BAND; MAX_BAND+1 is an unoccupied sentinel"
        ),
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help=(
            "write k-resolved, shortest-periodic-image K-to-K' line, and "
            "(with transport) sigma plots"
        ),
    )
    parser.add_argument(
        "--plot-domain",
        choices=("fractional", "first-bz"),
        default="fractional",
        help=(
            "k-resolved map domain: legacy fractional q1/q2 primitive cell "
            "(default), or the Cartesian Wigner-Seitz first BZ"
        ),
    )
    parser.add_argument("--mu-min", type=float, default=0.40)
    parser.add_argument("--mu-max", type=float, default=0.55)
    parser.add_argument("--mu-num", type=int, default=151)
    parser.add_argument("--valley-k", type=parse_fractional_vector)
    parser.add_argument("--valley-kp", type=parse_fractional_vector)
    parser.add_argument(
        "--valley-radius",
        type=float,
        help="optional K/K' mask radius in inverse Angstrom; default is a Voronoi partition",
    )
    parser.add_argument("--min-link-sv", type=float, default=1.0e-3)
    parser.add_argument("--min-neighbor-gap-ev", type=float, default=1.0e-5)
    parser.add_argument(
        "--min-pw-coverage",
        type=float,
        default=0.90,
        help="minimum common plane-wave fraction for every nearest-neighbor link",
    )
    parser.add_argument("--max-abs-phi", type=float, default=0.8 * np.pi)
    parser.add_argument(
        "--allow-invalid-transport",
        action="store_true",
        help="write an explicitly INVALID diagnostic scan instead of refusing it",
    )
    args = parser.parse_args()

    full_only_no_legacy_maps = (
        args.transport_full_t0 is not None
        and args.maps is None
        and not args.plot
    )
    if not full_only_no_legacy_maps and args.energy_band < 2:
        parser.error("target energy band must be at least 2")
    if args.transport_t0 and args.transport_full_t0 is not None:
        parser.error("--transport-t0 and --transport-full-t0 are mutually exclusive")
    if args.transport_full_t0 is not None and args.transport_full_t0 < 1:
        parser.error("--transport-full-t0 MAX_BAND must be positive")
    if not (math.isfinite(args.mu_min) and math.isfinite(args.mu_max)):
        parser.error("--mu-min and --mu-max must be finite")
    if args.valley_radius is not None and not math.isfinite(args.valley_radius):
        parser.error("--valley-radius must be finite")
    for option, value in (
        ("--min-link-sv", args.min_link_sv),
        ("--min-neighbor-gap-ev", args.min_neighbor_gap_ev),
        ("--min-pw-coverage", args.min_pw_coverage),
        ("--max-abs-phi", args.max_abs_phi),
    ):
        if not math.isfinite(value):
            parser.error(f"{option} must be finite")
    if args.min_link_sv < 0.0 or args.min_neighbor_gap_ev < 0.0:
        parser.error("--min-link-sv and --min-neighbor-gap-ev must be nonnegative")
    if not (0.0 <= args.min_pw_coverage <= 1.0):
        parser.error("--min-pw-coverage must be in [0, 1]")
    if not (0.0 < args.max_abs_phi < np.pi):
        parser.error("--max-abs-phi must be in (0, pi)")
    try:
        if args.maps is not None:
            maps = args.maps
        elif full_only_no_legacy_maps:
            # Full transport builds its cumulative maps in one batched pass.
            # Avoid repeating the five legacy/raw maps unless explicitly asked.
            maps = []
        else:
            maps = default_map_specifications(args.energy_band)
        output_names = planned_output_names(
            maps, args.plot, args.transport_t0,
            args.transport_full_t0 is not None,
        )
        validate_input_output_collision(
            args.wavecar,
            args.output_dir,
            output_names,
        )
        remove_planned_outputs(args.output_dir, output_names)
    except ValueError as error:
        parser.error(str(error))
    wavecar = Wavecar(
        args.wavecar, spinor_components=args.spinor_components, spin=args.spin
    )
    grid = infer_uniform_grid(wavecar.kpoints, args.nx, args.ny)
    try:
        if not full_only_no_legacy_maps:
            validate_map_specifications(
                maps,
                args.energy_band,
                wavecar.header.nbands,
                args.plot or args.transport_t0,
                args.transport_t0,
            )
        if (
            args.transport_full_t0 is not None
            and args.transport_full_t0 >= wavecar.header.nbands
        ):
            raise ValueError(
                "--transport-full-t0 MAX_BAND must be below WAVECAR NBANDS "
                "so MAX_BAND+1 can serve as the upper occupation sentinel"
            )
        if args.plot or args.transport_t0 or args.transport_full_t0 is not None:
            if args.valley_k is None or args.valley_kp is None:
                raise ValueError(
                    "plot/transport requires --valley-k and --valley-kp"
                )
            validate_valley_definition(
                args.valley_k,
                args.valley_kp,
                wavecar.header.reciprocal,
                args.valley_radius,
            )
    except ValueError as error:
        parser.error(str(error))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    summaries: list[dict[str, object]] = []
    map_results: dict[str, tuple[list[int], np.ndarray, LinkSet, LinkSet]] = {}
    for label, bands in maps:
        print(f"computing {label}: bands {bands[0]}:{bands[-1]}", flush=True)
        flux, xlinks, ylinks = compute_band_map(wavecar, grid, bands)
        write_table(
            args.output_dir / f"fukui_{label}.csv",
            wavecar,
            grid,
            label,
            bands,
            flux,
            xlinks,
            ylinks,
            args.energy_band,
        )
        summaries.append(
            map_summary(
                label, bands, flux, xlinks, ylinks, args.min_pw_coverage
            )
        )
        map_results[label] = (bands, flux, xlinks, ylinks)

    energy_diagnostics: dict[str, object] | None = None
    if not full_only_no_legacy_maps:
        band = args.energy_band - 1
        gap_below = wavecar.energies[:, band] - wavecar.energies[:, band - 1]
        gap_above = wavecar.energies[:, band + 1] - wavecar.energies[:, band]
        energy_diagnostics = {
            "band": args.energy_band,
            "energy_min_ev": float(np.min(wavecar.energies[:, band])),
            "energy_max_ev": float(np.max(wavecar.energies[:, band])),
            "min_gap_below_ev": float(np.min(gap_below)),
            "min_gap_below_k_index": int(np.argmin(gap_below) + 1),
            "min_gap_above_ev": float(np.min(gap_above)),
            "min_gap_above_k_index": int(np.argmin(gap_above) + 1),
        }
    diagnostics = {
        # Keep shareable diagnostics independent of the caller's local
        # directory layout while retaining the input filename as provenance.
        "wavecar": args.wavecar.name,
        "logical_recl": wavecar.header.logical_recl,
        "physical_record_bytes": wavecar.header.stride_bytes,
        "legacy_four_byte_word_recl": (
            wavecar.header.stride_bytes == 4 * wavecar.header.logical_recl
        ),
        "ispin": wavecar.header.ispin,
        "spin_channel": wavecar.spin,
        "spinor_components": wavecar.spinor_components,
        "rtag": wavecar.header.rtag,
        "nkpoints": wavecar.header.nkpoints,
        "nbands": wavecar.header.nbands,
        "encut_ev": wavecar.header.encut_ev,
        "grid": [grid.nx, grid.ny],
        "grid_offset_fractional": grid.offset.tolist(),
        "plot_domain": args.plot_domain if args.plot else None,
        "lattice_angstrom": wavecar.header.lattice.tolist(),
        "reciprocal_inv_angstrom": wavecar.header.reciprocal.tolist(),
        "nbmax": list(wavecar.nbmax),
        "energy_diagnostics": energy_diagnostics,
        "maps": summaries,
        "plane_wave_coverage_validated": all(
            bool(item["plane_wave_coverage_pass"]) for item in summaries
        ),
        "conventions": {
            "loop": "k -> k+dk1 -> k+dk1+dk2 -> k+dk2 -> k",
            "flux": "-Arg(product of link determinants), principal branch",
            "curvature": "flux / (|b1 x b2|/(nx*ny))",
            "chern": "sum(flux)/(2*pi)",
        },
    }
    write_strict_json(args.output_dir / "diagnostics.json", diagnostics)
    if not diagnostics["plane_wave_coverage_validated"]:
        raise RuntimeError(
            "common plane-wave coverage is below --min-pw-coverage; raw files and "
            "diagnostics were written, but plotting and transport were refused"
        )
    if args.plot:
        if args.valley_k is None or args.valley_kp is None:
            parser.error("--plot requires --valley-k and --valley-kp")
        plot_k_resolved_maps(
            args.output_dir,
            wavecar,
            grid,
            map_results,
            args.energy_band,
            args.valley_k,
            args.valley_kp,
            args.plot_domain,
        )
        plot_k_to_kp_line(
            args.output_dir,
            wavecar,
            grid,
            map_results,
            args.energy_band,
            args.valley_k,
            args.valley_kp,
        )
    if args.transport_t0:
        if args.valley_k is None or args.valley_kp is None:
            parser.error("--transport-t0 requires --valley-k and --valley-kp")
        if args.mu_num < 2 or args.mu_max <= args.mu_min:
            parser.error("transport requires mu_num>=2 and mu_max>mu_min")
        write_t0_transport(
            args.output_dir,
            wavecar,
            grid,
            map_results,
            args.energy_band,
            args.mu_min,
            args.mu_max,
            args.mu_num,
            args.valley_k,
            args.valley_kp,
            args.valley_radius,
            args.min_link_sv,
            args.min_neighbor_gap_ev,
            args.max_abs_phi,
            args.allow_invalid_transport,
        )
        if args.plot:
            plot_transport_csv(
                args.output_dir / "transport_t0.csv",
                args.output_dir / "wavecar_fukui_sigma_mu.png",
            )
    if args.transport_full_t0 is not None:
        if args.valley_k is None or args.valley_kp is None:
            parser.error("--transport-full-t0 requires --valley-k and --valley-kp")
        if args.mu_num < 2 or args.mu_max <= args.mu_min:
            parser.error("full transport requires mu_num>=2 and mu_max>mu_min")
        write_full_t0_transport(
            args.output_dir,
            wavecar,
            grid,
            args.transport_full_t0,
            args.mu_min,
            args.mu_max,
            args.mu_num,
            args.valley_k,
            args.valley_kp,
            args.valley_radius,
            args.min_link_sv,
            args.min_neighbor_gap_ev,
            args.min_pw_coverage,
            args.max_abs_phi,
            args.allow_invalid_transport,
        )
        if args.plot:
            plot_full_transport_csv(
                args.output_dir / "transport_full_t0.csv",
                args.output_dir / "wavecar_fukui_sigma_full_mu.png",
                args.output_dir / "transport_full_t0_diagnostics.json",
            )
    print(json.dumps(json_safe(diagnostics), indent=2, allow_nan=False))


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, RuntimeError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(2)
