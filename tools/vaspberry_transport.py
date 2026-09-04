#!/usr/bin/env python3
"""Plot VASPBERRY curvature and calculate intrinsic 2D Hall transport.

The transport implemented here is the Berry-curvature (Fermi-sea) term

    C(mu,T) = sum_n integral_BZ f(E_nk-mu,T) Omega_n(k) d^2k / (2*pi)
    sigma_xy(mu,T)/(e^2/h) = -C(mu,T)

for a two-dimensional, uniformly sampled *full* Brillouin zone. For Fukui
plaquettes, ``Omega*dS`` is weighted by the average Fermi occupation of its four
vertices. A line-mode calculation, including legacy Kubo path output, is useful
for plotting but is deliberately rejected for transport.

This is a standalone validation/post-processing tool.  It does not modify the
VASPBERRY executable or its legacy output files.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


__version__ = "1.2.0"
KB_EV_PER_K = 8.617333262145e-5
TWO_PI = 2.0 * math.pi


@dataclass
class CurvatureData:
    cart: np.ndarray
    frac: np.ndarray
    omega: np.ndarray
    metadata: dict[str, object] = field(default_factory=dict)
    energy: np.ndarray | None = None
    vertex_energies: np.ndarray | None = None
    source: Path | None = None

    def __post_init__(self) -> None:
        self.cart = np.asarray(self.cart, dtype=float)
        self.frac = np.asarray(self.frac, dtype=float)
        self.omega = np.asarray(self.omega, dtype=float)
        if self.energy is not None:
            self.energy = np.asarray(self.energy, dtype=float)
        if self.vertex_energies is not None:
            self.vertex_energies = np.asarray(self.vertex_energies, dtype=float)
        n = len(self.omega)
        if n == 0:
            raise ValueError("curvature data must contain at least one point")
        if self.cart.shape != (n, 3) or self.frac.shape != (n, 3):
            raise ValueError("cart and frac coordinates must both have shape (N, 3)")
        if self.energy is not None and self.energy.shape != (n,):
            raise ValueError("energy must have shape (N,)")
        if self.vertex_energies is not None and self.vertex_energies.shape != (n, 4):
            raise ValueError("vertex_energies must have shape (N, 4)")
        arrays = {"cart coordinates": self.cart, "fractional coordinates": self.frac,
                  "Berry curvature": self.omega}
        if self.energy is not None:
            arrays["band energy"] = self.energy
        if self.vertex_energies is not None:
            arrays["Fukui vertex energies"] = self.vertex_energies
        for label, values in arrays.items():
            if not np.all(np.isfinite(values)):
                raise ValueError(f"{label} contains NaN or infinity")


@dataclass
class EigenvalData:
    frac: np.ndarray
    energies: np.ndarray
    occupations: np.ndarray
    weights: np.ndarray
    nelect: int
    nk: int
    nband: int
    source: Path

    def __post_init__(self) -> None:
        self.frac = np.asarray(self.frac, dtype=float)
        self.energies = np.asarray(self.energies, dtype=float)
        self.occupations = np.asarray(self.occupations, dtype=float)
        self.weights = np.asarray(self.weights, dtype=float)
        if self.nk <= 0 or self.nband <= 0:
            raise ValueError("EIGENVAL NK and NBANDS must be positive")
        if self.frac.shape != (self.nk, 3):
            raise ValueError("EIGENVAL fractional coordinates must have shape (NK, 3)")
        if self.energies.shape != (self.nk, self.nband):
            raise ValueError("EIGENVAL energies must have shape (NK, NBANDS)")
        if self.occupations.shape != (self.nk, self.nband):
            raise ValueError("EIGENVAL occupations must have shape (NK, NBANDS)")
        if self.weights.shape != (self.nk,):
            raise ValueError("EIGENVAL weights must have shape (NK,)")
        for label, values in {
            "EIGENVAL coordinates": self.frac,
            "EIGENVAL energies": self.energies,
            "EIGENVAL occupations": self.occupations,
            "EIGENVAL weights": self.weights,
        }.items():
            if not np.all(np.isfinite(values)):
                raise ValueError(f"{label} contains NaN or infinity")


def validate_finite_curvature(data: CurvatureData) -> None:
    """Defense-in-depth for arrays mutated after dataclass construction."""

    arrays = {"cart coordinates": data.cart, "fractional coordinates": data.frac,
              "Berry curvature": data.omega}
    if data.energy is not None:
        arrays["band energy"] = data.energy
    if data.vertex_energies is not None:
        arrays["Fukui vertex energies"] = data.vertex_energies
    for label, values in arrays.items():
        if not np.all(np.isfinite(values)):
            raise ValueError(f"{label} contains NaN or infinity")


def _header_lines(path: Path) -> list[str]:
    lines: list[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.lstrip().startswith("#"):
                lines.append(line.rstrip())
    return lines


def _parse_metadata(headers: Sequence[str], source: Path | None = None) -> dict[str, object]:
    metadata: dict[str, object] = {"headers": list(headers)}
    patterns = {
        "nk": r"NKPOINT\s*:\s*(\d+)",
        "nbands": r"NBANDS\s*:\s*(\d+)",
        "dk_area": r"dk\^2\s*=.*?=\s*([-+0-9.Ee]+)",
        "chern_reported": r"Chern Number\s*=\s*([-+0-9.Ee]+)",
    }
    joined = "\n".join(headers)
    for key, pattern in patterns.items():
        matches = re.findall(pattern, joined)
        if matches:
            value = matches[-1]
            metadata[key] = int(value) if key in {"nk", "nbands"} else float(value)

    grid_match = re.search(r"K-GRID\s*:\s*(\d+)\s*X\s*(\d+)", joined)
    if grid_match:
        metadata["k_grid"] = (int(grid_match.group(1)), int(grid_match.group(2)))

    ispin_match = re.search(r"ISPIN\s*:\s*(\d+)", joined)
    if ispin_match:
        metadata["ispin"] = int(ispin_match.group(1))
    lsorbit_match = re.search(r"LSORBIT\s*=\s*\.(TRUE|FALSE)\.", joined, re.IGNORECASE)
    if lsorbit_match:
        metadata["lsorbit"] = lsorbit_match.group(1).upper() == "TRUE"

    explicit_spin = re.search(r"(?:^|\n)#\s*SPIN\s*:\s*(UP|DN|DOWN)\b", joined, re.IGNORECASE)
    channel = None
    if explicit_spin:
        channel = 1 if explicit_spin.group(1).upper() == "UP" else 2
    if source is not None:
        filename_match = re.search(r"\.(UP|DN|DOWN)\.DAT$", source.name, re.IGNORECASE)
        if filename_match:
            filename_channel = 1 if filename_match.group(1).upper() == "UP" else 2
            if channel is not None and channel != filename_channel:
                raise ValueError(f"{source}: spin label in header conflicts with filename")
            channel = filename_channel
    if channel is not None:
        if metadata.get("ispin") not in (None, 2):
            raise ValueError(f"{source or '<data>'}: UP/DN channel requires ISPIN=2")
        metadata["spin_channel"] = channel
        metadata["spin_label"] = "UP" if channel == 1 else "DN"

    for index in (1, 2, 3):
        match = re.search(
            rf"RECIVEC B{index}(?:\s*\(A\^-1\))?\s*:\s*"
            r"([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)",
            joined,
        )
        if match:
            metadata[f"b{index}"] = np.array([float(v) for v in match.groups()])

    single = re.search(r"Chern Number for the BAND\s*:\s*(\d+)", joined)
    multi = re.search(r"Chern Number for the BANDS\s*:\s*(\d+)\s*-\s*(\d+)", joined)
    kubo = re.search(r"Berry curvature using kubo for BAND\s+index n\s+(\d+)", joined)
    if kubo:
        metadata["band"] = int(kubo.group(1))
        metadata["band_mode"] = "single-kubo"
    elif single:
        metadata["band"] = int(single.group(1))
        metadata["band_mode"] = "single-fukui"
    elif multi:
        metadata["band_range"] = (int(multi.group(1)), int(multi.group(2)))
        metadata["band_mode"] = "manifold"
    return metadata


def read_vaspberry(path: str | Path, energy_column: int | None = None) -> CurvatureData:
    """Read legacy seven-column VASPBERRY curvature output.

    Expected data columns are cartesian kx, ky, kz, Omega, reciprocal/fractional
    k1, k2, k3.  An optional energy column can be selected with a zero-based
    index; otherwise an eighth column is treated as energy when present.
    """

    source = Path(path)
    raw = np.loadtxt(source, comments="#", ndmin=2)
    if raw.shape[1] < 7:
        raise ValueError(f"{source}: expected at least 7 numeric columns, got {raw.shape[1]}")
    if energy_column is None and raw.shape[1] >= 8:
        energy_column = 7
    if energy_column is not None:
        if not isinstance(energy_column, (int, np.integer)):
            raise ValueError("energy_column must be a zero-based integer column index")
        if energy_column < 7 or energy_column >= raw.shape[1]:
            available = (
                f"7..{raw.shape[1] - 1}" if raw.shape[1] >= 8
                else "none (the file has only the seven required columns)"
            )
            raise ValueError(
                f"{source}: energy_column {energy_column} is not an available extra "
                f"zero-based energy column; valid range: {available}"
            )
    energy = None if energy_column is None else raw[:, energy_column]
    return CurvatureData(
        cart=raw[:, 0:3],
        omega=raw[:, 3],
        frac=raw[:, 4:7],
        energy=energy,
        metadata=_parse_metadata(_header_lines(source), source),
        source=source,
    )


def read_eigenval(path: str | Path, spin: int = 1) -> EigenvalData:
    """Read standard VASP EIGENVAL energies and occupations.

    ``spin`` is one-based.  For a non-spin-polarized/SOC EIGENVAL use spin=1.
    For collinear ISPIN=2 EIGENVAL, spin=1/2 selects the corresponding column.
    """

    source = Path(path)
    lines = source.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 7:
        raise ValueError(f"{source}: too short to be an EIGENVAL file")
    try:
        nelect, nk, nband = (int(float(v)) for v in lines[5].split()[:3])
    except (ValueError, IndexError) as exc:
        raise ValueError(f"{source}: cannot parse NELECT/NK/NBANDS on line 6") from exc

    coords = np.zeros((nk, 3), dtype=float)
    energies = np.zeros((nk, nband), dtype=float)
    occupations = np.zeros((nk, nband), dtype=float)
    weights = np.zeros(nk, dtype=float)
    cursor = 6
    for ik in range(nk):
        while cursor < len(lines) and not lines[cursor].strip():
            cursor += 1
        if cursor >= len(lines):
            raise ValueError(f"{source}: missing k point {ik + 1}")
        k_tokens = lines[cursor].split()
        cursor += 1
        if len(k_tokens) < 4:
            raise ValueError(f"{source}: malformed k-point record {ik + 1}")
        coords[ik] = [float(v) for v in k_tokens[:3]]
        weights[ik] = float(k_tokens[3])

        for ib in range(nband):
            if cursor >= len(lines):
                raise ValueError(f"{source}: missing band {ib + 1} at k point {ik + 1}")
            tokens = lines[cursor].split()
            cursor += 1
            values = [float(v) for v in tokens[1:]]
            if len(values) >= 4:  # collinear: E_up E_dn occ_up occ_dn
                if spin not in (1, 2):
                    raise ValueError("spin must be 1 or 2 for a collinear EIGENVAL")
                energies[ik, ib] = values[spin - 1]
                occupations[ik, ib] = values[spin + 1]
            elif len(values) >= 2:  # non-spin-polarized or SOC: E occ
                if spin != 1:
                    raise ValueError("spin=2 requested for a single-energy-column EIGENVAL")
                energies[ik, ib], occupations[ik, ib] = values[:2]
            else:
                raise ValueError(f"{source}: malformed band record at k={ik + 1}, band={ib + 1}")

    return EigenvalData(coords, energies, occupations, weights, nelect, nk, nband, source)


def _wrap_fractional(delta: np.ndarray) -> np.ndarray:
    return delta - np.rint(delta)


def validate_spin_channels(
    datasets: Sequence[CurvatureData],
    requested_spin: int,
    require_collinear_label: bool = True,
) -> None:
    """Require curvature spin provenance to agree with one EIGENVAL channel."""

    if requested_spin not in (1, 2):
        raise ValueError("spin must be 1 (UP) or 2 (DN)")
    channels = [data.metadata.get("spin_channel") for data in datasets]
    labeled = {int(channel) for channel in channels if channel is not None}
    if len(labeled) > 1:
        raise ValueError("UP and DN curvature files cannot be mixed in one sigma calculation")
    if labeled and any(channel is None for channel in channels):
        raise ValueError("spin-labeled and unlabeled curvature files cannot be mixed")
    if labeled and next(iter(labeled)) != requested_spin:
        label = "UP" if next(iter(labeled)) == 1 else "DN"
        raise ValueError(
            f"curvature file is {label}, but global --spin selects channel {requested_spin}"
        )
    for data, channel in zip(datasets, channels, strict=True):
        ispin = data.metadata.get("ispin")
        if ispin is not None and int(ispin) not in (1, 2):
            raise ValueError(f"{data.source or '<data>'}: unsupported ISPIN={ispin}")
        if ispin == 1 and requested_spin != 1:
            raise ValueError(
                f"{data.source or '<data>'}: ISPIN=1 curvature cannot use global --spin 2"
            )
        if require_collinear_label and ispin == 2 and channel is None:
            raise ValueError(
                f"{data.source or '<data>'}: ISPIN=2 curvature channel is ambiguous; "
                "retain a .UP/.DN filename or an explicit SPIN header"
            )


def validate_requested_band(curvature: CurvatureData, band: int) -> None:
    """Prevent a CLI band number from silently overriding file metadata."""

    header_band = curvature.metadata.get("band")
    if header_band is not None and int(header_band) != band:
        raise ValueError(
            f"{curvature.source or '<data>'}: requested band {band} conflicts with "
            f"BERRYCURV header band {int(header_band)}"
        )


def validate_eigenval_band_count(curvature: CurvatureData, eigenval: EigenvalData) -> None:
    header_nbands = curvature.metadata.get("nbands")
    if header_nbands is not None and int(header_nbands) != eigenval.nband:
        raise ValueError(
            f"{curvature.source or '<data>'}: BERRYCURV header NBANDS={int(header_nbands)} "
            f"does not match EIGENVAL NBANDS={eigenval.nband}"
        )


def validate_unique_active_bands(bands: Sequence[int]) -> None:
    """Reject repeated active-band inputs before they can be double-counted."""

    seen: set[int] = set()
    duplicates: set[int] = set()
    for band in bands:
        if band in seen:
            duplicates.add(band)
        seen.add(band)
    if duplicates:
        raise ValueError(
            f"duplicate active band inputs {sorted(duplicates)} would double-count transport"
        )


def validate_core_chern(core_chern: float, active_bands: Sequence[int]) -> int:
    """Validate the integer Chern baseline assigned to omitted occupied bands."""

    if not np.isfinite(core_chern):
        raise ValueError("core Chern baseline must be finite")
    nearest = int(np.rint(core_chern))
    if abs(core_chern - nearest) > 5.0e-3:
        raise ValueError(
            f"core Chern baseline {core_chern:.9g} is not within 0.005 of an integer "
            f"(nearest {nearest})"
        )
    if active_bands and min(active_bands) == 1 and core_chern != 0.0:
        raise ValueError(
            "a nonzero core Chern baseline is invalid when the active window starts at band 1"
        )
    return nearest


def validate_legacy_sigma_energy_window(
    eigenval: EigenvalData,
    active_bands: Sequence[int],
    mu_values: np.ndarray,
    temperature: float,
    occupation_tolerance: float = 1.0e-8,
    core_chern: float = 0.0,
) -> dict[str, object]:
    """Require a complete finite band window for legacy EIGENVAL transport.

    Bands below the active window are represented only by ``core_chern`` and
    therefore must be fully occupied throughout the requested (mu, T) range.
    The first band above the active window must remain unoccupied. Active bands
    may be partially occupied, but must form one consecutive window.
    """

    bands = [int(band) for band in active_bands]
    if not bands:
        raise ValueError("at least one active band is required")
    validate_unique_active_bands(bands)
    ordered = sorted(bands)
    if ordered[0] < 1 or ordered[-1] > eigenval.nband:
        raise ValueError(
            f"active bands {ordered} are outside EIGENVAL range 1..{eigenval.nband}"
        )
    expected = list(range(ordered[0], ordered[-1] + 1))
    if ordered != expected:
        raise ValueError(
            f"active bands must be consecutive; got {ordered}, expected the complete "
            f"window {expected}"
        )
    if not np.isfinite(occupation_tolerance) or not 0.0 <= occupation_tolerance < 0.5:
        raise ValueError("occupation tolerance must be finite and in the range [0, 0.5)")
    mu = np.asarray(mu_values, dtype=float)
    if mu.ndim != 1 or len(mu) == 0 or not np.all(np.isfinite(mu)):
        raise ValueError("chemical-potential grid must be a nonempty finite one-dimensional array")
    if not np.isfinite(temperature) or temperature < 0.0:
        raise ValueError("temperature must be finite and nonnegative")

    nearest_core_chern = validate_core_chern(core_chern, ordered)
    first, last = ordered[0], ordered[-1]
    if last >= eigenval.nband:
        raise ValueError(
            f"active window ends at EIGENVAL's final band {last}; band {last + 1} is "
            "required to verify that omitted higher bands remain unoccupied"
        )

    mu_min = float(np.min(mu))
    mu_max = float(np.max(mu))
    core_bands = list(range(1, first))
    minimum_core_occupation: float | None = None
    if core_bands:
        core_occupation = fermi_dirac(
            eigenval.energies[:, : first - 1], mu_min, temperature
        )
        minimum_core_occupation = float(np.min(core_occupation))
        if minimum_core_occupation < 1.0 - occupation_tolerance:
            location = np.unravel_index(int(np.argmin(core_occupation)), core_occupation.shape)
            failed_band = int(location[1]) + 1
            failed_energy = float(eigenval.energies[location[0], location[1]])
            raise ValueError(
                f"omitted core band {failed_band} is not fully occupied over the requested "
                f"window: min f={minimum_core_occupation:.9g} at E={failed_energy:.9g} eV, "
                f"mu_min={mu_min:.9g} eV, T={temperature:.9g} K; require f >= "
                f"{1.0 - occupation_tolerance:.9g}"
            )

    next_band = last + 1
    next_occupation = fermi_dirac(
        eigenval.energies[:, next_band - 1], mu_max, temperature
    )
    maximum_next_occupation = float(np.max(next_occupation))
    if maximum_next_occupation > occupation_tolerance:
        failed_k = int(np.argmax(next_occupation))
        failed_energy = float(eigenval.energies[failed_k, next_band - 1])
        raise ValueError(
            f"unrepresented band {next_band} is occupied over the requested window: "
            f"max f={maximum_next_occupation:.9g} at E={failed_energy:.9g} eV, "
            f"mu_max={mu_max:.9g} eV, T={temperature:.9g} K; require f <= "
            f"{occupation_tolerance:.9g}"
        )

    return {
        "validated": True,
        "active_bands": ordered,
        "active_window": [first, last],
        "core_bands": core_bands,
        "next_unrepresented_band": next_band,
        "occupation_tolerance": occupation_tolerance,
        "minimum_core_occupation": minimum_core_occupation,
        "maximum_next_band_occupation": maximum_next_occupation,
        "core_chern_input": core_chern,
        "core_chern_nearest_integer": nearest_core_chern,
        "mu_range_eV": [mu_min, mu_max],
        "temperature_K": temperature,
        "zero_temperature_convention": "E <= mu is occupied",
    }


def _canonical_xy(frac: np.ndarray, tolerance: float) -> np.ndarray:
    canonical = np.mod(np.asarray(frac, dtype=float)[:, :2], 1.0)
    canonical[np.isclose(canonical, 1.0, atol=tolerance, rtol=0.0)] = 0.0
    canonical[np.isclose(canonical, 0.0, atol=tolerance, rtol=0.0)] = 0.0
    return canonical


def _xy_keys(frac: np.ndarray, tolerance: float) -> np.ndarray:
    return np.rint(_canonical_xy(frac, tolerance) / tolerance).astype(np.int64)


@dataclass
class UniformGrid:
    nx: int
    ny: int
    keys: np.ndarray
    x_keys: np.ndarray
    y_keys: np.ndarray
    index_by_xy: dict[tuple[int, int], int]


def uniform_grid(
    frac: np.ndarray,
    tolerance: float = 2.0e-5,
    label: str = "k mesh",
) -> UniformGrid:
    """Validate and index a complete uniform product mesh on a 2D torus."""

    frac_array = np.asarray(frac, dtype=float)
    if np.max(np.abs(_wrap_fractional(frac_array[:, 2] - frac_array[0, 2]))) > tolerance:
        raise ValueError(f"{label} spans multiple kz planes; a two-dimensional mesh is required")
    keys = _xy_keys(frac_array, tolerance)
    if len(np.unique(keys, axis=0)) != len(keys):
        raise ValueError(f"{label} contains duplicate points modulo reciprocal lattice vectors")
    x_keys = np.unique(keys[:, 0])
    y_keys = np.unique(keys[:, 1])
    nx, ny = len(x_keys), len(y_keys)
    if nx * ny != len(keys):
        raise ValueError(
            f"{label} is not a complete full-BZ product mesh "
            f"(N={len(keys)}, unique grid={nx}x{ny})"
        )
    for values, axis in ((x_keys, "kx"), (y_keys, "ky")):
        coordinates = np.sort(values.astype(float) * tolerance)
        gaps = np.diff(np.r_[coordinates, coordinates[0] + 1.0])
        expected = 1.0 / len(values)
        if not np.allclose(gaps, expected, atol=5.0 * tolerance, rtol=0.0):
            raise ValueError(f"{label} {axis} coordinates are not uniformly spaced")
    mapping = {(int(key[0]), int(key[1])): i for i, key in enumerate(keys)}
    return UniformGrid(nx, ny, keys, x_keys, y_keys, mapping)


def deduplicate_reciprocal_replicas(
    data: CurvatureData,
    tolerance: float = 2.0e-5,
    omega_atol: float = 5.0e-5,
    omega_rtol: float = 1.0e-10,
) -> CurvatureData:
    """Collapse extended-BZ translations after verifying replica agreement.

    Legacy VASPBERRY normally prints nine reciprocal translations. Grouping is
    performed modulo integer reciprocal lattice vectors, not by row blocks,
    because legacy output may be globally sorted after extension.
    """

    kz_delta = _wrap_fractional(data.frac[:, 2] - data.frac[0, 2])
    if np.max(np.abs(kz_delta)) > tolerance:
        raise ValueError("curvature data contains multiple inequivalent kz planes")
    canonical_kz = float(np.mod(data.frac[0, 2], 1.0))
    if math.isclose(canonical_kz, 1.0, abs_tol=tolerance):
        canonical_kz = 0.0
    canonical = _canonical_xy(data.frac, tolerance)
    keys = _xy_keys(data.frac, tolerance)
    unique_keys, inverse, counts = np.unique(keys, axis=0, return_inverse=True, return_counts=True)
    if len(set(int(value) for value in counts)) != 1:
        raise ValueError(
            f"{data.source or '<data>'}: reciprocal replicas have inconsistent multiplicities "
            f"{sorted(set(int(value) for value in counts))}"
        )
    expected_nk = data.metadata.get("nk")
    if isinstance(expected_nk, int) and expected_nk != len(unique_keys):
        raise ValueError(
            f"{data.source or '<data>'}: header NKPOINT={expected_nk}, but "
            f"{len(unique_keys)} unique points remain modulo reciprocal vectors"
        )

    omega = np.empty(len(unique_keys), dtype=float)
    frac = np.zeros((len(unique_keys), 3), dtype=float)
    energy = None if data.energy is None else np.empty(len(unique_keys), dtype=float)
    max_spread = 0.0
    for group in range(len(unique_keys)):
        members = np.flatnonzero(inverse == group)
        values = data.omega[members]
        spread = float(np.ptp(values))
        max_spread = max(max_spread, spread)
        if not np.allclose(values, values[0], atol=omega_atol, rtol=omega_rtol):
            raise ValueError(
                f"{data.source or '<data>'}: Berry curvature differs among reciprocal "
                f"replicas at key {tuple(unique_keys[group])}; spread={spread:.6g}"
            )
        omega[group] = float(np.mean(values))
        frac[group, :2] = np.mean(canonical[members], axis=0)
        frac[group, 2] = canonical_kz
        if energy is not None:
            e_values = data.energy[members]
            if not np.allclose(e_values, e_values[0], atol=1.0e-7, rtol=1.0e-10):
                raise ValueError("band energy differs among reciprocal replicas")
            energy[group] = float(np.mean(e_values))

    basis = infer_reciprocal_basis(data)
    cart = frac[:, :2] @ basis
    if "b3" in data.metadata:
        cart += frac[:, 2, None] * np.asarray(data.metadata["b3"])[None, :]
    metadata = dict(data.metadata)
    metadata["reciprocal_replica_count"] = int(counts[0])
    metadata["omega_replica_max_spread"] = max_spread
    metadata["raw_row_count"] = len(data.omega)
    return CurvatureData(
        cart=cart,
        frac=frac,
        omega=omega,
        energy=energy,
        metadata=metadata,
        source=data.source,
    )


def attach_band_energy(
    curvature: CurvatureData,
    eigenval: EigenvalData,
    band: int,
    tolerance: float = 2.0e-4,
) -> CurvatureData:
    """Attach E_nk by periodic reciprocal-coordinate matching."""

    if not 1 <= band <= eigenval.nband:
        raise ValueError(f"band {band} is outside EIGENVAL range 1..{eigenval.nband}")
    validate_requested_band(curvature, band)
    validate_eigenval_band_count(curvature, eigenval)
    matched = np.empty(len(curvature.frac), dtype=int)
    max_distance = 0.0
    for i, point in enumerate(curvature.frac):
        distance = np.linalg.norm(_wrap_fractional(eigenval.frac - point), axis=1)
        matched[i] = int(np.argmin(distance))
        max_distance = max(max_distance, float(distance[matched[i]]))
    if max_distance > tolerance:
        raise ValueError(
            f"curvature/EIGENVAL k meshes do not match periodically: max fractional "
            f"distance {max_distance:.3g} > tolerance {tolerance:.3g}"
        )
    curvature.energy = eigenval.energies[matched, band - 1].copy()
    other = np.delete(eigenval.energies[matched], band - 1, axis=1)
    direct_gap = np.min(np.abs(other - curvature.energy[:, None]), axis=1)
    curvature.metadata["band"] = band
    curvature.metadata["energy_source"] = eigenval.source.name
    curvature.metadata["max_k_match_error"] = max_distance
    curvature.metadata["min_direct_band_gap_eV"] = float(np.min(direct_gap))
    curvature.metadata["median_direct_band_gap_eV"] = float(np.median(direct_gap))
    return curvature


def attach_fukui_vertex_energies(
    curvature: CurvatureData,
    eigenval: EigenvalData,
    band: int,
    tolerance: float = 2.0e-5,
) -> CurvatureData:
    """Attach the four EIGENVAL vertex energies of every Fukui plaquette.

    VASPBERRY reports the plaquette center. Its four WAVECAR/EIGENVAL vertices
    are center +/- dk1/2 +/- dk2/2. Fermi weighting therefore uses the
    arithmetic mean of the four *occupations*, rather than an energy or
    occupation evaluated only at the center.
    """

    if curvature.metadata.get("band_mode") != "single-fukui":
        raise ValueError("Fukui vertex attachment requires a single-band Fukui curvature map")
    if not 1 <= band <= eigenval.nband:
        raise ValueError(f"band {band} is outside EIGENVAL range 1..{eigenval.nband}")
    validate_requested_band(curvature, band)
    validate_eigenval_band_count(curvature, eigenval)

    curvature_grid = uniform_grid(curvature.frac, tolerance, "curvature mesh")
    eigen_grid = uniform_grid(eigenval.frac, tolerance, "EIGENVAL mesh")
    if curvature_grid.nx < 2 or curvature_grid.ny < 2:
        raise ValueError(
            f"Fukui transport requires a 2D mesh with Nx,Ny >= 2; "
            f"got {curvature_grid.nx}x{curvature_grid.ny}"
        )
    kz_delta = abs(float(_wrap_fractional(
        np.array([curvature.frac[0, 2] - eigenval.frac[0, 2]])
    )[0]))
    if kz_delta > tolerance:
        raise ValueError(
            f"curvature and EIGENVAL kz planes do not match modulo G "
            f"(fractional difference {kz_delta:.6g})"
        )
    if eigenval.nk != len(curvature.omega):
        raise ValueError(
            f"EIGENVAL is not the matching full mesh: NK={eigenval.nk}, "
            f"curvature unique NK={len(curvature.omega)}"
        )
    if (curvature_grid.nx, curvature_grid.ny) != (eigen_grid.nx, eigen_grid.ny):
        raise ValueError("curvature and EIGENVAL full meshes have different shapes")
    if not np.isclose(np.sum(eigenval.weights), 1.0, atol=2.0e-5, rtol=0.0):
        raise ValueError(f"EIGENVAL k weights sum to {np.sum(eigenval.weights):.8g}, not 1")
    expected_weight = 1.0 / eigenval.nk
    if not np.allclose(eigenval.weights, expected_weight, atol=2.0e-6, rtol=0.0):
        raise ValueError(
            "EIGENVAL contains nonuniform k weights; a full ISYM=-1 uniform mesh is required"
        )

    vertices = np.empty((len(curvature.omega), 4), dtype=int)
    center_xy = _canonical_xy(curvature.frac, tolerance)
    for i, center in enumerate(center_xy):
        x_minus = np.mod(center[0] - 0.5 / curvature_grid.nx, 1.0)
        x_plus = np.mod(center[0] + 0.5 / curvature_grid.nx, 1.0)
        y_minus = np.mod(center[1] - 0.5 / curvature_grid.ny, 1.0)
        y_plus = np.mod(center[1] + 0.5 / curvature_grid.ny, 1.0)
        target_frac = np.array(
            [[x_minus, y_minus, 0.0], [x_plus, y_minus, 0.0],
             [x_plus, y_plus, 0.0], [x_minus, y_plus, 0.0]]
        )
        target_keys = _xy_keys(target_frac, tolerance)
        try:
            vertices[i] = [
                eigen_grid.index_by_xy[(int(key[0]), int(key[1]))] for key in target_keys
            ]
        except KeyError as exc:
            raise ValueError(
                "EIGENVAL mesh does not contain the four vertices implied by the "
                f"Fukui center {center.tolist()}"
            ) from exc

    curvature.vertex_energies = eigenval.energies[vertices, band - 1]
    # Center energy is for visualization only. Transport uses mean[f(E_vertex)].
    curvature.energy = np.mean(curvature.vertex_energies, axis=1)
    other = np.delete(eigenval.energies, band - 1, axis=1)
    direct_gap = np.min(np.abs(other - eigenval.energies[:, band - 1, None]), axis=1)
    curvature.metadata["band"] = band
    curvature.metadata["energy_source"] = eigenval.source.name
    curvature.metadata["energy_quadrature"] = "mean Fermi occupation at four Fukui vertices"
    curvature.metadata["min_direct_band_gap_eV"] = float(np.min(direct_gap))
    curvature.metadata["median_direct_band_gap_eV"] = float(np.median(direct_gap))
    curvature.metadata["eigenval_full_mesh_validated"] = True
    return curvature


def infer_reciprocal_basis(data: CurvatureData) -> np.ndarray:
    """Return a 2x3 matrix whose rows are b1 and b2 in inverse Angstrom."""

    if "b1" in data.metadata and "b2" in data.metadata:
        basis = np.vstack([data.metadata["b1"], data.metadata["b2"]]).astype(float)
        if basis.shape != (2, 3) or not np.all(np.isfinite(basis)):
            raise ValueError("reciprocal basis contains NaN, infinity, or an invalid shape")
        return basis
    frac = data.frac[:, :2]
    design = np.column_stack([frac, np.ones(len(frac))])
    fit, *_ = np.linalg.lstsq(design, data.cart, rcond=None)
    basis = fit[:2]
    residual = np.max(np.linalg.norm(design @ fit - data.cart, axis=1))
    if residual > 2.0e-3:
        raise ValueError(
            "reciprocal vectors are absent from the header and cannot be reliably "
            f"inferred (maximum Cartesian residual {residual:.3g} A^-1)"
        )
    return basis


def validate_cartesian_coordinates(data: CurvatureData, tolerance: float = 2.0e-4) -> None:
    """Check printed Cartesian k coordinates against reciprocal header vectors."""

    validate_finite_curvature(data)
    if "b1" not in data.metadata or "b2" not in data.metadata:
        return
    basis = infer_reciprocal_basis(data)
    expected = data.frac[:, :2] @ basis
    if "b3" in data.metadata:
        b3 = np.asarray(data.metadata["b3"], dtype=float)
        if b3.shape != (3,) or not np.all(np.isfinite(b3)):
            raise ValueError("reciprocal vector b3 contains NaN, infinity, or an invalid shape")
        expected += data.frac[:, 2, None] * b3[None, :]
    residual = np.max(np.linalg.norm(expected - data.cart, axis=1))
    if residual > tolerance:
        raise ValueError(
            f"{data.source or '<data>'}: Cartesian/fractional k coordinates disagree with "
            f"the reciprocal header (max residual {residual:.6g} A^-1)"
        )


def uniform_full_bz_shape(data: CurvatureData, tolerance: float = 2.0e-5) -> tuple[int, int]:
    """Validate a uniform full-BZ product mesh in reciprocal coordinates."""

    grid = uniform_grid(data.frac, tolerance, "curvature mesh")
    if grid.nx < 2 or grid.ny < 2:
        raise ValueError(
            f"transport/map requires a genuinely two-dimensional full-BZ mesh; "
            f"got {grid.nx}x{grid.ny}"
        )
    return grid.nx, grid.ny


def reciprocal_area(data: CurvatureData) -> float:
    basis = infer_reciprocal_basis(data)
    area = float(np.linalg.norm(np.cross(basis[0], basis[1])))
    if not np.isfinite(area) or area <= 0.0:
        raise ValueError("invalid two-dimensional reciprocal-cell area")
    return area


def validate_fukui_geometry(
    data: CurvatureData,
    area_rtol: float = 2.0e-4,
    chern_atol: float = 1.0e-3,
) -> tuple[float, float]:
    """Validate dS and the full-band Fukui Chern sum on the unique mesh."""

    nx, ny = uniform_full_bz_shape(data)
    header_nk = data.metadata.get("nk")
    if header_nk is not None and int(header_nk) != nx * ny:
        raise ValueError(f"header NKPOINT={header_nk} does not match inferred grid {nx}x{ny}")
    header_grid = data.metadata.get("k_grid")
    if header_grid is not None and tuple(int(v) for v in header_grid) != (nx, ny):
        raise ValueError(f"header K-GRID={header_grid} does not match inferred grid {(nx, ny)}")
    dsk = reciprocal_area(data) / (nx * ny)
    header_dsk = data.metadata.get("dk_area")
    if header_dsk is not None and not math.isclose(
        float(header_dsk), dsk, rel_tol=area_rtol, abs_tol=2.0e-6
    ):
        raise ValueError(
            f"{data.source or '<data>'}: header dS={float(header_dsk):.9g} A^-2 "
            f"does not match |b1 x b2|/(Nx Ny)={dsk:.9g} A^-2"
        )
    chern = float(np.sum(data.omega * dsk) / TWO_PI)
    max_flux = float(np.max(np.abs(data.omega * dsk)))
    if max_flux >= 0.8 * math.pi:
        raise ValueError(
            f"{data.source or '<data>'}: maximum Fukui plaquette flux {max_flux:.6g} rad "
            "is too close to the principal-log branch; use a denser k mesh"
        )
    nearest_chern = int(np.rint(chern))
    if abs(chern - nearest_chern) > 5.0e-3:
        raise ValueError(
            f"{data.source or '<data>'}: fully occupied single-Fukui Chern {chern:.9g} "
            f"is not within 0.005 of an integer (nearest {nearest_chern}); check mesh, "
            "band isolation, and branch resolution"
        )
    header_chern = data.metadata.get("chern_reported")
    if header_chern is not None and not math.isclose(
        float(header_chern), chern, rel_tol=2.0e-4, abs_tol=chern_atol
    ):
        raise ValueError(
            f"{data.source or '<data>'}: header Chern={float(header_chern):.9g} "
            f"does not match the unique-grid Fukui sum {chern:.9g}"
        )
    data.metadata["validated_dk_area_A^-2"] = dsk
    data.metadata["validated_unique_grid_chern"] = chern
    data.metadata["nearest_integer_chern"] = nearest_chern
    data.metadata["maximum_abs_plaquette_flux_rad"] = max_flux
    data.metadata["fukui_geometry_validated"] = True
    return dsk, chern


def fermi_dirac(energy: np.ndarray, mu: float, temperature: float) -> np.ndarray:
    energy = np.asarray(energy, dtype=float)
    if not np.all(np.isfinite(energy)):
        raise ValueError("Fermi-Dirac energies contain NaN or infinity")
    if not np.isfinite(mu):
        raise ValueError("chemical potential must be finite")
    if not np.isfinite(temperature):
        raise ValueError("temperature must be finite")
    if temperature < 0.0:
        raise ValueError("temperature must be nonnegative")
    if temperature == 0.0:
        return np.asarray(energy <= mu, dtype=float)
    x = (energy - mu) / (KB_EV_PER_K * temperature)
    out = np.empty_like(x)
    out[x > 40.0] = np.exp(-x[x > 40.0])
    out[x < -40.0] = 1.0 - np.exp(x[x < -40.0])
    middle = np.abs(x) <= 40.0
    out[middle] = 1.0 / (np.exp(x[middle]) + 1.0)
    return out


def periodic_cartesian_distance(
    frac: np.ndarray, center: Sequence[float], basis: np.ndarray
) -> np.ndarray:
    center3 = np.zeros(3, dtype=float)
    center_values = np.asarray(center, dtype=float)
    center3[: len(center_values)] = center_values
    delta = np.asarray(frac, dtype=float) - center3
    return np.asarray([
        np.linalg.norm(_shortest_lattice_image(row[:2], basis) @ basis)
        for row in delta
    ])


def _shortest_lattice_image(delta: Sequence[float], basis: np.ndarray) -> np.ndarray:
    """Solve the two-dimensional closest reciprocal-image problem exactly.

    An initial component-wrapped image supplies an upper bound. The smallest
    singular value of the 2x3 basis then gives a finite coefficient-space ball
    containing every image that can improve that bound. This remains correct
    for oblique, non-reduced bases where the nearest image may require an
    integer shift outside the surrounding 3x3 copies.
    """

    base = np.asarray(delta, dtype=float)
    reciprocal_basis = np.asarray(basis, dtype=float)
    if base.shape != (2,) or reciprocal_basis.shape != (2, 3):
        raise ValueError("reciprocal-image search requires a 2-vector and a 2x3 basis")
    if not np.all(np.isfinite(base)) or not np.all(np.isfinite(reciprocal_basis)):
        raise ValueError("reciprocal-image search inputs must be finite")
    singular_values = np.linalg.svd(reciprocal_basis, compute_uv=False)
    scale = float(singular_values[0])
    smallest = float(singular_values[-1])
    if scale <= 0.0 or smallest <= 1.0e-12 * scale:
        raise ValueError("reciprocal in-plane basis is singular or numerically ill-conditioned")

    initial_shift = -np.rint(base).astype(int)
    best = base + initial_shift
    best_norm = float(np.linalg.norm(best @ reciprocal_basis))
    coefficient_radius = best_norm / smallest + 32.0 * np.finfo(float).eps
    lower = np.ceil(-base - coefficient_radius).astype(int)
    upper = np.floor(-base + coefficient_radius).astype(int)
    counts = upper - lower + 1
    candidate_count = int(counts[0]) * int(counts[1])
    if np.any(counts <= 0):
        return best
    if candidate_count > 1_000_000:
        raise ValueError(
            "reciprocal basis is too ill-conditioned for a safe nearest-image search; "
            "use a reduced in-plane reciprocal basis"
        )

    for shift_x in range(int(lower[0]), int(upper[0]) + 1):
        for shift_y in range(int(lower[1]), int(upper[1]) + 1):
            candidate = base + np.array([shift_x, shift_y])
            norm = float(np.linalg.norm(candidate @ reciprocal_basis))
            if norm < best_norm:
                best = candidate
                best_norm = norm
    return best


def shortest_reciprocal_delta(
    start: Sequence[float], end: Sequence[float], basis: np.ndarray
) -> np.ndarray:
    """Return the 2D reciprocal-image displacement with shortest Cartesian norm."""

    base = np.asarray(end, dtype=float)[:2] - np.asarray(start, dtype=float)[:2]
    return _shortest_lattice_image(base, basis)


def valley_masks(
    data: CurvatureData,
    k_center: Sequence[float] | None,
    kp_center: Sequence[float] | None,
    radius: float | None,
) -> dict[str, np.ndarray]:
    n = len(data.omega)
    if k_center is None and kp_center is None:
        return {"total": np.ones(n, dtype=bool)}
    if k_center is None or kp_center is None or radius is None:
        raise ValueError("K, K' and a Cartesian valley radius must be supplied together")
    if not np.all(np.isfinite(np.asarray(k_center, dtype=float))) or not np.all(
        np.isfinite(np.asarray(kp_center, dtype=float))
    ):
        raise ValueError("K and K' centers must be finite")
    if not np.isfinite(radius) or radius <= 0.0:
        raise ValueError("valley radius must be positive and finite")
    basis = infer_reciprocal_basis(data)
    center_probe = np.asarray([list(k_center) + [0.0]])
    center_separation = periodic_cartesian_distance(center_probe, kp_center, basis)[0]
    if 2.0 * radius >= center_separation:
        raise ValueError(
            f"K and K' disks overlap (2r={2.0 * radius:.6g} A^-1 >= "
            f"separation={center_separation:.6g} A^-1); reduce --valley-radius"
        )
    # Fukui outputs are already located at plaquette centers.
    k_distance = periodic_cartesian_distance(data.frac, k_center, basis)
    kp_distance = periodic_cartesian_distance(data.frac, kp_center, basis)
    mask_k = (k_distance <= radius) & (k_distance <= kp_distance)
    mask_kp = (kp_distance <= radius) & (kp_distance < k_distance)
    if not np.any(mask_k) or not np.any(mask_kp):
        raise ValueError(
            f"valley mask is empty on this mesh (K cells={int(np.sum(mask_k))}, "
            f"K' cells={int(np.sum(mask_kp))}); increase the radius or k-mesh density"
        )
    return {
        "total": np.ones(n, dtype=bool),
        "K": mask_k,
        "Kp": mask_kp,
        "outside": ~(mask_k | mask_kp),
    }


def integrate_sigma(
    datasets: Sequence[CurvatureData],
    mu_values: np.ndarray,
    temperature: float,
    k_center: Sequence[float] | None = None,
    kp_center: Sequence[float] | None = None,
    valley_radius: float | None = None,
    isolation_tolerance: float = 1.0e-5,
    core_chern: float = 0.0,
) -> dict[str, np.ndarray]:
    """Integrate band-resolved curvature on a uniform full 2D BZ.

    Returns both the occupation-weighted Chern integral and the conventional
    electron Hall sign, sigma_xy/(e^2/h) = -C_occ.
    """

    if not datasets:
        raise ValueError("at least one curvature dataset is required")
    mu_values = np.asarray(mu_values, dtype=float)
    if mu_values.ndim != 1 or len(mu_values) == 0 or not np.all(np.isfinite(mu_values)):
        raise ValueError("chemical-potential grid must be a nonempty finite one-dimensional array")
    if not np.isfinite(temperature) or temperature < 0.0:
        raise ValueError("temperature must be finite and nonnegative")
    if not np.isfinite(isolation_tolerance) or isolation_tolerance < 0.0:
        raise ValueError("isolation tolerance must be finite and nonnegative")
    if not np.isfinite(core_chern):
        raise ValueError("core Chern baseline must be finite")
    active_bands = [
        int(data.metadata["band"])
        for data in datasets
        if data.metadata.get("band") is not None
    ]
    if len(active_bands) == len(datasets):
        validate_core_chern(core_chern, active_bands)
    elif abs(core_chern - np.rint(core_chern)) > 5.0e-3:
        raise ValueError("core Chern baseline must be within 0.005 of an integer")
    reference = datasets[0]
    for data in datasets:
        validate_finite_curvature(data)
    shape = uniform_full_bz_shape(reference)
    reference_basis = infer_reciprocal_basis(reference)
    validate_cartesian_coordinates(reference)
    reference_dsk, _ = validate_fukui_geometry(reference)
    masks = valley_masks(reference, k_center, kp_center, valley_radius)
    for data in datasets:
        validate_cartesian_coordinates(data)
        if not np.allclose(
            infer_reciprocal_basis(data), reference_basis, atol=2.0e-6, rtol=2.0e-6
        ):
            raise ValueError("band datasets have inconsistent reciprocal basis vectors")
        mode = data.metadata.get("band_mode")
        if mode == "manifold":
            raise ValueError(
                f"{data.source or '<data>'}: an entire band-manifold curvature cannot be "
                "weighted by one band energy; provide band-resolved curvature instead"
            )
        if mode == "single-kubo":
            raise ValueError(
                f"{data.source or '<data>'}: legacy Kubo output is supported for plotting "
                "only and is not accepted as validated transport input"
            )
        if mode != "single-fukui":
            raise ValueError(f"{data.source or '<data>'}: unsupported curvature mode {mode!r}")
        if data.vertex_energies is None:
            raise ValueError(
                f"{data.source or '<data>'}: four Fukui-plaquette vertex energies are "
                "required; attach a matching full-mesh EIGENVAL"
            )
        min_gap = data.metadata.get("min_direct_band_gap_eV")
        if min_gap is None:
            raise ValueError(
                f"{data.source or '<data>'}: missing min_direct_band_gap_eV provenance; "
                "validated single-band transport requires an isolation check"
            )
        if not np.isfinite(float(min_gap)) or float(min_gap) < 0.0:
            raise ValueError("minimum direct band-gap provenance must be finite and nonnegative")
        if float(min_gap) <= isolation_tolerance:
            raise ValueError(
                f"{data.source or '<data>'}: band {data.metadata.get('band', '?')} is not "
                f"numerically isolated (minimum direct gap {float(min_gap):.6g} eV <= "
                f"{isolation_tolerance:.6g} eV); use a gauge-covariant multiband treatment"
            )
        if uniform_full_bz_shape(data) != shape:
            raise ValueError("all band datasets must use the same full-BZ grid shape")
        if len(data.omega) != len(reference.omega):
            raise ValueError("all band datasets must have the same number of k points")
        mismatch = np.max(np.linalg.norm(_wrap_fractional(data.frac - reference.frac), axis=1))
        if mismatch > 2.0e-4:
            raise ValueError("band datasets are not in the same periodic k-point order")
        data_dsk, _ = validate_fukui_geometry(data)
        if not math.isclose(data_dsk, reference_dsk, rel_tol=2.0e-4, abs_tol=1.0e-10):
            raise ValueError("band datasets have inconsistent Fukui plaquette areas")

    result: dict[str, np.ndarray] = {"mu_eV": np.asarray(mu_values, dtype=float)}
    for region, mask in masks.items():
        c_values = np.zeros(len(mu_values), dtype=float)
        if region == "total":
            c_values += core_chern
        for index, mu in enumerate(mu_values):
            for data in datasets:
                vertex_occupation = fermi_dirac(data.vertex_energies, float(mu), temperature)
                plaquette_occupation = np.mean(vertex_occupation, axis=1)
                plaquette_flux = data.omega * reference_dsk
                c_values[index] += (
                    np.sum(plaquette_flux[mask] * plaquette_occupation[mask]) / TWO_PI
                )
        result[f"chern_weight_{region}"] = c_values
        result[f"sigma_xy_{region}_e2_over_h"] = -c_values
    if "sigma_xy_K_e2_over_h" in result:
        result["sigma_xy_valley_contrast_e2_over_h"] = (
            result["sigma_xy_K_e2_over_h"] - result["sigma_xy_Kp_e2_over_h"]
        )
        result["sigma_xy_valley_sum_e2_over_h"] = (
            result["sigma_xy_K_e2_over_h"] + result["sigma_xy_Kp_e2_over_h"]
        )
    result["temperature_K"] = np.full(len(mu_values), temperature)
    result["core_chern"] = np.full(len(mu_values), core_chern)
    return result


def write_csv(path: str | Path, columns: dict[str, np.ndarray]) -> None:
    output = Path(path)
    names = list(columns)
    length = len(columns[names[0]])
    if any(len(columns[name]) != length for name in names):
        raise ValueError("all CSV columns must have the same length")
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(names)
        for row in zip(*(columns[name] for name in names), strict=True):
            writer.writerow([f"{float(value):.12g}" for value in row])


def remove_stale_sigma_outputs(
    csv_path: str | Path,
    plot_path: str | Path,
    summary_path: str | Path | None,
) -> None:
    """Remove only the requested sigma artifacts before a new validation run."""

    for target in (csv_path, plot_path, summary_path):
        if target is not None:
            Path(target).unlink(missing_ok=True)


def validate_sigma_output_targets(
    input_paths: Sequence[str | Path],
    eigenval_path: str | Path,
    csv_path: str | Path,
    plot_path: str | Path,
    summary_path: str | Path | None,
) -> None:
    """Reject aliased sigma targets before any stale artifact is removed."""

    outputs = [("--csv", csv_path), ("--plot", plot_path)]
    if summary_path is not None:
        outputs.append(("--summary", summary_path))
    protected = [("--input", path) for path in input_paths]
    protected.append(("--eigenval", eigenval_path))

    def resolved(label: str, path: str | Path) -> Path:
        try:
            return Path(path).resolve(strict=False)
        except (OSError, RuntimeError) as exc:
            raise ValueError(f"cannot safely resolve {label} path {path!s}: {exc}") from exc

    resolved_outputs = [(label, path, resolved(label, path)) for label, path in outputs]
    for index, (label, path, canonical) in enumerate(resolved_outputs):
        for other_label, other_path, other_canonical in resolved_outputs[:index]:
            if canonical == other_canonical:
                raise ValueError(
                    f"sigma output targets must be distinct: {other_label}={other_path} "
                    f"and {label}={path} resolve to {canonical}"
                )

    resolved_protected = [
        (label, path, resolved(label, path)) for label, path in protected
    ]
    for label, path, canonical in resolved_outputs:
        for protected_label, protected_path, protected_canonical in resolved_protected:
            if canonical == protected_canonical:
                raise ValueError(
                    f"refusing to overwrite transport input: {label}={path} and "
                    f"{protected_label}={protected_path} resolve to {canonical}"
                )


def validate_plot_output_target(
    input_path: str | Path,
    output_path: str | Path,
    eigenval_path: str | Path | None = None,
) -> None:
    """Reject a plot target that aliases any input before reading or writing."""

    try:
        canonical_output = Path(output_path).resolve(strict=False)
        protected = [("--input", Path(input_path).resolve(strict=False))]
        if eigenval_path is not None:
            protected.append(
                ("--eigenval", Path(eigenval_path).resolve(strict=False))
            )
    except (OSError, RuntimeError) as exc:
        raise ValueError(f"cannot safely resolve plot input/output paths: {exc}") from exc

    for label, canonical_input in protected:
        if canonical_output == canonical_input:
            raise ValueError(
                f"refusing to overwrite plot input: --output={output_path} and "
                f"{label} resolve to {canonical_output}"
            )


def plot_map(
    data: CurvatureData,
    output: str | Path,
    k_center: Sequence[float] | None = None,
    kp_center: Sequence[float] | None = None,
    title: str | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(7.1, 5.8), constrained_layout=True)
    vmax = float(np.nanpercentile(np.abs(data.omega), 99.0))
    if vmax == 0.0:
        vmax = 1.0
    if len(data.omega) >= 4:
        artist = ax.tricontourf(
            data.cart[:, 0], data.cart[:, 1], data.omega, levels=101,
            cmap="RdBu_r", vmin=-vmax, vmax=vmax, extend="both",
        )
    else:
        artist = ax.scatter(data.cart[:, 0], data.cart[:, 1], c=data.omega, cmap="RdBu_r")
    colorbar = fig.colorbar(artist, ax=ax)
    colorbar.set_label(r"Berry curvature $\Omega_z$ ($\AA^2$)")
    basis = infer_reciprocal_basis(data)
    for center, label, marker in ((k_center, "K", "o"), (kp_center, "K'", "s")):
        if center is not None:
            frac = np.asarray(center[:2], dtype=float)
            cart = frac @ basis
            ax.scatter(cart[0], cart[1], s=65, marker=marker, facecolors="none", edgecolors="black")
            ax.annotate(label, cart[:2], xytext=(5, 5), textcoords="offset points")
    ax.set_xlabel(r"$k_x$ ($\AA^{-1}$)")
    ax.set_ylabel(r"$k_y$ ($\AA^{-1}$)")
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(title or "k-resolved Berry curvature")
    fig.savefig(output, dpi=220)
    plt.close(fig)


def _path_distance(cart: np.ndarray) -> np.ndarray:
    step = np.linalg.norm(np.diff(cart, axis=0), axis=1)
    return np.r_[0.0, np.cumsum(step)]


def select_fundamental_copy(data: CurvatureData) -> CurvatureData:
    """Select one translated copy when legacy ``-kp`` output repeats data.

    VASPBERRY can print a 3x3 set of reciprocal translations. Those copies must
    not be connected into one path or counted nine times in a BZ integral. The
    header's NKPOINT count identifies the length of one copy; the copy closest
    to the Cartesian origin is selected only when the row count is an integer
    multiple of NKPOINT.
    """

    nk = data.metadata.get("nk")
    if not isinstance(nk, int) or nk <= 0 or len(data.omega) <= nk or len(data.omega) % nk:
        return data
    blocks = len(data.omega) // nk
    scores = [
        float(np.mean(np.linalg.norm(data.cart[i * nk:(i + 1) * nk, :2], axis=1)))
        for i in range(blocks)
    ]
    block = int(np.argmin(scores))
    section = slice(block * nk, (block + 1) * nk)
    metadata = dict(data.metadata)
    metadata["selected_periodic_block"] = block
    return CurvatureData(
        cart=data.cart[section],
        frac=data.frac[section],
        omega=data.omega[section],
        energy=None if data.energy is None else data.energy[section],
        metadata=metadata,
        source=data.source,
    )


def select_fundamental_path(data: CurvatureData) -> CurvatureData:
    """Backward-compatible name for selecting one repeated path copy."""

    return select_fundamental_copy(data)


def prepare_full_bz_map(data: CurvatureData) -> CurvatureData:
    """Collapse reciprocal replicas and require a genuine 2D full-BZ mesh."""

    try:
        validate_cartesian_coordinates(data)
        unique = deduplicate_reciprocal_replicas(data)
        uniform_full_bz_shape(unique)
    except ValueError as exc:
        raise ValueError(f"map/cut input must be a 2D full-BZ mesh: {exc}") from exc
    return unique


def prepare_k_path(data: CurvatureData) -> CurvatureData:
    """Select one periodic path copy and reject full 2D maps."""

    validate_cartesian_coordinates(data)
    header_grid = data.metadata.get("k_grid")
    header_nk = data.metadata.get("nk")
    header_declares_map = (
        header_grid is not None
        and header_nk is not None
        and int(header_grid[0]) >= 2
        and int(header_grid[1]) >= 2
        and int(header_grid[0]) * int(header_grid[1]) == int(header_nk)
    )
    if header_declares_map:
        # This also preserves replica-consistency failures instead of turning a
        # corrupted 2D map into a plausible-looking line.
        prepare_full_bz_map(data)
        raise ValueError("line input is a 2D full-BZ map; use the map or cut subcommand")
    path = select_fundamental_path(data)
    if len(path.omega) < 2:
        raise ValueError("line input must contain at least two ordered path points")
    try:
        uniform_full_bz_shape(path)
    except ValueError:
        pass
    else:
        raise ValueError("line input is a 2D full-BZ map; use the map or cut subcommand")
    if not np.isfinite(_path_distance(path.cart)[-1]) or _path_distance(path.cart)[-1] <= 0.0:
        raise ValueError("line input does not define a nonzero finite k path")
    return path


def plot_path(data: CurvatureData, output: str | Path, title: str | None = None) -> None:
    data = select_fundamental_path(data)
    distance = _path_distance(data.cart)
    rows = 2 if data.energy is not None else 1
    fig, axes = plt.subplots(rows, 1, figsize=(7.4, 3.2 + 2.3 * (rows - 1)), sharex=True,
                             constrained_layout=True)
    axes_array = np.atleast_1d(axes)
    axes_array[0].plot(distance, data.omega, color="#a3132a", lw=1.6)
    axes_array[0].axhline(0.0, color="0.55", lw=0.7)
    axes_array[0].set_ylabel(r"$\Omega_z$ ($\AA^2$)")
    if data.energy is not None:
        axes_array[1].plot(distance, data.energy, color="#174f78", lw=1.6)
        axes_array[1].set_ylabel("Energy (eV)")
    axes_array[-1].set_xlabel(r"Cumulative path length ($\AA^{-1}$)")
    axes_array[0].set_title(title or "Berry curvature along k path")
    fig.savefig(output, dpi=220)
    plt.close(fig)


def periodic_bilinear_interpolate(
    data: CurvatureData,
    values: np.ndarray,
    sample_frac: np.ndarray,
    tolerance: float = 2.0e-5,
) -> np.ndarray:
    """Interpolate a scalar field on a validated uniform 2D reciprocal torus."""

    grid = uniform_grid(data.frac, tolerance, "cut interpolation mesh")
    if grid.nx < 2 or grid.ny < 2:
        raise ValueError(
            f"cut interpolation requires Nx,Ny >= 2; got {grid.nx}x{grid.ny}"
        )
    field = np.asarray(values, dtype=float)
    samples = np.asarray(sample_frac, dtype=float)
    if field.shape != (len(data.omega),):
        raise ValueError("cut interpolation values must have shape (N,)")
    if samples.ndim != 2 or samples.shape[1] != 2:
        raise ValueError("cut interpolation samples must have shape (M, 2)")
    if not np.all(np.isfinite(field)) or not np.all(np.isfinite(samples)):
        raise ValueError("cut interpolation inputs contain NaN or infinity")

    x_keys = np.sort(grid.x_keys)
    y_keys = np.sort(grid.y_keys)
    x_rank = {int(key): index for index, key in enumerate(x_keys)}
    y_rank = {int(key): index for index, key in enumerate(y_keys)}
    value_grid = np.empty((grid.nx, grid.ny), dtype=float)
    for row, key in enumerate(grid.keys):
        value_grid[x_rank[int(key[0])], y_rank[int(key[1])]] = field[row]

    canonical = _canonical_xy(data.frac, tolerance)
    x0 = float(np.mean(canonical[grid.keys[:, 0] == x_keys[0], 0]))
    y0 = float(np.mean(canonical[grid.keys[:, 1] == y_keys[0], 1]))
    scaled_x = np.mod(samples[:, 0] - x0, 1.0) * grid.nx
    scaled_y = np.mod(samples[:, 1] - y0, 1.0) * grid.ny
    ix0_unwrapped = np.floor(scaled_x).astype(int)
    iy0_unwrapped = np.floor(scaled_y).astype(int)
    tx = scaled_x - ix0_unwrapped
    ty = scaled_y - iy0_unwrapped
    ix0 = np.mod(ix0_unwrapped, grid.nx)
    iy0 = np.mod(iy0_unwrapped, grid.ny)
    ix1 = (ix0 + 1) % grid.nx
    iy1 = (iy0 + 1) % grid.ny
    interpolated = (
        (1.0 - tx) * (1.0 - ty) * value_grid[ix0, iy0]
        + tx * (1.0 - ty) * value_grid[ix1, iy0]
        + (1.0 - tx) * ty * value_grid[ix0, iy1]
        + tx * ty * value_grid[ix1, iy1]
    )
    if not np.all(np.isfinite(interpolated)):
        raise ValueError("periodic cut interpolation produced NaN or infinity")
    return interpolated


def plot_valley_cut(
    data: CurvatureData,
    k_center: Sequence[float],
    kp_center: Sequence[float],
    output: str | Path,
    samples: int = 401,
) -> None:
    """Interpolate Omega (and optionally E) on a straight K-to-K' cut."""

    basis = infer_reciprocal_basis(data)
    k0 = np.asarray(k_center[:2], dtype=float)
    delta = shortest_reciprocal_delta(k0, kp_center, basis)
    t = np.linspace(0.0, 1.0, samples)
    cut_frac = k0[None, :] + t[:, None] * delta[None, :]

    omega_cut = periodic_bilinear_interpolate(data, data.omega, cut_frac)
    cart_cut = cut_frac @ basis
    distance = _path_distance(cart_cut)

    rows = 2 if data.energy is not None else 1
    fig, axes = plt.subplots(rows, 1, figsize=(7.4, 3.2 + 2.3 * (rows - 1)), sharex=True,
                             constrained_layout=True)
    axes_array = np.atleast_1d(axes)
    axes_array[0].plot(distance, omega_cut, color="#a3132a", lw=1.7)
    axes_array[0].axhline(0.0, color="0.55", lw=0.7)
    axes_array[0].set_ylabel(r"$\Omega_z$ ($\AA^2$)")
    if data.energy is not None:
        energy_cut = periodic_bilinear_interpolate(data, data.energy, cut_frac)
        axes_array[1].plot(distance, energy_cut, color="#174f78", lw=1.7)
        axes_array[1].set_ylabel("Energy (eV)")
    axes_array[-1].set_xlabel(r"K $\rightarrow$ K' path length ($\AA^{-1}$)")
    axes_array[0].set_title("Interpolated K–K' valley cut")
    for ax in axes_array:
        ax.axvline(distance[0], color="0.25", lw=0.8, ls="--")
        ax.axvline(distance[-1], color="0.25", lw=0.8, ls="--")
    axes_array[-1].set_xticks([distance[0], distance[-1]], ["K", "K'"])
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_sigma(columns: dict[str, np.ndarray], output: str | Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    mu = columns["mu_eV"]
    preferred = [
        ("sigma_xy_total_e2_over_h", r"total $\sigma_{xy}$", "black"),
        ("sigma_xy_K_e2_over_h", r"K contribution", "#b31b34"),
        ("sigma_xy_Kp_e2_over_h", r"K' contribution", "#1668a5"),
        ("sigma_xy_valley_contrast_e2_over_h", r"K$-$K' contrast", "#6b3fa0"),
    ]
    for key, label, color in preferred:
        if key in columns:
            ax.plot(mu, columns[key], label=label, color=color, lw=1.8)
    ax.axhline(0.0, color="0.6", lw=0.7)
    ax.set_xlabel(r"Chemical potential $\mu$ (eV)")
    ax.set_ylabel(r"$\sigma_{xy}$ ($e^2/h$)")
    ax.set_title("Intrinsic Hall conductivity from Berry curvature")
    ax.legend(frameon=False)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def _load_with_optional_energy(args: argparse.Namespace, sampling: str) -> CurvatureData:
    data = read_vaspberry(args.input, getattr(args, "energy_column", None))
    eigenval_path = getattr(args, "eigenval", None)
    validate_spin_channels(
        [data], args.spin, require_collinear_label=eigenval_path is not None
    )
    if sampling == "map":
        data = prepare_full_bz_map(data)
    elif sampling == "path":
        data = prepare_k_path(data)
    else:
        raise ValueError(f"unknown sampling mode {sampling!r}")
    band = getattr(args, "band", None)
    if eigenval_path:
        if band is None:
            inferred = data.metadata.get("band")
            if inferred is None:
                raise ValueError("--band is required because the curvature header has no single band index")
            band = int(inferred)
        eigenval = read_eigenval(eigenval_path, args.spin)
        if sampling == "map":
            attach_fukui_vertex_energies(data, eigenval, band)
        else:
            attach_band_energy(data, eigenval, band)
    return data


def _pair_values(values: Iterable[Sequence[float]] | None) -> Sequence[float] | None:
    if not values:
        return None
    value = list(values)
    if len(value) != 1:
        raise ValueError("a valley center must be specified once")
    return value[0]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--input", required=True, help="VASPBERRY curvature data file")
    common.add_argument("--eigenval", help="matching EIGENVAL used to attach band energies")
    common.add_argument("--band", type=int, help="one-based band index in EIGENVAL")
    common.add_argument("--spin", type=int, default=1, choices=(1, 2), help="EIGENVAL spin column")
    common.add_argument("--energy-column", type=int, help="zero-based energy column already in --input")

    map_parser = subparsers.add_parser("map", parents=[common], help="plot a k-resolved heat map")
    map_parser.add_argument("--output", required=True)
    map_parser.add_argument("--k-center", nargs=2, type=float, metavar=("K1", "K2"))
    map_parser.add_argument("--kp-center", nargs=2, type=float, metavar=("KP1", "KP2"))
    map_parser.add_argument("--title")

    line_parser = subparsers.add_parser("line", parents=[common], help="plot an ordered k path")
    line_parser.add_argument("--output", required=True)
    line_parser.add_argument("--title")

    cut_parser = subparsers.add_parser(
        "cut", parents=[common], help="interpolate a K-to-K' cut from a 2D full-BZ map"
    )
    cut_parser.add_argument("--output", required=True)
    cut_parser.add_argument("--k-center", nargs=2, type=float, required=True,
                            metavar=("K1", "K2"))
    cut_parser.add_argument("--kp-center", nargs=2, type=float, required=True,
                            metavar=("KP1", "KP2"))
    cut_parser.add_argument("--samples", type=int, default=401)

    sigma_parser = subparsers.add_parser("sigma", help="integrate sigma_xy(mu) on a full BZ mesh")
    sigma_parser.add_argument("--input", action="append", required=True,
                              help="band-resolved curvature file; repeat for multiple bands")
    sigma_parser.add_argument("--bands", nargs="+", type=int,
                              help="one band index per input (otherwise infer from headers)")
    sigma_parser.add_argument("--eigenval", required=True)
    sigma_parser.add_argument("--spin", type=int, default=1, choices=(1, 2))
    sigma_parser.add_argument("--mu-min", type=float, required=True)
    sigma_parser.add_argument("--mu-max", type=float, required=True)
    sigma_parser.add_argument("--mu-points", type=int, default=401)
    sigma_parser.add_argument("--temperature", type=float, default=0.0)
    sigma_parser.add_argument(
        "--isolation-tolerance", type=float, default=1.0e-5,
        help="reject a nominal single band when its direct gap to another band is below this value (eV)",
    )
    sigma_parser.add_argument(
        "--core-chern", type=float,
        help=(
            "integer Chern baseline of fully occupied bands below the active "
            "window; must be supplied explicitly when the first active band is >1"
        ),
    )
    sigma_parser.add_argument(
        "--occupation-tolerance", type=float, default=1.0e-8,
        help=(
            "maximum allowed Fermi-Dirac occupation tail in omitted bands "
            "(default: 1e-8)"
        ),
    )
    sigma_parser.add_argument("--k-center", nargs=2, type=float, metavar=("K1", "K2"))
    sigma_parser.add_argument("--kp-center", nargs=2, type=float, metavar=("KP1", "KP2"))
    sigma_parser.add_argument("--valley-radius", type=float, help="K/K' disk radius in inverse Angstrom")
    sigma_parser.add_argument("--csv", required=True)
    sigma_parser.add_argument("--plot", required=True)
    sigma_parser.add_argument("--summary", help="optional JSON validation summary")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "map":
            validate_plot_output_target(args.input, args.output, args.eigenval)
            for label, center in (("--k-center", args.k_center), ("--kp-center", args.kp_center)):
                if center is not None and not np.all(np.isfinite(center)):
                    raise ValueError(f"{label} must contain finite values")
            data = _load_with_optional_energy(args, "map")
            plot_map(data, args.output, args.k_center, args.kp_center, args.title)
        elif args.command == "line":
            validate_plot_output_target(args.input, args.output, args.eigenval)
            data = _load_with_optional_energy(args, "path")
            plot_path(data, args.output, args.title)
        elif args.command == "cut":
            validate_plot_output_target(args.input, args.output, args.eigenval)
            if args.samples < 2:
                raise ValueError("--samples must be at least 2")
            if not np.all(np.isfinite(args.k_center)) or not np.all(np.isfinite(args.kp_center)):
                raise ValueError("cut valley centers must contain finite values")
            data = _load_with_optional_energy(args, "map")
            plot_valley_cut(data, args.k_center, args.kp_center, args.output, args.samples)
        elif args.command == "sigma":
            # A failed validation must not leave an older result at the exact
            # paths requested for this run. Never delete a parent or glob.
            validate_sigma_output_targets(
                args.input, args.eigenval, args.csv, args.plot, args.summary
            )
            remove_stale_sigma_outputs(args.csv, args.plot, args.summary)
            if not np.isfinite(args.isolation_tolerance) or args.isolation_tolerance < 0.0:
                raise ValueError("--isolation-tolerance must be finite and nonnegative")
            if args.mu_points < 1:
                raise ValueError("--mu-points must be at least 1")
            if not np.isfinite(args.mu_min) or not np.isfinite(args.mu_max):
                raise ValueError("chemical-potential bounds must be finite")
            if args.mu_max < args.mu_min:
                raise ValueError("--mu-max must be greater than or equal to --mu-min")
            if (
                not np.isfinite(args.occupation_tolerance)
                or not 0.0 <= args.occupation_tolerance < 0.5
            ):
                raise ValueError(
                    "--occupation-tolerance must be finite and in the range [0, 0.5)"
                )
            eigenval = read_eigenval(args.eigenval, args.spin)
            raw_datasets = [read_vaspberry(path) for path in args.input]
            for data in raw_datasets:
                validate_cartesian_coordinates(data)
            datasets = [deduplicate_reciprocal_replicas(data) for data in raw_datasets]
            validate_spin_channels(datasets, args.spin)
            if args.bands is not None and len(args.bands) != len(datasets):
                raise ValueError("--bands must contain exactly one index per --input")
            bands: list[int] = []
            for index, data in enumerate(datasets):
                band = args.bands[index] if args.bands is not None else data.metadata.get("band")
                if band is None:
                    raise ValueError(f"cannot infer band for {data.source}; provide --bands")
                bands.append(int(band))
            validate_unique_active_bands(bands)
            if args.core_chern is None:
                if min(bands) > 1:
                    raise ValueError(
                        "--core-chern must be supplied explicitly when active bands "
                        "start above band 1 (use --core-chern 0.0 when zero has been "
                        "established)"
                    )
                core_chern = 0.0
            else:
                core_chern = float(args.core_chern)
            mu = np.linspace(args.mu_min, args.mu_max, args.mu_points)
            energy_window_validation = validate_legacy_sigma_energy_window(
                eigenval,
                bands,
                mu,
                args.temperature,
                args.occupation_tolerance,
                core_chern,
            )
            for data, band in zip(datasets, bands, strict=True):
                attach_fukui_vertex_energies(data, eigenval, int(band))
            columns = integrate_sigma(
                datasets, mu, args.temperature, args.k_center, args.kp_center,
                args.valley_radius, args.isolation_tolerance, core_chern,
            )
            write_csv(args.csv, columns)
            plot_sigma(columns, args.plot)
            if args.summary:
                masks = valley_masks(
                    datasets[0], args.k_center, args.kp_center, args.valley_radius
                )
                summary = {
                    "inputs": [data.source.name for data in datasets],
                    "eigenval": Path(args.eigenval).name,
                    "spin": args.spin,
                    "bands": bands,
                    "grid": uniform_full_bz_shape(datasets[0]),
                    "reciprocal_area_A^-2": reciprocal_area(datasets[0]),
                    "plaquette_area_A^-2": datasets[0].metadata.get("validated_dk_area_A^-2"),
                    "unique_grid_chern": [
                        data.metadata.get("validated_unique_grid_chern") for data in datasets
                    ],
                    "replica_count": [
                        data.metadata.get("reciprocal_replica_count") for data in datasets
                    ],
                    "temperature_K": args.temperature,
                    "mu_range_eV": [args.mu_min, args.mu_max],
                    "mu_grid": {
                        "points": args.mu_points,
                        "values_eV": mu.tolist(),
                    },
                    "valleys": {
                        "K_center_fractional": args.k_center,
                        "Kp_center_fractional": args.kp_center,
                        "radius_A^-1": args.valley_radius,
                        "mask_counts": {name: int(np.sum(mask)) for name, mask in masks.items()},
                    },
                    "isolation_tolerance_eV": args.isolation_tolerance,
                    "minimum_direct_band_gaps_eV": [
                        data.metadata.get("min_direct_band_gap_eV") for data in datasets
                    ],
                    "core_chern": core_chern,
                    "core_chern_provenance": "explicit CLI input" if args.core_chern is not None else "no omitted lower bands",
                    "occupation_tolerance": args.occupation_tolerance,
                    "energy_window_validation": energy_window_validation,
                    "sign_convention": "sigma_xy/(e^2/h) = - occupied_Chern_weight",
                    "zero_temperature_occupation": "E <= mu is occupied",
                    "energy_quadrature": "mean Fermi occupation at four Fukui plaquette vertices",
                    "model_scope": "intrinsic 2D Berry-curvature term; rigid-band occupations",
                }
                Path(args.summary).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        return 0
    except (OSError, ValueError) as exc:
        parser.exit(2, f"error: {exc}\n")


if __name__ == "__main__":
    sys.exit(main())
