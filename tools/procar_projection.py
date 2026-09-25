"""Strict VASP 5.4.4 LORBIT=11 noncollinear PROCAR projections.

Raw component order is (q, m1, m2, m3), where m is projected Pauli weight
in the SAXIS frame, not spin angular momentum. No weights are normalized or
clipped. Group joint-spin weights are (q_group +/- axis.dot(m_group))/2;
they are not group charge times the whole-state spin polarization.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
from typing import Mapping, Sequence

import numpy as np


_FLOAT = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
# VASP writes coordinates in adjacent fixed-width fields. A minus sign may
# therefore follow the previous coordinate with no whitespace (including -0).
_COORD_SEPARATOR = r"(?:\s+|(?=[+-]))"
_KPOINT_HEADER = re.compile(
    rf"k-point\s+(\d+)\s*:\s*({_FLOAT}){_COORD_SEPARATOR}"
    rf"({_FLOAT}){_COORD_SEPARATOR}({_FLOAT})\s+weight\s*=\s*({_FLOAT})")


class ProcarFormatError(ValueError):
    """Missing, ambiguous, truncated, or inconsistent projection data."""


@dataclass(frozen=True)
class ProcarProjection:
    source: Path
    kpoint_ids: np.ndarray
    band_ids: np.ndarray
    ion_ids: np.ndarray
    kpoints: np.ndarray                 # (k, 3), fractional
    weights: np.ndarray                 # (k,), as printed; not renormalized
    energies_eV: np.ndarray              # (k, band)
    occupations: np.ndarray              # (k, band), not folded into projections
    orbital_names: tuple[str, ...]
    ion_orbital: np.ndarray              # (k, band, ion, orbital, component)
    ion_totals: np.ndarray               # (k, band, ion, component), printed tot
    orbital_totals: np.ndarray           # (k, band, orbital, component), tot row
    state_totals: np.ndarray             # (k, band, component), final tot column
    spin_to_cartesian: np.ndarray        # m_cart = R @ m_SAXIS
    rounding_diagnostics: dict[str, float]

    @property
    def nkpoints(self) -> int:
        return len(self.kpoint_ids)

    @property
    def nbands(self) -> int:
        return len(self.band_ids)

    @property
    def nions(self) -> int:
        return len(self.ion_ids)

    @property
    def kpoint_weights(self) -> np.ndarray:
        return self.weights

    def pauli_cartesian(self, raw_components: np.ndarray) -> np.ndarray:
        """Rotate a raw (...,4) q/m array; retain raw Pauli normalization."""
        raw = np.asarray(raw_components, dtype=float)
        if raw.ndim < 1 or raw.shape[-1] != 4 or not np.isfinite(raw).all():
            raise ValueError("finite raw q,m1,m2,m3 components required")
        return np.einsum("ij,...j->...i", self.spin_to_cartesian, raw[..., 1:])


@dataclass(frozen=True)
class GroupProjection:
    name: str
    ion_ids: tuple[int, ...]             # user-supplied 1-based atom IDs
    q: np.ndarray                       # (k, band), sum of printed ion totals
    m_cartesian: np.ndarray             # (k, band, 3), Pauli weights
    m_axis: np.ndarray                   # (k, band), axis.dot(m_cartesian)
    p_plus: np.ndarray                   # (q + m_axis)/2, not clipped
    p_minus: np.ndarray                  # (q - m_axis)/2, not clipped
    axis_cartesian: np.ndarray


def _rotation(matrix: np.ndarray) -> np.ndarray:
    if matrix is None:
        raise ValueError("explicit SAXIS-to-Cartesian matrix required; supply identity explicitly")
    rotation = np.asarray(matrix, dtype=float)
    if rotation.shape != (3, 3) or not np.isfinite(rotation).all():
        raise ValueError("spin_to_cartesian must be a finite 3x3 matrix")
    # OUTCAR's printed Euler matrix has finite decimal precision. Do not repair
    # it silently, but admit the rounding error of its seven-place output.
    if (not np.allclose(rotation @ rotation.T, np.eye(3), atol=2e-6, rtol=0)
            or not np.isclose(np.linalg.det(rotation), 1, atol=2e-6, rtol=0)):
        raise ValueError("spin_to_cartesian must be a proper orthonormal rotation")
    return rotation.copy()


def _numbers(cells: Sequence[str], context: str) -> np.ndarray:
    try:
        values = np.array([float(x.replace("D", "E").replace("d", "e")) for x in cells])
    except ValueError as exc:
        raise ProcarFormatError(f"{context}: invalid numeric value") from exc
    if not np.isfinite(values).all():
        raise ProcarFormatError(f"{context}: nonfinite numeric value")
    return values


def parse_procar(path: str | Path, *, noncollinear: bool,
                 spin_to_cartesian: np.ndarray) -> ProcarProjection:
    """Read the four-block 5.4.4 ncl format; unsupported modes fail explicitly.

    Both ``noncollinear=True`` and an explicit SAXIS-to-Cartesian rotation are
    required. The format has one orbital header per band followed by q/m1/m2/m3
    blocks. Collinear/spin-component sections and LORBIT>=12 phase blocks are
    intentionally unsupported. Arrays are indexed by canonical 1-based IDs,
    converted to zero-based NumPy positions, even if file IDs are shuffled.
    """
    if noncollinear is not True:
        raise ValueError("this parser requires explicit noncollinear=True and four projection blocks")
    rotation = _rotation(spin_to_cartesian)
    path = Path(path)
    lines = [(number, line.strip()) for number, line in enumerate(path.read_text().splitlines(), 1)
             if line.strip()]
    if len(lines) < 2 or not lines[0][1].startswith("PROCAR"):
        raise ProcarFormatError("missing PROCAR title or dimension header")
    dimensions = re.fullmatch(
        r"# of k-points:\s*(\d+)\s+# of bands:\s*(\d+)\s+# of ions:\s*(\d+)",
        lines[1][1])
    if dimensions is None:
        raise ProcarFormatError("invalid PROCAR k-point/band/ion dimension header")
    nk, nb, ni = map(int, dimensions.groups())
    if min(nk, nb, ni) <= 0:
        raise ProcarFormatError("PROCAR dimensions must be positive")
    required_lines = 2 + nk * (1 + nb * (2 + 4 * (ni + 1)))
    if len(lines) < required_lines:
        raise ProcarFormatError("truncated PROCAR: expected four complete q/m blocks per state")

    def orbital_header(text: str) -> tuple[str, ...]:
        cells = text.split()
        if len(cells) < 3 or cells[0] != "ion" or cells[-1] != "tot":
            raise ProcarFormatError("expected one 'ion ... tot' orbital header per band")
        names = tuple(cells[1:-1])
        if len(set(names)) != len(names) or "tot" in names or "ion" in names:
            raise ProcarFormatError("duplicate or reserved orbital names")
        return names

    orbitals = orbital_header(lines[4][1])
    no = len(orbitals)
    kpoints = np.full((nk, 3), np.nan)
    weights = np.full(nk, np.nan)
    energies = np.full((nk, nb), np.nan)
    occupations = np.full((nk, nb), np.nan)
    ion_orbital = np.full((nk, nb, ni, no, 4), np.nan)
    ion_totals = np.full((nk, nb, ni, 4), np.nan)
    orbital_totals = np.full((nk, nb, no, 4), np.nan)
    state_totals = np.full((nk, nb, 4), np.nan)
    cursor = 2

    def take(context: str) -> tuple[int, str]:
        nonlocal cursor
        if cursor >= len(lines):
            raise ProcarFormatError(f"truncated PROCAR while reading {context}")
        result = lines[cursor]
        cursor += 1
        return result

    seen_k = set()
    for _ in range(nk):
        line_number, text = take("k-point header")
        match = _KPOINT_HEADER.fullmatch(text)
        if match is None:
            raise ProcarFormatError(f"line {line_number}: expected k-point header; collinear sections unsupported")
        ik = int(match[1]) - 1
        if not 0 <= ik < nk or ik in seen_k:
            raise ProcarFormatError(f"line {line_number}: duplicate or out-of-range k-point ID")
        seen_k.add(ik)
        values = _numbers(match.groups()[1:], f"line {line_number} k-point")
        kpoints[ik], weights[ik] = values[:3], values[3]
        seen_b = set()
        for _ in range(nb):
            line_number, text = take("band header")
            match = re.fullmatch(r"band\s+(\d+)\s+# energy\s+(\S+)\s+# occ\.\s+(\S+)", text)
            if match is None:
                raise ProcarFormatError(f"line {line_number}: expected band energy/occupation header")
            band = int(match[1]) - 1
            if not 0 <= band < nb or band in seen_b:
                raise ProcarFormatError(f"line {line_number}: duplicate or out-of-range band ID")
            seen_b.add(band)
            energies[ik, band], occupations[ik, band] = _numbers(
                match.groups()[1:], f"line {line_number} band")
            _, text = take("orbital header")
            if orbital_header(text) != orbitals:
                raise ProcarFormatError("orbital names/order changed between bands")
            for component in range(4):
                seen_ions = set()
                for _ in range(ni):
                    line_number, text = take("ion projection")
                    cells = text.split()
                    if len(cells) != no + 2 or not cells[0].isdigit():
                        raise ProcarFormatError(f"line {line_number}: expected ion row in q/m block {component}")
                    ion = int(cells[0]) - 1
                    if not 0 <= ion < ni or ion in seen_ions:
                        raise ProcarFormatError(f"line {line_number}: duplicate or out-of-range ion ID")
                    seen_ions.add(ion)
                    values = _numbers(cells[1:], f"line {line_number} projection")
                    ion_orbital[ik, band, ion, :, component] = values[:-1]
                    ion_totals[ik, band, ion, component] = values[-1]
                line_number, text = take("tot row")
                cells = text.split()
                if len(cells) != no + 2 or cells[0] != "tot":
                    raise ProcarFormatError(f"line {line_number}: missing tot row in q/m block {component}")
                values = _numbers(cells[1:], f"line {line_number} totals")
                orbital_totals[ik, band, :, component] = values[:-1]
                state_totals[ik, band, component] = values[-1]
    if cursor != len(lines):
        raise ProcarFormatError("unexpected trailing data: duplicate states or unsupported phase/spin blocks")

    arrays = (kpoints, weights, energies, occupations, ion_orbital, ion_totals,
              orbital_totals, state_totals)
    if not all(np.isfinite(a).all() for a in arrays):
        raise ProcarFormatError("missing or nonfinite k/band/ion data")
    # Independently rounded entries and totals can disagree legitimately.
    row_error = float(np.max(np.abs(ion_orbital.sum(axis=3) - ion_totals)))
    ion_error = float(np.max(np.abs(ion_orbital.sum(axis=2) - orbital_totals)))
    state_error = float(np.max(np.abs(ion_totals.sum(axis=2) - state_totals)))
    orbital_error = float(np.max(np.abs(orbital_totals.sum(axis=2) - state_totals)))
    orbital_tolerance = (no + 1) * 0.0005 + 1e-12
    ion_tolerance = (ni + 1) * 0.0005 + 1e-12
    if max(row_error, orbital_error) > orbital_tolerance or max(ion_error, state_error) > ion_tolerance:
        raise ProcarFormatError("projection totals disagree beyond three-decimal rounding bounds")
    return ProcarProjection(
        source=path, kpoint_ids=np.arange(1, nk+1), band_ids=np.arange(1, nb+1),
        ion_ids=np.arange(1, ni+1), kpoints=kpoints, weights=weights,
        energies_eV=energies, occupations=occupations, orbital_names=orbitals,
        ion_orbital=ion_orbital, ion_totals=ion_totals,
        orbital_totals=orbital_totals, state_totals=state_totals,
        spin_to_cartesian=rotation,
        rounding_diagnostics={"max_orbital_sum_vs_ion_total": row_error,
                              "max_ion_sum_vs_orbital_total": ion_error,
                              "max_ion_sum_vs_state_total": state_error,
                              "max_orbital_sum_vs_state_total": orbital_error,
                              "orbital_rounding_tolerance": orbital_tolerance,
                              "ion_rounding_tolerance": ion_tolerance})


def aggregate_groups(data: ProcarProjection, groups: Mapping[str, Sequence[int]], *,
                     axis_cartesian: Sequence[float]) -> dict[str, GroupProjection]:
    """Return raw joint ion-group/spin projections along an explicit unit axis.

    Ion IDs are 1-based. Duplicate IDs inside a group fail. Distinct named groups
    may overlap (e.g. all atoms and one layer); this API does not claim they form
    a partition. Printed ion ``tot`` values are summed without renormalization,
    and small negative p values caused by projection/rounding are preserved.
    """
    axis = np.asarray(axis_cartesian, dtype=float)
    if axis.shape != (3,) or not np.isfinite(axis).all() or not np.isclose(np.linalg.norm(axis), 1, atol=1e-10, rtol=0):
        raise ValueError("axis_cartesian must be an explicit finite unit 3-vector")
    if not isinstance(groups, Mapping) or not groups:
        raise ValueError("at least one explicitly named ion group is required")
    result = {}
    for name, ids in groups.items():
        if not isinstance(name, str) or not name.strip():
            raise ValueError("group names must be nonempty strings")
        try:
            raw_ids = tuple(ids)
        except TypeError as exc:
            raise ValueError(f"group {name}: a sequence of 1-based integer ion IDs required") from exc
        if not raw_ids or any(isinstance(i, (bool, np.bool_)) or not isinstance(i, (int, np.integer)) for i in raw_ids):
            raise ValueError(f"group {name}: nonempty 1-based integer ion IDs required")
        ion_ids = tuple(int(i) for i in raw_ids)
        if len(set(ion_ids)) != len(ion_ids) or min(ion_ids) < 1 or max(ion_ids) > data.nions:
            raise ValueError(f"group {name}: duplicate or out-of-range ion IDs")
        raw = data.ion_totals[:, :, np.array(ion_ids)-1, :].sum(axis=2)
        q = raw[..., 0]
        m = data.pauli_cartesian(raw)
        along = np.einsum("...i,i->...", m, axis)
        result[name] = GroupProjection(name, ion_ids, q, m, along,
                                       (q+along)/2, (q-along)/2, axis.copy())
    return result


def validate_state_alignment(data: ProcarProjection, kpoints: np.ndarray,
                             energies_eV: np.ndarray, *, kpoint_ids=None,
                             band_ids=None, weights=None,
                             energy_atol_eV: float = 5.1e-9,
                             k_atol: float = 5.1e-9,
                             weight_atol: float = 5.1e-9) -> dict[str, float]:
    """Reject mismatched states; no nearest-energy matching/reordering is done.

    References without explicit IDs use canonical ascending 1-based IDs. The
    reference must include every represented k point and band. Fractional k
    coordinates may differ by an exact reciprocal-lattice vector. Eight-place
    PROCAR header rounding sets the default tolerances. Matching alone cannot
    prove equal Hamiltonians or equal degenerate-state gauges; run provenance
    must also agree before combining projection and curvature data.
    """
    for label, value in (("energy", energy_atol_eV), ("k", k_atol), ("weight", weight_atol)):
        if not np.isfinite(value) or value < 0:
            raise ValueError(f"{label} matching tolerance must be finite and nonnegative")
    points = np.asarray(kpoints, dtype=float)
    energies = np.asarray(energies_eV, dtype=float)
    if points.shape != data.kpoints.shape or energies.shape != data.energies_eV.shape:
        raise ValueError("PROCAR/reference k-point or complete band shape mismatch")
    if not np.isfinite(points).all() or not np.isfinite(energies).all():
        raise ValueError("reference k points and band energies must be finite")
    for label, expected, supplied in (("k-point", data.kpoint_ids, kpoint_ids),
                                      ("band", data.band_ids, band_ids)):
        if supplied is not None and not np.array_equal(np.asarray(supplied), expected):
            raise ValueError(f"PROCAR/reference {label} IDs/order mismatch")
    difference = data.kpoints - points
    difference -= np.rint(difference)
    kd = float(np.max(np.abs(difference)))
    ed = float(np.max(np.abs(data.energies_eV-energies)))
    if kd > k_atol or ed > energy_atol_eV:
        raise ValueError(f"PROCAR/reference states mismatch: k delta {kd:.3g}, energy delta {ed:.3g} eV")
    result = {"max_k_fractional_delta": kd, "max_energy_delta_eV": ed}
    if weights is not None:
        reference = np.asarray(weights, dtype=float)
        if reference.shape != data.weights.shape or not np.isfinite(reference).all():
            raise ValueError("reference k weights must be complete and finite")
        wd = float(np.max(np.abs(data.weights-reference)))
        if wd > weight_atol:
            raise ValueError("PROCAR/reference k weights mismatch")
        result["max_weight_delta"] = wd
    return result
