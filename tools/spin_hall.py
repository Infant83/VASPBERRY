#!/usr/bin/env python3
"""Conventional intrinsic spin Hall response of a gapped 2D occupied bundle.

Input matrices are in the same orthonormal eigenstate basis. ``velocity_eVA``
means D = hbar*v, in eV Angstrom, not v in m/s. ``spin_pauli`` is dimensionless
sigma, so physical spin is hbar*sigma/2. A supplied ``spin_current_eVA`` is
K = {sigma, D}/2, with axes [k, spin, current_direction, band, band].

The optional product of finite-band spin and velocity matrices is explicitly
a projected-product approximation. It misses P s Q v P and P v Q s P terms.
Providing K avoids that extra product approximation; its physical provenance
and the convergence of the represented intermediate-band sum remain external
requirements. This module never infers PAW completeness from a matrix's shape.

References: Ryoo et al., PRB 99, 235113 (2019), Eqs. (1)-(3), Sec. II D;
Qiao et al., PRB 98, 214402 (2018). No metal, finite-T or lifetime model here.
"""
from __future__ import annotations

import numpy as np

ELEMENTARY_CHARGE_C = 1.602176634e-19
PLANCK_J_S = 6.62607015e-34
HBAR_EV_S = PLANCK_J_S / (2*np.pi*ELEMENTARY_CHARGE_C)
CONDUCTANCE_QUANTUM_S = ELEMENTARY_CHARGE_C**2 / PLANCK_J_S


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _hermitian(value, shape, label, atol=1e-8, rtol=1e-10):
    result = np.asarray(value, dtype=np.complex128)
    _require(result.shape == shape and np.isfinite(result).all(),
             f'{label}: expected finite matrix array with shape {shape}')
    defect = np.max(np.abs(result-result.conj().swapaxes(-1, -2)), axis=(-2, -1))
    scale = np.max(np.abs(result), axis=(-2, -1))
    _require(np.all(defect <= atol + rtol*scale), f'{label}: non-Hermitian input')
    return result, float(defect.max())


def velocity_from_m_per_s(velocity_m_per_s):
    """Convert physical velocity to hbar*v in eV Angstrom, without normalization."""
    value = np.asarray(velocity_m_per_s, dtype=np.complex128)
    _require(np.isfinite(value).all(), 'nonfinite velocity')
    return value * (HBAR_EV_S*1e10)


def projected_spin_current(spin_pauli, velocity_eVA):
    """Return K={sigma,D}/2 within the supplied finite band space.

    Band multiplication must precede selecting occupied/empty pairs. This is
    not automatically the full physical P{sigma,D}P/2 operator.
    """
    d = np.asarray(velocity_eVA)
    _require(d.ndim == 4 and d.shape[1] == 3 and d.shape[2] == d.shape[3]
             and d.shape[0] > 0 and d.shape[2] > 1,
             'velocity_eVA must have shape [nk,3,nb,nb]')
    d, _ = _hermitian(d, d.shape, 'velocity_eVA')
    s, _ = _hermitian(spin_pauli, d.shape, 'spin_pauli')
    return .5*(s[:, :, None] @ d[:, None, :] + d[:, None, :] @ s[:, :, None])


def occupied_spin_curvature(energies_eV, velocity_eVA, occupied, *,
                            input_method, spin_pauli=None, spin_current_eVA=None,
                            basis_overlap=None, gap_threshold_eV=1e-8,
                            mu_eV=None, basis_overlap_tolerance=1e-7):
    """Compute the T=0 insulating conventional spin-Kubo tensor.

    Returns omega_spin_A2[k,spin,current,electric_field]. Only occupied-empty
    denominators are evaluated: equal-occupation terms cancel pairwise before
    division, including exactly degenerate internal states. A touching selected
    boundary or nonpositive sampled indirect gap is rejected. ``mu_eV``, when
    supplied, must lie strictly inside the sampled gap.

    D=hbar*v; K={sigma,D}/2 (or its independently supplied full matrix):
        omega_spin = -2 Im sum_(n occupied,m empty) K_nm D_mn/(Em-En)^2.
    Spin sheet response in (hbar/e)*(e^2/h) is +A_BZ/(4*pi)*sum(w*omega_spin).
    Charge response in e^2/h uses -A_BZ/(2*pi)*sum(w*omega_charge).
    Thus for conserved sigma=+1, the normalized spin response is minus one
    half of the normalized charge response. No extra spin multiplicity enters.
    """
    _require(not np.iscomplexobj(energies_eV), 'energies_eV must be real')
    energy = np.asarray(energies_eV, dtype=float)
    _require(energy.ndim == 2 and energy.shape[0] > 0 and energy.shape[1] > 1
             and np.isfinite(energy).all(), 'finite energies [nk,nb] required')
    nk, nb = energy.shape
    _require(type(occupied) is int and 0 < occupied < nb,
             'occupied count must lie strictly inside the supplied band space')
    _require(isinstance(input_method, str) and bool(input_method.strip()),
             'input_method must describe the operator producer/approximation')
    _require(np.isfinite(gap_threshold_eV) and gap_threshold_eV > 0,
             'positive finite gap_threshold_eV required')
    _require(np.all(np.diff(energy, axis=1) >= 0), 'energies must be sorted ascending')
    _require((spin_pauli is None) != (spin_current_eVA is None),
             'supply exactly one of spin_pauli or spin_current_eVA')
    d, d_defect = _hermitian(velocity_eVA, (nk, 3, nb, nb), 'velocity_eVA')
    _require(np.ndim(basis_overlap_tolerance) == 0 and np.isfinite(basis_overlap_tolerance)
             and basis_overlap_tolerance > 0, 'positive finite basis_overlap_tolerance required')
    overlap_defect = None
    if basis_overlap is not None:
        overlap, _ = _hermitian(basis_overlap, (nk, nb, nb), 'basis_overlap')
        overlap_defect = float(np.max(np.abs(overlap-np.eye(nb))))
        _require(overlap_defect <= basis_overlap_tolerance,
                 'orthonormal eigenstate basis required; do not normalize pseudo matrices silently')
    gap = energy[:, occupied]-energy[:, occupied-1]
    _require(np.all(gap > gap_threshold_eV), 'occupied and empty states touch at the selected boundary')
    vbm, cbm = float(energy[:, occupied-1].max()), float(energy[:, occupied].min())
    _require(cbm-vbm > gap_threshold_eV,
             'positive sampled global gap required; metallic occupations are unsupported')
    if mu_eV is not None:
        mu = np.asarray(mu_eV, dtype=float)
        _require(mu.size > 0 and np.isfinite(mu).all() and np.all((mu > vbm) & (mu < cbm)),
                 'all chemical potentials must lie strictly inside the sampled gap')
    if spin_pauli is not None:
        s, s_defect = _hermitian(spin_pauli, (nk, 3, nb, nb), 'spin_pauli')
        current = projected_spin_current(s, d)
        method = 'finite_band_projected_product'
        limitation = 'P{sigma,D}P/2 differs in general from {P sigma P,P D P}/2; source-band closure must converge.'
    else:
        current, s_defect = _hermitian(spin_current_eVA, (nk, 3, 3, nb, nb), 'spin_current_eVA')
        method = 'provided_spin_current'
        limitation = 'Full spin-current accuracy is a producer claim to validate; finite intermediate-band convergence remains required.'
    de = energy[:, None, occupied:]-energy[:, :occupied, None]
    inv_de2 = 1/de**2
    cross_d = d[:, :, :occupied, occupied:]
    cross_j = current[:, :, :, :occupied, occupied:]
    omega = -2*np.imag(np.einsum('ksinm,kjnm,knm->ksij', cross_j,
                                 cross_d.conj(), inv_de2, optimize=True))
    charge = -2*np.imag(np.einsum('kinm,kjnm,knm->kij', cross_d,
                                  cross_d.conj(), inv_de2, optimize=True))
    _require(np.isfinite(omega).all() and np.isfinite(charge).all(), 'nonfinite Kubo response')
    return dict(omega_spin_A2=omega, omega_charge_A2=charge,
                minimum_cross_gap_eV=float(gap.min()), global_gap_eV=cbm-vbm,
                valence_max_eV=vbm, conduction_min_eV=cbm,
                occupied=occupied, source_band_count=nb,
                spin_current_method=method, input_method=input_method.strip(),
                spin_current_limitation=limitation,
                spin_operator_convention='physical s=(hbar/2)*spin_pauli',
                velocity_convention='velocity_eVA=hbar*physical_velocity in eV Angstrom',
                spin_current_convention='spin_current_eVA={spin_pauli,velocity_eVA}/2',
                tensor_axis_order=['spin', 'current', 'electric_field'],
                basis_contract='same orthonormal eigenstate basis; producer association is external',
                basis_overlap_checked=basis_overlap is not None,
                basis_overlap_max_defect=overlap_defect,
                basis_overlap_tolerance=float(basis_overlap_tolerance),
                max_input_hermiticity_defect=max(d_defect, s_defect),
                temperature_K=0.0, metallic_occupations_supported=False,
                integer_rounding_applied=False)


def integrate_sheet(omega_spin_A2, weights, bz_area_inv_A2):
    """Integrate all Cartesian tensor components over a normalized 2D BZ.

    Weights must be nonnegative and sum to one; they are never renormalized.
    ``sigma_hbar_over_e_e2_over_h`` is dimensionless; ``sigma_hbar_over_e_S``
    is its value multiplied by e^2/h. These are the SAME physical coefficient,
    expressed in units (hbar/e)*(e^2/h) and (hbar/e)*S, respectively.
    Sampling completeness and k-point/band association belong to the caller.
    """
    omega = np.asarray(omega_spin_A2, dtype=float)
    w = np.asarray(weights, dtype=float)
    _require(omega.ndim == 4 and omega.shape[1:] == (3, 3, 3)
             and len(omega) > 0 and np.isfinite(omega).all(),
             'finite omega_spin_A2[nk,3,3,3] required')
    _require(w.shape == (len(omega),) and np.isfinite(w).all()
             and np.all(w >= 0), 'finite nonnegative weights[nk] required')
    _require(abs(float(w.sum())-1) <= 1e-10, 'weights must sum to one; no automatic normalization')
    _require(np.ndim(bz_area_inv_A2) == 0 and np.isfinite(bz_area_inv_A2)
             and bz_area_inv_A2 > 0, 'positive finite 2D BZ area required')
    response = float(bz_area_inv_A2)/(4*np.pi)*np.einsum('k,ksij->sij', w, omega)
    return dict(sigma_hbar_over_e_e2_over_h=response,
                sigma_hbar_over_e_S=response*CONDUCTANCE_QUANTUM_S,
                bz_area_inv_A2=float(bz_area_inv_A2), weight_sum=float(w.sum()),
                tensor_axis_order=['spin', 'current', 'electric_field'],
                conductance_quantum_S=CONDUCTANCE_QUANTUM_S,
                formula='sigma/[(hbar/e)*(e^2/h)]=+A_BZ/(4*pi)*sum_k w_k Omega_spin',
                spin_multiplicity=1, dimensionality=2, integer_rounding_applied=False)


def tensor_component(tensor, spin_direction, current_direction, field_direction):
    """Contract a [...,spin,current,field] tensor along three Cartesian unit axes.

    Explicit directed axes preserve signs under axis reversal. Directions must
    already be unit vectors; this helper does not rescale or choose orientations.
    """
    value = np.asarray(tensor, dtype=float)
    _require(value.ndim >= 3 and value.shape[-3:] == (3, 3, 3)
             and np.isfinite(value).all(), 'finite tensor with trailing shape [3,3,3] required')
    directions = [np.asarray(v, dtype=float) for v in
                  (spin_direction, current_direction, field_direction)]
    _require(all(v.shape == (3,) and np.isfinite(v).all()
                 and abs(np.linalg.norm(v)-1) <= 1e-10 for v in directions),
             'explicit Cartesian unit directions required')
    return np.einsum('...sij,s,i,j->...', value, *directions)
