#!/usr/bin/env python3
"""Full-connection occupied-bundle Berry curvature from real-space operators.

The Hamiltonian and Cartesian position matrices must describe the same finite
Wannier subspace. This module evaluates the response itself, without postw90.
Wang et al., Phys. Rev. B 74, 195118 (2006); Lopez et al., Phys. Rev. B
85, 014435 (2012), Eq. (51).
Only occupied-to-empty energy denominators are needed for an insulating bundle.
"""
from __future__ import annotations

import numpy as np

from berry_data import require

AXIAL = ((1, 2), (2, 0), (0, 1))


def prepare_fourier(data):
    """Cartesian derivatives use i R_cart; fractional phases use 2 pi q.R."""
    h = np.asarray(data.hamiltonian_eV)
    a = np.asarray(data.position_A)
    r = np.asarray(data.irvec)
    lattice = np.asarray(data.lattice_A)
    require(h.ndim == 3 and h.shape[1] == h.shape[2] and h.shape[1] > 1,
            'square Hamiltonian with at least two states required')
    nr, nb, _ = h.shape
    require(a.shape == (nr, 3, nb, nb) and r.shape == (nr, 3)
            and lattice.shape == (3, 3) and abs(np.linalg.det(lattice)) > 1e-12,
            'inconsistent position/lattice/translation dimensions')
    require(all(np.isfinite(v).all() for v in (h, a, r, lattice)), 'nonfinite real-space operator')
    require(np.array_equal(r, np.rint(r)), 'integer lattice translations required')
    rc = r @ lattice
    fields = np.empty((nr, 10, nb, nb), dtype=np.complex128)
    fields[:, 0] = h
    fields[:, 1:4] = a
    fields[:, 4:7] = 1j * rc[:, :, None, None] * h[:, None]
    for c, (alpha, beta) in enumerate(AXIAL):
        fields[:, 7+c] = 1j * (rc[:, alpha, None, None]*a[:, beta]
                                - rc[:, beta, None, None]*a[:, alpha])
    return fields


def occupied_curvature(data, q, occupied, gap_threshold=1e-8, fields=None):
    """Return the gauge-invariant trace, retaining full J0, J1 and J2 terms.

    D_a = <occupied|dH/dk_a|empty> / (E_empty-E_occupied).
    Rotations within degenerate occupied or empty spaces leave these traces
    unchanged. Internal denominators and gauge-dependent band curvatures are
    never evaluated. A gap closure across the selected bundle is an error.
    """
    q = np.asarray(q, dtype=float)
    nb = data.hamiltonian_eV.shape[-1]
    require(q.ndim == 2 and q.shape[1] == 3 and len(q) and np.isfinite(q).all(),
            'nonempty finite fractional k points required')
    require(type(occupied) is int and 0 < occupied < nb, 'occupied count must lie strictly inside the model')
    require(np.isfinite(gap_threshold) and gap_threshold > 0, 'positive finite cross-gap threshold required')
    fields = prepare_fourier(data) if fields is None else fields
    require(fields.shape == (len(data.irvec), 10, nb, nb), 'inconsistent prepared Fourier fields')
    phase = np.exp(2j*np.pi*(q @ data.irvec.T))
    value = (phase @ fields.reshape(len(data.irvec), -1)).reshape(len(q), 10, nb, nb)
    defect = np.max(np.abs(value-value.conj().swapaxes(-1, -2)), axis=(-2, -1))
    scale = np.max(np.abs(value), axis=(-2, -1))
    require(np.all(defect <= 1e-8 + 1e-10*scale), 'Fourier operators are not Hermitian; check real-space input')
    energy, unitary = np.linalg.eigh(value[:, 0])
    gap = energy[:, occupied]-energy[:, occupied-1]
    require(np.all(gap > gap_threshold), 'selected occupied bundle touches excluded states')
    occ, emp = unitary[:, :, :occupied], unitary[:, :, occupied:]
    occ_h = occ.conj().swapaxes(-1, -2)
    cross_a = occ_h[:, None] @ (value[:, 1:4] @ emp[:, None])
    cross_dh = occ_h[:, None] @ (value[:, 4:7] @ emp[:, None])
    denominator = energy[:, None, occupied:]-energy[:, :occupied, None]
    cross_d = cross_dh/denominator[:, None]
    terms = np.empty((len(q), 3, 3), dtype=float)
    # J0 is the ordinary curl of A^W, with no extra connection commutator.
    terms[:, 0] = np.einsum('kni,kanm,kmi->ka', occ.conj(), value[:, 7:10],
                             occ, optimize=True).real
    for c, (alpha, beta) in enumerate(AXIAL):
        terms[:, 1, c] = 2*np.real(np.sum(
            cross_a[:, alpha]*cross_d[:, beta].conj()
            - cross_d[:, alpha]*cross_a[:, beta].conj(), axis=(1, 2)))
        terms[:, 2, c] = -2*np.imag(np.sum(
            cross_d[:, alpha]*cross_d[:, beta].conj(), axis=(1, 2)))
    require(np.isfinite(energy).all() and np.isfinite(terms).all(), 'nonfinite full-connection response')
    return dict(energies_eV=energy, omega_terms_A2=terms, omega_A2=terms.sum(axis=1),
                minimum_gap_eV=float(gap.min()), max_hermiticity_defect=float(defect.max()))
