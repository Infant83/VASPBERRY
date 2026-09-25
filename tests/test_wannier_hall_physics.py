"""Independent spinor/projector checks of full-connection Wannier curvature.

The finite-plaquette oracle uses physical three-component wavefunctions,
not the solver's derivative or J0/J1/J2 formulas.
"""
from pathlib import Path
from types import SimpleNamespace
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from wannier_hall import occupied_curvature, prepare_fourier

SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.diag([1., -1.]).astype(complex)


def model_from_records(h, a, lattice=None):
    keys = sorted(set(h) | set(a))
    n = next(iter(h.values())).shape[0]
    return SimpleNamespace(
        irvec=np.array(keys, dtype=int),
        hamiltonian_eV=np.array([h.get(r, np.zeros((n, n), complex)) for r in keys]),
        position_A=np.array([a.get(r, np.zeros((3, n, n), complex)) for r in keys]),
        lattice_A=np.eye(3) if lattice is None else np.asarray(lattice, float),
    )


def qwz(mass=-1.2, axes=(0, 1), lattice=None, embedded=False):
    """H=sin(k_a)σx+sin(k_b)σy+(m+cos(k_a)+cos(k_b))σz."""
    origin = (0, 0, 0)
    h = {origin: mass*SZ}
    for axis, pauli in zip(axes, (SX, SY)):
        r = np.eye(3, dtype=int)[axis]
        h[tuple(r)] = (SZ-1j*pauli)/2
        h[tuple(-r)] = (SZ+1j*pauli)/2
    a = {}
    if embedded:
        # Basis columns are v1=e^(ikx/2)(cos(kx/2),0,e^iky sin(kx/2))
        # and v2=(0,1,0). Their connection is Ax_11=-1/2,
        # Ay_11=(cos(kx)-1)/2; this is not a pure basis-gauge shift.
        a[origin] = np.zeros((3, 2, 2), complex)
        a[origin][0, 0, 0] = -.5
        a[origin][1, 0, 0] = -.5
        for sign in (-1, 1):
            a[(sign, 0, 0)] = np.zeros((3, 2, 2), complex)
            a[(sign, 0, 0)][1, 0, 0] = .25
    return model_from_records(h, a, lattice)


def analytic_lower_curvature(kx, ky, mass=-1.2):
    d = np.array([np.sin(kx), np.sin(ky), mass+np.cos(kx)+np.cos(ky)])
    # With A=i<u|du>, the lower eigenstate has +d·(dx d × dy d)/(2|d|³).
    numerator = np.cos(kx)+np.cos(ky)+mass*np.cos(kx)*np.cos(ky)
    return numerator/(2*np.linalg.norm(d)**3)


def embedded_lower_state(kx, ky, mass=-1.2):
    h = np.sin(kx)*SX + np.sin(ky)*SY + (mass+np.cos(kx)+np.cos(ky))*SZ
    _, u = np.linalg.eigh(h)
    basis = np.array([[np.exp(.5j*kx)*np.cos(.5*kx), 0],
                      [0, 1],
                      [np.exp(.5j*kx+1j*ky)*np.sin(.5*kx), 0]], complex)
    return basis @ u[:, 0]


def plaquette_curvature(state, kx, ky, step=2e-4):
    # Centered physical overlaps with counter-clockwise boundary; A=i<u|du>
    # gives minus the overlap-product phase divided by Cartesian area.
    offsets = [(-.5, -.5), (.5, -.5), (.5, .5), (-.5, .5)]
    states = [state(kx+step*x, ky+step*y) for x, y in offsets]
    product = np.prod([np.vdot(states[i], states[(i+1) % 4]) for i in range(4)])
    return -np.angle(product)/step**2


def integer_basis_shift(model, shifts):
    """Apply W(q)=diag exp(+2πiq·t): H'=W†HW, A'=W†AW-t_cart."""
    shifts = np.asarray(shifts, int)
    n = model.hamiltonian_eV.shape[1]
    h, a = {}, {}
    for r, hr, ar in zip(model.irvec, model.hamiltonian_eV, model.position_A):
        for i in range(n):
            for j in range(n):
                moved = tuple(r + shifts[j] - shifts[i])
                h.setdefault(moved, np.zeros((n, n), complex))[i, j] += hr[i, j]
                a.setdefault(moved, np.zeros((3, n, n), complex))[:, i, j] += ar[:, i, j]
    zero = a.setdefault((0, 0, 0), np.zeros((3, n, n), complex))
    cart = shifts @ model.lattice_A
    for i in range(n):
        zero[:, i, i] -= cart[i]
    return model_from_records(h, a, model.lattice_A)


class WannierHallPhysicsTests(unittest.TestCase):
    def test_qwz_analytic_curvature_and_cartesian_axis_signs(self):
        q = np.array([[.07, .13, .19], [.24, .31, .11], [.41, .08, .23]])
        for axes, component, sign in [((0, 1), 2, 1), ((1, 2), 0, 1), ((0, 2), 1, -1)]:
            with self.subTest(axes=axes):
                out = occupied_curvature(qwz(axes=axes), q, 1)
                expected = np.zeros((len(q), 3))
                expected[:, component] = [sign*analytic_lower_curvature(*(2*np.pi*p[list(axes)])) for p in q]
                np.testing.assert_allclose(out["omega_A2"], expected, rtol=2e-12, atol=2e-12)
                np.testing.assert_allclose(out["omega_terms_A2"][:, :2], 0, atol=1e-14)

    def test_cartesian_derivatives_use_real_lattice_not_fractional_q(self):
        q = np.array([[.13, .27, 0], [.31, .41, 0]])
        lattice = np.diag([2., 3., 5.])
        out = occupied_curvature(qwz(lattice=lattice), q, 1)
        expected = [6*analytic_lower_curvature(*(2*np.pi*p[:2])) for p in q]
        np.testing.assert_allclose(out["omega_A2"][:, 2], expected, rtol=2e-12, atol=2e-12)

    def test_full_connection_matches_physical_wavefunction_plaquettes(self):
        points = np.array([[.71, -.46], [1.2, .8], [-.9, 1.3]])
        q = np.column_stack([points/(2*np.pi), np.zeros(len(points))])
        out = occupied_curvature(qwz(embedded=True), q, 1)
        expected = [plaquette_curvature(embedded_lower_state, *k) for k in points]
        np.testing.assert_allclose(out["omega_A2"][:, 2], expected, rtol=3e-7, atol=3e-8)
        # This physical embedding activates the connection curl and the mixed
        # term; an H-only backend cannot accidentally satisfy this oracle.
        self.assertGreater(np.max(np.abs(out["omega_terms_A2"][:, 0, 2])), .05)
        self.assertGreater(np.max(np.abs(out["omega_terms_A2"][:, 1, 2])), .05)
        self.assertGreater(np.max(np.abs(out["omega_terms_A2"][:, 2, 2])), .05)
        h_only = occupied_curvature(qwz(), q, 1)["omega_A2"][:, 2]
        self.assertGreater(np.max(np.abs(h_only-expected)), .05)

    def test_k_dependent_basis_shift_preserves_physical_curvature(self):
        model = qwz(embedded=True)
        moved = integer_basis_shift(model, [[1, -1, 0], [-1, 2, 0]])
        q = np.array([[.07, .13, 0], [.24, .31, 0], [.41, .08, 0]])
        baseline = occupied_curvature(model, q, 1)
        out = occupied_curvature(moved, q, 1)
        np.testing.assert_allclose(out["energies_eV"], baseline["energies_eV"], atol=2e-14)
        np.testing.assert_allclose(out["omega_A2"], baseline["omega_A2"], rtol=1e-11, atol=1e-11)
        # Individual J terms depend on the chosen representation.
        self.assertGreater(np.max(np.abs(out["omega_terms_A2"]-baseline["omega_terms_A2"])), .05)

    def test_exact_internal_degeneracy_is_regular_for_occupied_bundle(self):
        original = qwz()
        rng = np.random.default_rng(29)
        z = rng.normal(size=(4, 4))+1j*rng.normal(size=(4, 4))
        rotation, _ = np.linalg.qr(z)
        h = np.array([rotation.conj().T @ np.kron(np.eye(2), r) @ rotation
                      for r in original.hamiltonian_eV])
        model = SimpleNamespace(irvec=original.irvec, hamiltonian_eV=h,
                                position_A=np.zeros((len(h), 3, 4, 4), complex), lattice_A=np.eye(3))
        q = np.array([[.09, .17, 0], [.31, .28, 0]])
        out = occupied_curvature(model, q, 2)
        expected = [2*analytic_lower_curvature(*(2*np.pi*p[:2])) for p in q]
        np.testing.assert_allclose(out["omega_A2"][:, 2], expected, rtol=2e-12, atol=2e-12)
        np.testing.assert_allclose(out["energies_eV"][:, 0], out["energies_eV"][:, 1], atol=2e-14)

    def test_closed_occupied_empty_gap_is_rejected(self):
        with self.assertRaises(ValueError):
            occupied_curvature(qwz(mass=-2), np.zeros((1, 3)), 1, gap_threshold=1e-8)

    def test_periodic_images_and_precomputed_fields_match_analytic_values(self):
        model = qwz()
        fields = prepare_fourier(model)
        q = np.array([[.13, .27, .11], [1.13, -1.73, 2.11]])
        out = occupied_curvature(model, q, 1, fields=fields)
        expected = analytic_lower_curvature(*(2*np.pi*q[0, :2]))
        np.testing.assert_allclose(out["omega_A2"][:, 2], expected, atol=2e-12)


if __name__ == "__main__":
    unittest.main()
