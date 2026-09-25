"""Independent physical checks of the conventional occupied spin-Kubo response."""
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))
from spin_hall import (occupied_spin_curvature, projected_spin_current,
                       integrate_sheet, tensor_component, velocity_from_m_per_s,
                       HBAR_EV_S, CONDUCTANCE_QUANTUM_S)

X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1., -1.]).astype(complex)
I = np.eye(2, dtype=complex)


def orbital(kx, ky, mass=-1.):
    h = np.sin(kx)*X + np.sin(ky)*Y + (mass+np.cos(kx)+np.cos(ky))*Z
    d = np.array([np.cos(kx)*X-np.sin(kx)*Z,
                  np.cos(ky)*Y-np.sin(ky)*Z, np.zeros((2, 2))])
    return h, d


def bands_and_operators(h, d, s):
    energy, u = np.linalg.eigh(h)
    uh = u.conj().swapaxes(-1, -2)
    return energy, uh[:, None] @ d @ u[:, None], uh[:, None] @ s @ u[:, None]


def qsh(points, same_spin_blocks=False):
    hs, ds = [], []
    for kx, ky in points:
        hu, du = orbital(kx, ky)
        hm, dm = orbital(-kx, -ky)
        hd, dd = (hu, du) if same_spin_blocks else (hm.conj(), -dm.conj())
        h = np.zeros((4, 4), complex)
        d = np.zeros((3, 4, 4), complex)
        h[:2, :2], h[2:, 2:] = hu, hd
        d[:, :2, :2], d[:, 2:, 2:] = du, dd
        hs.append(h); ds.append(d)
    spin = np.array([np.kron(pauli, I) for pauli in [X, Y, Z]])
    return bands_and_operators(np.array(hs), np.array(ds), spin[None])


def lower_analytic(kx, ky):
    radius = np.sqrt(np.sin(kx)**2+np.sin(ky)**2+(-1+np.cos(kx)+np.cos(ky))**2)
    return (np.cos(kx)+np.cos(ky)-np.cos(kx)*np.cos(ky))/(2*radius**3)


def lower_plaquette(kx, ky, step=2e-4):
    states = [np.linalg.eigh(orbital(kx+dx*step, ky+dy*step)[0])[1][:, 0]
              for dx, dy in [(-.5, -.5), (.5, -.5), (.5, .5), (-.5, .5)]]
    links = np.prod([np.vdot(states[i], states[(i+1) % 4]) for i in range(4)])
    return -np.angle(links)/step**2


def mixed_three_state(occupied_energy=(-2., -1.)):
    # Every diagonal spin is zero. Spin mixing within the occupied space
    # couples through Dx to the empty state; mean-spin*Berry misses it.
    energy = np.array([[*occupied_energy, 1.]])
    d = np.zeros((1, 3, 3, 3), complex)
    d[0, 0, 1, 2] = d[0, 0, 2, 1] = 1
    d[0, 1, 0, 2], d[0, 1, 2, 0] = 1j, -1j
    s = np.zeros_like(d)
    s[0, 2, 0, 1] = s[0, 2, 1, 0] = 1
    return energy, d, s


def calculate(e, d, s, occupied=2, **kwargs):
    return occupied_spin_curvature(e, d, occupied, spin_pauli=s,
                                    input_method='exact finite test model', **kwargs)


class SpinHallPhysicsTests(unittest.TestCase):
    def test_spin_conserving_qsh_matches_analytic_and_physical_links(self):
        points = np.array([[.21, -.37], [1.12, .44], [-1.24, 1.1]])
        e, d, s = qsh(points)
        result = calculate(e, d, s)
        expected = np.array([2*lower_analytic(*p) for p in points])
        np.testing.assert_allclose(result['omega_spin_A2'][:, 2, 0, 1], expected, atol=2e-14)
        plaquettes = np.array([2*lower_plaquette(*p) for p in points])
        np.testing.assert_allclose(expected, plaquettes, rtol=2e-7, atol=3e-8)
        np.testing.assert_allclose(result['omega_charge_A2'], 0, atol=2e-14)
        np.testing.assert_allclose(result['omega_spin_A2'][:, 2, 1, 0], -expected, atol=2e-14)

    def test_qsh_sheet_integral_sign_and_physical_spin_half(self):
        n = 41
        a = 2*np.pi*(np.arange(n)+.5)/n
        points = np.array(np.meshgrid(a, a)).reshape(2, -1).T
        e, d, s = qsh(points)
        result = calculate(e, d, s)
        sheet = integrate_sheet(result['omega_spin_A2'], np.full(n*n, 1/(n*n)), (2*np.pi)**2)
        # C_up=+1 and C_down=-1, so charge cancels while conventional spin
        # sigma is +1 in (hbar/e)*(e^2/h), not +2 or -1.
        self.assertAlmostEqual(sheet['sigma_hbar_over_e_e2_over_h'][2, 0, 1], 1., places=9)
        self.assertAlmostEqual(sheet['sigma_hbar_over_e_S'][2, 0, 1], CONDUCTANCE_QUANTUM_S, places=13)

    def test_single_up_spin_current_has_opposite_sign_to_charge_current(self):
        h, d = orbital(.4, .8)
        spin = np.array([np.zeros((2, 2)), np.zeros((2, 2)), I])
        e, d, s = bands_and_operators(h[None], d[None], spin[None])
        result = calculate(e, d, s, occupied=1)
        sheet = integrate_sheet(result['omega_spin_A2'], [1.], 3.)
        charge = -3/(2*np.pi)*result['omega_charge_A2'][0, 0, 1]
        self.assertAlmostEqual(sheet['sigma_hbar_over_e_e2_over_h'][2, 0, 1], -.5*charge, places=14)

    def test_zero_soc_spin_degenerate_blocks_have_zero_spin_response(self):
        e, d, s = qsh([[.3, .8], [.7, -.2]], same_spin_blocks=True)
        result = calculate(e, d, s)
        np.testing.assert_allclose(result['omega_spin_A2'], 0., atol=2e-14)
        self.assertGreater(np.max(np.abs(result['omega_charge_A2'])), .1)

    def test_offdiagonal_spin_and_all_intermediate_product_states_are_required(self):
        e, d, s = mixed_three_state()
        result = calculate(e, d, s)
        self.assertTrue(np.all(np.diagonal(s, axis1=-2, axis2=-1) == 0))
        # K_zx[0,2]=1/2 from the occupied intermediate state 1.
        # Direct two-amplitude algebra gives -2 Im[(1/2)*(-i)]/3^2=1/9.
        self.assertAlmostEqual(result['omega_spin_A2'][0, 2, 0, 1], 1/9, places=15)
        # Spin current need not produce an antisymmetric current/field tensor.
        self.assertAlmostEqual(result['omega_spin_A2'][0, 2, 1, 0], -1/4, places=15)
        np.testing.assert_allclose(result['omega_charge_A2'], 0, atol=1e-15)
        self.assertEqual(result['spin_current_method'], 'finite_band_projected_product')

    def test_supplied_current_preserves_term_lost_by_subspace_product(self):
        e, d, s = mixed_three_state()
        full_current = projected_spin_current(s, d)
        keep = [0, 2]
        er = e[:, keep]
        dr = d[:, :, keep, :][:, :, :, keep]
        sr = s[:, :, keep, :][:, :, :, keep]
        jr = full_current[:, :, :, keep, :][:, :, :, :, keep]
        direct = occupied_spin_curvature(er, dr, 1, spin_current_eVA=jr,
                                         input_method='current projected from complete three-state toy')
        product = calculate(er, dr, sr, occupied=1)
        self.assertAlmostEqual(direct['omega_spin_A2'][0, 2, 0, 1], 1/9, places=15)
        self.assertEqual(product['omega_spin_A2'][0, 2, 0, 1], 0.)
        self.assertEqual(direct['spin_current_method'], 'provided_spin_current')

    def test_internal_exact_degeneracy_and_complex_gauge_rotations(self):
        e, d, s = qsh([[.3, .7], [.9, -.4]])
        result = calculate(e, d, s)
        rng = np.random.default_rng(12)
        unitary = np.zeros((2, 4, 4), complex)
        for k in range(2):
            for offset in [0, 2]:
                z = rng.normal(size=(2, 2))+1j*rng.normal(size=(2, 2))
                unitary[k, offset:offset+2, offset:offset+2] = np.linalg.qr(z)[0]
        uh = unitary.conj().swapaxes(-1, -2)
        changed = calculate(e, uh[:, None] @ d @ unitary[:, None],
                             uh[:, None] @ s @ unitary[:, None])
        np.testing.assert_allclose(changed['omega_spin_A2'], result['omega_spin_A2'], atol=4e-14)
        # Nondegenerate band phases also preserve all complex off-diagonal data.
        e, d, s = mixed_three_state()
        u = np.diag(np.exp(1j*np.array([.2, 1.4, -.7])))
        changed = calculate(e, u.conj().T @ d @ u, u.conj().T @ s @ u)
        np.testing.assert_allclose(changed['omega_spin_A2'], calculate(e, d, s)['omega_spin_A2'], atol=1e-15)

    def test_time_reversal_spin_even_charge_odd_pair(self):
        e, d, s = mixed_three_state()
        original = calculate(e, d, s)
        # In the partner eigenbasis, T-odd Hermitian operators transform
        # as -O*. Spin current is even because both factors are T odd.
        partner = calculate(e, -d.conj(), -s.conj())
        np.testing.assert_allclose(partner['omega_spin_A2'], original['omega_spin_A2'], atol=1e-15)
        # Add a same-pair velocity to make a nonzero charge reference.
        d[0, 0, 0, 2] = d[0, 0, 2, 0] = .4
        first = calculate(e, d, s)
        second = calculate(e, -d.conj(), -s.conj())
        self.assertGreater(abs(first['omega_charge_A2'][0, 0, 1]), .01)
        np.testing.assert_allclose(second['omega_charge_A2'], -first['omega_charge_A2'], atol=1e-15)
        np.testing.assert_allclose(second['omega_spin_A2'], first['omega_spin_A2'], atol=1e-15)

    def test_units_normalized_weights_and_directed_axis_contraction(self):
        e, d, s = mixed_three_state()
        physical = d/(HBAR_EV_S*1e10)
        np.testing.assert_allclose(velocity_from_m_per_s(physical), d, atol=1e-15)
        result = calculate(e, d, s)
        sheet = integrate_sheet(result['omega_spin_A2'], [1.], 2.5)
        tensor = sheet['sigma_hbar_over_e_e2_over_h']
        x, y, z = np.eye(3)
        self.assertEqual(tensor_component(tensor, z, x, y), tensor[2, 0, 1])
        self.assertEqual(tensor_component(tensor, -z, x, y), -tensor[2, 0, 1])
        self.assertEqual(tensor_component(tensor, z, -x, y), -tensor[2, 0, 1])
        with self.assertRaisesRegex(ValueError, 'sum to one'):
            integrate_sheet(result['omega_spin_A2'], [2.], 2.5)
        with self.assertRaisesRegex(ValueError, 'unit directions'):
            tensor_component(tensor, z, 2*x, y)

    def test_gap_metal_mu_and_nonorthonormal_input_are_rejected(self):
        e, d, s = mixed_three_state()
        with self.assertRaisesRegex(ValueError, 'touch'):
            calculate(np.array([[-2., 1., 1.]]), d, s)
        e2 = np.concatenate([e, e+5], axis=0)
        with self.assertRaisesRegex(ValueError, 'global gap'):
            calculate(e2, np.repeat(d, 2, axis=0), np.repeat(s, 2, axis=0))
        for mu in [-1., 1., np.nan]:
            with self.assertRaisesRegex(ValueError, 'chemical potentials'):
                calculate(e, d, s, mu_eV=mu)
        with self.assertRaisesRegex(ValueError, 'orthonormal'):
            calculate(e, d, s, basis_overlap=.8*np.eye(3)[None])
        self.assertTrue(calculate(e, d, s, basis_overlap=np.eye(3)[None])['basis_overlap_checked'])
        defect = 2e-6
        overlap = (1+defect)*np.eye(3)[None]
        with self.assertRaisesRegex(ValueError, 'orthonormal'):
            calculate(e, d, s, basis_overlap=overlap)
        record = calculate(e, d, s, basis_overlap=overlap, basis_overlap_tolerance=5e-5)
        self.assertAlmostEqual(record['basis_overlap_max_defect'], defect, places=14)
        self.assertEqual(record['basis_overlap_tolerance'], 5e-5)

    def test_missing_provenance_and_nonhermitian_matrices_are_rejected(self):
        e, d, s = mixed_three_state()
        with self.assertRaisesRegex(ValueError, 'input_method'):
            occupied_spin_curvature(e, d, 2, spin_pauli=s, input_method='')
        with self.assertRaisesRegex(ValueError, 'must be real'):
            calculate(e+0j, d, s)
        d[0, 0, 1, 2] = 2j
        with self.assertRaisesRegex(ValueError, 'non-Hermitian'):
            calculate(e, d, s)


if __name__ == '__main__':
    unittest.main()
