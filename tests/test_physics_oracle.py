"""Independent physics oracles for the valley-transport postprocessor."""

import math
import unittest

import numpy as np

import vaspberry_transport as vt


TWO_PI = 2.0 * math.pi
B1 = np.array([TWO_PI, 0.0, 0.0])
B2 = np.array([0.0, TWO_PI, 0.0])


def _qwz_lower_energy(kx, ky, mass):
    norm = np.sqrt(
        np.sin(kx) ** 2
        + np.sin(ky) ** 2
        + (mass + np.cos(kx) + np.cos(ky)) ** 2
    )
    return -norm


def _qwz_lower_curvature(kx, ky, mass):
    """QWZ lower-band curvature for A=i<u|grad u>."""

    zero = np.zeros_like(kx)
    d = np.stack(
        [np.sin(kx), np.sin(ky), mass + np.cos(kx) + np.cos(ky)], axis=-1
    )
    dkx = np.stack([np.cos(kx), zero, -np.sin(kx)], axis=-1)
    dky = np.stack([zero, np.cos(ky), -np.sin(ky)], axis=-1)
    triple = np.einsum("...i,...i->...", d, np.cross(dkx, dky))
    return 0.5 * triple / np.linalg.norm(d, axis=-1) ** 3


def _qwz_center_mesh(n, mass, omega, band_mode="single-fukui"):
    center = (np.arange(n, dtype=float) + 0.5) / n
    fx, fy = np.meshgrid(center, center, indexing="ij")
    frac = np.column_stack([fx.ravel(), fy.ravel(), np.zeros(n * n)])
    cart = frac[:, 0, None] * B1 + frac[:, 1, None] * B2

    vertex = np.arange(n, dtype=float) * TWO_PI / n
    kx, ky = np.meshgrid(vertex, vertex, indexing="ij")
    energy = _qwz_lower_energy(kx, ky, mass)
    vertices = np.empty((n, n, 4), dtype=float)
    vertices[:, :, 0] = energy
    vertices[:, :, 1] = np.roll(energy, -1, axis=0)
    vertices[:, :, 2] = np.roll(np.roll(energy, -1, axis=0), -1, axis=1)
    vertices[:, :, 3] = np.roll(energy, -1, axis=1)

    area = np.linalg.norm(np.cross(B1, B2))
    return vt.CurvatureData(
        cart=cart,
        frac=frac,
        omega=np.asarray(omega, dtype=float).ravel(),
        energy=vertices[:, :, 0].ravel(),
        vertex_energies=vertices.reshape(n * n, 4),
        metadata={
            "b1": B1,
            "b2": B2,
            "band": 1,
            "band_mode": band_mode,
            "nk": n * n,
            "k_grid": (n, n),
            "dk_area": area / (n * n),
            "min_direct_band_gap_eV": 2.0 * np.min(np.abs(energy)),
        },
    )


def qwz_analytic_data(n, mass=1.0):
    center = (np.arange(n, dtype=float) + 0.5) * TWO_PI / n
    kx, ky = np.meshgrid(center, center, indexing="ij")
    return _qwz_center_mesh(
        n, mass, _qwz_lower_curvature(kx, ky, mass),
        band_mode="analytic-point-quadrature",
    )


def qwz_fukui_data(n, mass=1.0):
    """Construct VASPBERRY-sign plaquette fluxes from independent eigenvectors."""

    coordinate = np.arange(n, dtype=float) * TWO_PI / n
    state = np.empty((n, n, 2), dtype=complex)
    for ix, kx in enumerate(coordinate):
        for iy, ky in enumerate(coordinate):
            dx = math.sin(kx)
            dy = math.sin(ky)
            dz = mass + math.cos(kx) + math.cos(ky)
            hamiltonian = np.array(
                [[dz, dx - 1j * dy], [dx + 1j * dy, -dz]], dtype=complex
            )
            _, eigenvectors = np.linalg.eigh(hamiltonian)
            state[ix, iy] = eigenvectors[:, 0]

    flux = np.empty((n, n), dtype=float)
    for ix in range(n):
        for iy in range(n):
            u00 = state[ix, iy]
            u10 = state[(ix + 1) % n, iy]
            u11 = state[(ix + 1) % n, (iy + 1) % n]
            u01 = state[ix, (iy + 1) % n]
            overlap_loop = (
                np.vdot(u00, u10)
                * np.vdot(u10, u11)
                * np.vdot(u11, u01)
                * np.vdot(u01, u00)
            )
            # This is the sign used by vaspberry.f.
            flux[ix, iy] = -np.angle(overlap_loop)

    dsk = np.linalg.norm(np.cross(B1, B2)) / (n * n)
    return _qwz_center_mesh(n, mass, flux / dsk)


def _periodic_square_distance(frac, center):
    delta = frac[:, :2] - np.asarray(center)
    delta -= np.rint(delta)
    return np.linalg.norm(delta @ np.vstack([B1, B2]), axis=1)


def split_valley_data(n=60, splitting_eV=0.20):
    """Opposite, discretely normalized +/-1/2 valley-Chern packets."""

    center = (np.arange(n, dtype=float) + 0.5) / n
    fx, fy = np.meshgrid(center, center, indexing="ij")
    frac = np.column_stack([fx.ravel(), fy.ravel(), np.zeros(n * n)])
    cart = frac[:, 0, None] * B1 + frac[:, 1, None] * B2
    k_center = (1.0 / 3.0, 1.0 / 3.0)
    kp_center = (2.0 / 3.0, 2.0 / 3.0)
    radius = 0.85
    width = 0.24
    d_k = _periodic_square_distance(frac, k_center)
    d_kp = _periodic_square_distance(frac, kp_center)
    mask_k = (d_k <= radius) & (d_k <= d_kp)
    mask_kp = (d_kp <= radius) & (d_kp < d_k)
    packet_k = np.exp(-0.5 * (d_k / width) ** 2) * mask_k
    packet_kp = np.exp(-0.5 * (d_kp / width) ** 2) * mask_kp

    area = np.linalg.norm(np.cross(B1, B2))
    chern_weight = area / (n * n) / TWO_PI
    omega = (
        0.5 * packet_k / (np.sum(packet_k) * chern_weight)
        - 0.5 * packet_kp / (np.sum(packet_kp) * chern_weight)
    )
    energy = np.ones(n * n)
    energy[mask_k] = 0.0
    energy[mask_kp] = splitting_eV
    data = vt.CurvatureData(
        cart=cart,
        frac=frac,
        omega=omega,
        energy=energy,
        vertex_energies=np.repeat(energy[:, None], 4, axis=1),
        metadata={
            "b1": B1,
            "b2": B2,
            "band": 1,
            "band_mode": "single-fukui",
            "nk": n * n,
            "k_grid": (n, n),
            "dk_area": area / (n * n),
            "chern_reported": 0.0,
            "min_direct_band_gap_eV": 1.0,
        },
    )
    return data, k_center, kp_center, radius


class PhysicsOracleTests(unittest.TestCase):
    def test_qwz_fukui_chern_and_hall_sign(self):
        # Right-handed b1 x b2 and -Im(log loop): m=+1 has C=-1.
        for mass, expected_chern in ((1.0, -1.0), (-1.0, 1.0), (2.5, 0.0)):
            with self.subTest(mass=mass):
                data = qwz_fukui_data(16, mass)
                _, chern = vt.validate_fukui_geometry(data)
                self.assertAlmostEqual(chern, expected_chern, places=12)
                result = vt.integrate_sigma([data], np.array([-6.0, 0.0]), 0.0)
                np.testing.assert_allclose(
                    result["chern_weight_total"], [0.0, expected_chern], atol=1.0e-12
                )
                np.testing.assert_allclose(
                    result["sigma_xy_total_e2_over_h"], [0.0, -expected_chern],
                    atol=1.0e-12,
                )

    def test_qwz_analytic_grid_convergence(self):
        bounds = {8: 1.7e-2, 12: 1.1e-3, 16: 7.0e-5, 24: 3.0e-7, 32: 1.1e-9}
        previous_error = float("inf")
        for n, bound in bounds.items():
            with self.subTest(n=n):
                # This is an independent point-quadrature oracle, not a Fukui
                # plaquette-flux input. Its coarse-grid sum is not exactly
                # quantized and must not weaken the validated transport API's
                # nearest-integer Fukui gate.
                data = qwz_analytic_data(n)
                dsk = vt.reciprocal_area(data) / (n * n)
                chern = np.sum(data.omega * dsk) / TWO_PI
                error = abs(chern + 1.0)
                self.assertLess(error, bound)
                self.assertLess(error, previous_error)
                previous_error = error

    def test_degenerate_opposite_valleys_cancel(self):
        data, k, kp, radius = split_valley_data(splitting_eV=0.0)
        result = vt.integrate_sigma(
            [data], np.array([-0.1, 0.1]), 0.0,
            k_center=k, kp_center=kp, valley_radius=radius,
        )
        np.testing.assert_allclose(result["sigma_xy_total_e2_over_h"], [0.0, 0.0], atol=5e-14)
        np.testing.assert_allclose(
            result["sigma_xy_valley_contrast_e2_over_h"], [0.0, -1.0], atol=5e-14
        )

    def test_valley_splitting_unmasks_then_recancels(self):
        data, k, kp, radius = split_valley_data(splitting_eV=0.20)
        result = vt.integrate_sigma(
            [data], np.array([-0.1, 0.1, 0.3]), 0.0,
            k_center=k, kp_center=kp, valley_radius=radius,
        )
        np.testing.assert_allclose(
            result["sigma_xy_total_e2_over_h"], [0.0, -0.5, 0.0], atol=5e-14
        )
        np.testing.assert_allclose(
            result["sigma_xy_K_e2_over_h"], [0.0, -0.5, -0.5], atol=5e-14
        )
        np.testing.assert_allclose(
            result["sigma_xy_Kp_e2_over_h"], [0.0, 0.0, 0.5], atol=5e-14
        )
        np.testing.assert_allclose(
            result["sigma_xy_valley_contrast_e2_over_h"], [0.0, -0.5, -1.0],
            atol=5e-14,
        )


if __name__ == "__main__":
    unittest.main()
