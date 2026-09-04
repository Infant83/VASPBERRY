import unittest
import json
from contextlib import redirect_stderr
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

import vaspberry_transport as vt


REPO = Path(__file__).resolve().parents[1]


def synthetic_uniform_band(nx=12, ny=10, energy=0.0, chern=1.0, shift=(0.0, 0.0)):
    # VASPBERRY reports Fukui plaquette centers, not WAVECAR vertices.
    x = np.mod((np.arange(nx) + 0.5) / nx + shift[0], 1.0)
    y = np.mod((np.arange(ny) + 0.5) / ny + shift[1], 1.0)
    xx, yy = np.meshgrid(x, y, indexing="ij")
    frac = np.column_stack([xx.ravel(), yy.ravel(), np.zeros(nx * ny)])
    b1 = np.array([2.0 * np.pi, 0.0, 0.0])
    b2 = np.array([0.0, 2.0 * np.pi, 0.0])
    cart = frac[:, 0, None] * b1 + frac[:, 1, None] * b2
    area = np.linalg.norm(np.cross(b1, b2))
    omega = np.full(nx * ny, chern * 2.0 * np.pi / area)
    return vt.CurvatureData(
        cart,
        frac,
        omega,
        metadata={
            "b1": b1,
            "b2": b2,
            "band": 1,
            "band_mode": "single-fukui",
            "nk": nx * ny,
            "dk_area": area / (nx * ny),
            "chern_reported": chern,
            "min_direct_band_gap_eV": 2.0,
        },
        energy=np.full(nx * ny, energy),
        vertex_energies=np.full((nx * ny, 4), energy),
    )


def synthetic_eigenval(nx=12, ny=10, band_energy=0.0, weights=None, shift=(0.0, 0.0)):
    # WAVECAR/EIGENVAL mesh gives the vertices of the centered Fukui plaquettes.
    x = np.mod(np.arange(nx) / nx + shift[0], 1.0)
    y = np.mod(np.arange(ny) / ny + shift[1], 1.0)
    xx, yy = np.meshgrid(x, y, indexing="ij")
    frac = np.column_stack([xx.ravel(), yy.ravel(), np.zeros(nx * ny)])
    nk = nx * ny
    energies = np.column_stack([np.full(nk, band_energy), np.full(nk, band_energy + 2.0)])
    occupations = np.column_stack([np.ones(nk), np.zeros(nk)])
    if weights is None:
        weights = np.full(nk, 1.0 / nk)
    return vt.EigenvalData(
        frac=frac,
        energies=energies,
        occupations=occupations,
        weights=np.asarray(weights),
        nelect=1,
        nk=nk,
        nband=2,
        source=Path("synthetic-EIGENVAL"),
    )


def synthetic_window_eigenval(energies):
    energy_array = np.atleast_2d(np.asarray(energies, dtype=float))
    nk, nband = energy_array.shape
    frac = np.zeros((nk, 3), dtype=float)
    if nk > 1:
        frac[:, 0] = np.arange(nk) / nk
    return vt.EigenvalData(
        frac=frac,
        energies=energy_array,
        occupations=np.zeros_like(energy_array),
        weights=np.full(nk, 1.0 / nk),
        nelect=0,
        nk=nk,
        nband=nband,
        source=Path("synthetic-window-EIGENVAL"),
    )


def nine_reciprocal_copies(data):
    frac_blocks = []
    cart_blocks = []
    omega_blocks = []
    b1 = np.asarray(data.metadata["b1"])
    b2 = np.asarray(data.metadata["b2"])
    for sx in (-1, 0, 1):
        for sy in (-1, 0, 1):
            frac_blocks.append(data.frac + np.array([sx, sy, 0.0]))
            cart_blocks.append(data.cart + sx * b1 + sy * b2)
            omega_blocks.append(data.omega.copy())
    frac = np.vstack(frac_blocks)
    cart = np.vstack(cart_blocks)
    omega = np.concatenate(omega_blocks)
    permutation = np.random.default_rng(7).permutation(len(omega))
    return vt.CurvatureData(
        cart=cart[permutation],
        frac=frac[permutation],
        omega=omega[permutation],
        metadata=dict(data.metadata),
    )


class TransportTests(unittest.TestCase):
    def test_zero_temperature_uses_closed_occupied_side(self):
        mu = 0.5
        energies = np.array([mu - 5.0e-13, mu, mu + 5.0e-13])
        np.testing.assert_array_equal(
            vt.fermi_dirac(energies, mu, temperature=0.0),
            [1.0, 1.0, 0.0],
        )

    def test_energy_column_accepts_only_available_extra_columns(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / "curvature.dat"
            path.write_text("0 0 0 1 0 0 0 2 3\n", encoding="utf-8")
            self.assertEqual(vt.read_vaspberry(path, energy_column=7).energy[0], 2.0)
            self.assertEqual(vt.read_vaspberry(path, energy_column=8).energy[0], 3.0)
            self.assertEqual(vt.read_vaspberry(path).energy[0], 2.0)
            for column in (-1, 6, 9):
                with self.subTest(column=column), self.assertRaisesRegex(
                    ValueError, "available extra.*energy column"
                ):
                    vt.read_vaspberry(path, energy_column=column)

            stderr = StringIO()
            with self.assertRaises(SystemExit), redirect_stderr(stderr):
                vt.main([
                    "line", "--input", str(path), "--energy-column", "9",
                    "--output", str(Path(directory) / "line.png"),
                ])
            self.assertNotIn("Traceback", stderr.getvalue())
            self.assertIn("error:", stderr.getvalue())

    def test_legacy_energy_window_completeness_and_provenance(self):
        eig = synthetic_window_eigenval([
            [0.10, 0.45, 1.00, 1.50],
            [0.15, 0.55, 1.10, 1.60],
        ])
        validation = vt.validate_legacy_sigma_energy_window(
            eig, [2], np.array([0.2, 0.8]), 0.0,
            occupation_tolerance=1.0e-8, core_chern=-1.0,
        )
        self.assertTrue(validation["validated"])
        self.assertEqual(validation["active_window"], [2, 2])
        self.assertEqual(validation["core_bands"], [1])
        self.assertEqual(validation["next_unrepresented_band"], 3)
        self.assertEqual(validation["minimum_core_occupation"], 1.0)
        self.assertEqual(validation["maximum_next_band_occupation"], 0.0)
        self.assertEqual(validation["occupation_tolerance"], 1.0e-8)
        self.assertEqual(validation["zero_temperature_convention"], "E <= mu is occupied")

    def test_legacy_energy_window_rejects_incomplete_band_sets_and_tails(self):
        eig = synthetic_window_eigenval([[0.10, 0.45, 1.00, 1.50, 2.00]])
        with self.assertRaisesRegex(ValueError, "must be consecutive"):
            vt.validate_legacy_sigma_energy_window(eig, [2, 4], np.array([0.2]), 0.0)
        with self.assertRaisesRegex(ValueError, "core band 1 is not fully occupied"):
            vt.validate_legacy_sigma_energy_window(eig, [2], np.array([0.05, 0.8]), 0.0)
        with self.assertRaisesRegex(ValueError, "unrepresented band 3 is occupied"):
            vt.validate_legacy_sigma_energy_window(eig, [2], np.array([0.2, 1.0]), 0.0)
        with self.assertRaisesRegex(ValueError, "final band"):
            vt.validate_legacy_sigma_energy_window(eig, [5], np.array([1.9]), 0.0)

        finite_temperature = synthetic_window_eigenval([[-1.0, 0.0, 0.2]])
        with self.assertRaisesRegex(ValueError, "unrepresented band 3 is occupied"):
            vt.validate_legacy_sigma_energy_window(
                finite_temperature, [2], np.array([0.0]), 300.0,
                occupation_tolerance=1.0e-8,
            )
        accepted = vt.validate_legacy_sigma_energy_window(
            finite_temperature, [2], np.array([0.0]), 300.0,
            occupation_tolerance=1.0e-3,
        )
        self.assertLessEqual(accepted["maximum_next_band_occupation"], 1.0e-3)

    def test_core_chern_must_match_omitted_integer_manifold(self):
        eig = synthetic_window_eigenval([[-1.0, 0.0, 1.0]])
        with self.assertRaisesRegex(ValueError, "not within 0.005 of an integer"):
            vt.validate_legacy_sigma_energy_window(
                eig, [2], np.array([0.0]), 0.0, core_chern=0.01
            )
        with self.assertRaisesRegex(ValueError, "starts at band 1"):
            vt.validate_legacy_sigma_energy_window(
                eig, [1], np.array([-0.5]), 0.0, core_chern=1.0
            )

    def test_sigma_cli_records_energy_window_validation_in_summary(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            nx = ny = 4
            band = synthetic_uniform_band(nx=nx, ny=ny, chern=1.0)
            band.metadata["band"] = 2
            curvature_path = root / "BERRYCURV.dat"
            rows = np.column_stack([band.cart, band.omega, band.frac])
            header = (
                f"# NKPOINT : {nx * ny}\n"
                "# NBANDS : 3\n"
                f"# K-GRID : {nx} X {ny}\n"
                "# RECIVEC B1 (A^-1) : 6.283185307179586 0 0\n"
                "# RECIVEC B2 (A^-1) : 0 6.283185307179586 0\n"
                "# Chern Number for the BAND : 2\n"
            )
            curvature_path.write_text(
                header + "\n".join(" ".join(f"{value:.16g}" for value in row) for row in rows) + "\n",
                encoding="utf-8",
            )

            eig = synthetic_eigenval(nx=nx, ny=ny, band_energy=0.0)
            eig.energies = np.column_stack([
                np.full(eig.nk, -1.0), np.zeros(eig.nk), np.full(eig.nk, 2.0)
            ])
            eig.occupations = np.column_stack([
                np.ones(eig.nk), np.zeros(eig.nk), np.zeros(eig.nk)
            ])
            eig.nband = 3
            eigenval_path = root / "EIGENVAL"
            eigenval_lines = ["header"] * 5 + [f"1 {eig.nk} {eig.nband}"]
            for kpoint, weight, energies, occupations in zip(
                eig.frac, eig.weights, eig.energies, eig.occupations, strict=True
            ):
                eigenval_lines.append("")
                eigenval_lines.append(
                    " ".join(str(value) for value in (*kpoint, weight))
                )
                eigenval_lines.extend(
                    f"{index} {energy} {occupation}"
                    for index, (energy, occupation) in enumerate(
                        zip(energies, occupations, strict=True), start=1
                    )
                )
            eigenval_path.write_text("\n".join(eigenval_lines) + "\n", encoding="utf-8")

            summary_path = root / "summary.json"
            exit_code = vt.main([
                "sigma", "--input", str(curvature_path), "--bands", "2",
                "--eigenval", str(eigenval_path), "--mu-min", "-0.5", "--mu-max", "1.0",
                "--mu-points", "3", "--occupation-tolerance", "1e-7",
                "--core-chern", "0.0",
                "--csv", str(root / "sigma.csv"), "--plot", str(root / "sigma.png"),
                "--summary", str(summary_path),
            ])
            self.assertEqual(exit_code, 0)
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            self.assertEqual(summary["occupation_tolerance"], 1.0e-7)
            self.assertTrue(summary["energy_window_validation"]["validated"])
            self.assertEqual(summary["energy_window_validation"]["active_window"], [2, 2])
            self.assertEqual(summary["zero_temperature_occupation"], "E <= mu is occupied")
            self.assertEqual(summary["core_chern_provenance"], "explicit CLI input")
            self.assertEqual(summary["inputs"], ["BERRYCURV.dat"])
            self.assertEqual(summary["eigenval"], "EIGENVAL")

            stderr = StringIO()
            with self.assertRaises(SystemExit), redirect_stderr(stderr):
                vt.main([
                    "sigma", "--input", str(curvature_path), "--bands", "2",
                    "--eigenval", str(eigenval_path), "--mu-min", "-0.5",
                    "--mu-max", "1.0", "--mu-points", "3",
                    "--csv", str(root / "missing-core.csv"),
                    "--plot", str(root / "missing-core.png"),
                ])
            self.assertIn("--core-chern must be supplied explicitly", stderr.getvalue())

    def test_sigma_mesh_requires_two_dimensions(self):
        for shape in ((1, 1), (1, 4), (4, 1)):
            with self.subTest(shape=shape), self.assertRaisesRegex(ValueError, "genuinely two-dimensional"):
                vt.uniform_full_bz_shape(synthetic_uniform_band(*shape))
        self.assertEqual(vt.uniform_full_bz_shape(synthetic_uniform_band(2, 2)), (2, 2))

    def test_constant_curvature_recovers_integer_chern_and_hall_sign(self):
        data = synthetic_uniform_band(chern=1.0)
        result = vt.integrate_sigma([data], np.array([-1.0, 1.0]), temperature=0.0)
        np.testing.assert_allclose(result["chern_weight_total"], [0.0, 1.0], atol=1e-13)
        np.testing.assert_allclose(result["sigma_xy_total_e2_over_h"], [0.0, -1.0], atol=1e-13)

    def test_equal_opposite_bands_cancel(self):
        plus = synthetic_uniform_band(energy=0.0, chern=1.0)
        minus = synthetic_uniform_band(energy=0.0, chern=-1.0)
        result = vt.integrate_sigma([plus, minus], np.array([1.0]), temperature=0.0)
        self.assertAlmostEqual(result["sigma_xy_total_e2_over_h"][0], 0.0, places=13)

    def test_nine_copy_deduplication_preserves_transport(self):
        unique = synthetic_uniform_band(nx=6, ny=4, energy=0.0, chern=1.0)
        repeated = nine_reciprocal_copies(unique)
        collapsed = vt.deduplicate_reciprocal_replicas(repeated)
        self.assertEqual(len(collapsed.omega), 24)
        self.assertEqual(collapsed.metadata["reciprocal_replica_count"], 9)
        self.assertEqual(collapsed.metadata["omega_replica_max_spread"], 0.0)
        vt.attach_fukui_vertex_energies(collapsed, synthetic_eigenval(6, 4), 1)
        reference = vt.integrate_sigma([unique], np.array([1.0]), 0.0)
        actual = vt.integrate_sigma([collapsed], np.array([1.0]), 0.0)
        np.testing.assert_allclose(
            actual["sigma_xy_total_e2_over_h"], reference["sigma_xy_total_e2_over_h"], atol=1e-13
        )

    def test_inconsistent_reciprocal_copy_is_rejected(self):
        repeated = nine_reciprocal_copies(synthetic_uniform_band(nx=4, ny=4))
        repeated.omega[0] += 0.1
        with self.assertRaisesRegex(ValueError, "differs among reciprocal replicas"):
            vt.deduplicate_reciprocal_replicas(repeated)

    def test_fukui_uses_mean_of_four_vertex_occupations(self):
        data = synthetic_uniform_band(nx=2, ny=2, energy=1.0, chern=0.0)
        dsk = vt.reciprocal_area(data) / 4.0
        data.omega[:] = 0.0
        data.omega[0] = 0.5 * np.pi / dsk  # safely below the principal-log branch
        data.omega[1] = -data.omega[0]  # keep the fully occupied Chern integer (zero)
        data.metadata.pop("chern_reported")
        data.vertex_energies[:] = 1.0
        data.vertex_energies[0] = [-1.0, -1.0, 1.0, 1.0]
        result = vt.integrate_sigma([data], np.array([0.0]), temperature=0.0)
        self.assertAlmostEqual(result["chern_weight_total"][0], 0.125, places=13)
        self.assertAlmostEqual(result["sigma_xy_total_e2_over_h"][0], -0.125, places=13)

    def test_finite_temperature_uses_mean_f_not_f_of_mean_energy(self):
        data = synthetic_uniform_band(nx=2, ny=2, chern=0.0)
        dsk = vt.reciprocal_area(data) / 4.0
        data.omega[:] = 0.0
        data.omega[0] = 0.5 * np.pi / dsk
        data.omega[1] = -data.omega[0]
        data.metadata.pop("chern_reported")
        energies = np.array([-0.1, 0.3, 0.3, 0.3])
        data.vertex_energies[:] = 1.0
        data.vertex_energies[0] = energies
        temperature = 300.0
        result = vt.integrate_sigma([data], np.array([0.0]), temperature)
        expected = 0.25 * np.mean(vt.fermi_dirac(energies, 0.0, temperature))
        wrong = 0.25 * vt.fermi_dirac(np.array([np.mean(energies)]), 0.0, temperature)[0]
        self.assertAlmostEqual(result["chern_weight_total"][0], expected, places=13)
        self.assertGreater(abs(expected - wrong), 1.0e-3)

    def test_full_eigenval_mesh_attaches_centered_plaquette_vertices(self):
        data = synthetic_uniform_band(nx=6, ny=4)
        data.vertex_energies = None
        eig = synthetic_eigenval(6, 4, band_energy=-0.3)
        vt.attach_fukui_vertex_energies(data, eig, 1)
        self.assertEqual(data.vertex_energies.shape, (24, 4))
        np.testing.assert_allclose(data.vertex_energies, -0.3)
        self.assertTrue(data.metadata["eigenval_full_mesh_validated"])

    def test_attach_rejects_one_dimensional_fukui_mesh(self):
        data = synthetic_uniform_band(nx=1, ny=4)
        data.vertex_energies = None
        with self.assertRaisesRegex(ValueError, "Nx,Ny >= 2"):
            vt.attach_fukui_vertex_energies(data, synthetic_eigenval(1, 4), 1)

    def test_curvature_and_eigenval_kz_planes_must_match_periodically(self):
        data = synthetic_uniform_band(nx=4, ny=4)
        data.vertex_energies = None
        data.frac[:, 2] = 0.25
        mismatch = synthetic_eigenval(4, 4)
        with self.assertRaisesRegex(ValueError, "kz planes do not match"):
            vt.attach_fukui_vertex_energies(data, mismatch, 1)
        periodic = synthetic_eigenval(4, 4)
        periodic.frac[:, 2] = 1.25
        vt.attach_fukui_vertex_energies(data, periodic, 1)

    def test_dedup_rejects_multiple_kz_planes_but_accepts_periodic_images(self):
        base = synthetic_uniform_band(nx=4, ny=4)
        multiple = vt.CurvatureData(
            cart=np.vstack([base.cart, base.cart]),
            frac=np.vstack([base.frac, base.frac + np.array([0.0, 0.0, 0.25])]),
            omega=np.r_[base.omega, base.omega], metadata=dict(base.metadata),
        )
        with self.assertRaisesRegex(ValueError, "multiple inequivalent kz planes"):
            vt.deduplicate_reciprocal_replicas(multiple)
        periodic = vt.CurvatureData(
            cart=np.vstack([base.cart, base.cart]),
            frac=np.vstack([base.frac, base.frac + np.array([0.0, 0.0, 1.0])]),
            omega=np.r_[base.omega, base.omega], metadata=dict(base.metadata),
        )
        collapsed = vt.deduplicate_reciprocal_replicas(periodic)
        np.testing.assert_allclose(collapsed.frac[:, 2], 0.0)

    def test_shifted_mesh_and_shuffled_eigenval_order_are_accepted(self):
        shift = (0.03125, 0.0625)
        data = synthetic_uniform_band(nx=6, ny=4, shift=shift)
        data.vertex_energies = None
        eig = synthetic_eigenval(6, 4, band_energy=-0.2, shift=shift)
        order = np.random.default_rng(19).permutation(eig.nk)
        eig.frac = eig.frac[order]
        eig.energies = eig.energies[order]
        eig.occupations = eig.occupations[order]
        eig.weights = eig.weights[order]
        vt.attach_fukui_vertex_energies(data, eig, 1)
        np.testing.assert_allclose(data.vertex_energies, -0.2)

    def test_vertex_mapping_preserves_tuple_order_across_periodic_boundaries(self):
        nx, ny = 4, 3
        shift = (0.071, 0.113)
        data = synthetic_uniform_band(nx=nx, ny=ny, shift=shift)
        data.vertex_energies = None
        eig = synthetic_eigenval(nx=nx, ny=ny, shift=shift)
        ix, iy = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
        labels = (100.0 * ix + iy).ravel()
        eig.energies[:, 0] = labels
        eig.energies[:, 1] = labels + 1000.0

        order = np.random.default_rng(23).permutation(eig.nk)
        eig.frac = eig.frac[order]
        eig.energies = eig.energies[order]
        eig.occupations = eig.occupations[order]
        eig.weights = eig.weights[order]
        vt.attach_fukui_vertex_energies(data, eig, 1)

        # Vertex order is (00, 10, 11, 01). Check one interior cell and the
        # corner cell that wraps across both reciprocal boundaries.
        np.testing.assert_array_equal(
            data.vertex_energies[1 * ny + 1], [101.0, 201.0, 202.0, 102.0]
        )
        np.testing.assert_array_equal(
            data.vertex_energies[3 * ny + 2], [302.0, 2.0, 0.0, 300.0]
        )

    def test_duplicate_or_missing_eigenval_point_is_rejected(self):
        data = synthetic_uniform_band(nx=4, ny=4)
        eig = synthetic_eigenval(4, 4)
        eig.frac[-1] = eig.frac[0]
        with self.assertRaisesRegex(ValueError, "duplicate points"):
            vt.attach_fukui_vertex_energies(data, eig, 1)

    def test_requested_band_and_total_band_count_must_match_headers(self):
        data = synthetic_uniform_band(nx=4, ny=4)
        data.vertex_energies = None
        data.metadata["band"] = 2
        with self.assertRaisesRegex(ValueError, "requested band 1 conflicts"):
            vt.attach_fukui_vertex_energies(data, synthetic_eigenval(4, 4), 1)
        data.metadata["band"] = 1
        data.metadata["nbands"] = 3
        with self.assertRaisesRegex(ValueError, "header NBANDS=3"):
            vt.attach_fukui_vertex_energies(data, synthetic_eigenval(4, 4), 1)
        vt.validate_unique_active_bands([1, 2, 3])
        with self.assertRaisesRegex(ValueError, "duplicate active band inputs"):
            vt.validate_unique_active_bands([1, 2, 1])

    def test_up_dn_filename_parsing_and_global_spin_validation(self):
        content = """# ISPIN : 2 (LSORBIT = .FALSE.)
# Chern Number for the BAND : 1
0 0 0 0 0 0 0
"""
        with TemporaryDirectory() as directory:
            up_path = Path(directory) / "BERRYCURV.UP.dat"
            dn_path = Path(directory) / "BERRYCURV.DN.dat"
            plain_path = Path(directory) / "BERRYCURV.dat"
            up_path.write_text(content, encoding="utf-8")
            dn_path.write_text(content, encoding="utf-8")
            plain_path.write_text(content, encoding="utf-8")
            up = vt.read_vaspberry(up_path)
            dn = vt.read_vaspberry(dn_path)
            plain = vt.read_vaspberry(plain_path)
            self.assertEqual(up.metadata["spin_channel"], 1)
            self.assertEqual(dn.metadata["spin_channel"], 2)
            vt.validate_spin_channels([up], 1)
            vt.validate_spin_channels([dn], 2)
            # A curvature-only plot may show a renamed collinear file, but an
            # energy-resolved calculation must know which channel it contains.
            vt.validate_spin_channels([plain], 1, require_collinear_label=False)
            with self.assertRaisesRegex(ValueError, "channel is ambiguous"):
                vt.validate_spin_channels([plain], 1)
            for datasets, spin, pattern in (
                ([up], 2, "global --spin"), ([dn], 1, "global --spin"),
                ([up, dn], 1, "cannot be mixed"),
                ([up, plain], 1, "labeled and unlabeled"),
            ):
                with self.subTest(spin=spin, pattern=pattern), self.assertRaisesRegex(ValueError, pattern):
                    vt.validate_spin_channels(datasets, spin)

    def test_conflicting_or_invalid_spin_metadata_is_rejected(self):
        with TemporaryDirectory() as directory:
            conflict = Path(directory) / "BERRYCURV.UP.dat"
            conflict.write_text(
                "# ISPIN: 2 (LSORBIT=.FALSE.)\n# SPIN: DN\n0 0 0 0 0 0 0\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "conflicts with filename"):
                vt.read_vaspberry(conflict)
            invalid = Path(directory) / "BERRYCURV.DN.dat"
            invalid.write_text("# ISPIN: 1 (LSORBIT=.FALSE.)\n0 0 0 0 0 0 0\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "requires ISPIN=2"):
                vt.read_vaspberry(invalid)

            scalar = Path(directory) / "BERRYCURV.scalar.dat"
            scalar.write_text("# ISPIN: 1 (LSORBIT=.TRUE.)\n0 0 0 0 0 0 0\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "ISPIN=1 curvature"):
                vt.validate_spin_channels([vt.read_vaspberry(scalar)], 2)

    def test_reduced_or_nonuniform_eigenval_weights_are_rejected(self):
        data = synthetic_uniform_band(nx=4, ny=4)
        weights = np.full(16, 1.0 / 16)
        weights[0] += 0.01
        weights[1] -= 0.01
        with self.assertRaisesRegex(ValueError, "nonuniform k weights"):
            vt.attach_fukui_vertex_energies(data, synthetic_eigenval(4, 4, weights=weights), 1)

    def test_legacy_periodic_map_collapses_and_validates_chern(self):
        repeated = vt.read_vaspberry(
            REPO / "examples/1H-MoS2/BERRYCURV.dat"
        )
        self.assertEqual(len(repeated.omega), 9 * 144)
        fundamental = vt.deduplicate_reciprocal_replicas(repeated)
        self.assertEqual(len(fundamental.omega), 144)
        self.assertEqual(vt.uniform_full_bz_shape(fundamental), (12, 12))
        _, chern = vt.validate_fukui_geometry(fundamental)
        self.assertAlmostEqual(chern, 0.0, places=3)

    def test_map_and_line_sampling_modes_are_not_interchangeable(self):
        raw_map = vt.read_vaspberry(
            REPO / "examples/1H-MoS2/BERRYCURV.dat"
        )
        kubo_path = vt.read_vaspberry(
            REPO
            / "examples/1H-MoS2/KPATH/3.BC_kubo/BERRYCURV_KUBO.EIG-10.dat"
        )
        with self.assertRaisesRegex(ValueError, "2D full-BZ map"):
            vt.prepare_k_path(raw_map)
        with self.assertRaisesRegex(ValueError, "map/cut input must be"):
            vt.prepare_full_bz_map(kubo_path)
        unique_map = vt.prepare_full_bz_map(raw_map)
        self.assertEqual(len(unique_map.omega), 144)

    def test_corrupt_extended_map_is_not_reinterpreted_as_line(self):
        raw_map = vt.read_vaspberry(
            REPO / "examples/1H-MoS2/BERRYCURV.dat"
        )
        raw_map.omega[0] += 0.2
        with self.assertRaisesRegex(ValueError, "differs among reciprocal replicas"):
            vt.prepare_k_path(raw_map)

    def test_kubo_path_is_plot_only_and_not_a_full_bz(self):
        curvature = vt.read_vaspberry(
            REPO
            / "examples/1H-MoS2/KPATH/3.BC_kubo/BERRYCURV_KUBO.EIG-10.dat"
        )
        eig = vt.read_eigenval(
            REPO / "examples/1H-MoS2/KPATH/2.band/EIGENVAL"
        )
        vt.attach_band_energy(curvature, eig, 10)
        self.assertEqual(curvature.metadata["min_direct_band_gap_eV"], 0.0)
        fundamental = vt.select_fundamental_path(curvature)
        self.assertEqual(len(fundamental.energy), 48)
        with self.assertRaisesRegex(ValueError, "duplicate points|full-BZ"):
            vt.integrate_sigma([curvature], np.array([0.0]), temperature=0.0)

    def test_full_grid_kubo_is_explicitly_rejected_for_transport(self):
        data = synthetic_uniform_band(nx=4, ny=4)
        data.metadata["band_mode"] = "single-kubo"
        with self.assertRaisesRegex(ValueError, "plotting only"):
            vt.integrate_sigma([data], np.array([1.0]), temperature=0.0)

    def test_manifold_curvature_is_rejected_for_single_energy_weight(self):
        data = synthetic_uniform_band()
        data.metadata["band_mode"] = "manifold"
        with self.assertRaisesRegex(ValueError, "band-manifold"):
            vt.integrate_sigma([data], np.array([1.0]), temperature=0.0)

    def test_headerless_mode_is_not_accepted_as_validated_fukui_transport(self):
        data = synthetic_uniform_band(nx=4, ny=4)
        data.metadata.pop("band_mode")
        with self.assertRaisesRegex(ValueError, "unsupported curvature mode"):
            vt.integrate_sigma([data], np.array([1.0]), temperature=0.0)
        data.vertex_energies = None
        with self.assertRaisesRegex(ValueError, "single-band Fukui"):
            vt.attach_fukui_vertex_energies(data, synthetic_eigenval(4, 4), 1)

    def test_valley_partition_sums_to_total(self):
        data = synthetic_uniform_band(nx=24, ny=20, energy=0.0, chern=1.0)
        result = vt.integrate_sigma(
            [data], np.array([1.0]), 0.0,
            k_center=(0.25, 0.5), kp_center=(0.75, 0.5), valley_radius=0.45,
        )
        partition = (
            result["sigma_xy_K_e2_over_h"]
            + result["sigma_xy_Kp_e2_over_h"]
            + result["sigma_xy_outside_e2_over_h"]
        )
        np.testing.assert_allclose(partition, result["sigma_xy_total_e2_over_h"], atol=1e-13)

    def test_nonpositive_valley_radius_is_rejected(self):
        data = synthetic_uniform_band(nx=4, ny=4)
        with self.assertRaisesRegex(ValueError, "must be positive"):
            vt.integrate_sigma(
                [data], np.array([1.0]), 0.0,
                k_center=(0.25, 0.5), kp_center=(0.75, 0.5), valley_radius=0.0,
            )

    def test_oblique_shortest_cut_uses_cartesian_minimum_over_three_by_three_images(self):
        basis = np.array([[1.0, 0.0, 0.0], [0.8, 0.6, 0.0]])
        start = np.array([0.0, 0.0])
        end = np.array([-0.9, -0.55])
        delta = vt.shortest_reciprocal_delta(start, end, basis)
        np.testing.assert_allclose(delta, [0.1, -0.55])
        component_wrapped = vt._wrap_fractional(np.r_[end - start, 0.0])[:2]
        self.assertLess(
            np.linalg.norm(delta @ basis), np.linalg.norm(component_wrapped @ basis)
        )

        # A deliberately non-reduced but nonsingular basis needs a shift of
        # two reciprocal vectors, beyond the surrounding 3x3 copies.
        nonreduced = np.array([[1.0, 0.0, 0.0], [1.9, 0.1, 0.0]])
        extended = vt.shortest_reciprocal_delta(
            [0.0, 0.0], [-0.99, -0.94], nonreduced
        )
        np.testing.assert_allclose(extended, [-1.99, 1.06])

    def test_periodic_cut_interpolation_is_finite_for_highly_skew_basis(self):
        data = synthetic_uniform_band(nx=8, ny=6)
        basis = np.array([[1.0, 0.0, 0.0], [1.9, 0.1, 0.0]])
        data.metadata["b1"] = basis[0]
        data.metadata["b2"] = basis[1]
        data.cart = data.frac[:, :2] @ basis
        data.omega = np.sin(2.0 * np.pi * data.frac[:, 0]) + np.cos(
            2.0 * np.pi * data.frac[:, 1]
        )
        data.energy = 0.2 * data.omega

        start = np.array([0.0, 0.0])
        end = np.array([-0.99, -0.94])
        delta = vt.shortest_reciprocal_delta(start, end, basis)
        self.assertGreater(np.max(np.abs(delta)), 1.0)
        samples = start + np.linspace(0.0, 1.0, 401)[:, None] * delta
        omega_cut = vt.periodic_bilinear_interpolate(data, data.omega, samples)
        energy_cut = vt.periodic_bilinear_interpolate(data, data.energy, samples)
        self.assertTrue(np.all(np.isfinite(omega_cut)))
        self.assertTrue(np.all(np.isfinite(energy_cut)))

        translated_mesh = data.frac[:, :2] + np.array([10.0, -7.0])
        np.testing.assert_allclose(
            vt.periodic_bilinear_interpolate(data, data.omega, translated_mesh),
            data.omega,
            atol=2.0e-14,
        )
        with TemporaryDirectory() as directory:
            output = Path(directory) / "skew-cut.png"
            vt.plot_valley_cut(data, start, end, output, samples=401)
            self.assertTrue(output.is_file())
            self.assertGreater(output.stat().st_size, 0)

    def test_nonfinite_inputs_are_rejected_early_and_after_mutation(self):
        base = synthetic_uniform_band(nx=4, ny=4)
        constructor_cases = {
            "cart": (base.cart.copy(), base.frac.copy(), base.omega.copy()),
            "frac": (base.cart.copy(), base.frac.copy(), base.omega.copy()),
            "omega": (base.cart.copy(), base.frac.copy(), base.omega.copy()),
        }
        constructor_cases["cart"][0][0, 0] = np.nan
        constructor_cases["frac"][1][0, 0] = np.inf
        constructor_cases["omega"][2][0] = -np.inf
        for label, values in constructor_cases.items():
            with self.subTest(label=label), self.assertRaisesRegex(ValueError, "NaN or infinity"):
                vt.CurvatureData(cart=values[0], frac=values[1], omega=values[2])
        mutated = synthetic_uniform_band(nx=4, ny=4)
        mutated.vertex_energies[0, 0] = np.nan
        with self.assertRaisesRegex(ValueError, "NaN or infinity"):
            vt.integrate_sigma([mutated], np.array([0.0]), 0.0)
        with self.assertRaisesRegex(ValueError, "Fermi-Dirac energies"):
            vt.fermi_dirac(np.array([np.inf]), 0.0, 10.0)
        with self.assertRaisesRegex(ValueError, "chemical potential"):
            vt.fermi_dirac(np.array([0.0]), np.nan, 10.0)
        with self.assertRaisesRegex(ValueError, "temperature"):
            vt.fermi_dirac(np.array([0.0]), 0.0, np.inf)

    def test_isolation_provenance_integer_chern_and_finite_options_are_gated(self):
        missing = synthetic_uniform_band(nx=4, ny=4)
        missing.metadata.pop("min_direct_band_gap_eV")
        with self.assertRaisesRegex(ValueError, "missing min_direct_band_gap_eV"):
            vt.integrate_sigma([missing], np.array([0.0]), 0.0)
        noninteger = synthetic_uniform_band(nx=4, ny=4, chern=0.25)
        with self.assertRaisesRegex(ValueError, "not within 0.005 of an integer"):
            vt.integrate_sigma([noninteger], np.array([0.0]), 0.0)
        valid = synthetic_uniform_band(nx=4, ny=4)
        with self.assertRaisesRegex(ValueError, "isolation tolerance"):
            vt.integrate_sigma([valid], np.array([0.0]), 0.0, isolation_tolerance=-1.0)
        with self.assertRaisesRegex(ValueError, "core Chern"):
            vt.integrate_sigma([valid], np.array([0.0]), 0.0, core_chern=np.nan)

    def test_failed_sigma_run_removes_only_requested_stale_outputs(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            csv_path = root / "sigma.csv"
            plot_path = root / "sigma.png"
            summary_path = root / "sigma.json"
            unrelated = root / "keep.txt"
            for path in (csv_path, plot_path, summary_path, unrelated):
                path.write_text("stale\n", encoding="utf-8")

            with self.assertRaises(SystemExit), redirect_stderr(StringIO()):
                vt.main([
                    "sigma",
                    "--input", str(root / "missing-curvature.dat"),
                    "--eigenval", str(root / "missing-EIGENVAL"),
                    "--mu-min", "0", "--mu-max", "1",
                    "--isolation-tolerance", "-1",
                    "--csv", str(csv_path),
                    "--plot", str(plot_path),
                    "--summary", str(summary_path),
                ])

            for path in (csv_path, plot_path, summary_path):
                self.assertFalse(path.exists())
            self.assertEqual(unrelated.read_text(encoding="utf-8"), "stale\n")

    def test_sigma_path_collisions_are_rejected_before_cleanup(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            curvature = root / "curvature.dat"
            eigenval = root / "EIGENVAL"
            plot_path = root / "sigma.png"
            summary_path = root / "sigma.json"
            for path in (curvature, eigenval, plot_path, summary_path):
                path.write_text(f"preserve {path.name}\n", encoding="utf-8")

            common = [
                "sigma", "--input", str(curvature), "--eigenval", str(eigenval),
                "--mu-min", "0", "--mu-max", "1",
                "--isolation-tolerance", "-1",
            ]
            with self.assertRaises(SystemExit), redirect_stderr(StringIO()):
                vt.main(common + [
                    "--csv", str(curvature),
                    "--plot", str(plot_path),
                    "--summary", str(summary_path),
                ])
            self.assertEqual(
                curvature.read_text(encoding="utf-8"), "preserve curvature.dat\n"
            )
            self.assertTrue(plot_path.exists())
            self.assertTrue(summary_path.exists())

            shared_output = root / "shared.out"
            shared_output.write_text("preserve shared\n", encoding="utf-8")
            with self.assertRaises(SystemExit), redirect_stderr(StringIO()):
                vt.main(common + [
                    "--csv", str(shared_output),
                    "--plot", str(shared_output),
                    "--summary", str(summary_path),
                ])
            self.assertEqual(shared_output.read_text(encoding="utf-8"), "preserve shared\n")
            self.assertTrue(summary_path.exists())

    def test_plot_subcommands_refuse_to_overwrite_their_input(self):
        with TemporaryDirectory() as directory:
            source = Path(directory) / "curvature.dat"
            original = "input must survive\n"
            commands = {
                "map": [],
                "line": [],
                "cut": [
                    "--k-center", "0.3333333333", "0.3333333333",
                    "--kp-center", "0.6666666667", "0.6666666667",
                ],
            }
            for command, extra in commands.items():
                with self.subTest(command=command):
                    source.write_text(original, encoding="utf-8")
                    stderr = StringIO()
                    with self.assertRaises(SystemExit), redirect_stderr(stderr):
                        vt.main([
                            command,
                            "--input", str(source),
                            "--output", str(source),
                            *extra,
                        ])
                    self.assertIn("refusing to overwrite plot input", stderr.getvalue())
                    self.assertEqual(source.read_text(encoding="utf-8"), original)

    def test_core_chern_adds_constant_total_baseline_only(self):
        data = synthetic_uniform_band(energy=2.0, chern=1.0)
        data.metadata["band"] = 2
        result = vt.integrate_sigma([data], np.array([0.0]), 0.0, core_chern=-1.0)
        self.assertAlmostEqual(result["chern_weight_total"][0], -1.0)
        self.assertAlmostEqual(result["sigma_xy_total_e2_over_h"][0], 1.0)


if __name__ == "__main__":
    unittest.main()
