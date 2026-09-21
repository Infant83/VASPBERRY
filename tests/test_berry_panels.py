"""Independent numerical and data-integrity checks for map/path/band panels."""
from pathlib import Path
import csv
import hashlib
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import plot_berry_panels as panels
from plot_berry_curvature import polygon_area


def grid(n=4, centered=True):
    axis = (np.arange(n) + (.5 if centered else 0.)) / n
    return np.array([[x, y, 0.] for x in axis for y in axis])


def band_rows(q, reciprocal):
    distance = np.r_[0., np.cumsum(np.linalg.norm(np.diff(q @ reciprocal, axis=0), axis=1))]
    return [dict(k_index=i + 1, q1=p[0], q2=p[1], q3=p[2], path_inv_A=s,
                 band=b, energy_eV=e)
            for i, (p, s) in enumerate(zip(q, distance))
            for b, e in [(1, -1. + p[0] ** 2), (2, 1. + p[0] ** 2)]]


def write_bands(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_kubo(path, q, values, gaps=None, normalization="STANDARD_MINUS_TWO_IM", intermediate=2):
    if gaps is None:
        gaps = np.full(len(q), 2.)
    with path.open("w", newline="") as handle:
        handle.write("# schema=VASPBERRY_BARE_MOMENTUM_KUBO_V2\n")
        handle.write(f"# normalization={normalization}\n")
        handle.write("# operator=WAVECAR_BARE_MOMENTUM_NO_PAW_NONLOCAL_VELOCITY\n")
        handle.write("# berry_connection=A_i=i<u|d/dk_i u>\n")
        handle.write(f"# intermediate_bands=1:{intermediate}\n")
        writer = csv.writer(handle)
        writer.writerow(["spin", "k_index", "band", "kx_frac", "ky_frac", "kz_frac", "energy_eV", "omega_z_A2", "min_gap_eV"])
        for i, (p, v, gap) in enumerate(zip(q, values, gaps)):
            writer.writerow([1, i + 1, 1, *p, -1. + p[0] ** 2, v, gap])


def write_bundle(path, q, values, gaps=None, metadata=None, spins=None):
    if gaps is None:
        gaps = np.full(len(q), 2.)
    if spins is None:
        spins = np.ones(len(q), dtype=int)
    headers = {
        "schema": "VASPBERRY_BARE_MOMENTUM_KUBO_BUNDLE_V1",
        "normalization": "STANDARD_MINUS_TWO_IM",
        "operator": "WAVECAR_BARE_MOMENTUM_NO_PAW_NONLOCAL_VELOCITY",
        "berry_connection": "A_i=i<u|d/dk_i u>",
        "intermediate_bands": "EXTERNAL_TO_SELECTED_BUNDLE_WITHIN_SOURCE_NBANDS",
        "result_status": "PASS", "result_kind": "ISOLATED_BUNDLE_TRACE",
        "curvature": "TRACE_OF_SELECTED_BUNDLE", "internal_transitions": "EXCLUDED_ANALYTICALLY",
        "occupation_weighting": "NONE", "band_min": "1", "band_max": "1", "band_rank": "1",
        "source_nbands": "2", "no_external_states": "false", "gap_threshold_eV": "1e-5",
    }
    headers.update(metadata or {})
    with path.open("w", newline="") as handle:
        for key, value in headers.items():
            handle.write(f"# {key}={value}\n")
        writer = csv.writer(handle)
        writer.writerow(["spin", "k_index", "kx_frac", "ky_frac", "kz_frac", "omega_z_A2", "min_external_gap_eV"])
        for i, (p, v, gap, spin) in enumerate(zip(q, values, gaps, spins)):
            writer.writerow([spin, i + 1, *p, v, gap])


class PeriodicDisplayCutTests(unittest.TestCase):
    def test_constant_field_is_unchanged_far_outside_primitive_cell(self):
        q = grid()
        target = np.array([[.04, -.31, 0.], [18.55, -21.3, 0.], [-7., 14., 0.]])
        np.testing.assert_allclose(panels.periodic_bilinear(q, np.full(len(q), 7.25), target), 7.25,
                                   atol=1e-13, rtol=0)

    def test_linear_field_is_exact_inside_a_cell(self):
        q = grid()
        values = 4 + 2 * q[:, 0] + 3 * q[:, 1]
        target = np.array([[.45, .60, 0.], [.51, .39, 0.], [.38, .49, 0.]])
        expected = 4 + 2 * target[:, 0] + 3 * target[:, 1]
        np.testing.assert_allclose(panels.periodic_bilinear(q, values, target), expected, atol=1e-13, rtol=0)

    def test_periodic_sinusoidal_grid_wraps_at_both_boundaries(self):
        q = grid()
        values = np.cos(2 * np.pi * q[:, 0]) + .5 * np.sin(2 * np.pi * q[:, 1])
        # Four surrounding samples at (±1/8, ±1/8): sine cancels and cosines equal sqrt(1/2).
        target = np.array([[0., 0., 0.], [1., -1., 0.], [-3., 7., 0.]])
        np.testing.assert_allclose(panels.periodic_bilinear(q, values, target), np.sqrt(.5), atol=1e-13, rtol=0)
        np.testing.assert_allclose(panels.periodic_bilinear(q, values, q), values, atol=1e-13, rtol=0)
        point = np.array([[.98, .01, 0.], [.23, -.02, 0.]])
        np.testing.assert_allclose(panels.periodic_bilinear(q, values, point),
                                   panels.periodic_bilinear(q, values, point + [3., -2., 0.]), atol=1e-13, rtol=0)

    def test_missing_native_grid_sample_is_rejected(self):
        q = grid()[:-1]
        with self.assertRaises(ValueError):
            panels.periodic_bilinear(q, np.ones(len(q)), np.array([[0., 0., 0.]]))

    def test_unknown_neighbour_propagates_without_contaminating_exact_known_node(self):
        q = grid(centered=False)
        values = np.full(len(q), 7.)
        values[0] = np.nan
        # The cell interior uses the undefined origin, while its opposite known node does not.
        target = np.array([[.125, .125, 0.], [.25, .25, 0.], [.25, 0., 0.], [0., 0., 0.]])
        actual = panels.periodic_bilinear(q, values, target)
        np.testing.assert_array_equal(np.isnan(actual), [True, False, False, True])
        np.testing.assert_allclose(actual[[1, 2]], [7., 7.], atol=1e-13, rtol=0)
        np.testing.assert_allclose(actual, panels.periodic_bilinear(q, values, target + [2., -3., 0.]),
                                   atol=1e-13, rtol=0, equal_nan=True)

    def test_unknown_mask_tracks_values_when_grid_is_reordered(self):
        q = grid(centered=False)
        values = 2 + q[:, 0] + q[:, 1]
        values[5] = np.nan
        order = np.random.default_rng(103).permutation(len(q))
        target = np.array([[.3, .3, 0.], [.6, .6, 0.]])
        np.testing.assert_allclose(panels.periodic_bilinear(q, values, target),
                                   panels.periodic_bilinear(q[order], values[order], target),
                                   atol=1e-13, rtol=0, equal_nan=True)

    def test_infinite_values_are_rejected_and_unknown_field_remains_unknown(self):
        q = grid()
        with self.assertRaisesRegex(ValueError, "infinite"):
            panels.periodic_bilinear(q, np.full(len(q), np.inf), q[:1])
        self.assertTrue(np.isnan(panels.periodic_bilinear(q, np.full(len(q), np.nan), q[:1])).all())


class PanelDataIntegrityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name)
        self.poscar = self.path / "POSCAR"
        self.poscar.write_text("square\n1\n1 0 0\n0 1 0\n0 0 5\nH\n1\nDirect\n0 0 0\n")
        self.reciprocal = np.diag([2 * np.pi, 2 * np.pi, 2 * np.pi / 5])
        self.qpath = np.array([[-.25, 0., 0.], [0., 0., 0.], [.25, 0., 0.]])
        self.bands = self.path / "bands.csv"
        self.rows = band_rows(self.qpath, self.reciprocal)
        write_bands(self.bands, self.rows)
        self.qmesh = grid(centered=False)
        self.native_values = np.arange(1., len(self.qmesh) + 1.)
        gaps = np.full(len(self.qmesh), 2.)
        gaps[0] = 1e-9
        self.mesh = self.path / "mesh.csv"
        self.path_kubo = self.path / "path.csv"
        write_kubo(self.mesh, self.qmesh, self.native_values, gaps)
        write_kubo(self.path_kubo, self.qpath, [-2., 1e6, 2.], [2., 1e-9, 2.])

    def call_plot(self, **updates):
        args = dict(method="kubo", input_path=self.mesh, path_input=self.path_kubo,
                    poscar_path=self.poscar, bands_csv=self.bands, output_path=self.path / "figure.png",
                    occupied=1, band=1, node_indices=(1, 2, 3), node_labels=("A", "Gamma", "B"),
                    curve_output=self.path / "curve.csv")
        args.update(updates)
        return panels.plot_panels(**args)

    def test_band_table_preserves_energies_and_uses_cartesian_path_length(self):
        q, distance, energies = panels.read_bands(self.bands, self.reciprocal)
        np.testing.assert_allclose(q, self.qpath)
        np.testing.assert_allclose(distance, [0., np.pi / 2, np.pi])
        np.testing.assert_allclose(energies, [[-.9375, 1.0625], [-1., 1.], [-.9375, 1.0625]])

    def test_duplicate_band_row_is_rejected_even_if_identical(self):
        write_bands(self.bands, self.rows + [dict(self.rows[0])])
        with self.assertRaisesRegex(ValueError, "unique"):
            panels.read_bands(self.bands, self.reciprocal)

    def test_nonfirst_band_nonfinite_distance_is_rejected(self):
        self.rows[1]["path_inv_A"] = float("nan")
        write_bands(self.bands, self.rows)
        with self.assertRaises(ValueError):
            panels.read_bands(self.bands, self.reciprocal)

    def test_band_distance_from_another_reciprocal_metric_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "distance"):
            panels.read_bands(self.bands, np.eye(3))

    def test_wrong_or_missing_kubo_normalization_is_rejected(self):
        for normalization in ["LEGACY_MINUS_IM", "", "STANDARD_MINUS_FOUR_IM"]:
            with self.subTest(normalization=normalization):
                write_kubo(self.path_kubo, self.qpath, [-1., 0., 1.], normalization=normalization)
                with self.assertRaisesRegex(ValueError, "normalization"):
                    panels.read_kubo(self.path_kubo, 1, 1e-5)

    def test_invalid_degeneracy_threshold_is_rejected(self):
        for threshold in [-1., float("nan"), float("inf")]:
            with self.subTest(threshold=threshold), self.assertRaises(ValueError):
                panels.read_kubo(self.path_kubo, 1, threshold)

    def test_kubo_path_energy_must_match_displayed_band(self):
        self.rows[2]["energy_eV"] += .02
        write_bands(self.bands, self.rows)
        with self.assertRaisesRegex(ValueError, "energies disagree"):
            self.call_plot()
        self.assertFalse((self.path / "figure.png").exists())

    def test_kubo_path_coordinates_must_match_displayed_band(self):
        shifted = self.qpath + [.01, 0., 0.]
        write_kubo(self.path_kubo, shifted, [-2., 0., 2.])
        with self.assertRaisesRegex(ValueError, "coordinates disagree"):
            self.call_plot()

    def test_kubo_map_and_path_must_use_the_same_intermediate_band_window(self):
        write_kubo(self.path_kubo, self.qpath, [-2., 0., 2.], intermediate=3)
        with self.assertRaisesRegex(ValueError, "intermediate_bands"):
            self.call_plot()

    def test_native_fukui_group_cannot_be_relabelled_as_another_band_group(self):
        source = ROOT / "examples/features/fukui-berry-curvature/reference/BERRYCURV.dat"
        with self.assertRaisesRegex(ValueError, "band range"):
            self.call_plot(method="fukui", input_path=source, path_input=None, occupied=1)

    def test_native_values_and_degeneracy_mask_survive_map_clipping(self):
        sources = [self.mesh, self.path_kubo, self.bands]
        before = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
        with mock.patch.object(panels, "PolyCollection", wraps=panels.PolyCollection) as collection:
            metadata = self.call_plot()
        displayed = collection.call_args.kwargs["array"]
        polygons = collection.call_args.args[0]
        self.assertTrue(np.ma.getmaskarray(displayed).any())
        self.assertEqual(set(displayed.compressed()), set(self.native_values[1:]))
        areas = np.array([polygon_area(p) for p in polygons])
        self.assertAlmostEqual(float(areas.sum()), (2 * np.pi) ** 2, places=9)
        # The single undefined native node occupies exactly one of the 16 equal sampling cells.
        self.assertAlmostEqual(float(areas[np.ma.getmaskarray(displayed)].sum()), (2 * np.pi) ** 2 / 16, places=9)
        self.assertEqual(metadata["masked_map_points"], 1)
        self.assertEqual(metadata["masked_path_points"], 1)
        with (self.path / "curve.csv").open() as handle:
            output = list(csv.DictReader(handle))
        self.assertEqual([r["valid"] for r in output], ["1", "0", "1"])
        self.assertEqual(output[1]["omega_z_A2"], "")
        self.assertEqual([float(output[i]["omega_z_A2"]) for i in [0, 2]], [-2., 2.])
        self.assertEqual(before, {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})

    def test_smooth_map_preserves_input_files_and_direct_path_curve(self):
        sources = [self.mesh, self.path_kubo, self.bands]
        before = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
        captured = []
        original = panels.plt.Axes.pcolormesh

        def capture(axis, *args, **kwargs):
            if len(args) >= 3 and np.shape(args[2]) == (25, 25):
                captured.append(np.ma.asarray(args[2]).copy())
            return original(axis, *args, **kwargs)

        with mock.patch.object(panels.plt.Axes, "pcolormesh", new=capture):
            metadata = self.call_plot(map_style="smooth", display_grid=25)
        self.assertEqual(len(captured), 1)
        displayed = captured[0]
        self.assertTrue(np.ma.getmaskarray(displayed).any())
        self.assertGreaterEqual(displayed.min(), self.native_values[1:].min())
        self.assertLessEqual(displayed.max(), self.native_values.max())
        self.assertEqual(metadata["map_interpolation"], "periodic_bilinear")
        self.assertEqual(metadata["display_grid"], [25, 25])
        self.assertEqual(metadata["masked_map_points"], 1)
        with (self.path / "curve.csv").open() as handle:
            output = list(csv.DictReader(handle))
        self.assertEqual([r["valid"] for r in output], ["1", "0", "1"])
        self.assertEqual(output[1]["omega_z_A2"], "")
        self.assertEqual([float(output[i]["omega_z_A2"]) for i in [0, 2]], [-2., 2.])
        self.assertEqual(before, {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})

    def test_bundle_reader_retains_finite_trace_and_checks_occupied_rank(self):
        write_bundle(self.mesh, self.qmesh, self.native_values)
        q, values, gaps = panels.read_bundle(self.mesh, 1, 1e-5)
        np.testing.assert_allclose(q, self.qmesh)
        np.testing.assert_array_equal(values, self.native_values)
        np.testing.assert_array_equal(gaps, np.full(len(q), 2.))
        for change in [{"band_min": "2"}, {"band_max": "2"}, {"band_rank": "2"}]:
            with self.subTest(change=change):
                write_bundle(self.mesh, self.qmesh, self.native_values, metadata=change)
                with self.assertRaisesRegex(ValueError, "range"):
                    panels.read_bundle(self.mesh, 1, 1e-5)

    def test_bundle_reader_rejects_failed_or_noncurrent_result(self):
        for change in [{"schema": "VASPBERRY_BARE_MOMENTUM_KUBO_V2"},
                       {"normalization": "STANDARD_MINUS_FOUR_IM"}, {"result_status": "FAIL"},
                       {"operator": "UNKNOWN_OPERATOR"}, {"berry_connection": "A=-i<u|dk u>"},
                       {"intermediate_bands": "1:2"}, {"result_kind": "SUM_OF_MASKED_BAND_ROWS"}]:
            with self.subTest(change=change):
                write_bundle(self.mesh, self.qmesh, self.native_values, metadata=change)
                with self.assertRaises(ValueError):
                    panels.read_bundle(self.mesh, 1, 1e-5)

    def test_bundle_reader_requires_consistent_soc_spin_one_rows(self):
        cases = [np.full(len(self.qmesh), 2), np.zeros(len(self.qmesh), dtype=int),
                 np.array([1] * (len(self.qmesh) - 1) + [2])]
        for spins in cases:
            with self.subTest(spins=spins.tolist()):
                write_bundle(self.mesh, self.qmesh, self.native_values, spins=spins)
                with self.assertRaises(ValueError):
                    panels.read_bundle(self.mesh, 1, 1e-5)

    def test_bundle_reader_rejects_nonisolated_or_nonfinite_external_gap(self):
        for gap in [0., 1e-5, -1., float("nan"), float("inf")]:
            with self.subTest(gap=gap):
                gaps = np.full(len(self.qmesh), 2.)
                gaps[3] = gap
                write_bundle(self.mesh, self.qmesh, self.native_values, gaps)
                with self.assertRaises(ValueError):
                    panels.read_bundle(self.mesh, 1, 1e-5)

    def test_bundle_reader_rejects_duplicate_k_index(self):
        write_bundle(self.mesh, self.qmesh, self.native_values)
        lines = self.mesh.read_text().splitlines()
        with self.mesh.open("a") as handle:
            handle.write(lines[-1] + "\n")
        with self.assertRaisesRegex(ValueError, "one selected-spin row"):
            panels.read_bundle(self.mesh, 1, 1e-5)

    def test_bundle_panels_preserve_trace_and_matching_direct_path_values(self):
        write_bundle(self.mesh, self.qmesh, self.native_values)
        write_bundle(self.path_kubo, self.qpath, [-3., .01, 3.])
        before = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in [self.mesh, self.path_kubo, self.bands]}
        metadata = self.call_plot(method="kubo-bundle", map_style="smooth", display_grid=25)
        self.assertEqual(metadata["map_quantity"], "pointwise_Kubo_bundle_trace")
        self.assertEqual(metadata["masked_path_points"], 0)
        self.assertIsNone(metadata["kubo_band"])
        with (self.path / "curve.csv").open() as handle:
            output = list(csv.DictReader(handle))
        np.testing.assert_array_equal([float(r["omega_z_A2"]) for r in output], [-3., .01, 3.])
        self.assertEqual(before, {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in before})

    def test_bundle_path_gap_must_match_displayed_band_energies(self):
        write_bundle(self.mesh, self.qmesh, self.native_values)
        write_bundle(self.path_kubo, self.qpath, [-3., 0., 3.], gaps=[2., 2.1, 2.])
        with self.assertRaisesRegex(ValueError, "gap"):
            self.call_plot(method="kubo-bundle")
        self.assertFalse((self.path / "figure.png").exists())

    def test_small_bundle_ignores_unused_default_band_but_single_band_checks_bounds(self):
        write_bundle(self.mesh, self.qmesh, self.native_values)
        write_bundle(self.path_kubo, self.qpath, [-3., 0., 3.])
        # Omit band deliberately: the default is 18, while this fixture has only two bands.
        metadata = panels.plot_panels(
            method="kubo-bundle", input_path=self.mesh, path_input=self.path_kubo,
            poscar_path=self.poscar, bands_csv=self.bands,
            output_path=self.path / "small_bundle.png", occupied=1,
            node_indices=(1, 2, 3), node_labels=("A", "Gamma", "B"),
        )
        self.assertIsNone(metadata["kubo_band"])
        self.assertEqual(metadata["occupied_bands"], [1, 1])
        with self.assertRaisesRegex(ValueError, "outside band table"):
            self.call_plot(method="kubo", band=18)


class RealMoS2PanelInputTests(unittest.TestCase):
    def test_real_fukui_display_cut_respects_time_reversal_and_periodicity(self):
        feature = ROOT / "examples/features/fukui-berry-curvature"
        reciprocal = panels.reciprocal_from_poscar(feature / "inputs/POSCAR")
        q, values = panels.read_native_curvature(feature / "reference/BERRYCURV.dat", reciprocal)
        q, values, mesh = panels.uniform_plaquettes(q, values)
        self.assertEqual(mesh, (12, 12))
        path = np.column_stack([np.linspace(1 / 3, -1 / 3, 49), np.linspace(2 / 3, -2 / 3, 49), np.zeros(49)])
        cut = panels.periodic_bilinear(q, values, path)
        np.testing.assert_allclose(cut, -cut[::-1], atol=2.1e-4, rtol=0)
        np.testing.assert_allclose(cut, panels.periodic_bilinear(q, values, path + [1, -2, 0]), atol=1e-11, rtol=0)
        self.assertGreater(np.max(np.abs(cut)), 1.)
        self.assertLessEqual(np.max(np.abs(cut)), np.max(np.abs(values)))

    def test_real_archived_kubo_has_valid_valleys_and_masked_gamma(self):
        source = ROOT / "examples/features/kubo-curvature/reference/KUBO.csv"
        q, values, energies = panels.read_kubo(source, 18, 1e-5)
        gamma = np.max(np.abs(q), axis=1) < 1e-12
        self.assertTrue(gamma.any())
        self.assertTrue(np.isnan(values[gamma]).all())
        self.assertTrue(np.isfinite(values[[0, -1]]).all())
        self.assertAlmostEqual(float(values[0]), -float(values[-1]), delta=2e-5)
        self.assertTrue(np.isfinite(energies).all())


if __name__ == "__main__":
    unittest.main()
