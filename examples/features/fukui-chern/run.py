#!/usr/bin/env python3
"""Recalculate QWZ eigenstates and production Fukui plaquette/Chern results."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT/'tools'))
from wavecar_fukui import LinkSet, infer_uniform_grid, plaquette_flux
from vaspberry_transport import CurvatureData, validate_fukui_geometry, integrate_sigma


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def calculate(nx, ny, mass):
    x, y = np.meshgrid(np.arange(nx)/nx, np.arange(ny)/ny, indexing='ij')
    q = np.column_stack((x.ravel(), y.ravel(), np.zeros(nx*ny)))
    kx, ky = 2*np.pi*x, 2*np.pi*y
    h = np.zeros((nx, ny, 2, 2), dtype=np.complex128)
    h[:, :, 0, 0] = mass+np.cos(kx)+np.cos(ky)
    h[:, :, 1, 1] = -h[:, :, 0, 0]
    h[:, :, 0, 1] = np.sin(kx)-1j*np.sin(ky)
    h[:, :, 1, 0] = h[:, :, 0, 1].conj()
    energies, vectors = np.linalg.eigh(h)
    states = vectors[:, :, :, 0]
    links = []
    for axis in (0, 1):
        overlap = np.sum(states.conj()*np.roll(states, -1, axis=axis), axis=-1)
        magnitude = abs(overlap)
        if not np.isfinite(magnitude).all() or magnitude.min() <= 1e-10:
            raise ValueError('singular analytic-model link; refine the mesh or change the model')
        links.append(LinkSet(np.angle(overlap), np.log(magnitude), magnitude, magnitude,
                             np.ones_like(magnitude), np.ones_like(magnitude)))
    flux = plaquette_flux(*links, infer_uniform_grid(q, nx, ny))
    center = q+np.array([.5/nx, .5/ny, 0.])
    area = (2*np.pi)**2/(nx*ny)
    lower = energies[:, :, 0]
    vertex_energies = np.stack((lower, np.roll(lower, -1, axis=0),
        np.roll(np.roll(lower, -1, axis=0), -1, axis=1), np.roll(lower, -1, axis=1)), axis=-1)
    gap = float(np.min(energies[:, :, 1]-lower))
    data = CurvatureData(cart=center*2*np.pi, frac=center, omega=flux.ravel()/area,
        vertex_energies=vertex_energies.reshape(nx*ny, 4),
        metadata={'b1': np.array([2*np.pi, 0., 0.]), 'b2': np.array([0., 2*np.pi, 0.]),
                  'band': 1, 'band_mode': 'single-fukui', 'nk': nx*ny, 'k_grid': (nx, ny),
                  'dk_area': area, 'min_direct_band_gap_eV': gap})
    _, chern = validate_fukui_geometry(data)
    empty_mu = float(lower.min()-1.)
    transport = integrate_sigma([data], np.array([empty_mu, 0.]), 0.)
    sigma = transport['sigma_xy_total_e2_over_h']
    if not np.isclose(sigma[0], 0., atol=1e-12, rtol=0) or not np.isclose(sigma[1], -chern, atol=1e-12, rtol=0):
        raise ValueError('empty/full lower-band Hall limits do not reproduce 0 and -C')
    expected = 1 if -2 < mass < 0 else -1 if 0 < mass < 2 else 0
    if not np.isclose(chern, expected, atol=1e-10, rtol=0):
        raise ValueError('Fukui result disagrees with the specified gapped QWZ phase')
    summary = {'mass': mass, 'mesh_nx': nx, 'mesh_ny': ny, 'band': 1,
               'chern': chern, 'expected_chern': expected,
               'absolute_error': abs(chern-expected),
               'minimum_sampled_gap_eV': gap,
               'empty_reference_mu_eV': empty_mu,
               'lower_band_sigma_empty_e2h': float(transport['sigma_xy_total_e2_over_h'][0]),
               'lower_band_sigma_mu0_e2h': float(transport['sigma_xy_total_e2_over_h'][1]),
               'minimum_link_magnitude': float(min(link.min_singular.min() for link in links)),
               'maximum_abs_plaquette_flux_rad': float(abs(flux).max())}
    return summary, center, flux


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=HERE/'input.json')
    parser.add_argument('--output-dir', type=Path, required=True, help='fresh result directory')
    args = parser.parse_args()
    if args.output_dir.exists():
        parser.error('output directory exists; use a fresh path')
    config = json.loads(args.input.read_text())
    mesh, masses = config['mesh'], config['masses']
    if not (isinstance(mesh, list) and len(mesh) == 2 and all(type(v) is int and v >= 4 for v in mesh)
            and isinstance(masses, list) and masses and np.isfinite(masses).all()
            and not any(m in (-2., 0., 2.) for m in masses)):
        parser.error('mesh needs two integers >=4; masses must be finite and avoid -2, 0, 2')
    if config.get('selected_band') != 1 or config.get('lattice_constant_A') != 1.0 or config.get('energy_unit_eV') != 1.0:
        parser.error('this explicit fixture uses lower band 1, a=1 Angstrom and the 1 eV model scale')
    cases = [calculate(*mesh, float(mass)) for mass in masses]
    declared = config.get('expected_lower_band_chern')
    if declared is not None and declared != [case[0]['expected_chern'] for case in cases]:
        parser.error('input expected_lower_band_chern contradicts its gapped QWZ masses')
    args.output_dir.mkdir(parents=True)
    with (args.output_dir/'summary.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(cases[0][0])); writer.writeheader()
        writer.writerows(case[0] for case in cases)
    with (args.output_dir/'plaquettes.csv').open('w', newline='') as f:
        writer = csv.writer(f); writer.writerow(['mass', 'cell_id', 'q1_center', 'q2_center', 'flux_rad'])
        for summary, centers, flux in cases:
            for index, (center, value) in enumerate(zip(centers, flux.ravel()), 1):
                writer.writerow([summary['mass'], index, *center[:2], value])
    fig, axes = plt.subplots(1, len(cases), figsize=(4*len(cases), 3.8), layout='constrained', squeeze=False)
    limit = max(float(abs(case[2]).max()) for case in cases)
    for ax, (summary, _, flux) in zip(axes[0], cases):
        plot = ax.imshow(flux.T, origin='lower', extent=(0, 1, 0, 1), interpolation='none',
                         cmap='RdBu_r', vmin=-limit, vmax=limit)
        ax.set(xlabel='fractional k1', ylabel='fractional k2',
               title=f"mass = {summary['mass']:g}; C = {summary['chern']:.6f}")
    fig.colorbar(plot, ax=axes.ravel().tolist(), label='Fukui plaquette flux (rad)', shrink=.85)
    fig.suptitle('Analytic QWZ eigenstates → production Fukui plaquettes')
    fig.savefig(args.output_dir/'figure.png', dpi=180); plt.close(fig)
    result = {'schema_version': 1, 'feature_id': 'fukui-chern', 'status': 'PASS',
              'workflow_mode': 'new_analytic_model_calculation', 'mode': 'new_analytic_model_calculation',
              'software_version': (ROOT/'VERSION').read_text().strip(), 'input': config,
              'cases': [case[0] for case in cases],
              'numerical_checks': {'status': 'PASS', 'expected_chern': 'PASS', 'isolated_links_and_flux_geometry': 'PASS',
                                   'empty_band_sigma_equals_zero': 'PASS', 'mu0_sigma_equals_minus_chern': 'PASS'},
              'production_functions': ['wavecar_fukui.plaquette_flux', 'vaspberry_transport.validate_fukui_geometry',
                                       'vaspberry_transport.integrate_sigma'],
              'provenance': {'input_sha256': sha(args.input), 'runner_sha256': sha(__file__),
                             'wavecar_fukui_sha256': sha(ROOT/'tools/wavecar_fukui.py'),
                             'transport_sha256': sha(ROOT/'tools/vaspberry_transport.py')},
              'limitations': ['Analytic orthonormal two-orbital model, not a VASP or material calculation.',
                              'Model eigenvector overlaps are generated here; no WAVECAR parser or PAW overlap is exercised.',
                              'Flux is cell-integrated radians, not Kubo point curvature in Angstrom squared.',
                              'Hall values are the lower-band contribution with four-vertex occupations at T=0; the upper band is empty at mu=0 in this model.',
                              'A resolved known gapped model is checked; integer output alone is not a material convergence proof.'],
              'output_sha256': {name: sha(args.output_dir/name) for name in ('summary.csv', 'plaquettes.csv', 'figure.png')}}
    (args.output_dir/'result.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status': 'PASS', 'chern': [case[0]['chern'] for case in cases]}))


if __name__ == '__main__':
    main()
