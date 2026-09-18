"""Physical oracles, composable CLI, migration and malformed-data regression."""
import copy
from dataclasses import replace
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tools'))
import vaspberry_kubo as cli
from berry_data import (hall_spectrum, read_curvature, regions_from_spec, validate,
                       write_curvature, NORMALIZATION)
from vaspberry_transport import fermi_dirac


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.run_cli('demo', '--mesh', '16', '--output-dir', self.root/'model')
        self.run_cli('matrix', '--matrices', self.root/'model/matrix.npz',
            '--metadata', self.root/'model/matrix.json', '--n-bands', '1:2', '--m-bands', '1:2',
            '--degeneracy-threshold-eV', '1e-9', '--mesh', '16', '16',
            '--spin-multiplicity', '1', '--energy-reference', 'model zero', '--output-dir', self.root/'curv')
        self.data = read_curvature(self.root/'curv')

    def run_cli(self, *args):
        cli.main([str(a) for a in args])

    def test_qwz_independent_d_vector_oracle_and_units(self):
        x,y=(2*np.pi*self.data.kpoints_fractional[:,:2]).T
        d=np.column_stack((np.sin(x),np.sin(y),-1+np.cos(x)+np.cos(y)))
        dx=np.column_stack((np.cos(x),np.zeros_like(x),-np.sin(x)))
        dy=np.column_stack((np.zeros_like(x),np.cos(y),-np.sin(y)))
        lower=np.einsum('ka,ka->k',d,np.cross(dx,dy))/(2*np.linalg.norm(d,axis=1)**3)
        np.testing.assert_allclose(self.data.omega_A2[:,0,2],lower,atol=1e-14)
        rows,meta=hall_spectrum(self.data,[0.],[0.],mu_reference=0.)
        self.assertAlmostEqual(rows[0]['sigma_e2_over_h'],-1.,delta=1e-4)
        self.assertEqual(meta['zero_T_equal_energy_occupation'],1.)

    def test_zero_T_and_finite_T_match_direct_quadrature_at_boundaries(self):
        mus=np.unique(np.r_[-.5,0.,.5,self.data.energies_eV[3,0]])
        rows,_=hall_spectrum(self.data,mus,[0.,100.,300.],mu_reference=0.,mu_chunk=2)
        for row in rows:
            if row['region']!='total':continue
            f=fermi_dirac(self.data.energies_eV,row['mu_eV'],row['temperature_K'])
            expected=-2*np.pi*np.mean(np.sum(f*self.data.omega_A2[:,:,2],axis=1))
            self.assertAlmostEqual(row['sigma_e2_over_h'],expected,places=13)
            self.assertAlmostEqual(row['electrons_per_cell'],np.mean(f.sum(axis=1)),places=13)

    def test_partition_difference_band_sum_and_multiplicity(self):
        spec={'regions':[{'name':'A','k_ids':list(range(1,65))},
                         {'name':'B','k_ids':list(range(65,129))}]}
        rows,meta=hall_spectrum(self.data,[0.],[0.],mu_reference=-.5,region_spec=spec,
                               differences=[('contrast','A','B')],band_resolved=True)
        values={(r['region'],r['band_id']):r['sigma_e2_over_h'] for r in rows}
        self.assertAlmostEqual(values['total',0],sum(values[name,0] for name in ('A','B','rest')),places=13)
        self.assertAlmostEqual(values['contrast',0],values['A',0]-values['B',0],places=13)
        self.assertAlmostEqual(values['total',0],values['total',1]+values['total',2],places=13)
        doubled=replace(self.data,metadata=dict(self.data.metadata,spin_multiplicity=2))
        rows2,_=hall_spectrum(doubled,[0.],[0.],mu_reference=0.)
        self.assertAlmostEqual(rows2[0]['sigma_e2_over_h'],2*values['total',0],places=13)

    def test_plane_rotation_and_orientation(self):
        # Rotate physical x->y->z, and keep fractional sampling axes unchanged.
        rot=np.array([[0,1,0],[0,0,1],[1,0,0.]])
        d=replace(self.data,lattice_A=self.data.lattice_A@rot,
                  reciprocal_inv_A=self.data.reciprocal_inv_A@rot,omega_A2=self.data.omega_A2@rot)
        before,_=hall_spectrum(self.data,[0.],[0.],mu_reference=0.)
        after,_=hall_spectrum(d,[0.],[0.],mu_reference=0.)
        self.assertAlmostEqual(before[0]['sigma_e2_over_h'],after[0]['sigma_e2_over_h'],places=13)
        m=copy.deepcopy(d.metadata);m['sampling']['plane_axes']=[1,0]
        reversed_rows,_=hall_spectrum(replace(d,metadata=m),[0.],[0.],mu_reference=0.)
        self.assertAlmostEqual(before[0]['sigma_e2_over_h'],-reversed_rows[0]['sigma_e2_over_h'],places=13)

    def test_partial_band_selection_requires_explicit_scope(self):
        d=replace(self.data,band_ids=self.data.band_ids[:1],energies_eV=self.data.energies_eV[:,:1],
                  omega_A2=self.data.omega_A2[:,:1],min_gap_eV=self.data.min_gap_eV[:,:1],
                  valid_nondegenerate=self.data.valid_nondegenerate[:,:1])
        with self.assertRaisesRegex(ValueError,'incomplete occupied-band'):
            hall_spectrum(d,[0.],[0.],mu_reference=0.)
        _,m=hall_spectrum(d,[0.],[0.],mu_reference=0.,allow_partial_bands=True)
        self.assertEqual(m['scope'],'partial_band_contribution')

    def test_guard_occupied_top_nonfinite_and_negative_T(self):
        for mus,temps in [([3.],[0.]),([np.nan],[0.]),([0.],[-1.])]:
            with self.subTest(mus=mus,temps=temps),self.assertRaises(ValueError):
                hall_spectrum(self.data,mus,temps,mu_reference=0.)

    def test_small_doping_response_survives_huge_filled_baseline(self):
        energies=np.broadcast_to(np.array([-10.,.25]),self.data.energies_eV.shape).copy()
        omega=np.zeros_like(self.data.omega_A2);omega[:,0,2]=1e20;omega[:,1,2]=1.
        data=replace(self.data,energies_eV=energies,omega_A2=omega,
                     min_gap_eV=np.full_like(energies,10.25))
        rows,_=hall_spectrum(data,[.5],[0.,100.],mu_reference=0.,allow_partial_bands=True)
        for r in rows:
            if r['region']!='total':continue
            expected=-2*np.pi*(fermi_dirac(np.array([.25]),.5,r['temperature_K'])[0]
                              -fermi_dirac(np.array([.25]),0.,r['temperature_K'])[0])
            self.assertAlmostEqual(r['delta_sigma_e2_over_h'],expected,places=12)

    def test_grid_duplicate_path_nonuniform_weights_rejected(self):
        points=self.data.kpoints_fractional.copy();points[0]=points[1]
        with self.assertRaisesRegex(ValueError,'duplicate'):
            validate(replace(self.data,kpoints_fractional=points))
        weights=self.data.weights.copy();weights[:2]+=[.0001,-.0001]
        with self.assertRaisesRegex(ValueError,'uniform'):
            validate(replace(self.data,weights=weights))
        meta=copy.deepcopy(self.data.metadata);meta['sampling']={'kind':'points'}
        with self.assertRaisesRegex(ValueError,'full 2D mesh'):
            hall_spectrum(replace(self.data,metadata=meta),[0.],[0.],mu_reference=0.)

    def test_nan_masks_reject_integration_and_checksum_rejects_corruption(self):
        valid=self.data.valid_nondegenerate.copy();valid[0,0]=False
        omega=self.data.omega_A2.copy();omega[0,0]=np.nan
        with self.assertRaisesRegex(ValueError,'invalid/degenerate'):
            hall_spectrum(replace(self.data,valid_nondegenerate=valid,omega_A2=omega),[0.],[0.],mu_reference=0.)
        with (self.root/'curv/curvature.npz').open('ab') as f:f.write(b'corrupt')
        with self.assertRaisesRegex(ValueError,'SHA256'):
            read_curvature(self.root/'curv')

    def test_region_validation_periodic_translation(self):
        def spec(center):return {'regions':[{'name':'pocket','center_fractional':center,'radius_inv_A':.5}]}
        a=regions_from_spec(self.data,spec([0,0,0]))
        b=regions_from_spec(self.data,spec([1,-1,0]))
        np.testing.assert_array_equal(a['pocket'],b['pocket'])
        with self.assertRaisesRegex(ValueError,'overlap'):
            regions_from_spec(self.data,{'regions':[{'name':'a','k_ids':[1]},{'name':'b','k_ids':[1]}]})
        with self.assertRaisesRegex(ValueError,'sampled reciprocal plane'):
            regions_from_spec(self.data,spec([0,0,.1]))

    def test_hall_cli_serialization_no_pickle_no_overwrite(self):
        self.run_cli('hall','--curvature',self.root/'curv','--mu-min','-.5','--mu-max','.5',
                     '--mu-num','5','--mu-reference','0','--temperatures','0','100',
                     '--band-resolved','--output-dir',self.root/'hall')
        with np.load(self.root/'hall/conductivity.npz',allow_pickle=False) as a:
            self.assertEqual(len(a['sigma_e2_over_h']),60)
            self.assertEqual(a['region'].dtype.kind,'U')
        with self.assertRaisesRegex(ValueError,'exists'):
            write_curvature(self.root/'curv',self.data)

    def test_legacy_conversion_once_and_physical_values_unchanged(self):
        class Header:
            lattice=np.eye(3);reciprocal=2*np.pi*np.eye(3);ispin=1
        class FakeWavecar:
            header=Header();energies=self.data.energies_eV;kpoints=self.data.kpoints_fractional
            def coefficients(self,*args):return None
        wave=self.root/'WAVECAR';wave.write_bytes(b'test fixture')
        for norm,col,factor in [('legacy-double','omega_legacy_A2',2.),('physical','omega_z_A2',1.)]:
            path=self.root/(norm+'.csv')
            import csv
            with path.open('w') as f:
                if norm=='physical':f.write('# normalization=STANDARD_MINUS_TWO_IM\n')
                w=csv.writer(f);w.writerow(['k_index','band','kx_frac','ky_frac','kz_frac','energy_eV',col,'min_gap_eV'])
                for k,n in np.ndindex(self.data.energies_eV.shape):
                    w.writerow([k+1,n+1,*self.data.kpoints_fractional[k],self.data.energies_eV[k,n],
                                factor*self.data.omega_A2[k,n,2],self.data.min_gap_eV[k,n]])
            with patch.object(cli,'Wavecar',return_value=FakeWavecar()):
                self.run_cli('import-legacy','--csv',path,'--wavecar',wave,'--normalization',norm,
                    '--mesh','16','16','--spinor-components','2','--spin-multiplicity','1',
                    '--energy-reference','model zero','--degeneracy-threshold-eV','1e-9',
                    '--output-dir',self.root/norm)
            data=read_curvature(self.root/norm)
            np.testing.assert_array_equal(data.omega_A2[:,:,2],self.data.omega_A2[:,:,2])
            self.assertTrue(np.isnan(data.omega_A2[:,:,:2]).all())
            self.assertEqual(data.metadata['normalization'],NORMALIZATION)
        args=cli.parser().parse_args(['import-legacy','--csv',str(self.root/'physical.csv'),
             '--wavecar',str(wave),'--normalization','legacy-double','--spinor-components','2',
             '--spin-multiplicity','1','--energy-reference','model','--degeneracy-threshold-eV','1e-9',
             '--output-dir',str(self.root/'bad')])
        with patch.object(cli,'Wavecar',return_value=FakeWavecar()),self.assertRaisesRegex(ValueError,'divided by two again'):
            cli.import_legacy(args)


if __name__=='__main__':unittest.main()
