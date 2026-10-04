"""Independent finite-band pair oracles and public selected optical workflows."""
import contextlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from waveder_selected import (parse_bands, optical_pairs, geometric_trace,
                              occupation_response, logistic_difference)
from vaspberry_transport import fermi_dirac
import waveder_hall as wh
import vaspberry_kubo as cli
from postprocess_config import load_settings
import vaspberry_post as post
import test_waveder_hall as occupied_fixture


def hermitian_fixture(nk=4, nb=4, seed=96):
    rng = np.random.default_rng(seed)
    a = rng.normal(size=(nk, 3, nb, nb))+1j*rng.normal(size=(nk, 3, nb, nb))
    return ((a+a.conj().swapaxes(-1,-2))/2).astype(np.complex64)


def band_oracle(c, energies, bands, mu, temperature):
    result = np.zeros((len(energies),3))
    occ = fermi_dirac(energies, mu, temperature)
    for k in range(len(energies)):
        for n in bands:
            n -= 1
            for m in range(c.shape[2]):
                if m == n: continue
                z = c[k,:,m,n].astype(complex)
                result[k] += occ[k,n]*np.array([-2*(z[a].conjugate()*z[b]).imag
                                          for a,b in ((1,2),(2,0),(0,1))])
    return result


class SelectedPairTests(unittest.TestCase):
    def setUp(self):
        self.c = hermitian_fixture()
        self.e = np.tile([-2.,-.2,.3,2.], (4,1))

    def test_band_parser_and_bad_selections(self):
        self.assertEqual(parse_bands('31,33:34'), [31,33,34])
        self.assertEqual(parse_bands('3,1:2'), [1,2,3])
        for text in ('0','2:1','1,1','1:2,2','1,,2','1 2','-1','1:',''):
            with self.subTest(text=text), self.assertRaises(ValueError): parse_bands(text)

    def test_weighted_finiteT_and_delta_match_resolved_band_oracle(self):
        p = optical_pairs(self.c, self.e)
        for bands in ([1],[2,3],[1,3,4]):
            for t in (0., 500.):
                for mu in (-.2,.25,1.):
                    value,delta,_,_,_ = occupation_response(p,bands,mu,t,-.6)
                    expected=band_oracle(self.c,self.e,bands,mu,t)
                    reference=band_oracle(self.c,self.e,bands,-.6,t)
                    np.testing.assert_allclose(value,expected,atol=2e-14)
                    np.testing.assert_allclose(delta,expected-reference,atol=2e-14)

    def test_internal_unequal_weights_do_not_use_bundle_mean_occupation(self):
        p=optical_pairs(self.c,self.e)
        value,*_=occupation_response(p,[2,3],0.,0.,-1.)
        expected=band_oracle(self.c,self.e,[2],0.,0.)
        np.testing.assert_allclose(value,expected,atol=1e-14)
        self.assertGreater(np.max(abs(value-geometric_trace(p,[2,3]))),.1)

    def test_reverse_orientation_high_complete_bundle_and_missing_internal_weights(self):
        p=optical_pairs(self.c[:,:,:,:2],self.e)
        expected=geometric_trace(optical_pairs(self.c,self.e),[3,4])
        np.testing.assert_allclose(geometric_trace(p,[3,4]),expected,atol=1e-14)
        value,*_=occupation_response(p,[3,4],3.,0.,3.)
        np.testing.assert_allclose(value,expected,atol=1e-14)
        with self.assertRaisesRegex(ValueError,'missing rectangular'):
            geometric_trace(p,[3])
        with self.assertRaisesRegex(ValueError,'missing rectangular'):
            occupation_response(p,[3,4],1.,0.,3.)
        with self.assertRaisesRegex(ValueError,'missing rectangular'):
            occupation_response(p,[3,4],-100.,1.,-100.)  # Underflow does not certify equal weights.

    def test_transitive_full_source_cluster_and_complete_group(self):
        e=np.tile([0.,.0015,.003,2.],(4,1));p=optical_pairs(self.c,e)
        self.assertTrue(np.all(p.clusters[:,:3]==0))
        with self.assertRaisesRegex(ValueError,'producer-erased'):
            geometric_trace(p,[1,3])
        trace=geometric_trace(p,[1,2,3])
        value,*_=occupation_response(p,[1,2,3],1.,0.,1.)
        np.testing.assert_allclose(value,trace,atol=1e-14)
        with self.assertRaisesRegex(ValueError,'producer-erased'):
            occupation_response(p,[1,2,3],-100.,1.,-100.)
        with self.assertRaisesRegex(ValueError,'producer-erased'):
            occupation_response(p,[1,2,3],.002,0.,1.)

    def test_exact_degenerate_group_and_unitary_rotation(self):
        e=np.tile([-1.,-1.,1.,2.],(4,1));p=optical_pairs(self.c,e)
        value,*_=occupation_response(p,[1,2],-1.,300.,-2.)
        np.testing.assert_allclose(value,.5*geometric_trace(p,[1,2]),atol=1e-14)
        rng=np.random.default_rng(166)
        u=np.eye(4,dtype=complex);u[:2,:2]=np.linalg.qr(rng.normal(size=(2,2))+1j*rng.normal(size=(2,2)))[0]
        rotated=u.conj().T@self.c.astype(complex)@u
        rotated_value,*_=occupation_response(optical_pairs(rotated,e),[1,2],-1.,300.,-2.)
        np.testing.assert_allclose(rotated_value,value,atol=1e-13)

    def test_required_quality_only_and_gauge_invariance(self):
        c=self.c.copy();c[:,:,3,2]+=50  # Unselected high-high pair need not be globally valid.
        p=optical_pairs(c,self.e)
        occupation_response(p,[1],0.,0.,0.)
        with self.assertRaisesRegex(ValueError,'Hermitian consistency'):
            occupation_response(p,[3],1.,0.,0.)
        c[:,:,2,0]=np.nan
        with self.assertRaisesRegex(ValueError,'nonfinite required'):
            occupation_response(optical_pairs(c,self.e),[1],0.,0.,0.)
        phase=np.exp(1j*np.arange(4)*.731)
        gauged=self.c.astype(complex)*phase[None,None,None,:]*phase.conj()[None,None,:,None]
        p=optical_pairs(self.c,self.e);g=optical_pairs(gauged,self.e)
        np.testing.assert_allclose(geometric_trace(g,[1,3]),geometric_trace(p,[1,3]),atol=1e-13)

    def test_stable_fermi_tail_and_empty_complement(self):
        self.assertGreater(float(logistic_difference(np.array(-100.),np.array(-101.))),-1e-40)
        self.assertNotEqual(float(logistic_difference(np.array(-100.),np.array(-101.))),0.)
        np.testing.assert_array_equal(geometric_trace(optical_pairs(self.c,self.e),[1,2,3,4]),0.)


class SelectedWorkflowTests(unittest.TestCase):
    # Reuse input fixture construction, not its occupied-mode test cases.
    setUp_base=occupied_fixture.WaveDerHallTests.setUp
    fake_wavecar=occupied_fixture.WaveDerHallTests.fake_wavecar
    write_run=occupied_fixture.WaveDerHallTests.write_run
    make_chunks=occupied_fixture.WaveDerHallTests.make_chunks
    cli_arguments=occupied_fixture.WaveDerHallTests.cli_arguments

    def setUp(self):
        self.setUp_base()
        self.c=hermitian_fixture(nb=3)[None]
        self.energies[:]=[-1.,0.,2.]
        self.write_run()
        self.options.pop('occupied')
        self.options.update(bands=[2],temperatures=[0.,300.])

    def spectrum(self,**opts):
        return wh.waveder_hall_spectrum(self.run,[-.3,.3],**dict(self.options,**opts))

    def test_selected_mesh_sign_regions_and_metadata(self):
        rows,meta=self.spectrum(region_spec={'regions':[{'name':'left','k_ids':[1,2]},
                                                     {'name':'right','k_ids':[3,4]}]},
                                differences=[('contrast','left','right')])
        self.assertEqual(meta['scope'],'selected_band_contribution')
        self.assertEqual(meta['selected_band_ids'],[2]);self.assertEqual(meta['intermediate_band_ids'],[1,2,3])
        self.assertFalse(meta['total_AHC_certified']);self.assertFalse(meta['total_delta_AHC_certified'])
        for t in (0.,300.):
            for mu in (-.3,.3):
                r={x['region']:x for x in rows if x['temperature_K']==t and x['mu_eV']==mu}
                expected=-2*np.pi*np.mean(band_oracle(self.c[0],self.energies[0],[2],mu,t)[:,2])
                self.assertAlmostEqual(r['total']['sigma_e2_over_h'],expected,places=12)
                self.assertAlmostEqual(r['total']['sigma_e2_over_h'],r['left']['sigma_e2_over_h']+r['right']['sigma_e2_over_h'])
                self.assertAlmostEqual(r['contrast']['sigma_e2_over_h'],r['left']['sigma_e2_over_h']-r['right']['sigma_e2_over_h'])
        reverse,_=self.spectrum(sampling={'kind':'uniform_full_2d','mesh':[2,2],'plane_axes':[1,0]})
        forward,_=self.spectrum()
        np.testing.assert_allclose([r['sigma_e2_over_h'] for r in reverse],[-r['sigma_e2_over_h'] for r in forward])
        with self.assertRaises(ValueError): self.spectrum(sampling={'kind':'uniform_full_2d','mesh':[3,3],'plane_axes':[0,1]})

    def test_true_chunks_and_source_binding(self):
        single,_=self.spectrum()
        chunks=self.make_chunks([[0,1],[2,3]])
        rows,meta=wh.waveder_hall_spectrum(chunks,[-.3,.3],**self.options)
        self.assertEqual(rows,single);self.assertEqual(meta['source_metadata']['source_run_count'],2)
        (chunks[1]/'CHGCAR').write_bytes(b'different')
        with self.assertRaisesRegex(ValueError,'different CHGCAR'):
            wh.waveder_hall_spectrum(chunks,[-.3,.3],**self.options)

    def test_standalone_and_unified_cli(self):
        args=self.cli_arguments('standalone');i=args.index('--occupied');args[i:i+2]=['--bands','2']
        args+=['--temperatures','0','300']
        with contextlib.redirect_stdout(io.StringIO()): wh.main(args)
        m=json.loads((self.root/'standalone/conductivity.json').read_text())
        self.assertEqual(m['selected_band_ids'],[2])
        self.assertEqual(m['provenance']['vaspberry_version'],(ROOT/'VERSION').read_text().strip())
        for name in ('waveder_selected.py','band_selection.py'):
            self.assertEqual(m['provenance']['implementation_sha256'][name],wh.sha256(ROOT/'tools'/name))
        args[args.index(str(self.root/'standalone'))]=str(self.root/'unified')
        with contextlib.redirect_stdout(io.StringIO()): cli.main(['kubo-hall',*args])
        args[args.index(str(self.root/'unified'))]=str(self.root/'rejected')
        with contextlib.redirect_stderr(io.StringIO()),self.assertRaises(SystemExit):
            cli.main(['kubo-hall',*args,'--kubo-source','wavecar'])
        self.assertFalse((self.root/'rejected').exists())
        with contextlib.redirect_stderr(io.StringIO()),self.assertRaises(SystemExit):
            wh.main([*args,'--occupied','1'])

    def config(self,extra=''):
        p=self.root/'study.ini'
        p.write_text('[run]\ninput_dir=run\noutput=study-output\nmesh=2 2\nspin_mode=spinor\nenergy_reference=fixture zero\n'
                     '[hall]\nbands=2\nmu=-.3 .3 3\nreference=.25\ntemperatures=0 300\n'+extra)
        return p

    def test_realistic_system_title_semicolon_preserves_physical_assignments(self):
        path=self.run/'INCAR';old=path.read_text()
        path.write_text('SYSTEM=probe; no PAW exporter patch\n'+old.replace('SYSTEM = fixture\n',''))
        self.spectrum()
        path.write_text(old.replace('SYSTEM = fixture','SYSTEM=probe; title; LOPTICS=.FALSE.'))
        with self.assertRaises(ValueError): self.spectrum()
        self.assertEqual(wh.incar_values('SYSTEM=a; title; LOPTICS=F')['LOPTICS'],'F')
        for text in ('LOPTICS=T; invalid words','SYSTEM=a; malformed assignment = x',
                     'SYSTEM=a; title; LOPTICS=T; invalid words'):
            with self.subTest(text=text),self.assertRaisesRegex(ValueError,'ambiguous INCAR'):
                wh.incar_values(text)

    def test_ini_selected_route_and_invalid_options(self):
        settings=load_settings(self.config());self.assertEqual(settings['hall']['bands'],[2])
        inputs,cache=post.preflight(settings)
        self.assertEqual(cache['waveder_hall'][1]['scope'],'selected_band_contribution')
        self.assertEqual(set(inputs),{'wavecar','waveder','incar','outcar'})
        with self.assertRaisesRegex(ValueError,'mutually exclusive'):load_settings(self.config('occupied=1\n'))
        p=self.config();p.write_text(p.read_text().replace('input_dir=run','input_dir=run\nkubo_source=wavecar'))
        with self.assertRaisesRegex(ValueError,'only to kubo_source'):load_settings(p)


if __name__=='__main__': unittest.main()
