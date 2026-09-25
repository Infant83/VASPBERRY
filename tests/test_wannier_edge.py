"""Analytic finite-chain and topological-edge oracles, plus output guards."""
from pathlib import Path
from types import SimpleNamespace
import json
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tools'))
import wannier_edge as edge
from wannier_operators import WannierOperators, SCHEMA, UNITS, PHASE, write_operators

SX=np.array([[0,1],[1,0]],complex)
SY=np.array([[0,-1j],[1j,0]],complex)
SZ=np.diag([1.,-1.]).astype(complex)


def model(records):
    keys=sorted(records); n=next(iter(records.values())).shape[-1]
    metadata=dict(schema=SCHEMA,version=1,complete=True,source_format='wannier90_effective_HH_R_AA_R',
        units=UNITS,fourier_phase=PHASE,real_space_degeneracy='already_absorbed',
        position_components=['x','y','z'],energy_reference='analytic test only',
        source_operator={'kind':'Hamiltonian_and_position_connection','scope':'represented finite Wannier subspace'},
        source_sha256={'HH_R':'0'*64,'AA_R':'0'*64},spinor_components=1,spin_multiplicity=1,
        num_wann=n,nrpts=len(keys))
    return WannierOperators(np.array(keys,dtype=np.int64),np.array([records[k] for k in keys],complex),
        np.zeros((len(keys),3,n,n),complex),np.eye(3,dtype=float),metadata)


def qwz(mass=-1.2):
    return model({(0,0,0):mass*SZ,(1,0,0):(SZ-1j*SX)/2,(-1,0,0):(SZ+1j*SX)/2,
                  (0,1,0):(SZ-1j*SY)/2,(0,-1,0):(SZ+1j*SY)/2})


class WannierEdgeTests(unittest.TestCase):
    def test_exact_finite_chain_spectrum_and_edge_sum_rules(self):
        onsite=np.diag([-.7,1.3]); hopping=.2*np.eye(2); along=.15*np.eye(2)
        data=model({(0,0,0):onsite,(0,1,0):hopping,(0,-1,0):hopping,
                    (1,0,0):along,(-1,0,0):along})
        q=np.array([-.31,0.,.27]);width=9
        result=edge.strip_spectrum(data,q,width,edge_cells=2)
        for k,point in enumerate(q):
            modes=.4*np.cos(np.arange(1,width+1)*np.pi/(width+1))
            expected=np.sort((np.diag(onsite)[:,None]+modes+.3*np.cos(2*np.pi*point)).ravel())
            np.testing.assert_allclose(result['energies_eV'][k],expected,atol=3e-15)
        np.testing.assert_allclose(result['left_edge_weight'].sum(axis=1),4.,atol=1e-13)
        np.testing.assert_allclose(result['right_edge_weight'].sum(axis=1),4.,atol=1e-13)

    def test_chern_edge_branches_have_opposite_velocities_and_boundaries(self):
        q=np.array([-.025,.025]);out=edge.strip_spectrum(qwz(),q,24,3)
        for i,point in enumerate(q):
            candidates=np.argsort(abs(out['energies_eV'][i]))[:2]
            values=out['energies_eV'][i,candidates]
            np.testing.assert_allclose(np.sort(values),[-abs(np.sin(2*np.pi*point)),abs(np.sin(2*np.pi*point))],atol=1e-10)
            l=candidates[np.argmax(out['left_edge_weight'][i,candidates])]
            r=candidates[np.argmax(out['right_edge_weight'][i,candidates])]
            self.assertGreater(out['left_edge_weight'][i,l],.999)
            self.assertGreater(out['right_edge_weight'][i,r],.999)
            self.assertLess(out['energies_eV'][i,l]*out['energies_eV'][i,r],0)
        left_values=[out['energies_eV'][i,np.argmax(out['left_edge_weight'][i])] for i in range(2)]
        self.assertLess(left_values[0]*left_values[1],0)

    def test_trivial_insulator_has_no_midgap_edge_and_width_gap_decays(self):
        trivial=edge.strip_spectrum(qwz(-3),[0.],20)
        self.assertGreater(np.min(abs(trivial['energies_eV'])),1.)
        narrow=edge.strip_spectrum(qwz(),[0.],4,1)
        wide=edge.strip_spectrum(qwz(),[0.],12,1)
        self.assertLess(np.min(abs(wide['energies_eV'])),np.min(abs(narrow['energies_eV']))*1e-4)

    def test_orbital_rotation_preserves_spectrum_and_edge_projectors(self):
        data=qwz();u=np.array([[1,1j],[1j,1]],complex)/np.sqrt(2)
        rotated=model({tuple(r):u.conj().T@h@u for r,h in zip(data.irvec,data.hamiltonian_eV)})
        a=edge.strip_spectrum(data,[.13,.27],10)
        b=edge.strip_spectrum(rotated,[.13,.27],10)
        for key in ['energies_eV','left_edge_weight','right_edge_weight']:
            np.testing.assert_allclose(a[key],b[key],atol=1e-11)
        np.testing.assert_allclose(edge.strip_hamiltonian(data,.13,10),edge.strip_hamiltonian(data,1.13,10),atol=2e-15)

    def test_nonhermitian_invalid_regions_and_resource_bounds(self):
        data=qwz()
        for kwargs in [{'edge_cells':3},{'memory_limit_mib':.001},{'time_limit':-1}]:
            with self.assertRaises(ValueError):edge.strip_spectrum(data,[0.],4,**kwargs)
        with self.assertRaises(ValueError):edge.strip_hamiltonian(data,0.,4,0,0)
        data.hamiltonian_eV[0,0,1]+=1j
        with self.assertRaises(ValueError):edge.strip_spectrum(data,[0.],4)

    def test_command_formats_integrity_and_failed_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);cache=root/'cache';write_operators(cache,qwz())
            args=SimpleNamespace(operators=cache,output_dir=root/'edge',width=8,edge_cells=2,
                periodic_axis=0,open_axis=1,q_min=-.5,q_max=.5,kpoints=5,memory_limit_mib=64.,
                time_limit=30.,formats=['csv','dat','npz'])
            meta=edge.command(args);self.assertEqual(meta['schema'],'vaspberry.wannier-edge')
            csv=np.loadtxt(args.output_dir/'edge.csv',delimiter=',',skiprows=1)
            dat=np.loadtxt(args.output_dir/'edge.dat',skiprows=1)
            np.testing.assert_array_equal(csv,dat)
            with np.load(args.output_dir/'edge.npz') as z:
                np.testing.assert_array_equal(csv[:,3],z['energies_eV'].ravel())
            with self.assertRaises(ValueError):edge.command(args)
            args.output_dir=root/'failed'
            with patch.object(edge.np.linalg,'eigh',side_effect=ValueError('independent failure fixture')):
                with self.assertRaises(ValueError):edge.command(args)
            self.assertFalse(args.output_dir.exists())
            self.assertEqual(json.loads((root/'failed.partial/run.json').read_text())['status'],'FAILED')


if __name__=='__main__':unittest.main()
