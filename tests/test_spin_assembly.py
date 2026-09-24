"""Independent record, mesh, and provenance tests for spin chunk assembly.

Only audit_producer is mocked: synthetic binary WAVECAR records still pass the
real reader and the assembled physical bundle passes the real matrix contract.
These zero-spin fixtures do not represent a VASP material or PAW calculation.
"""
from pathlib import Path
import json
import struct
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
import spin_assembly as assembly
from exported_matrix_kubo import sha256
from spin_hall_workflow import read_matrices
from wavecar_fukui import Wavecar

Q=np.array([[0.,0.,0.],[0.,.5,0.],[.5,0.,0.],[.5,.5,0.]])


def write_wavecar(path,q,energies,occupations,stride,seed,fermi):
    """Byte-RECL45200 fixture with nonzero unused bytes to detect record loss."""
    nk,nb=energies.shape
    raw=bytearray([seed]*(stride*(2+nk*(nb+1))))
    struct.pack_into('<3d',raw,0,stride,1,45200)
    struct.pack_into('<13d',raw,stride,nk,nb,120.,*np.eye(3).ravel(),fermi)
    for ik,point in enumerate(q):
        # Independent G enumeration for a cubic1A cell with120eV cutoff.
        count=sum(np.dot(2*np.pi*(point+g),2*np.pi*(point+g))/.262465831 < 120.
            for g in np.array(np.meshgrid(range(-3,4),range(-3,4),range(-3,4))).T.reshape(-1,3))
        nplane=2*count
        assert nplane*8<=stride
        record=2+ik*(nb+1)
        values=[nplane,*point]
        for e,f in zip(energies[ik],occupations[ik]):values.extend([e,0.,f])
        struct.pack_into('<'+str(len(values))+'d',raw,record*stride,*values)
        for band in range(nb):
            c=(np.arange(nplane)+1j*(seed+band+ik+1)).astype('<c8')
            offset=(record+1+band)*stride;raw[offset:offset+c.nbytes]=c.tobytes()
    path.write_bytes(raw)


class SpinAssemblyTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.path=Path(self.tmp.name);self.data={};self.reports={}
        self.cases=self.make_cases('two',[[3,1],[2,0]])

    def make_cases(self,name,groups,nb=4):
        cases=[]
        for index,group in enumerate(groups):
            case=self.path/name/str(index);case.mkdir(parents=True)
            q=Q[group].copy();nk=len(q)
            e=np.tile([-3.,-1.,1.,4.],(nk,1))[:,:nb]+.1*q.sum(axis=1)[:,None]
            occ=np.zeros((nk,nb));occ[:,:2]=1
            stride=256+256*index
            write_wavecar(case/'WAVECAR',q,e,occ,stride,11+index,.25+index)
            for n in ['CHGCAR','POTCAR','POSCAR']:(case/n).write_text('Synthetic same-Hamiltonian test input '+n+'\n')
            (case/'KPOINTS').write_text(json.dumps(q.tolist()))
            (case/'INCAR').write_text('ICHARG=11\nENCUT=120\nLSORBIT=T\nSAXIS=0 0 1\nNELM='+str(30+index)+'\n')
            (case/'OUTCAR').write_text('Synthetic independent assembly fixture; producer validation mocked\n')
            (case/'SPIN_VELOCITY.bin').write_bytes(b'fixture'+bytes(16384*nk))
            pseudo=np.zeros((nk,4,nb,nb),complex);pseudo[:,0]=np.eye(nb)
            velocity=np.zeros((nk,3,nb,nb),complex)
            velocity[:,0,0,2]=q[:,0]+1j*q[:,1];velocity[:,0,2,0]=velocity[:,0,0,2].conj()
            d=dict(kpoints_fractional=q,energies_eV=e,occupations=occ,weights=np.full(nk,1/nk),
                lattice_A=np.eye(3),pseudo=pseudo,delta=np.zeros_like(pseudo),velocity_eVA=velocity)
            self.data[case.resolve()]=d
            self.refresh(case);cases.append(case)
        return cases

    def refresh(self,case):
        source={n:sha256(case/n) for n in ['WAVECAR','INCAR','POSCAR','KPOINTS','OUTCAR','SPIN_VELOCITY.bin']}
        run=dict(status='FINISHED',input_sha256={n:sha256(case/n) for n in ['CHGCAR','POTCAR','POSCAR','KPOINTS']})
        (case/'run.json').write_text(json.dumps(run))
        source['run.json']=sha256(case/'run.json')
        self.reports[case.resolve()]=dict(status='PASS',source_sha256=source,producer_binary_sha256='a'*64,
            producer='Synthetic assembly fixture, not a material calculation')

    def audit(self,case):
        return self.data[case.resolve()],self.reports[case.resolve()]

    def merge(self,name='result',cases=None,**kwargs):
        with patch.object(assembly,'audit_producer',side_effect=self.audit):
            return assembly.merge_runs(cases or self.cases,self.path/name,[2,2],2,**kwargs)

    def test_two_and_four_chunks_reorder_every_raw_record_and_preserve_gauge(self):
        for name,cases in [('two-result',self.cases),('four-result',self.make_cases('four',[[3],[1],[2],[0]]))]:
            with self.subTest(chunks=len(cases)):
                meta=self.merge(name,cases);target=self.path/name
                wave=Wavecar(target/'WAVECAR',spinor_components=2)
                np.testing.assert_array_equal(wave.kpoints,Q)
                source_waves=[Wavecar(c/'WAVECAR',spinor_components=2) for c in cases]
                self.assertEqual(wave.header.stride_bytes,max(w.header.stride_bytes for w in source_waves))
                actual=(target/'WAVECAR').read_bytes();stride=wave.header.stride_bytes
                first=(cases[0]/'WAVECAR').read_bytes();oldstride=source_waves[0].header.stride_bytes
                self.assertEqual(actual[8:oldstride],first[8:oldstride])
                self.assertEqual(actual[stride+8:stride+oldstride],first[oldstride+8:2*oldstride])
                self.assertEqual(struct.unpack_from('<d',actual,stride+96)[0],.25)
                for ik,q in enumerate(Q):
                    for case,old in zip(cases,source_waves):
                        match=np.flatnonzero(np.all(old.kpoints==q,axis=1))
                        if not len(match):continue
                        oldk=int(match[0]);raw=(case/'WAVECAR').read_bytes();s=old.header.stride_bytes
                        for band in range(5):
                            start=(wave.record_number(ik)+band-1)*stride
                            previous=(old.record_number(oldk)+band-1)*s
                            self.assertEqual(actual[start:start+s],raw[previous:previous+s])
                            self.assertEqual(actual[start+s:start+stride],bytes(stride-s))
                        np.testing.assert_array_equal(wave.coefficients(ik,[1,2,3,4]),old.coefficients(oldk,[1,2,3,4]))
                arrays,m,_=read_matrices(target/'physical-matrices.npz',target/'physical-matrices.json')
                np.testing.assert_array_equal(arrays['kpoints_fractional'],Q)
                np.testing.assert_array_equal(arrays['weights'],np.full(4,.25))
                np.testing.assert_array_equal(arrays['spin_pauli'],np.zeros((4,3,4,4)))
                np.testing.assert_array_equal(arrays['velocity_eVA'][:,0,0,2],Q[:,0]+1j*Q[:,1])
                self.assertEqual(m['source_wavecar_sha256'],sha256(target/'WAVECAR'))
                self.assertEqual(m['source_sha256']['WAVECAR'],sha256(target/'WAVECAR'))
                self.assertNotIn(m['source_wavecar_sha256'],[sha256(c/'WAVECAR') for c in cases])
                self.assertEqual(meta['assembly']['source_chunk_point_counts'],[len(self.data[c.resolve()]['weights']) for c in cases])
                self.assertEqual(json.loads((target/'run.json').read_text())['status'],'PASS')
                self.assertFalse((self.path/(name+'.partial')).exists())

    def test_missing_or_duplicate_grid_reject_before_output(self):
        for name,groups in [('missing',[[0],[1,2]]),('duplicate',[[0,1],[2,2]])]:
            with self.subTest(name=name):
                cases=self.make_cases(name,groups)
                with self.assertRaises(ValueError):self.merge(name+'-out',cases)
                self.assertFalse((self.path/(name+'-out')).exists());self.assertFalse((self.path/(name+'-out.partial')).exists())

    def test_same_hamiltonian_and_input_association_reject(self):
        for file in ['CHGCAR','POTCAR','POSCAR','INCAR']:
            cases=self.make_cases('mismatch-'+file,[[0,1],[2,3]])
            original=(cases[1]/file).read_text()
            (cases[1]/file).write_text(original.replace('ENCUT=120','ENCUT=121') if file=='INCAR' else original+'changed\n')
            self.refresh(cases[1])
            with self.subTest(file=file),self.assertRaisesRegex(ValueError,'same fixed Hamiltonian'):
                self.merge('bad-'+file,cases)
        cases=self.make_cases('association',[[0,1],[2,3]])
        (cases[0]/'CHGCAR').write_text('changed after producer')
        with self.assertRaisesRegex(ValueError,'source association'):self.merge('bad-association',cases)

    def test_band_coverage_producer_and_occupation_guards(self):
        cases=self.make_cases('bands',[[0,1],[2,3]])
        d=self.data[cases[1].resolve()]
        write_wavecar(cases[1]/'WAVECAR',d['kpoints_fractional'],d['energies_eV'][:,:3],d['occupations'][:,:3],512,12,1.25)
        self.refresh(cases[1])
        with self.assertRaisesRegex(ValueError,'band coverage'):self.merge('bad-bands',cases)
        cases=self.make_cases('producer',[[0,1],[2,3]]);self.reports[cases[1].resolve()]['producer_binary_sha256']='b'*64
        with self.assertRaisesRegex(ValueError,'different producer'):self.merge('bad-producer',cases)
        self.data[self.cases[0].resolve()]['occupations'][0,1]=.5
        with self.assertRaisesRegex(ValueError,'insulating occupied'):self.merge('bad-occupations')

    def test_memory_and_duplicate_directory_preflight_never_calls_producer(self):
        with patch.object(assembly,'audit_producer') as audit:
            for limit in [.000001,float('nan'),-1]:
                with self.assertRaisesRegex(ValueError,'memory'):assembly.merge_runs(self.cases,self.path/'memory',[2,2],2,memory_limit_mib=limit)
            with self.assertRaisesRegex(ValueError,'distinct'):assembly.merge_runs([self.cases[0]]*2,self.path/'duplicate-case',[2,2],2)
            audit.assert_not_called()
        self.assertFalse((self.path/'memory').exists())

    def test_real_cli_memory_preflight_rejects_before_mock_producer_data(self):
        target=self.path/'cli-memory'
        command=[sys.executable,str(ROOT/'tools/vaspberry_kubo.py'),'spin-merge',
            '--run-dirs',*[str(c) for c in self.cases],'--mesh','2','2','--occupied','2',
            '--memory-limit-mib','0.000001','--output-dir',str(target)]
        result=subprocess.run(command,capture_output=True,text=True,timeout=30)
        self.assertNotEqual(result.returncode,0)
        self.assertIn('memory estimate exceeds limit',result.stderr)
        self.assertFalse(target.exists());self.assertFalse(target.with_name(target.name+'.partial').exists())

    def test_output_guard_and_failed_partial_keep_source_unchanged(self):
        before={str(c):sha256(c/'WAVECAR') for c in self.cases}
        self.merge('existing')
        with self.assertRaisesRegex(ValueError,'exists'):self.merge('existing')
        with patch.object(assembly,'assemble_wavecar',side_effect=ValueError('deliberate copy failure')):
            with self.assertRaisesRegex(ValueError,'copy failure'):self.merge('failed')
        self.assertFalse((self.path/'failed').exists())
        self.assertEqual(json.loads((self.path/'failed.partial/run.json').read_text())['status'],'FAILED')
        self.assertEqual(before,{str(c):sha256(c/'WAVECAR') for c in self.cases})

    def test_legacy_word_recl_and_coefficient_record_overflow_rejected(self):
        original=(self.cases[0]/'WAVECAR').read_bytes()
        for name,replacement in [('word-recl',None),('overflow',128.)]:
            raw=bytearray(original)
            if replacement is None:struct.pack_into('<d',raw,0,64.)
            else:struct.pack_into('<d',raw,2*256,replacement)
            (self.cases[0]/'WAVECAR').write_bytes(raw)
            waves=[Wavecar(c/'WAVECAR',spinor_components=2) for c in self.cases]
            target=self.path/(name+'-WAVECAR')
            with self.subTest(name=name),self.assertRaises(ValueError):
                assembly.assemble_wavecar([c/'WAVECAR' for c in self.cases],waves,[3,1,2,0],target)
            self.assertFalse(target.exists())
        (self.cases[0]/'WAVECAR').write_bytes(original)

    def test_source_change_during_assembly_prevents_final_publication(self):
        original=assembly.assemble_wavecar
        def change(*args,**kwargs):
            wave=original(*args,**kwargs);(self.cases[0]/'INCAR').write_text('changed after copy');return wave
        with patch.object(assembly,'assemble_wavecar',side_effect=change):
            with self.assertRaisesRegex(ValueError,'changed during assembly'):self.merge('changed')
        self.assertFalse((self.path/'changed').exists())
        self.assertEqual(json.loads((self.path/'changed.partial/run.json').read_text())['status'],'FAILED')


if __name__=='__main__':unittest.main()
