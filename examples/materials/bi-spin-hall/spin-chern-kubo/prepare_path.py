#!/usr/bin/env python3
"""Prepare a same-density ordinary-VASP Bi Gamma-K-M-Gamma SOC path."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--potcar',type=Path,required=True,help='Locally licensed potential matching the parent Bi example')
 p.add_argument('--charge',type=Path,help='Own matching SCF charge; default is the checksum-verified supplied reference density')
 p.add_argument('--nbands',type=int,default=48)
 p.add_argument('--points-per-segment',type=int,default=16,help='Intervals in each of the three path segments; default gives49points')
 p.add_argument('--output-dir',type=Path,required=True)
 a=p.parse_args()
 if a.points_per_segment<2:p.error('--points-per-segment must be at least2')
 if a.nbands<=10:p.error('--nbands must exceed the10occupied spinor bands')
 command=[sys.executable,str(HERE.parent/'prepare_vasp.py'),'--stage','wavecar','--potcar',str(a.potcar),'--mesh','6','6','--nbands',str(a.nbands),'--output-dir',str(a.output_dir)]
 if a.charge:command+=['--charge',str(a.charge)]
 subprocess.run(command,check=True)
 original=json.loads((a.output_dir/'input_manifest.json').read_text())
 vertices=[(0.,0.,0.),(1/3,2/3,0.),(0.,.5,0.),(0.,0.,0.)]
 points=[];n=a.points_per_segment
 for first,last in zip(vertices[:-1],vertices[1:]):
  points.extend(tuple((1-i/n)*x+i/n*y for x,y in zip(first,last)) for i in range(n))
 points.append(vertices[-1])
 kpoints='Bi SOC Gamma-K-M-Gamma path; exact same fixed density as full mesh\n'+str(len(points))+'\nReciprocal\n'
 kpoints+=''.join(' '.join(f'{x:.16f}' for x in q)+' 1.0\n' for q in points)
 (a.output_dir/'KPOINTS').write_text(kpoints)
 record=dict(status='PREPARED',stage='path',parent_preparation=original,path_labels=['Gamma','K','M','Gamma'],path_node_indices=[1,n+1,2*n+1,3*n+1],path_nodes_fractional=vertices,npoints=len(points),source_bands=a.nbands,occupied_bands=10,physics='Same parent Bi fixed-charge INCAR/POSCAR/POTCAR/CHGCAR; only KPOINTS replaced with explicit path',inputs_sha256={name:sha(a.output_dir/name) for name in ['INCAR','POSCAR','POTCAR','CHGCAR','KPOINTS']},preparation_sha256=sha(Path(__file__)))
 (a.output_dir/'input_manifest.json').write_text(json.dumps(record,indent=2)+'\n')
 print('Prepared path '+str(a.output_dir))
if __name__=='__main__':main()
