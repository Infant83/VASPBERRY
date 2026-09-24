"""Restore checked generated Bi Wannier matrices into a new directory."""
from pathlib import Path
import argparse,gzip,hashlib,json

def digest(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for chunk in iter(lambda:f.read(1024**2),b''):h.update(chunk)
 return h.hexdigest()
def restore(kind,output):
 base=Path(__file__).resolve().parent;m=json.loads((base/'manifest.json').read_text())
 assert m['schema']=='vaspberry.bi16-wannier-reference' and m['complete'] is True
 if output.exists():raise ValueError('Output exists; choose a new directory')
 selected=[r for r in m['files'] if r['file'].startswith(kind+'/')]
 expected={'raw':{'wannier90.amn','wannier90.mmn','wannier90.eig'},'operators':{'wannier90_HH_R.dat','wannier90_AA_R.dat'}}[kind]
 assert {r['restored_name'] for r in selected}==expected
 for r in selected:
  assert r['file']==kind+'/'+r['restored_name']+'.gz';source=base/r['file'];assert source.stat().st_size==r['bytes'] and digest(source)==r['sha256']
 output.mkdir(parents=True);log={'status':'RESTORING','kind':kind,'outputs':{}}
 try:
  for r in selected:
   target=output/r['restored_name'];size=0
   with gzip.open(base/r['file'],'rb') as f,target.open('xb') as g:
    while chunk:=f.read(1024**2):
     size+=len(chunk)
     if size>r['uncompressed_bytes']:raise ValueError('Unexpected decompressed size')
     g.write(chunk)
   assert size==r['uncompressed_bytes'] and digest(target)==r['uncompressed_sha256']
   log['outputs'][target.name]=digest(target)
  log['status']='PASS'
 except BaseException as exc:
  log.update(status='FAILED',error=str(exc));raise
 finally:(output/'restore.json').write_text(json.dumps(log,indent=2)+'\n')
 return log
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--kind',choices=['raw','operators'],required=True);p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args();print(json.dumps(restore(a.kind,a.output_dir),indent=2))
