"""Reproduce the fixed Bi16 localization sequence using supplied Wannier90 3.1.0."""
from pathlib import Path
import argparse,hashlib,json,os,re,shutil,subprocess,time
from restore import restore,digest
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--wannier90',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--timeout',type=float,default=300);a=p.parse_args()
base=Path(__file__).resolve().parent;exe=a.wannier90.resolve()
if a.output_dir.exists():raise ValueError('Output exists; choose a new directory')
if a.timeout<=0:raise ValueError('Positive timeout required')
a.output_dir.mkdir(parents=True);restore('raw',a.output_dir/'raw')
m=json.loads((base/'manifest.json').read_text());assert digest(base/'wannier90.win')==m['plain_files']['wannier90.win']['sha256']
env=dict(os.environ)
for name in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):env[name]='1'
records=[];previous=None
for stage,count in enumerate([1000,5000,5000,5000]):
 case=a.output_dir/f'stage{stage:02d}';case.mkdir()
 for ext in ('amn','mmn','eig'):shutil.copyfile(a.output_dir/'raw'/f'wannier90.{ext}',case/f'wannier90.{ext}')
 text=re.sub(r'(?m)^num_iter\s*=.*$',f'num_iter = {count}',(base/'wannier90.win').read_text())
 if previous:
  shutil.copyfile(previous/'wannier90.chk',case/'wannier90.chk');text+='restart = wannierise\n'
 (case/'wannier90.win').write_text(text)
 record=dict(status='RUNNING',stage=stage,iteration_cap=count,command=[str(exe),'wannier90'],binary_sha256=digest(exe),input_sha256={f.name:digest(f) for f in case.iterdir() if f.is_file()});start=time.monotonic()
 try:
  with (case/'stdout.log').open('w') as out,(case/'stderr.log').open('w') as err:r=subprocess.run([str(exe),'wannier90'],cwd=case,env=env,stdout=out,stderr=err,timeout=a.timeout)
  record['returncode']=r.returncode
  wout=(case/'wannier90.wout').read_text(errors='replace') if (case/'wannier90.wout').exists() else ''
  if r.returncode or 'All done: wannier90 exiting' not in wout or '3.1.0' not in wout:raise RuntimeError('Wannier90 3.1.0 normal completion required; inspect retained logs')
  record.update(status='PASS',output_checkpoint_sha256=digest(case/'wannier90.chk'))
 except BaseException as exc:
  record.update(status='FAILED',error=str(exc));raise
 finally:
  record['wall_seconds']=time.monotonic()-start;(case/'run.json').write_text(json.dumps(record,indent=2)+'\n')
 records.append(record);previous=case
(a.output_dir/'result.json').write_text(json.dumps(dict(status='PASS_EXECUTION',stages=records,scope='Fixed iteration sequence reproduced. Numerical tolerances, raw-source band fit and actual spread history must still be checked; completion does not assert strict spread convergence.'),indent=2)+'\n')
print(previous)
