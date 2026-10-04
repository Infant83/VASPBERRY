"""A/B Intel Hydra SIGPIPE handling without modifying binaries or verdict rules.

Compile mpi_abort_signal_probe.f90 first with the same MPI Fortran wrapper.
This standalone diagnostic only uses Python's standard library. It preserves
all failing attempts and distinguishes a supported mitigation from an
inconclusive non-reproduction. No networking, publication or source editing.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import signal
import shutil
import subprocess
import time


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def signal_policy(command, variant):
    if variant == 'default':
        return command, True
    if variant == 'rank-files':
        if command[1:3] == ['-n', '2']:
            command = [*command[:3], '/bin/sh', '-c',
                       'exec "$@" >"mpi-process-$$.stdout.log" 2>"mpi-process-$$.stderr.log"',
                       'mpi-direct-log', *command[3:]]
        return command, True
    if variant == 'ignore-pipe':
        # A literal script and a separate argv array; no interpolated shell code.
        # exec preserves the actual launcher's exit/signal result. Only SIGPIPE
        # is ignored, leaving Python's default handling of other signals intact.
        return ['/bin/sh', '-c', 'trap "" PIPE; exec "$@"', 'mpi-sigpipe-wrapper', *command], True
    if variant == 'inherit-python':
        return command, False
    raise ValueError(variant)


def run_case(command, directory, variant, expected, diagnostic, timeout):
    directory.mkdir()
    argv, restore = signal_policy(command, variant)
    row = dict(command=argv, original_command=command, cwd=str(directory),
               variant=variant, restore_signals=restore, expected_exit=expected,
               required_diagnostic=diagnostic, status='RUNNING')
    start = time.monotonic()
    try:
        with (directory/'stdout.log').open('wb') as out, (directory/'stderr.log').open('wb') as err:
            process = subprocess.Popen(argv, cwd=directory, stdout=out, stderr=err,
                                       restore_signals=restore, start_new_session=True)
            try:
                exit_code = process.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid,signal.SIGKILL)
                process.wait()
                raise
        row['exit_code'] = exit_code
        text = ' '.join('\n'.join(p.read_text(errors='replace') for p in sorted(directory.glob('*.log'))).lower().split())
        row['diagnostic_present'] = not diagnostic or ' '.join(diagnostic.lower().split()) in text
        row['completed_csv'] = [p.name for p in directory.glob('*.csv*')
                                if p.read_text(errors='replace').strip().endswith('# result_status=PASS')]
        row['any_csv_files'] = [p.name for p in directory.glob('*.csv*')]
        row['status'] = ('PASS' if exit_code == expected and row['diagnostic_present']
                         and not row['any_csv_files'] else 'FAIL')
    except subprocess.TimeoutExpired:
        row.update(status='TIMEOUT', exit_code=None)
    except Exception as error:
        row.update(status='ERROR', error=f'{type(error).__name__}: {error}')
    finally:
        row['elapsed_seconds'] = time.monotonic()-start
        row['logs_sha256'] = {p.name:sha(p) for p in directory.glob('*.log')}
        (directory/'command.json').write_text(json.dumps(row,indent=2)+'\n')
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--launcher','--mpiexec',dest='launcher',required=True)
    parser.add_argument('--probe', type=Path)
    parser.add_argument('--compiler', help='MPI compiler wrapper, required unless --probe already built')
    parser.add_argument('--product','--native',dest='product',type=Path,required=True)
    parser.add_argument('--repo', type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--repeats', type=int, default=20)
    parser.add_argument('--timeout', type=float, default=12)
    parser.add_argument('--variants', nargs='+', choices=['default','ignore-pipe','inherit-python','rank-files'],
                        default=['default','ignore-pipe'])
    args = parser.parse_args()
    if args.repeats < 1 or args.timeout <= 0:
        parser.error('positive repeats/timeout required')
    launcher = Path(shutil.which(args.launcher) or args.launcher).absolute()
    if not launcher.is_file():
        parser.error('MPI launcher does not exist')
    product,repo = args.product.resolve(strict=True),args.repo.resolve(strict=True)
    output = args.output_dir.resolve()
    output.mkdir(parents=True,exist_ok=False)
    if args.probe:
        probe=args.probe.resolve(strict=True)
    else:
        if not args.compiler:
            parser.error('--compiler is required without --probe')
        compiler=Path(shutil.which(args.compiler) or args.compiler).absolute()
        if not compiler.is_file():
            parser.error('MPI compiler does not exist')
        probe=output/'mpi-abort-probe'
        source=Path(__file__).with_name('mpi_abort_signal_probe.f90')
        command=[str(compiler),'-O0','-o',str(probe),str(source)]
        with (output/'build.stdout.log').open('wb') as out,(output/'build.stderr.log').open('wb') as err:
            built=subprocess.run(command,stdout=out,stderr=err,timeout=90)
        (output/'build.json').write_text(json.dumps(dict(command=command,exit_code=built.returncode,
            source_sha256=sha(source)),indent=2)+'\n')
        if built.returncode:
            raise RuntimeError('probe compilation failed; see build.stderr.log')
    spec = importlib.util.spec_from_file_location('original_intel_validation',repo/'tests/run_intel_mpi_validation.py')
    original = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(original)
    inputs = output/'inputs'; inputs.mkdir()
    wave = inputs/'WAVECAR'
    truncated = inputs/'TRUNCATED'
    original.write_synthetic_wavecar(wave)
    original.write_synthetic_wavecar(truncated, corruption='truncated')
    prefix = [str(launcher),'-n','2']
    tasks = [
        ('all', [*prefix,str(probe),'all'],1,'PROBE_ABORT'),
        ('rank0', [*prefix,str(probe),'rank0'],1,'PROBE_ABORT rank=0 code=1'),
        ('rank1', [*prefix,str(probe),'rank1'],1,'PROBE_ABORT rank=1 code=1'),
        ('directory', [*prefix,str(product),'--input-dir',str(wave),'--wavecar',str(wave),
                       '--spinor','1','--task','kubo','--kubo-source','wavecar','--bands','1:2'],
                      1,'must be an existing directory'),
        ('truncated-per-band', [*prefix,str(product),'--wavecar',str(truncated),'--spinor','1',
                       '--task','kubo','--kubo-source','wavecar','--bands','1:2','--per-band','1',
                       '--curvature-csv','BAND.csv'],1,'invalid/truncated WAVECAR layout'),
    ]
    summary = dict(started_utc=datetime.now(timezone.utc).isoformat(), status='RUNNING',
                   input_sha256={str(p):sha(p) for p in (probe,product,wave,truncated)},
                   launcher=str(launcher), python_parent_sigpipe=str(signal.getsignal(signal.SIGPIPE)),
                   environment={key:os.environ.get(key) for key in ('I_MPI_FABRICS','I_MPI_PIN')},
                   variants={}, limitations='A diagnostic pass is not the complete release CI gate.')
    try:
        for variant in args.variants:
            directory=output/variant; directory.mkdir()
            rows=[]
            controls = [
                ('signal-policy', ['/bin/sh','-c','kill -PIPE $$; exit 17'],
                 -signal.SIGPIPE if variant in ('default','rank-files') else 17,None),
                ('mpi-success', [*prefix,str(probe),'success'],0,'PROBE_SUCCESS'),
                ('mpi-abort7', [*prefix,str(probe),'abort7'],7,'code=7'),
            ]
            for name,command,expected,diagnostic in controls:
                row=run_case(command,directory/name,variant,expected,diagnostic,args.timeout)
                row['case']=name; rows.append(row)
            for iteration in range(args.repeats):
                for name,command,expected,diagnostic in tasks:
                    row=run_case(command,directory/f'{name}-{iteration:03d}',variant,expected,diagnostic,args.timeout)
                    row['case']=name; rows.append(row)
                print(json.dumps({'variant':variant,'iteration':iteration,'failures':sum(r['status']!='PASS' for r in rows)}),flush=True)
            summary['variants'][variant] = dict(rows=rows,
                status='PASS' if all(r['status']=='PASS' for r in rows) else 'FAIL',
                exit_counts=dict(Counter(str(r.get('exit_code')) for r in rows)),
                failed_cases=[r['case'] for r in rows if r['status']!='PASS'])
            (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
        baseline = summary['variants'].get('default',{})
        baseline_sigpipe = any(r.get('exit_code') == -signal.SIGPIPE and r['case']!='signal-policy'
                               for r in baseline.get('rows',[]))
        candidates = [v for name,v in summary['variants'].items() if name!='default']
        good = bool(candidates) and all(v['status']=='PASS' for v in candidates)
        summary['status'] = ('REJECT_CANDIDATE' if not good else
                             'SUPPORTED_SIGNAL_POLICY_MITIGATION' if baseline_sigpipe else
                             'CANDIDATE_PASS_BASELINE_NOT_REPRODUCED')
        summary['baseline_mpi_sigpipe_reproduced']=baseline_sigpipe
    finally:
        summary['completed_utc']=datetime.now(timezone.utc).isoformat()
        (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({'status':summary['status']}))
    return 1 if summary['status']=='REJECT_CANDIDATE' else 0


if __name__=='__main__':
    raise SystemExit(main())
