"""Run one prepared, instrumented VASP case with bounded resources and records.

The VASP executable and PAW datasets remain locally licensed inputs. This
wrapper does not prepare a material, copy potentials, or overwrite prior logs.
It always uses one scientific thread; parallelize independent k-point chunks
only within the available CPU and memory budget.
"""
from __future__ import annotations
import argparse
import datetime
import json
import os
from pathlib import Path
import resource
import signal
import subprocess
import time

from berry_data import require
from exported_matrix_kubo import sha256
from vasp_optical_export import FOOTER as OPTICAL_FOOTER
from vasp_spin_export import FOOTER as SPIN_FOOTER


def _complete(path, footer):
    if not path.is_file() or path.stat().st_size < len(footer):
        return False
    with path.open('rb') as stream:
        stream.seek(-len(footer), 2)
        return stream.read() == footer


def _run_case(directory, binary, producer_manifest, *, timeout, memory_limit_mib, state):
    case, binary, producer_manifest = Path(directory).resolve(), Path(binary).resolve(), Path(producer_manifest).resolve()
    require(case.is_dir() and binary.is_file() and producer_manifest.is_file(), 'run directory, binary and producer manifest must exist')
    require(0 < timeout <= 86400 and 0 < memory_limit_mib <= 1048576, 'positive bounded timeout and memory limit required')
    reserved = ('run.json', 'stdout.log', 'stderr.log', 'OUTCAR', 'SPIN_VELOCITY.bin', 'BERRY_CONNECTION.bin')
    require(not any((case/n).exists() for n in reserved), 'run directory already contains output; prepare a new directory')
    inputs = ['INCAR', 'POSCAR', 'KPOINTS', 'POTCAR', 'CHGCAR']
    require(all((case/n).is_file() for n in inputs), 'prepared fixed-charge VASP inputs including CHGCAR required')
    if (case/'WAVECAR').exists():
        inputs.append('WAVECAR')
    manifest = json.loads(producer_manifest.read_text())
    require(manifest.get('schema') == 'vaspberry.vasp544-spin-producer' and manifest.get('version') == 1,
            'compatible instrumentation manifest required')
    soft, hard = resource.getrlimit(resource.RLIMIT_STACK)
    desired = min(hard, 64*1024**2) if hard != resource.RLIM_INFINITY else 64*1024**2
    try:
        if soft != resource.RLIM_INFINITY:
            resource.setrlimit(resource.RLIMIT_STACK, (max(soft, desired), hard))
    except (ValueError, OSError):
        # Some hosts advertise an unlimited hard limit but enforce a lower cap.
        resource.setrlimit(resource.RLIMIT_STACK, (soft, hard))
    now = lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
    record = {'status': 'RUNNING', 'started_utc': now(), 'command': [str(binary)],
              'binary_sha256': sha256(binary), 'producer_manifest_sha256': sha256(producer_manifest),
              'wrapper_sha256': sha256(Path(__file__)),
              'stack_limit_bytes': resource.getrlimit(resource.RLIMIT_STACK)[0],
              'producer_manifest': manifest, 'threads': 1, 'timeout_seconds': timeout,
              'memory_limit_bytes': int(memory_limit_mib*1024**2),
              'input_sha256': {n:sha256(case/n) for n in inputs}}

    def save():
        temporary = case/'run.json.tmp'
        temporary.write_text(json.dumps(record, indent=2, allow_nan=False)+'\n')
        temporary.replace(case/'run.json')

    env = dict(os.environ)
    for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        env[name] = '1'
    save(); begin = time.monotonic(); peak = 0; reason = None
    with (case/'stdout.log').open('x') as stdout, (case/'stderr.log').open('x') as stderr:
        process = subprocess.Popen([str(binary)], cwd=case, env=env, stdout=stdout, stderr=stderr, start_new_session=True)
        state['process'] = process
        record['pid'] = process.pid; save()
        while process.poll() is None:
            elapsed = time.monotonic()-begin
            rss = subprocess.run(['ps', '-o', 'rss=', '-p', str(process.pid)], capture_output=True, text=True)
            value = rss.stdout.strip()
            if value.isdigit():
                peak = max(peak, int(value)*1024)
            if elapsed > timeout:
                reason = 'TIMEOUT'
            elif peak > record['memory_limit_bytes']:
                reason = 'MEMORY_LIMIT'
            if reason:
                os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL); process.wait()
                break
            record.update(elapsed_s=elapsed, sampled_peak_rss_bytes=peak); save()
            time.sleep(1)
    text = (case/'OUTCAR').read_text(errors='replace') if (case/'OUTCAR').is_file() else ''
    normal = 'General timing and accounting' in text
    converged = 'EDIFF is reached' in text
    exports = _complete(case/'SPIN_VELOCITY.bin', SPIN_FOOTER) and _complete(case/'BERRY_CONNECTION.bin', OPTICAL_FOOTER)
    unchanged = all(sha256(case/n) == record['input_sha256'][n] for n in ('INCAR', 'POSCAR', 'KPOINTS', 'POTCAR', 'CHGCAR'))
    unchanged &= sha256(binary) == record['binary_sha256'] and sha256(producer_manifest) == record['producer_manifest_sha256']
    success = process.returncode == 0 and normal and converged and exports and unchanged and reason is None
    outputs = ('OUTCAR', 'OSZICAR', 'EIGENVAL', 'WAVECAR', 'SPIN_VELOCITY.bin', 'BERRY_CONNECTION.bin', 'WAVEDER', 'stdout.log', 'stderr.log')
    record.update(status=reason or ('FINISHED' if success else 'FAILED_OUTPUT'), returncode=process.returncode,
                  normal_footer=normal, converged=converged, completed_export=exports, unchanged_inputs_and_binary=unchanged,
                  finished_utc=now(), wall_seconds=time.monotonic()-begin, sampled_peak_rss_bytes=peak,
                  output_sha256={n:sha256(case/n) for n in outputs if (case/n).is_file()})
    save()
    return record


def run_case(directory, binary, producer_manifest, *, timeout=3600., memory_limit_mib=4096.):
    """Retain a terminal failure record even if launch or validation raises."""
    path = Path(directory).resolve()/'run.json'
    existed = path.exists()
    state = {}
    try:
        return _run_case(directory, binary, producer_manifest, timeout=timeout,
                         memory_limit_mib=memory_limit_mib, state=state)
    except BaseException as exc:
        process = state.get('process')
        if process is not None and process.poll() is None:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL); process.wait()
        if not existed and path.is_file():
            record = json.loads(path.read_text())
            record.update(status='INTERRUPTED' if isinstance(exc, KeyboardInterrupt) else
                          ('FAILED_LAUNCH' if process is None else 'FAILED_VALIDATION'),
                          error=type(exc).__name__+': '+str(exc),
                          finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
            temporary = path.with_name('run.json.tmp')
            temporary.write_text(json.dumps(record, indent=2)+'\n'); temporary.replace(path)
        raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', required=True, type=Path)
    parser.add_argument('--binary', required=True, type=Path)
    parser.add_argument('--producer-manifest', required=True, type=Path)
    parser.add_argument('--timeout', type=float, default=3600.)
    parser.add_argument('--memory-limit-mib', type=float, default=4096.)
    args = parser.parse_args()
    result = run_case(args.run_dir, args.binary, args.producer_manifest, timeout=args.timeout, memory_limit_mib=args.memory_limit_mib)
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result['status'] == 'FINISHED' else 1)
