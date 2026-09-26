"""CI adaptation of the public INIs; runs current native source and records results."""
import argparse
import configparser
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            result.update(block)
    return result.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wavecar', type=Path, required=True)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    receipt = dict(status='RUNNING', commands=[], inputs={
        'wavecar': sha(args.wavecar), 'binary': sha(args.binary)})
    try:
        for name in ('bi', 'bi-rescan'):
            config = configparser.ConfigParser(interpolation=None)
            source = ROOT / f'examples/features/simple-postprocess/{name}.ini'
            config.read(source)
            receipt['inputs'][name+'.ini'] = sha(source)
            config['run'].update(wavecar=str(args.wavecar.resolve()),
                                 binary=str(args.binary.resolve()),
                                 output=str(output/name), mpi_procs='1')
            settings = output/(name+'.ini')
            with settings.open('w') as stream:
                config.write(stream)
            reuse = ['--reuse', str(output/'bi')] if name == 'bi-rescan' else []
            for stage, arguments in (('check', [str(settings), *reuse]),
                                     ('run', [str(settings), *reuse]),
                                     ('plot', [str(output/name)])):
                command = [sys.executable, str(ROOT/'tools/vaspberry_post.py'), stage, *arguments]
                entry = dict(command=command, status='RUNNING')
                receipt['commands'].append(entry)
                with (output/(name+'-'+stage+'.stdout.log')).open('w') as out, \
                     (output/(name+'-'+stage+'.stderr.log')).open('w') as err:
                    result = subprocess.run(command, cwd=ROOT, stdout=out, stderr=err)
                entry.update(exit_code=result.returncode, status='PASS' if result.returncode == 0 else 'FAIL')
                if result.returncode:
                    raise RuntimeError(f'{name} {stage} failed; see retained logs')
            with np.load(output/name/'hall/conductivity.npz', allow_pickle=False) as data:
                mask = data['region'] == 'total'
                assert mask.sum() == (9 if name == 'bi' else 13)
                # Finite-input regression reference, not a convergence/zero-Hall claim.
                np.testing.assert_allclose(data['sigma_e2_over_h'][mask], -6.27668023e-6,
                                           atol=1e-9, rtol=0)
                np.testing.assert_allclose(data['delta_sigma_e2_over_h'][mask], 0, atol=1e-14, rtol=0)
                np.testing.assert_allclose(data['electrons_per_cell'][mask], 10, atol=1e-12, rtol=0)
                np.testing.assert_allclose(data['delta_electrons_per_cell'][mask], 0, atol=1e-14, rtol=0)
            for extension in ('png', 'pdf', 'svg'):
                assert (output/name/f'figures/charge-hall/hall.{extension}').stat().st_size > 100
        for name in ('pairs.npz', 'pairs.json'):
            assert sha(output/'bi'/'pairs'/name) == sha(output/'bi-rescan'/'pairs'/name)
        with (output/'bi/native/PAIRS.csv').open() as stream:
            assert sum(1 for line in stream if line.strip() and not line.startswith('#'))-1 == 22032
        reused = json.loads((output/'bi-rescan/run.json').read_text())
        assert [row['name'] for row in reused['stages']] == ['charge-hall']
        assert not (output/'bi-rescan/native').exists()
        receipt['status'] = 'PASS'
    except BaseException as exc:
        receipt.update(status='FAIL', error=repr(exc))
        raise
    finally:
        receipt['outputs'] = {str(p.relative_to(output)): sha(p)
                              for p in sorted(output.rglob('*')) if p.is_file()}
        (output/'validation.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print('Public settings example: native pairs, Hall scans, reuse and figures PASS')


if __name__ == '__main__':
    main()
