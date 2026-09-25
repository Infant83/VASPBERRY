"""Expand exact official-postw90 H and position operators into its effective format.

Input is the export-only v3.1.0 hook's native little-endian stream. No physical
matrix is fitted or thresholded. Pair-dependent Wigner-Seitz translations and
both degeneracy factors are absorbed before using use_ws_distance=false.
"""
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import numpy as np


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def read_exact(path):
    with path.open('rb') as f:
        def read(dtype, shape):
            a = np.fromfile(f, dtype=dtype, count=int(np.prod(shape)))
            if a.size != np.prod(shape):
                raise ValueError('Truncated exact operator stream')
            return a.reshape(shape, order='F')
        n, nr, nd = read('<i4', (3,)).tolist()
        assert 1 <= n <= 10000 and 1 <= nr <= 100000 and 1 <= nd <= 1000
        lattice = read('<f8', (3, 3))
        irvec = read('<i4', (3, nr))
        ndegen = read('<i4', (nr,))
        hh = read('<c16', (n, n, nr))
        aa = read('<c16', (n, n, nr, 3))
        wdeg = read('<i4', (n, n, nr))
        shifts = read('<i4', (3, nd, n, n, nr))
        assert f.read(1) == b''
    assert np.all(ndegen > 0) and np.all((wdeg > 0) & (wdeg <= nd))
    assert np.isfinite(lattice).all() and np.isfinite(hh).all() and np.isfinite(aa).all()
    return lattice, irvec, ndegen, hh, aa, wdeg, shifts


def field(value):
    # An explicit decimal point overrides F12.6's implied decimal. A scientific
    # exponent is accepted by the standard Fortran F input descriptor as well.
    for precision in range(11, 4, -1):
        s = format(float(value), f'.{precision}g')
        if '.' not in s.split('e')[0]:
            s = s.replace('e', '.e') if 'e' in s else s + '.'
        s = s.replace('e-0', 'e-').replace('e+0', 'e+')
        if len(s) <= 12:
            return s.rjust(12)
    raise ValueError(f'Value cannot be represented in a 12-character field: {value}')


def write_effective(path, vectors, operators, header):
    # Readers sum repeated entries. Two short fields preserve the input double
    # to roundoff, rather than reducing it to the format's usual six decimals.
    max_error = 0.
    rows = 0
    with path.open('w', newline='\n') as f:
        f.write(header)
        for r in vectors:
            block = operators[r]
            mask = np.any(block != 0, axis=2)
            indices = np.argwhere(mask)
            if not len(indices):
                indices = np.array([[0, 0]])  # preserve identical R order
            for j, i in indices:
                values = block[j, i].view(np.float64)
                prefix = ''.join(f'{v:5d}' for v in (*r, j+1, i+1))
                first = [field(v) for v in values]
                parsed = np.array([float(v) for v in first])
                residual = values - parsed
                f.write(prefix + ''.join(first) + '\n')
                rows += 1
                if np.any(residual != 0):
                    second = [field(v) for v in residual]
                    f.write(prefix + ''.join(second) + '\n')
                    restored = parsed + np.array([float(v) for v in second])
                    max_error = max(max_error, float(np.max(np.abs(restored-values))))
                    rows += 1
    return dict(rows=rows, max_serialization_absolute_error=max_error,
                bytes=path.stat().st_size, sha256=sha(path))


def compress_parts(path, max_uncompressed=24*1024*1024):
    records = []
    # Byte chunks are concatenated after decompression, so no matrix rows or
    # digits are dropped even if a split occurs inside a text line.
    with path.open('rb') as src:
        index = 0
        while chunk := src.read(max_uncompressed):
            out = path.with_name(path.name + f'.part{index:03d}.gz')
            with out.open('wb') as f:
                with gzip.GzipFile(fileobj=f, mode='wb', filename='', mtime=0, compresslevel=6) as z:
                    z.write(chunk)
            records.append(dict(file=out.name, bytes=out.stat().st_size,
                                uncompressed_bytes=len(chunk), sha256=sha(out)))
            assert out.stat().st_size < 15*1024*1024
            index += 1
    return records


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(exist_ok=False)
    lattice, irvec, ndegen, hh, aa, wdeg, shifts = read_exact(args.input)
    n, _, nr = hh.shape
    hh_out, aa_out = {}, {}
    for ir in range(nr):
        for j in range(n):
            for i in range(n):
                count = int(wdeg[j, i, ir])
                denominator = int(ndegen[ir]) * count
                for ideg in range(count):
                    r = tuple(shifts[:, ideg, j, i, ir].tolist())
                    if r not in hh_out:
                        hh_out[r] = np.zeros((n, n, 1), complex)
                        aa_out[r] = np.zeros((n, n, 3), complex)
                    hh_out[r][j, i, 0] += hh[j, i, ir] / denominator
                    aa_out[r][j, i] += aa[j, i, ir] / denominator
    vectors = sorted(hh_out)
    hh_path = args.output/'wannier90_HH_R.dat'
    aa_path = args.output/'wannier90_AA_R.dat'
    h_record = write_effective(hh_path, vectors, hh_out,
                               f'Exact full HH_R with true WS translations expanded\n{n}\n{len(vectors)}\n')
    a_record = write_effective(aa_path, vectors, aa_out,
                               'Exact Hermitized AA_R; transl_inv=false; true WS translations expanded\n')
    h_record['compressed_parts'] = compress_parts(hh_path)
    a_record['compressed_parts'] = compress_parts(aa_path)
    complete = hashlib.sha256()
    archive = hashlib.sha256()
    for path, record in [(hh_path, h_record), (aa_path, a_record)]:
        with path.open('rb') as f:
            for block in iter(lambda: f.read(1024*1024), b''):
                complete.update(block)
        for part in record['compressed_parts']:
            with (args.output/part['file']).open('rb') as f:
                for block in iter(lambda: f.read(1024*1024), b''):
                    archive.update(block)
    np.savez_compressed(args.output/'expanded-operators.npz',
                        lattice_A=lattice, irvec=np.array(vectors),
                        HH=np.stack([hh_out[r][..., 0] for r in vectors]),
                        AA=np.stack([aa_out[r] for r in vectors]))
    result = dict(status='EXPORTED_PENDING_PHYSICAL_PARITY', source=str(args.input),
                  source_sha256=sha(args.input), num_wann=n, original_nrpts=nr,
                  expanded_nrpts=len(vectors), source_transl_inv=False,
                  source_use_ws_distance=True, effective_use_ws_distance=False,
                  source_operator='Full official postw90 HH_R and Hermitized AA_R',
                  filtering='Only exact complex zeros omitted; no magnitude cutoff',
                  effective_ndegen=1, lattice_A=lattice.tolist(),
                  HH=h_record, AA=a_record,
                  complete_payload_order=['HH', 'AA'],
                  complete_uncompressed_payload_sha256=complete.hexdigest(),
                  complete_ordered_gzip_archive_sha256=archive.hexdigest(),
                  complete_archive_definition='Concatenate compressed_parts in listed order, HH then AA; uncompressed payload consists of full HH then full AA, boundaries fixed by each operator byte count.',
                  reconstruction='Concatenate decompressed .partNNN.gz files in numeric order.')
    (args.output/'export.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
