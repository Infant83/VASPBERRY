"""Restore complete Hamiltonian and position matrices from ordered gzip parts.

Usage: python restore_compact_operators.py INPUT_DIRECTORY NEW_OUTPUT_DIRECTORY
The destination must not exist. Both full operators are checked before success.
"""
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import os


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def restore(source, destination):
    if destination.exists():
        raise FileExistsError(f'Destination already exists: {destination}')
    manifest = json.loads((source/'export.json').read_text())
    archive = hashlib.sha256()
    for key in ('HH', 'AA'):
        parts = manifest[key]['compressed_parts']
        assert parts and len({p['file'] for p in parts}) == len(parts)
        for index, part in enumerate(parts):
            expected = f'wannier90_{key}_R.dat.part{index:03d}.gz'
            if part['file'] != expected:
                raise ValueError(f'Incorrect or unsafe part order: {part["file"]}')
            path = source/expected
            if path.stat().st_size != part['bytes'] or sha(path) != part['sha256']:
                raise ValueError(f'Input part failed verification: {expected}')
            with path.open('rb') as f:
                for block in iter(lambda: f.read(1024*1024), b''):
                    archive.update(block)
    complete = hashlib.sha256()
    # Validate decompression and all complete-file checks before creating any
    # user output. A damaged or inconsistent package leaves no output directory.
    for key in ('HH', 'AA'):
        h = hashlib.sha256()
        size = 0
        for part in manifest[key]['compressed_parts']:
            part_size = 0
            with gzip.open(source/part['file'], 'rb') as z:
                for block in iter(lambda: z.read(1024*1024), b''):
                    h.update(block)
                    complete.update(block)
                    size += len(block)
                    part_size += len(block)
            if part_size != part['uncompressed_bytes']:
                raise ValueError(f'Decompressed part size mismatch: {part["file"]}')
        if size != manifest[key]['bytes'] or h.hexdigest() != manifest[key]['sha256']:
            raise ValueError(f'Complete operator prevalidation failed: {key}')
    for name, actual in [('complete_uncompressed_payload_sha256', complete.hexdigest()),
                         ('complete_ordered_gzip_archive_sha256', archive.hexdigest())]:
        if name in manifest and manifest[name] != actual:
            raise ValueError(f'Complete archive prevalidation failed: {name}')
    destination.mkdir(parents=False, exist_ok=False)
    record = dict(status='RESTORING', operators={})
    log = destination/'restoration.json'
    log.write_text(json.dumps(record, indent=2)+'\n')
    combined = hashlib.sha256()
    for key in ('HH', 'AA'):
        name = f'wannier90_{key}_R.dat'
        staging = destination/(name+'.partial')
        h = hashlib.sha256()
        size = 0
        with staging.open('xb') as out:
            for part in manifest[key]['compressed_parts']:
                with gzip.open(source/part['file'], 'rb') as z:
                    for block in iter(lambda: z.read(1024*1024), b''):
                        out.write(block)
                        h.update(block)
                        combined.update(block)
                        size += len(block)
        if h.hexdigest() != manifest[key]['sha256'] or size != manifest[key]['bytes']:
            raise ValueError(f'Reconstructed operator failed verification: {name}')
        # Link is atomic and fails if a target appeared after initial checks.
        os.link(staging, destination/name)
        staging.unlink()
        record['operators'][name] = dict(bytes=size, sha256=h.hexdigest())
    expected = manifest.get('complete_uncompressed_payload_sha256')
    if expected is not None and combined.hexdigest() != expected:
        raise ValueError('Complete ordered payload failed verification')
    record.update(status='PASS', complete_uncompressed_payload_sha256=combined.hexdigest())
    log.write_text(json.dumps(record, indent=2)+'\n')
    return record


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('source', type=Path)
    p.add_argument('destination', type=Path)
    a = p.parse_args()
    print(json.dumps(restore(a.source, a.destination), indent=2))
