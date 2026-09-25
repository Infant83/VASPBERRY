#!/usr/bin/env python3
"""Restore the actual Bi PAW operator matrices supplied in this example."""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import shutil

HERE=Path(__file__).resolve().parent

def restore(mesh, output):
    source=HERE/'inputs'/f'matrices-{mesh}x{mesh}-b48'
    manifest=json.loads((source/'manifest.json').read_text())
    if output.exists():raise ValueError('Output exists; choose a new directory')
    payload=[]
    for part in manifest['parts']:
        path=source/part['name']
        if path.parent!=source or not path.is_file():raise ValueError('Invalid matrix part')
        raw=path.read_bytes()
        if len(raw)!=part['bytes'] or hashlib.sha256(raw).hexdigest()!=part['sha256']:
            raise ValueError('Matrix part failed integrity check: '+path.name)
        payload.append(raw)
    raw=b''.join(payload)
    if len(raw)!=manifest['bytes'] or hashlib.sha256(raw).hexdigest()!=manifest['sha256']:
        raise ValueError('Assembled matrix archive failed integrity check')
    metadata=source/'physical-matrices.json'
    if hashlib.sha256(metadata.read_bytes()).hexdigest()!=manifest['metadata_sha256']:
        raise ValueError('Matrix metadata failed integrity check')
    output.mkdir(parents=True)
    (output/'physical-matrices.npz').write_bytes(raw)
    shutil.copy2(metadata,output/metadata.name)
    return output

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',type=int,choices=[6,12],default=12)
    p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args()
    try:restore(a.mesh,a.output_dir)
    except (OSError,ValueError,KeyError) as e:p.error(str(e))
    print(a.output_dir)

if __name__=='__main__':main()
