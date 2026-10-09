#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Copy the maintained VASPBERRY skill to a local agent skill directory."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys


def file_hashes(directory: Path) -> dict[str, str]:
    return {
        str(path.relative_to(directory)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(directory.rglob('*'))
        if path.is_file() and '__pycache__' not in path.parts and path.suffix != '.pyc'
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--destination', type=Path,
        default=Path.home() / '.agents' / 'skills' / 'vaspberry',
        help='Skill folder to create (default: ~/.agents/skills/vaspberry).',
    )
    parser.add_argument('--dry-run', action='store_true', help='Print paths and hashes without writing files.')
    args = parser.parse_args()
    source = Path(__file__).resolve().parents[1] / 'plugins' / 'vaspberry' / 'skills' / 'vaspberry'
    destination = args.destination.expanduser().absolute()
    try:
        required = ('SKILL.md', 'LICENSE', 'scripts/vaspberry_agent.py', 'references/setup.md')
        for name in required:
            if not (source / name).is_file():
                raise ValueError(f'Missing skill file: {source / name}. Fetch plugins/vaspberry first.')
        if destination.exists() or destination.is_symlink():
            raise FileExistsError(f'Destination already exists: {destination}. Preserve it and choose a new destination.')
        hashes = file_hashes(source)
        receipt = {'source': str(source), 'destination': str(destination),
                   'dry_run': args.dry_run, 'installed': False, 'sha256': hashes}
        if not args.dry_run:
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copytree(source, destination, ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
            if file_hashes(destination) != hashes:
                raise RuntimeError(f'Copied file hashes differ. Inspect the partial installation: {destination}')
            receipt['installed'] = True
        print(json.dumps(receipt, indent=2))
        return 0
    except (OSError, ValueError, RuntimeError) as exc:
        print(json.dumps({'installed': False, 'destination': str(destination), 'error': str(exc)}), file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
