#!/usr/bin/env python3
"""Fetch the checksum-pinned public Bi VASP WAVECAR used by the tutorials."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import urllib.request

BI_SHA256 = "a8d81854f2efc561938e478dde1be29a17ccc95d5d37325c122b9d9e82fa0838"
BI_BYTES = 200421600
BI_URL = ("https://media.githubusercontent.com/media/Infant83/VASPBERRY/"
          "a692e21482d24767859d02b7885dd70b234c6a2c/examples/Bi_Z2/WAVECAR")


def verify(path: Path) -> str:
    if path.stat().st_size != BI_BYTES:
        raise ValueError(f"Expected the {BI_BYTES}-byte Bi WAVECAR payload, not a pointer or different file")
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    if digest.hexdigest() != BI_SHA256:
        raise ValueError("Bi WAVECAR checksum mismatch; the tutorial requires the named public input")
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset", choices=["bi"])
    parser.add_argument("--output-dir", type=Path, required=True, help="New directory to contain the actual WAVECAR")
    parser.add_argument("--source", type=Path, help="Copy an already downloaded, identical public Bi WAVECAR")
    args = parser.parse_args()
    if args.output_dir.exists():
        parser.error("output directory exists; choose a new path")
    if args.source is not None:
        try:
            verify(args.source)
        except (OSError, ValueError) as exc:
            parser.error(str(exc))
    args.output_dir.mkdir(parents=True, exist_ok=False)
    part = args.output_dir / "WAVECAR.part"
    record = {"dataset": "Bi SOC 12x12 WAVECAR", "url": BI_URL,
              "expected_sha256": BI_SHA256, "expected_bytes": BI_BYTES,
              "started_utc": datetime.now(timezone.utc).isoformat(), "status": "RUNNING"}
    report = args.output_dir / "download.json"
    report.write_text(json.dumps(record, indent=2) + "\n")
    try:
        if args.source is not None:
            shutil.copyfile(args.source, part)
            record["copied_from"] = str(args.source.resolve())
        else:
            request = urllib.request.Request(BI_URL, headers={"User-Agent": "VASPBERRY-public-input-fetcher"})
            with urllib.request.urlopen(request, timeout=120) as source, part.open("xb") as output:
                shutil.copyfileobj(source, output, length=1024 * 1024)
        record["sha256"] = verify(part)
        part.rename(args.output_dir / "WAVECAR")
        record.update(status="PASS", exit_code=0)
    except Exception as exc:
        record.update(status="FAIL", exit_code=1, error=str(exc))
    record["finished_utc"] = datetime.now(timezone.utc).isoformat()
    report.write_text(json.dumps(record, indent=2) + "\n")
    if record["status"] != "PASS":
        parser.exit(1, f"Input fetch failed; see {report}\n")
    print(f"Verified VASP WAVECAR: {args.output_dir / 'WAVECAR'}")
    print("Use this path as --wavecar for an individual tutorial or --bi-wavecar for the batch runner.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
