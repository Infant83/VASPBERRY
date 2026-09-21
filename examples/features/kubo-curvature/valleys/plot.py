#!/usr/bin/env python3
"""Redraw the published MoS2 valley figures from the exported native samples."""
from pathlib import Path
import argparse
import csv
import json
from run import plot, sha


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",type=Path,default=Path(__file__).parent/"reference/summary.csv")
    parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args()
    with args.input.open() as handle:rows=list(csv.DictReader(handle))
    if len(rows)!=162:parser.error("expected 162 exported valley samples")
    for row in rows:
        for key in row:
            if key!="valley":row[key]=float(row[key])
    if args.output_dir.exists():parser.error("output exists; choose a new directory")
    args.output_dir.mkdir(parents=True)
    record=plot(rows,args.output_dir)
    record["source_sha256"]=sha(args.input)
    (args.output_dir/"plot.json").write_text(json.dumps(record,indent=2)+"\n")
    print(args.output_dir/"figure.png")


if __name__=="__main__":main()
