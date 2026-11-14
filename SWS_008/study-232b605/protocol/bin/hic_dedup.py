#!/usr/bin/env python3
"""Remove duplicate pairs from a .pairs file.

Usage:
    hic_dedup.py --settings assay_settings.json --input IN.pairs.gz --output OUT.pairs.gz

The first occurrence of each pair is kept. Header lines are copied unchanged.
"""
import argparse
import gzip
import json

RULES = {"exact_coordinates_strand": ("chrom1", "pos1", "strand1", "chrom2", "pos2", "strand2")}


def write_gz(path, data, level=3):
    with open(path, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=level) as gz:
            gz.write(data)


def main():
    ap = argparse.ArgumentParser(description="Remove duplicate pairs.")
    ap.add_argument("--settings", required=True)
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    with open(args.settings) as fh:
        settings = json.load(fh)
    key_cols = RULES[settings["dedup_rule"]]

    with gzip.open(args.input, "rb") as fh:
        lines = fh.read().split(b"\n")
    if lines and lines[-1] == b"":
        lines.pop()
    header = [ln for ln in lines if ln.startswith(b"#")]
    body = [ln for ln in lines if not ln.startswith(b"#")]
    cols = [ln for ln in header if ln.startswith(b"#columns:")][0].split()[1:]
    i0, i1, i2, i3, i4, i5 = [cols.index(c.encode()) for c in key_cols]
    seen = set()
    kept = []
    for ln in body:
        f = ln.split(b"\t")
        key = (f[i0], f[i1], f[i2], f[i3], f[i4], f[i5])
        if key in seen:
            continue
        seen.add(key)
        kept.append(ln)
    out = b"\n".join(header + kept) + b"\n"
    write_gz(args.output, out)


if __name__ == "__main__":
    main()
