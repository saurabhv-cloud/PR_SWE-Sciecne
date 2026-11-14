#!/usr/bin/env python3
"""Apply the assay settings to a .pairs file.

Usage:
    hic_filter.py --settings assay_settings.json --input IN.pairs.gz --output OUT.pairs.gz
                  [--baitmap BAITS.baitmap] [--classes DIR]

Pairs that pass are written in their original order. Header lines are copied unchanged.
--baitmap keeps only pairs with at least one end in a baited fragment (fragment names in column 4 of the baitmap).
--classes DIR also writes the rejected pairs of each class to DIR, one file per class, named after OUTPUT
(<name>.low_mapq.pairs.gz, <name>.same_fragment.pairs.gz, <name>.short_cis.pairs.gz, and <name>.off_target.pairs.gz
with --baitmap).
"""
import argparse
import gzip
import json
import os


def write_gz(path, data, level=3):
    with open(path, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=level) as gz:
            gz.write(data)


def read_baits(path):
    baits = set()
    with open(path, "rb") as fh:
        for ln in fh:
            f = ln.split()
            if len(f) >= 4 and not ln.startswith(b"#"):
                baits.add(f[3])
    return baits


def main():
    ap = argparse.ArgumentParser(description="Filter pairs by mapping quality, fragment and cis distance.")
    ap.add_argument("--settings", required=True)
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--baitmap")
    ap.add_argument("--classes")
    args = ap.parse_args()
    with open(args.settings) as fh:
        settings = json.load(fh)
    mapq_min = int(settings["mapq_min"])
    min_dist = int(settings["min_cis_distance"])
    baits = read_baits(args.baitmap) if args.baitmap else None

    with gzip.open(args.input, "rb") as fh:
        lines = fh.read().split(b"\n")
    if lines and lines[-1] == b"":
        lines.pop()
    header = [ln for ln in lines if ln.startswith(b"#")]
    body = [ln for ln in lines if not ln.startswith(b"#")]
    cols = [ln for ln in header if ln.startswith(b"#columns:")][0].split()[1:]
    ic1, ip1 = cols.index(b"chrom1"), cols.index(b"pos1")
    ic2, ip2 = cols.index(b"chrom2"), cols.index(b"pos2")
    iq1, iq2 = cols.index(b"mapq1"), cols.index(b"mapq2")
    if1, if2 = cols.index(b"frag1"), cols.index(b"frag2")
    kept = []
    rejected = {"low_mapq": [], "same_fragment": [], "short_cis": [], "off_target": []}
    for ln in body:
        f = ln.split(b"\t")
        if int(f[iq1]) < mapq_min or int(f[iq2]) < mapq_min:
            rejected["low_mapq"].append(ln)
            continue
        if f[if1] == f[if2]:
            rejected["same_fragment"].append(ln)
            continue
        if f[ic1] == f[ic2] and abs(int(f[ip2]) - int(f[ip1])) < min_dist:
            rejected["short_cis"].append(ln)
            continue
        if baits is not None and f[if1] not in baits and f[if2] not in baits:
            rejected["off_target"].append(ln)
            continue
        kept.append(ln)
    write_gz(args.output, b"\n".join(header + kept) + b"\n")

    if args.classes:
        os.makedirs(args.classes, exist_ok=True)
        name = os.path.basename(args.output)
        stem = name[: -len(".pairs.gz")] if name.endswith(".pairs.gz") else name
        for cls in ("low_mapq", "same_fragment", "short_cis", "off_target"):
            if cls == "off_target" and baits is None:
                continue
            path = os.path.join(args.classes, "%s.%s.pairs.gz" % (stem, cls))
            write_gz(path, b"\n".join(header + rejected[cls]) + b"\n")


if __name__ == "__main__":
    main()
