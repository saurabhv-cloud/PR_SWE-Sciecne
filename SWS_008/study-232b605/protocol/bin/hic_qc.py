#!/usr/bin/env python3
"""Write per-library counts for the three pair files of a library.

Usage:
    hic_qc.py --settings assay_settings.json --library ID --parsed P.pairs.gz --dedup D.pairs.gz
              --filtered F.pairs.gz --output ID.qc.json [--verdict ID.verdict]

The counts include `passed`: true when the deduplicated pairs number at least `min_unique_pairs` of the settings.
--verdict writes the word pass or fail to a file.
"""
import argparse
import gzip
import json


def count_pairs(path):
    n = 0
    with gzip.open(path, "rb") as fh:
        for line in fh:
            if not line.startswith(b"#"):
                n += 1
    return n


def rate(k, n):
    return round(k / n, 6) if n else 0.0


def main():
    ap = argparse.ArgumentParser(description="Per-library pair counts.")
    ap.add_argument("--settings", required=True)
    ap.add_argument("--library", required=True)
    ap.add_argument("--parsed", required=True)
    ap.add_argument("--dedup", required=True)
    ap.add_argument("--filtered", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--verdict")
    args = ap.parse_args()
    with open(args.settings) as fh:
        settings = json.load(fh)
    mapq_min = int(settings["mapq_min"])
    min_dist = int(settings["min_cis_distance"])
    min_unique = int(settings["min_unique_pairs"])

    n_parsed = count_pairs(args.parsed)
    n_filtered = count_pairs(args.filtered)
    with gzip.open(args.dedup, "rb") as fh:
        lines = fh.read().split(b"\n")
    header = [ln for ln in lines if ln.startswith(b"#")]
    body = [ln for ln in lines if ln and not ln.startswith(b"#")]
    cols = [ln for ln in header if ln.startswith(b"#columns:")][0].split()[1:]
    ic1, ip1 = cols.index(b"chrom1"), cols.index(b"pos1")
    ic2, ip2 = cols.index(b"chrom2"), cols.index(b"pos2")
    iq1, iq2 = cols.index(b"mapq1"), cols.index(b"mapq2")
    if1, if2 = cols.index(b"frag1"), cols.index(b"frag2")
    n_dedup = len(body)
    low = same = short = cis = 0
    for ln in body:
        f = ln.split(b"\t")
        if int(f[iq1]) < mapq_min or int(f[iq2]) < mapq_min:
            low += 1
        if f[if1] == f[if2]:
            same += 1
        elif f[ic1] == f[ic2] and abs(int(f[ip2]) - int(f[ip1])) < min_dist:
            short += 1
        if f[ic1] == f[ic2]:
            cis += 1

    qc = {
        "library": args.library,
        "n_parsed": n_parsed,
        "n_dedup": n_dedup,
        "n_filtered": n_filtered,
        "n_duplicates": n_parsed - n_dedup,
        "duplicate_rate": rate(n_parsed - n_dedup, n_parsed),
        "n_low_mapq": low,
        "low_mapq_rate": rate(low, n_dedup),
        "n_same_fragment": same,
        "same_fragment_rate": rate(same, n_dedup),
        "n_short_cis": short,
        "short_cis_rate": rate(short, n_dedup),
        "short_range_rate": rate(same + short, n_dedup),
        "cis_fraction_dedup": rate(cis, n_dedup),
        "retained_fraction": rate(n_filtered, n_parsed),
        "passed": n_dedup >= min_unique,
    }
    with open(args.output, "wb") as fh:
        fh.write((json.dumps(qc, indent=2, sort_keys=True) + "\n").encode("ascii"))
    if args.verdict:
        with open(args.verdict, "wb") as fh:
            fh.write(b"pass\n" if qc["passed"] else b"fail\n")


if __name__ == "__main__":
    main()
