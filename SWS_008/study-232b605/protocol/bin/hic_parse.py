#!/usr/bin/env python3
"""Convert a raw alignment-pair file into a 4DN-style .pairs file with restriction-fragment annotation.

Usage:
    hic_parse.py --settings assay_settings.json --annotations DIR --input RAW.tsv.gz --output OUT.pairs.gz

Input columns (tab separated, no header):
    read_id chrom1 pos1 strand1 mapq1 chrom2 pos2 strand2 mapq2
The fragment map named in the settings file is read from DIR.
Output is gzip text (mtime 0), sorted by chromosome order and position.
"""
import argparse
import gzip
import json
import os

import numpy as np


def read_bytes(path):
    with gzip.open(path, "rb") as fh:
        return fh.read()


def write_gz(path, data, level=3):
    with open(path, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=level) as gz:
            gz.write(data)


def read_fragments(path):
    chroms, starts, ends, names = [], [], [], []
    with open(path, "rb") as fh:
        for line in fh:
            f = line.split()
            if len(f) < 4:
                continue
            chroms.append(f[0]); starts.append(int(f[1])); ends.append(int(f[2])); names.append(f[3])
    order = []
    for c in chroms:
        if c not in order:
            order.append(c)
    chroms = np.array(chroms)
    starts = np.array(starts, dtype=np.int64)
    ends = np.array(ends, dtype=np.int64)
    names = np.array(names)
    per = {}
    for c in order:
        m = chroms == c
        per[c] = (starts[m], names[m], int(ends[m].max()))
    return order, per


def main():
    ap = argparse.ArgumentParser(description="Annotate alignment pairs with restriction fragments.")
    ap.add_argument("--settings", required=True)
    ap.add_argument("--annotations", required=True)
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    with open(args.settings) as fh:
        settings = json.load(fh)
    order, per = read_fragments(os.path.join(args.annotations, settings["fragment_map"]))
    cidx = {c: i for i, c in enumerate(order)}

    tokens = read_bytes(args.input).split()
    n = len(tokens) // 9
    arr = np.array(tokens[: n * 9]).reshape(n, 9)
    rid = arr[:, 0]
    uc1, inv1 = np.unique(arr[:, 1], return_inverse=True)
    uc2, inv2 = np.unique(arr[:, 5], return_inverse=True)
    c1 = np.array([cidx[x] for x in uc1.tolist()], dtype=np.int64)[inv1]
    c2 = np.array([cidx[x] for x in uc2.tolist()], dtype=np.int64)[inv2]
    p1 = arr[:, 2].astype(np.int64)
    p2 = arr[:, 6].astype(np.int64)
    s1, q1 = arr[:, 3], arr[:, 4]
    s2, q2 = arr[:, 7], arr[:, 8]

    swap = (c1 > c2) | ((c1 == c2) & (p1 > p2))
    c1, c2 = np.where(swap, c2, c1), np.where(swap, c1, c2)
    p1, p2 = np.where(swap, p2, p1), np.where(swap, p1, p2)
    s1, s2 = np.where(swap, s2, s1), np.where(swap, s1, s2)
    q1, q2 = np.where(swap, q2, q1), np.where(swap, q1, q2)

    sel = np.lexsort((p2, p1, c2, c1))
    rid, c1, c2, p1, p2, s1, s2, q1, q2 = (x[sel] for x in (rid, c1, c2, p1, p2, s1, s2, q1, q2))

    def frag_names(ci, pos):
        out = np.empty(len(pos), dtype="S12")
        for k, c in enumerate(order):
            m = ci == k
            if m.any():
                starts, names, _ = per[c]
                out[m] = names[np.searchsorted(starts, pos[m], side="right") - 1]
        return out

    f1 = frag_names(c1, p1)
    f2 = frag_names(c2, p2)
    cnames = np.array(order)
    header = ["## pairs format v1.0", "#sorted: chr1-chr2-pos1-pos2", "#shape: upper triangle",
              "#genome_assembly: synthetic1"]
    for c in order:
        header.append("#chromsize: %s %d" % (c.decode(), per[c][2]))
    header.append("#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type mapq1 mapq2 frag1 frag2")
    ut = np.full(len(p1), b"UU", dtype="S2")
    cols = [rid, cnames[c1], p1.astype("S12"), cnames[c2], p2.astype("S12"), s1, s2, ut, q1, q2, f1, f2]
    rows = zip(*[x.tolist() for x in cols])
    body = b"\n".join(b"\t".join(r) for r in rows)
    out = ("\n".join(header) + "\n").encode("ascii") + body + (b"\n" if len(p1) else b"")
    write_gz(args.output, out)


if __name__ == "__main__":
    main()
