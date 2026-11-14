#!/usr/bin/env python3
"""Collect the per-library QC files named by a QC index into one table.

Usage:
    qc_table.py QC_INDEX.csv OUT.tsv

QC_INDEX.csv has a header row with `library`, `condition`, `passed` and `qc` columns. Relative `qc` paths are taken
relative to the index file. OUT.tsv has one row per library, sorted by library id.
"""
import csv
import json
import os
import sys

FIELDS = (
    "n_parsed", "n_dedup", "n_filtered", "n_duplicates", "duplicate_rate", "n_low_mapq", "low_mapq_rate",
    "n_same_fragment", "same_fragment_rate", "n_short_cis", "short_cis_rate", "short_range_rate",
    "cis_fraction_dedup", "retained_fraction",
)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 2:
        sys.stderr.write(__doc__)
        return 2
    index, out = argv
    index_dir = os.path.dirname(os.path.abspath(index))
    with open(index, newline="") as fh:
        rows = list(csv.DictReader(fh))
    lines = ["\t".join(("library", "condition", "passed") + FIELDS)]
    for row in sorted(rows, key=lambda r: r["library"]):
        path = row["qc"]
        if not os.path.isabs(path):
            path = os.path.join(index_dir, path)
        with open(path) as fh:
            qc = json.load(fh)
        cells = [row["library"], row["condition"], "true" if qc["passed"] else "false"]
        cells += [str(qc[k]) for k in FIELDS]
        lines.append("\t".join(cells))
    with open(out, "wb") as fh:
        fh.write(("\n".join(lines) + "\n").encode("ascii"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
