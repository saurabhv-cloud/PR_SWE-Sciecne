"""Command line: python -m hicparams INDEX.csv COLUMN OUT.tsv [--annotations DIR]

INDEX.csv has a header row with a `library` column and one path column per pair file.
Relative paths are taken relative to the index file. For every row the pair file named in COLUMN is read
and one row of parameters is written to OUT.tsv. The annotation directory is taken from --annotations,
then the HICPARAMS_ANNOTATIONS environment variable, then an `annotations` directory next to the index.
"""
import argparse
import csv
import os

from . import PARAMETERS, compute_parameters, load_annotations, read_pairs


def main(argv=None):
    ap = argparse.ArgumentParser(prog="hicparams", description="Estimate Hi-C parameters per library.")
    ap.add_argument("index")
    ap.add_argument("column")
    ap.add_argument("out")
    ap.add_argument("--annotations")
    args = ap.parse_args(argv)

    index_dir = os.path.dirname(os.path.abspath(args.index))
    ann_dir = args.annotations or os.environ.get("HICPARAMS_ANNOTATIONS") or os.path.join(index_dir, "annotations")
    ann = load_annotations(ann_dir)
    with open(args.index, newline="") as fh:
        rows = list(csv.DictReader(fh))
    lines = ["\t".join(("library",) + PARAMETERS)]
    for row in rows:
        path = row[args.column]
        if not os.path.isabs(path):
            path = os.path.join(index_dir, path)
        table = read_pairs(path, ann.chrom_names)
        res = compute_parameters(table, ann)
        lines.append("\t".join([row["library"]] + ["%.6f" % res[p] for p in PARAMETERS]))
    with open(args.out, "wb") as fh:
        fh.write(("\n".join(lines) + "\n").encode("ascii"))


if __name__ == "__main__":
    main()
