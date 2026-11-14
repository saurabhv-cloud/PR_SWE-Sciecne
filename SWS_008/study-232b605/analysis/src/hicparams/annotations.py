"""Reading annotation files."""
import os
from dataclasses import dataclass

import numpy as np


@dataclass
class Annotations:
    chrom_names: list
    chrom_sizes: np.ndarray          # bp, in chrom_names order
    ep_pairs: list                   # (chrom_index, start1, start2) per annotated pair
    domains: list                    # (chrom_index, start, end) per annotated domain


def load_annotations(directory):
    """Read chrom.sizes, ep_pairs.bedpe and domains.bed from `directory`."""
    names, sizes = [], []
    with open(os.path.join(directory, "chrom.sizes")) as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 2:
                names.append(f[0])
                sizes.append(int(f[1]))
    cidx = {c: i for i, c in enumerate(names)}
    ep = []
    with open(os.path.join(directory, "ep_pairs.bedpe")) as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 6 and f[0] == f[3]:
                ep.append((cidx[f[0]], int(f[1]), int(f[4])))
    dom = []
    with open(os.path.join(directory, "domains.bed")) as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 3:
                dom.append((cidx[f[0]], int(f[1]), int(f[2])))
    return Annotations(names, np.array(sizes, dtype=np.int64), ep, dom)
