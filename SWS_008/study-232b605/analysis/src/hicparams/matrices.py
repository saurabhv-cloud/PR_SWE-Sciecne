"""Binned cis contact matrices."""
import numpy as np


def cis_matrices(table, chrom_sizes, bin_size):
    """Return a list of dense upper-triangle count matrices, one per chromosome."""
    t = table.normalised()
    cis = t.c1 == t.c2
    out = []
    for ci, length in enumerate(chrom_sizes.tolist()):
        n = -(-length // bin_size)
        m = cis & (t.c1 == ci)
        b1 = t.p1[m] // bin_size
        b2 = t.p2[m] // bin_size
        flat = np.bincount(b1 * n + b2, minlength=n * n)
        out.append(flat.reshape(n, n))
    return out
