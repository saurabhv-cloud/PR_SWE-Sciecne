"""Reading .pairs files into arrays."""
import gzip
from dataclasses import dataclass

import numpy as np


@dataclass
class PairTable:
    """Contact pairs as integer arrays; chromosome indices follow the chromosome-size file order."""

    c1: np.ndarray
    p1: np.ndarray
    c2: np.ndarray
    p2: np.ndarray

    def __len__(self):
        return len(self.p1)

    def take(self, idx):
        return PairTable(self.c1[idx], self.p1[idx], self.c2[idx], self.p2[idx])

    def normalised(self):
        """Return a table with the lower chromosome (then lower position) as end 1."""
        swap = (self.c1 > self.c2) | ((self.c1 == self.c2) & (self.p1 > self.p2))
        return PairTable(
            np.where(swap, self.c2, self.c1),
            np.where(swap, self.p2, self.p1),
            np.where(swap, self.c1, self.c2),
            np.where(swap, self.p1, self.p2),
        )


def read_pairs(path, chrom_order):
    """Read a (gzip) .pairs file. `chrom_order` is the list of chromosome names."""
    cidx = {c.encode(): i for i, c in enumerate(chrom_order)}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rb") as fh:
        data = fh.read()
    lines = data.split(b"\n")
    cols = None
    body = []
    for ln in lines:
        if not ln:
            continue
        if ln.startswith(b"#"):
            if ln.startswith(b"#columns:"):
                cols = ln.split()[1:]
            continue
        body.append(ln)
    if cols is None:
        raise ValueError("no #columns header in %s" % path)
    ic1, ip1 = cols.index(b"chrom1"), cols.index(b"pos1")
    ic2, ip2 = cols.index(b"chrom2"), cols.index(b"pos2")
    n = len(body)
    if n == 0:
        z = np.zeros(0, dtype=np.int64)
        return PairTable(z, z, z, z)
    rows = [ln.split(b"\t") for ln in body]
    arr = np.array(rows)
    names1, inv1 = np.unique(arr[:, ic1], return_inverse=True)
    names2, inv2 = np.unique(arr[:, ic2], return_inverse=True)
    c1 = np.array([cidx[x] for x in names1.tolist()], dtype=np.int64)[inv1]
    c2 = np.array([cidx[x] for x in names2.tolist()], dtype=np.int64)[inv2]
    p1 = arr[:, ip1].astype(np.int64)
    p2 = arr[:, ip2].astype(np.int64)
    return PairTable(c1, p1, c2, p2)
