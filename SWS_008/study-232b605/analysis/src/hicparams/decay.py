"""Contact decay exponent."""
import numpy as np


def decay_exponent(table, chrom_sizes, lo=2000, hi=4_000_000, n_bins=20):
    """Return the decay exponent of cis contact frequency over log-spaced distance bins in [lo, hi)."""
    t = table.normalised()
    cis = t.c1 == t.c2
    s = (t.p2 - t.p1)[cis]
    edges = np.logspace(np.log10(lo), np.log10(hi), n_bins + 1)
    counts, _ = np.histogram(s, bins=edges)
    mid = np.sqrt(edges[:-1] * edges[1:])
    width = np.diff(edges)
    sizes = chrom_sizes.astype(float)
    positions = np.array([np.clip(sizes - m, 0, None).sum() for m in mid])
    freq = counts / (positions * width)
    ok = counts > 0
    slope = np.polyfit(np.log10(mid[ok]), np.log10(freq[ok]), 1)[0]
    return float(-slope)
