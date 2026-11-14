"""Domain score from annotated domains."""
import numpy as np


def domain_score(table, domains):
    """Return the mean domain score over `domains`, a list of (chrom_index, start, end) tuples."""
    domains = sorted(domains)
    keys = np.array([c * (1 << 32) + a for c, a, _ in domains], dtype=np.int64)
    t = table.normalised()
    cis = t.c1 == t.c2
    d1 = np.searchsorted(keys, t.c1[cis] * (1 << 32) + t.p1[cis], side="right") - 1
    d2 = np.searchsorted(keys, t.c2[cis] * (1 << 32) + t.p2[cis], side="right") - 1
    within = d1 == d2
    n = len(keys)
    n_within = np.bincount(d1[within], minlength=n)
    n_touching = np.bincount(d1, minlength=n) + np.bincount(d2, minlength=n) - n_within
    ok = n_touching > 0
    return float(np.mean(n_within[ok] / n_touching[ok]))
