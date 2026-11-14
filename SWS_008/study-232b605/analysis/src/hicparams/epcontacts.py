"""Contact frequency at annotated enhancer-promoter pairs."""
from .matrices import cis_matrices


def ep_contact_frequency(table, chrom_sizes, ep_pairs, bin_size=10_000, radius=1):
    """Return the contacts in the windows around annotated enhancer-promoter pairs per million pairs in `table`."""
    mats = cis_matrices(table, chrom_sizes, bin_size)
    obs = 0.0
    for ci, start1, start2 in ep_pairs:
        m = mats[ci]
        n = m.shape[0]
        i, j = start1 // bin_size, start2 // bin_size
        if i - radius < 0 or j + radius >= n:
            continue
        obs += m[i - radius:i + radius + 1, j - radius:j + radius + 1].sum()
    return float(obs / (len(table) / 1e6))
