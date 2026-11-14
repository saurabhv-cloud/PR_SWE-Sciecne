"""Cis-to-trans contact ratio."""


def cis_trans_ratio(table):
    """Return the number of cis pairs divided by the number of trans pairs."""
    cis = int((table.c1 == table.c2).sum())
    trans = len(table) - cis
    if trans == 0:
        return float("inf")
    return cis / trans
