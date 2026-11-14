"""Per-library parameter estimates from Hi-C pair files."""
from .annotations import Annotations, load_annotations
from .cistrans import cis_trans_ratio
from .decay import decay_exponent
from .domains import domain_score
from .epcontacts import ep_contact_frequency
from .reader import PairTable, read_pairs

__version__ = "0.1.0"

PARAMETERS = ("decay_exponent", "cis_trans_ratio", "domain_score", "ep_contact_frequency")


def compute_parameters(table, annotations):
    """Return a dict with the four parameters for a PairTable."""
    return {
        "decay_exponent": decay_exponent(table, annotations.chrom_sizes),
        "cis_trans_ratio": cis_trans_ratio(table),
        "domain_score": domain_score(table, annotations.domains),
        "ep_contact_frequency": ep_contact_frequency(table, annotations.chrom_sizes, annotations.ep_pairs),
    }
