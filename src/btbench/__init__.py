"""Research reference implementation; see docs/METHODS.md for model limits."""

from .core import Acquisition, FilterResult
from .filters import svd_filter, aer_filter, local_filter, rpca_filter
from .axes import motion_packet_filter, multilag_filter

__version__ = "0.5.0"
__all__ = ["Acquisition", "FilterResult", "svd_filter", "aer_filter",
           "local_filter", "rpca_filter", "motion_packet_filter", "multilag_filter"]
