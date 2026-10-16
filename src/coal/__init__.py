"""COAL: codon alignment from a protein alignment and unaligned coding DNA."""
from .core import CoalError, codon_align

__version__ = "0.1.0"
__all__ = ["CoalError", "codon_align", "__version__"]
