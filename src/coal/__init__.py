"""COAL: codon alignment from a protein alignment and unaligned coding DNA.

Public API (docs/python-api.md): codon_align (in memory), align_files (files in, file out, as the `coal` command)
and CoalError. Everything else in the package is internal.
"""
from .core import CoalError, codon_align

__version__ = "0.4.2"

from .pipeline import run as align_files  # noqa: E402  (pipeline needs __version__)

__all__ = ["CoalError", "codon_align", "align_files", "__version__"]
