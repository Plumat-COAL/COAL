"""COAL: codon alignments from protein alignments, Ka/Ks (dN/dS) and the one-step pipeline.

Public API (docs/python-api.md):

  * `import coal`: codon_align (in memory), align_files (files in, file out, as the `coal` command), CoalError and
    __version__. These four names are what the package exports (__all__).
  * `coal.dnds_pipeline`: run_pairs, run_clusters and their PipelineResult, codon_dnds, the column constants
    (PAIR_COLUMNS, SUMMARY_COLUMNS, FAILED_COLUMNS, SCHEMA), find_diamond / install_diamond.
  * `coal.citations`: what to cite for what was run. `coal.align`, `coal.msa`: protein alignments.

The supported names of each module are its __all__; python-api.md has their signatures and the versioning promise.
Everything else in the package (coal.core, coal.pipeline, coal.report, coal.seqio, names with a leading
underscore) is internal. The package is typed (py.typed).
"""
from .core import CoalError, codon_align

__version__ = "0.6.0.dev0"

from .pipeline import run as align_files  # noqa: E402  (pipeline needs __version__)

__all__ = ["CoalError", "codon_align", "align_files", "__version__"]
