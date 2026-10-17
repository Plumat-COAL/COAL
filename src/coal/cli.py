"""Command-line interface.

    coal PEP_ALN NUC_FASTA [NUC_FASTA ...] [-codontable N] [-nogap] [-nomismatch] [-o OUT]

PAL2NAL-style single-dash options are accepted, so COAL can be dropped into
existing pipelines; -output clustal|paml|fasta|codon and -blockonly as in PAL2NAL v14.
"""
from __future__ import annotations

import argparse
import os
import signal
import sys

from . import __version__
from .core import CoalError
from .pipeline import FORMATS, run
from .seqio import InputFormatError
from .tables import SUPPORTED


def _on_sigterm(signum, frame):
    raise SystemExit(128 + signum)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        prog="coal",
        description="Convert a protein alignment + unaligned coding DNA into a codon alignment.",
    )
    ap.add_argument("pep_aln", help="protein alignment (FASTA or CLUSTAL)")
    ap.add_argument("nuc_fasta", nargs="+", help="coding DNA sequences (FASTA; one or more files)")
    ap.add_argument("-codontable", "--codontable", type=int, default=1, choices=SUPPORTED,
                    metavar="N", help=f"NCBI genetic code (default 1); one of {', '.join(map(str, SUPPORTED))}")
    ap.add_argument("-nogap", "--nogap", action="store_true",
                    help="remove codon columns containing gaps or stop codons (stop codons of the -codontable "
                         "code, any case; PAL2NAL uses the standard code's TAA/TAG/TGA, upper case only)")
    # CHANGED(FZ-5): help said "stop codons" (the standard code's, upper case, as PAL2NAL); now the selected code's
    ap.add_argument("-nomismatch", "--nomismatch", action="store_true",
                    help="remove alignment columns where a codon does not match its residue")
    # CHANGED(DEC-4): was choices=["fasta"]; now PAL2NAL v14's four formats, byte-identical to PAL2NAL. The default
    # stays fasta (deliberate difference; the drop-in `pal2nal` command defaults to clustal like PAL2NAL)
    ap.add_argument("-output", "--output-format", default="fasta", choices=list(FORMATS),
                    help="output format, as PAL2NAL v14 writes it (default fasta; PAL2NAL's default is clustal)")
    ap.add_argument("-blockonly", "--blockonly", action="store_true",
                    help="keep only the columns marked '#' in the CLUSTAL input's conservation lines (PAL2NAL's "
                         "-blockonly; no effect if the input has no '#' marks)")
    ap.add_argument("-o", "--out", help="output file (default: stdout)")
    ap.add_argument("-j", "--jobs", type=int, default=0, metavar="N",
                    help="maximum worker processes (default 0 = automatic: all CPU cores for large "
                         "inputs, one process for small ones; 1 = never use extra processes)")
    # CHANGED(R2-E4): help was "do not print warnings" (notices are printed anyway)
    # CHANGED(DEC-2): was "notices about probably-wrong data are still printed"; now -q silences them too
    ap.add_argument("-q", "--quiet", action="store_true",
                    help="do not print warnings or notices (errors are always printed)")
    # CHANGED(DEC-1): new. Was: these three rules were always errors (ENG-c, R3-E1) or a warning (ENG-d). Now
    # PAL2NAL's behaviour by default (output + notice) and errors only with --strict. See overnight/integration/CHANGES.md
    ap.add_argument("-strict", "--strict", action="store_true",
                    help="refuse probably-wrong data instead of writing it as PAL2NAL does (PAL2NAL never refuses): "
                         "exit 1 with no output if (1) more than half of a sequence's checkable codons (A/C/G/T "
                         "only) do not match its protein; (2) no 10-residue segment of a protein is found in its "
                         "DNA, the DNA is longer than the gene and fewer than 10 codons are checkable (where the gene "
                         "starts can't be told); (3) only some protein IDs are DNA IDs (pairing by order would "
                         "follow). Without --strict these are notices on stderr")
    ap.add_argument("--version", action="version", version=f"coal {__version__}")
    args = ap.parse_args(argv)

    warn = None if args.quiet else (lambda w: print(f"coal: warning: {w}", file=sys.stderr))
    # CHANGED(ENG-c): new; notices are printed even with -q. See overnight/engine-fixes/CHANGES.md
    # CHANGED(DEC-2): was printed even with -q; now -q silences notices too (GNU / clig.dev "quiet", and PAL2NAL's
    # -nostderr). Errors are still printed. See overnight/integration/CHANGES.md
    notice = warn
    # CHANGED(ENG-f): was default signal handling (SIGTERM left .coal_*.tmp files; Ctrl-C printed a
    # traceback; a closed stdout pipe gave "[Errno 32] Broken pipe", exit 1); now SIGTERM -> SystemExit
    # (cleanup runs, exit 143), Ctrl-C -> "coal: interrupted", exit 130, broken pipe -> silent exit 141.
    # See overnight/engine-fixes/CHANGES.md
    signal.signal(signal.SIGTERM, _on_sigterm)
    if hasattr(signal, "SIGHUP"):  # CHANGED(R2-E8): terminal closed -> same cleanup, exit 129
        signal.signal(signal.SIGHUP, _on_sigterm)
    try:
        run(args.pep_aln, args.nuc_fasta, out_path=args.out, table=args.codontable,
            nogap=args.nogap, nomismatch=args.nomismatch, jobs=args.jobs, warn=warn, notice=notice,
            strict=args.strict, out_format=args.output_format, blockonly=args.blockonly)
    except KeyboardInterrupt:
        print("coal: interrupted", file=sys.stderr)
        return 130
    except BrokenPipeError:
        try:  # stop Python from reporting the pipe again when it flushes stdout at exit
            os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno())
        except (OSError, ValueError):
            pass
        return 141
    except (CoalError, InputFormatError, UnicodeDecodeError, OSError) as e:
        print(f"coal: error: {e}", file=sys.stderr)
        return 1
    return 0
