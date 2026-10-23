"""``pal2nal``: accepts PAL2NAL v14's command line and runs COAL instead.

    pal2nal pep.aln nuc.fasta [nuc.fasta ...] [-output clustal|paml|fasta|codon]
            [-codontable N] [-nogap] [-nomismatch] [-blockonly] [-nostderr] [-h]

An unchanged PAL2NAL command line works: the default output format is clustal, as in PAL2NAL, and all four
formats, -nogap, -nomismatch, -blockonly and every -codontable give PAL2NAL v14's output byte for byte, except for
COAL's documented deliberate differences (docs/differences.md: newer NCBI genetic-code tables, -nogap using the chosen code's
stop codons, the chosen code's start codons in the fallback, and inputs where PAL2NAL v14 writes malformed output or
fails). The command line is read as PAL2NAL v14 reads it (a dangling -output or -codontable is ignored,
"-codontable 2x" is code 2). The codon alignment goes to standard output, messages to standard error.
``-`` as the protein alignment or one DNA file reads standard input (as for ``coal``; at most one ``-``).
``-nostderr`` is accepted and mapped to COAL's ``-q``: warnings and notices are silenced, errors still go to
standard error (as with PAL2NAL's ``-nostderr``, which silences its warnings but not its ERROR messages).
``-html`` (PAL2NAL's web-server mode) gives a clear error.
COAL's own options (``--strict``, ``-j``, ``-o`` ...) are not accepted here; use ``coal`` for those.

Differences from PAL2NAL: the exit code is non-zero on errors (PAL2NAL always exits 0): 1 for a usage error or
bad input (``coal`` itself uses 2 for a usage error), 130 / 141 / 143 as for ``coal``.
"""
from __future__ import annotations

import errno
import os
import re
import sys

USAGE = __doc__

FLAG_OPTS = {"-nogap", "-nomismatch", "-blockonly"}
UNSUPPORTED = {
    "-html": ("-html (PAL2NAL's web-server mode: HTML with <pre> and coloured mismatches on standard output) is not "
              "supported by COAL; leave it out to get PAL2NAL's plain output."),
}
QUIET_OPTS = {"-nostderr"}  # mapped to coal -q
FORMATS = ("clustal", "paml", "fasta", "codon")


class UsageError(Exception):
    pass


_PERL_NUM = re.compile(r"\s*([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)")
PAL2NAL_TABLES = (1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23)


def _perl_number(v: str) -> float:
    """The number Perl reads from a string in a numeric comparison: its leading numeric part, else 0."""
    m = _PERL_NUM.match(v)
    return float(m.group(1)) if m else 0.0


def translate(argv: list[str], supported_formats=FORMATS) -> list[str]:
    """Turn PAL2NAL arguments into ``coal`` arguments (raises UsageError).

    The command line is read as PAL2NAL v14 reads it (observed by running it), not strictly: a dangling -output /
    -codontable, "-codontable 1.0", "-codontable 2x" and "-codontable -nogap" all run in PAL2NAL v14. -output takes
    the next argument whatever it is; after -codontable the flags -h, -output, -nogap, -nomismatch, -blockonly and
    -codontable are still read as options and the next other argument is the table, compared as a Perl number
    (leading numeric part: "1.0" = 1, "2x" = "+2" = "2e0" = 2; "2.5" and "0x2" are invalid); a value never given
    leaves the default."""
    files: list[str] = []
    fmt = "clustal"  # PAL2NAL's default
    table = None
    flags: list[str] = []
    want_out = want_table = False
    for a in argv:
        if a == "-h":
            print(USAGE, file=sys.stderr)  # PAL2NAL prints its help and carries on
            continue
        if a == "-output":
            want_out = True
            continue
        if want_out:
            if a not in FORMATS:
                raise UsageError(f"unknown -output format '{a}' (PAL2NAL formats: {', '.join(FORMATS)})")
            fmt, want_out = a, False
            continue
        if a in FLAG_OPTS:
            if a not in flags:
                flags.append(a)
            continue
        if a == "-codontable":
            want_table = True
            continue
        if want_table:
            num = _perl_number(a)
            if num not in PAL2NAL_TABLES:
                raise UsageError(f"invalid codontable number, {a} (PAL2NAL codes: "
                                 f"{', '.join(map(str, PAL2NAL_TABLES))})")
            table, want_table = int(num), False
            continue
        if a in UNSUPPORTED:
            raise UsageError(UNSUPPORTED[a])
        if a in QUIET_OPTS:
            if "-q" not in flags:
                flags.append("-q")
        elif a.startswith("-") and len(a) > 1:
            raise UsageError(f"unknown option {a} (PAL2NAL options: -output, -codontable, -nogap, -nomismatch, "
                             f"-blockonly, -nostderr, -h)")
        else:
            files.append(a)
    if len(files) < 2:
        raise UsageError("give a protein alignment and at least one DNA file")
    if files.count("-") > 1:
        from coal.pipeline import STDIN_TWICE
        raise UsageError(STDIN_TWICE)
    if fmt not in supported_formats:
        raise UsageError(f"this COAL cannot write '{fmt}' output. Supported: {', '.join(supported_formats)}")
    out = ["-output", fmt]
    if table is not None:
        out += ["-codontable", str(table)]
    return out + flags + ["--"] + files


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    import signal
    from coal.cli import _on_sigterm
    signal.signal(signal.SIGTERM, _on_sigterm)  # SIGTERM before the engine starts -> exit 143 too
    if not argv or argv[0] in ("-h", "-help", "--help"):
        # `pal2nal | head`: a closed stdout is a silent exit 141, as for `coal` (Windows reports it as EINVAL)
        try:
            print(USAGE)
            sys.stdout.flush()
        except OSError as e:
            if isinstance(e, BrokenPipeError) or (os.name == "nt" and e.errno == errno.EINVAL):
                from coal.cli import closed_stdout_exit
                return closed_stdout_exit()
            raise
        return 0
    # The engine runs in this process, so its own list of formats is used.
    from coal.pipeline import FORMATS as ENGINE_FORMATS
    try:
        coal_args = translate(argv, ENGINE_FORMATS)
    except UsageError as e:
        print(f"pal2nal (COAL): error: {e}", file=sys.stderr)
        return 1
    # The engine runs in this process (not as a subprocess, which a SIGTERM to this wrapper would orphan, with its
    # workers and temp files), so its own signal handling applies (Ctrl-C -> "coal: interrupted", 130;
    # SIGTERM -> cleanup, 143).
    from coal.cli import main as coal_main
    # `coal` prints at most 20 per-codon warnings by default; PAL2NAL prints them all
    return coal_main(["--max-warnings", "0", *coal_args])
