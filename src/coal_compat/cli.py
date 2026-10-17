"""``pal2nal``: accepts PAL2NAL v14's command line and runs COAL instead.

    pal2nal pep.aln nuc.fasta [nuc.fasta ...] [-output clustal|paml|fasta|codon]
            [-codontable N] [-nogap] [-nomismatch] [-blockonly] [-nostderr] [-h]

An unchanged PAL2NAL command line works: the default output format is clustal, as in PAL2NAL, and all four
formats, -nogap, -nomismatch, -blockonly and every -codontable give PAL2NAL v14's output byte for byte on
well-formed input (COAL's documented deliberate differences, where PAL2NAL silently writes wrong output, are listed
in README.md). The codon alignment goes to standard output,
messages to standard error. ``-nostderr`` is accepted and mapped to COAL's ``-q``: warnings and
notices are silenced, errors still go to standard error (as with PAL2NAL's ``-nostderr``, which silences its
warnings but not its ERROR messages). CHANGED(R2-A10): was "(COAL's ``-q``)", which implied a silent stderr.
CHANGED(DEC-2): was "COAL's notices about probably-wrong data ... still go to standard error".
``-html`` (PAL2NAL's web-server mode: HTML on standard output, messages included) gives a clear error.
COAL's own options (``--strict``, ``-j``, ``-o`` ...) are not accepted here; use ``coal`` for those.

Differences from PAL2NAL: the exit code is non-zero on errors (PAL2NAL always exits 0).
CHANGED(DEC-4): was "the default output format is PAL2NAL's (clustal), so when this COAL can only write FASTA you
must pass -output fasta"; now clustal/paml/codon are written and the default is clustal.
See overnight/integration/CHANGES.md
"""
from __future__ import annotations

import sys

USAGE = __doc__

VALUE_OPTS = {"-output", "-codontable"}
FLAG_OPTS = {"-nogap", "-nomismatch", "-blockonly"}  # CHANGED(DEC-4): -blockonly was rejected as unsupported
UNSUPPORTED = {
    # CHANGED(APP-C2): was "use -output fasta or clustal" (clustal is not supported yet); now only what works
    # CHANGED(DEC-4): wording; all four formats are supported now
    "-html": ("-html (PAL2NAL's web-server mode: HTML with <pre> and coloured mismatches on standard output) is not "
              "supported by COAL; leave it out to get PAL2NAL's plain output."),
}
QUIET_OPTS = {"-nostderr"}  # CHANGED(APP-C2): was rejected as unknown; now mapped to coal -q
FORMATS = ("clustal", "paml", "fasta", "codon")


class UsageError(Exception):
    pass


def translate(argv: list[str], supported_formats=FORMATS) -> list[str]:
    """Turn PAL2NAL arguments into ``coal`` arguments (raises UsageError)."""
    files: list[str] = []
    fmt = "clustal"  # PAL2NAL's default
    table = None
    flags: list[str] = []
    i = 0
    while i < len(argv):
        a = argv[i]
        if a in UNSUPPORTED:
            raise UsageError(UNSUPPORTED[a])
        if a in VALUE_OPTS:
            if i + 1 >= len(argv):
                raise UsageError(f"{a} needs a value")
            v = argv[i + 1]
            if a == "-output":
                if v not in FORMATS:
                    raise UsageError(f"unknown -output format '{v}' (PAL2NAL formats: {', '.join(FORMATS)})")
                fmt = v
            else:
                try:
                    table = int(v)
                except ValueError:
                    raise UsageError(f"-codontable needs a number, got '{v}'") from None
            i += 2
            continue
        if a in FLAG_OPTS:
            if a not in flags:
                flags.append(a)
        elif a in QUIET_OPTS:
            if "-q" not in flags:
                flags.append("-q")
        elif a.startswith("-") and len(a) > 1:
            raise UsageError(f"unknown option {a} (PAL2NAL options: -output, -codontable, -nogap, -nomismatch, "
                             f"-blockonly, -nostderr, -h)")
        else:
            files.append(a)
        i += 1
    if len(files) < 2:
        raise UsageError("give a protein alignment and at least one DNA file")
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
    signal.signal(signal.SIGTERM, _on_sigterm)  # CHANGED(R2-A4): SIGTERM before the engine starts -> exit 143 too
    if not argv or argv[0] in ("-h", "-help", "--help"):
        print(USAGE)
        return 0
    # CHANGED(DEC-4): was detect_features() (a `python -m coal --help` subprocess) to find the formats; the engine
    # runs in this process (R2-A4), so its own list is used. See overnight/integration/CHANGES.md
    from coal.pipeline import FORMATS as ENGINE_FORMATS
    try:
        coal_args = translate(argv, ENGINE_FORMATS)
    except UsageError as e:
        print(f"pal2nal (COAL): error: {e}", file=sys.stderr)
        return 1
    # CHANGED(R2-A4): was subprocess.call(python -m coal ...): Ctrl-C gave a traceback and exit -2, SIGTERM
    # stopped only this wrapper (the engine kept running and writing), and the engine's workers and temp
    # files were orphaned. Now the engine runs in this process, so its own signal handling applies
    # (Ctrl-C -> "coal: interrupted", 130; SIGTERM -> cleanup, 143). See overnight/integration/CHANGES.md
    from coal.cli import main as coal_main
    return coal_main(coal_args)
