"""``pal2nal``: accepts PAL2NAL v14's command line and runs COAL instead.

    pal2nal pep.aln nuc.fasta [nuc.fasta ...] [-output clustal|paml|fasta|codon]
            [-codontable N] [-nogap] [-nomismatch] [-nostderr] [-blockonly] [-html] [-h]

Pipelines that already call PAL2NAL with ``-output fasta`` can swap ``pal2nal.pl`` for
``pal2nal`` (or ``python -m coal_compat``): the codon alignment goes to standard output,
messages to standard error. Pipelines relying on PAL2NAL's default (clustal) output must add
``-output fasta`` until COAL writes clustal. ``-nostderr`` is accepted and mapped to COAL's ``-q``: per-codon warnings are silenced, but
COAL's notices about probably-wrong data and its errors still go to standard error (PAL2NAL's ``-nostderr``
silences everything). CHANGED(R2-A10): was "(COAL's ``-q``)", which implied a silent stderr.
Options COAL cannot honour (including ``-verbose``) give a clear error instead of being
silently ignored.

Differences from PAL2NAL: the exit code is non-zero on errors (PAL2NAL always exits 0),
and the default output format is PAL2NAL's (clustal), so when this COAL can only write
FASTA you must pass ``-output fasta`` (we will not silently change the format).
"""
from __future__ import annotations

import sys

from coal_app.core import detect_features

USAGE = __doc__

VALUE_OPTS = {"-output", "-codontable"}
FLAG_OPTS = {"-nogap", "-nomismatch"}
UNSUPPORTED = {
    "-blockonly": "-blockonly (Gblocks-style blocks) is not supported by COAL yet.",
    # CHANGED(APP-C2): was "use -output fasta or clustal" (clustal is not supported yet); now only what works
    "-html": "-html output is not supported by COAL; use -output fasta.",
}
QUIET_OPTS = {"-nostderr"}  # CHANGED(APP-C2): was rejected as unknown; now mapped to coal -q
FORMATS = ("clustal", "paml", "fasta", "codon")


class UsageError(Exception):
    pass


def translate(argv: list[str], supported_formats=("fasta",)) -> list[str]:
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
            flags.append(a)
        elif a in QUIET_OPTS:
            if "-q" not in flags:
                flags.append("-q")
        elif a.startswith("-") and len(a) > 1:
            raise UsageError(f"unknown option {a} (PAL2NAL options: -output, -codontable, -nogap, -nomismatch, "
                             f"-nostderr, -blockonly, -html, -h)")
        else:
            files.append(a)
        i += 1
    if len(files) < 2:
        raise UsageError("give a protein alignment and at least one DNA file")
    if fmt not in supported_formats:
        hint = " (PAL2NAL's default is clustal; add '-output fasta')" if fmt == "clustal" else ""
        raise UsageError(f"this COAL cannot write '{fmt}' output yet{hint}. "
                         f"Supported: {', '.join(supported_formats)}")
    out = files[:]
    if fmt != "fasta":
        out += ["-output", fmt]
    if table is not None:
        out += ["-codontable", str(table)]
    return out + flags


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    import signal
    from coal.cli import _on_sigterm
    signal.signal(signal.SIGTERM, _on_sigterm)  # CHANGED(R2-A4): SIGTERM before the engine starts -> exit 143 too
    if not argv or argv[0] in ("-h", "-help", "--help"):
        print(USAGE)
        return 0
    try:
        coal_args = translate(argv, detect_features().output_formats)
    except UsageError as e:
        print(f"pal2nal (COAL): error: {e}", file=sys.stderr)
        return 1
    # CHANGED(R2-A4): was subprocess.call(python -m coal ...): Ctrl-C gave a traceback and exit -2, SIGTERM
    # stopped only this wrapper (the engine kept running and writing), and the engine's workers and temp
    # files were orphaned. Now the engine runs in this process, so its own signal handling applies
    # (Ctrl-C -> "coal: interrupted", 130; SIGTERM -> cleanup, 143). See overnight/integration/CHANGES.md
    from coal.cli import main as coal_main
    return coal_main(coal_args)
