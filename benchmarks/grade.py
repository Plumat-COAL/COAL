"""Grade one tool output for one test case into exactly one class.

Classes: TRUTH | REF | KNOWN_DIV | MISMATCH | ERROR | TIMEOUT | OOM.

Normalisation (the only one applied): line endings (CRLF/CR -> LF) and line wrapping (sequence lines of
a record are joined; trailing whitespace on each line dropped). Headers, record order, case and every
sequence character are compared exactly.

KNOWN_DIV rules (each tied to a documented divergence; see KNOWN_DIV_RULES) are applied only to an output
that is not TRUTH/REF. A rule must explain *every* differing record, otherwise the output is a MISMATCH.
"""
from __future__ import annotations

import difflib
import re
from pathlib import Path

WARN_RE = re.compile(rb"(\S+) pepAlnPos (\d+): (\S+) does not correspond to (\S+)")


def read_fasta(path) -> list[tuple[str, str]] | None:
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return None
    data = p.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n").decode("latin-1")
    recs, head, seq = [], None, []
    for line in data.split("\n"):
        line = line.rstrip()
        if line.startswith(">"):
            if head is not None:
                recs.append((head, "".join(seq)))
            head, seq = line[1:], []
        elif line:
            if head is None:
                return None  # not FASTA
            seq.append(line.strip())
    if head is not None:
        recs.append((head, "".join(seq)))
    return recs


def records_hash(recs) -> str:
    """Hash of the normalised records (what grading compares), to re-check later runs cheaply."""
    import hashlib
    h = hashlib.sha256()
    for head, seq in recs:
        h.update(head.encode("latin-1") + b"\n" + seq.encode("latin-1") + b"\n")
    return h.hexdigest()


def warnings_from_stderr(path) -> list[tuple]:
    p = Path(path)
    if not p.exists():
        return []
    return [tuple(x.decode("latin-1") for x in m.groups()) for m in WARN_RE.finditer(p.read_bytes())]


def expected_warning_set(case_dir: Path):
    f = case_dir / "expected_warnings.tsv"
    if not f.exists():
        return None
    rows = [l.split("\t") for l in f.read_text().splitlines()[1:] if l.strip()]
    return {(r[0], r[1], r[2], r[3]) for r in rows}


def write_diff(a, b, a_name, b_name, path: Path, max_lines=400) -> str:
    def lines(recs):
        return [f">{h}\n{s}\n" for h, s in recs] if recs else ["<no output>\n"]
    d = list(difflib.unified_diff(lines(a), lines(b), a_name, b_name, n=0))
    path.parent.mkdir(parents=True, exist_ok=True)
    text = "".join(d[:max_lines]) + (f"... ({len(d) - max_lines} more diff lines)\n" if len(d) > max_lines else "")
    path.write_text(text)
    return str(path)


# ---- KNOWN_DIV rules -------------------------------------------------------------------------------------
# A rule gets (tool, case meta, output records, truth records or None, p2n14 records or None) and returns a
# reason string if it explains the difference, else None. Rules are added only after triage, each citing
# its source (p2n16 tests/divergences.tsv case name, or COAL README "Deliberate differences" bullet).
KNOWN_DIV_RULES: list = []


def classify(tool, meta, case_dir: Path, run_status, out_path, ref_recs, diff_path: Path):
    """Return (class, reason, diff_path or '')."""
    if run_status.status in ("TIMEOUT", "OOM"):
        return run_status.status, run_status.reason, ""
    out = read_fasta(out_path) if run_status.status == "OK" else None
    truth = read_fasta(case_dir / meta["truth"]) if meta.get("truth") else None
    if run_status.status != "OK" or out is None:
        reason = run_status.reason or f"exit={run_status.exit_code}: empty/unparseable output"
        return "ERROR", reason, ""
    if truth is not None:
        if out == truth:
            return "TRUTH", "", ""
        target, tname = truth, "truth"
    else:
        if tool != "p2n14" and ref_recs is not None and out == ref_recs:
            return "REF", "", ""
        if tool == "p2n14":
            return "REF", "is the reference", ""
        if ref_recs is None:
            return "MISMATCH", "no truth and p2n14 produced no output", write_diff(None, out, "p2n14", tool, diff_path)
        target, tname = ref_recs, "p2n14"
    for rule in KNOWN_DIV_RULES:
        why = rule(tool, meta, out, truth, ref_recs)
        if why:
            return "KNOWN_DIV", why, write_diff(target, out, tname, tool, diff_path)
    return "MISMATCH", f"differs from {tname}", write_diff(target, out, tname, tool, diff_path)
