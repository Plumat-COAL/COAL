#!/usr/bin/env python3
"""Check an installed platform wheel of coal-align: are both binaries the wheel's own, and do they work?

    <venv>/bin/python packaging/wheels/check_wheel_binaries.py

Run it with the Python of the environment the wheel was installed into (release.yml does, after
build_platform_wheel.py; tests/test_wheels.py does on the local build). Self-contained: its sequences are in this
file, and it takes COAL_CPP, COAL_CPP_ALIGN, COAL_ENGINE and PYTHONPATH out of the environment and runs with a PATH
that holds no tool, so only what the wheel installed can be found:

1. `coal.dispatch.find_cpp()` and `coal.align.find_align_binary()` are coal/bin/coal-cpp and coal/bin/coal-cpp-align
   inside the installed package, executable, and `--version` is the package's version;
2. `coal.align.align_pairs(..., python_fallback=False)` (coal-cpp-align) gives the alignments of the pure-Python
   kernel;
3. `coal pipeline pairs --pairs ...` (translation, pairwise alignment, codon alignment and Ka/Ks by the bundled
   coal-cpp; no DIAMOND) writes a row per pair and method.

Exit 0 and one line per check, or exit 1 with the reason.
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile

# two "species" of three short genes (the second with substitutions and one codon deleted)
QUERY = {
    "q1": "ATGGCTAAAGGTCTGTTTGATCCGGAATGGCGTAAACTGATTGAAGCTCTGCGTGATCATGGTTATACCGTTTAA",
    "q2": "ATGAGCGAACAGCTGAAAGCGATTCTGGAAAAACATCCGGATGTGAAAGTTCTGAGCGGTGAAGCGTGGTAA",
    "q3": "ATGACCCCGATTCGTCATGATCTGGTTAACGCGTATAAAGGTCAGGAACTGCCGTTTGGTATTGCGTAA",
}
SUBJECT = {
    "s1": "ATGGCAAAGGGTTTGTTCGACCCGGAGTGGCGCAAATTGATCGAAGCACTGCGTGACCACGGCTACACCGTCTAA",
    "s2": "ATGAGTGAACAACTGAAGGCGATCCTGGAAAAACACCCGGACGTTAAAGTGCTGTCTGAAGCTTGGTAA",
    "s3": "ATGACGCCGATCCGTCACGATTTGGTGAATGCCTATAAAGGCCAAGAGCTGCCATTCGGCATCGCGTAA",
}
PAIRS = [("q1", "s1"), ("q2", "s2"), ("q3", "s3")]
PROTEINS = [("MAKGLFDPEWRKLIEALRDHGYTV", "MAKGLFDPEWRKLIEALRDHGYTV"), ("MSEQLKAILEKHPDVKVLSGEAW", "MSEQLKAILEKHPDVKVLSEAW"),
            ("MTPIRHDLVNAYKGQELPFGIA", "MTPIRHDLVNAYKGQELPFGIA")]


def fail(msg: str) -> "None":
    print(f"check_wheel_binaries: {msg}", file=sys.stderr)
    sys.exit(1)


def main() -> int:
    for k in ("COAL_CPP", "COAL_CPP_ALIGN", "COAL_ENGINE", "PYTHONPATH", "COAL_DIAMOND"):
        os.environ.pop(k, None)
    start = os.getcwd()
    with tempfile.TemporaryDirectory(prefix="coal_wheel_check_") as tmp:
        try:
            return check(tmp)
        finally:
            os.chdir(start)   # (Windows cannot remove the folder a process is in)


def check(tmp: str) -> int:
    if True:
        # an empty folder as PATH (kept: the system folder on Windows, where DLL lookup uses PATH)
        empty = os.path.join(tmp, "empty-path")
        os.mkdir(empty)
        os.environ["PATH"] = empty + (os.pathsep + os.path.join(os.environ.get("SystemRoot", r"C:\Windows"), "System32")
                                      if os.name == "nt" else "")
        os.chdir(tmp)   # not the checkout: `import coal` must be the installed package
        here = os.path.dirname(os.path.abspath(__file__))
        sys.path[:] = [p for p in sys.path if p and os.path.abspath(p) != here]
        import coal
        import coal.align as al
        from coal import dispatch

        pkg_bin = os.path.join(os.path.dirname(os.path.abspath(coal.__file__)), "bin")
        ext = ".exe" if os.name == "nt" else ""
        found = {"coal-cpp": dispatch.find_cpp(), "coal-cpp-align": al.find_align_binary()}
        for name, path in found.items():
            want = os.path.join(pkg_bin, name + ext)
            if path is None or os.path.normcase(os.path.realpath(path)) != os.path.normcase(os.path.realpath(want)):
                fail(f"{name}: found {path!r}, expected the wheel's {want}")
            if not os.access(path, os.X_OK):
                fail(f"{path} is not executable")
            r = subprocess.run([path, "--version"], capture_output=True, text=True, timeout=60,
                               stdin=subprocess.DEVNULL)
            words = r.stdout.split()
            if r.returncode != 0 or len(words) < 2 or words[1] != coal.__version__:
                fail(f"{path} --version: exit {r.returncode}, {r.stdout.strip()!r} (package {coal.__version__})")
            print(f"{name}: {path} ({r.stdout.strip()})")

        got = al.align_pairs(PROTEINS, jobs=2, python_fallback=False)
        for (a, b), g in zip(PROTEINS, got):
            w = al.align_pair(a, b)
            if g.error or (g.a, g.b, g.score, g.identities) != (w.a, w.b, w.score, w.identities):
                fail(f"coal-cpp-align and the Python kernel differ for {a} / {b}: {g} against {w}")
        if "-" not in got[1].a + got[1].b:
            fail("the second pair should need a gap")
        print(f"coal.align.align_pairs: {len(got)} pairs by coal-cpp-align, the Python kernel's alignments")

        for fn, recs in (("query.fasta", QUERY), ("subject.fasta", SUBJECT)):
            with open(fn, "w", encoding="utf-8", newline="\n") as fh:
                fh.writelines(f">{k}\n{v}\n" for k, v in recs.items())
        with open("pairs.tsv", "w", encoding="utf-8", newline="\n") as fh:
            fh.writelines(f"{a}\t{b}\n" for a, b in PAIRS)
        r = subprocess.run([sys.executable, "-m", "coal", "pipeline", "pairs", "--query", "query.fasta", "--subject",
                            "subject.fasta", "--pairs", "pairs.tsv", "-o", "out", "-j", "2", "-m", "YN", "-m",
                            "Comeron"], capture_output=True, text=True, timeout=600, stdin=subprocess.DEVNULL)
        if r.returncode != 0:
            fail(f"coal pipeline pairs: exit {r.returncode}\n{r.stdout}{r.stderr}")
        with open(os.path.join("out", "dnds.tsv"), encoding="utf-8") as fh:
            rows = [ln.rstrip("\n").split("\t") for ln in fh]
        head, rows = rows[0], rows[1:]
        st, dn, cc = head.index("status"), head.index("dN"), head.index("codons_compared")
        good = (len(rows) == 2 * len(PAIRS) and all(row[st] in ("ok", "saturated") and int(row[cc]) >= 20 for row in rows)
                and all(row[st] == "ok" and float(row[dn]) > 0 for row in rows if row[0] == "q2"))
        if not good:
            fail(f"coal pipeline pairs: expected {2 * len(PAIRS)} rows of >= 20 codons with status ok or saturated, and "
                 "dN > 0 for q2, got\n" + "\n".join("\t".join(row) for row in [head] + rows))
        print(f"coal pipeline pairs: {len(rows)} rows (YN and Comeron for {len(PAIRS)} pairs), no tool on PATH")
    return 0


if __name__ == "__main__":
    sys.exit(main())
