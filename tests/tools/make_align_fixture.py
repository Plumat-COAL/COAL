"""Writes tests/fixtures/align/reference_scores.json: protein pairs and their optimal global-alignment scores from
Biopython's PairwiseAligner and, where EMBOSS `needle` is on PATH, from needle, under the scoring settings
tests/test_align.py checks coal-cpp-align and coal.align against (so CI checks them without either tool).

    python tests/tools/make_align_fixture.py            # needs Biopython; EMBOSS needle optional

The pairs are random but fixed (seeded): related pairs (substitutions, indels, extensions at both ends), unrelated
pairs, short and empty-adjacent cases, lower case, B/Z/X and '*'.
"""
from __future__ import annotations

import json
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "fixtures", "align", "reference_scores.json")
AA = "ARNDCQEGHILKMFPSTWYV"
# (matrix, gap_open, gap_extend, free_end_gaps)
CONFIGS = [("BLOSUM62", 10, 0.5, True), ("BLOSUM62", 10, 0.5, False), ("BLOSUM62", 14, 4, False),
           ("BLOSUM45", 12, 2, True), ("BLOSUM80", 9, 1, False), ("PAM250", 11, 1, True), ("BLOSUM62", 11, 1, True)]


def mutate(s, rng, sub=0.3, indel=0.05):
    out = []
    for c in s:
        r = rng.random()
        if r < indel / 2:
            continue
        if r < indel:
            out.append("".join(rng.choice(AA) for _ in range(rng.randint(1, 6))))
        out.append(rng.choice(AA) if rng.random() < sub else c)
    return "".join(out)


def make_pairs():
    rng = random.Random(20261111)
    pairs = []
    for k in range(48):
        n = rng.randint(1, 4) if k % 12 == 0 else rng.randint(20, 260)
        alpha = AA + "BZX*" if k % 7 == 0 else AA
        a = "".join(rng.choice(alpha) for _ in range(n))
        b = mutate(a, rng) if k % 3 else "".join(rng.choice(AA) for _ in range(rng.randint(1, 200)))
        if k % 5 == 0:
            b = "".join(rng.choice(AA) for _ in range(rng.randint(0, 40))) + b + \
                "".join(rng.choice(AA) for _ in range(rng.randint(0, 40)))
        if not b:
            b = "M"
        if k % 11 == 0:
            a = a.lower()
        pairs.append([a, b])
    return pairs


def needle_score(a, b, cfg, needle, tmp):
    mat, o, e, free = cfg
    if "*" in a + b or any(c in "BZX" for c in a.upper() + b.upper()):
        return None  # EMBOSS reads these differently; Biopython is the reference for them
    fa, fb, out = (os.path.join(tmp, x) for x in ("a.fa", "b.fa", "o.txt"))
    with open(fa, "w") as fh:
        fh.write(f">a\n{a.upper()}\n")
    with open(fb, "w") as fh:
        fh.write(f">b\n{b.upper()}\n")
    cmd = [needle, "-asequence", "fasta::" + fa, "-bsequence", "fasta::" + fb, "-gapopen", str(o), "-gapextend", str(e), "-datafile",
           "E" + mat, "-outfile", out, "-sprotein1", "-sprotein2", "-auto"]
    if not free:
        cmd += ["-endweight", "-endopen", str(o), "-endextend", str(e)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(f"needle failed on {a!r} {b!r}: {r.stderr}")
    with open(out) as fh:
        return float(re.search(r"# Score: ([-\d.]+)", fh.read()).group(1))


def main():
    from Bio import __version__ as bio_version
    from Bio.Align import PairwiseAligner, substitution_matrices
    pairs = make_pairs()
    needle = shutil.which("needle")
    needle_version = ""
    if needle:
        r = subprocess.run([needle, "-version"], capture_output=True, text=True)
        needle_version = (r.stdout + r.stderr).strip().splitlines()[-1]
    data = {"generated_by": "tests/tools/make_align_fixture.py", "biopython": bio_version,
            "emboss_needle": needle_version, "pairs": pairs, "configs": []}
    with tempfile.TemporaryDirectory() as tmp:
        for cfg in CONFIGS:
            mat, o, e, free = cfg
            al = PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load(mat),
                                 open_gap_score=-o, extend_gap_score=-e)
            if free:
                al.end_gap_score = 0
            bio = [al.score(a.upper(), b.upper()) for a, b in pairs]
            nd = [needle_score(a, b, cfg, needle, tmp) for a, b in pairs] if needle and mat != "PAM250" else None
            if nd is not None:
                for (a, b), x, y in zip(pairs, bio, nd):
                    if y is not None and abs(x - y) > 0.051:
                        print(f"Biopython and needle disagree under {cfg}: {x} vs {y} ({a!r}, {b!r})")
            data["configs"].append({"matrix": mat, "gap_open": o, "gap_extend": e, "free_end_gaps": free,
                                    "biopython": bio, "needle": nd})
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(data, fh, indent=1)
        fh.write("\n")
    print(f"wrote {OUT}: {len(pairs)} pairs x {len(CONFIGS)} settings"
          + (f", needle {needle_version}" if needle else ", no needle"))


if __name__ == "__main__":
    main()
