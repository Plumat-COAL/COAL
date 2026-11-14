"""The synthetic protein pairs of tests/fixtures/pwalign (pairs.fasta). The expected alignments (expected.fasta)
are pwalign::pairwiseAlignment(type = "global", substitutionMatrix = "BLOSUM62")'s, written by make_expected.R.

    python tests/fixtures/pwalign/make_pairs.py && Rscript tests/fixtures/pwalign/make_expected.R

Repeats, runs of one residue and indels make most pairs have several co-optimal alignments."""
import pathlib
import random

AA = "ARNDCQEGHILKMFPSTWYV"
rng = random.Random(20261114)


def mutate(s, sub, indel):
    out = []
    for c in s:
        r = rng.random()
        if r < indel / 2:
            continue                                      # deletion
        if r < indel:
            out.append("".join(rng.choice(AA) for _ in range(rng.randint(1, 4))))   # insertion
        out.append(rng.choice(AA) if rng.random() < sub else c)
    return "".join(out) or "M"


pairs = []
for k in range(80):
    n = rng.randint(5, 160)
    unit = "".join(rng.choice(AA) for _ in range(rng.randint(1, 4)))
    s = "".join(rng.choice(AA) if rng.random() < 0.6 else unit for _ in range(n))[:n]
    t = mutate(s, rng.choice((0.05, 0.2, 0.4)), rng.choice((0.02, 0.08, 0.2)))
    if k % 10 == 3:
        t = t + "*"
    if k % 10 == 7:
        s, t = s.replace("A", "X", 1), t.replace("G", "B", 1)
    pairs.append((s, t))
pairs += [("MAAAKAAA", "MAAKAAAA"), ("GGGGAGGG", "GGGAGGGG"), ("MKVLAAGG", "MKVLAGG"), ("M", "MK"), ("WW", "W")]
with open(pathlib.Path(__file__).with_name("pairs.fasta"), "w", encoding="utf-8", newline="\n") as fh:
    for i, (a, b) in enumerate(pairs, 1):
        fh.write(f">p{i}_a\n{a}\n>p{i}_b\n{b}\n")
