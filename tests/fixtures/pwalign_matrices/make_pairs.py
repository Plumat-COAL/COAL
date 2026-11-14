"""The protein pairs of tests/fixtures/pwalign_matrices (pairs.fasta): real orthologous proteins, from the protein
alignments of COAL pipeline runs (--keep-protein-alignments) on public genomes, ungapped again:

    S. cerevisiae / C. glabrata, S. cerevisiae / K. lactis, A. thaliana / A. lyrata, human / mouse,
    E. coli K-12 / Sakai, Salmonella / Shigella, S. cerevisiae / C. glabrata mitochondria

Pairs of 20-260 residues, spread over the range of protein identity of each run; every pair with an X or an internal
'*' in the runs is taken. Then, so that every letter pwalign reads differently from the NCBI files is in the
fixture, some residues are changed (deterministically): B, Z, X, J, U, O and '*' (pwalign's matrices have rows for
B, J, Z, X and '*'; U and O are scored as X by COAL; see make_expected.R for what pwalign does with them).

    python make_pairs.py RUN_DIR... > pairs.fasta      (each RUN_DIR holds protein_alignments.fasta)
"""
import random
import sys

rng = random.Random(20261114)
EDITS = "BZXJUO*"


def read_fasta(path):
    out = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                out.append([line[1:].split()[0], ""])
            elif out:
                out[-1][1] += line.strip()
    return out


def identity(a, b):
    cols = [(x, y) for x, y in zip(a, b) if x != "-" or y != "-"]
    return sum(x == y for x, y in cols) / max(1, len(cols))


pairs = []
for run in sys.argv[1:]:
    rec = read_fasta(run.rstrip("/") + "/protein_alignments.fasta")
    cand, special = [], []
    for k in range(0, len(rec) - 1, 2):
        (na, ga), (nb, gb) = rec[k], rec[k + 1]
        a, b = ga.replace("-", ""), gb.replace("-", "")
        inner = a[:-1] + b[:-1]
        if ("X" in a + b or "*" in inner) and max(len(a), len(b)) <= 600:
            special.append((na, nb, a, b))
        elif 20 <= min(len(a), len(b)) and max(len(a), len(b)) <= 260:
            cand.append((identity(ga, gb), na, nb, a, b))
    cand.sort()
    take = min(14, len(cand))
    idx = sorted({round(i * (len(cand) - 1) / max(1, take - 1)) for i in range(take)})
    pairs += special[:4] + [c[1:] for i, c in enumerate(cand) if i in idx]

for k, (na, nb, a, b) in enumerate(pairs):
    if k % 3 == 1:
        a, b = list(a), list(b)
        for _ in range(rng.randint(1, 3)):
            s = rng.choice((a, b))
            s[rng.randrange(len(s))] = rng.choice(EDITS)
        a, b = "".join(a), "".join(b)
        na, nb = na + "_edited", nb + "_edited"
    sys.stdout.write(f">p{k + 1}_a {na}\n{a}\n>p{k + 1}_b {nb}\n{b}\n")
