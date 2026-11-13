"""Generate examples/dnds: two small "species" of coding sequences with known orthologs, for `coal pipeline`,
`coal dnds` and `coal kaks` (fixed seed; synthetic, not real genes).

    python examples/make_dnds_example.py

The simulator. For each of 12 genes an ancestral coding sequence is drawn (codons by the amino-acid and codon
frequencies below, which are typical of eukaryotic proteomes; ATG first, no stop inside). Each species' copy evolves
from it independently, codon by codon, by single-nucleotide substitutions in continuous time (Gillespie's
algorithm) with the rates of the Muse-Gaut / Goldman-Yang family of codon models: 1 for a transversion, KAPPA for
a transition, times the gene's OMEGA when the amino acid changes, 0 into a stop codon. The time on each branch is
the gene's DS / 2, scaled so that one unit is one synonymous substitution per synonymous site at the ancestral
sequence, so the two copies are about DS apart at synonymous sites and about OMEGA * DS at nonsynonymous ones. A
few genes also get an in-frame insertion or deletion in one species, so that the protein alignment has gaps.
The genes differ in length (45 to 480 codons), OMEGA (0.04 to 1.2) and DS (0.08 to 0.9), so a run shows a range of
dN/dS and, with `--min-codons 50`, one flagged pair. Each species also has one gene without an ortholog.

This script does not import coal: the example is an independent input. It writes species_a_cds.fasta,
species_b_cds.fasta and orthologs.tsv (the true pairs, with the parameters each was simulated under). The codon
alignments in codon.fasta are COAL's own output for these files (tests/test_examples.py checks that a run still
writes them; examples/dnds/README.md has the command).
"""
from __future__ import annotations

import os
import random
from itertools import product

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "dnds")
SEED = 20261113
KAPPA = 2.5

BASES = "TCAG"
CODONS = ["".join(c) for c in product(BASES, repeat=3)]
AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = dict(zip(CODONS, AAS))
SENSE = [c for c in CODONS if CODE[c] != "*"]
TRANSITIONS = {("T", "C"), ("C", "T"), ("A", "G"), ("G", "A")}
# amino-acid frequencies (percent, rounded; typical of eukaryotic proteomes)
AA_FREQ = {"A": 7.0, "R": 5.4, "N": 4.4, "D": 5.4, "C": 1.8, "Q": 3.9, "E": 6.5, "G": 6.4, "H": 2.4, "I": 5.3,
           "L": 9.5, "K": 6.3, "M": 2.4, "F": 4.2, "P": 4.8, "S": 9.0, "T": 5.1, "W": 1.2, "Y": 2.9, "V": 6.7}
# (name, codons, omega, dS, in-frame indel in species B: +n insertion / -n deletion of n codons, or 0)
GENES = [("g01", 320, 0.04, 0.35, 0), ("g02", 210, 0.08, 0.60, -3), ("g03", 480, 0.12, 0.25, 0),
         ("g04", 150, 0.18, 0.45, 2), ("g05", 260, 0.25, 0.90, 0), ("g06", 390, 0.33, 0.30, -6),
         ("g07", 180, 0.45, 0.20, 0), ("g08", 300, 0.60, 0.50, 4), ("g09", 120, 0.80, 0.15, 0),
         ("g10", 240, 1.00, 0.40, 0), ("g11", 200, 1.20, 0.08, 0), ("g12", 45, 0.30, 0.30, 0)]


def codon_weights(rng: random.Random) -> dict[str, float]:
    """Codon frequencies: the amino-acid frequency shared between its codons with a mild, fixed codon bias."""
    w = {}
    for aa, f in AA_FREQ.items():
        cs = [c for c in SENSE if CODE[c] == aa]
        bias = [rng.uniform(0.5, 1.5) for _ in cs]
        for c, b in zip(cs, bias):
            w[c] = f * b / sum(bias)
    return w


def neighbours(codon: str, omega: float) -> list[tuple[str, float, bool]]:
    """(codon one substitution away, rate, synonymous?) for every sense neighbour."""
    out = []
    for pos in range(3):
        for b in BASES:
            if b == codon[pos]:
                continue
            new = codon[:pos] + b + codon[pos + 1:]
            if CODE[new] == "*":
                continue
            rate = KAPPA if (codon[pos], b) in TRANSITIONS else 1.0
            syn = CODE[new] == CODE[codon]
            out.append((new, rate if syn else rate * omega, syn))
    return out


def synonymous_scale(codons: list[str]) -> float:
    """Time per unit of dS: with it, one unit of time is one synonymous substitution per synonymous site, where a
    codon's synonymous sites are counted as its synonymous share of the neutral (omega = 1) rate, times 3."""
    syn_rate = site = 0.0
    for c in codons:
        nb = neighbours(c, 1.0)
        total = sum(r for _, r, _ in nb)
        s = sum(r for _, r, syn in nb if syn)
        syn_rate += s
        site += 3.0 * s / total
    return site / syn_rate


def evolve(rng: random.Random, codons: list[str], omega: float, time: float) -> list[str]:
    out = [codons[0]]                      # the start codon stays
    for c in codons[1:]:
        t = 0.0
        while True:
            nb = neighbours(c, omega)
            total = sum(r for _, r, _ in nb)
            t += rng.expovariate(total)
            if t > time:
                break
            x = rng.uniform(0, total)
            for new, r, _ in nb:
                x -= r
                if x <= 0:
                    c = new
                    break
        out.append(c)
    return out


def indel(rng: random.Random, codons: list[str], n: int, weights: dict[str, float]) -> list[str]:
    if n == 0:
        return codons
    at = rng.randint(len(codons) // 4, 3 * len(codons) // 4)
    if n < 0:
        return codons[:at] + codons[at - n:]
    return codons[:at] + rng.choices(list(weights), list(weights.values()), k=n) + codons[at:]


def write_fasta(path: str, recs) -> None:
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        for name, desc, seq in recs:
            fh.write(f">{name} {desc}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


def main() -> None:
    rng = random.Random(SEED)
    weights = codon_weights(rng)
    names, w = list(weights), list(weights.values())
    a_recs, b_recs, truth = [], [], []
    for gene, n, omega, ds, gap in GENES:
        anc = ["ATG"] + rng.choices(names, w, k=n - 1)
        time = ds / 2 * synonymous_scale(anc)
        a = evolve(rng, anc, omega, time)
        b = indel(rng, evolve(rng, anc, omega, time), gap, weights)
        a_recs.append((f"spA_{gene}", f"synthetic gene {gene}, species A", "".join(a) + "TAA"))
        b_recs.append((f"spB_{gene}", f"synthetic gene {gene}, species B", "".join(b) + rng.choice(["TAA", "TGA", "TAG"])))
        truth.append((f"spA_{gene}", f"spB_{gene}", n, omega, ds, gap))
    for recs, sp in ((a_recs, "A"), (b_recs, "B")):   # one gene per species that the other does not have
        orphan = ["ATG"] + rng.choices(names, w, k=160)
        recs.append((f"sp{sp}_orphan", f"synthetic gene without an ortholog, species {sp}", "".join(orphan) + "TAA"))
    rng.shuffle(b_recs)                                # the two files are not in the same order
    os.makedirs(OUT, exist_ok=True)
    write_fasta(os.path.join(OUT, "species_a_cds.fasta"), a_recs)
    write_fasta(os.path.join(OUT, "species_b_cds.fasta"), b_recs)
    with open(os.path.join(OUT, "orthologs.tsv"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("query\tsubject\tcodons\tsimulated_omega\tsimulated_dS\tindel_codons\n")
        for row in truth:
            fh.write("\t".join(str(x) for x in row) + "\n")
    print(f"{OUT}: {len(a_recs)} + {len(b_recs)} CDS, {len(truth)} ortholog pairs")


if __name__ == "__main__":
    main()
