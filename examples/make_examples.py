"""Generate the small example datasets in examples/ (fixed seed, no PAL2NAL needed).

    python examples/make_examples.py

Every dataset is correct by construction: we invent an ancestral protein, derive
related proteins from it (substitutions and deletions), write the protein
alignment directly from that history, and back-translate each protein into DNA
with randomly chosen synonymous codons. So we know the right codon alignment
without running any aligner. The data are synthetic, not real genes.

Datasets (each folder has protein.aln.fasta, dna.fasta and README.md):
  clean            nuclear genes, standard genetic code (1)
  mito             vertebrate-mitochondrial-style genes that only fit code 2
  messy_ids        protein IDs like AT1G01010.1, DNA IDs like AT1G01010, DNA in another order
  one_broken_gene  19 good genes plus one with a frameshift (1 base deleted)

This script is self-contained on purpose (it does not import coal), so the
examples are an independent check of the engine.
"""
from __future__ import annotations

import os
import random
from itertools import product

HERE = os.path.dirname(os.path.abspath(__file__))
SEED = 20261017

BASES = "TCAG"
CODONS = ["".join(c) for c in product(BASES, repeat=3)]
# NCBI translation tables (same order as CODONS).
CODE1 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE2 = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG"
AMINO = "ACDEFGHIKLMNPQRSTVWY"


def codon_map(table: str) -> dict[str, list[str]]:
    m: dict[str, list[str]] = {}
    for c, a in zip(CODONS, table):
        m.setdefault(a, []).append(c)
    return m


def wrap(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def write_fasta(path: str, records: list[tuple[str, str]]) -> None:
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        for name, seq in records:
            fh.write(f">{name}\n{wrap(seq)}\n")


def make_family(rng: random.Random, n: int, length: int, sub=0.10, dele=0.05, alphabet=AMINO):
    """Return n aligned protein rows (with '-' gaps) derived from one ancestor."""
    anc = "M" + "".join(rng.choice(alphabet) for _ in range(length - 1))
    rows = []
    for _ in range(n):
        row = [anc[0]]
        for a in anc[1:]:
            r = rng.random()
            if r < dele:
                row.append("-")
            elif r < dele + sub:
                row.append(rng.choice(alphabet))
            else:
                row.append(a)
        rows.append("".join(row))
    # drop columns that are gaps in every row (a real aligner would not emit them)
    keep = [i for i in range(length) if any(r[i] != "-" for r in rows)]
    return ["".join(r[i] for i in keep) for r in rows]


def back_translate(rng: random.Random, protein: str, cmap: dict[str, list[str]], stop=True,
                   prefer: dict[str, str] | None = None) -> str:
    """Random synonymous codons; `prefer` maps aa -> codon used half the time."""
    out = []
    for a in protein:
        if prefer and a in prefer and rng.random() < 0.5:
            out.append(prefer[a])
        else:
            out.append(rng.choice(cmap[a]))
    if stop:
        out.append(rng.choice(cmap["*"]))
    return "".join(out)


def ungap(s: str) -> str:
    return s.replace("-", "")


def write_dataset(name: str, proteins: list[tuple[str, str]], dnas: list[tuple[str, str]], readme: str):
    d = os.path.join(HERE, name)
    os.makedirs(d, exist_ok=True)
    write_fasta(os.path.join(d, "protein.aln.fasta"), proteins)
    write_fasta(os.path.join(d, "dna.fasta"), dnas)
    with open(os.path.join(d, "README.md"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(readme.strip() + "\n")


def clean(rng):
    cmap = codon_map(CODE1)
    rows = make_family(rng, 8, 150)
    names = [f"gene{i + 1:02d}" for i in range(len(rows))]
    prots = list(zip(names, rows))
    dnas = [(n, back_translate(rng, ungap(r), cmap)) for n, r in prots]
    write_dataset("clean", prots, dnas, """
# clean: nuclear genes, standard genetic code

Eight related protein-coding genes (synthetic, 140-144 codons each, aligned in 150 columns) with a
ready-made protein alignment and their unaligned coding DNA, each ending in a stop codon. Everything fits the
standard genetic code (NCBI code 1), the IDs match exactly, and there are no frameshifts,
so COAL should report every sequence as OK. Use it as a first test:
`python -m coal examples/clean/protein.aln.fasta examples/clean/dna.fasta`.
""")


def mito(rng):
    cmap = codon_map(CODE2)
    # Vertebrate mitochondrial code: TGA = Trp (a stop in code 1), ATA = Met (Ile in code 1).
    # AGA/AGG are stops in code 2, so they never appear inside a gene here.
    genes = ["COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND4", "ATP6"]
    species = ["Homo", "Pan", "Mus", "Bos", "Gallus", "Danio"]
    prots, dnas = [], []
    for g in genes[:4]:
        rows = make_family(rng, len(species), 110)
        for sp, r in zip(species, rows):
            if "W" not in r[1:]:
                i = rng.randrange(5, len(r) - 5)
                r = r[:i] + "W" + r[i + 1:]
            name = f"{sp}_{g}"
            prots.append((name, r))
            dna = back_translate(rng, ungap(r), cmap, prefer={"W": "TGA", "M": "ATA"})
            # make sure every gene really needs code 2: at least one internal TGA for Trp
            p = ungap(r)
            if not any(dna[3 * k:3 * k + 3] == "TGA" for k in range(1, len(p)) if p[k] == "W"):
                k = next(k for k in range(1, len(p)) if p[k] == "W")
                dna = dna[:3 * k] + "TGA" + dna[3 * k + 3:]
            dnas.append((name, dna))
    write_dataset("mito", prots, dnas, """
# mito: vertebrate mitochondrial genes (genetic code 2)

Twenty-four synthetic mitochondrial-style genes (COX1, COX2, COX3 and CYTB from six
vertebrates, about 110 codons each). Mitochondria read DNA slightly differently from the
nucleus: TGA means tryptophan (W) instead of "stop", and ATA means methionine (M). Every gene
here uses TGA for at least one W, so with the default standard code (1) COAL finds codons
that don't match their amino acid. Re-running with the vertebrate mitochondrial code fixes
it: `python -m coal examples/mito/protein.aln.fasta examples/mito/dna.fasta -codontable 2`.
""")


def messy_ids(rng):
    cmap = codon_map(CODE1)
    rows = make_family(rng, 10, 120)
    loci = [f"AT{c}G{rng.randrange(1, 99) * 1000 + 10:05d}" for c in (1, 1, 2, 2, 3, 3, 4, 4, 5, 5)]
    prots = [(f"{locus}.1", r) for locus, r in zip(loci, rows)]
    dnas = [(locus, back_translate(rng, ungap(r), cmap)) for locus, r in zip(loci, rows)]
    rng.shuffle(dnas)
    write_dataset("messy_ids", prots, dnas, """
# messy_ids: protein and DNA IDs that don't quite match

Ten synthetic Arabidopsis-style genes. The protein alignment uses transcript IDs such as
`AT1G01010.1`, while the DNA file uses gene IDs such as `AT1G01010`, and the DNA records are
in a different order. This is very common when proteins and CDS come from different
downloads. A tool that pairs sequences by exact ID (or by position) will pair the wrong
protein with the wrong DNA. The fix is to match IDs after removing the `.1` version
suffix; once paired correctly, every gene fits the standard genetic code (1).
""")


def one_broken_gene(rng):
    cmap = codon_map(CODE1)
    rows = make_family(rng, 20, 200)
    names = [f"OG0001_sp{i + 1:02d}" for i in range(len(rows))]
    prots = list(zip(names, rows))
    dnas = [(n, back_translate(rng, ungap(r), cmap)) for n, r in prots]
    broken = 13  # the 14th gene
    n, dna = dnas[broken]
    cut = (len(dna) // 2 // 3) * 3 + 1  # delete one base in the middle of a codon
    dnas[broken] = (n, dna[:cut] + dna[cut + 1:])
    write_dataset("one_broken_gene", prots, dnas, f"""
# one_broken_gene: 19 good genes and one frameshift

Twenty synthetic orthologous genes (200 codons). Nineteen are fine. In `{n}` one DNA base
has been deleted halfway along the gene (a frameshift, as caused by a sequencing or
assembly error), so from that point on its codons no longer spell the protein. A good tool
should flag that one gene clearly, explain what's wrong in plain words, and (if asked)
still produce the alignment for the other nineteen.
""")


def main():
    rng = random.Random(SEED)
    clean(rng)
    mito(rng)
    messy_ids(rng)
    one_broken_gene(rng)
    for name in ("clean", "mito", "messy_ids", "one_broken_gene"):
        d = os.path.join(HERE, name)
        size = sum(os.path.getsize(os.path.join(d, f)) for f in os.listdir(d))
        print(f"{name:16s} {size / 1024:5.1f} KB")


if __name__ == "__main__":
    main()
