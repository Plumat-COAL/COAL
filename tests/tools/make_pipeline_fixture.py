"""Synthetic coding sequences for the `coal pipeline` tests and validation (tests/test_dnds_pipeline.py,
docs/pipeline.md "Validation"). Everything is generated from a seed: no real sequence data.

  python tests/tools/make_pipeline_fixture.py                 # rewrites tests/fixtures/pipeline/ (committed)
  python tests/tools/make_pipeline_fixture.py --big OUT_DIR   # the validation set: adds a 23,000-codon gene to the
                                                              # pairs and families of 30 and 300 members

Pairs (query_cds.fasta / subject_cds.fasta): ordinary orthologs g01..gNN / h01..hNN, plus the edge cases
  stop_q/stop_s   an internal stop codon in the query        frame_q/frame_s  subject CDS not a multiple of 3
  amb_q/amb_s     ambiguous bases (N, R, Y) in both          same_q/same_s    identical sequences
  rna_q/rna_s     lower case and U in the query              nohit_q, nohit_s unrelated genes (no hit)
  bad_q           a character that is no nucleotide code     empty_q          an empty record
  g01 (again)     a repeated ID in the query file            short_s          two bases
Families (family_cds.fasta + family_clusters.tsv, DIAMOND's two-column format: representative, member):
  famA (6 members), famB (12), famC (2), famD (1, a singleton), and famB holds one member without a usable CDS.
"""
from __future__ import annotations

import argparse
import os
import random

BASES = "TCAG"
AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
CODE = dict(zip(CODONS, AAS))
SENSE = [c for c in CODONS if CODE[c] != "*"]
SYN = {c: [d for d in SENSE if CODE[d] == CODE[c] and d != c] for c in SENSE}


def random_cds(rng: random.Random, n_codons: int) -> list[str]:
    return ["ATG"] + [rng.choice(SENSE) for _ in range(n_codons - 1)]


def evolve(rng: random.Random, codons: list[str], p_syn: float, p_non: float, p_indel: float = 0.0) -> list[str]:
    out = []
    for k, c in enumerate(codons):
        if k and p_indel and rng.random() < p_indel:
            if rng.random() < 0.5:
                continue                                                   # deletion
            out += [rng.choice(SENSE) for _ in range(rng.randint(1, 3))]   # insertion
        r = rng.random()
        if k and r < p_syn and SYN[c]:
            c = rng.choice(SYN[c])
        elif k and r < p_syn + p_non:
            c = rng.choice(SENSE)
        out.append(c)
    return out


def cds(codons: list[str], stop: str = "TAA") -> str:
    return "".join(codons) + stop


def write_fasta(path: str, recs) -> None:
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        for name, seq in recs:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")
            if not seq:
                fh.write("\n")


def make_pairs(rng: random.Random, n_plain: int, big: bool):
    q, s = [], []
    for k in range(1, n_plain + 1):
        anc = random_cds(rng, rng.randint(60, 180))
        q.append((f"g{k:02d}", cds(evolve(rng, anc, 0.10, 0.03, 0.01))))
        s.append((f"h{k:02d}", cds(evolve(rng, anc, 0.10, 0.03, 0.01), rng.choice(["TAA", "TAG", "TGA"]))))
    anc = random_cds(rng, 120)
    a, b = evolve(rng, anc, 0.1, 0.03), evolve(rng, anc, 0.1, 0.03)
    a[40] = "TAG"
    q.append(("stop_q", cds(a)))
    s.append(("stop_s", cds(b)))
    anc = random_cds(rng, 100)
    q.append(("frame_q", cds(evolve(rng, anc, 0.1, 0.03))))
    s.append(("frame_s", cds(evolve(rng, anc, 0.1, 0.03))[:-2]))            # 3n + 1 bases
    anc = random_cds(rng, 110)
    a, b = "".join(evolve(rng, anc, 0.1, 0.03)), "".join(evolve(rng, anc, 0.1, 0.03))
    a = a[:30] + "NNN" + a[33:60] + "R" + a[61:]
    b = b[:91] + "Y" + b[92:150] + "N" + b[151:]
    q.append(("amb_q", a + "TAA"))
    s.append(("amb_s", b + "TGA"))
    same = cds(random_cds(rng, 90))
    q.append(("same_q", same))
    s.append(("same_s", same))
    anc = random_cds(rng, 80)
    q.append(("rna_q", cds(evolve(rng, anc, 0.1, 0.03)).lower().replace("t", "u")))
    s.append(("rna_s", cds(evolve(rng, anc, 0.1, 0.03))))
    q.append(("nohit_q", cds(random_cds(rng, 100))))
    s.append(("nohit_s", cds(random_cds(rng, 100))))
    q.append(("bad_q", "ATGGCT-GCTAAAGCT!GCTTAA"))
    q.append(("empty_q", ""))
    q.append(("g01", cds(random_cds(rng, 50))))                             # a repeated ID
    s.append(("short_s", "AT"))
    if big:
        anc = random_cds(rng, 23000)
        q.append(("big_q", cds(evolve(rng, anc, 0.08, 0.02, 0.0005))))
        s.append(("big_s", cds(evolve(rng, anc, 0.08, 0.02, 0.0005))))
    return q, s


def make_families(rng: random.Random, sizes: dict[str, int]):
    recs, clusters = [], []
    for fam, n in sizes.items():
        anc = random_cds(rng, rng.randint(90, 160))
        for k in range(1, n + 1):
            name = f"{fam}_{k:03d}"
            recs.append((name, cds(evolve(rng, anc, 0.12, 0.04, 0.004))))
            clusters.append((f"{fam}_001", name))
    return recs, clusters


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--big", metavar="OUT_DIR")
    ap.add_argument("--seed", type=int, default=20261111)
    a = ap.parse_args()
    out = a.big or os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "fixtures", "pipeline")
    os.makedirs(out, exist_ok=True)
    rng = random.Random(a.seed)
    q, s = make_pairs(rng, 12, bool(a.big))
    write_fasta(os.path.join(out, "query_cds.fasta"), q)
    write_fasta(os.path.join(out, "subject_cds.fasta"), s)
    sizes = {"famA": 6, "famB": 12, "famC": 2, "famD": 1}
    if a.big:
        sizes.update({"famE": 30, "famF": 300})
    recs, clusters = make_families(rng, sizes)
    recs = [(n, "ATGNNQQ" if n == "famB_007" else sq) for n, sq in recs]   # one member without a usable CDS
    clusters.append(("famB_001", "famB_ghost"))                            # and one that has no CDS record at all
    write_fasta(os.path.join(out, "family_cds.fasta"), recs)
    with open(os.path.join(out, "family_clusters.tsv"), "w", encoding="utf-8", newline="\n") as fh:
        for rep, mem in clusters:
            fh.write(f"{rep}\t{mem}\n")
    print(f"{len(q)} query, {len(s)} subject, {len(recs)} family CDS -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
