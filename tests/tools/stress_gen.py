#!/usr/bin/env python3
"""Generate stress-test inputs for COAL together with the true codon alignment.

Writes into OUTDIR: the protein alignment (`pep.fasta` or `pep.aln`), the coding DNA (`nuc.fasta`, or
`nuc_0001.fasta` ... with --files), and `truth.fasta`: the codon alignment COAL must write with `-output fasta`
(60 columns per line), built from the same codons. Every DNA sequence is built from its protein row (or the other way
round), so the truth is known without running any aligner. The standard genetic code is a copy of the NCBI table
kept in this script; other codes are taken from coal.tables (they are checked against PAL2NAL by the stored goldens).

    python tests/tools/stress_gen.py /tmp/c1 --nseq 100000 --ncodons 400 --gaps
    python tests/tools/stress_gen.py /tmp/c2 --nseq 1 --ncodons 1000000 --frameshift 999900:1 --utr5 50 --stop
    python tests/tools/stress_gen.py /tmp/c3 --nseq 3 --ncodons 200000 --fmt clustal --marks --gaps

Many sequences come from a pool of --pool distinct variants (record k uses variant k mod pool), so a million
sequences are written in seconds; IDs are always unique. See --help for the odd-format and mismatch options.
Not collected by unittest; standard library only (plus this repo's src/ for the non-standard genetic codes).
"""
from __future__ import annotations

import argparse
import gzip
import os
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent

# NCBI standard code (table 1), bases in TCAG order, kept here on purpose (independent of coal.tables).
_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
STANDARD = {a + b + c: _AAS[16 * i + 4 * j + k]
            for i, a in enumerate(_BASES) for j, b in enumerate(_BASES) for k, c in enumerate(_BASES)}


def code_table(table: int) -> dict[str, str]:
    """codon -> amino acid ('*' for stops) of NCBI table `table`."""
    if table == 1:
        return dict(STANDARD)
    sys.path.insert(0, str(ROOT / "src"))
    from coal.tables import get_code  # noqa: E402
    return dict(get_code(table).code)


def wrap60(s: str) -> str:
    return "\n".join(s[i:i + 60] for i in range(0, len(s), 60)) + "\n" if s else "\n"


class Variant:
    __slots__ = ("prow", "dna", "truth")

    def __init__(self, prow: str, dna: str, truth: str):
        self.prow, self.dna, self.truth = prow, dna, truth


def make_variant(rng: random.Random, a) -> Variant:
    """One aligned row: (protein row, coding DNA, true codon row). Columns are width 3 (frameshift digits 1-2)."""
    table = code_table(a.table)
    by_aa: dict[str, list[str]] = {}
    for cod, aa in table.items():
        by_aa.setdefault(aa, []).append(cod)
    sense = sorted(c for c, aa in table.items() if aa != "*")
    L = a.ncodons
    if a.repeat == "poly":
        prot = ["A"] * L
    elif a.repeat == "period3":
        prot = [("G", "P", "P")[i % 3] for i in range(L)]
    else:
        codons = rng.choices(sense, k=L)
        prot = [table[c] for c in codons]
        if L and a.start_met:
            prot[0] = "M"
    if a.repeat != "none" or a.start_met or a.mismatch in ("unrelated",) or a.lower or a.rna or a.iupac:
        codons = [rng.choice(by_aa[x]) for x in prot]
    if L and a.start_met and a.repeat == "none":
        codons[0] = "ATG" if "ATG" in by_aa.get("M", []) else by_aa["M"][0]
    if a.allx:
        prot = ["X"] * L
    # gaps: a set of columns (never the first) where this row has no residue and no codon
    gapcols: set[int] = set()
    if a.gaps and L > 1:
        k = int(round((L - 1) * a.gap_frac))
        gapcols = set(rng.sample(range(1, L), k)) if k else set()
    # frameshift digits: column -> digit
    shifts: dict[int, int] = {}
    for spec in a.frameshift or ():
        col, d = spec.split(":")
        col, d = int(col), int(d)
        if 0 <= col < L and d in (1, 2):
            shifts[col] = d
            gapcols.discard(col)
    if a.mismatch == "late" and L:
        last = L - 1
        while last in gapcols and last > 0:
            last -= 1
        aa = prot[last]
        prot[last] = "W" if aa != "W" else "C"  # the codon stays: one wrong residue at the end
    if a.repeat != "none" and a.repeat_prefix and L:
        # the gene itself ends in a residue outside the repeat (and its codon), so it fits exactly only once
        prot[L - 1] = "W"
        codons[L - 1] = by_aa["W"][0]
    if a.mismatch == "unrelated":
        codons = rng.choices(sense, k=L)
    if a.iupac:
        # third base N where every reading gives the same residue (4-fold degenerate families)
        for i, c in enumerate(codons):
            if i in gapcols or i in shifts:
                continue
            fam = {table[c[:2] + b] for b in _BASES}
            if len(fam) == 1 and rng.random() < a.iupac:
                codons[i] = c[:2] + "N"
    prow, dna, truth = [], [], []
    gapch = "." if a.dot_gap else "-"
    for i in range(L):
        if i in gapcols:
            prow.append(gapch)
            truth.append("---")
            continue
        if i in shifts:
            d = shifts[i]
            bases = "".join(rng.choice("ACGT") for _ in range(d))
            prow.append(str(d))
            dna.append(bases)
            truth.append(bases + "-" * (3 - d))
            continue
        prow.append(prot[i])
        dna.append(codons[i])
        truth.append(codons[i])
    seq = "".join(dna)
    if a.repeat != "none" and a.repeat_prefix:
        # extra copies of the repeat BEFORE the gene (the first residue's codons, so a terminal mismatch residue
        # is not copied): the gene still fits exactly, further right, after a long low-complexity stretch
        unit = [prot[0]] * L
        seq = "".join("".join(rng.choice(by_aa[x]) for x in unit) for _ in range(a.repeat_prefix)) + seq
    if a.utr5:
        seq = "".join(rng.choice("ACGT") for _ in range(a.utr5)) + seq
    if a.stop:
        seq += "TAA"
    if a.utr3:
        seq += "".join(rng.choice("ACGT") for _ in range(a.utr3))
    t = "".join(truth)
    if a.lower:
        seq, t = seq.lower(), t.lower()
    if a.rna:
        seq, t = seq.replace("T", "U").replace("t", "u"), t.replace("T", "U").replace("t", "u")
    return Variant("".join(prow), seq, t)


def ids_for(k: int, a) -> tuple[str, str, str, str]:
    """(protein ID, protein header rest, DNA ID, DNA header rest) of record k."""
    base = f"seq{k}"
    m = a.ids
    if m == "same":
        return base, "", base, ""
    if m == "version":
        return f"g{k}", "", f"g{k}.1", ""
    if m == "case":
        return base, "", base.upper(), ""
    if m == "half":
        return base, "", (base if k % 2 else f"other{k}"), ""
    if m == "none":
        return base, "", f"d{k}", ""
    if m == "long":
        tag = f"{k:08d}"
        return tag + "x" * (a.id_len - len(tag)), "", tag + "x" * (a.id_len - len(tag)), ""
    if m == "spaces":
        return base, " some description words", base, " other description here"
    if m == "tabs":
        return base, "\tdesc\tmore", base, "\tdesc"
    if m == "unicode-desc":
        return base, " Krüppel-like factor αβ", base, " Krüppel αβ"
    if m == "nonascii-id":
        return f"gène{k}", "", f"gène{k}", ""
    if m == "empty-name":
        return ("" if k == 1 else base), "", base, ""
    raise SystemExit(f"unknown --ids {m}")


def eol(a) -> str:
    return "\r\n" if a.crlf else "\r" if a.cr else "\n"


def open_text(path, a):
    fh = open(path, "w", newline="", encoding="utf-8")
    if a.bom:
        fh.write("﻿")
    return fh


def write_fasta_pep(path, a, pool, nseq, nl):
    w = a.width
    with open_text(path, a) as fh:
        for k in range(nseq):
            pid, prest, _, _ = ids_for(k, a)
            v = pool[k % len(pool)]
            row = v.prow
            body = (row + nl) if not w else nl.join(row[i:i + w] for i in range(0, len(row), w)) + nl
            if a.blank_lines:
                body = body.replace(nl, nl + nl, 1) + nl
            fh.write(f">{pid}{prest}{nl}{body}")


def write_clustal_pep(path, a, pool, nseq, nl, mask):
    names = [ids_for(k, a)[0] for k in range(nseq)]
    idw = max(len(n) for n in names) + 2
    ncol = len(pool[0].prow)
    blk = a.clustal_block or 60
    with open_text(path, a) as fh:
        fh.write("CLUSTAL W (1.83) multiple sequence alignment" + nl + nl + nl)
        for b in range(0, ncol, blk):
            for k in range(nseq):
                fh.write(names[k].ljust(idw) + pool[k % len(pool)].prow[b:b + blk] + nl)
            if mask is not None:
                fh.write(" " * idw + mask[b:b + blk].rstrip() + nl)
            fh.write(nl)


def dna_body(seq: str, a, nl: str) -> str:
    w = a.width
    if a.blank_lines:
        w = w or 60
    body = (seq + nl) if not w else nl.join(seq[i:i + w] for i in range(0, len(seq), w)) + nl
    if a.blank_lines:
        body = body.replace(nl, nl + nl, 2) + nl
    return body


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--nseq", type=int, default=10)
    ap.add_argument("--ncodons", type=int, default=400)
    ap.add_argument("--pool", type=int, default=2000, help="distinct sequence variants (default 2000)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--table", type=int, default=1, help="genetic code used to build the DNA (default 1)")
    ap.add_argument("--no-start-met", dest="start_met", action="store_false", help="do not force a first Met/ATG")
    ap.add_argument("--gaps", action="store_true", help="about --gap-frac of the columns of each row are gaps")
    ap.add_argument("--gap-frac", type=float, default=0.05)
    ap.add_argument("--fmt", choices=("fasta", "clustal"), default="fasta", help="protein alignment format")
    ap.add_argument("--clustal-block", type=int, default=0, help="CLUSTAL block width (default 60)")
    ap.add_argument("--marks", action="store_true", help="CLUSTAL '#' block marks under about half of the columns")
    ap.add_argument("--width", type=int, default=60, help="sequence line width in FASTA files (0 = one line)")
    ap.add_argument("--files", type=int, default=1, help="split the DNA over this many files (round robin)")
    ap.add_argument("--shuffle-dna", action="store_true", help="DNA records in random order (pairing by ID)")
    ap.add_argument("--reverse-dna", action="store_true", help="DNA records in reverse order")
    ap.add_argument("--ids", default="same",
                    choices=("same", "version", "case", "half", "none", "long", "spaces", "tabs", "unicode-desc",
                             "nonascii-id", "empty-name"))
    ap.add_argument("--id-len", type=int, default=10000, help="ID length for --ids long")
    ap.add_argument("--utr5", type=int, default=0, help="random bases before the gene")
    ap.add_argument("--utr3", type=int, default=0, help="random bases after the gene (after the stop if --stop)")
    ap.add_argument("--stop", action="store_true", help="TAA after the gene")
    ap.add_argument("--frameshift", action="append", metavar="COL:D", help="frameshift digit D (1|2) at column COL (0-based), every row")
    ap.add_argument("--mismatch", choices=("none", "late", "unrelated"), default="none",
                    help="late: the last residue is wrong (codon kept); unrelated: random DNA of the same length")
    ap.add_argument("--repeat", choices=("none", "poly", "period3"), default="none",
                    help="low-complexity protein: poly-A, or (GPP)n")
    ap.add_argument("--repeat-prefix", type=int, default=0,
                    help="with --repeat: this many extra copies of the repeat before the gene in the DNA")
    ap.add_argument("--crlf", action="store_true")
    ap.add_argument("--cr", action="store_true", help="old Mac line endings (lone CR)")
    ap.add_argument("--bom", action="store_true", help="UTF-8 BOM at the start of every file")
    ap.add_argument("--blank-lines", action="store_true", help="blank lines inside and between records")
    ap.add_argument("--lower", action="store_true", help="lower-case DNA")
    ap.add_argument("--rna", action="store_true", help="RNA (U for T)")
    ap.add_argument("--iupac", type=float, default=0.0, help="fraction of 4-fold degenerate codons whose third base is N")
    ap.add_argument("--allx", action="store_true", help="protein rows made entirely of X")
    ap.add_argument("--dot-gap", action="store_true", help="'.' as the gap character")
    ap.add_argument("--empty-record", type=int, default=-1, help="DNA record k is empty (header only)")
    ap.add_argument("--gzip", action="store_true", help="also write nuc.fasta.gz")
    ap.add_argument("--no-truth", dest="truth", action="store_false")
    a = ap.parse_args(argv)
    if a.ncodons < 1 or a.nseq < 1:
        raise SystemExit("--nseq and --ncodons must be >= 1")
    rng = random.Random(a.seed)
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    npool = min(a.pool, a.nseq)
    pool = [make_variant(rng, a) for _ in range(npool)]
    nl = eol(a)
    ncol = len(pool[0].prow)
    mask = None
    if a.fmt == "clustal" and a.marks:
        mask = "".join("#" if rng.random() < 0.5 else " " for _ in range(ncol))
        (out / "marks.txt").write_text(mask + "\n")
    # protein alignment
    if a.fmt == "clustal":
        pep = out / "pep.aln"
        write_clustal_pep(pep, a, pool, a.nseq, nl, mask)
    else:
        pep = out / "pep.fasta"
        write_fasta_pep(pep, a, pool, a.nseq, nl)
    # DNA
    order = list(range(a.nseq))
    if a.shuffle_dna:
        rng.shuffle(order)
    if a.reverse_dna:
        order.reverse()
    nuc_paths = []

    def record(k):
        _, _, did, drest = ids_for(k, a)
        seq = "" if k == a.empty_record else pool[k % npool].dna
        return f">{did}{drest}{nl}{dna_body(seq, a, nl) if seq else ''}"

    gz = gzip.open(out / "nuc.fasta.gz", "wt", newline="", encoding="utf-8") if a.gzip else None
    if a.files == 1:
        nuc_paths.append(out / "nuc.fasta")
        with open_text(out / "nuc.fasta", a) as fh:
            for k in order:
                rec = record(k)
                fh.write(rec)
                if gz is not None:
                    gz.write(rec)
    else:
        # one file at a time (round robin over the records), so many thousand files do not need as many open handles
        for i in range(a.files):
            p = out / f"nuc_{i:05d}.fasta"
            nuc_paths.append(p)
            with open_text(p, a) as fh:
                for n in range(i, len(order), a.files):
                    rec = record(order[n])
                    fh.write(rec)
                    if gz is not None:
                        gz.write(rec)
    if gz is not None:
        gz.close()
    # truth: COAL -output fasta (60 columns), protein order and IDs
    if a.truth:
        with open(out / "truth.fasta", "w", newline="\n", encoding="utf-8") as fh:
            for k in range(a.nseq):
                pid = ids_for(k, a)[0]
                fh.write(f">{pid}\n{wrap60(pool[k % npool].truth)}")
    sizes = {p.name: os.path.getsize(p) for p in [pep, *nuc_paths]}
    total = sum(sizes.values())
    print(f"{out}: {a.nseq} x {a.ncodons} codons, {ncol} columns, pool {npool}; input {total / 1e6:.1f} MB "
          f"({', '.join(f'{k} {v / 1e6:.1f} MB' for k, v in list(sizes.items())[:3])}{', ...' if len(sizes) > 3 else ''})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
