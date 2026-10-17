"""Turn a simulated Family into a COAL test case: inputs, the true codon alignment, expected warnings.

Edits ("features") applied on top of the simulated truth:

  utr            random 5'/3' flanking DNA (0-300 nt) on every CDS        -> truth unchanged (trimmed)
  stop3          stop codon appended to every CDS (before any 3' UTR)     -> truth unchanged (trimmed)
  protein_stop   extra last column: protein '*', CDS stop codon           -> truth includes the stop
  internal_stop  one in-frame stop ('*') in ~30% of rows                  -> truth includes it
  alt_start      first codon is a non-ATG start codon of the table (still 'M')
  lower / mixed  lowercase / random-case DNA                               -> truth keeps the case
  rna            T -> U                                                     -> truth keeps U
  iupac          ambiguous bases: 3rd-position N/R/Y where synonymous (residue kept), or NNN = 'X'
  frameshift     ~20% of rows get a digit d in {1,2,4,5} for d bases (PAL2NAL convention:
                 the column is widened to the next multiple of 3 for every row)
  mismatch       ~25% of rows get 1-3 codons that don't encode their residue (spaced >= 40
                 residues apart, away from the ends) -> truth = codon placed anyway + a warning
  shuffle        CDS records in random order (pairing by ID)
  rename         CDS IDs differ from protein IDs, same order (pairing by position)
  split          CDS split over two files
  headers        CDS headers carry descriptions after the ID
  wrap           CDS wrapped at 70 columns
  crlf           Windows line endings in both inputs
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from simulate import Family, code, translate

BASES = "ACGT"


def _rand_dna(rng, n):
    return "".join(rng.choice(list(BASES), n)) if n else ""


def _fourfold_prefixes(table):
    c = code(table)
    out = set()
    for a in BASES:
        for b in BASES:
            aas = {c.aa[a + b + x] for x in BASES}
            if len(aas) == 1 and "*" not in aas:
                out.add(a + b)
    return out


def make_case(fam: Family, rng: np.random.Generator, features: set[str]) -> dict:
    table = fam.table
    c = code(table)
    rows = fam.rows()                       # per row: list of codon | None
    n, ncol = len(rows), fam.ncol
    prot = [[translate(x, table) if x else "-" for x in r] for r in rows]
    tok = [list(r) for r in rows]           # DNA token placed in each column (None = gap)
    warnings = []

    def residue_cols(i, lo=1, hi=None):
        hi = ncol if hi is None else hi
        return [j for j in range(lo, hi) if tok[i][j] is not None]

    if "alt_start" in features:
        alts = [s for s in c.starts if s != "ATG" and s not in ("ATT", "ATC", "ATA")]
        for i in range(n):
            tok[i][0] = alts[int(rng.integers(len(alts)))]

    if "internal_stop" in features:
        for i in range(n):
            if rng.random() < 0.3:
                cols = residue_cols(i, 1, ncol - 1)
                if cols:
                    j = cols[int(rng.integers(len(cols)))]
                    tok[i][j] = c.stops[int(rng.integers(len(c.stops)))]
                    prot[i][j] = "*"

    if "iupac" in features:
        ff = _fourfold_prefixes(table)
        for i in range(n):
            for j in residue_cols(i):
                if rng.random() >= 0.03:
                    continue
                cod = tok[i][j]
                k = rng.random()
                if k < 0.4 and cod[:2] in ff:
                    new = cod[:2] + "N"
                elif k < 0.8:
                    amb = "R" if cod[2] in "AG" else "Y"
                    new = cod[:2] + amb
                    if translate(new, table) != prot[i][j]:
                        continue
                else:
                    new = "NNN"
                tok[i][j] = new
                prot[i][j] = translate(new, table)

    if "mismatch" in features:
        sense = c.sense
        for i in range(n):
            if rng.random() >= 0.25:
                continue
            cols = residue_cols(i)
            if len(cols) < 120:
                continue
            k = int(rng.integers(1, 4))
            chosen = []
            for idx in rng.permutation(np.arange(20, len(cols) - 20)):
                if all(abs(int(idx) - x) >= 40 for x in chosen):
                    chosen.append(int(idx))
                if len(chosen) == k:
                    break
            for idx in chosen:
                j = cols[idx]
                aa = prot[i][j]
                if aa in "X*":
                    continue
                choices = [s for s in sense if c.aa[s] != aa]
                tok[i][j] = choices[int(rng.integers(len(choices)))]
                warnings.append({"id": fam.ids[i], "pep_col": j + 1, "residue": aa, "codon": tok[i][j]})

    widths = [3] * ncol
    if "frameshift" in features:
        for i in range(n):
            if rng.random() >= 0.2:
                continue
            cols = residue_cols(i, 2, ncol - 2)
            for j in rng.choice(cols, min(len(cols), int(rng.integers(1, 3))), replace=False) if cols else []:
                j = int(j)
                if prot[i][j] in "X*":
                    continue
                d = int(rng.choice([1, 2, 4, 5]))
                cod = tok[i][j]
                tok[i][j] = cod[:d] if d < 3 else cod + _rand_dna(rng, d - 3)
                prot[i][j] = str(d)
                widths[j] = max(widths[j], -(-d // 3) * 3)

    if "protein_stop" in features:
        widths.append(3)
        for i in range(n):
            tok[i].append(c.stops[int(rng.integers(len(c.stops)))])
            prot[i].append("*")

    # case / RNA apply to every DNA token, so truth and input agree
    def recase(s: str) -> str:
        if "lower" in features:
            s = s.lower()
        elif "mixed" in features:
            s = "".join(ch.lower() if rng.random() < 0.5 else ch for ch in s)
        if "rna" in features:
            s = s.replace("T", "U").replace("t", "u")
        return s

    for i in range(n):
        tok[i] = [recase(t) if t else t for t in tok[i]]

    truth = []
    cds = []
    for i in range(n):
        parts = []
        for j, w in enumerate(widths):
            t = tok[i][j] or ""
            parts.append(t + "-" * (w - len(t)))
        truth.append("".join(parts))
        core = "".join(t for t in tok[i] if t)
        if "stop3" in features:
            core += recase(c.stops[int(rng.integers(len(c.stops)))])
        if "utr" in features:
            core = recase(_rand_dna(rng, int(rng.integers(0, 301)))) + core + recase(_rand_dna(rng, int(rng.integers(0, 301))))
        cds.append(core)

    pep = ["".join(p) for p in prot]
    return {"ids": fam.ids, "pep": pep, "cds": cds, "truth": truth, "warnings": warnings, "table": table}


def write_fasta(path: Path, records, wrap: int | None = 60, eol: str = "\n"):
    with open(path, "w", newline="") as fh:
        for h, s in records:
            fh.write(f">{h}{eol}")
            if wrap:
                for k in range(0, len(s), wrap):
                    fh.write(s[k:k + wrap] + eol)
            else:
                fh.write(s + eol)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def write_case(outdir: Path, case: dict, features: set[str], rng, meta: dict) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    eol = "\r\n" if "crlf" in features else "\n"
    ids = case["ids"]
    write_fasta(outdir / "pep.fasta", zip(ids, case["pep"]), 60, eol)
    nuc_ids = [f"cds_{k + 1}" for k in range(len(ids))] if "rename" in features else list(ids)
    heads = [f"{x} simulated CDS len={len(s)}" if "headers" in features else x for x, s in zip(nuc_ids, case["cds"])]
    recs = list(zip(heads, case["cds"]))
    if "shuffle" in features:
        recs = [recs[k] for k in rng.permutation(len(recs))]
    wrap = 70 if "wrap" in features else None
    nuc_files = ["nuc.fasta"]
    if "split" in features:
        h = len(recs) // 2
        nuc_files = ["nuc1.fasta", "nuc2.fasta"]
        write_fasta(outdir / nuc_files[0], recs[:h], wrap, eol)
        write_fasta(outdir / nuc_files[1], recs[h:], wrap, eol)
    else:
        write_fasta(outdir / "nuc.fasta", recs, wrap, eol)
    write_fasta(outdir / "truth.fasta", zip(ids, case["truth"]), 60)
    if case["warnings"]:
        with open(outdir / "expected_warnings.tsv", "w") as fh:
            fh.write("id\tpep_col\tresidue\tcodon\n")
            for w in case["warnings"]:
                fh.write(f"{w['id']}\t{w['pep_col']}\t{w['residue']}\t{w['codon']}\n")
    opts = ["-codontable", str(case["table"])] if case["table"] != 1 else []
    rec = {
        **meta,
        "pep": "pep.fasta",
        "nuc": nuc_files,
        "truth": "truth.fasta",
        "opts": opts,
        "features": sorted(features),
        "n_seqs": len(ids),
        "aln_cols": len(case["pep"][0]) if ids else 0,
        "cds_nt": sum(len(s) for s in case["cds"]),
        "expected_warnings": len(case["warnings"]),
    }
    (outdir / "meta.json").write_text(json.dumps(rec, indent=1) + "\n")
    return rec
