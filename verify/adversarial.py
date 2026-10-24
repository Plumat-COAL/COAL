#!/usr/bin/env python3
"""Phase 2 of the verification: hand-built adversarial cases through PAL2NAL v14, the Python engine and coal-cpp.

    python verify/adversarial.py OUTDIR            # writes OUTDIR/<case>/..., OUTDIR/results.json, prints a table

Every case is run with `-output fasta` (PAL2NAL supports it) unless the case says otherwise. Recorded per tool:
exit status, stdout bytes (sha256 + first bytes), stderr (first lines). Then: Python vs C++ (stdout, stderr, rc
must be identical) and COAL vs PAL2NAL (identical stdout, or both fail, or a deliberate difference to classify).
"""
from __future__ import annotations

import hashlib
import json
import os
import random
import shutil
import stat
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PY = [sys.executable, "-m", "coal"]
CPP = os.environ.get("COAL_CPP_BIN") or str(ROOT / "build" / "coal-cpp")
PAL2NAL = os.environ.get("PAL2NAL") or str(ROOT / "PAL2NAL" / "pal2nal.v14" / "pal2nal.pl")
ENV = dict(os.environ, COAL_ENGINE="python", PYTHONPATH=str(ROOT / "src"), PYTHONIOENCODING="utf-8", PYTHONUTF8="1")
CODE1 = {"A": "GCT", "K": "AAA", "L": "CTG", "D": "GAC", "V": "GTT", "E": "GAA", "G": "GGC", "P": "CCA", "S": "TCT",
         "T": "ACC", "N": "AAC", "Q": "CAG", "R": "CGT", "I": "ATC", "F": "TTC", "Y": "TAC", "H": "CAC", "C": "TGC",
         "W": "TGG", "M": "ATG"}
RC = {"A": "T", "C": "G", "G": "C", "T": "A"}


def gene(rng, n, aas=CODE1):
    picks = [rng.choice(sorted(aas)) for _ in range(n - 1)]
    return "M" + "".join(picks), "ATG" + "".join(aas[a] for a in picks)


def fasta(recs, wrap=0):
    out = []
    for i, s in recs:
        out.append(f">{i}\n")
        out.append(s + "\n" if not wrap else "\n".join(s[k:k + wrap] for k in range(0, len(s), wrap)) + "\n")
    return "".join(out)


def clustal(recs, marks=None):
    w = max(len(i) for i, _ in recs) + 4
    s = "CLUSTAL W (1.83) multiple sequence alignment\n\n\n"
    n = len(recs[0][1])
    for b in range(0, n, 60):
        for i, r in recs:
            s += i.ljust(w) + r[b:b + 60] + "\n"
        if marks is not None:
            s += " " * w + marks[b:b + 60] + "\n"
        s += "\n"
    return s


def write(d: Path, name: str, data, binary=False):
    d.mkdir(parents=True, exist_ok=True)
    p = d / name
    if binary or isinstance(data, bytes):
        p.write_bytes(data if isinstance(data, bytes) else data.encode("utf-8"))
    else:
        p.write_text(data, encoding="utf-8", newline="\n")
    return p


def genes(rng, n, ncod):
    peps, nucs = [], []
    for k in range(n):
        p, s = gene(rng, ncod)
        peps.append((f"g{k}", p))
        nucs.append((f"g{k}", s))
    return peps, nucs


CASES = {}


def case(name, note=""):
    def deco(fn):
        CASES[name] = (fn, note)
        return fn
    return deco


# ---------------------------------------------------------------- the cases: each returns (args, stdin bytes|None)
@case("empty_files", "both files empty")
def _(d, rng):
    write(d, "pep.fasta", ""); write(d, "nuc.fasta", "")
    return ["pep.fasta", "nuc.fasta"], None


@case("empty_dna", "protein ok, DNA file empty")
def _(d, rng):
    p, n = genes(rng, 3, 20); write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", "")
    return ["pep.fasta", "nuc.fasta"], None


@case("header_only", "headers without sequences")
def _(d, rng):
    write(d, "pep.fasta", ">a\n>b\n"); write(d, "nuc.fasta", ">a\n>b\n")
    return ["pep.fasta", "nuc.fasta"], None


@case("one_sequence", "a single gene")
def _(d, rng):
    p, n = genes(rng, 1, 30); write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    return ["pep.fasta", "nuc.fasta"], None


@case("very_long_gene", "one gene of 1,000,000 codons (PAL2NAL's regex limit is about 22,000 residues)")
def _(d, rng):
    p, n = gene(rng, 1_000_000); write(d, "pep.fasta", fasta([("long", p)], 60)); write(d, "nuc.fasta", fasta([("long", n)], 60))
    return ["pep.fasta", "nuc.fasta"], None


@case("protein_all_gaps", "one row is only gaps")
def _(d, rng):
    p, n = genes(rng, 3, 20)
    p[1] = (p[1][0], "-" * 20)
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    return ["pep.fasta", "nuc.fasta"], None


@case("all_gap_columns", "two columns are gaps in every row")
def _(d, rng):
    p, n = genes(rng, 3, 20)
    p = [(i, r[:5] + "--" + r[5:]) for i, r in p]
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    return ["pep.fasta", "nuc.fasta"], None


@case("internal_stop_star", "an internal stop codon with '*' in the protein")
def _(d, rng):
    p, n = genes(rng, 3, 20)
    i, r = p[0]; p[0] = (i, r[:10] + "*" + r[11:])
    j, s = n[0]; n[0] = (j, s[:30] + "TAA" + s[33:])
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    return ["pep.fasta", "nuc.fasta"], None


@case("internal_stop_no_star", "an internal stop codon where the protein has a residue")
def _(d, rng):
    p, n = genes(rng, 3, 20)
    j, s = n[0]; n[0] = (j, s[:30] + "TAA" + s[33:])
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    return ["pep.fasta", "nuc.fasta"], None


@case("terminal_stop", "DNA ends with a stop codon the protein does not have")
def _(d, rng):
    p, n = genes(rng, 3, 20)
    n = [(i, s + "TGA") for i, s in n]
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    return ["pep.fasta", "nuc.fasta"], None


@case("terminal_stop_star", "protein ends with '*' and the DNA with a stop")
def _(d, rng):
    p, n = genes(rng, 3, 20)
    p = [(i, r + "*") for i, r in p]; n = [(i, s + "TAG") for i, s in n]
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    return ["pep.fasta", "nuc.fasta"], None


@case("lowercase_mixed", "lower-case protein and mixed-case DNA")
def _(d, rng):
    p, n = genes(rng, 3, 20)
    p = [(i, r.lower()) for i, r in p]
    n = [(i, "".join(c.lower() if k % 2 else c for k, c in enumerate(s))) for i, s in n]
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    return ["pep.fasta", "nuc.fasta"], None


@case("iupac_and_n_runs", "IUPAC codes and a run of N in the DNA")
def _(d, rng):
    p, n = genes(rng, 3, 30)
    j, s = n[0]; n[0] = (j, s[:9] + "NNNNNNNNN" + s[18:30] + "GCN" + s[33:])
    j, s = n[1]; n[1] = (j, s[:6] + "CTY" + s[9:])  # Y = C/T: CTC/CTT both L
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    return ["pep.fasta", "nuc.fasta"], None


@case("xuobzj_residues", "X, U, O, B, Z and J residues")
def _(d, rng):
    p, n = genes(rng, 2, 30)
    i, r = p[0]; p[0] = (i, r[:3] + "X" + r[4:8] + "U" + r[9:12] + "O" + r[13:16] + "B" + r[17:20] + "Z" + r[21:24] + "J" + r[25:])
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    return ["pep.fasta", "nuc.fasta"], None


@case("frameshift_digits", "frameshift digits 1, 2, 4 in the protein alignment")
def _(d, rng):
    p, n = genes(rng, 3, 20)
    i, r = p[0]; p[0] = (i, r[:5] + "2" + r[6:])
    j, s = n[0]; n[0] = (j, s[:15] + "GT" + s[18:])
    i, r = p[1]; p[1] = (i, r[:5] + "1" + r[6:])
    j, s = n[1]; n[1] = (j, s[:15] + "G" + s[18:])
    i, r = p[2]; p[2] = (i, r[:5] + "4" + r[6:])
    j, s = n[2]; n[2] = (j, s[:15] + "GTAC" + s[18:])
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    return ["pep.fasta", "nuc.fasta"], None


@case("tfastx_notation", "tfastx '\\' (1 base) and '/' (gap then 2 bases)")
def _(d, rng):
    p, n = genes(rng, 2, 20)
    i, r = p[0]; p[0] = (i, r[:5] + "\\" + r[6:])
    j, s = n[0]; n[0] = (j, s[:15] + "G" + s[18:])
    i, r = p[1]; p[1] = (i, r[:5] + "/" + r[6:])
    j, s = n[1]; n[1] = (j, s[:15] + "GT" + s[18:])
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    return ["pep.fasta", "nuc.fasta"], None


@case("utr_before_after", "extra DNA before and after every gene")
def _(d, rng):
    p, n = genes(rng, 3, 30)
    n = [(i, "".join(rng.choice("ACGT") for _ in range(50)) + s + "TAA" + "A" * 30) for i, s in n]
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    return ["pep.fasta", "nuc.fasta"], None


@case("reverse_complement", "the DNA is the reverse complement of the gene: must fail clearly, never align silently")
def _(d, rng):
    p, n = genes(rng, 3, 30)
    j, s = n[1]; n[1] = (j, "".join(RC[c] for c in reversed(s)))
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    return ["pep.fasta", "nuc.fasta"], None


@case("wrong_gene_dna", "one protein gets another gene's DNA")
def _(d, rng):
    p, n = genes(rng, 4, 40)
    other, _ = gene(rng, 40)
    _, s = gene(rng, 40)
    n[2] = (n[2][0], s)
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    return ["pep.fasta", "nuc.fasta"], None


@case("wrong_genetic_code", "mitochondrial (code 2) genes run with the standard code")
def _(d, rng):
    aas = dict(CODE1, W="TGA", M="ATA")
    peps, nucs = [], []
    for k in range(6):
        pp, ss = gene(rng, 40, aas)
        peps.append((f"g{k}", pp)); nucs.append((f"g{k}", ss))
    write(d, "pep.fasta", fasta(peps)); write(d, "nuc.fasta", fasta(nucs))
    return ["pep.fasta", "nuc.fasta"], None


@case("right_genetic_code_2", "the same genes with -codontable 2")
def _(d, rng):
    aas = dict(CODE1, W="TGA", M="ATA")
    peps, nucs = [], []
    for k in range(6):
        pp, ss = gene(rng, 40, aas)
        peps.append((f"g{k}", pp)); nucs.append((f"g{k}", ss))
    write(d, "pep.fasta", fasta(peps)); write(d, "nuc.fasta", fasta(nucs))
    return ["pep.fasta", "nuc.fasta", "-codontable", "2"], None


@case("duplicate_ids", "two proteins with the same ID")
def _(d, rng):
    p, n = genes(rng, 3, 20); p[2] = ("g0", p[2][1]); n[2] = ("g0", n[2][1])
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    return ["pep.fasta", "nuc.fasta"], None


@case("empty_id", "a '>' header without a name")
def _(d, rng):
    p, n = genes(rng, 3, 20); p[1] = ("", p[1][1]); n[1] = ("", n[1][1])
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    return ["pep.fasta", "nuc.fasta"], None


@case("pipe_style_ids", "sp|P12345|NAME_HUMAN style IDs in both files")
def _(d, rng):
    p, n = genes(rng, 3, 20)
    p = [(f"sp|P1234{k}|G{k}_HUMAN", r) for k, (i, r) in enumerate(p)]
    n = [(f"sp|P1234{k}|G{k}_HUMAN", s) for k, (i, s) in enumerate(n)]
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    return ["pep.fasta", "nuc.fasta"], None


@case("versioned_ids_reordered", "protein IDs with .1 suffix, DNA without, DNA in another order")
def _(d, rng):
    p, n = genes(rng, 4, 20)
    p = [(f"AT1G0{k}.1", r) for k, (i, r) in enumerate(p)]
    n = [(f"AT1G0{k}", s) for k, (i, s) in enumerate(n)]
    n = n[::-1]
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    return ["pep.fasta", "nuc.fasta"], None


@case("whitespace_ids", "IDs with descriptions after a space and a tab")
def _(d, rng):
    p, n = genes(rng, 3, 20)
    p = [(f"{i} some description here", r) for i, r in p]
    n = [(f"{i}\tother description", s) for i, s in n]
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    return ["pep.fasta", "nuc.fasta"], None


@case("different_counts", "DNA has one record more")
def _(d, rng):
    p, n = genes(rng, 3, 20); extra = gene(rng, 20); n.append(("extra", extra[1]))
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    return ["pep.fasta", "nuc.fasta"], None


@case("ids_other_order", "same IDs, DNA records in another order")
def _(d, rng):
    p, n = genes(rng, 5, 20); n = n[::-1]
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    return ["pep.fasta", "nuc.fasta"], None


@case("crlf", "CRLF line ends in both files")
def _(d, rng):
    p, n = genes(rng, 3, 20)
    write(d, "pep.fasta", fasta(p).replace("\n", "\r\n").encode()); write(d, "nuc.fasta", fasta(n).replace("\n", "\r\n").encode())
    return ["pep.fasta", "nuc.fasta"], None


@case("lone_cr", "lone CR line ends")
def _(d, rng):
    p, n = genes(rng, 3, 20)
    write(d, "pep.fasta", fasta(p).replace("\n", "\r").encode()); write(d, "nuc.fasta", fasta(n).replace("\n", "\r").encode())
    return ["pep.fasta", "nuc.fasta"], None


@case("utf8_bom", "UTF-8 BOM at the start of both files")
def _(d, rng):
    p, n = genes(rng, 3, 20)
    write(d, "pep.fasta", b"\xef\xbb\xbf" + fasta(p).encode()); write(d, "nuc.fasta", b"\xef\xbb\xbf" + fasta(n).encode())
    return ["pep.fasta", "nuc.fasta"], None


@case("non_ascii_headers", "non-ASCII characters in the headers")
def _(d, rng):
    p, n = genes(rng, 3, 20)
    p = [(f"gène_{k}_ü", r) for k, (i, r) in enumerate(p)]; n = [(f"gène_{k}_ü", s) for k, (i, s) in enumerate(n)]
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    return ["pep.fasta", "nuc.fasta"], None


@case("clustal_input", "CLUSTAL protein alignment")
def _(d, rng):
    p, n = genes(rng, 4, 70)
    write(d, "pep.aln", clustal(p)); write(d, "nuc.fasta", fasta(n))
    return ["pep.aln", "nuc.fasta"], None


@case("clustal_marks_blockonly", "CLUSTAL with '#' marks and -blockonly")
def _(d, rng):
    p, n = genes(rng, 3, 70)
    marks = "".join("#" if 10 <= c < 50 else " " for c in range(70))
    write(d, "pep.aln", clustal(p, marks)); write(d, "nuc.fasta", fasta(n))
    return ["pep.aln", "nuc.fasta", "-blockonly"], None


@case("clustal_no_header_line", "CLUSTAL without the 'CLUSTAL W' first line (variant)")
def _(d, rng):
    p, n = genes(rng, 3, 30)
    write(d, "pep.aln", clustal(p).split("\n", 1)[1]); write(d, "nuc.fasta", fasta(n))
    return ["pep.aln", "nuc.fasta"], None


@case("path_with_spaces_and_unicode", "input paths with spaces and non-ASCII characters")
def _(d, rng):
    p, n = genes(rng, 3, 20)
    sub = d / "my data é"
    write(sub, "pep file.fasta", fasta(p)); write(sub, "nuc file.fasta", fasta(n))
    return ["my data é/pep file.fasta", "my data é/nuc file.fasta"], None


@case("readonly_output_folder", "-o into a read-only folder")
def _(d, rng):
    p, n = genes(rng, 3, 20)
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    ro = d / "ro"; ro.mkdir(exist_ok=True); os.chmod(ro, stat.S_IRUSR | stat.S_IXUSR)
    return ["pep.fasta", "nuc.fasta", "-o", "ro/out.fa"], None


@case("output_equals_input", "-o names the DNA input: must never destroy it")
def _(d, rng):
    p, n = genes(rng, 3, 20)
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    return ["pep.fasta", "nuc.fasta", "-o", "nuc.fasta"], None


@case("stdin_protein", "protein alignment from standard input")
def _(d, rng):
    p, n = genes(rng, 3, 20)
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    return ["-", "nuc.fasta"], (d / "pep.fasta").read_bytes()


@case("stdin_dna", "DNA from standard input")
def _(d, rng):
    p, n = genes(rng, 3, 20)
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    return ["pep.fasta", "-"], (d / "nuc.fasta").read_bytes()


@case("two_dna_files", "DNA split over two files")
def _(d, rng):
    p, n = genes(rng, 4, 20)
    write(d, "pep.fasta", fasta(p)); write(d, "nuc0.fasta", fasta(n[:2])); write(d, "nuc1.fasta", fasta(n[2:]))
    return ["pep.fasta", "nuc0.fasta", "nuc1.fasta"], None


@case("nogap_code2", "-nogap with code 2 (TGA = Trp kept, AGA/AGG stops removed)")
def _(d, rng):
    aas = dict(CODE1, W="TGA")
    peps, nucs = [], []
    for k in range(3):
        pp, ss = gene(rng, 30, aas)
        pp = pp[:10] + "*" + pp[11:]; ss = ss[:30] + "AGA" + ss[33:]
        peps.append((f"g{k}", pp)); nucs.append((f"g{k}", ss))
    write(d, "pep.fasta", fasta(peps)); write(d, "nuc.fasta", fasta(nucs))
    return ["pep.fasta", "nuc.fasta", "-codontable", "2", "-nogap"], None


@case("nomismatch_mixed", "-nomismatch on genes with a few wrong codons")
def _(d, rng):
    p, n = genes(rng, 4, 30)
    j, s = n[1]; n[1] = (j, s[:15] + CODE1["W" if p[1][1][5] != "W" else "A"] + s[18:])
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    return ["pep.fasta", "nuc.fasta", "-nomismatch"], None


@case("codon_format", "-output codon")
def _(d, rng):
    p, n = genes(rng, 3, 25)
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    return ["pep.fasta", "nuc.fasta", "-output", "codon"], None


@case("paml_format", "-output paml")
def _(d, rng):
    p, n = genes(rng, 3, 25)
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    return ["pep.fasta", "nuc.fasta", "-output", "paml"], None


@case("clustal_format", "-output clustal")
def _(d, rng):
    p, n = genes(rng, 3, 65)
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    return ["pep.fasta", "nuc.fasta", "-output", "clustal"], None


@case("dna_too_short", "DNA truncated before the end of the protein")
def _(d, rng):
    p, n = genes(rng, 3, 30)
    j, s = n[0]; n[0] = (j, s[:60])
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    return ["pep.fasta", "nuc.fasta"], None


@case("files_swapped", "DNA given as the protein alignment and vice versa")
def _(d, rng):
    p, n = genes(rng, 3, 30)
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    return ["nuc.fasta", "pep.fasta"], None


@case("binary_file", "a binary (NUL bytes) file as the DNA")
def _(d, rng):
    p, n = genes(rng, 3, 30)
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", b">g0\n" + bytes(range(256)))
    return ["pep.fasta", "nuc.fasta"], None


@case("rows_different_lengths", "protein rows of different lengths (not aligned)")
def _(d, rng):
    p, n = genes(rng, 3, 30)
    i, r = p[1]; p[1] = (i, r + "AA")
    write(d, "pep.fasta", fasta(p)); write(d, "nuc.fasta", fasta(n))
    return ["pep.fasta", "nuc.fasta"], None


# ---------------------------------------------------------------- running
def run(cmd, cwd, stdin, env=None, timeout=900):
    try:
        p = subprocess.run(cmd, cwd=cwd, input=stdin, capture_output=True, timeout=timeout, env=env)
        return {"rc": p.returncode, "stdout": p.stdout, "stderr": p.stderr, "timeout": False}
    except subprocess.TimeoutExpired as e:
        return {"rc": None, "stdout": e.stdout or b"", "stderr": e.stderr or b"", "timeout": True}


def digest(b: bytes):
    return {"sha256": hashlib.sha256(b).hexdigest(), "bytes": len(b), "head": b[:160].decode("utf-8", "replace")}


def pal2nal_args(args):
    """The PAL2NAL subset of a coal command line (-o and COAL-only options dropped)."""
    out, it = [], iter(args)
    for a in it:
        if a == "-o":
            next(it, None)
        elif a in ("-codontable", "-output"):
            out += [a, next(it)]
        elif a in ("-nogap", "-nomismatch", "-blockonly"):
            out.append(a)
        elif a.startswith("-") and a != "-":
            continue
        else:
            out.append(a)
    return out


def main(outdir):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    results = {}
    for name, (fn, note) in CASES.items():
        d = out / name
        if d.exists():
            for root, dirs, files in os.walk(d):
                for x in dirs:
                    os.chmod(Path(root) / x, 0o755)
            shutil.rmtree(d)
        d.mkdir(parents=True)
        rng = random.Random(name)
        args, stdin = fn(d, rng)
        fmt_args = args if "-output" in args else args + ["-output", "fasta"]
        rec = {"note": note, "args": fmt_args, "tools": {}}
        for tool, cmd, env in (("python", PY + fmt_args, ENV), ("cpp", [CPP] + fmt_args, ENV),
                               ("pal2nal", ["perl", PAL2NAL] + pal2nal_args(fmt_args), None)):
            if tool == "pal2nal" and "-" in fmt_args:
                r = {"rc": None, "stdout": b"", "stderr": b"(not run: PAL2NAL has no standard-input mode)", "timeout": False}
            elif tool == "pal2nal" and name == "very_long_gene":
                r = run(cmd, d, stdin, env, timeout=600)
            else:
                r = run(cmd, d, stdin, env)
            rec["tools"][tool] = {"rc": r["rc"], "timeout": r["timeout"], "stdout": digest(r["stdout"]),
                                  "stderr_lines": r["stderr"].decode("utf-8", "replace").splitlines()[:6],
                                  "stderr_sha": hashlib.sha256(r["stderr"]).hexdigest()}
            (d / f"stdout.{tool}").write_bytes(r["stdout"])
            (d / f"stderr.{tool}").write_bytes(r["stderr"])
            # a produced -o file: keep its digest and check the inputs survived
            if "-o" in fmt_args:
                o = d / fmt_args[fmt_args.index("-o") + 1]
                rec["tools"][tool]["outfile"] = digest(o.read_bytes()) if o.exists() else None
        py, cpp, pal = rec["tools"]["python"], rec["tools"]["cpp"], rec["tools"]["pal2nal"]
        rec["engines_identical"] = (py["rc"] == cpp["rc"] and py["stdout"]["sha256"] == cpp["stdout"]["sha256"]
                                    and py["stderr_sha"] == cpp["stderr_sha"])
        pal_ok = pal["rc"] == 0 and pal["stdout"]["bytes"] > 0
        if pal["rc"] is None:
            rec["vs_pal2nal"] = "not run"
        elif pal_ok and py["rc"] == 0 and py["stdout"]["sha256"] == pal["stdout"]["sha256"]:
            rec["vs_pal2nal"] = "identical"
        elif not pal_ok and py["rc"] != 0:
            rec["vs_pal2nal"] = "both fail"
        elif pal_ok and py["rc"] != 0:
            rec["vs_pal2nal"] = "PAL2NAL output, COAL error"
        elif not pal_ok and py["rc"] == 0:
            rec["vs_pal2nal"] = "COAL output, PAL2NAL failed"
        else:
            rec["vs_pal2nal"] = "DIFFERENT output"
        # the leftovers: anything but the inputs and our stdout/stderr files
        rec["leftovers"] = sorted(x.name for x in d.iterdir() if x.name.startswith((".coal", "out.")) and
                                  not x.name.startswith(("stdout.", "stderr.")) and x.name not in fmt_args)
        results[name] = rec
        print(f"{name:32s} py rc={py['rc']} cpp rc={cpp['rc']} pal rc={pal['rc']} engines={'=' if rec['engines_identical'] else 'DIFF'} "
              f"vs PAL2NAL: {rec['vs_pal2nal']}" + (f" leftovers={rec['leftovers']}" if rec["leftovers"] else ""), flush=True)
    (out / "results.json").write_text(json.dumps(results, indent=1), encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "adversarial_out"))
