"""Differential fuzz: COAL's ``pal2nal`` drop-in vs the real PAL2NAL v14, on random inputs.

Needs Perl and PAL2NAL v14 (fetch it with ``python scripts/fetch_pal2nal.py``; it is GPL and never committed):

    PAL2NAL=PAL2NAL/pal2nal.v14/pal2nal.pl python tests/fuzz_vs_pal2nal.py          # FUZZ_N cases (default 3000)
    PAL2NAL=... FUZZ_N=500 FUZZ_SEED=7 python tests/fuzz_vs_pal2nal.py
    PAL2NAL=... COAL_FUZZ_ENGINE=cpp COAL_CPP_BIN=build/coal-cpp python tests/fuzz_vs_pal2nal.py   # coal-cpp instead

Without PAL2NAL set the test is skipped. By default COAL is the Python ``pal2nal`` command run in-process;
``COAL_FUZZ_ENGINE=cpp`` runs coal-cpp (``COAL_CPP_BIN``, default ``build/coal-cpp``) through the same comparison,
with the PAL2NAL command line translated by the drop-in's own ``coal_compat.cli.translate`` (what
tests/test_cpp_vs_reference.py does). ``FUZZ_KEEP=DIR`` copies every case with an unexplained difference to
``DIR/c<N>/`` with its command line in ``args.txt``. Every case is a random protein alignment (FASTA or CLUSTAL, sometimes with
a '#'-marked block line) plus coding DNA (UTRs, stop codons, mismatching / unrelated / masked / lower-case / RNA DNA,
frameshift digits, other genetic codes, shuffled / different / partly matching IDs ...) and a random option set
(-codontable, -nogap, -nomismatch, -blockonly). Each case runs through both tools with exactly the same command
line, once for each output format (-output clustal / paml / fasta / codon), so 4 comparisons per case.

Rule per comparison: PAL2NAL always exits 0; it succeeded if its stdout is non-empty. Then COAL (the ``pal2nal``
command, run in-process) must exit 0 with byte-identical stdout; if PAL2NAL failed, COAL must exit non-zero with
empty stdout. stderr is not compared (COAL's messages are its own). The only exceptions are the deliberate
differences in DELIBERATE below (each is in docs/differences.md).
A difference is accepted only if it is exactly the documented one. The generator
also writes a neutralised copy of each case (U residues -> B, RNA -> DNA, alternative start codons -> ATG, table-10
TGA codons -> TGC, each IUPAC codon -> one of its readings chosen so that it matches the residue / is a stop exactly
when COAL says so; substituted codons are written in the other letter case so they can be recognised), and the real
PAL2NAL is run on that copy without -nogap. COAL's output must then equal that reference output, with -nogap applied
by COAL's documented rule (gap or stop codon of the selected code, any case) and only the substituted codons allowed
to differ (and only by the documented substitution); an invalid protein character must be a COAL error. Anything
else is an unexplained difference.
"""
from __future__ import annotations

import io
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from coal.codon import ANCHOR_STARTS, IUPAC, accepted  # noqa: E402
from coal.core import nogap_bad_triplets  # noqa: E402
from coal.tables import SUPPORTED, get_code  # noqa: E402
from coal_compat.cli import main as pal2nal_main  # noqa: E402

FORMATS = ("clustal", "paml", "fasta", "codon")

# Deliberate differences (docs/differences.md): feature -> why COAL may differ from PAL2NAL.
# A case is allowed to differ only if the generator gave it one of these features.
DELIBERATE = {
    "invalid_char": "invalid protein characters (e.g. '@', '#'): PAL2NAL drops a codon silently, COAL refuses (exit 1)",
    "genetic_code": "COAL uses the current NCBI tables, PAL2NAL v14 older ones. Measured against the real PAL2NAL "
                    "(every codon x residue, every first-Met codon): table 10 TGA (NCBI Cys, PAL2NAL stop) and the "
                    "first-Met start codons GTG in table 3, ATT/ATC/ATA in table 11 and ATA in table 13 (NCBI "
                    "starts PAL2NAL lacks). Nothing else differs",
    "iupac": "IUPAC codes other than N in the DNA: COAL accepts a codon if every reading gives the residue (and "
             "-nogap removes it if every reading is a stop, e.g. TRA); PAL2NAL has a fixed list",
    "u_residue": "U (selenocysteine) in the protein: COAL treats it as X (any codon, with a warning) like B/Z/J/O; "
                 "PAL2NAL's search pattern gives U no codon at all, so the row is misaligned or the gene fails",
    "nogap_code": "-nogap removes the stop codons of the selected genetic code, in any case (COAL, deliberate); "
                  "PAL2NAL always removes the standard code's TAA/TAG/TGA, upper case only",
    "first_met_code": "the first Met of a gene that needs the fallback (one N, one mismatch, ...): COAL accepts the "
                      "selected code's start and Met codons (e.g. ATT and ATA in code 2, TTG in code 11), as in its "
                      "exact-fit path; PAL2NAL accepts only ATG/CTG/GTG there, whatever -codontable says (FIN-E1)",
    "normalised_ids": "DNA IDs equal to the protein IDs only after normalising (version suffix, case, database "
                      "fields such as sp|...|) and in another order: PAL2NAL pairs by file order (each protein gets "
                      "another gene's DNA); COAL pairs by the normalised IDs. Exact check: COAL's output is "
                      "byte-identical to PAL2NAL's on the same case with the DNA records put in protein order",
    "rna_first_met": "RNA (U) input with a protein starting with M: when PAL2NAL needs its fallback, it does not "
                     "accept AUG/CUG/GUG as the first Met (its fallback pattern has T only), so it flags that "
                     "codon (dropped by -nomismatch) or places the gene elsewhere; COAL reads U as T everywhere",
}
# Documented deliberate differences the generator does not produce (README.md lists them with repros): protein rows
# of different lengths (COAL error); FASTA starting with a blank line, BOM, "> id" (COAL reads them, PAL2NAL fails);
# duplicate IDs, empty DNA header, headerless second DNA file, '_' as stop (COAL error); B/Z/J/O fallback warnings
# (stderr only, not compared); partial-N codons such as NTG (FZ-D8).
# Codons behind "genetic_code" (run table -> codons): TGA anywhere in the DNA; a first Met on these codons.
TABLE_CODONS = {10: {"TGA"}}
START_CODONS = {3: {"GTG"}, 11: {"ATT", "ATC", "ATA"}, 13: {"ATA"}}

AA = "ACDEFGHIKLMNPQRSTVWY"


def _codons(table):
    code = get_code(table)
    by = {a: sorted(c) for a, c in code.codons.items()}
    return code, by


class Case:
    def __init__(self):
        self.files = {}      # name -> text
        self.ref_files = {}  # name -> neutralised text (see the module docstring)
        self.ref_files_b = {}  # the same without the first-Met substitution (for a gene placed off its own frame)
        self.ref_files_p = {}  # the case itself with the DNA records in protein order (normalised_ids)
        self.pep = None
        self.nucs = []
        self.opts = []
        self.features = set()
        self.info = set()  # what the case exercises (reported as coverage counts)


def _readings(trip):
    return ["".join(x) for x in __import__("itertools").product(*(IUPAC.get(b, b) for b in trip.upper()))]


def _reference_reading(table, aa, trip):
    """A concrete reading of an IUPAC codon that matches residue `aa` exactly when COAL accepts
    the codon (every reading matches) and is a stop exactly when COAL's -nogap removes it (every reading is a stop);
    None if no reading has both properties."""
    code = get_code(table).code
    reads = _readings(trip)

    def enc(r):
        return True if aa in "XBZJOU" else r in accepted(table, aa)
    all_enc = trip.upper() in accepted(table, aa) if aa not in "XBZJOU" else True
    all_stop = all(code[r] == "*" for r in reads)
    for r in reads:
        if enc(r) == all_enc and (code[r] == "*") == all_stop:
            return r
    return None


def _flip(s, marks):
    """The marked positions in the other letter case."""
    t = list(s)
    for i in marks:
        if 0 <= i < len(t):
            t[i] = t[i].lower() if t[i].isupper() else t[i].upper()
    return "".join(t)


def _wrap(seq, width):
    if not width or not seq:
        return seq + "\n"
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width)) + "\n"


def make_case(rng: random.Random) -> Case:
    c = Case()
    gen_table = 1 if rng.random() < 0.6 else rng.choice(SUPPORTED)
    run_table = gen_table if rng.random() < 0.85 else rng.choice(SUPPORTED)
    code, by = _codons(gen_table)
    nseq = rng.choice((1, 2, 2, 3, 3, 4, 5, 6))
    r = rng.random()
    length = rng.randint(1, 12) if r < 0.15 else rng.randint(13, 70) if r < 0.85 else rng.randint(71, 260)
    anc = [rng.choice(AA) for _ in range(length)]
    if rng.random() < 0.7:
        anc[0] = "M"
    ncol_extra = rng.randint(0, 4)
    cols = list(range(length))
    ins_at = sorted(rng.randint(0, length) for _ in range(ncol_extra))
    p_gap = rng.choice((0, 0.03, 0.08, 0.2))
    p_sub = rng.choice((0, 0.05, 0.15))
    rows = []
    for k in range(nseq):
        row = []
        ins = list(ins_at)
        for i in cols:
            while ins and ins[0] == i:
                ins.pop(0)
                row.append(rng.choice(AA) if rng.random() < 0.4 else "-")
            if rng.random() < p_gap and i > 0:
                row.append("-")
            else:
                row.append(rng.choice(AA) if rng.random() < p_sub else anc[i])
        while ins:
            ins.pop(0)
            row.append(rng.choice(AA) if rng.random() < 0.4 else "-")
        rows.append(row)
    # rare residue kinds
    if rng.random() < 0.06:
        k, i = rng.randrange(nseq), rng.randrange(len(rows[0]))
        if rows[k][i] != "-":
            rows[k][i] = rng.choice("BZJOUX")
            if rows[k][i] == "U":
                c.features.add("u_residue")
    if rng.random() < 0.05:
        k, i = rng.randrange(nseq), rng.randrange(len(rows[0]))
        if rows[k][i] != "-" and i > 0:
            rows[k][i] = "*"
    if rng.random() < 0.05:  # frameshift digit
        k, i = rng.randrange(nseq), rng.randrange(len(rows[0]))
        if rows[k][i] != "-" and i > 0:
            rows[k][i] = str(rng.choice((1, 2, 4, 5)))
    if rng.random() < 0.04 and len(rows[0]) > 1:  # tfastx/tfasty frameshift notation (needs 2+ columns): '\\' = 1 base, '/' + residue = gap + 2 bases
        k, i = rng.randrange(nseq), rng.randrange(1, len(rows[0]))
        j = next((j for j in range(i + 1, len(rows[0])) if rows[k][j] != "-"), None)
        if rows[k][i].isalpha():
            if rng.random() < 0.5:
                rows[k][i] = "\\"
            elif j is not None and rows[k][j].isalpha():
                rows[k][i] = "/"
    if rng.random() < 0.01:
        k, i = rng.randrange(nseq), rng.randrange(len(rows[0]))
        rows[k][i] = rng.choice("@#")
        c.features.add("invalid_char")
    term_stop = rng.random() < 0.15
    if term_stop:
        for row in rows:
            row.append("*")

    # IDs
    def name(k):
        if rng.random() < 0.8:
            return f"s{k}"
        return "".join(rng.choice("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.|")
                       for _ in range(rng.randint(1, 24)))
    pids = []
    while len(pids) < nseq:
        n = name(len(pids))
        if n not in pids:
            pids.append(n)

    # DNA
    dnas = []
    mode = rng.choices(("clean", "mismatch", "unrelated", "masked", "wrongtable"),
                       (0.55, 0.22, 0.06, 0.12, 0.05))[0]
    ref_dnas = []  # neutralised DNA (see the module docstring)
    for k, row in enumerate(rows):
        out = []
        out_res = []  # residue of each entry of `out` (None for frameshift bases)
        first = True
        start_idx = None  # index in `out` of the codon of the first Met (an M before any non-X residue)
        slash = False
        for ch in row:
            if ch == "-":
                continue
            if ch == "/":
                slash = True
                continue
            if ch.isdigit() or ch == "\\" or slash:
                n = int(ch) if ch.isdigit() else 1 if ch == "\\" else 2
                slash = False
                out.append("".join(rng.choice("ACGT") for _ in range(n)))
                out_res.append(None)
                continue
            if ch == "M" and first and rng.random() < 0.3 and code.starts:
                cod = rng.choice(sorted(code.starts))
            elif ch in by:
                cod = rng.choice(by[ch])
            else:
                cod = "".join(rng.choice("ACGT") for _ in range(3))
            if first and ch == "M":
                start_idx = len(out)
            out.append(cod)
            out_res.append(ch)
            if ch not in "XBZJOU":
                first = False
        cds = out
        if mode == "mismatch":
            frac = rng.choice((0.02, 0.1, 0.3, 0.5, 0.6, 0.9))
            cds = [("".join(rng.choice("ACGT") for _ in range(3)) if len(x) == 3 and rng.random() < frac else x)
                   for x in cds]
        elif mode == "unrelated":
            cds = ["".join(rng.choice("ACGT") for _ in range(len(x))) for x in cds]
        elif mode == "masked":
            frac = rng.choice((0.3, 0.6, 0.85, 0.95, 1.0))
            cds = [("NNN" if len(x) == 3 and rng.random() < frac else x) for x in cds]
        cds = list(cds)
        ref = list(cds)
        marked = set()  # indices into cds/ref whose reference codon is a substitution
        # IUPAC codes in the gene's codons only, each with a
        # concrete reading for the reference run (codons with no reading of the right kind are left unchanged)
        if rng.random() < 0.04:
            for i, (x, aa) in enumerate(zip(cds, out_res)):
                if aa is None or i == start_idx or len(x) != 3 or x.strip("ACGT"):
                    continue
                trip = "".join(rng.choice("RYSWKMBDHV") if rng.random() < 0.05 else b for b in x)
                if trip == x:
                    continue
                r = _reference_reading(run_table, aa, trip)
                if r is None:
                    continue
                cds[i], ref[i] = trip, r
                marked.add(i)
                c.features.add("iupac")
        keep_tga = set()  # table 10: TGA is Cys (NCBI), a stop for PAL2NAL v14 (neutralised as TGC for Cys codons)
        if run_table in TABLE_CODONS:
            pos = 0
            for i, x in enumerate(ref):
                if x.upper() == "TGA" and out_res[i] == "C":
                    ref[i] = "TGC"
                    marked.add(i)
                    keep_tga.add(pos)
                pos += len(x)
        first_met = cds[start_idx].upper() if start_idx is not None else None
        if first_met in START_CODONS.get(run_table, ()):
            c.features.add("genetic_code")
        elif first_met in FIRST_MET(run_table) and first_met not in ANCHOR_STARTS:
            c.features.add("first_met_code")  # matters only when the gene needs the fallback
        start_marked = set()
        if first_met in START_CODONS.get(run_table, ()) or (
                first_met in FIRST_MET(run_table) and first_met not in ANCHOR_STARTS):
            ref[start_idx] = "ATG"
            start_marked.add(start_idx)
        seq = "".join(cds)
        rseq = "".join(ref)
        rseq_b = "".join(cds[i] if i in start_marked else x for i, x in enumerate(ref))  # without the start fix
        pos, mpos, mpos_b = 0, set(), set()
        for i, x in enumerate(ref):
            if i in marked or i in start_marked:
                mpos.update(range(pos, pos + len(x)))
            if i in marked:
                mpos_b.update(range(pos, pos + len(x)))
            pos += len(x)
        if rng.random() < 0.3:
            utr = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 30)))
            seq, rseq, rseq_b = utr + seq, utr + rseq, utr + rseq_b
            mpos, mpos_b = {p + len(utr) for p in mpos}, {p + len(utr) for p in mpos_b}
            keep_tga = {p + len(utr) for p in keep_tga}
        if rng.random() < 0.3:
            tail = rng.choice(("TAA", "TAG", "TGA")) + "".join(rng.choice("ACGT") for _ in range(rng.randint(0, 20)))
            seq, rseq, rseq_b = seq + tail, rseq + tail, rseq_b + tail
        if run_table in TABLE_CODONS:
            # any other TGA (UTR, out of frame, a random codon) would change PAL2NAL's placement search in a
            # way the reference cannot neutralise, so it becomes TGG in both copies
            sl, rl, rb = list(seq), list(rseq), list(rseq_b)
            while True:
                hits = [m for m in range(len(sl) - 2) if "".join(sl[m:m + 3]) == "TGA" and m not in keep_tga]
                if not hits:
                    break
                for m in hits:
                    sl[m + 2] = rl[m + 2] = rb[m + 2] = "G"
            seq, rseq, rseq_b = "".join(sl), "".join(rl), "".join(rb)
        if rng.random() < 0.03 and len(seq) > 3:
            cut = rng.randint(1, 3)
            seq, rseq, rseq_b = seq[:-cut], rseq[:-cut], rseq_b[:-cut]
        if rng.random() < 0.05:
            seq = seq.replace("T", "U")  # the reference keeps T
            if start_idx is not None:
                c.features.add("rna_first_met")
        if rng.random() < 0.1:
            seq, rseq, rseq_b = seq.lower(), rseq.lower(), rseq_b.lower()
        dnas.append(seq)
        ref_dnas.append((_flip(rseq, mpos), _flip(rseq_b, mpos_b)))

    if any(t in d.upper().replace("U", "T") for d in dnas for t in TABLE_CODONS.get(run_table, ())):
        c.features.add("genetic_code")

    # DNA IDs and order
    did_mode = rng.choices(("same", "shuffled", "different", "partial", "count", "normalised"),
                           (0.57, 0.17, 0.1, 0.05, 0.03, 0.08))[0]
    order = list(range(nseq))
    dids = list(pids)
    if did_mode == "normalised":
        def variant(p):
            r = rng.random()
            if r < 0.15:
                return p
            if r < 0.45:
                return f"{p}.{rng.randint(1, 9)}"
            if r < 0.7 and p.swapcase() != p:
                return p.swapcase()
            return f"{rng.choice(('sp', 'tr', 'gi'))}|{p}|{p.upper()}_HUMAN"
        dids = [variant(p) for p in pids]
        if len(set(dids)) < nseq:
            dids = list(pids)  # a duplicate DNA ID is a COAL error (FZ-D7): keep the case plain
        else:
            c.features.add("normalised_ids")
        if rng.random() < 0.75:
            rng.shuffle(order)
    elif did_mode == "shuffled":
        rng.shuffle(order)
    elif did_mode == "different":
        dids = [f"d{k}x" for k in range(nseq)]
    elif did_mode == "partial" and nseq >= 2:
        dids = [p if k == 0 else f"other{k}" for k, p in enumerate(pids)]
    recs = [(dids[j], dnas[j], ref_dnas[j][0], ref_dnas[j][1]) for j in order]
    if did_mode == "count":
        recs = recs[:-1] if (len(recs) > 1 and rng.random() < 0.5) else recs + [("extra",) + ("ATGAAA",) * 3]
    # with normalised IDs every reference has the DNA in protein order (what COAL pairs)
    ref_recs = [(dids[j], dnas[j], ref_dnas[j][0], ref_dnas[j][1]) for j in range(nseq)] \
        if "normalised_ids" in c.features else recs
    width = rng.choice((0, 60, 60, 70, 13))
    split = len(recs) > 1 and rng.random() < 0.1
    cut = rng.randint(1, len(recs) - 1) if split else None
    nuc_parts, ref_parts, ref_parts_b, ref_parts_p = [], [], [], []
    for col, parts, rr in ((1, nuc_parts, recs), (2, ref_parts, ref_recs), (3, ref_parts_b, ref_recs),
                           (1, ref_parts_p, ref_recs)):
        if split:
            parts += ["".join(f">{r[0]}\n" + _wrap(r[col], width) for r in rr[:cut]),
                      "".join(f">{r[0]}\n" + _wrap(r[col], width) for r in rr[cut:])]
        else:
            parts.append("".join(f">{r[0]} desc\n" + _wrap(r[col], width) for r in rr))

    # protein file
    lower = rng.random() < 0.1
    srows = ["".join(r) for r in rows]
    if lower:
        srows = [s.lower() for s in srows]
    # the reference protein has B where the case has U (PAL2NAL treats B as COAL treats U)
    rrows = [s.replace("U", "B").replace("u", "b") for s in srows]
    ncol = len(srows[0])
    mask = None
    if rng.random() < 0.35:
        idw = max(len(p) for p in pids) + rng.choice((1, 2, 6))
        blk = rng.choice((60, 50, 7))
        hashed = rng.random() < 0.6
        if hashed:
            mask = []
            while len(mask) < ncol:
                n = rng.randint(1, 15)
                mask += [rng.choice("# ")] * n
            mask = "".join(mask[:ncol])
        text = rtext = "CLUSTAL W (1.83) multiple sequence alignment\n\n\n"
        for b in range(0, ncol, blk):
            for p, s, r in zip(pids, srows, rrows):
                text += p.ljust(idw) + s[b:b + blk] + "\n"
                rtext += p.ljust(idw) + r[b:b + blk] + "\n"
            last = b + blk >= ncol
            if mask is not None:
                ml = " " * idw + mask[b:b + blk].rstrip() + "\n"
            else:
                ml = " " * idw + "".join(rng.choice(" *:.") for _ in s[b:b + blk]).rstrip() + "\n"
            text, rtext = text + ml, rtext + ml
            if not last or rng.random() < 0.8:
                text, rtext = text + "\n", rtext + "\n"
        c.pep = "pep.aln"
    else:
        text = rtext = ""
        for p, s, r in zip(pids, srows, rrows):
            w = rng.choice((0, 60, 10))
            text += f">{p}\n" + _wrap(s, w)
            rtext += f">{p}\n" + _wrap(r, w)
        c.pep = "pep.fasta"
    c.files[c.pep] = text
    c.ref_files[c.pep] = rtext
    c.ref_files_b[c.pep] = rtext
    c.ref_files_p[c.pep] = text
    for i, part in enumerate(nuc_parts):
        c.files[f"nuc{i}.fasta"] = part
        c.ref_files[f"nuc{i}.fasta"] = ref_parts[i]
        c.ref_files_b[f"nuc{i}.fasta"] = ref_parts_b[i]
        c.ref_files_p[f"nuc{i}.fasta"] = ref_parts_p[i]
        c.nucs.append(f"nuc{i}.fasta")

    opts = []
    if run_table != 1 or rng.random() < 0.05:
        opts += ["-codontable", str(run_table)]
    if rng.random() < 0.25:
        opts.append("-nogap")
    if rng.random() < 0.25:
        opts.append("-nomismatch")
    if rng.random() < (0.5 if mask else 0.05):
        opts.append("-blockonly")
    stops = {k for k, v in get_code(run_table).code.items() if v == "*"}
    if "-nogap" in opts and (stops != {"TAA", "TAG", "TGA"} or any(ch.islower() for d in dnas for ch in d)):
        c.features.add("nogap_code")
    c.opts = opts
    c.info.add(f"mode {mode}")
    c.info.add(f"ids {did_mode}")
    if run_table != 1:
        c.info.add("codontable != 1")
    if any(ch in "\\/" for r in rows for ch in r):
        c.info.add("tfastx notation")
    if any(ch.isdigit() for r in rows for ch in r):
        c.info.add("frameshift digit")
    if any(ch.isdigit() and int(ch) >= 4 for r in rows for ch in r):
        c.info.add("frameshift digit >= 4")
    if 3 * ncol > 60:
        c.info.add("more than one 60-nt block")
    if c.pep == "pep.aln":
        c.info.add("CLUSTAL input")
    if mask and "#" in mask:
        c.info.add("'#' marks + -blockonly" if "-blockonly" in opts else "'#' marks shown (no -blockonly)")
    for o in ("-nogap", "-nomismatch"):
        if o in opts:
            c.info.add(o)
    return c


def parse_output(fmt, text):
    """[(id, aligned sequence)] from PAL2NAL/COAL output in any of the four formats
    (amino-acid, mark and header lines skipped)."""
    lines = text.split("\n")
    recs = []
    if fmt == "fasta":
        for ln in lines:
            if ln.startswith(">"):
                recs.append([ln[1:].strip(), ""])
            elif recs:
                recs[-1][1] += ln.strip()
    elif fmt == "paml":
        n, width = (int(x) for x in lines[0].split())
        i = 1
        for _ in range(n):
            while not lines[i].strip():
                i += 1
            sid, seq = lines[i].strip(), ""
            i += 1
            while len(seq) < width:
                seq += lines[i].strip()
                i += 1
            recs.append([sid, seq])
    else:  # clustal, codon
        order = {}
        for ln in lines[1:] if fmt == "clustal" else lines:
            if not ln.strip() or ln[0] in " #" or ln.startswith("CLUSTAL"):
                continue
            sid, _, rest = ln.partition(" ")
            if sid not in order:
                order[sid] = len(recs)
                recs.append([sid, ""])
            recs[order[sid]][1] += rest.replace(" ", "")
    return [tuple(r) for r in recs]


def _allowed_codon(table, coal3, ref3):
    """True if a substituted reference codon `ref3` is the documented neutralisation of
    COAL's codon `coal3`: one of its IUPAC readings, TGC for table 10's TGA, or ATG for a first-Met codon that is a
    start or Met codon of the selected code."""
    c = coal3.upper().replace("U", "T")
    r = ref3.upper()
    if len(c) == len(r) and all(a == b or b in IUPAC.get(a, a) for a, b in zip(c, r)):
        return True  # base by base: equal, or the reference has one reading of COAL's IUPAC base
    reads = _readings(c)
    if table in TABLE_CODONS and r == "TGC" and "TGA" in reads:
        return True
    return r == "ATG" and any(x in FIRST_MET(table) for x in reads)


def FIRST_MET(table):
    code = get_code(table)
    return code.starts | code.codons.get("M", frozenset())


def _nogap_rule(fmt, text, table):
    """COAL's documented -nogap rule applied to an output without -nogap: [(id, row)] with every 3-nt column that
    holds a gap or a stop codon of `table` (any case, U as T, IUPAC triplets whose every reading is a stop) removed."""
    rows = parse_output(fmt, text)
    bad = set().union(*(nogap_bad_triplets(sq, table) for _, sq in rows)) if rows else set()
    return [(i, "".join(sq[k:k + 3] for k in range(0, len(sq), 3) if k // 3 not in bad)) for i, sq in rows]


def check_deliberate(case, fmt, rc, out, ref_out, coal_without_nogap=None):
    """None if COAL's result is exactly the documented deliberate difference for this case,
    else the reason it is not. `ref_out`: PAL2NAL's stdout on the neutralised case (without -nogap).
    With -nogap, `coal_without_nogap()` gives COAL's (rc, stdout) for the same case without -nogap: that must match
    the reference, and COAL's -nogap output must be the documented -nogap rule applied to it (so an IUPAC codon's
    stop status is judged on COAL's own letters, e.g. TRA, and not on the reference's reading)."""
    table = int(case.opts[case.opts.index("-codontable") + 1]) if "-codontable" in case.opts else 1
    if "-nogap" in case.opts and "invalid_char" not in case.features:
        rc0, out0 = coal_without_nogap()
        why = check_deliberate(_without_nogap(case), fmt, rc0, out0, ref_out)
        if why is not None:
            return "without -nogap: " + why
        if (rc0 == 0) != (rc == 0):
            return "COAL succeeds with or without -nogap, not both"
        if rc != 0:
            return None if out == "" else "COAL failed but wrote output"
        try:
            return None if parse_output(fmt, out) == _nogap_rule(fmt, out0, table) else (
                "-nogap output is not the documented rule applied to the output without -nogap")
        except (ValueError, IndexError) as e:
            return f"could not parse output ({e})"
    if "invalid_char" in case.features:
        return None if rc != 0 and out == "" else "invalid protein character but COAL did not refuse"
    if ref_out == "":
        return None if rc != 0 and out == "" else "PAL2NAL fails on the neutralised case but COAL succeeded"
    if rc != 0:
        return f"COAL failed (exit {rc}) but PAL2NAL succeeds on the neutralised case"
    try:
        want, got = parse_output(fmt, ref_out), parse_output(fmt, out)
    except (ValueError, IndexError) as e:
        return f"could not parse output ({e})"
    if [i for i, _ in want] != [i for i, _ in got]:
        return "different IDs or order"
    for (sid, w), (_, g) in zip(want, got):
        if len(w) != len(g):
            return f"{sid}: length {len(g)}, expected {len(w)}"
        for k in range(0, len(w), 3):
            c3, r3 = g[k:k + 3], w[k:k + 3]
            if c3.replace("U", "T").replace("u", "t") == r3:
                continue
            marked = [a.isupper() != b.isupper() for a, b in zip(c3, r3) if a.isalpha() and b.isalpha()]
            if any(marked) and _allowed_codon(table, c3, r3):
                continue
            return f"{sid} at nt {k + 1}: COAL {c3}, expected {r3}"
    return None


def _fasta_records(text):
    """[(header line, sequence)] of a FASTA text (lines joined)."""
    recs = []
    for ln in text.split("\n"):
        if ln.startswith(">"):
            recs.append([ln, ""])
        elif recs:
            recs[-1][1] += ln.strip()
    return recs


def placed_first_met_files(case, coal_fasta):
    """A third neutralised copy (from the copy without first-Met substitutions), built
    after the run: for every gene that starts with M,
    find where COAL placed it (its ungapped output row is a contiguous piece of its DNA) and, if that first codon is
    a start or Met codon of the selected code other than ATG, write a marked ATG there. Needed when a gene is placed
    off the generator's own frame (e.g. an unrelated DNA), where the first Met falls on another codon. None if
    nothing would change."""
    table = int(case.opts[case.opts.index("-codontable") + 1]) if "-codontable" in case.opts else 1
    pfmt = "clustal" if case.pep.endswith(".aln") else "fasta"
    prot = {i: r.upper() for i, r in parse_output(pfmt, case.files[case.pep])}
    files = dict(case.ref_files_b)  # no first-Met substitution in the generator's frame: only where COAL placed it
    changed = False
    for sid, row in parse_output("fasta", coal_fasta):
        prow = prot.get(sid, "")
        first = next((ch for ch in prow if ch != "-"), "")
        ungapped = row.replace("-", "")
        if first != "M" or len(ungapped) < 3:
            continue
        cod = ungapped[:3].upper().replace("U", "T")
        if cod == "ATG" or cod not in FIRST_MET(table):
            continue
        for name in case.nucs:
            real = _fasta_records(case.files[name])
            ref = _fasta_records(files[name])
            for k, (hdr, seq) in enumerate(real):
                pos = seq.find(ungapped)
                if pos < 0:
                    continue
                r = ref[k][1]
                atg = "atg" if seq[pos:pos + 3].isupper() else "ATG"  # the other case than the real DNA
                ref[k][1] = r[:pos] + atg + r[pos + 3:]
                files[name] = "".join(f"{h}\n{q}\n" for h, q in ref)
                changed = True
                break
            else:
                continue
            break
    return files if changed else None


def _without_nogap(case):
    c = Case()
    c.__dict__.update(case.__dict__)
    c.opts = [o for o in case.opts if o != "-nogap"]
    return c


def run_pal2nal(pal2nal, d, args):
    p = subprocess.run(["perl", pal2nal, *args], cwd=d, capture_output=True, timeout=120)
    return p.returncode, p.stdout.decode("latin-1"), p.stderr.decode("latin-1", "replace")


def run_coal(d, args):
    out, err = io.StringIO(), io.StringIO()
    cwd = os.getcwd()
    os.chdir(d)
    try:
        with redirect_stdout(out), redirect_stderr(err):
            try:
                rc = pal2nal_main(list(args))
            except SystemExit as e:
                rc = e.code
    finally:
        os.chdir(cwd)
    return rc, out.getvalue(), err.getvalue()


def cpp_binary():
    """coal-cpp from COAL_CPP_BIN (default build/coal-cpp; .exe and build/Release/ tried too), or None."""
    p = Path(os.environ.get("COAL_CPP_BIN", ROOT / "build" / "coal-cpp"))
    if not p.is_absolute():
        p = ROOT / p
    for cand in (p, p.with_suffix(".exe"), p.parent / "Release" / (p.name + ".exe")):
        if cand.is_file():
            return cand
    return None


def cpp_run_coal(binary):
    """A `run_coal` that runs coal-cpp with the PAL2NAL command line translated as the drop-in does."""
    from coal.pipeline import FORMATS as ENGINE_FORMATS
    from coal_compat.cli import UsageError, translate

    def run(d, args):
        try:
            cargs = translate(list(args), ENGINE_FORMATS)
        except UsageError as e:
            return 2, "", f"pal2nal: error: {e}\n"
        r = subprocess.run([str(binary), "--max-warnings", "0", *cargs], cwd=d, capture_output=True, timeout=120)
        err = r.stderr.decode("utf-8", "replace").replace("\r\n", "\n")
        err = "".join(ln for ln in err.splitlines(True) if not ln.startswith("coal: summary:"))
        return r.returncode, r.stdout.decode("latin-1"), err
    return run


def select_engine(name=None):
    """Point `run_coal` at the engine named by COAL_FUZZ_ENGINE ("python", the default, or "cpp"); returns the
    name. "cpp" without a binary raises FileNotFoundError."""
    global run_coal
    name = (name or os.environ.get("COAL_FUZZ_ENGINE") or "python").lower()
    if name == "python":
        run_coal = _python_run_coal
    elif name == "cpp":
        binary = cpp_binary()
        if binary is None:
            raise FileNotFoundError("COAL_FUZZ_ENGINE=cpp but coal-cpp is not built (build/coal-cpp or COAL_CPP_BIN)")
        run_coal = cpp_run_coal(binary)
    else:
        raise ValueError(f"COAL_FUZZ_ENGINE={name!r}: expected 'python' or 'cpp'")
    return name


_python_run_coal = run_coal


def fuzz(pal2nal, n, seed, report=print, keep_dir=None, formats=FORMATS):
    """Run n cases; return (counts per format, list of unexplained differences)."""
    rng = random.Random(seed)
    counts = {f: {"cases": 0, "identical": 0, "both_fail": 0, "deliberate": 0, "differ": 0} for f in formats}
    by_feature = {}
    coverage = {}
    bad = []
    base = tempfile.mkdtemp(prefix="coal_fuzz_")
    try:
        for i in range(n):
            case = make_case(rng)
            for t in case.info:
                coverage[t] = coverage.get(t, 0) + 1
            d = os.path.join(base, f"c{i}")
            os.mkdir(d)
            for name, text in case.files.items():
                with open(os.path.join(d, name), "w", newline="\n") as fh:
                    fh.write(text)
            procs = {}
            for fmt in formats:
                args = [case.pep, *case.nucs, "-output", fmt, *case.opts]
                procs[fmt] = (args, subprocess.Popen(["perl", pal2nal, *args], cwd=d, stdout=subprocess.PIPE,
                                                     stderr=subprocess.PIPE))
            keep = False
            refs = {}

            def reference(fmt, variant, case=case, d=d, refs=refs):
                """PAL2NAL's stdout on a neutralised copy of the case, without -nogap; variant 3: the case
                itself with only the DNA records reordered, with every option."""
                if (fmt, variant) not in refs:
                    rd = os.path.join(d, f"ref{variant}")
                    os.makedirs(rd, exist_ok=True)
                    for name, text in (case.ref_files, case.ref_files_b, getattr(case, "ref_files_c", None),
                                       case.ref_files_p)[variant].items():  # variant 3: ref_files_p
                        with open(os.path.join(rd, name), "w", newline="\n") as fh:
                            fh.write(text)
                    opts = [o for o in case.opts if o != "-nogap" or variant == 3]
                    refs[fmt, variant] = run_pal2nal(pal2nal, rd, [case.pep, *case.nucs, "-output", fmt, *opts])[1]
                return refs[fmt, variant]

            def verified(fmt, rc, out, case=case, d=d):
                """None if the difference is exactly a documented one (either reference), else why not."""
                if "normalised_ids" in case.features:  # exact, byte for byte
                    ref = reference(fmt, 3)
                    if (rc == 0 and out == ref != "") or (ref == "" and rc != 0 and out == ""):
                        return None
                    if case.features & set(DELIBERATE) == {"normalised_ids"}:
                        return "normalised IDs: COAL's output is not PAL2NAL's on the DNA in protein order"
                    # with other deliberate features too: the neutralised references below (DNA also in protein order)
                def base():
                    a = [case.pep, *case.nucs, "-output", fmt, *[o for o in case.opts if o != "-nogap"]]
                    r = run_coal(d, a)
                    return r[0], r[1]
                why = "no deliberate feature"
                if case.features & set(DELIBERATE):  # a case without one must be byte-identical (or both fail)
                    why = check_deliberate(case, fmt, rc, out, reference(fmt, 0), base)
                if why is not None and case.features & set(DELIBERATE) and case.ref_files_b != case.ref_files:
                    if check_deliberate(case, fmt, rc, out, reference(fmt, 1), base) is None:
                        return None
                if why is not None:  # the first Met where COAL placed the gene
                    if not hasattr(case, "ref_files_c"):
                        a = [case.pep, *case.nucs, "-output", "fasta", *[o for o in case.opts
                                                                        if o not in ("-nogap", "-nomismatch")]]
                        r = run_coal(d, a)
                        case.ref_files_c = placed_first_met_files(case, r[1]) if r[0] == 0 else None
                    if case.ref_files_c is not None:
                        if check_deliberate(case, fmt, rc, out, reference(fmt, 2), base) is None:
                            case.features.add("first_met_code")
                            return None
                return why
            why = None
            for fmt in formats:
                args, p = procs[fmt]
                rc, out, err = run_coal(d, args)
                try:
                    pout, perr = p.communicate(timeout=120)
                except subprocess.TimeoutExpired:
                    p.kill()
                    pout, perr = p.communicate()
                pout = pout.decode("latin-1")
                cnt = counts[fmt]
                cnt["cases"] += 1
                pal_ok = pout != ""
                if pal_ok and rc == 0 and out == pout:
                    cnt["identical"] += 1
                elif not pal_ok and rc != 0 and out == "":
                    cnt["both_fail"] += 1
                elif (why := verified(fmt, rc, out)) is None:  # every difference is checked
                    cnt["deliberate"] += 1
                    for f in case.features & set(DELIBERATE):
                        by_feature[f] = by_feature.get(f, 0) + 1
                else:
                    cnt["differ"] += 1
                    keep = True
                    bad.append({"case": i, "dir": d, "args": args, "pal2nal_ok": pal_ok, "coal_rc": rc,
                                "features": sorted(case.features),
                                "why": why,
                                "pal2nal_stdout": pout[:400], "coal_stdout": out[:400],
                                "coal_stderr": err[-400:], "pal2nal_stderr": perr.decode("latin-1")[-400:]})
            if keep and keep_dir:
                shutil.copytree(d, os.path.join(keep_dir, f"c{i}"), dirs_exist_ok=True)
                with open(os.path.join(keep_dir, f"c{i}", "args.txt"), "w") as fh:
                    fh.write(" ".join(bad[-1]["args"]) + "\n")
            shutil.rmtree(d, ignore_errors=True)
    finally:
        shutil.rmtree(base, ignore_errors=True)
    for fmt in formats:
        c = counts[fmt]
        report(f"fuzz vs PAL2NAL v14, -output {fmt}: {c['cases']} cases, identical {c['identical']}, both fail "
               f"{c['both_fail']}, deliberate differences {c['deliberate']}, unexplained differences {c['differ']}")
    if by_feature:
        report("deliberate differences by feature: " + ", ".join(f"{k} {v}" for k, v in sorted(by_feature.items())))
    report("cases exercising: " + ", ".join(f"{k} {v}" for k, v in sorted(coverage.items())))
    return counts, bad


@unittest.skipUnless(os.environ.get("PAL2NAL"), "set PAL2NAL=/path/to/pal2nal.pl (scripts/fetch_pal2nal.py)")
class TestFuzzVsPal2nal(unittest.TestCase):
    def test_random_inputs_identical(self):
        pal2nal = os.path.abspath(os.environ["PAL2NAL"])
        self.assertTrue(os.path.isfile(pal2nal), pal2nal)
        self.assertIsNotNone(shutil.which("perl"), "perl is needed to run PAL2NAL")
        n = int(os.environ.get("FUZZ_N", "3000"))
        seed = int(os.environ.get("FUZZ_SEED", "20261017"))
        formats = tuple(os.environ.get("FUZZ_FORMATS", ",".join(FORMATS)).split(","))  # e.g. fasta,paml
        engine = select_engine()  # COAL_FUZZ_ENGINE: python (default) or cpp
        print(f"fuzz engine {engine}, seed {seed}, {n} cases", file=sys.stderr)
        counts, bad = fuzz(pal2nal, n, seed, report=lambda s: print(s, file=sys.stderr),
                           keep_dir=os.environ.get("FUZZ_KEEP"), formats=formats)
        for b in bad[:5]:
            print(b, file=sys.stderr)
        self.assertEqual(bad, [], f"{len(bad)} unexplained differences (first: {bad[:1]})")
        self.assertEqual(sum(c["cases"] for c in counts.values()), len(formats) * n)


if __name__ == "__main__":
    unittest.main()
