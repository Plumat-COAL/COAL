"""Differential fuzz: COAL's ``pal2nal`` drop-in vs the real PAL2NAL v14, on random inputs.

Needs Perl and PAL2NAL v14 (fetch it with ``python scripts/fetch_pal2nal.py``; it is GPL and never committed):

    PAL2NAL=PAL2NAL/pal2nal.v14/pal2nal.pl python tests/fuzz_vs_pal2nal.py          # FUZZ_N cases (default 3000)
    PAL2NAL=... FUZZ_N=500 FUZZ_SEED=7 python tests/fuzz_vs_pal2nal.py

Without PAL2NAL set the test is skipped. Every case is a random protein alignment (FASTA or CLUSTAL, sometimes with
a '#'-marked block line) plus coding DNA (UTRs, stop codons, mismatching / unrelated / masked / lower-case / RNA DNA,
frameshift digits, other genetic codes, shuffled / different / partly matching IDs ...) and a random option set
(-codontable, -nogap, -nomismatch, -blockonly). Each case runs through both tools with exactly the same command
line, once for each output format (-output clustal / paml / fasta / codon), so 4 comparisons per case.

Rule per comparison: PAL2NAL always exits 0; it succeeded if its stdout is non-empty. Then COAL (the ``pal2nal``
command, run in-process) must exit 0 with byte-identical stdout; if PAL2NAL failed, COAL must exit non-zero with
empty stdout. stderr is not compared (COAL's messages are its own). The only exceptions are the deliberate
differences in DELIBERATE below (each is in README.md's list): a case that has one of those features may differ,
and is counted under that feature instead.
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
from coal.tables import SUPPORTED, get_code  # noqa: E402
from coal_compat.cli import main as pal2nal_main  # noqa: E402

FORMATS = ("clustal", "paml", "fasta", "codon")

# Deliberate differences (README.md, "Deliberate differences"): feature -> why COAL may differ from PAL2NAL.
# A case is allowed to differ only if the generator gave it one of these features.
DELIBERATE = {
    "invalid_char": "invalid protein characters (e.g. '@', '#'): PAL2NAL drops a codon silently, COAL refuses (exit 1)",
    "genetic_code": "COAL uses the current NCBI tables, PAL2NAL v14 older ones. Measured against the real PAL2NAL "
                    "(every codon x residue, every first-Met codon): table 10 TGA (NCBI Cys, PAL2NAL stop) and the "
                    "first-Met start codons GTG in table 3, ATT/ATC/ATA in table 11 and ATA in table 13 (NCBI "
                    "starts PAL2NAL lacks). Nothing else differs",
    "iupac": "IUPAC codes other than N in the DNA: COAL accepts a codon if every reading gives the residue; "
             "PAL2NAL has a fixed list",
    "u_residue": "U (selenocysteine) in the protein: COAL treats it as X (any codon, with a warning) like B/Z/J/O; "
                 "PAL2NAL's search pattern gives U no codon at all, so the row is misaligned or the gene fails",
    "nogap_code": "-nogap removes the stop codons of the selected genetic code, in any case (COAL, deliberate); "
                  "PAL2NAL always removes the standard code's TAA/TAG/TGA, upper case only",
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
        self.pep = None
        self.nucs = []
        self.opts = []
        self.features = set()
        self.info = set()  # what the case exercises (reported as coverage counts)


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
    if rng.random() < 0.04:  # tfastx/tfasty frameshift notation: '\\' = 1 base, '/' + residue = gap + 2 bases
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
    for k, row in enumerate(rows):
        out = []
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
        if start_idx is not None and cds[start_idx].upper().replace("U", "T") in START_CODONS.get(run_table, ()):
            c.features.add("genetic_code")
        seq = "".join(cds)
        if rng.random() < 0.3:
            seq = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 30))) + seq
        if rng.random() < 0.3:
            seq += rng.choice(("TAA", "TAG", "TGA")) + "".join(rng.choice("ACGT") for _ in range(rng.randint(0, 20)))
        if rng.random() < 0.03 and len(seq) > 3:
            seq = seq[:-rng.randint(1, 3)]
        if rng.random() < 0.04:
            seq = "".join(rng.choice("RYSWKMBDHV") if b in "ACGT" and rng.random() < 0.05 else b for b in seq)
            c.features.add("iupac")
        if rng.random() < 0.05:
            seq = seq.replace("T", "U")
            if start_idx is not None:
                c.features.add("rna_first_met")
        if rng.random() < 0.1:
            seq = seq.lower()
        dnas.append(seq)

    if any(t in d.upper().replace("U", "T") for d in dnas for t in TABLE_CODONS.get(run_table, ())):
        c.features.add("genetic_code")

    # DNA IDs and order
    did_mode = rng.choices(("same", "shuffled", "different", "partial", "count"), (0.65, 0.17, 0.1, 0.05, 0.03))[0]
    order = list(range(nseq))
    dids = list(pids)
    if did_mode == "shuffled":
        rng.shuffle(order)
    elif did_mode == "different":
        dids = [f"d{k}x" for k in range(nseq)]
    elif did_mode == "partial" and nseq >= 2:
        dids = [p if k == 0 else f"other{k}" for k, p in enumerate(pids)]
    recs = [(dids[j], dnas[j]) for j in order]
    if did_mode == "count":
        recs = recs[:-1] if (len(recs) > 1 and rng.random() < 0.5) else recs + [("extra", "ATGAAA")]
    width = rng.choice((0, 60, 60, 70, 13))
    nuc_parts = ["".join(f">{i} desc\n" + _wrap(s, width) for i, s in recs)]
    if len(recs) > 1 and rng.random() < 0.1:
        cut = rng.randint(1, len(recs) - 1)
        nuc_parts = ["".join(f">{i}\n" + _wrap(s, width) for i, s in recs[:cut]),
                     "".join(f">{i}\n" + _wrap(s, width) for i, s in recs[cut:])]

    # protein file
    lower = rng.random() < 0.1
    srows = ["".join(r) for r in rows]
    if lower:
        srows = [s.lower() for s in srows]
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
        text = "CLUSTAL W (1.83) multiple sequence alignment\n\n\n"
        for b in range(0, ncol, blk):
            for p, s in zip(pids, srows):
                text += p.ljust(idw) + s[b:b + blk] + "\n"
            last = b + blk >= ncol
            if mask is not None:
                text += " " * idw + mask[b:b + blk].rstrip() + "\n"
            else:
                text += " " * idw + "".join(rng.choice(" *:.") for _ in s[b:b + blk]).rstrip() + "\n"
            if not last or rng.random() < 0.8:
                text += "\n"
        c.pep = "pep.aln"
    else:
        text = "".join(f">{p}\n" + _wrap(s, rng.choice((0, 60, 10))) for p, s in zip(pids, srows))
        c.pep = "pep.fasta"
    c.files[c.pep] = text
    for i, part in enumerate(nuc_parts):
        c.files[f"nuc{i}.fasta"] = part
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
                elif case.features & set(DELIBERATE):
                    cnt["deliberate"] += 1
                    for f in case.features & set(DELIBERATE):
                        by_feature[f] = by_feature.get(f, 0) + 1
                else:
                    cnt["differ"] += 1
                    keep = True
                    bad.append({"case": i, "dir": d, "args": args, "pal2nal_ok": pal_ok, "coal_rc": rc,
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
        counts, bad = fuzz(pal2nal, n, seed, report=lambda s: print(s, file=sys.stderr),
                           keep_dir=os.environ.get("FUZZ_KEEP"), formats=formats)
        for b in bad[:5]:
            print(b, file=sys.stderr)
        self.assertEqual(bad, [], f"{len(bad)} unexplained differences (first: {bad[:1]})")
        self.assertEqual(sum(c["cases"] for c in counts.values()), len(formats) * n)


if __name__ == "__main__":
    unittest.main()
