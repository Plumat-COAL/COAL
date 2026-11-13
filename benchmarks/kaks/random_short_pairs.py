"""Random short pairs against KaKs_Calculator 2.0: the cases real genes rarely reach (docs/kaks.md, "Short and
saturated pairs").

  python benchmarks/kaks/random_short_pairs.py --coal-cpp build/coal-cpp --kaks-hp KaKs_Calculator_15digits \
      [--mode single|multi] [--n 4000] [--seed 5] [--methods YN,MYN] [--kaks-fixed KaKs_Calculator_corrected] \
      [--yn00 yn00_12digits] [--one-file] --out work/

Two generators:

  --mode single (default; seed 5, 4,000 pairs)   2 to 8 codons over full or reduced alphabets; the second
      sequence is the first with at most ONE changed position per codon. Reaches: sites of one kind only, base
      frequencies of zero, a kappa that cannot be estimated. No pathway is ever weighted.
  --mode multi (seed 20261112, 3,000 pairs)      1 to 12 codons from ten codon pools; the second sequence is
      drawn independently, codon by codon, so codons differ at two or three positions. Reaches in addition:
      pathway weighting with NaN inside, saturation of one kind of site, YN / MYN iterations that do not
      converge. (The generator of the independent review of 2026-11-13.)

Both tools are run with 15 significant digits. KaKs_Calculator is run ONCE PER PAIR, each in a process of its
own: its MYN carries state from one pair to the next, so in one file its output for a pair can depend on the
pairs before it (--one-file runs it the old way and reports how many rows that changes).

Every field (Ka, Ks, Ka/Ks, S-Sites, N-Sites) where the two outputs differ is put in a class:

  nan -> NA               KaKs_Calculator prints nan, COAL NA
  nan -> number           KaKs_Calculator loses a defined value to NaN arithmetic (differences 8 to 10)
  NA -> number            KaKs_Calculator prints NA, COAL a number (difference 8)
  number -> NA            allowed only with COAL's status `saturated` (difference 1)
  number != number        two numbers that differ by more than 1e-9 (relative). Each must be explained:
      one kind of site      COAL status no_data (S = 0 or N = 0): KaKs_Calculator drops the differences of codons
                            that differ at several positions (difference 11); with --yn00, COAL's Ka (Ks) must be
                            yn00's with weighting = 0
      reference nan row     the KaKs_Calculator row has nan / NA in Ka or Ks: a NaN got into its iteration and
                            the number it still prints lacks the same codons (difference 11)
      not converged         COAL's note says the iteration did not converge (difference 12)
      saturated, other path COAL status saturated and the difference is below 1e-6 (difference 13)
      unexplained           anything else: the script exits 1

With --kaks-fixed (the same source with the two corrected lines of docs/kaks.md) the same table is made against
that build.
"""
from __future__ import annotations

import argparse
import math
import random
import subprocess
import sys
import tempfile
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

BASES = "TCAG"
CODE1 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
AA = {BASES[i // 16] + BASES[(i // 4) % 4] + BASES[i % 4]: CODE1[i] for i in range(64)}
SENSE = [c for c in AA if AA[c] != "*"]
FIELDS = {2: "Ka", 3: "Ks", 4: "Ka/Ks", 7: "S-Sites", 8: "N-Sites"}


def make_pairs(n: int, seed: int):
    """One changed position per codon at most (the generator of 2026-11-12)."""
    rng = random.Random(seed)
    out = []
    for k in range(n):
        ncod = rng.randint(2, 8)
        alpha = rng.choice(["TCAG", "TCA", "TCG", "CAG", "TAG", "TC", "AG"])
        pool = [c for c in SENSE if all(b in alpha for b in c)]
        if len(pool) < 2:
            continue
        a = [rng.choice(pool) for _ in range(ncod)]
        b = list(a)
        for i in range(ncod):
            if rng.random() < 0.5:
                pos = rng.randrange(3)
                c = a[i][:pos] + rng.choice([x for x in alpha if x != a[i][pos]]) + a[i][pos + 1:]
                if AA[c] != "*":
                    b[i] = c
        out.append((f"r{k}", "".join(a), "".join(b)))
    return out


def make_pairs_multi(n: int, seed: int):
    """Codons of the second sequence drawn independently: differences at two and three positions."""
    rng = random.Random(seed)
    pools = [SENSE, ["ATG", "TGG"], ["AAA", "AAG"], ["GCT", "GCA", "GCC", "GCG"],
             [c for c in SENSE if set(c) <= set("TC")], [c for c in SENSE if set(c) <= set("AG")],
             [c for c in SENSE if set(c) <= set("TCA")], [c for c in SENSE if set(c) <= set("CAG")],
             ["ATG", "TGG", "AAA", "AAG", "GCT", "GCA"], ["TTT", "TTC", "CTT", "CTC"]]
    out = []
    for k in range(n):
        ncod = rng.randint(1, 12)
        pool = rng.choice(pools)
        a = [rng.choice(pool) for _ in range(ncod)]
        mode = rng.random()
        if mode < 0.1:
            b = list(a)
        elif mode < 0.6:
            b = [c if rng.random() < 0.6 else rng.choice(pool) for c in a]
        else:
            b = [rng.choice(pool) for _ in range(ncod)]
        out.append((f"p{k}", "".join(a), "".join(b)))
    return out


def num(x: str):
    try:
        v = float(x)
    except ValueError:
        return None
    return v if math.isfinite(v) else None


def rows(text: str) -> dict:
    out = {}
    for line in text.splitlines():
        f = line.split("\t")
        if len(f) > 8 and f[0] != "Sequence":
            out[(f[0], f[1])] = f
    return out


def run_reference(exe: str, pairs, method: str, one_file: bool, workdir: Path) -> str:
    """KaKs_Calculator's rows for `method`: one process per pair (default), or one run over the whole file."""
    if one_file:
        axt, out = workdir / "all.axt", workdir / f"all.{method}.kaks"
        axt.write_text("".join(f"{n}\n{s1}\n{s2}\n\n" for n, s1, s2 in pairs))
        subprocess.run([exe, "-i", str(axt), "-o", str(out), "-m", method], capture_output=True, check=True)
        return out.read_text()

    def one(k: int) -> str:
        n, s1, s2 = pairs[k]
        axt, out = workdir / f"{k}.axt", workdir / f"{k}.kaks"
        axt.write_text(f"{n}\n{s1}\n{s2}\n\n")
        subprocess.run([exe, "-i", str(axt), "-o", str(out), "-m", method], capture_output=True)
        text = out.read_text() if out.exists() else ""
        axt.unlink()
        out.unlink(missing_ok=True)
        return "\n".join(text.splitlines()[1:])
    with ThreadPoolExecutor(2) as ex:
        return "\n".join(ex.map(one, range(len(pairs))))


def yn00_w0(yn00: str, s1: str, s2: str):
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from validate_kaks import yn00_one
    return yn00_one(yn00, s1, s2, 0)


def classify(ref: dict, ours: dict, status: dict, seqs: dict, yn00: str | None):
    counts, example, bad = Counter(), {}, []
    for key, rf in ref.items():
        of, (st, note) = ours[key], status[key]
        ref_nan = any(num(rf[i]) is None for i in (2, 3))
        for i, name in FIELDS.items():
            x, y = rf[i], of[i]
            if x == y:
                continue
            nx, ny = num(x), num(y)
            if nx is not None and ny is not None:
                rel = abs(nx - ny) / max(1.0, abs(nx), abs(ny))
                if rel <= 1e-9:
                    continue
                if "did not converge" in note:
                    kind = "number != number: not converged"
                elif st == "no_data":
                    kind = "number != number: one kind of site"
                    if yn00 and name in ("Ka", "Ks"):
                        y0 = yn00_w0(yn00, *seqs[key[0]])
                        v = y0 and y0["dN" if name == "Ka" else "dS"]
                        if not (v is not None and math.isfinite(v) and abs(v - ny) <= 1e-9 * max(1.0, abs(ny))):
                            kind = "number != number: unexplained (one kind of site, but not yn00 weighting 0)"
                elif st == "saturated" and rel <= 1e-6:
                    kind = "number != number: saturated, other path"
                elif ref_nan:
                    kind = "number != number: reference nan row"
                else:
                    kind = "number != number: unexplained"
            elif nx is not None:
                kind = "number -> NA"
                if st != "saturated" and nx != 0:
                    kind = "number -> NA: unexplained"
            elif ny is not None:
                kind = "nan -> number" if "nan" in x or "inf" in x else f"{x} -> number"
            else:
                kind = "nan -> NA"
            if "unexplained" in kind:
                bad.append((key, name, x, y, st))
            counts[(key[1], name, kind, st)] += 1
            example.setdefault((key[1], name, kind, st), f"{key[0]}: {x} / {y}")
    return counts, example, bad


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--coal-cpp", required=True)
    ap.add_argument("--kaks-hp", required=True, help="KaKs_Calculator 2.0 built to print 15 significant digits")
    ap.add_argument("--kaks-fixed", help="the same with the two corrected lines of docs/kaks.md (optional)")
    ap.add_argument("--yn00", help="PAML yn00 printing 12 digits (optional): checks the `one kind of site` class")
    ap.add_argument("--mode", choices=("single", "multi"), default="single")
    ap.add_argument("--n", type=int)
    ap.add_argument("--seed", type=int)
    ap.add_argument("--methods", default="NG,LWL,MLWL,LPB,MLPB,YN,MYN")
    ap.add_argument("--one-file", action="store_true",
                    help="also run KaKs_Calculator once over the whole file and count the rows that differ")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    methods = a.methods.split(",")
    if a.mode == "single":
        pairs = make_pairs(a.n or 4000, 5 if a.seed is None else a.seed)
    else:
        pairs = make_pairs_multi(a.n or 3000, 20261112 if a.seed is None else a.seed)
    seqs = {n: (s1, s2) for n, s1, s2 in pairs}
    axt = out / "random_short.axt"
    axt.write_text("".join(f"{n}\n{s1}\n{s2}\n\n" for n, s1, s2 in pairs))
    mflags = sum((["-m", m] for m in methods), [])
    subprocess.run([a.coal_cpp, "kaks", "-i", str(axt), "-o", str(out / "coal.kaks"), "--digits", "15", "-q",
                    *mflags], check=True)
    subprocess.run([a.coal_cpp, "dnds", str(axt), "-o", str(out / "coal.dnds.tsv"), "--digits", "15", "-q",
                    *mflags], check=True)
    status = {}
    for line in (out / "coal.dnds.tsv").read_text().splitlines()[1:]:
        f = line.split("\t")
        status[(f[0], f[2])] = (f[13], f[14] if len(f) > 14 else "")
    ours = rows((out / "coal.kaks").read_text())
    total_bad = 0
    builds = [("KaKs_Calculator", a.kaks_hp)] + ([("KaKs_Calculator_corrected", a.kaks_fixed)] if a.kaks_fixed else [])
    with tempfile.TemporaryDirectory() as td:
        for label, exe in builds:
            ref = {}
            for m in methods:
                text = run_reference(exe, pairs, m, False, Path(td))
                (out / f"{label}.{m}.kaks").write_text(text + "\n")
                ref.update(rows(text))
                if a.one_file:
                    whole = rows(run_reference(exe, pairs, m, True, Path(td)))
                    moved = sum(1 for k, f in rows(text).items() if k in whole and whole[k][2:5] != f[2:5])
                    print(f"{label} {m}: {moved} of {len(whole)} rows differ (Ka, Ks or Ka/Ks) between one run "
                          "over the file and one process per pair")
            counts, example, bad = classify(ref, ours, status, seqs, a.yn00)
            table = out / (f"random_short_pairs.{a.mode}.tsv" if label == "KaKs_Calculator"
                           else f"random_short_pairs.{a.mode}.corrected.tsv")
            with open(table, "w", newline="\n") as fh:
                fh.write("method\tfield\tdifference\tcoal_status\tcount\texample (pair: KaKs_Calculator / COAL)\n")
                for k in sorted(counts):
                    fh.write("\t".join(k) + f"\t{counts[k]}\t{example[k]}\n")
            nn = Counter()
            for (m, name, kind, st), c in counts.items():
                if kind.startswith("number != number"):
                    nn[kind.split(": ", 1)[1]] += c
            print(f"{label}: {len(pairs)} pairs ({a.mode}), methods {','.join(methods)}: {sum(counts.values())} "
                  f"differing fields in {len(counts)} classes; number != number: "
                  + (", ".join(f"{c} {k}" for k, c in sorted(nn.items())) or "none")
                  + f"; {len(bad)} unexplained -> {table}")
            for b in bad[:10]:
                print("  unexplained:", b)
            total_bad += len(bad)
    return 1 if total_bad else 0


if __name__ == "__main__":
    sys.exit(main())
