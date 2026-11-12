"""Random short pairs against KaKs_Calculator 2.0: the edge cases real genes rarely reach (docs/kaks.md).

  python benchmarks/kaks/random_short_pairs.py --coal-cpp build/coal-cpp --kaks-hp KaKs_Calculator_15digits \
      [--n 4000] [--seed 5] [--methods YN,MYN] --out work/

Makes N pairs of 2 to 8 codons over full or reduced alphabets (so that sites of one kind, base frequencies of
zero and kappa that cannot be estimated all occur), runs both tools with 15 significant digits and classifies
every difference in Ka, Ks, Ka/Ks, S-Sites and N-Sites:

  number != number     never allowed (the script exits 1)
  number -> NA         allowed only with COAL's status `saturated` (difference 1 of docs/kaks.md)
  nan -> NA            KaKs_Calculator prints nan, COAL NA
  nan -> number        KaKs_Calculator loses a defined value to NaN arithmetic (differences 8 to 10)

With a KaKs_Calculator built with the two corrections named in docs/kaks.md (after difference 10), the
`nan -> number` class is empty for MYN (run with --methods MYN and that build as --kaks-hp).
"""
from __future__ import annotations

import argparse
import math
import random
import subprocess
import sys
from collections import Counter
from pathlib import Path

BASES = "TCAG"
CODE1 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
AA = {BASES[i // 16] + BASES[(i // 4) % 4] + BASES[i % 4]: CODE1[i] for i in range(64)}
SENSE = [c for c in AA if AA[c] != "*"]
FIELDS = {2: "Ka", 3: "Ks", 4: "Ka/Ks", 7: "S-Sites", 8: "N-Sites"}


def make_pairs(n: int, seed: int):
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


def num(x: str):
    try:
        v = float(x)
    except ValueError:
        return None
    return v if math.isfinite(v) else None


def rows(path: Path) -> dict:
    out = {}
    for line in path.read_text().splitlines()[1:]:
        f = line.split("\t")
        if len(f) > 8:
            out[(f[0], f[1])] = f
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--coal-cpp", required=True)
    ap.add_argument("--kaks-hp", required=True)
    ap.add_argument("--n", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=5)
    ap.add_argument("--methods", default="NG,LWL,MLWL,LPB,MLPB,YN,MYN")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    methods = a.methods.split(",")
    axt = out / "random_short.axt"
    pairs = make_pairs(a.n, a.seed)
    axt.write_text("".join(f"{n}\n{s1}\n{s2}\n\n" for n, s1, s2 in pairs))
    mflags = sum((["-m", m] for m in methods), [])
    subprocess.run([a.coal_cpp, "kaks", "-i", str(axt), "-o", str(out / "coal.kaks"), "--digits", "15", "-q",
                    *mflags], check=True)
    subprocess.run([a.coal_cpp, "dnds", str(axt), "-o", str(out / "coal.dnds.tsv"), "--digits", "15", "-q",
                    *mflags], check=True)
    status = {}
    for line in (out / "coal.dnds.tsv").read_text().splitlines()[1:]:
        f = line.split("\t")
        status[(f[0], f[2])] = f[13]
    ours = rows(out / "coal.kaks")
    counts, example, bad = Counter(), {}, 0
    for m in methods:
        ref = out / f"KaKs_Calculator.{m}.kaks"
        subprocess.run([a.kaks_hp, "-i", str(axt), "-o", str(ref), "-m", m], capture_output=True, check=True)
        for key, rf in rows(ref).items():
            of, st = ours[key], status[key]
            for i, name in FIELDS.items():
                x, y = rf[i], of[i]
                if x == y:
                    continue
                nx, ny = num(x), num(y)
                if nx is not None and ny is not None:
                    if abs(nx - ny) <= 1e-9 * max(1.0, abs(nx), abs(ny)):
                        continue
                    kind = "number != number"
                elif nx is not None:
                    kind = "number -> NA"
                elif ny is not None:
                    kind = "nan -> number" if "nan" in x or "inf" in x else f"{x} -> number"
                else:
                    kind = "nan -> NA"
                if kind == "number != number" or (kind == "number -> NA" and st != "saturated" and nx != 0):
                    bad += 1
                counts[(m, name, kind, st)] += 1
                example.setdefault((m, name, kind, st), f"{key[0]}: {x} / {y}")
    with open(out / "random_short_pairs.tsv", "w", newline="\n") as fh:
        fh.write("method\tfield\tdifference\tcoal_status\tcount\texample (pair: KaKs_Calculator / COAL)\n")
        for k in sorted(counts):
            fh.write("\t".join(k) + f"\t{counts[k]}\t{example[k]}\n")
    print(f"{len(pairs)} pairs, methods {','.join(methods)}: {sum(counts.values())} differing fields in "
          f"{len(counts)} classes, {bad} not allowed -> {out / 'random_short_pairs.tsv'}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
