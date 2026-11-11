"""Make the Ka/Ks test inputs and golden files (tests/test_kaks.py; docs/kaks.md).

  python tests/kaks/make_fixtures.py                        # inputs only: edge_cases.axt, sim_small.axt
  python tests/kaks/make_fixtures.py --coal-cpp build/coal-cpp \
      --kaks KaKs_Calculator_15digits --yn00 yn00_12digits    # + reference and COAL golden files

All sequences are made here (edge cases by hand, the rest by a seeded simulator), so the data are COAL's own.
Reference files: golden/<set>.KaKs_Calculator.<M>.kaks (KaKs_Calculator 2.0 built with 15 significant digits,
see benchmarks/kaks/validate_kaks.py) and golden/sim_small.yn00_w<0|1>.tsv (PAML yn00 4.10.10, rst line at 12
digits). COAL's own outputs for regression: golden/<set>.coal.kaks and golden/<set>.coal.dnds.tsv.
"""
from __future__ import annotations

import argparse
import random
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
METHODS = ["NG", "LWL", "MLWL", "LPB", "MLPB", "YN", "MYN"]
BASES = "TCAG"
CODE1 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
SENSE = [BASES[i // 16] + BASES[(i // 4) % 4] + BASES[i % 4] for i in range(64) if CODE1[i] != "*"]
AA = {BASES[i // 16] + BASES[(i // 4) % 4] + BASES[i % 4]: CODE1[i] for i in range(64)}


def rand_cds(rng: random.Random, n: int) -> str:
    return "ATG" + "".join(rng.choice(SENSE) for _ in range(n - 1))


def mutate(rng: random.Random, s: str, frac: float, synonymous_only: bool = False) -> str:
    cod = [s[i:i + 3] for i in range(0, len(s), 3)]
    for k in range(1, len(cod)):
        if rng.random() >= frac:
            continue
        for _ in range(20):
            pos = rng.randrange(3)
            c = cod[k][:pos] + rng.choice([b for b in "ACGT" if b != cod[k][pos]]) + cod[k][pos + 1:]
            if AA[c] == "*":
                continue
            if synonymous_only and AA[c] != AA[cod[k]]:
                continue
            cod[k] = c
            break
    return "".join(cod)


def edge_cases() -> list[tuple[str, str, str]]:
    rng = random.Random(20261111)
    base = rand_cds(rng, 120)
    out = []
    out.append(("identical", base, base))
    out.append(("all_synonymous", base, mutate(rng, base, 0.3, synonymous_only=True)))
    out.append(("typical", base, mutate(rng, base, 0.15)))
    out.append(("saturated_unrelated", base, rand_cds(random.Random(7), 120)))
    m = mutate(rng, base, 0.15)
    g1 = base[:30] + "---" * 5 + base[45:]
    g2 = m[:60] + "A-G" + m[63:90] + "------" + m[96:]
    out.append(("gaps", g1, g2))
    a2 = list(mutate(rng, base, 0.15))
    for p in (10, 50, 51, 200):
        a2[p] = "N"
    a2[120] = "R"
    out.append(("ambiguous", base, "".join(a2)))
    m = mutate(rng, base, 0.15)
    out.append(("internal_stop", base, m[:90] + "TAA" + m[93:150] + "TGA" + m[153:]))
    out.append(("terminal_stop", base + "TAA", mutate(rng, base, 0.15) + "TGA"))
    out.append(("terminal_stop_one_side", base + "TAG", mutate(rng, base, 0.15) + "---"))
    out.append(("lowercase_rna", base.lower().replace("t", "u"), mutate(rng, base, 0.15).lower()))
    out.append(("length_not_multiple_of_3", base[:299], mutate(rng, base, 0.1)[:299]))
    out.append(("all_gap", base, "-" * len(base)))
    mw = "".join(rng.choice(["ATG", "TGG"]) for _ in range(40))
    out.append(("met_trp_only", mw, "".join(rng.choice(["ATG", "TGG"]) for _ in range(40))))
    # codon pairs differing at 2 or 3 positions, some pathways through stop codons
    multi1 = "TGGTACCAATTACGAAGAAGGATAATGTTT" * 3
    multi2 = "TAATAGTGACTAAGGCGGAGAATGATATTA" * 3
    a, b = [], []
    for i in range(0, len(multi1), 3):
        c1, c2 = multi1[i:i + 3], multi2[i:i + 3]
        if AA[c1] != "*" and AA[c2] != "*":
            a.append(c1)
            b.append(c2)
    out.append(("multi_hit_codons", base + "".join(a), base + "".join(b)))
    out.append(("short_one_codon", "ATGAAA", "ATGAAG"))
    return out


def write_axt(path: Path, pairs) -> None:
    with open(path, "w", newline="\n") as f:
        for name, s1, s2 in pairs:
            f.write(f"{name}\n{s1}\n{s2}\n\n")


def sim_small() -> list[tuple[str, str, str]]:
    sys.path.insert(0, str(ROOT / "benchmarks" / "kaks"))
    import numpy as np
    from simulate_pairs import simulate
    rng = np.random.default_rng(20261111)
    out = []
    for w in (0.1, 0.5, 1.0, 2.0):
        for t in (0.3, 1.0, 2.0):
            for r in range(4):
                s1, s2 = simulate(rng, 300, 2.0, w, t)
                out.append((f"sim_w{w:g}_t{t:g}_k2_{r}&b", s1, s2))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--coal-cpp")
    ap.add_argument("--kaks")
    ap.add_argument("--yn00")
    a = ap.parse_args()
    write_axt(HERE / "edge_cases.axt", edge_cases())
    if not (HERE / "sim_small.axt").exists():
        write_axt(HERE / "sim_small.axt", sim_small())
    gold = HERE / "golden"
    gold.mkdir(exist_ok=True)
    mflags = sum((["-m", m] for m in METHODS), [])
    for name in ("edge_cases", "sim_small"):
        axt = HERE / f"{name}.axt"
        if a.kaks:
            for m in METHODS:
                with tempfile.TemporaryDirectory() as td:
                    subprocess.run([a.kaks, "-i", str(axt), "-o", str(Path(td) / "o"), "-m", m],
                                   capture_output=True, check=True)
                    (gold / f"{name}.KaKs_Calculator.{m}.kaks").write_text((Path(td) / "o").read_text())
        if a.coal_cpp:
            subprocess.run([a.coal_cpp, "kaks", "-i", str(axt), "-o", str(gold / f"{name}.coal.kaks"), "-q",
                            *mflags], check=True)
            subprocess.run([a.coal_cpp, "dnds", str(axt), "-o", str(gold / f"{name}.coal.dnds.tsv"), "-q",
                            *mflags], check=True)
    if a.yn00:
        sys.path.insert(0, str(ROOT / "benchmarks" / "kaks"))
        from validate_kaks import read_axt, yn00_one
        for w in (1, 0):
            with open(gold / f"sim_small.yn00_w{w}.tsv", "w", newline="\n") as f:
                f.write("pair\tdN\tdS\tS\tN\tSEdN\tSEdS\n")
                for p in read_axt(HERE / "sim_small.axt"):
                    r = yn00_one(a.yn00, p[1], p[2], w)
                    f.write(f"{p[0]}\t{r['dN']!r}\t{r['dS']!r}\t{r['S']!r}\t{r['N']!r}\t{r['SEdN']!r}\t{r['SEdS']!r}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
