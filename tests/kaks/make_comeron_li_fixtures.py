"""Make the test input and the golden files for the methods Comeron and Li (tests/test_kaks_comeron_li.py).

  python tests/kaks/make_comeron_li_fixtures.py                 # the input only: comeron_li_cases.axt
  python tests/kaks/make_comeron_li_fixtures.py --coal-cpp build/coal-cpp \
      --gestimator-ref work/gestimator_ref --rscript Rscript      # + reference and COAL golden files

The sequences are made here (by hand, and by seeded generators), so the data are COAL's own. For each of the
sets comeron_li_cases, edge_cases and sim_small:
  golden/<set>.gestimator_ref.tsv   orthologr 0.4.3's gestimator code, 17 digits (benchmarks/kaks/ref/
                                    gestimator_ref.cpp; build line in that file)
  golden/<set>.seqinr_ref.tsv       seqinr::kaks, 17 digits (benchmarks/kaks/ref/seqinr_ref.R; the file's last
                                    column names are in that script)
  golden/<set>.coal.comeron_li.tsv  `coal dnds -m Comeron -m Li`, the default flavours (regression)
  golden/<set>.coal.comeron_li.kaks `coal kaks -m Comeron -m Li` (regression of the 22-column format)
"""
from __future__ import annotations

import argparse
import importlib.util
import random
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
BASES = "TCAG"
CODE1 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
SENSE = [BASES[i // 16] + BASES[(i // 4) % 4] + BASES[i % 4] for i in range(64) if CODE1[i] != "*"]


def _bench(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "benchmarks" / "kaks" / (name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def mutate(rng, cds, p):
    out = []
    for i in range(0, len(cds), 3):
        c = cds[i:i + 3]
        if rng.random() < p:
            c = rng.choice(SENSE) if rng.random() < 0.3 else c[:2] + rng.choice(BASES)
            if CODE1[BASES.index(c[0]) * 16 + BASES.index(c[1]) * 4 + BASES.index(c[2])] == "*":
                c = "GCT"
        out.append(c)
    return "".join(out)


def cases():
    rng = random.Random(19950601)
    g = "".join(rng.choice(SENSE) for _ in range(40))   # 120 nucleotides
    h = mutate(rng, g, 0.25)
    far = "".join(rng.choice(SENSE) for _ in range(40))
    out = [
        ("identical", g, g),
        ("typical", g, h),
        ("unrelated_saturated", g, far),
        ("one_codon_identical", "ATG", "ATG"),
        ("one_codon_synonymous", "GCT", "GCA"),
        ("one_codon_nonsynonymous", "GCT", "GAT"),
        ("one_codon_two_changes", "GCT", "AAT"),
        ("one_codon_three_changes", "GCT", "AAG"),
        ("two_codons", "ATGGCT", "ATGGCC"),
        ("three_codons", "ATGGCTAAA", "ATGGCAAAG"),
        ("three_codons_all_differ", "TTTCCCAAA", "GGGAAATTT"),
        ("gap_codons", g[:60] + "---" + g[63:120], h[:30] + "------" + h[36:120]),
        ("gap_inside_codon", g[:61] + "-" + g[62:120], h),
        ("dot_gaps", g[:60] + "..." + g[63:120], h),
        ("n_in_codon", g[:61] + "N" + g[62:120], h),
        ("iupac_letters", g[:30] + "RAY" + g[33:60] + "ACN" + g[63:120], h[:90] + "KWS" + h[93:120]),
        ("internal_stop_one_sequence", g[:30] + "TAA" + g[33:120], h),
        ("internal_stop_both", g[:30] + "TAA" + g[33:120], h[:30] + "TAA" + h[33:120]),
        ("internal_stops_differ", g[:30] + "TAA" + g[33:120], h[:30] + "TGA" + h[33:120]),
        ("three_internal_stops", g[:30] + "TAG" + g[33:60] + "TGA" + g[63:120], h[:90] + "TAA" + h[93:120]),
        ("terminal_stop_both", g + "TAA", h + "TAA"),
        ("terminal_stops_differ", g + "TAA", h + "TGA"),
        ("terminal_stop_one_side", g + "TAG", h + "GCT"),
        ("terminal_stop_then_gap", g + "TAA---", h + "GCTTGA"),
        ("lower_case", g.lower(), h.lower()),
        ("mixed_case", g[:60].lower() + g[60:], h),
        ("rna", g.replace("T", "U"), h),
        ("rna_lower_both", g.replace("T", "U").lower(), h.replace("T", "U").lower()),
        ("all_gap", "---------", "---------"),
        ("one_sequence_all_gap", "ATGGCTAAA", "---------"),
        ("only_stop_codons", "TAATGA", "TAGTAA"),
        ("stop_against_sense_only", "TAA", "TAC"),
        ("met_trp_only", "ATGTGGATGTGG", "ATGTGGTGGATG"),
        ("fourfold_only", "GCTGCTGCTGCT", "GCAGCCGCGGCT"),
        ("transitions_only", "".join(c for c in ["GCT", "AAA", "GAT", "TTT", "CAT"] * 6),
         "".join(c for c in ["GCC", "AAG", "GAC", "TTC", "CAC"] * 6)),
        ("transversions_only", "".join(c for c in ["GCT", "ACA", "GGT", "CCT", "GTA"] * 6),
         "".join(c for c in ["GCA", "ACC", "GGA", "CCG", "GTT"] * 6)),
        ("synonymous_saturated", "".join(["GCT", "GGT", "CCT", "ACT", "GTT"] * 6),
         "".join(["GCA", "GGC", "CCG", "ACA", "GTG"] * 6)),
        ("nonsynonymous_saturated", "AAACCCGGGTTT" * 4, "CCCAAATTTGGG" * 4),
        ("length_not_multiple_of_3", "ATGGCTAA", "ATGGCTAA"),
    ]
    vcl = _bench("validate_comeron_li")
    out += [("stops_" + n, a, b) for n, a, b in vcl.make_pairs_stops(260, 20261113)]
    rsp = _bench("random_short_pairs")
    out += [("multi_" + n, a, b) for n, a, b in rsp.make_pairs_multi(160, 20261112)]
    longer = random.Random(7)
    for k in range(40):   # genes of 150 codons at four divergences, with a terminal stop in most
        a = "".join(longer.choice(SENSE) for _ in range(150))
        b = mutate(longer, a, (0.05, 0.2, 0.5, 0.9)[k % 4])
        stop = longer.choice(["TAA", "TAG", "TGA"])
        if k % 5:
            a, b = a + stop, b + longer.choice(["TAA", "TAG", "TGA"])
        out.append((f"gene{k}", a, b))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--coal-cpp")
    ap.add_argument("--gestimator-ref")
    ap.add_argument("--rscript")
    a = ap.parse_args()
    axt = HERE / "comeron_li_cases.axt"
    with open(axt, "w", newline="\n") as f:
        for name, s1, s2 in cases():
            f.write(f"{name}\n{s1}\n{s2}\n\n")
    gold = HERE / "golden"
    for s in ("comeron_li_cases", "edge_cases", "sim_small"):
        src = HERE / (s + ".axt")
        if a.gestimator_ref:
            with open(gold / f"{s}.gestimator_ref.tsv", "w", newline="\n") as f:
                subprocess.run([a.gestimator_ref, str(src)], check=True, stdout=f)
        if a.rscript:
            subprocess.run([a.rscript, str(ROOT / "benchmarks" / "kaks" / "ref" / "seqinr_ref.R"), str(src),
                            str(gold / f"{s}.seqinr_ref.tsv")], check=True, stdout=subprocess.DEVNULL)
        if a.coal_cpp:
            subprocess.run([a.coal_cpp, "dnds", str(src), "--format", "axt", "-m", "Comeron", "-m", "Li", "-q",
                            "--digits", "15", "-o", str(gold / f"{s}.coal.comeron_li.tsv")], check=True)
            subprocess.run([a.coal_cpp, "kaks", "-i", str(src), "-m", "Comeron", "-m", "Li", "-q",
                            "-o", str(gold / f"{s}.coal.comeron_li.kaks")], check=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
