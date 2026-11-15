"""Validate COAL's method GY (maximum likelihood, Goldman & Yang 1994) against PAML's codeml and simulated truth.

  python3 benchmarks/kaks/validate_gy.py sim  --coal build/coal-cpp --codeml codeml --out work/gy
  python3 benchmarks/kaks/validate_gy.py real --coal build/coal-cpp --codeml codeml --out work/gy \
        --pairs name=in.axt[:N] [--pairs name=codon.fasta[:N]] ...

sim: pairs simulated under GY94 / F3x4 (sim_truth/codonsim.py; truth is the model's dN and dS as Yang & Nielsen
2000 define them): dS 0.1 ... 3, omega 0.2 / 1, kappa 2 / 4, equal and GC3-skewed codon usage, 300 codons.
The first --codeml-reps replicates of each cell go to codeml as well (runmode = -2, CodonFreq = 2, one data set
per pair); all --reps replicates are compared with the truth, for GY and YN.
real: pairs of a real data set (AXT, or an aligned FASTA of consecutive pairs). codeml gets the columns COAL
compares (no gap, no ambiguous base, no stop codon), COAL gets the rows as they are.

Reports (stdout and <out>/*.tsv): the difference in lnL (COAL - codeml; codeml prints 6 decimals), the relative
differences of dN, dS, omega, kappa and t, and the relative error against the truth by dS.
"""
import argparse
import os
import re
import statistics
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
STOPS = {"TAA", "TAG", "TGA"}


def read_axt(path, limit=None):
    out = []
    with open(path) as fh:
        block = []
        for line in fh:
            line = line.strip()
            if not line:
                if len(block) >= 3:
                    out.append((block[0], block[1], block[2]))
                block = []
                if limit and len(out) >= limit:
                    break
                continue
            block.append(line)
        if len(block) >= 3 and not (limit and len(out) >= limit):
            out.append((block[0], block[1], block[2]))
    return out


def read_fasta_pairs(path, limit=None):
    recs, name, seq = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(seq)))
                    if limit and len(recs) >= 2 * limit:
                        name = None
                        break
                name, seq = line[1:].split()[0], []
            elif line:
                seq.append(line)
        if name is not None:
            recs.append((name, "".join(seq)))
    return [(f"{recs[i][0]}&{recs[i + 1][0]}", recs[i][1], recs[i + 1][1]) for i in range(0, len(recs) - 1, 2)]


def clean_pair(a, b):
    """The codon columns COAL compares under the standard code: both codons ACGT only, neither a stop codon."""
    a, b = a.upper().replace("U", "T"), b.upper().replace("U", "T")
    oa, ob = [], []
    for i in range(0, min(len(a), len(b)) - 2, 3):
        x, y = a[i:i + 3], b[i:i + 3]
        if set(x) <= set("ACGT") and set(y) <= set("ACGT") and x not in STOPS and y not in STOPS:
            oa.append(x)
            ob.append(y)
    return "".join(oa), "".join(ob)


def run_coal(coal, pairs, workdir, methods=("GY",), table=1, jobs=1):
    path = os.path.join(workdir, "coal_in.axt")
    with open(path, "w") as fh:
        for i, (_, a, b) in enumerate(pairs):
            fh.write(f"p{i}&x\n{a}\n{b}\n\n")
    cmd = [coal, "dnds", path, "--ml-columns", "--digits", "12", "-q", "-j", str(jobs), "-c", str(table)]
    for m in methods:
        cmd += ["-m", m]
    t0 = time.time()
    txt = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
    dt = time.time() - t0
    lines = txt.splitlines()
    head = lines[0].split("\t")
    res = {m: [None] * len(pairs) for m in methods}
    for ln in lines[1:]:
        r = dict(zip(head, ln.split("\t")))
        res[r["method"]][int(r["seq1"][1:])] = r
    return res, dt


CTL = ("seqfile = seq.phy\noutfile = mlc\nnoisy = 0\nverbose = 0\nrunmode = -2\nseqtype = 1\nCodonFreq = 2\n"
       "model = 0\nNSsites = 0\nicode = 0\nfix_kappa = 0\nkappa = 2\nfix_omega = 0\nomega = 0.4\ncleandata = 0\n"
       "ndata = {n}\n")
NUM = r"(-?[0-9.]+(?:e[-+]?[0-9]+)?|-?nan|-?inf)"
PAT_LNL = re.compile(r"lnL\s*=\s*" + NUM + r"\s*\n\s*" + NUM + r"\s+" + NUM + r"\s+" + NUM)
PAT_D = re.compile(r"t=\s*" + NUM + r"\s+S=\s*" + NUM + r"\s+N=\s*" + NUM + r"\s+dN/dS=\s*" + NUM + r"\s+dN\s*=\s*" +
                   NUM + r"\s+dS\s*=\s*" + NUM)


def run_codeml(codeml, pairs, workdir):
    """One data set per pair. Returns a list of dicts (None where codeml printed nothing) and the seconds."""
    with tempfile.TemporaryDirectory(dir=workdir) as td:
        with open(f"{td}/seq.phy", "w") as fh:
            for _, a, b in pairs:
                fh.write(f"2 {len(a)}\nseqA  {a}\nseqB  {b}\n\n")
        with open(f"{td}/codeml.ctl", "w") as fh:
            fh.write(CTL.format(n=len(pairs)))
        t0 = time.time()
        subprocess.run(["nice", "-n", "10", codeml, "codeml.ctl"], cwd=td, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, check=False, stdin=subprocess.DEVNULL)
        dt = time.time() - t0
        txt = open(f"{td}/mlc").read() if os.path.exists(f"{td}/mlc") else ""
    # split by data set
    blocks = re.split(r"(?:^|\n)\s*Data ?set \d+\s*\n", txt)
    blocks = blocks[1:] if len(blocks) > 1 else [txt]
    out = []
    for b in blocks:
        m1, m2 = PAT_LNL.search(b), PAT_D.search(b)
        if not (m1 and m2):
            out.append(None)
            continue
        out.append(dict(lnL=float(m1.group(1)), t=float(m1.group(2)), kappa=float(m1.group(3)),
                        omega=float(m1.group(4)), S=float(m2.group(2)), N=float(m2.group(3)),
                        dN=float(m2.group(5)), dS=float(m2.group(6))))
    out += [None] * (len(pairs) - len(out))
    return out[:len(pairs)], dt


def fl(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-12)


def compare(tag, pairs, coal_rows, cm, out):
    """COAL GY rows against codeml's. Estimates are compared where the two found the same maximum (lnL within
    1e-3) and codeml is not at one of its bounds (t 50, kappa 0.4 / 999, omega 0.001 / 99): at a bound the two
    optimise over different regions, and where COAL's likelihood is higher codeml stopped at a lower maximum."""
    n = both = worse = bound = other = 0
    dl, rels = [], {k: [] for k in ("dN", "dS", "omega", "kappa", "t")}
    with open(os.path.join(out, f"{tag}_vs_codeml.tsv"), "w") as fh:
        fh.write("pair\tstatus\tlnL_coal\tlnL_codeml\tdiff\tdN_coal\tdN_codeml\tdS_coal\tdS_codeml\tw_coal\tw_codeml"
                 "\tkappa_coal\tkappa_codeml\tt_coal\tt_codeml\tcodeml_at_bound\n")
        for (name, _, _), r, c in zip(pairs, coal_rows, cm):
            n += 1
            if c is None or r is None:
                continue
            lc = fl(r["lnL"])
            atb = (c["t"] >= 49.99 or c["kappa"] <= 0.4001 or c["kappa"] >= 998.9 or c["omega"] <= 0.00101 or
                   c["omega"] >= 98.99)
            fh.write("\t".join(str(x) for x in (name, r["status"], r["lnL"], c["lnL"], lc - c["lnL"], r["dN"],
                                                c["dN"], r["dS"], c["dS"], r["dNdS"], c["omega"], r["kappa"],
                                                c["kappa"], r["t"], c["t"], int(atb))) + "\n")
            if lc != lc:
                continue
            both += 1
            dl.append(lc - c["lnL"])
            if lc < c["lnL"] - 1e-4:
                worse += 1
                print(f"  WORSE {tag} {name}: COAL {lc} codeml {c['lnL']}")
            if atb:
                bound += 1
                continue
            if abs(lc - c["lnL"]) > 1e-3:
                other += 1
                continue
            if r["status"] not in ("ok", "warning"):
                continue
            # codeml prints 4 or 5 decimals: compare where its number has at least 3 significant digits
            for k, a, b in (("dN", fl(r["dN"]), c["dN"]), ("dS", fl(r["dS"]), c["dS"]),
                            ("omega", fl(r["dNdS"]), c["omega"]), ("kappa", fl(r["kappa"]), c["kappa"]),
                            ("t", fl(r["t"]), c["t"])):
                if a == a and b >= 0.01:
                    rels[k].append(rel(a, b))
    print(f"{tag}: {n} pairs, {both} with lnL from both; COAL - codeml lnL: min {min(dl):.3g} median "
          f"{statistics.median(dl):.3g} max {max(dl):.3g}; COAL worse by > 1e-4: {worse}; COAL better by > 1e-4: "
          f"{sum(d > 1e-4 for d in dl)}; codeml at one of its bounds: {bound}; another maximum (|diff| > 1e-3, "
          f"codeml not at a bound): {other}")
    print("  where both found the same maximum (codeml prints 4 or 5 decimals):")
    for k, v in rels.items():
        if v:
            p95 = sorted(v)[int(0.95 * (len(v) - 1))]
            print(f"  {k}: n {len(v)} median rel diff {statistics.median(v):.2e}  95% {p95:.2e}  max {max(v):.2e}")
    return worse


def cmd_sim(a):
    import numpy as np
    sys.path.insert(0, os.path.join(HERE, "sim_truth"))
    from codonsim import FREQS, CodonModel, sample_pair_P, to_str
    os.makedirs(a.out, exist_ok=True)
    DS = [0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
    pairs, meta = [], []
    for fi, f in enumerate(("equal", "skewed")):
        for k in (2.0, 4.0):
            for w in (0.2, 1.0):
                m = CodonModel(k, w, FREQS[f])
                rng = np.random.default_rng([20261115, fi, int(k), int(w * 10), a.codons])
                for d in DS:
                    t = m.t_for_dS(d)
                    for rep in range(a.reps):
                        x, y = sample_pair_P(m, t, a.codons, rng)
                        pairs.append((f"{f}_k{k:g}_w{w:g}_d{d:g}_r{rep}", to_str(x), to_str(y)))
                        meta.append(dict(f=f, k=k, w=w, dS=d, dN=m.dN_of_t(t), rep=rep))
    print(f"sim: {len(pairs)} pairs of {a.codons} codons")
    res, dt = run_coal(a.coal, pairs, a.out, methods=("GY", "YN"), jobs=a.jobs)
    _, dt_gy = run_coal(a.coal, pairs, a.out, methods=("GY",), jobs=1)
    _, dt_yn = run_coal(a.coal, pairs, a.out, methods=("YN",), jobs=1)
    print(f"speed (1 thread, {a.codons} codons): GY {1000 * dt_gy / len(pairs):.1f} s per 1,000 pairs, "
          f"YN {1000 * dt_yn / len(pairs):.2f} s per 1,000 pairs")
    # truth
    print("accuracy against the truth: median relative error of dS / of dN/dS [fraction of pairs with status ok]")
    with open(os.path.join(a.out, "sim_truth.tsv"), "w") as fh:
        fh.write("method\tusage\tdS\tn\tok\tmedian_err_dS\tmedian_err_dNdS\tmedian_err_dN\n")
        for usage in ("equal", "skewed", "all"):
            for meth in ("GY", "YN"):
                line = []
                for d in DS:
                    eS, eW, eN, nall = [], [], [], 0
                    for mt, r in zip(meta, res[meth]):
                        if mt["dS"] != d or (usage != "all" and mt["f"] != usage):
                            continue
                        nall += 1
                        if r["status"] != "ok":
                            continue
                        eS.append(fl(r["dS"]) / d - 1)
                        eN.append(fl(r["dN"]) / mt["dN"] - 1)
                        if fl(r["dNdS"]) == fl(r["dNdS"]):
                            eW.append(fl(r["dNdS"]) / mt["w"] - 1)
                    ms = statistics.median(eS) if eS else float("nan")
                    mw = statistics.median(eW) if eW else float("nan")
                    mn = statistics.median(eN) if eN else float("nan")
                    fh.write(f"{meth}\t{usage}\t{d:g}\t{nall}\t{len(eS)}\t{ms:.4f}\t{mw:.4f}\t{mn:.4f}\n")
                    line.append(f"{100 * ms:+.0f}% / {100 * mw:+.0f}% [{len(eS) / nall:.2f}]")
                print(f"  {meth:3s} {usage:7s} " + " | ".join(line))
    print("  (columns: dS " + ", ".join(f"{d:g}" for d in DS) + ")")
    st = {}
    for r in res["GY"]:
        st[r["status"]] = st.get(r["status"], 0) + 1
    print("  GY statuses:", st)
    # codeml on the first replicates
    sub = [i for i, mt in enumerate(meta) if mt["rep"] < a.codeml_reps]
    sp = [pairs[i] for i in sub]
    cm, dtc = run_codeml(a.codeml, sp, a.out)
    print(f"codeml: {len(sp)} pairs in {dtc:.0f} s ({1000 * dtc / len(sp):.0f} s per 1,000 pairs)")
    return compare("sim", sp, [res["GY"][i] for i in sub], cm, a.out)


def cmd_real(a):
    os.makedirs(a.out, exist_ok=True)
    worse = 0
    for spec in a.pairs:
        name, path = spec.split("=", 1)
        limit = None
        if ":" in path:
            path, lim = path.rsplit(":", 1)
            limit = int(lim)
        raw = read_fasta_pairs(path, limit) if open(path).read(1) == ">" else read_axt(path, limit)
        pairs, cleaned = [], []
        for nm, x, y in raw:
            cx, cy = clean_pair(x, y)
            if len(cx) >= 3 * a.min_codons and cx != cy:
                pairs.append((nm, x, y))
                cleaned.append((nm, cx, cy))
        res, dt = run_coal(a.coal, pairs, a.out, methods=("GY", "YN"), jobs=a.jobs)
        _, dt_gy = run_coal(a.coal, pairs, a.out, methods=("GY",), jobs=1)
        _, dt_yn = run_coal(a.coal, pairs, a.out, methods=("YN",), jobs=1)
        cod = statistics.median(len(c[1]) // 3 for c in cleaned)
        print(f"{name}: {len(pairs)} pairs, median {cod} codons; speed (1 thread): GY "
              f"{1000 * dt_gy / len(pairs):.1f} s per 1,000 pairs, YN {1000 * dt_yn / len(pairs):.2f}")
        st = {}
        for r in res["GY"]:
            st[r["status"]] = st.get(r["status"], 0) + 1
        print("  GY statuses:", st)
        cm, dtc = run_codeml(a.codeml, cleaned, a.out)
        print(f"  codeml: {dtc:.0f} s ({1000 * dtc / len(pairs):.0f} s per 1,000 pairs)")
        worse += compare(name, pairs, res["GY"], cm, a.out)
    return worse


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for nm in ("sim", "real"):
        p = sub.add_parser(nm)
        p.add_argument("--coal", required=True)
        p.add_argument("--codeml", required=True)
        p.add_argument("--out", required=True)
        p.add_argument("--jobs", type=int, default=1)
    s = sub.choices["sim"]
    s.add_argument("--reps", type=int, default=50)
    s.add_argument("--codeml-reps", type=int, default=6)
    s.add_argument("--codons", type=int, default=300)
    r = sub.choices["real"]
    r.add_argument("--pairs", action="append", required=True)
    r.add_argument("--min-codons", type=int, default=30)
    a = ap.parse_args()
    worse = cmd_sim(a) if a.cmd == "sim" else cmd_real(a)
    sys.exit(1 if worse else 0)


if __name__ == "__main__":
    main()
