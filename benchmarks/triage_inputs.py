"""Tier C triage of cases every tool rejects: inspects the inputs (duplicate IDs, protein vs CDS translation) and,
without --quiet, reruns coal-cpp/coal-py/p2n16 untimed to show their errors. Writes the INPUT_INVALID list.
usage: triage_inputs.py CASE_LIST_FILE [OUT.csv] [--quiet]   (case paths relative to the data dir)"""
import json, subprocess, sys, collections
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import tiers, grade

DATA = tiers.DATA
B = "TCAG"
CODONS = [a + b + c for a in B for b in B for c in B]
STD = dict(zip(CODONS, "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))


def table(n):
    t = dict(STD)
    if n == 4:
        t["TGA"] = "W"
    return t


def diagnose(d: Path):
    m = json.loads((d / "meta.json").read_text())
    tn = int(m["opts"][m["opts"].index("-codontable") + 1]) if "-codontable" in m["opts"] else 1
    pep = grade.read_fasta(d / m["pep"]) or []
    nuc = []
    for f in m["nuc"]:
        nuc += grade.read_fasta(d / f) or []
    out = {"case": m["dir"], "scenario": m["scenario"], "features": m.get("features"), "n_pep": len(pep),
           "n_nuc": len(nuc), "table": tn}
    lens = collections.Counter(len(s) for _, s in pep)
    out["pep_row_lengths"] = dict(lens.most_common(4)) | ({"...distinct": len(lens)} if len(lens) > 4 else {})
    pid = [h.split()[0] for h, _ in pep]
    nid = [h.split()[0] for h, _ in nuc]
    out["ids_match_in_order"] = pid == nid
    out["dup_pep_ids"] = len(pid) - len(set(pid))
    seqs = collections.defaultdict(set)
    for h, q in pep:
        seqs[h.split()[0]].add(q.replace("-", ""))
    dups = [k for k in seqs if pid.count(k) > 1]
    out["dup_ids_identical_seqs"] = sum(len(seqs[k]) == 1 for k in dups)
    out["dup_id_example"] = dups[0] if dups else ""
    probs = []
    t = table(tn)
    for (ph, ps), (nh, ns) in zip(pep, nuc):
        p = ps.replace("-", "").replace(".", "").upper()
        n = ns.replace("-", "").upper()
        tr = "".join(t.get(n[i:i + 3], "X") for i in range(0, len(n) - len(n) % 3, 3))
        trs = tr[:-1] if tr.endswith("*") else tr
        issue = []
        if len(n) % 3:
            issue.append(f"cds_len%3={len(n) % 3}")
        if len(trs) != len(p):
            issue.append(f"len cds/3={len(trs)} (stop stripped) vs prot={len(p)}")
        k = next((i for i, (a, b) in enumerate(zip(trs, p)) if a != b and b not in "XBZJOU" and not (i == 0 and b == "M")), None)
        if k is not None:
            issue.append(f"first aa mismatch at {k + 1}: cds {trs[k:k + 6]} vs prot {p[k:k + 6]}")
        if "*" in trs:
            issue.append(f"internal stop at aa {trs.index('*') + 1}")
        if issue:
            probs.append((ph.split()[0], issue))
    out["bad_seqs"] = len(probs)
    out["bad_examples"] = probs[:2]
    return m, out


def stderr_of(tool, jobs, m, d):
    r = subprocess.run(tiers.tool_cmd(tool, jobs, m, d), capture_output=True, text=True, errors="replace",
                       env=tiers.jobrun.clean_env(), timeout=600)
    lines = [l for l in r.stderr.splitlines() if l.strip()]
    return r.returncode, len(r.stdout), lines


if __name__ == "__main__":
    cases = [l.strip() for l in open(sys.argv[1]) if l.strip()]
    quiet = "--quiet" in sys.argv
    agg = collections.Counter(); rows_out = []
    for c in cases:
        d = DATA / c
        m, o = diagnose(d)
        if o["dup_pep_ids"]:
            cat = "duplicate_ids"
            why = (f"duplicate sequence IDs ({o['dup_pep_ids']} extra rows, e.g. {o['dup_id_example']}; "
                   f"{o['dup_ids_identical_seqs']} of the duplicated IDs carry identical proteins)")
        elif o["bad_seqs"]:
            cat = "cds_protein_mismatch"
            ex = o["bad_examples"][0]
            why = (f"protein is not the translation of its CDS in {o['bad_seqs']} of {o['n_pep']} seqs "
                   f"(features {';'.join(o['features']) or '-'}; e.g. {ex[0]}: {'; '.join(ex[1])})")
        else:
            cat, why = "no_input_problem_found", ""
        agg[cat] += 1
        feats = ";".join(sorted(o["features"])) or "-"
        agg[(cat, feats)] += 1
        rows_out.append({"case": c, "category": cat, "reason": why})
        if quiet:
            continue
        print("=" * 100)
        print(json.dumps(o, default=str))
        for tool, jobs in (("coal-cpp", 1), ("coal-py", 1), ("p2n16", None)):
            rc, nout, lines = stderr_of(tool, jobs, m, d)
            print(f"--- {tool}: exit={rc} stdout_bytes={nout}")
            for l in lines[:5] + (["  ..."] + lines[-2:] if len(lines) > 7 else lines[5:]):
                print("   ", l[:220])
    for k, v in sorted(agg.items(), key=str): print(v, k)
    if len(sys.argv) > 2 and sys.argv[2].endswith(".csv"):
        import csv
        with open(sys.argv[2], "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=["case", "category", "reason"]); w.writeheader(); w.writerows(rows_out)
