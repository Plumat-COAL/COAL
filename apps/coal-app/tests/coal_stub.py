"""Test stand-in for `python -m coal` (NOT COAL): answers `--cpp-path`, `dnds` and `pipeline pairs` the way the app
reads them (the same files, columns, progress lines and exit statuses, made-up numbers), so the app's Ka/Ks and
pipeline features are tested on machines without coal-cpp or DIAMOND (the CI's app jobs install the pure-Python
engine). Everything else is handed to the real `python -m coal`.

Environment: COAL_STUB_NO_CPP=1 (no coal-cpp), COAL_STUB_SLEEP=seconds (between pipeline stages, to test Stop),
COAL_STUB_PIDFILE=path (the pipeline writes its PID there), COAL_STUB_FAIL=1 (the pipeline fails: exit 1).
"""
import json
import os
import signal
import subprocess
import sys
import time

DNDS_COLUMNS = ("seq1", "seq2", "method", "dN", "dS", "dNdS", "SE_dN", "SE_dS", "N_sites", "S_sites", "Nd", "Sd",
                "codons_compared", "status", "note")
PIPE_COLUMNS = DNDS_COLUMNS[:13] + ("prot_aln_length", "prot_identity", "status", "note")


def read_fasta(path):
    recs, name = [], None
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                recs.append([line[1:].split()[0], ""])
            elif line and recs:
                recs[-1][1] += line
    return recs


def read_axt(path):
    recs, block = [], []
    with open(path, encoding="utf-8") as fh:
        for line in list(fh) + [""]:
            s = line.strip()
            if not s:
                if block:
                    name = block[0]
                    a, _, b = name.partition("&")
                    recs.append(((a, b or "NA"), block[1] if len(block) > 1 else "", block[2] if len(block) > 2
                                 else ""))
                block = []
            elif not (s.startswith("#") and not block):
                block.append(s)
    return recs


def opt(args, flag, default=None):
    return args[args.index(flag) + 1] if flag in args else default


def opts(args, flag):
    return [args[i + 1] for i, a in enumerate(args) if a == flag]


def row(a, b, sa, sb, method, extra=()):
    if len(sa) != len(sb) or len(sa) % 3:
        vals = dict(seq1=a, seq2=b, method=method, status="error",
                    note=f"aligned sequences differ in length ({len(sa)} and {len(sb)})")
    else:
        diff = sum(x != y for x, y in zip(sa, sb))
        dn = round(0.001 * (1 + diff % 7), 6)
        ds = round(0.1 * (1 + diff % 5), 6)
        vals = dict(seq1=a, seq2=b, method=method, dN=dn, dS=ds, dNdS=round(dn / ds, 6), codons_compared=len(sa) // 3,
                    status="ok", note="")
    vals.update(extra)
    return vals


def write_tsv(path, cols, rows):
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "NA")) for c in cols) + "\n")


def dnds(args):
    src, out = args[0], opt(args, "-o")
    methods = opts(args, "-m") or ["YN"]
    first = open(src, encoding="utf-8").read().lstrip()[:1]
    if first == ">":
        recs = read_fasta(src)
        if opt(args, "--pairs", "all") == "consecutive":
            pairs = [(recs[i], recs[i + 1]) for i in range(0, len(recs) - 1, 2)]
        else:
            pairs = [(recs[i], recs[j]) for i in range(len(recs)) for j in range(i + 1, len(recs))]
        pairs = [((x[0], y[0]), x[1], y[1]) for x, y in pairs]
    else:
        pairs = read_axt(src)
    rows = [row(n[0], n[1], sa, sb, m) for n, sa, sb in pairs for m in methods]
    write_tsv(out, DNDS_COLUMNS, rows)
    ok = sum(r["status"] == "ok" for r in rows)
    print(f"coal dnds: {len(pairs)} pair(s), {len(methods)} method(s): {ok} ok, 0 warning, 0 not_converged, "
          f"0 saturated, 0 no_data, {len(rows) - ok} error", file=sys.stderr)
    if len(rows) - ok:
        print(f"coal dnds: warning: {len(rows) - ok} result(s) could not be computed (status error)", file=sys.stderr)
    if not ok:
        print("coal dnds: no pair has an estimate: the two sequences of the pairs differ in length", file=sys.stderr)
        return 2
    return 0


def say(text):
    print(f"coal pipeline: {text}", file=sys.stderr, flush=True)


def pipeline(args):
    out = opt(args, "-o")
    os.makedirs(os.path.join(out, "work"), exist_ok=True)
    with open(os.path.join(out, "work", "state.json"), "w") as fh:
        fh.write("{}")
    if os.environ.get("COAL_STUB_PIDFILE"):
        with open(os.environ["COAL_STUB_PIDFILE"], "w") as fh:
            fh.write(str(os.getpid()))
    report = {"schema": "coal.pipeline/1", "mode": "pairs", "status": "running", "exit_code": None,
              "command": ["pipeline"] + args, "counts": {}, "failures_total": 0, "failures": [], "outputs": {},
              "cite": ["coal", "suyama2006", "yang2000", "wang2010", "yang2007", "buchfink2021"]}

    def write_report():
        with open(os.path.join(out, "report.json"), "w", encoding="utf-8") as fh:
            json.dump(report, fh)

    def on_term(*_):
        report.update(status="interrupted", exit_code=143, error="interrupted by a signal")
        write_report()
        say("interrupted (the report says so; nothing is left running)")
        os._exit(143)

    signal.signal(signal.SIGTERM, on_term)
    pause = float(os.environ.get("COAL_STUB_SLEEP") or 0)
    q, s = read_fasta(opt(args, "--query")), read_fasta(opt(args, "--subject"))
    say(f"translated {len(q)} query and {len(s)} subject records (0 not usable)")
    time.sleep(pause)
    if os.environ.get("COAL_STUB_FAIL"):
        report.update(status="failed", exit_code=1, error="DIAMOND failed: stub failure")
        write_report()
        say("error: DIAMOND failed: stub failure")
        return 1
    if opt(args, "--pairs"):
        with open(opt(args, "--pairs"), encoding="utf-8") as fh:
            want = [line.split("\t")[:2] for line in fh.read().splitlines() if line.strip()]
        qd, sd = dict(q), dict(s)
        pairs = [((a, b.strip()), qd.get(a, ""), sd.get(b.strip(), "")) for a, b in want]
    else:
        pairs = [((x[0], y[0]), x[1], y[1]) for x, y in zip(q, s)]
    say(f"{len(pairs)} pairs ({'given' if opt(args, '--pairs') else opt(args, '--ortho', 'rbh')})")
    time.sleep(pause)
    methods = opts(args, "-m") or ["YN"]
    rows = [row(n[0], n[1], a, b, m, {"prot_aln_length": len(a) // 3, "prot_identity": 90.0})
            for n, a, b in pairs for m in methods]
    write_tsv(os.path.join(out, "dnds.tsv"), PIPE_COLUMNS, rows)
    write_tsv(os.path.join(out, "failed.tsv"), ("level", "cluster_id", "seq1", "seq2", "code", "message"), [])
    write_tsv(os.path.join(out, "genes.tsv"), ("set", "id", "status"),
              [{"set": "query", "id": n, "status": "ok"} for n, _ in q])
    counts = {k: sum(r["status"] == k for r in rows) for k in ("ok", "warning", "not_converged", "saturated",
                                                                "no_data", "error")}
    report.update(status="completed", exit_code=0,
                  counts={"genes": {"total": len(q) + len(s), "ok": len(q) + len(s), "warning": 0, "failed": 0},
                          "orthologs": len(pairs), "pairs": dict(total=len(pairs), **counts)},
                  outputs={"dnds": "dnds.tsv", "failed": "failed.tsv", "genes": "genes.tsv"})
    write_report()
    say(f"summary: {len(q) + len(s)} records (0 not usable); {len(pairs)} pairs: {counts['ok']} ok")
    say(f"summary: results in {os.path.join(out, 'dnds.tsv')}")
    return 0


def main(args):
    if args == ["--cpp-path"]:
        if os.environ.get("COAL_STUB_NO_CPP"):
            print("coal: coal-cpp not found (stub)", file=sys.stderr)
            return 1
        print(os.path.join(os.sep, "stub", "coal-cpp"))
        return 0
    if args[:1] == ["dnds"]:
        return dnds(args[1:])
    if args[:2] == ["pipeline", "pairs"]:
        return pipeline(args[2:])
    return subprocess.call([sys.executable, "-m", "coal", *args])


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
