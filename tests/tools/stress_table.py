#!/usr/bin/env python3
"""Turn the JSON lines written by stress_run.py / stress_resume.py into the Markdown table of the stress report.

    python tests/tools/stress_table.py results/*.jsonl > table.md

One row per case and engine: shape, input size, engine, wall time, peak memory, exit status and whether the output
was correct (= identical to the simulated truth where the run was a plain fasta run, to PAL2NAL v14 where it ran,
and identical between the two engines). A PAL2NAL row is added where PAL2NAL ran. Resume cases get a short row.
"""
from __future__ import annotations

import json
import sys


def fmt_mb(x):
    return f"{x:,.0f}" if x is not None else ""


def correct(r, res, eng):
    """What the row's 'correct?' column says."""
    parts = []
    if r.get("timeout"):
        return "timeout"
    if "truth_equal" in r:
        parts.append("truth " + ("=" if r["truth_equal"] else "DIFF"))
    p = res.get("pal2nal")
    if p and p.get("out_bytes"):
        parts.append("PAL2NAL " + ("=" if p.get(f"same_as_{eng}") else "DIFF"))
    elif p:
        parts.append("PAL2NAL failed")
    ident = res.get("engines_identical")
    if ident:
        same = all(v in (True, None) for k, v in ident.items() if not isinstance(v, dict))
        parts.append("engines " + ("=" if same else "DIFF"))
    if r.get("stderr", {}).get("traceback"):
        parts.append("TRACEBACK")
    return ", ".join(parts) if parts else ("ran" if r.get("rc") == 0 else "error")


def main(paths):
    rows = []
    for path in paths:
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                if not line.strip():
                    continue
                res = json.loads(line)
                if "resume" in res:
                    rs = res["resume"]
                    n = sum(len(x) for x in rs["rounds"])
                    rows.append((res["name"], f"1,000,000 x 400, {' '.join(rs['args'])}", "", rs["engine"],
                                 f"{rs['ref_seconds']} (ref)", "", f"{n} runs", "OK: byte-identical, no leftovers"
                                 if rs["ok"] else f"PROBLEMS {len(rs['problems'])}"))
                    continue
                for eng, r in res["engines"].items():
                    rc = "timeout" if r.get("timeout") else str(r.get("rc"))
                    rows.append((res["name"], res.get("shape") or "", f"{res['input_mb']:,.1f}", eng,
                                 f"{r['wall_s']:.2f}", fmt_mb(r.get("peak_mb")), rc, correct(r, res, eng)))
                p = res.get("pal2nal")
                if p:
                    rc = "timeout" if p.get("timeout") else str(p.get("rc"))
                    rows.append((res["name"], res.get("shape") or "", f"{res['input_mb']:,.1f}", "PAL2NAL v14",
                                 f"{p['wall_s']:.2f}", fmt_mb(p.get("peak_mb")), rc,
                                 "output" if p.get("out_bytes") else "no output"))
    print("| case | input shape | input MB | engine | time s | peak MB | exit | correct? |")
    print("|---|---|---:|---|---:|---:|---:|---|")
    for row in rows:
        print("| " + " | ".join(str(x) for x in row) + " |")
    print(f"\n{len(rows)} rows", file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv[1:])
