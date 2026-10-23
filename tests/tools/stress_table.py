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


# What a non-zero exit or a difference from the truth means for these cases (by prefix of the case name).
NOTES = [
    ("fmt_empty_header", "error by design: '>' without a name"),
    ("fmt_ids_nonascii", "error by design: non-ASCII ID"),
    ("fmt_empty_dna_record_keepgoing", "by design: the empty gene is left out (exit 2)"),
    ("fmt_empty_dna_record", "error by design: a gene without DNA"),
    ("keepgoing_100k_one_empty", "by design: the empty gene is left out (exit 2)"),
    ("strict_100k_one_empty", "error by design: a gene without DNA"),
    ("fmt_mixed_lf_then_cr", "FINDING 4: PAL2NAL reads this file"),
    ("files_1000", "FINDING 3: open-file limit (Python engine with workers)"),
    ("files_5000", "FINDING 3: open-file limit (both engines, clean error)"),
    ("files_70000_ulimit_65536", "FINDING 2 (Python, fixed) / FINDING 3 (both)"),
    ("files_70000", "FINDING 3: open-file limit"),
    ("files_60000", "FINDING 3: 129 KB per open file (Python)"),
    ("repeat_polyA_prefix_1000000", "FINDING 10: quadratic scan, killed at 30 min"),
    ("repeat_polyA_prefix", "FINDING 10: quadratic scan in coal-cpp"),
    ("repeat_GPP_prefix", "FINDING 10: quadratic scan in coal-cpp"),
    ("unrelated_1M", "FINDING 9: killed by the 16 GB memory watchdog"),
    ("unrelated_100k", "FINDING 9: memory per warning"),
    ("wrongcode_1M", "FINDING 9: memory per warning"),
    ("repeat_polyA_late_", "first run: the generator's truth was wrong here (COAL's leftmost exact fit is right); superseded by repeat_polyA_prefix"),
    ("repeat_GPP_late_", "first run: the generator's truth was wrong here (COAL's leftmost exact fit is right); superseded by repeat_GPP_prefix"),
    ("fmt_dot_gap", "documented difference: '.' is a gap for COAL, PAL2NAL misaligns the row"),
    ("fmt_bom", "documented: PAL2NAL fails on a BOM"),
    ("fmt_ids_nonascii", "error by design: non-ASCII ID (PAL2NAL accepts it)"),
    ("opt_100k_progress", "--progress lines are timing-based (expected)"),
    ("opt_1M_progress", "--progress lines are timing-based (expected)"),
    ("unrelated_100k_strict_keepgoing", "by design: every gene refused with --strict"),
    ("half_ids_1M_strict", "error by design: partial ID match with --strict"),
]


def note_for(name: str) -> str:
    for prefix, text in NOTES:
        if name.startswith(prefix):
            return text
    return ""


def correct(r, res, eng):
    """What the row's 'correct?' column says."""
    parts = []
    if r.get("timeout"):
        return "did not finish in 30 min (killed); " + note_for(res["name"])
    if r.get("rc") == -9:
        return "killed by the memory watchdog at 16 GB; " + note_for(res["name"])
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
    text = ", ".join(parts) if parts else ("ran" if r.get("rc") == 0 else "error")
    note = note_for(res["name"])
    if note and (r.get("rc") not in (0, None) or "DIFF" in text or "TRACEBACK" in text or "timeout" in text
                 or note.startswith("FINDING") or "progress" in note):
        text += f"; {note}"
    return text


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
                    corrupt = any("corrupt" in x for rnd in rs["rounds"] for x in rnd)
                    if rs["ok"]:
                        verdict = ("OK: every damaged checkpoint refused as documented, output byte-identical" if corrupt
                                   else "OK: byte-identical after every kill/resume, no leftovers")
                    elif all(p.get("had_checkpoint") is False for p in rs["problems"]):
                        verdict = (f"OK: {len(rs['problems'])} corruption round(s) killed before the first checkpoint "
                                   f"(nothing to corrupt; normal run, output byte-identical)")
                    else:
                        verdict = f"PROBLEMS {len(rs['problems'])}"
                    what = "corruption rounds" if corrupt else "kill/resume rounds"
                    rows.append((res["name"], f"1,000,000 x 400, {' '.join(rs['args'])}, {what}", "", rs["engine"],
                                 f"{rs['ref_seconds']} (ref)", "", f"{n} runs", verdict))
                    continue
                for eng, r in res["engines"].items():
                    rc = "timeout" if r.get("timeout") else str(r.get("rc"))
                    rows.append((res["name"], res.get("shape") or "", f"{res['input_mb']:,.1f}", eng,
                                 f"{r['wall_s']:.2f}", fmt_mb(r.get("peak_mb")), rc, correct(r, res, eng)))
                p = res.get("pal2nal")
                if p:
                    rc = "timeout" if p.get("timeout") else str(p.get("rc"))
                    pal = "output" if p.get("out_bytes") else "no output"
                    if not p.get("out_bytes") and "50000" in res["name"]:
                        pal += " (fails above about 22,000 residues: documented)"
                    rows.append((res["name"], res.get("shape") or "", f"{res['input_mb']:,.1f}", "PAL2NAL v14",
                                 f"{p['wall_s']:.2f}", fmt_mb(p.get("peak_mb")), rc, pal))
    print("| case | input shape | input MB | engine | time s | peak MB | exit | correct? |")
    print("|---|---|---:|---|---:|---:|---:|---|")
    for row in rows:
        print("| " + " | ".join(str(x) for x in row) + " |")
    print(f"\n{len(rows)} rows", file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv[1:])
