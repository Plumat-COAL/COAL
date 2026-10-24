#!/usr/bin/env python3
"""Phase 13: can the tests catch real bugs? In a scratch copy OUTSIDE the repository, one realistic bug at a time, then
the relevant suites. Each bug the suites miss is a test gap.

    python verify/mutations.py SCRATCH_DIR OUTDIR
"""
from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PY = ROOT / ".venv" / "bin" / "python"

STD = "---M------**--*----M---------------M----------------------------"   # the standard code's start codons (NCBI)
NO_CTG = STD[:19] + "-" + STD[20:]
ONLY_ATG = STD[:3] + "-" + STD[4:19] + "-" + STD[20:]
T2 = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG"       # table 2 (TGA = W at index 14)
T2_BAD = T2[:14] + "*" + T2[15:]

# (name, file, old text, new text, suites, needs a C++ build)
MUTATIONS = [
    ("py: mismatch warning column off by one", "src/coal/core.py",
     "            res.add_mismatch(c, res.dna[3 * i:3 * i + 3])",
     "            res.add_mismatch(c + 1, res.dna[3 * i:3 * i + 3])", ["engine"], False),
    ("py: CTG no longer a standard-code start codon", "src/coal/tables.py", '"' + STD + '"', '"' + NO_CTG + '"', ["engine"], False),
    ("py: -nogap keeps TAG stop columns", "src/coal/core.py",
     "    bad.update(m.start() // 3 for m in _stop_lookahead(table).finditer(up) if m.start() % 3 == 0)",
     "    bad.update(m.start() // 3 for m in _stop_lookahead(table).finditer(up) if m.start() % 3 == 0 and up[m.start():m.start() + 3] != 'TAG')",
     ["engine"], False),
    ("py: table 2 TGA translated as a stop", "src/coal/tables.py", '"' + T2 + '"', '"' + T2_BAD + '"', ["engine"], False),
    ("py: exit code 0 after a CoalError", "src/coal/cli.py",
     '        say(f"error: {e}")\n        rc = 1\n', '        say(f"error: {e}")\n        rc = 0\n', ["engine"], False),
    ("py: report summary ok count off by one", "src/coal/report.py",
     '        summary = {"total": n, "ok": c["ok"],', '        summary = {"total": n, "ok": c["ok"] + 1,', ["engine"], False),
    ("py: -nomismatch drops the column after the mismatch", "src/coal/pipeline.py",
     "                    drop.update(c - 1 for c in mcols)", "                    drop.update(c for c in mcols)", ["engine"], False),
    ("py: genetic-code hint suggests the worst tied code", "src/coal/report.py",
     "        c = tied[0]\n        total, table = n, self.table", "        c = tied[-1]\n        total, table = n, self.table", ["engine"], False),
    ("py: paml header counts one sequence too many", "src/coal/pipeline.py",
     '    return " %3d %6d\\n" % (nseq, alilen)', '    return " %3d %6d\\n" % (nseq + 1, alilen)', ["engine"], False),
    ("py: --keep-going exit 0 instead of 2", "src/coal/cli.py",
     "        if skipped:\n            rc = 2", "        if skipped:\n            rc = 0", ["engine"], False),
    ("py: -blockonly keeps the column after each mark", "src/coal/pipeline.py",
     '        keep_col = (lambda c: mask[c - 1:c] == "#") if block_cols else None',
     '        keep_col = (lambda c: mask[c:c + 1] == "#") if block_cols else None', ["engine"], False),
    ("cpp: output rows emitted in reverse batch order with -j 4", "cpp/src/pool.hpp",
     "        emit(cur);\n        cur.swap(nxt);",
     "        if (nthreads > 1) std::reverse(cur.begin(), cur.end());\n        emit(cur);\n        cur.swap(nxt);", ["cpp"], True),
    ("cpp: only ATG accepted as the first Met (standard code)", "cpp/src/tables.cpp", '"' + STD + '"', '"' + ONLY_ATG + '"', ["cpp"], True),
    ("app: GUI download hands out the report instead of the alignment", "apps/coal-app/src/coal_app/ui.py",
     "        gr.update(value=output_path, visible=bool(output_path)),",
     "        gr.update(value=report_path, visible=bool(output_path)),", ["app"], False),
    ("app: MCP all_ok true with a failed gene", "apps/coal-app/src/coal_app/tools.py",
     "        all_ok=all_ok,", "        all_ok=True,", ["app"], False),
]


def suite_cmd(which, scratch):
    env = dict(os.environ, COAL_CPP_BIN=str(scratch / "build" / "coal-cpp"), PYTHONPATH=str(scratch / "src"))
    if which == "engine":
        return [str(PY), "-m", "unittest", "discover", "-s", "tests"], env
    if which == "cpp":
        env["COAL_IMPL"] = "cpp"
        return [str(PY), "-m", "unittest", "tests.test_differential", "tests.test_cpp_vs_python", "tests.test_cpp_pool",
                "tests.test_cpp_cli", "tests.test_decisions", "tests.test_cpp_vs_reference"], env
    if which == "app":
        return [str(PY), "-m", "unittest", "discover", "-s", "apps/coal-app/tests"], env
    raise ValueError(which)


def apply(name, scratch):
    """Apply one mutation to the scratch copy (an exact old/new text pair)."""
    for m in MUTATIONS:
        if m[0] == name:
            break
    _, file, old, new, suites, cpp = m
    path = scratch / file
    t = path.read_text(encoding="utf-8")
    if t.count(old) != 1:
        raise RuntimeError(f"{name}: {t.count(old)} matches of the old text in {file}")
    path.write_text(t.replace(old, new), encoding="utf-8")


def main(scratch_dir, outdir):
    scratch = Path(scratch_dir); out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    results = []
    for m in MUTATIONS:
        name, file, old, new, suites, cpp = m
        rec = {"name": name, "file": file, "suites": suites}
        t0 = time.monotonic()
        if scratch.exists():
            shutil.rmtree(scratch)
        subprocess.run(["rsync", "-a", "--exclude", ".git", "--exclude", ".venv", "--exclude", ".venv-mcp", "--exclude", "PAL2NAL",
                        "--exclude", "verify", "--exclude", "build", f"{ROOT}/", f"{scratch}/"], check=True)
        try:
            apply(name, scratch)
        except Exception as e:  # noqa: BLE001
            rec["applied"] = False; rec["error"] = repr(e); results.append(rec)
            print(f"NOT APPLIED {name}: {e!r}", flush=True); continue
        rec["applied"] = True
        if cpp:
            b = subprocess.run(["bash", "-lc", f"cd {scratch} && cmake -S . -B build -DCMAKE_BUILD_TYPE=Release > /dev/null && cmake --build build -j 16 2>&1 | tail -2"],
                               capture_output=True, text=True, timeout=900)
            rec["build"] = b.stdout[-300:]
            if not (scratch / "build" / "coal-cpp").exists():
                rec["caught"] = None; rec["note"] = "build failed"; results.append(rec); print(f"BUILD FAILED {name}", flush=True); continue
        caught = {}
        for which in suites:
            cmd, env = suite_cmd(which, scratch)
            try:
                r = subprocess.run(cmd, cwd=scratch, env=env, capture_output=True, text=True, timeout=1500)
                tail = r.stderr[-3000:]
                m2 = re.search(r"Ran (\d+) tests.*?\n\n(OK|FAILED \(([^)]*)\))", tail, re.S)
                summary = m2.group(2) if m2 else f"rc={r.returncode}"
                failed = [l for l in r.stderr.splitlines() if l.startswith(("FAIL:", "ERROR:"))]
                caught[which] = {"summary": summary, "failures": failed[:8], "n_failed": len(failed), "rc": r.returncode}
            except subprocess.TimeoutExpired:
                caught[which] = {"summary": "TIMEOUT", "failures": [], "n_failed": 0, "rc": None}
        rec["suites_result"] = caught
        rec["caught"] = any(v["rc"] != 0 for v in caught.values())
        rec["seconds"] = round(time.monotonic() - t0, 1)
        results.append(rec)
        print(f"{'CAUGHT     ' if rec['caught'] else 'NOT CAUGHT '} {name}: "
              + "; ".join(f"{k}: {v['summary']} ({v['n_failed']} failing: {', '.join(x.split(' ')[1] for x in v['failures'][:3])})" for k, v in caught.items()), flush=True)
    (out / "results.json").write_text(json.dumps(results, indent=1), encoding="utf-8")
    n = sum(1 for r in results if r.get("applied"))
    c = sum(1 for r in results if r.get("caught"))
    print(f"\n{n} bugs introduced, {c} caught, {n - c} test gaps; {len(results) - n} recipes could not be applied")
    if scratch.exists():
        shutil.rmtree(scratch)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
