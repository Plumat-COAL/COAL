#!/usr/bin/env python3
"""Phase 13 follow-up: the two app-level mutants again, with the app suite run from BOTH app environments
(.venv has the [gui] extra, .venv-mcp the [mcp] extra), reporting skips, so a 'not caught' is not just a skipped suite.

    python verify/mutations_app.py SCRATCH_DIR OUTDIR
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

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mutations  # noqa: E402

ROOT = mutations.ROOT
VENVS = {"gui": ROOT / ".venv" / "bin" / "python", "mcp": ROOT / ".venv-mcp" / "bin" / "python"}


def main(scratch_dir, outdir):
    scratch = Path(scratch_dir); out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    results = []
    for m in mutations.MUTATIONS:
        name, file, old, new, suites, cpp = m
        if suites != ["app"]:
            continue
        rec = {"name": name, "file": file}
        if scratch.exists():
            shutil.rmtree(scratch)
        subprocess.run(["rsync", "-a", "--exclude", ".git", "--exclude", ".venv", "--exclude", ".venv-mcp", "--exclude", "PAL2NAL",
                        "--exclude", "verify", "--exclude", "build", f"{ROOT}/", f"{scratch}/"], check=True)
        mutations.apply(name, scratch)
        for venv, py in VENVS.items():
            env = dict(os.environ, COAL_CPP_BIN=str(scratch / "build" / "coal-cpp"), PYTHONPATH=f"{scratch}/src:{scratch}/apps/coal-app/src")
            t0 = time.monotonic()
            r = subprocess.run([str(py), "-m", "unittest", "discover", "-s", "apps/coal-app/tests", "-v"], cwd=scratch, env=env,
                               capture_output=True, text=True, timeout=1500)
            tail = r.stderr
            ran = re.search(r"Ran (\d+) tests", tail); summ = re.search(r"\n\n(OK.*|FAILED.*)\n?$", tail)
            failed = [l for l in tail.splitlines() if l.startswith(("FAIL:", "ERROR:"))]
            skipped = re.findall(r"\.\.\. skipped '([^']*)'", tail)
            rec[venv] = {"ran": int(ran.group(1)) if ran else None, "summary": summ.group(1).strip() if summ else f"rc={r.returncode}",
                         "rc": r.returncode, "failures": failed[:10], "n_skipped": len(skipped), "skip_reasons": sorted(set(skipped))[:8],
                         "seconds": round(time.monotonic() - t0, 1)}
            print(f"{name} [{venv}]: ran {rec[venv]['ran']}, {rec[venv]['summary']}, {len(skipped)} skipped, failing: {', '.join(x.split(' ')[1] for x in failed[:4])}", flush=True)
        rec["caught"] = any(rec[v]["rc"] != 0 for v in VENVS)
        results.append(rec)
    (out / "results_app.json").write_text(json.dumps(results, indent=1), encoding="utf-8")
    if scratch.exists():
        shutil.rmtree(scratch)
    print("\n" + "; ".join(f"{r['name']}: {'CAUGHT' if r['caught'] else 'NOT CAUGHT'}" for r in results))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
