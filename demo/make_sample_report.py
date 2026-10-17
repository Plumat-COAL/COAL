"""Regenerate src/coal_app/sample_report.json: the GUI's "Show a sample report" / --demo report.

    python demo/make_sample_report.py

It is a real run of this app on examples/mito with genetic code 1 (the mistake a new user makes with
vertebrate mitochondrial genes), so every number in it comes from COAL itself. File paths are cut to file
names and coal_version is "demo" so it is never mistaken for a run of a particular release.
"""
from __future__ import annotations

import json
import os
import shlex
import shutil
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

from coal_app import core  # noqa: E402

DEMO_SOURCE = ("DEMO SAMPLE, not your data: examples/mito (24 vertebrate mitochondrial genes) run by the COAL app "
               "with genetic code 1, which is the wrong code for these genes")


def build() -> dict:
    p, d = core.example_paths("mito")
    t = tempfile.mkdtemp(prefix="coal_sample_")
    try:
        rep = core.run_coal(p, d, codontable=1, out_dir=t).report
    finally:
        shutil.rmtree(t, ignore_errors=True)
    names = {}
    for item in list(rep.get("inputs") or []) + [rep.get("output") or {}]:
        if item.get("path"):
            names[item["path"]] = os.path.basename(item["path"])
            item["path"] = names[item["path"]]
    # CHANGED(MRG-2): "command_line" is new (the engine's own report: the exact argv, with the run folder's paths)
    for k in ("command", "pal2nal_equivalent", "command_line"):
        if k not in rep:
            continue
        rep[k] = core.strip_folders(rep[k], names)  # CHANGED(SHIP-2): was replace(_q(full)), missing shlex quoting
        if k == "command_line":
            # CHANGED(SHIP-2): was a text replace of the run folder, which left Windows' quotes around paths that are
            # not inputs or outputs (--report). The engine writes this field with shlex.quote on every OS, so split
            # it, drop the run folder from each argument and quote it again.
            pre = os.path.join(t, "")
            rep[k] = " ".join(shlex.quote(a[len(pre):] if a.startswith(pre) else a) for a in shlex.split(rep[k]))
    rep["coal_version"] = "demo"
    # CHANGED(MRG-10): was kept; the merged engine now writes its own methods_text naming "COAL v0.4.0", which
    # contradicts coal_version "demo" and overrides the app's methods_text (APP-15, FIN-C11), so the sample drops it
    rep.pop("methods_text", None)
    rep["source"] = DEMO_SOURCE
    rep["demo"] = True
    rep["runtime_seconds"] = 0.0
    return rep


def main() -> int:
    rep = build()
    with open(core.SAMPLE_REPORT, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(rep, fh, indent=2)
        fh.write("\n")
    print(f"wrote {core.SAMPLE_REPORT}: {core.summary_sentence(rep)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
