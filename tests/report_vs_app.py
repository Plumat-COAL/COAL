"""Compare the engine's run report (coal --report) with the app's fallback report (src/coal_app/core.py, built from
the engine's stderr: build_fallback_report + mismatch_probe + locate_breaks + id_hints + genetic_code_hint).

CHANGED(MRG-2): new (merge reconciliation item 2). Run:  python tests/report_vs_app.py   (prints every difference)
Used by tests/test_merge.py. Needs src/coal_app (skips without it). Every engine call has a timeout (the app's).
"""
from __future__ import annotations

import json
import os
import sys
import tempfile
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

OPTION_SETS = [
    {},
    {"nomismatch": True},
    {"nogap": True},
    {"codontable": 2},
    {"codontable": 11},
    {"codontable": 13},
    {"codontable": 2, "nomismatch": True},
]
# What may legitimately differ, and why (each documented in overnight/integration/CHANGES.md, MRG-2):
#   hints "engine_warning"/"nomismatch_unknown": stderr-parsing artefacts of the fallback (the engine knows)
#   genetic_code hint: the engine tests every gene and all 17 codes (the fallback: 9 codes on the problem genes)
#   "not_checked": the fallback marks every gene of an aborted run "not checked"; the engine only those it never
#   reached (the others keep their status: "checked, but not written")


def inputs():
    """(name, protein, dna) for every example and every single-DNA-file fixture."""
    out = []
    for d in sorted((ROOT / "examples").iterdir()):
        if (d / "protein.aln.fasta").is_file():
            out.append((f"examples/{d.name}", d / "protein.aln.fasta", d / "dna.fasta"))
    cases = json.loads((ROOT / "tests" / "cases.json").read_text())["cases"]
    for name, spec in sorted(cases.items()):
        nuc = spec.get("nuc", ["nuc.fasta"])
        if len(nuc) == 1:
            out.append((f"fixtures/{name}", ROOT / "tests" / "fixtures" / name / spec["pep"],
                        ROOT / "tests" / "fixtures" / name / nuc[0]))
    return out


def run_both(app, protein, dna, opts, work):
    from coal.pipeline import FORMATS
    base = dict(output_formats=FORMATS, json=True, keep_going=False, version=_version())
    res = {}
    for kind, has_report in (("app", False), ("engine", True)):
        feats = app.Features(report=has_report, **base)
        d = os.path.join(work, kind)
        os.makedirs(d, exist_ok=True)
        with mock.patch.object(app, "detect_features", lambda: feats):
            r = app.run_coal(str(protein), str(dna), codontable=opts.get("codontable", 1),
                             nogap=opts.get("nogap", False), nomismatch=opts.get("nomismatch", False),
                             keep_going=False, out_dir=d, probe_codes=True)
            res[kind] = (r.report, app.methods_text(r.report), app.explain(r.report),
                         [app.sequence_sentence(s) for s in r.report.get("sequences", [])])
    return res


def _version():
    from coal import __version__
    return __version__


def _fits(protein, dna, table, nogap=False) -> int:
    """Genes that COAL aligns with no definite mismatch under `table` (an independent check of a hint)."""
    import subprocess
    from coal.report import is_ambiguous_codon
    r = subprocess.run([sys.executable, "-m", "coal", "--keep-going", "-codontable", str(table), "--", str(protein),
                        str(dna)], capture_output=True, text=True, timeout=120)
    from coal_app.core import read_ids
    bad = set()
    for line in r.stderr.splitlines():
        parts = line.split()
        if line.startswith("coal: warning: skipped "):
            bad.add(parts[3].rstrip(":"))
        elif len(parts) > 3 and parts[3] == "pepAlnPos" and not is_ambiguous_codon(parts[-1]):
            bad.add(parts[2])
        elif ": protein and DNA are probably unrelated" in line or ": position of the gene" in line:
            bad.add(parts[2].rstrip(":"))
    return sum(1 for i in read_ids(str(protein)) if i not in bad)


def differences(res, protein=None, dna=None, table=1) -> list[str]:
    """Differences between the two reports, except the documented ones (MRG-2, MRG-1):
    pairing: the fallback assumes order pairing; the engine pairs by normalised IDs where order would mispair;
    genetic-code hint: accepted when the engine's suggestion really fits more genes than the code used (checked by
    running it) and the fallback's suggestion, if any, fits no more genes than the engine's."""
    a, ma, _, sa = res["app"]
    e, me, _, se = res["engine"]
    diffs = []
    norm = e["summary"].get("pairing") == "normalised" and a["summary"].get("pairing") == "order"
    if norm:
        a = dict(a, summary=dict(a["summary"], pairing="normalised"),
                 hints=[h for h in a.get("hints") or [] if h["code"] != "id_mismatch"],
                 sequences=[dict(x, dna_id=y.get("dna_id")) for x, y in zip(a["sequences"], e["sequences"])])
        e = dict(e, hints=[h for h in e.get("hints") or [] if h["code"] != "normalised_ids"])
    skip_hints = {"engine_warning", "nomismatch_unknown"}
    sa_ = {k: a["summary"].get(k) for k in ("total", "ok", "warning", "failed", "pairing", "columns_removed",
                                              "columns_removed_ambiguous", "not_checked")}
    se_ = {k: e["summary"].get(k) for k in sa_}
    aborted = a.get("exit_code") not in (0,)
    if aborted:  # see the comment on "not_checked" above
        for k in ("ok", "warning", "not_checked", "failed"):
            sa_.pop(k), se_.pop(k)
    if sa_ != se_:
        diffs.append(f"summary {sa_} != {se_}")
    for x, y, sx, sy in zip(a["sequences"], e["sequences"], sa, se):
        if aborted and "not_checked" in (x["status"], y["status"]):
            continue
        cx = (x["id"], x["status"], x.get("dna_id"), [m.get("code") for m in x["messages"]])
        cy = (y["id"], y["status"], y.get("dna_id"), [m.get("code") for m in y["messages"]])
        if cx != cy:
            diffs.append(f"sequence {cx} != {cy}")
        elif sx != sy:
            diffs.append(f"sentence {x['id']}: {sx!r} != {sy!r}")
    if len(a["sequences"]) != len(e["sequences"]):
        diffs.append(f"{len(a['sequences'])} sequences != {len(e['sequences'])}")
    ha = [(h["code"], h.get("text")) for h in a.get("hints") or [] if h["code"] not in skip_hints | {"genetic_code"}]
    he = [(h["code"], h.get("text")) for h in e.get("hints") or [] if h["code"] not in skip_hints | {"genetic_code"}]
    if ha != he:
        diffs.append(f"hints {ha} != {he}")
    ga = next((h.get("suggested_codontable") for h in a.get("hints") or [] if h["code"] == "genetic_code"), None)
    ge = next((h.get("suggested_codontable") for h in e.get("hints") or [] if h["code"] == "genetic_code"), None)
    if ga != ge:
        ok = False
        if ge is not None and protein is not None:
            cur, fe = _fits(protein, dna, table), _fits(protein, dna, ge)
            fa = _fits(protein, dna, ga) if ga is not None else -1
            ok = fe > cur and fa <= fe
        if not ok:
            diffs.append(f"genetic-code hint {ga} != {ge}")
    if bool(a.get("output")) != bool(e.get("output")):
        diffs.append(f"output {a.get('output')} != {e.get('output')}")
    # CHANGED(ENG-UX-9): the engine puts a "fix the warnings first" note before the methods text when fewer than half
    # the genes are OK; the app's fallback text (built from stderr) has no such note. Compare the text after it.
    from coal.report import AS_IT_STANDS, NOT_READY
    if me and me.startswith(NOT_READY) and not (ma or "").startswith(NOT_READY):
        me = me.split(AS_IT_STANDS + " ", 1)[1]
    if ma != me:
        diffs.append(f"methods text\n    app:    {ma}\n    engine: {me}")
    return diffs


def compare_all(verbose=True, only=None):
    try:
        from coal_app import core as app
    except ImportError:
        return None
    total, bad = 0, {}
    with tempfile.TemporaryDirectory() as work:
        for name, protein, dna in inputs():
            if only and name not in only:
                continue
            for opts in OPTION_SETS:
                total += 1
                w = os.path.join(work, f"run{total}")
                os.makedirs(w)
                d = differences(run_both(app, protein, dna, opts, w), protein, dna, opts.get("codontable", 1))
                if d:
                    bad[(name, json.dumps(opts))] = d
                    if verbose:
                        print(f"== {name} {opts}")
                        for x in d:
                            print("  ", x)
    return total, bad


if __name__ == "__main__":
    r = compare_all()
    if r is None:
        print("coal_app not importable: skipped")
        sys.exit(0)
    total, bad = r
    print(f"{total} comparisons, {len(bad)} with differences")
    sys.exit(1 if bad else 0)
