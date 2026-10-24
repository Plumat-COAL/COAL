#!/usr/bin/env python3
"""Phase 10: docs and website. Every command example that names a file in examples/ is run and its shown output
compared; the live site is fetched (version, download links, checksums, broken links); numbers are traced.

    COAL_CPP_BIN=build/coal-cpp python verify/phase10.py OUTDIR
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import subprocess
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SITE = "https://plumat-coal.github.io/COAL/"
CHECKS = []
ENV = dict(os.environ, COAL_ENGINE="python", PYTHONPATH=str(ROOT / "src"))


def check(name, ok, detail=""):
    CHECKS.append({"name": name, "ok": bool(ok), "detail": str(detail)[:700]})
    print(("OK   " if ok else "FINDING-CANDIDATE ") + name + ("" if ok else f": {str(detail)[:300]}"), flush=True)


def fetch(url, timeout=60):
    req = urllib.request.Request(url, headers={"User-Agent": "coal-verify/1"}, method="GET")
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.status, r.read(), r.headers


def head(url):
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "coal-verify/1"}, method="HEAD")
        with urllib.request.urlopen(req, timeout=60) as r:
            return r.status
    except urllib.error.HTTPError as e:
        return e.code
    except Exception as e:  # noqa: BLE001
        return str(e)[:60]


def run(cmd, cwd=ROOT):
    return subprocess.run(cmd, cwd=cwd, shell=True, capture_output=True, timeout=600, env=ENV, executable="/bin/bash")


def command_examples():
    """(source file, command line) for every `coal ...` / `pal2nal ...` / `python -m coal ...` line that names
    an examples/ file in the docs, README, skill and site."""
    out = []
    files = [ROOT / "README.md", ROOT / ".claude/skills/coal/SKILL.md", ROOT / "apps/coal-app/desktop/README.md",
             *sorted((ROOT / "docs").glob("*.md")), *sorted((ROOT / "site").glob("*.html"))]
    for f in files:
        text = f.read_text(encoding="utf-8", errors="replace")
        if f.suffix == ".html":
            import html as h
            text = h.unescape(re.sub(r"<[^>]+>", "", text))
        for line in text.splitlines():
            s = line.strip().lstrip("$ ").strip()
            if re.match(r"^(coal|pal2nal|python -m coal)\b", s) and "examples/" in s and "|" not in s and ">" not in s.split("#")[0]:
                out.append((str(f.relative_to(ROOT)), s.split("#")[0].strip()))
    return out


def main(outdir):
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    # 1. command examples with examples/ files
    seen = set()
    results = {}
    for src, cmd in command_examples():
        if cmd in seen:
            continue
        seen.add(cmd)
        c = cmd.replace("-o codon.aln.fasta", f"-o {out}/codon.aln.fasta").replace("-o out.", f"-o {out}/out.")
        c = re.sub(r"--report (\S+)", lambda m: f"--report {out}/{Path(m.group(1)).name}", c)
        r = run(c)
        results[cmd] = {"source": src, "rc": r.returncode, "stdout_head": r.stdout[:200].decode("utf-8", "replace"),
                        "stderr": r.stderr.decode("utf-8", "replace").splitlines()[:40]}
        print(f"   [{src}] {cmd} -> rc {r.returncode}", flush=True)
    (out / "command_examples.json").write_text(json.dumps(results, indent=1))
    # expectations the docs state next to the commands
    mito_err = results.get("coal examples/mito/protein.aln.fasta examples/mito/dna.fasta -o codon.aln.fasta", {}).get("stderr", [])
    if mito_err:
        check("site/guide: the mito example prints 20 warnings then '... and 159 more warnings' (179 in all) and the summary with code 2",
              any("159 more warnings" in l for l in mito_err) and any("Under genetic code 2 (Vertebrate mitochondrial), all 24 sequences match" in l for l in mito_err),
              mito_err[-3:])
    r = run("coal examples/one_broken_gene/protein.aln.fasta examples/one_broken_gene/dna.fasta -o /dev/null 2>&1; echo rc=$?")
    t = r.stdout.decode("utf-8", "replace")
    check("site/index+guide: one_broken_gene without --keep-going: error for OG0001_sp14, frameshift at column 102, 6 not checked, exit 1",
          "OG0001_sp14: protein and DNA are inconsistent (no match found)" in t and "Frameshift at column 102: one base is missing" in t and "6 sequences were not checked" in t and "rc=1" in t, t[-600:])
    r = run(f"coal examples/one_broken_gene/protein.aln.fasta examples/one_broken_gene/dna.fasta --keep-going -o {out}/kg.fa 2>&1; echo rc=$?; grep -c '>' {out}/kg.fa")
    t = r.stdout.decode()
    check("site/guide: with --keep-going 19 genes written, exit 2", "rc=2" in t and t.strip().endswith("19"), t[-300:])
    r = run("coal examples/clean/protein.aln.fasta examples/clean/dna.fasta | grep -c '>'; coal examples/clean/protein.aln.fasta examples/clean/dna.fasta -q -o /dev/null; echo rc=$?")
    t = r.stdout.decode()
    check("site/guide: the clean example has 8 genes, exit 0", t.split()[0] == "8" and "rc=0" in t, t)
    r = run(f"coal examples/mito/protein.aln.fasta examples/mito/dna.fasta -codontable 2 -o {out}/mito2.fa 2>&1; echo rc=$?")
    t = r.stdout.decode()
    check("site/guide: the mito example with -codontable 2 is clean (no warnings, exit 0)", "warning" not in t and "rc=0" in t, t[-300:])
    # README's Python example
    r = subprocess.run([sys.executable, "-c", 'import coal\nrows, w = coal.codon_align([("a", "MK-V"), ("b", "MKLV")], [("a", "ATGAAAGTT"), ("b", "ATGAAACTGGTT")], table=1)\nprint(rows == [("a", "ATGAAA---GTT"), ("b", "ATGAAACTGGTT")], w)'], capture_output=True, text=True, env=ENV, cwd=ROOT)
    check("README: codon_align example gives the shown rows", r.stdout.startswith("True"), r.stdout + r.stderr[-200:])
    # README exit-code list vs usage.md
    readme = (ROOT / "README.md").read_text()
    check("README exit codes list 0/1/2/141/130/143/129", all(x in readme for x in ("0 success", "141 closed output pipe", "130 / 143 / 129 interrupted")), "")
    # 2. the live website
    pages = {}
    for page in ("", "install.html", "guide.html", "differences.html", "cite.html", "explainer.html"):
        try:
            st, body, hdr = fetch(SITE + page)
            pages[page or "index.html"] = body.decode("utf-8", "replace")
            check(f"site fetch {page or 'index.html'}: 200", st == 200, st)
        except Exception as e:  # noqa: BLE001
            check(f"site fetch {page or 'index.html'}: 200", False, repr(e))
    idx = pages.get("index.html", "")
    check("site index shows version 0.5.0", 'data-version="0.5.0"' in idx and ">0.5.0<" in idx, re.findall(r'data-version="[^"]*"', idx)[:2])
    local = (ROOT / "site" / "index.html").read_text(encoding="utf-8")
    check("site index equals the repository's site/index.html (deployed from main)", idx.strip() == local.strip(), f"live {len(idx)} bytes vs repo {len(local)} bytes")
    # links
    links = set()
    for page, body in pages.items():
        for href in re.findall(r'href="([^"#]+)"', body):
            if href.startswith("mailto:"):
                continue
            url = href if href.startswith("http") else SITE + href.lstrip("./")
            links.add(url)
    broken = {}
    for url in sorted(links):
        st = head(url)
        if st not in (200, 301, 302):
            broken[url] = st
    check(f"site: no broken links ({len(links)} checked)", not broken, broken)
    (out / "site_links.json").write_text(json.dumps({"checked": sorted(links), "broken": broken}, indent=1))
    # download links resolve to release files; checksums
    dl = [u for u in links if "releases/" in u and ("download" in u)]
    check("site: download links present", bool(dl), dl[:5])
    sums = {}
    try:
        st, body, _ = fetch("https://github.com/Plumat-COAL/COAL/releases/latest/download/SHA256SUMS.txt")
        for line in body.decode().splitlines():
            h_, name = line.split()
            sums[name.lstrip("*")] = h_
        check("SHA256SUMS.txt downloadable from the latest release", bool(sums), len(sums))
    except Exception as e:  # noqa: BLE001
        check("SHA256SUMS.txt downloadable from the latest release", False, repr(e))
    for u in dl:
        name = u.rsplit("/", 1)[-1]
        if name in ("SHA256SUMS.txt",) or name not in sums:
            continue
        if name.endswith((".exe", ".zip")):
            continue  # the installer and the zip were verified on Windows (phase 8 notes); sdists are small
        try:
            st, body, _ = fetch(u, timeout=300)
            check(f"site download {name}: sha256 matches SHA256SUMS", hashlib.sha256(body).hexdigest() == sums[name], sums[name])
        except Exception as e:  # noqa: BLE001
            check(f"site download {name}", False, repr(e))
    win = [u for u in dl if u.endswith(".exe") or u.endswith(".zip")]
    check("site: Windows installer and zip links are release assets", all(head(u) in (200, 302) for u in win), win)
    # install steps per OS present
    inst = pages.get("install.html", "")
    check("site install: Windows / macOS / Linux steps, pip commands and checksum instructions", all(x in inst for x in ("Windows", "macOS", "Linux", "pip install", "SHA256SUMS")), "")
    # 3. claims the team decided not to make
    bad = {}
    for f in [ROOT / "README.md", *sorted((ROOT / "docs").glob("*.md")), *sorted((ROOT / "site").glob("*.html")), ROOT / ".claude/skills/coal/SKILL.md", ROOT / "apps/coal-app/desktop/README.md"]:
        t = f.read_text(encoding="utf-8", errors="replace")
        for claim, pat in (("NVIDIA Brev", r"Brev"), ("hosted web version", r"hosted (web )?version|try it online|live demo"), ("code signing", r"code[- ]signed(?! yet)|digitally signed(?! yet)"), ("DIAMOND on real genomes", r"DIAMOND .*real genome")):
            m = re.search(pat, t)
            if m:
                bad.setdefault(claim, []).append(f"{f.relative_to(ROOT)}: {t[max(0, m.start() - 60):m.end() + 60]!r}")
    check("claims absent: NVIDIA Brev, hosted web version, code signing, DIAMOND on real genomes", not bad, bad)
    pypi = re.findall(r"pip install coal-align\b(?! @)", readme + (ROOT / "docs/install.md").read_text())
    check("README/install: 'pip install coal-align' only as a future (PyPI) statement", all("Once" in readme[max(0, readme.find(p) - 80):readme.find(p)] or "Once" in (ROOT / "docs/install.md").read_text()[max(0, (ROOT / "docs/install.md").read_text().find(p) - 80):(ROOT / "docs/install.md").read_text().find(p)] for p in pypi) if pypi else True, pypi)
    (out / "results.json").write_text(json.dumps(CHECKS, indent=1))
    bad = [c for c in CHECKS if not c["ok"]]
    print(f"\n{len(CHECKS)} checks, {len(bad)} finding candidates")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "phase10_out")
