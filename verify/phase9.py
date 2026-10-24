#!/usr/bin/env python3
"""Phase 9: distribution and versions. Clean venvs on the oldest and newest supported Pythons, the release wheels and
sdists, `pip install git+https://...`, the platform wheel with coal-cpp inside, and every version string.

    python verify/phase9.py OUTDIR
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PY310 = Path.home() / "miniconda3" / "envs" / "py310" / "bin" / "python"
PY313 = Path.home() / "miniconda3" / "envs" / "py313" / "bin" / "python3.13"
CHECKS = []


def check(name, ok, detail=""):
    CHECKS.append({"name": name, "ok": bool(ok), "detail": str(detail)[:600]})
    print(("OK   " if ok else "FINDING-CANDIDATE ") + name + ("" if ok else f": {str(detail)[:300]}"), flush=True)


def sh(cmd, cwd=None, timeout=900, env=None):
    return subprocess.run(cmd, cwd=cwd, capture_output=True, timeout=timeout, env=env, text=True)


def venv(py, where):
    r = sh([str(py), "-m", "venv", str(where)])
    return where / "bin" / "python"


def main(outdir):
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    rel = out / "release"; rel.mkdir(exist_ok=True)
    # 1. release assets
    r = sh(["gh", "release", "download", "v0.5.0", "-R", "Plumat-COAL/COAL", "-p", "*.whl", "-p", "*.tar.gz", "-p", "SHA256SUMS.txt", "--clobber"], cwd=rel, timeout=600)
    check("release v0.5.0 assets downloaded", r.returncode == 0, r.stderr[-300:])
    sums = (rel / "SHA256SUMS.txt").read_text().splitlines()
    import hashlib
    bad = []
    for line in sums:
        h, name = line.split()
        name = name.lstrip("*")
        f = rel / name
        if f.exists() and hashlib.sha256(f.read_bytes()).hexdigest() != h:
            bad.append(name)
    check("SHA256SUMS.txt matches every downloaded asset", not bad, bad)
    pure = rel / "coal_align-0.5.0-py3-none-any.whl"
    app = rel / "coal_align_app-0.5.0-py3-none-any.whl"
    plat = rel / "coal_align-0.5.0-py3-none-manylinux_2_28_x86_64.whl"
    sdist = rel / "coal_align-0.5.0.tar.gz"
    app_sdist = rel / "coal_align_app-0.5.0.tar.gz"
    for name, py in (("3.10", PY310), ("3.13", PY313)):
        if not py.exists():
            check(f"python {name} available", False, py); continue
        v = venv(py, out / f"venv-wheel-{name}")
        r = sh([str(v), "-m", "pip", "install", "-q", str(pure), str(app)])
        check(f"py{name}: release pure wheel + app wheel install", r.returncode == 0, r.stderr[-300:])
        r = sh([str(v), "-m", "coal", "--version"])
        check(f"py{name}: coal --version == coal 0.5.0", r.stdout.strip() == "coal 0.5.0", r.stdout + r.stderr[-200:])
        r = sh([str(v), "-m", "coal", str(ROOT / "examples" / "clean" / "protein.aln.fasta"), str(ROOT / "examples" / "clean" / "dna.fasta")])
        check(f"py{name}: the clean example runs from the wheel", r.returncode == 0 and r.stdout.startswith(">"), r.stderr[-200:])
        r = sh([str(v / ".." / "pal2nal").resolve().__str__() if False else str(v.parent / "pal2nal"), "-h"])
        check(f"py{name}: pal2nal entry point installed", r.returncode == 0, r.stderr[-200:])
        r = sh([str(v), "-m", "coal_app", "--check"])
        check(f"py{name}: python -m coal_app --check works", r.returncode == 0, (r.stdout + r.stderr)[-300:])
        r = sh([str(v), "-m", "coal", "--cpp-path"])
        check(f"py{name}: pure wheel has no coal-cpp (--cpp-path exits 1 with a reason)", r.returncode == 1 and "not found" in r.stderr, (r.returncode, r.stderr[-200:]))
        # platform wheel: coal-cpp inside
        v2 = venv(py, out / f"venv-plat-{name}")
        r = sh([str(v2), "-m", "pip", "install", "-q", str(plat)])
        check(f"py{name}: manylinux platform wheel installs", r.returncode == 0, r.stderr[-300:])
        r = sh([str(v2), "-m", "coal", "--cpp-path"])
        check(f"py{name}: --cpp-path finds the bundled coal-cpp", r.returncode == 0 and "coal/bin/coal-cpp" in r.stdout, (r.returncode, r.stdout, r.stderr[-200:]))
        cpp = r.stdout.strip()
        if cpp:
            r = sh([cpp, "--version"])
            check(f"py{name}: bundled coal-cpp --version == coal 0.5.0 (C++)", r.stdout.strip() == "coal 0.5.0 (C++)", r.stdout + r.stderr[-200:])
            r = sh([cpp, str(ROOT / "examples" / "clean" / "protein.aln.fasta"), str(ROOT / "examples" / "clean" / "dna.fasta")])
            r2 = sh([str(v2), "-m", "coal", str(ROOT / "examples" / "clean" / "protein.aln.fasta"), str(ROOT / "examples" / "clean" / "dna.fasta")], env=dict(os.environ, COAL_ENGINE="python"))
            check(f"py{name}: bundled coal-cpp output equals the Python engine's on the clean example", r.returncode == 0 and r.stdout == r2.stdout, r.stderr[-200:])
        # sdist
        v3 = venv(py, out / f"venv-sdist-{name}")
        r = sh([str(v3), "-m", "pip", "install", "-q", str(sdist), str(app_sdist)])
        check(f"py{name}: release sdists install", r.returncode == 0, r.stderr[-300:])
        r = sh([str(v3), "-m", "coal", "--version"])
        check(f"py{name}: sdist coal --version", r.stdout.strip() == "coal 0.5.0", r.stdout + r.stderr[-200:])
    # 2. sdist built from the checkout
    v = venv(PY313, out / "venv-build")
    sh([str(v), "-m", "pip", "install", "-q", "build"])
    r = sh([str(v), "-m", "build", "--sdist", "--outdir", str(out / "dist")], cwd=ROOT, timeout=900)
    check("python -m build --sdist from the checkout", r.returncode == 0, r.stderr[-400:])
    built = sorted((out / "dist").glob("*.tar.gz"))
    if built:
        r = sh([str(v), "-m", "pip", "install", "-q", str(built[0])])
        check("built sdist installs", r.returncode == 0, r.stderr[-300:])
    # 3. git+https
    v = venv(PY313, out / "venv-git")
    r = sh([str(v), "-m", "pip", "install", "-q", "git+https://github.com/Plumat-COAL/COAL"], timeout=900)
    check("pip install git+https://github.com/Plumat-COAL/COAL", r.returncode == 0, r.stderr[-300:])
    r = sh([str(v), "-m", "coal", "--version"])
    check("git install: coal --version", r.stdout.strip() == "coal 0.5.0", r.stdout + r.stderr[-200:])
    # 4. version strings
    v = {}
    v["src/coal/__init__.py"] = re.search(r'__version__ = "([^"]+)"', (ROOT / "src/coal/__init__.py").read_text()).group(1)
    v["pyproject.toml"] = re.search(r'^version = "([^"]+)"', (ROOT / "pyproject.toml").read_text(), re.M).group(1)
    v["apps/coal-app/pyproject.toml"] = re.search(r'^version = "([^"]+)"', (ROOT / "apps/coal-app/pyproject.toml").read_text(), re.M).group(1)
    v["app pin coal-align=="] = re.search(r'coal-align==([0-9.]+)', (ROOT / "apps/coal-app/pyproject.toml").read_text()).group(1)
    v["cpp/src/main.cpp"] = re.search(r'VERSION = "([^"]+)"', (ROOT / "cpp/src/main.cpp").read_text()).group(1)
    v["cpp/src/report.cpp"] = re.search(r'VERSION = "([^"]+)"', (ROOT / "cpp/src/report.cpp").read_text()).group(1)
    v["CMakeLists.txt"] = re.search(r'project\(coal_cpp VERSION ([0-9.]+)', (ROOT / "CMakeLists.txt").read_text()).group(1)
    v["CITATION.cff"] = re.search(r'^version: ([0-9.]+)', (ROOT / "CITATION.cff").read_text(), re.M).group(1)
    v["site/index.html"] = re.search(r'data-version="([^"]+)"', (ROOT / "site/index.html").read_text()).group(1)
    v["site/cite.html"] = re.search(r'data-version-ref>([0-9.]+)', (ROOT / "site/cite.html").read_text()).group(1)
    v["tests/report_example.json coal_version"] = json.loads((ROOT / "tests/report_example.json").read_text())["coal_version"]
    v["CHANGELOG first release heading"] = re.search(r'^## ([0-9][0-9.]*)', (ROOT / "CHANGELOG.md").read_text(), re.M).group(1)
    r = sh([str(ROOT / ".venv" / "bin" / "python"), "-m", "coal", "--version"]); v["coal --version"] = r.stdout.strip().split()[-1]
    r = sh([str(ROOT / "build" / "coal-cpp"), "--version"]); v["coal-cpp --version"] = r.stdout.strip().split()[1]
    r = sh(["gh", "release", "view", "-R", "Plumat-COAL/COAL", "--json", "tagName", "-q", ".tagName"]); v["latest release tag"] = r.stdout.strip().lstrip("v")
    m = re.search(r"COAL-Setup-([0-9.]+)\.exe", (ROOT / "apps/coal-app/desktop/README.md").read_text()); v["desktop/README.md installer name"] = m.group(1) if m else "none"
    skill = (ROOT / ".claude/skills/coal/SKILL.md").read_text()
    v["SKILL.md"] = "no version string (not required)" if not re.search(r"\b0\.[0-9]\.[0-9]\b", skill) else re.search(r"\b(0\.[0-9]\.[0-9])\b", skill).group(1)
    (out / "versions.json").write_text(json.dumps(v, indent=1))
    vals = {k: x for k, x in v.items() if x not in ("no version string (not required)",)}
    distinct = sorted(set(vals.values()))
    check("every version string agrees (0.5.0)", distinct == ["0.5.0"], {k: x for k, x in vals.items() if x != "0.5.0"})
    (out / "results.json").write_text(json.dumps(CHECKS, indent=1))
    bad = [c for c in CHECKS if not c["ok"]]
    print(f"\n{len(CHECKS)} checks, {len(bad)} finding candidates")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "phase9_out")
