#!/usr/bin/env python3
"""Phase 8, the parts a script can drive on Windows: the installed app's command-line modes, its stdio MCP server,
--connect-claude / --disconnect-claude against a config with another server, the local app (--serve) driven with
Playwright (examples, a path with spaces and a OneDrive-style folder, a long path, a >= 1 MB input, dark theme,
binding to 127.0.0.1), the first launch of the real window (time to ready), and the portable zip.

    python phase8.py OUTDIR
Install / uninstall steps are in phase8_install.ps1 (run after this)."""
from __future__ import annotations

import json
import os
import queue
import re
import shutil
import socket
import subprocess
import sys
import threading
import time
import zipfile
from pathlib import Path

sys.path.insert(0, r"D:\coal-stress\verify")
import phase7  # noqa: E402
from phase7 import CHECKS, check, new_page, run_example, banner, health, click_button, download, wait_results  # noqa: E402

OUT = Path(sys.argv[1] if len(sys.argv) > 1 else r"D:\coal-stress\verify\p8")
APP = Path(os.environ["LOCALAPPDATA"]) / "Programs" / "COAL" / "COAL.exe"
DATA = Path(os.environ["LOCALAPPDATA"]) / "COAL"
CLAUDE = Path(os.environ["APPDATA"]) / "Claude" / "claude_desktop_config.json"
REPO = Path(r"D:\coal-handoff-kit\COAL")
EX = REPO / "examples"
BIG = Path(r"D:\coal-stress\verify\bench-data\n10000_L500")   # 13.7 MB of DNA: coal-cpp territory
RELEASE = Path(r"D:\coal-handoff-kit\release")
texts = {}


def run(cmd, timeout=120, **kw):
    try:
        r = subprocess.run(cmd, capture_output=True, timeout=timeout, **kw)
        return r.returncode, r.stdout.decode("utf-8", "replace"), r.stderr.decode("utf-8", "replace")
    except subprocess.TimeoutExpired as e:
        return None, (e.stdout or b"").decode("utf-8", "replace"), f"TIMEOUT after {timeout}s"


class Client:
    """A stdio MCP client for Windows (a reader thread instead of select)."""

    def __init__(self, cmd, cwd=None, timeout=90):
        self.p = subprocess.Popen(cmd, cwd=cwd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.q: queue.Queue = queue.Queue(); self.timeout = timeout; self.n = 0; self.stray = []
        threading.Thread(target=self._reader, daemon=True).start()

    def _reader(self):
        for line in iter(self.p.stdout.readline, b""):
            self.q.put(line)
        self.q.put(None)

    def send(self, method, params=None, notify=False):
        msg = {"jsonrpc": "2.0", "method": method}
        if params is not None:
            msg["params"] = params
        if not notify:
            self.n += 1; msg["id"] = self.n
        self.p.stdin.write((json.dumps(msg) + "\n").encode()); self.p.stdin.flush()
        return None if notify else self.recv(self.n)

    def recv(self, want_id):
        deadline = time.monotonic() + self.timeout
        while time.monotonic() < deadline:
            try:
                line = self.q.get(timeout=0.5)
            except queue.Empty:
                if self.p.poll() is not None:
                    raise RuntimeError(f"server exited {self.p.returncode}")
                continue
            if line is None:
                raise RuntimeError("server closed stdout")
            try:
                m = json.loads(line)
            except ValueError:
                self.stray.append(line[:200]); continue
            if m.get("id") == want_id:
                return m
        raise TimeoutError(f"no reply to {want_id} in {self.timeout}s")

    def call(self, name, args):
        return self.send("tools/call", {"name": name, "arguments": args})

    def close(self):
        try:
            self.p.stdin.close()
        except Exception:  # noqa: BLE001
            pass
        try:
            self.p.wait(timeout=20)
        except subprocess.TimeoutExpired:
            self.p.kill()
        return self.p.stderr.read()


def result_of(m):
    r = m.get("result", {})
    if "structuredContent" in r:
        d = dict(r["structuredContent"]); d.setdefault("isError", r.get("isError")); return d
    txt = "".join(c.get("text", "") for c in r.get("content", []) if isinstance(c, dict))
    try:
        d = json.loads(txt); d.setdefault("isError", r.get("isError")); return d
    except ValueError:
        return {"text": txt, "isError": r.get("isError", "error" in m)}


def mcp():
    c = Client([str(APP), "-m", "coal_app", "--mcp"], cwd=str(OUT))
    try:
        init = c.send("initialize", {"protocolVersion": "2025-06-18", "capabilities": {}, "clientInfo": {"name": "verify", "version": "0"}})
        check("installed COAL.exe -m coal_app --mcp: initialize answers", "result" in init, str(init)[:200])
        c.send("notifications/initialized", notify=True)
        tools = c.send("tools/list")
        names = sorted(t["name"] for t in tools.get("result", {}).get("tools", []))
        check("installed MCP: tools/list gives codon_align and explain_report", names == ["codon_align", "explain_report"], names)
        r = result_of(c.call("codon_align", {"protein_alignment_path": str(EX / "clean" / "protein.aln.fasta"), "dna_fasta_path": str(EX / "clean" / "dna.fasta")}))
        check("installed MCP: codon_align on the clean example completes, all_ok", r.get("completed") is True and r.get("all_ok") is True, str(r)[:300])
        texts["mcp_clean"] = r
        if r.get("report_path"):
            e = result_of(c.call("explain_report", {"report_path": r["report_path"]}))
            check("installed MCP: explain_report on that report is plain text", isinstance(e.get("text"), str) and len(e["text"]) > 40 and not e.get("isError"), str(e)[:200])
        r = result_of(c.call("codon_align", {"protein_alignment_path": str(BIG / "pep.fasta"), "dna_fasta_path": str(BIG / "nuc.fasta")}))
        rep = json.loads(Path(r["report_path"]).read_text(encoding="utf-8")) if r.get("report_path") else {}
        check("installed MCP: a 13.7 MB input runs in the bundled coal-cpp (report engine coal-cpp)", rep.get("engine") == "coal-cpp" and r.get("all_ok") is True, (rep.get("engine"), rep.get("command", "")[:120] if isinstance(rep.get("command"), str) else str(rep.get("command"))[:120]))
        texts["mcp_big_report_command"] = rep.get("command")
    except Exception as e:  # noqa: BLE001
        check("installed MCP server", False, repr(e)[:300])
    finally:
        err = c.close()
    check("installed MCP: nothing on stdout outside the protocol", not c.stray, c.stray[:3])
    check("installed MCP: no traceback on stderr", b"Traceback" not in err, err[-300:])
    (OUT / "mcp_stderr.txt").write_bytes(err)


def claude_config():
    orig = CLAUDE.read_bytes() if CLAUDE.exists() else None
    bak = CLAUDE.with_name(CLAUDE.name + ".bak-coal")
    try:
        other = {"mcpServers": {"other": {"command": "C:\\\\x\\\\other.exe", "args": ["--serve"]}}, "theme": "dark"}
        CLAUDE.write_text(json.dumps(other, indent=2), encoding="utf-8")
        if bak.exists():
            bak.unlink()
        rc, so, se = run([str(APP), "--connect-claude", "--quiet"], timeout=60)
        texts["connect_out"] = (rc, so[-400:], se[-400:])
        d = json.loads(CLAUDE.read_text(encoding="utf-8"))
        coal = d.get("mcpServers", {}).get("coal", {})
        check("--connect-claude adds a coal entry (COAL.exe -m coal_app --mcp)", coal.get("command", "").lower() == str(APP).lower() and coal.get("args") == ["-m", "coal_app", "--mcp"], coal)
        check("--connect-claude keeps the other server and other keys", d.get("mcpServers", {}).get("other") == other["mcpServers"]["other"] and d.get("theme") == "dark", str(d)[:300])
        check("--connect-claude writes claude_desktop_config.json.bak-coal, byte for byte the old file", bak.exists() and json.loads(bak.read_text(encoding="utf-8")) == other, bak.exists())
        check("--connect-claude exits 0", rc == 0, (rc, se[-200:]))
        rc2, so2, se2 = run([str(APP), "--connect-claude", "--quiet"], timeout=60)
        d2 = json.loads(CLAUDE.read_text(encoding="utf-8"))
        check("--connect-claude twice: one coal entry, file unchanged", d2 == d and rc2 == 0, (rc2, se2[-200:]))
        rc3, so3, se3 = run([str(APP), "--disconnect-claude", "--quiet"], timeout=60)
        d3 = json.loads(CLAUDE.read_text(encoding="utf-8"))
        check("--disconnect-claude removes coal, keeps the other server and keys", "coal" not in d3.get("mcpServers", {}) and d3.get("mcpServers", {}).get("other") == other["mcpServers"]["other"] and d3.get("theme") == "dark" and rc3 == 0, (rc3, str(d3)[:200], se3[-200:]))
        texts["disconnect_out"] = (rc3, so3[-400:], se3[-400:])
    except Exception as e:  # noqa: BLE001
        check("--connect-claude / --disconnect-claude", False, repr(e)[:300])
    finally:
        if orig is not None:
            CLAUDE.write_bytes(orig)
        if bak.exists():
            bak.unlink()


def wait_port(port, timeout=60):
    t0 = time.monotonic()
    while time.monotonic() - t0 < timeout:
        try:
            with socket.create_connection(("127.0.0.1", port), timeout=1):
                return time.monotonic() - t0
        except OSError:
            time.sleep(0.3)
    return None


def serve_and_drive(exe, port, label, full=True):
    p = subprocess.Popen([str(exe), "--serve", "--port", str(port)], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    ready = None; t0 = time.monotonic(); line = b""
    while time.monotonic() - t0 < 90:
        line = p.stdout.readline()
        if not line or line.startswith(b"COAL-READY"):
            ready = time.monotonic() - t0; break
    check(f"{label}: --serve prints COAL-READY <url>", bool(line and line.startswith(b"COAL-READY")), (line[:100], f"{ready}"))
    texts[f"{label}_ready_s"] = ready
    try:
        if not (line and line.startswith(b"COAL-READY")):
            return
        url = line.decode().split()[1].rstrip("/") + "/"
        rc, so, se = run(["netstat", "-ano", "-p", "tcp"], timeout=30)
        bound = [l.strip() for l in so.splitlines() if f":{port} " in l and "LISTEN" in l]
        check(f"{label}: the port is bound to 127.0.0.1 only", bound and all(l.split()[1].startswith("127.0.0.1:") for l in bound), bound[:3])
        phase7.URL = url
        from playwright.sync_api import sync_playwright
        with sync_playwright() as pw:
            browser = pw.chromium.launch(); ctx = browser.new_context()
            pg = new_page(ctx)
            for name, want in (("Clean nuclear genes", "8"), ("Mitochondrial (wrong code)", "24"), ("Messy IDs", ""), ("One broken gene", "19 of 20")):
                try:
                    t = run_example(pg, name); b = banner(pg)
                    check(f"{label}: example '{name}' runs ({want or 'a result'})", want in b and "Traceback" not in b, (f"{t:.1f}s", b[:120]))
                except Exception as e:  # noqa: BLE001
                    check(f"{label}: example '{name}' runs", False, repr(e)[:200])
            if not full:
                pg.close(); browser.close(); return
            # a path with spaces, a OneDrive-style folder name
            od = OUT / "OneDrive - Plumat Lab" / "my sequences (2026)"; od.mkdir(parents=True, exist_ok=True)
            shutil.copy2(EX / "clean" / "protein.aln.fasta", od / "protein aln.fasta"); shutil.copy2(EX / "clean" / "dna.fasta", od / "dna seqs.fasta")
            pg.close(); pg = new_page(ctx)
            click_button(pg, "Use file paths (large files)"); pg.wait_for_timeout(300)
            box = pg.locator("#coal-protein-path textarea, #coal-protein-path input").first; box.fill(str(od / "protein aln.fasta")); box.press("Tab")
            box = pg.locator("#coal-dna-path textarea, #coal-dna-path input").first; box.fill(str(od / "dna seqs.fasta")); box.press("Tab"); pg.wait_for_timeout(400)
            pg.locator("#coal-run").click(); wait_results(pg); b = banner(pg)
            check(f"{label}: path boxes with spaces and a 'OneDrive - ...' folder: 8 genes match", "8" in b and "match" in b.lower(), b[:200])
            # a long path (> 260 characters)
            longdir = OUT
            while len(str(longdir)) < 270:
                longdir = longdir / ("long_folder_name_" + "x" * 40)
            try:
                longdir.mkdir(parents=True, exist_ok=True)
                shutil.copy2(EX / "clean" / "protein.aln.fasta", longdir / "p.fasta"); shutil.copy2(EX / "clean" / "dna.fasta", longdir / "d.fasta")
                pg.close(); pg = new_page(ctx)
                click_button(pg, "Use file paths (large files)"); pg.wait_for_timeout(300)
                box = pg.locator("#coal-protein-path textarea, #coal-protein-path input").first; box.fill(str(longdir / "p.fasta")); box.press("Tab")
                box = pg.locator("#coal-dna-path textarea, #coal-dna-path input").first; box.fill(str(longdir / "d.fasta")); box.press("Tab"); pg.wait_for_timeout(400)
                pg.locator("#coal-run").click(); wait_results(pg); b = banner(pg); texts["long_path_banner"] = b[:300]
                check(f"{label}: a {len(str(longdir / 'p.fasta'))}-character path works or gives a plain message", "Traceback" not in b and len(b) > 5, b[:200])
                check(f"{label}: the long path aligns (8 genes)", "8" in b and "match" in b.lower(), b[:200])
            except Exception as e:  # noqa: BLE001
                check(f"{label}: long path", False, repr(e)[:200])
            # >= 1 MB through the path boxes: the bundled coal-cpp
            pg.close(); pg = new_page(ctx)
            click_button(pg, "Use file paths (large files)"); pg.wait_for_timeout(300)
            box = pg.locator("#coal-protein-path textarea, #coal-protein-path input").first; box.fill(str(BIG / "pep.fasta")); box.press("Tab")
            box = pg.locator("#coal-dna-path textarea, #coal-dna-path input").first; box.fill(str(BIG / "nuc.fasta")); box.press("Tab"); pg.wait_for_timeout(400)
            t0 = time.monotonic(); pg.locator("#coal-run").click(); wait_results(pg); t = time.monotonic() - t0; b = banner(pg)
            f = download(pg, "Download report")
            rep = json.loads(f.read_text(encoding="utf-8"))
            check(f"{label}: 10,000 x 500 through the path boxes runs in the bundled coal-cpp", rep.get("engine") == "coal-cpp" and "10" in b, (f"{t:.1f}s", rep.get("engine"), b[:100]))
            texts[f"{label}_big_command"] = rep.get("command")
            pg.screenshot(path=str(OUT / f"{label}_big.png"), full_page=False)
            pg.close()
            ctx3 = browser.new_context(viewport={"width": 1366, "height": 768}, color_scheme="dark")
            pgd = ctx3.new_page(); pgd.goto(url + "?__theme=dark", wait_until="domcontentloaded", timeout=60000); pgd.wait_for_selector("#coal-run", timeout=60000); pgd.wait_for_timeout(800)
            bg = pgd.evaluate("getComputedStyle(document.body).backgroundColor")
            run_example(pgd, "One broken gene"); pgd.screenshot(path=str(OUT / f"{label}_dark.png"), full_page=False)
            m = re.findall(r"\d+", bg); dark = m and sum(int(x) for x in m[:3]) < 200
            check(f"{label}: dark theme renders (dark body background)", bool(dark), bg)
            ctx3.close(); browser.close()
    except Exception as e:  # noqa: BLE001
        check(f"{label}: local app drive", False, repr(e)[:300])
    finally:
        try:
            p.stdin.close(); p.wait(timeout=20)
        except Exception:  # noqa: BLE001
            p.kill()
        check(f"{label}: --serve exits when its stdin closes", p.returncode == 0, p.returncode)


def first_launch():
    log = DATA / "logs" / "coal_desktop.log"
    before = log.read_text(encoding="utf-8", errors="replace") if log.exists() else ""
    t0 = time.monotonic()
    p = subprocess.Popen([str(APP)])
    ready = None; crashed = None
    while time.monotonic() - t0 < 90:
        if p.poll() is not None:
            crashed = p.returncode; break
        txt = log.read_text(encoding="utf-8", errors="replace") if log.exists() else ""
        new = txt[len(before):]
        if "ready at" in new:
            ready = time.monotonic() - t0; break
        time.sleep(0.5)
    texts["window_launch_s"] = ready; texts["window_log_tail"] = (log.read_text(encoding="utf-8", errors="replace")[len(before):])[-1500:] if log.exists() else ""
    check("COAL.exe window: the local app is ready within 90 s (log says 'ready at')", ready is not None and crashed is None, (ready, crashed, texts["window_log_tail"][-300:]))
    time.sleep(3)
    rc, so, se = run(["tasklist", "/FI", "IMAGENAME eq COAL.exe", "/FO", "CSV"], timeout=30)
    n = so.count("COAL.exe")
    check("COAL.exe window: the window and its local-app child are running (2 COAL.exe processes)", n >= 2, so[:300])
    run(["taskkill", "/PID", str(p.pid), "/T", "/F"], timeout=30)
    time.sleep(2)
    rc, so, se = run(["tasklist", "/FI", "IMAGENAME eq COAL.exe", "/FO", "CSV"], timeout=30)
    check("COAL.exe window: killing the window leaves no COAL.exe behind", "COAL.exe" not in so, so[:300])


def portable():
    z = RELEASE / "COAL-0.5.0-windows-x64.zip"; d = OUT / "portable"
    if d.exists():
        shutil.rmtree(d)
    with zipfile.ZipFile(z) as zf:
        zf.extractall(d)
    exes = list(d.rglob("COAL.exe"))
    check("portable zip: holds COAL.exe", len(exes) == 1, [str(e.relative_to(d)) for e in exes])
    if not exes:
        return
    exe = exes[0]
    rc, so, se = run([str(exe), "--version"], timeout=60)
    check("portable COAL.exe --version", rc == 0 and "0.5.0" in so, (rc, so.strip(), se[-100:]))
    rc, so, se = run([str(exe), "-m", "coal", "--cpp-path"], timeout=60)
    check("portable COAL.exe -m coal --cpp-path names its own bundled coal-cpp.exe", rc == 0 and so.strip().lower().startswith(str(d).lower()), (rc, so.strip()))
    rc, so, se = run([str(exe), "-m", "coal", str(BIG / "pep.fasta"), str(BIG / "nuc.fasta"), "-o", str(OUT / "portable_big.fa"), "--report", str(OUT / "portable_big.json")], timeout=300)
    rep = json.loads((OUT / "portable_big.json").read_text(encoding="utf-8")) if (OUT / "portable_big.json").exists() else {}
    check("portable COAL.exe -m coal on 13.7 MB: exit 0, engine coal-cpp", rc == 0 and rep.get("engine") == "coal-cpp", (rc, rep.get("engine"), se[-200:]))
    serve_and_drive(exe, 7871, "portable", full=False)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rc, so, se = run([str(APP), "--version"], timeout=60); check("installed COAL.exe --version", rc == 0 and "0.5.0" in so, so.strip())
    rc, so, se = run([str(APP), "-m", "coal", "--version"], timeout=60); check("installed COAL.exe -m coal --version is coal 0.5.0", rc == 0 and so.strip() == "coal 0.5.0", so.strip())
    rc, so, se = run([str(APP), "-m", "coal", "--cpp-path"], timeout=60); texts["cpp_path"] = so.strip()
    check("installed COAL.exe -m coal --cpp-path names the bundled coal-cpp.exe", rc == 0 and so.strip().lower().startswith(str(APP.parent).lower()), so.strip())
    rc, so, se = run([so.strip(), "--version"], timeout=60); check("bundled coal-cpp.exe --version is coal 0.5.0 (C++)", rc == 0 and "0.5.0" in so and "C++" in so, so.strip())
    rc, so, se = run([str(APP), "-m", "coal", str(BIG / "pep.fasta"), str(BIG / "nuc.fasta"), "-o", str(OUT / "installed_big.fa"), "--report", str(OUT / "installed_big.json")], timeout=300)
    rep = json.loads((OUT / "installed_big.json").read_text(encoding="utf-8")) if (OUT / "installed_big.json").exists() else {}
    check("installed COAL.exe -m coal on 13.7 MB: exit 0, engine coal-cpp (the bundled one)", rc == 0 and rep.get("engine") == "coal-cpp", (rc, rep.get("engine"), se[-200:]))
    texts["installed_big_command"] = rep.get("command")
    rc, so, se = run([str(APP), "-m", "coal_app", "--check"], timeout=120); check("installed COAL.exe -m coal_app --check", rc == 0, (rc, (so + se)[-200:]))
    rc, so, se = run([str(APP), "-m", "other"], timeout=60); check("installed COAL.exe -m other: refused with a plain message, non-zero", rc not in (0, None) and "only" in (so + se), (rc, (so + se)[:200]))
    for step in (mcp, claude_config, lambda: serve_and_drive(APP, 7870, "installed"), first_launch, portable):
        try:
            step()
        except Exception as e:  # noqa: BLE001
            check(f"step {getattr(step, '__name__', 'serve')}", False, repr(e)[:300])
    (OUT / "phase8_texts.json").write_text(json.dumps(texts, indent=1, default=str), encoding="utf-8")
    (OUT / "phase8_results.json").write_text(json.dumps(CHECKS, indent=1), encoding="utf-8")
    bad = [c for c in CHECKS if not c["ok"]]
    print(f"\n{len(CHECKS)} checks, {len(bad)} finding candidates")


if __name__ == "__main__":
    main()
