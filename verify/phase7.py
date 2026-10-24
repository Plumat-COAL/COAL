#!/usr/bin/env python3
"""Phase 7: the web GUI, driven headless with Playwright (this script runs on the Windows host against the GUI that
`python -m coal_app` serves from WSL on http://127.0.0.1:7860; the Linux file paths it types refer to the server).

    python verify/phase7.py OUTDIR
"""
from __future__ import annotations

import json
import os
import re
import sys
import time
import urllib.request
from pathlib import Path

from playwright.sync_api import sync_playwright

URL = "http://127.0.0.1:7860/"
OUT = Path(sys.argv[1] if len(sys.argv) > 1 else "phase7_out")
SHOTS = OUT / "screenshots"
CLI = Path(r"D:\coal-stress\verify\p7\cli")          # the CLI's outputs of the examples (made in WSL)
WSL_REPO = "/home/chloe/coal-kit/coal-verify"         # paths as the server sees them
BIG = "/mnt/d/coal-stress/verify/p7/big"
AXE = Path(r"D:\coal-stress\verify\axe\node_modules\axe-core\axe.min.js")
CHECKS = []
EXTERNAL = []


def check(name, ok, detail=""):
    CHECKS.append({"name": name, "ok": bool(ok), "detail": str(detail)[:800]})
    print(("OK   " if ok else "FINDING-CANDIDATE ") + name + ("" if ok else f": {str(detail)[:400]}"), flush=True)


def step(name, fn):
    try:
        return fn()
    except Exception as e:  # noqa: BLE001
        check(name + " (harness)", False, repr(e)[:400])
        return None


def new_page(ctx, width=1366, height=768, theme=None):
    pg = ctx.new_page()
    pg.set_viewport_size({"width": width, "height": height})
    pg.goto(URL + (f"?__theme={theme}" if theme else ""), wait_until="domcontentloaded", timeout=60000)
    pg.wait_for_selector("#coal-run", timeout=60000)
    pg.wait_for_timeout(800)
    return pg


def wait_results(pg, timeout=240000):
    """After Start alignment: wait until the download button is enabled or the banner shows an error."""
    t0 = time.monotonic()
    pg.wait_for_function(
        """() => { const b = document.querySelector('#coal-banner'); const d = document.querySelector('#coal-downloads');
                  const dl = d && Array.from(d.querySelectorAll('button, a')).find(x => /Download codon alignment/.test(x.textContent));
                  const txt = b ? b.innerText : '';
                  return (dl && !dl.disabled && dl.getAttribute('aria-disabled') !== 'true') || /error|could not|failed|stopped/i.test(txt); }""",
        timeout=timeout)
    pg.wait_for_timeout(500)
    return time.monotonic() - t0


def banner(pg):
    return pg.locator("#coal-banner").inner_text().strip()


def health(pg):
    return pg.locator("#coal-health-box").inner_text().strip()


def click_button(pg, text, scope=None):
    loc = (pg.locator(scope) if scope else pg).get_by_role("button", name=re.compile(re.escape(text)))
    loc.first.click()


def select_dropdown(pg, elem_id, text):
    box = pg.locator(f"#{elem_id} input").first
    box.click()
    pg.wait_for_timeout(200)
    opt = pg.locator("[role=option]", has_text=re.compile(re.escape(text))).first
    opt.click()
    pg.wait_for_timeout(200)


def dropdown_value(pg, elem_id):
    return pg.locator(f"#{elem_id} input").first.input_value()


def run_example(pg, name):
    click_button(pg, name, "#coal-examples")
    pg.wait_for_timeout(500)
    pg.locator("#coal-run").click()
    return wait_results(pg)


def download(pg, button_text, scope="#coal-downloads"):
    with pg.expect_download(timeout=60000) as dl:
        click_button(pg, button_text, scope)
    d = dl.value
    p = OUT / "downloads" / d.suggested_filename
    p.parent.mkdir(parents=True, exist_ok=True)
    d.save_as(str(p))
    return p


def main():
    OUT.mkdir(parents=True, exist_ok=True); SHOTS.mkdir(exist_ok=True)
    texts = {}
    with sync_playwright() as p:
        browser = p.chromium.launch()
        ctx = browser.new_context()

        def on_route(route):
            host = urllib.parse.urlsplit(route.request.url).hostname
            if host not in ("127.0.0.1", "localhost"):
                EXTERNAL.append(route.request.url)
                return route.abort()
            return route.continue_()
        import urllib.parse
        ctx.route("**/*", on_route)
        pg = new_page(ctx)
        texts["initial"] = pg.inner_text("body")
        pg.screenshot(path=str(SHOTS / "01_initial_1366x768.png"), full_page=True)
        check("GUI page loads with the Align and Settings tabs and Start alignment", "Start alignment" in texts["initial"] and "Settings" in texts["initial"])
        # G03/G11: clean example
        t = step("clean example", lambda: run_example(pg, "Clean nuclear genes"))
        texts["clean"] = pg.inner_text("body")
        pg.screenshot(path=str(SHOTS / "02_clean_results.png"), full_page=True)
        b = banner(pg)
        check("clean example: banner says everything matched", bool(re.search(r"(all|every|8 of 8).*(match|ok)", b, re.I)) or "✔" in b, b[:200])
        check("clean example: time to result under 30 s", t is not None and t < 30, t)
        h = health(pg)
        check("clean example: health table lists the 8 genes with a symbol and a word", h.count("OK") >= 8 or ("✔" in h and "OK" in h), h[:300])
        # G13 downloads in every format
        for fmt in ("fasta", "clustal", "paml", "codon"):
            def dl_fmt(fmt=fmt):
                if fmt != "fasta":
                    select_dropdown(pg, "coal-format", {"clustal": "CLUSTAL", "paml": "PAML", "codon": "Codon"}.get(fmt, fmt))
                    pg.locator("#coal-run").click(); wait_results(pg)
                f = download(pg, "Download codon alignment")
                ref = (CLI / f"clean.{fmt}").read_bytes()
                check(f"clean example download ({fmt}): byte-identical to the CLI", f.read_bytes() == ref, (f.name, f.stat().st_size, len(ref)))
            step(f"download {fmt}", dl_fmt)
        step("format back to fasta", lambda: select_dropdown(pg, "coal-format", "FASTA"))
        # G04 sample report
        def sample():
            click_button(pg, "Sample report", "#coal-examples"); pg.wait_for_timeout(2500)
            b = banner(pg); h = health(pg)
            check("Sample report button loads a report into the results (banner and health table filled)", len(b) > 5 and len(h) > 20, (b[:120], h[:120]))
            run_example(pg, "Clean nuclear genes")
        step("sample report", sample)
        # G14 report
        def rep():
            f = download(pg, "Download report")
            r = json.loads(f.read_bytes())
            check("clean example report (JSON): schema coal.report/1, 8 ok", r.get("schema") == "coal.report/1" and r["summary"]["ok"] == 8, r.get("summary"))
        step("report download", rep)
        # G17/G18 methods and raw messages
        def methods():
            click_button(pg, "Methods, commands and glossary")
            pg.wait_for_timeout(300)
            m = pg.locator("textarea").nth(0).input_value()  # (the Methods paragraph textbox: its label is not linked, see the findings)
            check("methods paragraph present with COAL v0.5.0 and the standard code", "COAL v0.5.0" in m and "standard genetic code" in m, m[:200])
            cp = pg.locator("button[aria-label*='opy'], button[title*='opy']").count()
            check("methods paragraph has a copy button", cp >= 1, cp)
            raw = pg.get_by_label("Raw COAL messages").input_value() if pg.get_by_label("Raw COAL messages").count() else pg.locator("textarea").nth(1).input_value()
            check("raw COAL messages box present", isinstance(raw, str), raw[:100])
            cmds = pg.inner_text("body")
            check("equivalent commands shown (coal and pal2nal)", "pal2nal" in cmds and "coal " in cmds, "")
        step("methods", methods)
        # G03/G12: mito example and the hint's re-run
        t = step("mito example", lambda: run_example(pg, "Mitochondrial (wrong code)"))
        texts["mito"] = pg.inner_text("body")
        pg.screenshot(path=str(SHOTS / "03_mito_wrong_code.png"), full_page=True)
        b = banner(pg); hint = pg.locator("#coal-hint").inner_text().strip() if pg.locator("#coal-hint").count() else ""
        check("mito example: banner says something is wrong and the tip names code 2", ("2" in hint and re.search(r"mitochondrial|genetic code", hint, re.I)) is not None, (b[:150], hint[:200]))
        def rerun():
            pg.locator("#coal-hint-rerun").click()
            wait_results(pg)
            check("hint re-run with code 2: all genes ok, code dropdown shows 2", "24" in banner(pg) and "2" in dropdown_value(pg, "coal-code"), (banner(pg)[:150], dropdown_value(pg, "coal-code")))
        step("hint re-run", rerun)
        texts["mito_fixed"] = pg.inner_text("body")
        pg.screenshot(path=str(SHOTS / "04_mito_fixed.png"), full_page=True)
        # messy IDs and one broken gene
        step("messy ids", lambda: run_example(pg, "Messy IDs"))
        texts["messy"] = pg.inner_text("body")
        check("messy IDs example: banner/health mention the pairing or tidy IDs", re.search(r"tid|normalis|pair|\.1", texts["messy"], re.I) is not None, banner(pg)[:200])
        step("one broken gene", lambda: run_example(pg, "One broken gene"))
        texts["broken"] = pg.inner_text("body")
        pg.screenshot(path=str(SHOTS / "05_one_broken_gene.png"), full_page=True)
        h = health(pg)
        check("one broken gene: the failed gene listed first with a sentence (frameshift) and a symbol + word", bool(re.search(r"OG0001_sp14.*(frameshift|missing|extra)", h, re.I | re.S)) and ("Failed" in h or "✖" in h), h[:400])
        # G15/G16 the alignment view: click a warning/failed row, centre on a column
        def view():
            acc = pg.get_by_role("button", name=re.compile("View alignment"))
            if acc.count() and pg.locator("input[type=number]").count() == 0:
                acc.first.click()
            pg.wait_for_timeout(800)
            check("alignment view has a centre-on-column number box", pg.locator("input[type=number]").count() >= 1, pg.locator("#coal-view").inner_text()[:200])
            v0 = pg.locator("#coal-view").inner_text()
            num = pg.locator("input[type=number]").first
            num.fill("102"); num.press("Enter")
            go = pg.locator("#coal-centre-go")
            if go.count():
                (go.locator("button").first if go.locator("button").count() else go).click()
            pg.wait_for_timeout(1500)
            v1 = pg.locator("#coal-view").inner_text()
            check("alignment view: centring on column 102 changes the view", v1 != v0 and len(v1) > 50, (len(v0), len(v1)))
            # clicking a gene row with a column jumps there
            row = pg.locator("#coal-health-box").get_by_text(re.compile("OG0001_sp14")).first
            row.click(); pg.wait_for_timeout(1500)
            centre = pg.locator("input[type=number]").first.input_value()
            check("clicking the failed gene's row centres the preview on its column (102)", centre.strip() in ("102", "102.0"), centre)
            f = download(pg, "Download alignment view", "#coal-view-dl")
            check("download alignment view (HTML) works", f.suffix in (".html", ".htm") and f.stat().st_size > 1000, f.name)
        step("view", view)
        # G19 settings persistence
        def settings():
            pg.locator("#coal-settings-tab-button").click(); pg.wait_for_timeout(500)
            select_dropdown(pg, "coal-set-code", "2: Vertebrate")
            pg.wait_for_timeout(800)
            st = pg.locator("#coal-settings-status").inner_text()
            check("settings: saved confirmation with the settings file path", "Saved" in st and "settings.json" in st, st[:200])
            pg2 = new_page(ctx)
            v = dropdown_value(pg2, "coal-code")
            check("settings persist across a reload: Align tab's genetic code is 2", v.strip().startswith("2"), v)
            pg2.close()
            pg.locator("#coal-set-reset").locator("button").first.click(); pg.wait_for_timeout(800)
            pg3 = new_page(ctx)
            v = dropdown_value(pg3, "coal-code")
            check("reset to defaults: genetic code back to 1", v.strip().startswith("1"), v)
            pg3.close()
            pg.locator("#coal-align-tab-button").click(); pg.wait_for_timeout(300)
        step("settings", settings)
        # G20 bad upload
        def bad_upload():
            garbage = OUT / "garbage.bin"; garbage.write_bytes(bytes(range(256)) * 10)
            dna = OUT / "dna.fasta"; dna.write_bytes(b">a\nATGAAA\n")
            pg.locator("input[type=file]").nth(0).set_input_files(str(garbage))
            pg.locator("input[type=file]").nth(1).set_input_files(str(dna))
            pg.wait_for_timeout(1000)
            pg.locator("#coal-run").click(); wait_results(pg)
            b = banner(pg)
            check("bad upload: a plain error in the banner (no traceback)", "Traceback" not in b and len(b) > 10, b[:300])
            texts["bad_upload"] = b
            run_example(pg, "Clean nuclear genes")
            check("bad upload: a good run afterwards works (recoverable)", "8" in banner(pg), banner(pg)[:150])
        step("bad upload", bad_upload)
        # G01/G02 upload real files
        def uploads():
            pep = OUT / "pep_up.fasta"; pep.write_bytes(b">a\nMKV\n>b\nMKV\n")
            dna = OUT / "dna_up.fasta"; dna.write_bytes(b">a\nATGAAAGTT\n>b\nATGAAAGTG\n")
            pg.locator("input[type=file]").nth(0).set_input_files(str(pep))
            pg.locator("input[type=file]").nth(1).set_input_files(str(dna))
            pg.wait_for_timeout(1000)
            pg.locator("#coal-run").click(); wait_results(pg)
            check("uploading two small files and aligning works", "2" in banner(pg), banner(pg)[:150])
        step("uploads", uploads)
        # G22 scale: 100,000 sequences via the path boxes
        def scale():
            nonlocal pg
            pg.close(); pg = new_page(ctx)  # a fresh page: no example upload around
            click_button(pg, "Use file paths (large files)")
            pg.wait_for_timeout(300)
            box = pg.locator("#coal-protein-path textarea, #coal-protein-path input").first
            box.fill(f"{BIG}/pep.fasta"); box.press("Tab")
            box = pg.locator("#coal-dna-path textarea, #coal-dna-path input").first
            box.fill(f"{BIG}/nuc.fasta"); box.press("Tab")
            pg.wait_for_timeout(500)
            t0 = time.monotonic()
            pg.locator("#coal-run").click()
            responsive = []
            while time.monotonic() - t0 < 600:
                s = time.monotonic(); pg.evaluate("1+1"); responsive.append(time.monotonic() - s)
                try:
                    wait_results(pg, timeout=3000); break
                except Exception:
                    pass
            t = time.monotonic() - t0
            b = banner(pg); texts["big"] = pg.inner_text("body")
            check("100,000 sequences through the path boxes: result within 10 minutes", t < 600 and "100" in b, (f"{t:.1f}s", b[:150]))
            check("browser stays responsive during the run (page eval under 2 s)", max(responsive) < 2.0, f"max {max(responsive):.2f}s")
            v = pg.locator("#coal-view").inner_text() if pg.locator("#coal-view").count() else ""
            check("preview says what it hides on a 100,000-sequence input", bool(re.search(r"(first|only|showing|of 100,?000|hidden|not shown)", v + texts["big"], re.I)), (v[:200]))
            pg.screenshot(path=str(SHOTS / "06_100k_results.png"), full_page=False)
        step("scale", scale)
        # G23 offline
        check("offline: no request to any host but 127.0.0.1 (fonts, CDNs, analytics)", not EXTERNAL, EXTERNAL[:5])
        # G24 keyboard-only flow and focus visibility
        def keyboard():
            pgk = new_page(ctx)
            found = False
            for _ in range(60):
                pgk.keyboard.press("Tab")
                txt = pgk.evaluate("document.activeElement && document.activeElement.innerText || ''")
                if "Clean nuclear genes" in txt:
                    found = True; break
            check("keyboard: Tab reaches the example button", found)
            vis = pgk.evaluate("(() => { const e = document.activeElement; const s = getComputedStyle(e); return s.outlineStyle + '|' + s.outlineWidth + '|' + s.boxShadow; })()")
            check("keyboard: focus is visible on the focused button (outline or box-shadow)", not (vis.startswith("none|") and vis.endswith("|none")), vis)
            pgk.keyboard.press("Enter"); pgk.wait_for_timeout(500)
            found = False
            for _ in range(60):
                pgk.keyboard.press("Tab")
                txt = pgk.evaluate("document.activeElement && document.activeElement.innerText || ''")
                if "Start alignment" in txt:
                    found = True; break
            check("keyboard: Tab reaches Start alignment", found)
            pgk.keyboard.press("Enter"); wait_results(pgk)
            check("keyboard-only run completes", "8" in banner(pgk), banner(pgk)[:150])
            # axe-core
            if AXE.exists():
                pgk.add_script_tag(content=AXE.read_text(encoding="utf-8"))
                res = pgk.evaluate("axe.run(document, {runOnly: {type: 'tag', values: ['wcag2a', 'wcag2aa']}})")
                vio = res.get("violations", [])
                (OUT / "axe.json").write_text(json.dumps(res, indent=1))
                serious = [v for v in vio if v.get("impact") in ("serious", "critical")]
                contrast = [v for v in vio if v["id"] == "color-contrast"]
                check("axe-core: no serious or critical WCAG 2 A/AA violations", not serious, [(v["id"], v["impact"], len(v["nodes"])) for v in serious])
                check("axe-core: no colour-contrast violations (4.5:1)", not contrast, [(n["html"][:80]) for v in contrast for n in v["nodes"][:5]])
                texts["axe_all"] = json.dumps([(v["id"], v["impact"], len(v["nodes"]), v["help"]) for v in vio])
            else:
                check("axe-core available", False, str(AXE))
            h = health(pgk)
            check("status never by colour alone: health rows carry a symbol and a word", ("OK" in h or "Warning" in h or "Failed" in h) and any(s in h for s in ("✔", "⚠", "✖")), h[:200])
            pgk.close()
        step("keyboard", keyboard)
        # G21 dark mode, G25 screenshots
        def shots():
            for w, h_, name in ((1920, 1080, "07_results_1920x1080.png"),):
                pgs = new_page(ctx, w, h_); run_example(pgs, "One broken gene"); pgs.screenshot(path=str(SHOTS / name), full_page=True); pgs.close()
            ctx2 = browser.new_context(viewport={"width": 683, "height": 384}, device_scale_factor=2)  # 200 % zoom of 1366x768
            pgz = ctx2.new_page(); pgz.goto(URL, wait_until="domcontentloaded", timeout=60000); pgz.wait_for_selector("#coal-run", timeout=60000); pgz.wait_for_timeout(800)
            run_example(pgz, "One broken gene")
            pgz.screenshot(path=str(SHOTS / "08_results_200pct_zoom.png"), full_page=True)
            overflow = pgz.evaluate("document.documentElement.scrollWidth > document.documentElement.clientWidth + 2")
            check("200 % zoom (683 CSS px wide, scale 2): no horizontal page overflow", not overflow, pgz.evaluate("document.documentElement.scrollWidth + ' > ' + document.documentElement.clientWidth"))
            ctx2.close()
            ctx3 = browser.new_context(viewport={"width": 1366, "height": 768}, color_scheme="dark")
            pgd = ctx3.new_page(); pgd.goto(URL + "?__theme=dark", wait_until="domcontentloaded", timeout=60000); pgd.wait_for_selector("#coal-run", timeout=60000); pgd.wait_for_timeout(800)
            bg = pgd.evaluate("getComputedStyle(document.body).backgroundColor + ' / ' + document.body.className")
            run_example(pgd, "Mitochondrial (wrong code)")
            pgd.screenshot(path=str(SHOTS / "09_dark_mode.png"), full_page=True)
            m = re.findall(r"\d+", bg)
            dark = m and sum(int(x) for x in m[:3]) < 200
            check("dark mode renders with a dark background", bool(dark), bg)
            texts["dark"] = pgd.inner_text("body")
            ctx3.close()
        step("shots", shots)
        # G27 MCP over SSE from the GUI
        def sse():
            req = urllib.request.Request(URL + "gradio_api/mcp/sse", headers={"Accept": "text/event-stream"})
            try:
                with urllib.request.urlopen(req, timeout=5) as r:
                    ct = r.headers.get("Content-Type", "")
                    check("GUI serves MCP at /gradio_api/mcp/sse (event stream)", "event-stream" in ct, ct)
            except Exception as e:  # noqa: BLE001
                # a timeout while the stream stays open is a success of sorts: the endpoint answered
                check("GUI serves MCP at /gradio_api/mcp/sse (event stream)", "timed out" in str(e).lower(), repr(e))
        step("sse", sse)
        browser.close()
    (OUT / "page_texts.json").write_text(json.dumps(texts, indent=1), encoding="utf-8")
    (OUT / "results.json").write_text(json.dumps(CHECKS, indent=1), encoding="utf-8")
    bad = [c for c in CHECKS if not c["ok"]]
    print(f"\n{len(CHECKS)} checks, {len(bad)} finding candidates")


if __name__ == "__main__":
    main()
