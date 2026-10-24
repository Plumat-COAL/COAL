#!/usr/bin/env python3
"""Phase 7 follow-up: the steps the main driver got wrong, each on a fresh page: settings reset and the settings-file
line, bad uploads (a refused extension, then garbage with a .fasta name) and a good run after them, two real uploads,
the 100,000-sequence run through the path boxes with the default code, the alignment view on that run, the centring
box on a 30,000-column preview, and the row click on the broken example.

Notes: the "View alignment" accordion is open by default (ui.py: open=True), so it is never clicked here; each file
input disappears after an upload, so the DNA file (the second input) is uploaded first."""
from __future__ import annotations

import json
import re
import sys
import time
from pathlib import Path

from playwright.sync_api import sync_playwright

sys.path.insert(0, str(Path(__file__).resolve().parent))
from phase7 import URL, BIG, CHECKS, check, new_page, wait_results, banner, health, click_button, run_example, download, select_dropdown, dropdown_value, SHOTS  # noqa: E402

OUT = Path(sys.argv[1] if len(sys.argv) > 1 else "phase7_out")
WIDE = "/mnt/d/coal-stress/verify/bench-data/n100_L10000"   # 100 sequences x 10,000 codons (30,000 columns)


def upload(pg, pep, dna):
    pg.locator("input[type=file]").nth(1).set_input_files(str(dna)); pg.wait_for_timeout(1200)
    pg.locator("input[type=file]").nth(0).set_input_files(str(pep)); pg.wait_for_timeout(1200)


def view_text(pg):
    return pg.locator("#coal-view").inner_text()


def toasts(pg):
    return " | ".join(t.strip() for t in pg.locator(".toast-body, .toast-text, [class*=toast]").all_inner_texts() if t.strip())[:400]


def main():
    OUT.mkdir(parents=True, exist_ok=True); SHOTS.mkdir(exist_ok=True)
    texts = {}
    with sync_playwright() as p:
        browser = p.chromium.launch(); ctx = browser.new_context()
        # settings: reset to defaults, the settings-file line
        pg = new_page(ctx)
        pg.locator("#coal-settings-tab-button").click(); pg.wait_for_timeout(500)
        tab = pg.locator("#coal-settings-tab").inner_text()
        check("settings tab says where the settings are stored (settings.json path)", "settings.json" in tab, tab[:400])
        pg.locator("#coal-set-reset").click(); pg.wait_for_timeout(1000)
        st = pg.locator("#coal-settings-status").inner_text()
        check("reset to defaults confirms", "default" in st.lower(), st)
        pg2 = new_page(ctx)
        v = dropdown_value(pg2, "coal-code")
        check("after reset the Align tab's genetic code is 1 (Standard)", v.strip().startswith("1"), v)
        pg2.close(); pg.close()
        # a refused extension
        pg = new_page(ctx)
        garbage = OUT / "garbage.bin"; garbage.write_bytes(bytes(range(256)) * 10)
        dna = OUT / "dna.fasta"; dna.write_bytes(b">a\nATGAAA\n")
        try:
            pg.locator("input[type=file]").nth(1).set_input_files(str(dna)); pg.wait_for_timeout(1200)
            pg.locator("input[type=file]").nth(0).set_input_files(str(garbage)); pg.wait_for_timeout(2500)
            texts["bin_toasts"] = toasts(pg)
            texts["bin_body_hint"] = " ".join(re.findall(r"[^\n]*(?:invalid|not allowed|accepted|file type)[^\n]*", pg.inner_text("body"), re.I))[:300]
            pg.locator("#coal-run").click(); pg.wait_for_timeout(8000)
            b = banner(pg); texts["bin_upload_banner"] = b
            pg.screenshot(path=str(SHOTS / "11a_bin_upload.png"), full_page=False)
            check("uploading a .bin as the protein file: the page says the file was refused (a toast or a banner line), not just 'Please add both files'",
                  "add both files" not in b.lower() or bool(texts["bin_toasts"]) or bool(texts["bin_body_hint"]), (b[:200], texts["bin_toasts"], texts["bin_body_hint"]))
        except Exception as e:  # noqa: BLE001
            check("uploading a .bin as the protein file: the page says the file was refused", False, repr(e)[:300])
        pg.close()
        # garbage with a .fasta name, then a good run
        pg = new_page(ctx)
        garbage2 = OUT / "garbage.fasta"; garbage2.write_bytes(bytes(range(256)) * 10)
        try:
            upload(pg, garbage2, dna)
            pg.locator("#coal-run").click(); pg.wait_for_timeout(15000)
            b = banner(pg); texts["bad_upload_banner"] = b
            pg.screenshot(path=str(SHOTS / "11_bad_upload.png"), full_page=False)
            check("garbage .fasta upload: a plain error in the banner, no traceback, says what to do",
                  "Traceback" not in b and len(b) > 10 and "add both files" not in b.lower() and ("not" in b.lower() or "error" in b.lower() or "could" in b.lower()), b[:400])
        except Exception as e:  # noqa: BLE001
            check("garbage .fasta upload: a plain error in the banner", False, repr(e)[:300])
        try:
            run_example(pg, "Clean nuclear genes")
            check("after a bad upload a good run works (recoverable)", "8" in banner(pg), banner(pg)[:150])
        except Exception as e:  # noqa: BLE001
            check("after a bad upload a good run works (recoverable)", False, repr(e)[:200])
        pg.close()
        # two real uploads
        pg = new_page(ctx)
        try:
            pep = OUT / "pep_up.fasta"; pep.write_bytes(b">a\nMKV\n>b\nMKV\n")
            dna2 = OUT / "dna_up.fasta"; dna2.write_bytes(b">a\nATGAAAGTT\n>b\nATGAAAGTG\n")
            upload(pg, pep, dna2)
            pg.locator("#coal-run").click(); wait_results(pg)
            b = banner(pg); texts["two_uploads_banner"] = b
            check("uploading two small files and aligning works (2 genes)", "2" in b and "error" not in b.lower(), b[:200])
            f = download(pg, "Download codon alignment")
            check("the uploaded run's alignment downloads (2 records)", f.read_text().count(">") == 2, f.read_text()[:120])
        except Exception as e:  # noqa: BLE001
            check("uploading two small files and aligning works (2 genes)", False, repr(e)[:300])
        pg.close()
        # the large input through the path boxes, default code
        pg = new_page(ctx)
        click_button(pg, "Use file paths (large files)"); pg.wait_for_timeout(300)
        box = pg.locator("#coal-protein-path textarea, #coal-protein-path input").first
        box.fill(f"{BIG}/pep.fasta"); box.press("Tab")
        box = pg.locator("#coal-dna-path textarea, #coal-dna-path input").first
        box.fill(f"{BIG}/nuc.fasta"); box.press("Tab"); pg.wait_for_timeout(500)
        t0 = time.monotonic(); pg.locator("#coal-run").click()
        lag = []
        while time.monotonic() - t0 < 900:
            s = time.monotonic(); pg.evaluate("1+1"); lag.append(time.monotonic() - s)
            try:
                wait_results(pg, timeout=3000); break
            except Exception:
                pass
        t = time.monotonic() - t0
        b = banner(pg); texts["big_banner"] = b; texts["big_health"] = health(pg)[:3000]
        check("100,000 clean sequences through the path boxes: everything matched, time to result", "100" in b.replace(",", "") and ("match" in b.lower() or "OK" in b), (f"{t:.1f}s", b[:200]))
        check("browser responsive during the 100,000-sequence run", max(lag) < 2.0, f"max page-eval lag {max(lag):.2f}s")
        texts["big_time_s"] = round(t, 1)
        pg.wait_for_timeout(1500)
        v = view_text(pg); texts["big_view_head"] = v[:800]
        check("large run: the alignment view says which sequences it shows (a preview of 100 of 100,000)", bool(re.search(r"(first|showing|of 100,?000|not shown|hidden|only|preview)", v, re.I)), v[:300].replace("\n", " "))
        check("large run: the centre-on-column box is hidden (all 300 columns fit the 1,000-column preview)", pg.locator("input[type=number]").count() == 0, pg.locator("input[type=number]").count())
        pg.screenshot(path=str(SHOTS / "10_100k_view.png"), full_page=True)
        pg.close()
        # a 30,000-column input: the preview leaves columns out, so the centring box shows
        pg = new_page(ctx)
        try:
            click_button(pg, "Use file paths (large files)"); pg.wait_for_timeout(300)
            box = pg.locator("#coal-protein-path textarea, #coal-protein-path input").first
            box.fill(f"{WIDE}/pep.fasta"); box.press("Tab")
            box = pg.locator("#coal-dna-path textarea, #coal-dna-path input").first
            box.fill(f"{WIDE}/nuc.fasta"); box.press("Tab"); pg.wait_for_timeout(500)
            pg.locator("#coal-run").click(); wait_results(pg); pg.wait_for_timeout(1500)
            v0 = view_text(pg); texts["wide_view_head"] = v0[:800]
            n = pg.locator("input[type=number]").count()
            check("30,000-column run: the centre-on-column box is shown", n >= 1, (n, v0[:200].replace("\n", " ")))
            check("30,000-column run: the preview says which columns it shows", bool(re.search(r"columns?\s+\d", v0, re.I)), v0[:300].replace("\n", " "))
            if n:
                num = pg.locator("input[type=number]").first
                num.fill("15000"); pg.locator("#coal-centre-go").click(); pg.wait_for_timeout(2500)
                v1 = view_text(pg); texts["wide_view_centred"] = v1[:800]
                check("centring the preview on column 15,000 redraws the view around it", v1 != v0 and bool(re.search(r"1[45],?\d{3}", v1)), (len(v0), len(v1), v1[:200].replace("\n", " ")))
                num.fill("99999999"); pg.locator("#coal-centre-go").click(); pg.wait_for_timeout(2500)
                v2 = view_text(pg); texts["wide_view_out_of_range"] = v2[:400]
                check("centring on a column past the end: no traceback, the view still draws", "Traceback" not in v2 and len(v2) > 100, v2[:200].replace("\n", " "))
                f = download(pg, "Download alignment view", "#coal-view-dl")
                check("download alignment view (HTML) works", f.suffix in (".html", ".htm") and f.stat().st_size > 1000, (f.name, f.stat().st_size))
            pg.screenshot(path=str(SHOTS / "10b_wide_view.png"), full_page=True)
        except Exception as e:  # noqa: BLE001
            check("30,000-column run: centring box and redraw", False, repr(e)[:300])
        pg.close()
        # row click on the broken example: the preview jumps to the failed gene's column
        pg = new_page(ctx)
        try:
            run_example(pg, "One broken gene"); pg.wait_for_timeout(1000)
            v0 = view_text(pg); texts["broken_view_head"] = v0[:800]
            row = pg.locator("#coal-health-box").get_by_text(re.compile("OG0001_sp14")).first
            row.click(); pg.wait_for_timeout(2000)
            v = view_text(pg); texts["broken_view_after_click"] = v[:800]
            marks = pg.locator("#coal-view").evaluate("e => Array.from(e.querySelectorAll('[class*=mark], [class*=frame], [class*=bad], [class*=issue]')).map(x => x.className).slice(0, 8)")
            check("clicking the failed gene's row keeps/brings its break column (102) in the view", "102" in v, (v0 == v, v[:200].replace("\n", " "), marks))
            pg.screenshot(path=str(SHOTS / "12_broken_row_click.png"), full_page=True)
        except Exception as e:  # noqa: BLE001
            check("clicking the failed gene's row", False, repr(e)[:300])
        pg.close()
        browser.close()
    (OUT / "followup_texts.json").write_text(json.dumps(texts, indent=1), encoding="utf-8")
    (OUT / "followup_results.json").write_text(json.dumps(CHECKS, indent=1), encoding="utf-8")
    bad = [c for c in CHECKS if not c["ok"]]
    print(f"\n{len(CHECKS)} checks, {len(bad)} finding candidates")


if __name__ == "__main__":
    main()
