import re, sys, json
from pathlib import Path
sys.path.insert(0, r"D:\coal-stress\verify")
from phase7 import URL, new_page, run_example, banner, health
from playwright.sync_api import sync_playwright
out = {}
with sync_playwright() as p:
    b = p.chromium.launch(); ctx = b.new_context(); pg = new_page(ctx)
    run_example(pg, "One broken gene"); pg.wait_for_timeout(800)
    out["dl_html"] = pg.locator("#coal-view-dl").evaluate("e => e.outerHTML.slice(0, 600)") if pg.locator("#coal-view-dl").count() else "NO #coal-view-dl"
    out["dl_buttons"] = pg.locator("#coal-view button, #coal-view a").all_inner_texts()[:10]
    try:
        with pg.expect_download(timeout=20000) as dl:
            pg.locator("#coal-view-dl").locator("button, a").first.click()
        f = Path(dl.value.path()); out["dl"] = (dl.value.suggested_filename, f.stat().st_size, f.read_text(errors="replace")[:120])
    except Exception as e:
        out["dl"] = repr(e)[:300]
    # the health table row for the failed gene: what is it (link?), what happens on click
    row = pg.locator("#coal-health-box").get_by_text(re.compile("OG0001_sp14")).first
    out["row_html"] = row.evaluate("e => (e.closest('a') || e.closest('tr') || e).outerHTML.slice(0, 500)")
    out["links_in_health"] = pg.locator("#coal-health-box a").count()
    out["view_scroll_before"] = pg.locator("#coal-view").evaluate("e => { const s = e.querySelector('.coalv, [style*=overflow], table'); return s ? [s.scrollLeft, s.scrollTop, window.scrollY] : [null, window.scrollY]; }")
    link = pg.locator("#coal-health-box a").first if pg.locator("#coal-health-box a").count() else row
    link.click(); pg.wait_for_timeout(1500)
    out["hash_after"] = pg.evaluate("location.hash")
    out["view_scroll_after"] = pg.locator("#coal-view").evaluate("e => { const s = e.querySelector('.coalv, [style*=overflow], table'); return s ? [s.scrollLeft, s.scrollTop, window.scrollY] : [null, window.scrollY]; }")
    out["target_exists"] = pg.evaluate("!!(location.hash && document.querySelector(location.hash.replace(/([:.])/g, '\\$1')))") if out["hash_after"] else None
    out["active_marked"] = pg.locator("#coal-view").evaluate("e => Array.from(e.querySelectorAll('[class*=target], [class*=active], [class*=sel], :target')).map(x => x.className).slice(0,5)")
    pg.screenshot(path=r"D:\coal-stress\verify\p7\screenshots\12_broken_row_click.png", full_page=False)
    # the garbage upload: the health table and raw messages
    pg2 = new_page(ctx)
    pg2.locator("input[type=file]").nth(1).set_input_files(r"D:\coal-stress\verify\p7\dna.fasta"); pg2.wait_for_timeout(1200)
    pg2.locator("input[type=file]").nth(0).set_input_files(r"D:\coal-stress\verify\p7\garbage.fasta"); pg2.wait_for_timeout(1200)
    pg2.locator("#coal-run").click(); pg2.wait_for_timeout(12000)
    out["bad_banner"] = banner(pg2); out["bad_health"] = health(pg2)[:800]
    out["bad_raw"] = " | ".join(t[:300] for t in pg2.locator("textarea").evaluate_all("els => els.map(e => e.value)")[:2])
    out["bad_body_error"] = " ".join(re.findall(r"[^\n]*(?:error|could not|not a|invalid|unreadable|decode|FASTA|CLUSTAL)[^\n]*", pg2.inner_text("body")))[:600]
    b.close()
Path(r"D:\coal-stress\verify\p7\probe2.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
for k, v in out.items(): print(k, ":", repr(v)[:500])
