import sys, json
from pathlib import Path
sys.path.insert(0, r"D:\coal-stress\verify")
from phase7 import new_page, run_example
from playwright.sync_api import sync_playwright
with sync_playwright() as p:
    b = p.chromium.launch(); ctx = b.new_context(); pg = new_page(ctx)
    run_example(pg, "One broken gene"); pg.wait_for_timeout(800)
    with pg.expect_download(timeout=30000) as dl:
        pg.locator("#coal-view-dl").click()
    f = Path(dl.value.path()); t = f.read_bytes()
    print("download:", dl.value.suggested_filename, len(t), "bytes; html:", t[:15], "; has coalv:", b"coalv" in t, "; has 102:", b"102" in t)
    b.close()
