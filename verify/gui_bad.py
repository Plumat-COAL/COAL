import sys, time
sys.path.insert(0, r"D:\coal-stress\verify")
from playwright.sync_api import sync_playwright
from phase7 import URL, new_page, banner, health
from pathlib import Path
OUT = Path(r"D:\coal-stress\verify\p7")
with sync_playwright() as p:
    b = p.chromium.launch(); ctx = b.new_context(); pg = new_page(ctx)
    garbage = OUT / "garbage.bin"; garbage.write_bytes(bytes(range(256)) * 10)
    dna = OUT / "dna.fasta"; dna.write_bytes(b">a\nATGAAA\n")
    pg.locator("input[type=file]").nth(0).set_input_files(str(garbage))
    pg.locator("input[type=file]").nth(1).set_input_files(str(dna))
    pg.wait_for_timeout(1500); pg.locator("#coal-run").click()
    for i in range(6):
        pg.wait_for_timeout(3000)
        print(f"--- t={3*(i+1)}s banner:", banner(pg)[:300].encode("ascii", "replace").decode())
    print("health:", health(pg)[:300].encode("ascii", "replace").decode())
    body = pg.inner_text("body")
    i = body.find("Files used"); print("files used:", body[i:i+200].encode("ascii","replace").decode())
    pg.screenshot(path=str(OUT / "screenshots" / "11_bad_upload.png"), full_page=True)
    b.close()
