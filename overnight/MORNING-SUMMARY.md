# COAL: overnight summary (Sat 17 Oct, ~05:00)

## Use this branch
**`keanu_overnight/integration`** has everything: the app work, the Windows desktop app, the engine fixes, and three rounds of independent review fixes. CI is green on every job (Windows, macOS and Linux; Python 3.10 and 3.13; real Gradio on Linux and Windows; MCP; the explainer in Chromium).

The other `keanu_overnight/*` branches (`app`, `windows-app`, `engine-fixes`) are earlier stages, kept for comparison. `main` is untouched.

## What was built overnight
- **GUI** (Gradio, `coal-gui`):
  - plain-language health table, genetic-code hint with one-click re-run, ready-to-paste methods sentence;
  - alignment viewer that jumps to the problem codon;
  - local file paths for large files (read in place, never copied);
  - MCP on by default (`--no-mcp` turns it off).
- **MCP tools** (`codon_align`, `explain_report`), plus a stdio server and **SKILL.md** for coding agents.
- **Windows desktop app**:
  - a native window around the local app;
  - "Open with" and drag-and-drop load files by path;
  - a per-user installer script and build scripts;
  - a smoke-test checklist.
- **Examples** (clean, mito, messy_ids, one_broken_gene), the **explainer animation**, **PITCH.md** and **QA.md** (all numbers tied to the test output by a guard test).
- **Demos**: race against PAL2NAL (per gene family), the DIAMOND pipeline preset, a scripted MCP stage demo.
- **Hugging Face Space**: prepared, not deployed (`space/DEPLOY.md`).
- **A drop-in `pal2nal` command.**

## Engine fixes (your teammate's engine code, changed carefully and marked)
- piped input works;
- `-o` no longer replaces shortcuts or devices, and keeps permissions and hard links;
- unrelated protein and DNA are an error (more than half of the *checkable* codons mismatch);
- masked genes are handled;
- the partial-ID-match warning;
- 50× faster fallback on long genes with N;
- clean Ctrl-C / SIGTERM / SIGHUP / broken pipe;
- table 3 GTG start;
- files from old Notepad and old Macs read correctly;
- error messages give `file:line`.
- **The test harness now checks what the README claims** and prints its own counts: 165 checks, 147 byte-identical (141 also with identical warnings), 12 both-fail, 6 deliberate differences.

## How every change is marked (for merging with your teammate's work)
- `overnight/*/CHANGES.md`: one row per change, with ID, where, Before, After, Why, the test, and a merge note.
- `CHANGED(<ID>)` comments at each behaviour change.
- The ID at the start of every commit message.
- Rows whose merge note says **"teammate's engine/C++ must mirror"** need copying into the C++ port.

## Review
Three rounds of independent review plus a final Windows CI round. Findings were fixed test-first, and every test was checked to fail without its fix.

| Suite | Result |
|---|---|
| Engine | 54 tests |
| App | 177 tests (5 skipped: need a non-root user, network for a wheel build, or Gradio) |
| Desktop | 63+ tests |

Several hundred thousand random known-answer cases across the rounds showed 0 wrong outputs on the final engine.

## Decisions waiting for you
1. **Mismatch threshold:** more than 50% of checkable codons mismatching is an error. A deliberate difference from PAL2NAL. Keep 50%?
2. **Notices and `-q`:** two notices (partial ID match, no anchor) still print with `-q`. Keep?
3. **Closed output pipe:** exits 141. Or 0?
4. **Drop-in `pal2nal` command:** without `-output fasta` it errors until CLUSTAL output exists. Or fall back to FASTA with a warning?
5. **Code signing** before any public Windows release (Smart App Control blocks unsigned apps).

## Needs a laptop / a human
- **Real-world runs:**
  - the race against real PAL2NAL (before quoting any number);
  - DIAMOND with the real binary;
  - deploying the Space;
  - the MCP demo in a real Claude client.
- **Windows desktop:** build the installer (`desktop/BUILD.md`) and run `desktop/SMOKE_TEST.md` on a Windows laptop. CI covers the logic, not the real window, WebView2 or the installer.
- **Rehearse the pitch on the presenting laptop;** record the 70-second backup video.
- **Fill the `[Team: …]` placeholders in `pitch/QA.md`.**
- **Person A (C++):** run `A2-engine-missed-steps.md` after the C++ run. Its task 0 now says to **merge `keanu_overnight/integration`**, not redo the engine fixes, and to mirror the marked rows in C++.

## Housekeeping
The `/dev/stdout` file in the cloud workspace was broken by a reviewer (it doesn't affect the repo or your computer). Say the word to restore it.
