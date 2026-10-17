"""Gradio web GUI for COAL, which also serves the MCP tools (``launch(mcp_server=True)``).

Written for Gradio 5.x (``pip install -e .[gui]``). All logic lives in ``core``
(tested without Gradio); this file only lays out components and wires events.
Offline: no CDN fonts (system fonts only) and analytics disabled.
"""
from __future__ import annotations

import os

os.environ.setdefault("GRADIO_ANALYTICS_ENABLED", "False")

import gradio as gr  # noqa: E402

from . import core, gui_html, tools, viewer  # noqa: E402


# CHANGED(GUI-3): was Gradio's default theme plus a few table rules; now a COAL theme (colours of the COAL logo and
# site: charcoal #2A2A2A, slate #4A4744, ember #FF6A13 under dark text #1E1E1E, amber #FFC24B, page #121110 in
# dark), system fonts only (no web fonts are requested), for light and for Gradio's dark theme (".dark" on the
# body). Text contrast is at least 4.5:1 in both (figures in overnight/integration/CHANGES.md GUI-3). Every status
# keeps an icon and a word, never colour alone.
SANS = ["Seravek", "Gill Sans Nova", "Ubuntu", "Calibri", "Segoe UI", "DejaVu Sans", "system-ui", "sans-serif"]
MONO = ["ui-monospace", "Cascadia Mono", "SF Mono", "Menlo", "Consolas", "DejaVu Sans Mono", "monospace"]

CSS = """
:root, .light { --coal-page:#F6F5F3; --coal-surface:#FFFFFF; --coal-ink:#1E1E1E; --coal-text:#2A2A2A;
  --coal-muted:#5E5853; --coal-line:#DEDAD6; --coal-ok:#4A4744; --coal-warn:#8A5D00; --coal-warn-fill:#FFC24B;
  --coal-fail:#B84400; --coal-fail-fill:#FF6A13; --coal-info:#0B5CAD; --coal-focus:#B84400;
  --coal-tile-ok:#4A4744; --coal-tile-line:#8C8782; --coal-m1:#2A2A2A; --coal-m2:#4A4744; }
.dark { --coal-page:#121110; --coal-surface:#1C1A19; --coal-ink:#F1EEEA; --coal-text:#E4E0DB;
  --coal-muted:#B3ADA6; --coal-line:#3A3633; --coal-ok:#B5AFA8; --coal-warn:#FFC24B; --coal-warn-fill:#FFC24B;
  --coal-fail:#FF8A45; --coal-fail-fill:#FF6A13; --coal-info:#8CC4FF; --coal-focus:#FFC24B;
  --coal-tile-ok:#6A655F; --coal-tile-line:#8C8782; --coal-m1:#8C8782; --coal-m2:#B5AFA8; }
.gradio-container { max-width:1080px !important; margin:0 auto !important; }
/* header */
.coal-head { display:flex; gap:16px; align-items:flex-start; padding:8px 0 4px; }
.coal-head .coal-mark { flex:none; margin-top:6px; }
.coal-mark .m1 { fill:var(--coal-m1); } .coal-mark .m2 { fill:var(--coal-m2); }
.coal-mark .m3 { fill:#FF6A13; } .coal-mark .m4 { fill:#FFC24B; }
.coal-head h1 { margin:0; font-size:1.7rem; line-height:1.2; color:var(--coal-ink); letter-spacing:0.01em; }
.coal-head .coal-head-sub { font-weight:400; color:var(--coal-muted); }
.coal-head p { margin:0.25em 0 0; max-width:44rem; color:var(--coal-text); line-height:1.5; }
.coal-head p.coal-head-note { font-size:0.9em; color:var(--coal-muted); }
/* example row: a compact line of small buttons, not a grid */
#coal-examples { flex-wrap:wrap !important; gap:8px !important; align-items:center; }
#coal-examples > * { flex:0 0 auto !important; min-width:0 !important; width:auto !important; }
#coal-examples .coal-ex-label p, #coal-examples .coal-ex-label { margin:0; color:var(--coal-muted); }
#coal-downloads { flex-wrap:wrap !important; gap:8px !important; }
#coal-downloads > * { flex:0 0 auto !important; min-width:0 !important; width:auto !important; }
#coal-run { font-size:1.1rem !important; min-height:52px; font-weight:700; }
/* result banner */
.coal-banner { border:1px solid var(--coal-line); border-left:6px solid var(--coal-muted); border-radius:8px;
  padding:14px 18px; background:var(--coal-surface); color:var(--coal-text); }
.coal-banner.coal-ok { border-left-color:var(--coal-ok); }
.coal-banner.coal-warning { border-left-color:var(--coal-warn-fill); }
.coal-banner.coal-failed { border-left-color:var(--coal-fail-fill); }
.coal-banner.coal-info { border-left-color:var(--coal-info); }
.coal-banner.coal-empty { border-left-color:var(--coal-line); border-style:dashed; border-left-style:solid; }
.coal-banner p { margin:0; }
.coal-status { font-size:0.9em; color:var(--coal-muted); }
.coal-headline { font-size:1.45rem; font-weight:700; line-height:1.25; color:var(--coal-ink); margin:2px 0 6px !important; }
.coal-icon { display:inline-block; min-width:1.1em; }
.coal-ok .coal-status .coal-icon { color:var(--coal-ok); }
.coal-warning .coal-status .coal-icon, .coal-warning > .coal-headline .coal-icon { color:var(--coal-warn); }
.coal-failed .coal-status .coal-icon, .coal-failed > .coal-headline .coal-icon { color:var(--coal-fail); }
.coal-info .coal-icon { color:var(--coal-info); }
.coal-chips { display:flex; flex-wrap:wrap; gap:6px 14px; margin:0 0 8px; }
.coal-chips:empty { display:none; }
.coal-chip { font-weight:600; }
.coal-chip.coal-failed { color:var(--coal-fail); } .coal-chip.coal-warning { color:var(--coal-warn); }
.coal-chip.coal-unchecked { color:var(--coal-muted); }
.coal-strip { display:flex; flex-wrap:wrap; gap:3px; margin:4px 0 10px; max-width:100%; }
.coal-tile { width:14px; height:14px; border-radius:3px; box-sizing:border-box; font:700 10px/14px system-ui,
  sans-serif; text-align:center; color:#1E1E1E; overflow:hidden; }
.coal-tile.coal-ok { background:var(--coal-tile-ok); }
.coal-tile.coal-warning { background:var(--coal-warn-fill); }
.coal-tile.coal-failed { background:var(--coal-fail-fill); }
.coal-tile.coal-unchecked { border:1.5px dashed var(--coal-tile-line); color:var(--coal-muted); line-height:11px; }
.coal-detail { color:var(--coal-text); line-height:1.5; }
.coal-files { margin-top:4px !important; font-size:0.88em; color:var(--coal-muted); overflow-wrap:anywhere; }
#coal-hint { border-left:6px solid var(--coal-warn-fill); border-radius:8px; padding:0.6em 0.9em; }
/* health table */
.coal-health { width:100%; border-collapse:collapse; font-size:1rem; color:var(--coal-text); }
.coal-health caption { text-align:left; font-weight:600; padding:0.4em 0; color:var(--coal-ink); }
.coal-health th, .coal-health td { text-align:left; padding:0.5em 0.7em; border:0;
  border-bottom:1px solid var(--coal-line); vertical-align:top; background:transparent; }
.coal-health th { color:var(--coal-muted); font-weight:600; font-size:0.9em; }
.coal-health tr.coal-ok td:first-child { border-left:6px solid var(--coal-tile-ok); }
.coal-health tr.coal-warning td:first-child { border-left:6px solid var(--coal-warn-fill); }
.coal-health tr.coal-failed td:first-child { border-left:6px solid var(--coal-fail-fill); }
.coal-health tr.coal-unchecked td:first-child { border-left:6px dashed var(--coal-tile-line); }
.coal-health td:first-child { white-space:nowrap; }
.coal-health summary { cursor:pointer; color:var(--coal-text); }
.coal-health .coal-oklist { margin:0.5em 0 0; line-height:1.9; }
.coal-health .coal-oklist code { margin-right:0.4em; white-space:nowrap; }
.coal-health a { color:var(--coal-info); }
.coal-badge { display:inline-block; min-width:1.2em; font-weight:700; }
.coal-badge.coal-ok { color:var(--coal-ok); } .coal-badge.coal-warning { color:var(--coal-warn); }
.coal-badge.coal-failed { color:var(--coal-fail); } .coal-badge.coal-unchecked { color:var(--coal-muted); }
.coal-small { font-size:0.88em; color:var(--coal-muted); }
abbr[title] { text-decoration:underline dotted; cursor:help; }
button:focus-visible, input:focus-visible, select:focus-visible, textarea:focus-visible, a:focus-visible,
summary:focus-visible, [tabindex]:focus-visible {
  outline:3px solid var(--coal-focus) !important; outline-offset:2px; }
@media (max-width: 640px) {
  .coal-head { gap:12px; } .coal-head .coal-mark { width:32px; height:25px; }
  .coal-head h1 { font-size:1.4rem; } .coal-headline { font-size:1.2rem; }
  .coal-health th, .coal-health td { padding:0.45em 0.45em; }
  .coal-health td:first-child { white-space:normal; } }
@media (prefers-reduced-motion: reduce) { * { transition:none !important; animation:none !important; } }
"""


def _theme():
    """CHANGED(GUI-3): COAL theme. Fonts are plain strings (never gr.themes.GoogleFont): nothing is fetched."""
    stone = gr.themes.Color(c50="#F6F5F3", c100="#EEECE9", c200="#DEDAD6", c300="#C4BEB8", c400="#A39C95",
                            c500="#7D7670", c600="#5E5853", c700="#4A4744", c800="#2E2B28", c900="#1C1A19",
                            c950="#121110", name="coal_stone")
    ember = gr.themes.Color(c50="#FFF3EB", c100="#FFE2CF", c200="#FFC4A0", c300="#FFA26B", c400="#FF8A45",
                            c500="#FF6A13", c600="#E25500", c700="#B84400", c800="#8F3500", c900="#662600",
                            c950="#401800", name="coal_ember")
    theme = gr.themes.Base(primary_hue=ember, secondary_hue=stone, neutral_hue=stone, font=SANS, font_mono=MONO,
                           text_size=gr.themes.sizes.text_lg, radius_size=gr.themes.sizes.radius_md)
    return theme.set(
        body_background_fill="#F6F5F3", body_background_fill_dark="#121110",
        body_text_color="#2A2A2A", body_text_color_dark="#E4E0DB",
        body_text_color_subdued="#5E5853", body_text_color_subdued_dark="#B3ADA6",
        background_fill_primary="#FFFFFF", background_fill_primary_dark="#1C1A19",
        background_fill_secondary="#F6F5F3", background_fill_secondary_dark="#161413",
        block_background_fill="#FFFFFF", block_background_fill_dark="#1C1A19",
        block_border_color="#DEDAD6", block_border_color_dark="#3A3633",
        border_color_primary="#DEDAD6", border_color_primary_dark="#3A3633",
        block_label_text_color="#4A4744", block_label_text_color_dark="#B3ADA6",
        block_title_text_color="#2A2A2A", block_title_text_color_dark="#E4E0DB",
        block_info_text_color="#5E5853", block_info_text_color_dark="#B3ADA6",
        input_background_fill="#FFFFFF", input_background_fill_dark="#121110",
        input_border_color="#C4BEB8", input_border_color_dark="#4A4744",
        input_border_color_focus="#B84400", input_border_color_focus_dark="#FFC24B",
        button_primary_background_fill="#FF6A13", button_primary_background_fill_dark="#FF6A13",
        button_primary_background_fill_hover="#FF8A45", button_primary_background_fill_hover_dark="#FF8A45",
        button_primary_text_color="#1E1E1E", button_primary_text_color_dark="#1E1E1E",
        button_primary_border_color="#E25500", button_primary_border_color_dark="#FF6A13",
        button_secondary_background_fill="#FFFFFF", button_secondary_background_fill_dark="#2E2B28",
        button_secondary_background_fill_hover="#EEECE9", button_secondary_background_fill_hover_dark="#3A3633",
        button_secondary_text_color="#2A2A2A", button_secondary_text_color_dark="#F1EEEA",
        button_secondary_border_color="#A39C95", button_secondary_border_color_dark="#5A534D",
        checkbox_background_color_selected="#FF6A13", checkbox_background_color_selected_dark="#FF6A13",
        checkbox_border_color_selected="#B84400", checkbox_border_color_selected_dark="#FF6A13",
        link_text_color="#0B5CAD", link_text_color_dark="#8CC4FF",
        color_accent="#FF6A13", color_accent_soft="#FFE2CF", color_accent_soft_dark="#3A2A20",
    )


# CHANGED(GUI-1): the title and paragraph are now the short gui_html.header_html (mark, name, one sentence).
# CHANGED(GUI-5): was the flags (-nogap, -nomismatch) as glossary terms; now the plain names used on the options,
# with the flag after them, plus the other two options
GLOSSARY = """
- **Codon**: three DNA letters that code for one amino acid (for example ATG = M, methionine).
- **Protein alignment**: protein sequences lined up so that matching positions are in the same column (from MAFFT, MUSCLE, Clustal...).
- **Genetic code**: the dictionary from codons to amino acids. Almost all nuclear genes use code 1; mitochondria and some microbes use slightly different ones.
- **Mismatch**: a codon that doesn't code for the amino acid at that position. A few are often a sign of the wrong genetic code.
- **Frameshift**: a missing or extra DNA letter, so every later codon is read out of step.
- **Remove gap columns** (`-nogap`): remove columns containing gaps or stop codons (some dN/dS tools need this).
- **Remove mismatched columns** (`-nomismatch`): remove columns where a codon doesn't match its amino acid.
- **Refuse probably-wrong genes** (`--strict`): genes whose DNA is probably not theirs are left out of the output and listed as failed.
- **Keep only marked columns** (`-blockonly`): keep only the columns marked `#` under a CLUSTAL alignment (e.g. Gblocks).
"""

# CHANGED(GUI-5): was "Stop with an error instead of writing genes whose DNA is probably wrong", but the app always
# runs COAL with --keep-going, so a refused gene is left out and listed and the run goes on (core.run_coal
# keep_going=True; MRG-3; tests_app/test_final.py TestFinA6StrictAndBlockonly). Only a partial ID match stops it
# (MRG-1). See overnight/integration/CHANGES.md
STRICT_INFO = ("Genes whose DNA is probably not theirs (unrelated, or where the gene starts is a guess) are left out "
               "of the alignment and listed as failed; the other genes are still aligned. If only some protein and "
               "DNA names match, the run stops instead. (--strict)")


def _code_choices():
    return [(f"{n}: {name}", n) for n, name in core.GENETIC_CODES.items()]


def _format_choices(features):
    return [(core.OUTPUT_FORMATS.get(f, f), f) for f in features.output_formats]


def _hint_parts(report):
    hints = report.get("hints") or []
    if not hints:
        return gr.update(value="", visible=False), gr.update(visible=False), None
    text = "\n\n".join(f"**Tip:** {h.get('text')}" for h in hints)
    code = next((h.get("suggested_codontable") for h in hints if h.get("suggested_codontable")), None)
    btn = gr.update(value=f"Re-run with genetic code {code}", visible=True) if code else gr.update(visible=False)
    return gr.update(value=text, visible=True), btn, code


def _results(report, output_path, report_path, coal_cmd, pal2nal_cmd, stderr, banner,
             view_html="", view_path=None, view_ok=False):
    """banner: HTML from gui_html (CHANGED(GUI-1): was markdown)."""
    hint_md, hint_btn, code = _hint_parts(report)
    # CHANGED(FIN-A5/FIN-C8): was "# Same run with PAL2NAL" (not the same run with -nogap and another code); now
    # "Closest PAL2NAL v14 command" with the caveats. See overnight/integration/CHANGES.md
    cmds = core.commands_text(coal_cmd, pal2nal_cmd, (report or {}).get("options") or {})
    return (
        gr.update(value=banner),
        hint_md,
        hint_btn,
        code,
        gr.update(value=output_path, visible=bool(output_path)),
        gr.update(value=report_path, visible=bool(report_path)),
        gr.update(value=core.health_table_html(report, links=view_ok)),  # CHANGED(FIN-A4): only if the view rendered
        gr.update(value=cmds),
        gr.update(value=stderr or "(no messages)"),
        gr.update(value=view_html or VIEW_EMPTY),
        gr.update(value=view_path, visible=bool(view_path)),
        gr.update(value=core.methods_text(report) if report else "",
                  visible=bool(report and core.methods_text(report))),
    )


VIEW_EMPTY = '<p class="coal-small">The alignment is drawn here after a run.</p>'  # CHANGED(GUI-6)


def run_ui(protein, dna, fmt, code, nogap, nomismatch, protein_path="", dna_path="", strict=False, blockonly=False):
    # CHANGED(APP-L1): was uploads only; now a typed/pasted path (local only) wins over an upload and is read
    # in place, never copied. See overnight/app/CHANGES.md
    try:
        protein, p_desc = core.resolve_input(protein, protein_path)
        dna, d_desc = core.resolve_input(dna, dna_path)
    except ValueError as e:
        return _results({}, None, None, "", "", str(e), gui_html.message_banner("failed", "Could not use that file",
                                                                                str(e)))
    if not protein or not dna:
        how = ("Add a protein alignment and a DNA FASTA file" +
               (" (or give their paths on this computer)" if core.path_inputs_enabled() else ""))
        return _results({}, None, None, "", "", "", gui_html.message_banner(
            "warning", "Please add both files", text_html=f"{how}, or pick one of the examples."))
    out_dir = core.new_run_dir()  # CHANGED(APP-12): runs older than 24 h are pruned
    try:
        # CHANGED(FIN-A6): strict and blockonly options (the app advised --strict but had no way to use it)
        r = core.run_coal(protein, dna, fmt, int(code), bool(nogap), bool(nomismatch), out_dir=out_dir,
                          strict=bool(strict), blockonly=bool(blockonly))
    except ValueError as e:
        return _results({}, None, None, "", "", str(e), gui_html.message_banner("failed", "Could not run", str(e)))
    head = core.banner_head(r.ok, r.report.get("summary", {}), r.hints)  # CHANGED(R2-A1): see core.banner_head
    # CHANGED(R3-A5): was no file names in the banner; now it says which files were used (path box or upload).
    # CHANGED(GUI-1): was markdown (heading, sentence, files); now the HTML summary banner (headline with counts,
    # one tile per gene, the same sentence and files line). See overnight/integration/CHANGES.md
    banner = gui_html.run_banner(head, r.report, core.summary_sentence(r.report),
                                 f"Files used: protein {p_desc}; DNA {d_desc}.")
    view_html, view_path, view_ok = "", None, False
    try:
        # CHANGED(R2-A3): was viewer.view twice (fragment + page); now rendered once
        view_html, view_path = viewer.render_run(protein, r.output_path, r.report, dna, int(code), out_dir)
        view_ok = True
    except (OSError, ValueError) as e:
        # CHANGED(FIN-A4): was the full server path in the message; now file names only
        import html as _html
        view_html = ("<p>Could not draw the alignment: "
                     f"{_html.escape(core.public_error(e, [protein, dna, r.output_path, out_dir]))}</p>")
    return _results(r.report, r.output_path, r.report_path, r.coal_cmd, r.pal2nal_cmd, r.stderr, banner,
                    view_html, view_path, view_ok)


def _clear_path(_upload):
    return ""


def show_sample():
    rep = core.load_sample_report()
    banner = gui_html.sample_banner(rep, core.summary_sentence(rep))  # CHANGED(GUI-1): was core.sample_banner (md)
    return _results(rep, None, core.SAMPLE_REPORT, rep.get("command", ""), rep.get("pal2nal_equivalent", ""),
                    "", banner)


def build(demo_mode: bool = False) -> gr.Blocks:
    """CHANGED(GUI-1): layout, top to bottom: short header; the two files side by side (path boxes under a
    "Use a file path instead" disclosure; examples as one line under them); genetic code and output format, the
    four flags under "More options"; one primary button; results: summary banner (with the genetic-code tip and
    its re-run button), downloads, health table (OK genes in one collapsed row), alignment view, and one
    "Details" area for the methods paragraph, commands, raw messages and the glossary. Every component, event
    and output of the previous layout is kept. See overnight/integration/CHANGES.md"""
    features = core.detect_features()
    with gr.Blocks(title="COAL codon alignment", theme=_theme(), css=CSS, analytics_enabled=False) as demo:
        gr.HTML(gui_html.header_html(core.intro_text()), elem_id="coal-header")
        local = core.path_inputs_enabled()  # CHANGED(APP-L1): path boxes only when running locally
        with gr.Row(equal_height=True):
            protein = gr.File(label="1. Protein alignment (FASTA or CLUSTAL)", type="filepath", height=150,
                              file_types=[".fasta", ".fa", ".faa", ".aln", ".fas", ".clustal", ".txt"],
                              elem_id="coal-protein")
            dna = gr.File(label="2. Coding DNA (FASTA, not aligned)", type="filepath", height=150,
                          file_types=[".fasta", ".fa", ".fna", ".ffn", ".fas", ".txt"], elem_id="coal-dna")
        # CHANGED(GUI-1): the path boxes were always shown under the drop zones; now under a closed disclosure
        # (opened by the desktop wrapper's ?protein=&dna= link). Labels CHANGED(GUI-1): were "or type/paste a path
        # on this computer (big files: read in place, not copied)" on both.
        with gr.Accordion("Use a file path instead (large files)", open=False, visible=local,
                          elem_id="coal-paths") as paths_box:
            gr.Markdown("Type or paste a path on this computer. The file is read in place, not copied, so this is "
                        "the way for big files. A path wins over an uploaded file.", elem_classes=["coal-small"])
            with gr.Row():
                with gr.Column():
                    protein_path = gr.Textbox(label="Protein alignment path", placeholder="/path/to/protein.aln.fasta",
                                              visible=local, max_lines=1, elem_id="coal-protein-path")
                    protein_info = gr.Markdown(visible=local)
                with gr.Column():
                    dna_path = gr.Textbox(label="Coding DNA path", placeholder="/path/to/dna.fasta",
                                          visible=local, max_lines=1, elem_id="coal-dna-path")
                    dna_info = gr.Markdown(visible=local)
        with gr.Row(elem_id="coal-examples"):
            gr.Markdown("No files to hand? Try an example:", elem_classes=["coal-ex-label"])
            ex_buttons = {n: gr.Button(core.EXAMPLES[n][0], size="sm", variant="secondary")
                          for n in core.available_examples()}
            sample_btn = gr.Button("Demo sample report", size="sm", variant="secondary")

        with gr.Row():
            code = gr.Dropdown(_code_choices(), value=1, label="Genetic code",
                               info="Code 1 for nuclear genes. Mitochondrial genes need their own code (2 for vertebrates).")
            fmt = gr.Dropdown(_format_choices(features), value="fasta", label="Output format",
                              info=("Which file format you want back." if len(features.output_formats) > 1 else
                                    "This COAL version writes FASTA only; more formats come with the next engine."))
        # CHANGED(GUI-5): the four flags were always shown, labelled with the flag; now under "More options",
        # plain label first, the flag at the end of the help text
        with gr.Accordion("More options", open=False, elem_id="coal-more"):
            with gr.Row():
                nogap = gr.Checkbox(label="Remove gap columns",
                                    info="Drop columns with a gap or stop codon in any sequence. (-nogap)")
                nomismatch = gr.Checkbox(label="Remove mismatched columns",
                                         info="Drop columns where a codon doesn't match its amino acid. (-nomismatch)")
            with gr.Row():  # CHANGED(FIN-A6): new options
                strict = gr.Checkbox(label="Refuse probably-wrong genes", info=STRICT_INFO)
                blockonly = gr.Checkbox(label="Keep only '#'-marked columns",
                                        info="CLUSTAL input with '#' marks under the blocks, e.g. Gblocks positions. "
                                             "(-blockonly)")
        run_btn = gr.Button("Make codon alignment", variant="primary", size="lg", elem_id="coal-run")

        with gr.Column(elem_id="coal-results"):
            banner = gr.HTML(gui_html.empty_banner(), elem_id="coal-banner")
            hint_md = gr.Markdown(visible=False, elem_id="coal-hint")
            hint_btn = gr.Button("Re-run", visible=False, variant="secondary", elem_id="coal-hint-rerun")
            hint_code = gr.State(None)
            with gr.Row(elem_id="coal-downloads"):
                dl_out = gr.DownloadButton("Download codon alignment", visible=False, variant="secondary")
                dl_rep = gr.DownloadButton("Download report (JSON)", visible=False, variant="secondary")
            table = gr.HTML(elem_id="coal-health-box")
            with gr.Accordion("View alignment (see exactly where a gene breaks)", open=True, elem_id="coal-view"):
                view = gr.HTML(VIEW_EMPTY)
                dl_view = gr.DownloadButton("Download alignment view (HTML)", visible=False, variant="secondary")
            # CHANGED(GUI-1): methods, commands, raw messages and the glossary were four separate areas; now one
            with gr.Accordion("Details: methods text, commands, raw messages, glossary", open=False,
                              elem_id="coal-details"):
                methods = gr.Textbox(label="Methods paragraph (ready to paste into a paper)", lines=4, max_lines=8,
                                     interactive=False, show_copy_button=True, visible=False)
                cmds = gr.Code(label="Equivalent commands (copy into a terminal or pipeline)", language="shell",
                               interactive=False)
                raw = gr.Textbox(label="Raw COAL messages", lines=6, max_lines=20, interactive=False)
                gr.Markdown("**What do these words mean?**\n" + GLOSSARY)

        outputs = [banner, hint_md, hint_btn, hint_code, dl_out, dl_rep, table, cmds, raw, view, dl_view, methods]
        inputs = [protein, dna, fmt, code, nogap, nomismatch, protein_path, dna_path, strict, blockonly]
        run_btn.click(run_ui, inputs, outputs, show_api=False)
        protein_path.change(core.path_status, protein_path, protein_info, show_api=False)
        dna_path.change(core.path_status, dna_path, dna_info, show_api=False)
        # CHANGED(R3-A5): was an upload silently ignored while its path box held a path ("Open with" fills the
        # boxes and a path wins); now uploading a file clears that box. See overnight/integration/CHANGES.md
        protein.upload(_clear_path, protein, protein_path, show_api=False)
        dna.upload(_clear_path, dna, dna_path, show_api=False)

        if local:
            # Desktop wrappers: open http://127.0.0.1:7860/?protein=<path>&dna=<path> (URL-encoded) to fill the
            # path boxes. See overnight/app/CHANGES.md (APP-L2).
            def fill_from_url(request: gr.Request):
                params = dict(request.query_params) if request is not None else {}
                p, d = core.paths_from_query(params)
                # CHANGED(GUI-1): also opens the path disclosure when the link filled a box
                return p, d, gr.update(open=bool(p or d))

            demo.load(fill_from_url, None, [protein_path, dna_path, paths_box], show_api=False)

        for name, btn in ex_buttons.items():
            p, d = core.example_paths(name)

            def load_example(p=p, d=d):
                return p, d, "", ""

            btn.click(load_example, None, [protein, dna, protein_path, dna_path], show_api=False).then(
                run_ui, inputs, outputs, show_api=False)

        sample_btn.click(show_sample, None, outputs, show_api=False)

        def apply_code(c):
            return gr.update(value=int(c)) if c else gr.update()

        hint_btn.click(apply_code, hint_code, code, show_api=False).then(run_ui, inputs, outputs, show_api=False)

        if demo_mode:
            demo.load(show_sample, None, outputs, show_api=False)

        # MCP / API: only these two tools are exposed to AI agents.
        if hasattr(gr, "api"):
            gr.api(tools.codon_align, api_name="codon_align")
            gr.api(tools.explain_report, api_name="explain_report")
    return demo


def launch(host="127.0.0.1", port=7860, mcp=True, demo_mode=False, prevent_thread_lock=False):
    os.makedirs(core.RUNS_DIR, exist_ok=True)
    core.prune_runs()  # CHANGED(APP-12): was never pruned; now runs older than 24 h are removed at startup
    demo = build(demo_mode)
    demo.launch(server_name=host, server_port=port, mcp_server=mcp, inbrowser=False,
                allowed_paths=[core.EXAMPLES_DIR, core.RUNS_DIR, os.path.dirname(core.SAMPLE_REPORT)],
                prevent_thread_lock=prevent_thread_lock, show_error=True)
    return demo
