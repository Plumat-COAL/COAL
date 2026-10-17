"""Gradio web GUI for COAL, which also serves the MCP tools (``launch(mcp_server=True)``).

Written for Gradio 5.x (``pip install -e .[gui]``). All logic lives in ``core``
(tested without Gradio); this file only lays out components and wires events.
Offline: no CDN fonts (system fonts only) and analytics disabled.
"""
from __future__ import annotations

import os
import tempfile

os.environ.setdefault("GRADIO_ANALYTICS_ENABLED", "False")

import gradio as gr  # noqa: E402

from . import core, tools, viewer  # noqa: E402

RUNS_DIR = os.path.join(tempfile.gettempdir(), "coal_app_runs")

# Okabe-Ito colour-blind-safe palette; every status also has an icon and a word.
CSS = """
:root { --coal-ok:#0072B2; --coal-warn:#E69F00; --coal-fail:#D55E00; --coal-muted:#666666; }
.gradio-container { font-size: 18px; }
.coal-health { width:100%; border-collapse:collapse; font-size:1rem; }
.coal-health caption { text-align:left; font-weight:600; padding:0.4em 0; }
.coal-health th, .coal-health td { text-align:left; padding:0.45em 0.6em; border-bottom:1px solid #8884; vertical-align:top; }
.coal-health tr.coal-ok td:first-child { border-left:6px solid var(--coal-ok); }
.coal-health tr.coal-warning td:first-child { border-left:6px solid var(--coal-warn); }
.coal-health tr.coal-failed td:first-child { border-left:6px solid var(--coal-fail); }
.coal-health tr.coal-unchecked td:first-child { border-left:6px dashed var(--coal-muted); }
.coal-badge { display:inline-block; min-width:1.4em; font-weight:700; }
.coal-badge.coal-ok { color:var(--coal-ok); } .coal-badge.coal-warning { color:#B07800; }
.coal-badge.coal-failed { color:var(--coal-fail); }
.coal-small { font-size:0.85em; opacity:0.8; }
.coal-hint { border:2px solid var(--coal-warn); border-radius:8px; padding:0.6em 0.9em; }
abbr[title] { text-decoration: underline dotted; cursor: help; }
button:focus-visible, input:focus-visible, select:focus-visible, a:focus-visible {
  outline:3px solid var(--coal-ok) !important; outline-offset:2px; }
"""

INTRO = """
# COAL: codon alignment made easy
Give COAL a **protein alignment** and the **DNA** of the same genes. It lines the DNA up
<abbr title="Three DNA letters that code for one amino acid">codon</abbr> by codon to match
the protein alignment, and tells you in plain words about any gene that doesn't fit.
Nothing leaves your computer.
"""

GLOSSARY = """
- **Codon**: three DNA letters that code for one amino acid (for example ATG = M, methionine).
- **Protein alignment**: protein sequences lined up so that matching positions are in the same column (from MAFFT, MUSCLE, Clustal...).
- **Genetic code**: the dictionary from codons to amino acids. Almost all nuclear genes use code 1; mitochondria and some microbes use slightly different ones.
- **Mismatch**: a codon that doesn't code for the amino acid at that position. A few are often a sign of the wrong genetic code.
- **Frameshift**: a missing or extra DNA letter, so every later codon is read out of step.
- **-nogap**: remove columns containing gaps or stop codons (some dN/dS tools need this).
- **-nomismatch**: remove columns where a codon doesn't match its amino acid.
"""


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
             view_html="", view_path=None):
    hint_md, hint_btn, code = _hint_parts(report)
    cmds = f"# COAL\n{coal_cmd}\n\n# Same run with PAL2NAL\n{pal2nal_cmd}"
    return (
        gr.update(value=banner),
        hint_md,
        hint_btn,
        code,
        gr.update(value=output_path, visible=bool(output_path)),
        gr.update(value=report_path, visible=bool(report_path)),
        gr.update(value=core.health_table_html(report, links=bool(view_html))),
        gr.update(value=cmds),
        gr.update(value=stderr or "(no messages)"),
        gr.update(value=view_html or "<p>Run COAL to see the alignment here.</p>"),
        gr.update(value=view_path, visible=bool(view_path)),
    )


def run_ui(protein, dna, fmt, code, nogap, nomismatch):
    if not protein or not dna:
        return _results({}, None, None, "", "", "",
                        "### ⚠ Please add both files\nUpload a protein alignment and a DNA FASTA file, "
                        "or press one of the **Try an example** buttons.")
    os.makedirs(RUNS_DIR, exist_ok=True)
    out_dir = tempfile.mkdtemp(prefix="run_", dir=RUNS_DIR)
    try:
        r = core.run_coal(protein, dna, fmt, int(code), bool(nogap), bool(nomismatch), out_dir=out_dir)
    except ValueError as e:
        return _results({}, None, None, "", "", str(e), f"### ✖ Could not run\n{e}")
    s = r.report.get("summary", {})
    if r.ok and not s.get("failed") and not s.get("warning"):
        head = "### ✔ Done: every sequence matched"
    elif r.ok:
        head = "### ⚠ Done, with some problems (see the table)"
    else:
        head = "### ✖ COAL stopped: no alignment was written"
    banner = f"{head}\n{core.summary_sentence(r.report)}"
    view_html, view_path = "", None
    try:
        view_html = viewer.view(protein, r.output_path, r.report, dna, int(code), page=False)
        view_path = os.path.join(out_dir, "alignment_view.html")
        with open(view_path, "w", encoding="utf-8") as fh:
            fh.write(viewer.view(protein, r.output_path, r.report, dna, int(code), page=True))
    except (OSError, ValueError) as e:
        view_html = f"<p>Could not draw the alignment: {e}</p>"
    return _results(r.report, r.output_path, r.report_path, r.coal_cmd, r.pal2nal_cmd, r.stderr, banner,
                    view_html, view_path)


def show_sample():
    rep = core.load_sample_report()
    banner = ("### Sample report (not your data)\nThis is what a run looks like with a newer COAL that writes "
              "reports. " + core.summary_sentence(rep))
    return _results(rep, None, core.SAMPLE_REPORT, rep.get("command", ""), rep.get("pal2nal_equivalent", ""),
                    "", banner)


def build(demo_mode: bool = False) -> gr.Blocks:
    features = core.detect_features()
    theme = gr.themes.Default(
        font=["system-ui", "Segoe UI", "Helvetica", "Arial", "sans-serif"],
        font_mono=["ui-monospace", "Consolas", "Menlo", "monospace"],
        text_size=gr.themes.sizes.text_lg,
    )
    with gr.Blocks(title="COAL codon alignment", theme=theme, css=CSS, analytics_enabled=False) as demo:
        gr.Markdown(INTRO)
        with gr.Row():
            protein = gr.File(label="1. Protein alignment (FASTA or CLUSTAL)", type="filepath",
                              file_types=[".fasta", ".fa", ".faa", ".aln", ".fas", ".clustal", ".txt"])
            dna = gr.File(label="2. Coding DNA (FASTA, not aligned)", type="filepath",
                          file_types=[".fasta", ".fa", ".fna", ".ffn", ".fas", ".txt"])
        with gr.Row():
            fmt = gr.Dropdown(_format_choices(features), value="fasta", label="Output format",
                              info=("Which file format you want back." if len(features.output_formats) > 1 else
                                    "This COAL version writes FASTA only; more formats come with the next engine."))
            code = gr.Dropdown(_code_choices(), value=1, label="Genetic code",
                               info="Code 1 for nuclear genes. Mitochondrial genes need their own code (2 for vertebrates).")
        with gr.Row():
            nogap = gr.Checkbox(label="Remove gap columns (-nogap)",
                                info="Drop columns with a gap or stop codon in any sequence.")
            nomismatch = gr.Checkbox(label="Remove mismatched columns (-nomismatch)",
                                     info="Drop columns where a codon doesn't match its amino acid.")
        run_btn = gr.Button("Make codon alignment", variant="primary")

        with gr.Row():
            gr.Markdown("**Try an example:**")
            ex_buttons = {n: gr.Button(core.EXAMPLES[n][0], size="sm") for n in core.available_examples()}
            sample_btn = gr.Button("Show a sample report", size="sm")

        banner = gr.Markdown("Results will appear here.")
        hint_md = gr.Markdown(visible=False, elem_classes=["coal-hint"])
        hint_btn = gr.Button("Re-run", visible=False, variant="primary")
        hint_code = gr.State(None)
        with gr.Row():
            dl_out = gr.DownloadButton("Download codon alignment", visible=False)
            dl_rep = gr.DownloadButton("Download report (JSON)", visible=False)
        table = gr.HTML()
        with gr.Accordion("View alignment (see exactly where a gene breaks)", open=True):
            view = gr.HTML("<p>Run COAL to see the alignment here.</p>")
            dl_view = gr.DownloadButton("Download alignment view (HTML)", visible=False)
        cmds = gr.Code(label="Equivalent commands (copy into a terminal or pipeline)", language="shell",
                       interactive=False)
        with gr.Accordion("Raw COAL messages", open=False):
            raw = gr.Textbox(show_label=False, lines=8, max_lines=20, interactive=False)
        with gr.Accordion("What do these words mean?", open=False):
            gr.Markdown(GLOSSARY)

        outputs = [banner, hint_md, hint_btn, hint_code, dl_out, dl_rep, table, cmds, raw, view, dl_view]
        inputs = [protein, dna, fmt, code, nogap, nomismatch]
        run_btn.click(run_ui, inputs, outputs, show_api=False)

        for name, btn in ex_buttons.items():
            p, d = core.example_paths(name)

            def load_example(p=p, d=d):
                return p, d

            btn.click(load_example, None, [protein, dna], show_api=False).then(
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
    os.makedirs(RUNS_DIR, exist_ok=True)
    demo = build(demo_mode)
    demo.launch(server_name=host, server_port=port, mcp_server=mcp, inbrowser=False,
                allowed_paths=[core.EXAMPLES_DIR, RUNS_DIR, os.path.dirname(core.SAMPLE_REPORT)],
                prevent_thread_lock=prevent_thread_lock, show_error=True)
    return demo
