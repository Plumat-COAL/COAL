"""Gradio web GUI for COAL, which also serves the MCP tools (``launch(mcp_server=True)``).

Written for Gradio 5.x (``pip install -e "apps/coal-app[gui]"``). All logic lives in ``core``
(tested without Gradio); this file only lays out components and wires events.
Offline: no CDN fonts (system fonts only) and analytics disabled.
"""
from __future__ import annotations

import os

os.environ.setdefault("GRADIO_ANALYTICS_ENABLED", "False")

import gradio as gr  # noqa: E402

from . import core, gui_html, settings, tools, viewer  # noqa: E402


# The COAL theme: colours of the COAL logo and site (charcoal #2A2A2A, slate #4A4744, ember #FF6A13 under dark text
# #1E1E1E, amber #FFC24B, page #121110 in dark), system fonts only (no web fonts are requested), for light and for
# Gradio's dark theme (".dark" on the body). Text contrast is at least 4.5:1 in both. Every status keeps an icon and
# a word, never colour alone.
# The accent (primary button, active tab, focus ring, link hover, drop-zone hover, "Show N more") is the logo's amber
# #FFC24B (button fill, dark text #1E1E1E, 10.4:1) and, where it is text or a thin line on the light theme, deep gold
# #8A5A00 (5.4:1 on the page). So that no status looks like the accent, warning is ember (fill #D45400 light /
# #FF6A13 dark, text #A83E00 / #FF8A45) and failed is crimson (fill #E0284F, text #B8143F / #FF8597).
SANS = ["Seravek", "Gill Sans Nova", "Ubuntu", "Calibri", "Segoe UI", "DejaVu Sans", "system-ui", "sans-serif"]
MONO = ["ui-monospace", "Cascadia Mono", "SF Mono", "Menlo", "Consolas", "DejaVu Sans Mono", "monospace"]

CSS = """
:root, .light { --coal-page:#F6F5F3; --coal-surface:#FFFFFF; --coal-ink:#1E1E1E; --coal-text:#2A2A2A;
  --coal-muted:#5E5853; --coal-line:#DEDAD6; --coal-ok:#4A4744; --coal-warn:#A83E00; --coal-warn-fill:#D45400;
  --coal-fail:#B8143F; --coal-fail-fill:#E0284F; --coal-info:#0B5CAD; --coal-focus:#8A5A00;
  --coal-accent:#FFC24B; --coal-accent-line:#8A5A00;
  --coal-tile-ok:#4A4744; --coal-tile-line:#8C8782; --coal-m1:#2A2A2A; --coal-m2:#4A4744; --coal-field:#8C8782; }
.dark { --coal-page:#121110; --coal-surface:#1C1A19; --coal-ink:#F1EEEA; --coal-text:#E4E0DB;
  --coal-muted:#B3ADA6; --coal-line:#3A3633; --coal-ok:#B5AFA8; --coal-warn:#FF8A45; --coal-warn-fill:#FF6A13;
  --coal-fail:#FF8597; --coal-fail-fill:#E0284F; --coal-info:#8CC4FF; --coal-focus:#FFC24B;
  --coal-accent:#FFC24B; --coal-accent-line:#FFC24B;
  --coal-tile-ok:#6A655F; --coal-tile-line:#8C8782; --coal-m1:#8C8782; --coal-m2:#B5AFA8; --coal-field:#6A655F; }
/* Gradio's theme has one --color-accent for both themes (active tab text and underline); the dark
   theme takes amber, the light one the deep gold */
.dark, .dark .gradio-container { --color-accent:#FFC24B; --color-accent-soft:#3A2E14; }
.gradio-container { max-width:1080px !important; margin:0 auto !important; }
/* header */
/* The header is the mark and the name only */
.coal-head { display:flex; gap:10px; align-items:center; padding:4px 0 0; }
.coal-head .coal-mark { flex:none; }
.coal-mark .m1 { fill:var(--coal-m1); } .coal-mark .m2 { fill:var(--coal-m2); }
.coal-mark .m3 { fill:#FF6A13; } .coal-mark .m4 { fill:#FFC24B; }
.coal-head h1 { margin:0; font-size:1.5rem; line-height:1.2; color:var(--coal-ink); letter-spacing:0.02em; }
.coal-foot { margin:0; font-size:0.85em; color:var(--coal-muted); }
/* example row: one quiet line of link-style buttons */
#coal-examples { flex-wrap:wrap !important; gap:2px 14px !important; align-items:baseline; }
#coal-examples > * { flex:0 0 auto !important; min-width:0 !important; width:auto !important; }
#coal-examples .coal-ex-label p, #coal-examples .coal-ex-label { margin:0; color:var(--coal-muted); }
#coal-examples button { background:none !important; border:0 !important; box-shadow:none !important;
  padding:2px 0 !important; min-height:0 !important; font-weight:400 !important; color:var(--coal-info) !important;
  text-decoration:underline; text-underline-offset:3px; }
#coal-examples button:hover { text-decoration-thickness:2px; color:var(--coal-accent-line) !important; }  /* Hover in the accent */
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
/* The genetic-code banner and its re-run button as one box (where :has() works; elsewhere the
   button sits right under the banner) */
#coal-fixbox:has(.coal-codefix) { border:1px solid var(--coal-line); border-left:6px solid var(--coal-info);
  border-radius:8px; padding:14px 18px; background:var(--coal-surface); gap:12px !important; }
#coal-fixbox:has(.coal-codefix) .coal-banner { border:0; border-radius:0; padding:0; background:none; }
.coal-codefix .coal-detail { margin-top:2px !important; }
.coal-codefix .coal-small { margin-top:6px !important; }
.coal-empty-note { margin:0; color:var(--coal-muted); }
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
.coal-tile.coal-failed { background:var(--coal-fail-fill); color:#FFFFFF; }  /* White ✖ on crimson 4.6:1 */
.coal-tile.coal-unchecked { border:1.5px dashed var(--coal-tile-line); color:var(--coal-muted); line-height:11px; }
.coal-detail { color:var(--coal-text); line-height:1.5; }
.coal-files { margin-top:4px !important; font-size:0.88em; color:var(--coal-muted); overflow-wrap:anywhere; }
#coal-hint { border-left:6px solid var(--coal-accent-line);
  border-radius:8px; padding:0.6em 0.9em; }
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
.coal-health a:hover, .coal-health summary:hover { text-decoration-color:var(--coal-accent-line); text-decoration-thickness:2px; }
.coal-badge { display:inline-block; min-width:1.2em; font-weight:700; }
.coal-badge.coal-ok { color:var(--coal-ok); } .coal-badge.coal-warning { color:var(--coal-warn); }
.coal-badge.coal-failed { color:var(--coal-fail); } .coal-badge.coal-unchecked { color:var(--coal-muted); }
.coal-small { font-size:0.88em; color:var(--coal-muted); }
abbr[title] { text-decoration:underline dotted; cursor:help; }
button:focus-visible, input:focus-visible, select:focus-visible, textarea:focus-visible, a:focus-visible,
summary:focus-visible, [tabindex]:focus-visible {
  outline:3px solid var(--coal-focus) !important; outline-offset:2px; }
/* The drop zones. Gradio's upload text has min-height 15rem; with height=150 it overflowed the
   zone and its icon covered the label. No fixed height now (the zone is as tall as its content, never clipped)
   and, where Gradio's markup is as in 5.49 (button[aria-dropeffect] > div.wrap), a compact text block. */
#coal-protein button[aria-dropeffect] > div, #coal-dna button[aria-dropeffect] > div {
  min-height:0 !important; padding:10px 0 !important; }
/* One flat form. No panels around fields (theme: block border 0, block fill = page); the drop
   zones read like the other fields: their label above (not a chip on the border), a dashed drop area under it with
   one line of text and the accepted formats as a quiet hint. Field borders #8C8782 / #6A655F (3:1 or more on the
   page and on the field). Where Gradio's markup differs, Gradio's own text and chip show instead. */
#coal-protein, #coal-dna { --block-border-width:0px; background:transparent !important; box-shadow:none !important;
  overflow:visible !important; }
#coal-protein label[data-testid="block-label"], #coal-dna label[data-testid="block-label"] {
  position:static !important; background:none !important; border:0 !important; box-shadow:none !important;
  padding:0 0 var(--spacing-lg) !important; color:var(--block-title-text-color) !important;
  font-size:var(--block-title-text-size) !important; font-weight:var(--block-title-text-weight) !important; }
#coal-protein label[data-testid="block-label"] > span, #coal-dna label[data-testid="block-label"] > span {
  display:none !important; }
#coal-protein button[aria-dropeffect], #coal-dna button[aria-dropeffect] { min-height:128px; height:auto !important;
  border:1.5px dashed var(--coal-field) !important; border-radius:8px; background:var(--coal-surface) !important; }
#coal-protein button[aria-dropeffect]:hover, #coal-dna button[aria-dropeffect]:hover { border-color:var(--coal-accent-line) !important; }
#coal-protein button[aria-dropeffect] > div, #coal-dna button[aria-dropeffect] > div { font-size:0 !important; gap:2px; }
#coal-protein button[aria-dropeffect] > div > .or, #coal-dna button[aria-dropeffect] > div > .or { display:none !important; }
#coal-protein button[aria-dropeffect] > div::before, #coal-dna button[aria-dropeffect] > div::before {
  order:1; content:"Drop a file here or click to choose"; font-size:var(--text-md); color:var(--coal-text); }
#coal-protein button[aria-dropeffect] > div::after, #coal-dna button[aria-dropeffect] > div::after {
  order:2; font-size:var(--text-sm); color:var(--coal-muted); }
#coal-protein button[aria-dropeffect] > div::after { content:"FASTA or CLUSTAL"; }
#coal-dna button[aria-dropeffect] > div::after { content:"FASTA, not aligned"; }
#coal-code .wrap, #coal-format .wrap, #coal-set-code .wrap, #coal-set-format .wrap { border-color:var(--coal-field) !important; }
#coal-more, #coal-paths, #coal-view, #coal-details { border-top:1px solid var(--coal-line) !important;
  border-radius:0 !important; }
/* A scrollbar showed under every collapsed accordion. Gradio 5.49.1 (js/accordion
   Accordion.svelte) turns the title's ▼ icon by rotate(90deg) when closed; the icon is a flex item taller than
   wide, so its turned box stuck out a few px past the right edge, and the block (js/atoms Block.svelte: inline
   overflow:auto) scrolled. It shows because the form's blocks have no side padding (the icon sits on the
   edge). The icon is a square now, so turning it changes nothing outside it. */
button.label-wrap > span.icon { display:inline-flex; align-items:center; justify-content:center; flex:none;
  width:1.4em; height:1.4em; line-height:1; transform-origin:50% 50%; }
#coal-tabs [role="tablist"] > button:first-child { padding-left:0; }
/* Gradio joins the fields of a row into one .form with no gap; keep a gap between the two boxes */
.form:has(#coal-code), .form:has(#coal-set-code) { gap:var(--layout-gap) !important; background:transparent !important; }
#coal-hint-rerun, #coal-set-reset, #coal-view-dl { align-self:flex-start; width:auto !important; flex:0 0 auto !important; }
.coal-health tr.coal-shared td { color:var(--coal-muted); border-bottom:0; padding-bottom:0.2em; }
/* The rest of a long group of rows, under "Show N more" (a nested table of the same columns) */
.coal-health tr.coal-morerows > td { white-space:normal; padding:0.45em 0; }
.coal-health tr.coal-morerows > td > details > summary { padding-left:calc(6px + 0.7em); color:var(--coal-accent-line); }
.coal-health .coal-nested { margin-top:0.45em; }
.coal-health .coal-nested tr:last-child td { border-bottom:0; }
.coal-sr { position:absolute !important; width:1px; height:1px; padding:0; margin:-1px; overflow:hidden;
  clip:rect(0 0 0 0); white-space:nowrap; border:0; }
/* in the app the health table above already lists every problem with its jump link: the view's own heading and
   problem list (kept in the downloadable page) are not repeated here */
#coal-view div.coalv > h2, #coal-view div.coalv > ul.issues, #coal-view div.coalv > details.hidden-issues {
  display:none !important; }
/* Settings tab */
#coal-settings-tab-button { margin-left:auto; }
#coal-settings-status p { margin:0; color:var(--coal-muted); }
#coal-centre { align-items:flex-end !important; gap:8px !important; flex-wrap:wrap !important; }
#coal-centre > * { flex:0 1 auto !important; min-width:0 !important; }
#coal-view-dl-note p { margin:0; }
/* Section headings, as the site's type scale (site/style.css h3: 1.1875rem, bold, ink): one per
   group (Files, Genetic code, Output, Results); more space above a heading than between the fields of a group;
   field labels lighter than headings (theme: block title #4A4744 / #B3ADA6, regular weight) */
.coal-sec { margin:0 !important; padding:0.9rem 0 0 !important; font-size:1.1875rem !important; font-weight:700 !important;
  line-height:1.25 !important; color:var(--coal-ink) !important; letter-spacing:0 !important; border:0 !important; }
.coal-sec.coal-sec-first { padding-top:0.25rem !important; }
/* The genetic code is its own section: its label is drawn as the section heading (the dropdown
   keeps it as its accessible name), a one-line hint under it, and a larger value than the other fields */
#coal-code span[data-testid="block-info"] { margin-top:0.9rem; font-size:1.1875rem; font-weight:700; line-height:1.25;
  color:var(--coal-ink); }
#coal-code input { font-size:1.125rem !important; font-weight:600; }
/* Gradio wraps each field in a .form (fill = border colour); the narrower boxes keep the page behind them */
.form:has(#coal-code) { max-width:640px; }
/* Output format is secondary: a narrower box under the "Output" heading */
.form:has(#coal-format) { max-width:320px; background:transparent !important; }
/* Gradio's own footer link "Settings" (js/core/src/Blocks.svelte in 5.49.1: footer >
   button.settings, after a div.divider) opens Gradio's panel (display language, theme), a second settings place
   whose language choice only translates Gradio's own labels. COAL has one: its Settings tab. Hidden with the "·"
   before it; "Use via API or MCP" and "Built with Gradio" stay. */
footer > button.settings, footer > div.divider:nth-last-child(2) { display:none !important; }
@media (max-width: 640px) {
  .coal-head { gap:12px; } .coal-head .coal-mark { width:32px; height:25px; }
  .coal-head h1 { font-size:1.4rem; } .coal-headline { font-size:1.2rem; }
  .coal-health th, .coal-health td { padding:0.45em 0.45em; }
  .coal-health td:first-child { white-space:normal; } }
@media (prefers-reduced-motion: reduce) { * { transition:none !important; animation:none !important; } }
"""


def _theme():
    """COAL theme. Fonts are plain strings (never gr.themes.GoogleFont): nothing is fetched."""
    stone = gr.themes.Color(c50="#F6F5F3", c100="#EEECE9", c200="#DEDAD6", c300="#C4BEB8", c400="#A39C95",
                            c500="#7D7670", c600="#5E5853", c700="#4A4744", c800="#2E2B28", c900="#1C1A19",
                            c950="#121110", name="coal_stone")
    # primary hue: the logo's amber #FFC24B (c500) down to deep gold
    amber = gr.themes.Color(c50="#FFF8EB", c100="#FFEFCF", c200="#FFE0A0", c300="#FFD27A", c400="#FFCB63",
                            c500="#FFC24B", c600="#E6A82E", c700="#B07A00", c800="#8A5A00", c900="#5E3D00",
                            c950="#3A2600", name="coal_amber")
    theme = gr.themes.Base(primary_hue=amber, secondary_hue=stone, neutral_hue=stone, font=SANS, font_mono=MONO,
                           text_size=gr.themes.sizes.text_lg, radius_size=gr.themes.sizes.radius_md)
    return theme.set(
        body_background_fill="#F6F5F3", body_background_fill_dark="#121110",
        body_text_color="#2A2A2A", body_text_color_dark="#E4E0DB",
        body_text_color_subdued="#5E5853", body_text_color_subdued_dark="#B3ADA6",
        background_fill_primary="#FFFFFF", background_fill_primary_dark="#1C1A19",
        background_fill_secondary="#F6F5F3", background_fill_secondary_dark="#161413",
        # the fields sit on the page (block fill = page colour, no border, no side padding so every edge lines up)
        block_background_fill="#F6F5F3", block_background_fill_dark="#121110",
        block_border_width="0px", block_border_width_dark="0px", block_shadow="none", block_shadow_dark="none",
        block_padding="*spacing_md 0px", layout_gap="*spacing_xl",
        block_border_color="#DEDAD6", block_border_color_dark="#3A3633",
        border_color_primary="#DEDAD6", border_color_primary_dark="#3A3633",
        block_label_text_color="#4A4744", block_label_text_color_dark="#B3ADA6",
        # field labels one step lighter than the section headings; #4A4744 on the page 8.6:1, #B3ADA6 on the dark
        # page 8.9:1
        block_title_text_color="#4A4744", block_title_text_color_dark="#B3ADA6",
        block_info_text_color="#5E5853", block_info_text_color_dark="#B3ADA6",
        # fields sit on the page, so they need their own edge: a 1px border (Base has none) and, in dark mode, the
        # surface colour (the page colour would hide the field)
        input_background_fill="#FFFFFF", input_background_fill_dark="#1C1A19", input_border_width="1px",
        input_border_color="#8C8782", input_border_color_dark="#6A655F",
        input_border_color_focus="#8A5A00", input_border_color_focus_dark="#FFC24B",
        # primary button amber #FFC24B (hover #FFD27A) with #1E1E1E text, 10.4:1 (11.7:1 hovered).
        # Light theme: a deep gold edge #8A5A00 (5.4:1 on the page) outlines it
        button_primary_background_fill="#FFC24B", button_primary_background_fill_dark="#FFC24B",
        button_primary_background_fill_hover="#FFD27A", button_primary_background_fill_hover_dark="#FFD27A",
        button_primary_text_color="#1E1E1E", button_primary_text_color_dark="#1E1E1E",
        button_primary_border_color="#8A5A00", button_primary_border_color_dark="#FFC24B",
        button_secondary_background_fill="#FFFFFF", button_secondary_background_fill_dark="#2E2B28",
        button_secondary_background_fill_hover="#EEECE9", button_secondary_background_fill_hover_dark="#3A3633",
        button_secondary_text_color="#2A2A2A", button_secondary_text_color_dark="#F1EEEA",
        button_secondary_border_color="#A39C95", button_secondary_border_color_dark="#5A534D",
        # deep gold in both themes so Gradio's white tick stays readable (5.9:1)
        checkbox_background_color_selected="#8A5A00", checkbox_background_color_selected_dark="#8A5A00",
        checkbox_border_color_selected="#8A5A00", checkbox_border_color_selected_dark="#FFC24B",
        link_text_color="#0B5CAD", link_text_color_dark="#8CC4FF",
        color_accent="#8A5A00", color_accent_soft="#FFEFCF", color_accent_soft_dark="#3A2E14",  # dark: CSS
    )


# the plain names used on the options, with the flag after them
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

# The app always runs COAL with --keep-going, so a refused gene is left out and listed and the run goes on
# (core.run_coal keep_going=True; tests/test_mcp_tools.py TestStrictAndBlockonlyOptions). Only a partial ID match
# stops it, whatever this option says (the glossary and the docs say so).
STRICT_INFO = "They are left out and listed as failed; the other genes are still aligned. (--strict)"


# Short names for the dropdowns (at most 38 characters; longer ones are cut off at about 40 in an 800 px window; the
# genetic-code box is the wider of the two). core.GENETIC_CODES has the full names, used in tips and the methods
# text.
CODE_SHORT: dict[int, str] = {
    1: "Standard (nuclear genes)",
    2: "Vertebrate mitochondrial",
    3: "Yeast mitochondrial",
    4: "Protozoan mitochondrial, Mycoplasma",
    5: "Invertebrate mitochondrial",
    6: "Ciliate and Hexamita nuclear",
    9: "Echinoderm, flatworm mitochondrial",
    10: "Euplotid nuclear",
    11: "Bacteria, archaea and plastids",
    12: "Alternative yeast nuclear",
    13: "Ascidian mitochondrial",
    14: "Alternative flatworm mitochondrial",
    15: "Blepharisma nuclear",
    16: "Chlorophycean mitochondrial",
    21: "Trematode mitochondrial",
    22: "Scenedesmus mitochondrial",
    23: "Thraustochytrium mitochondrial",
}
FORMAT_SHORT = {"fasta": "FASTA", "clustal": "CLUSTAL", "paml": "PAML (for codeml)", "codon": "Codon (spaced)"}
# not shown as a help line: the option names carry the meaning ("1: Standard (nuclear genes)", "2: Vertebrate
# mitochondrial", "11: Bacteria, archaea and plastids")
CODE_INFO = "1 for nuclear genes, 2 for vertebrate mitochondria, 11 for bacteria and plastids."


# the one line under the genetic code (it is the choice that most often makes results wrong)
CODE_HINT = "Mitochondrial, bacterial and plastid genes need their own code."


def _sec(text: str, first: bool = False) -> str:
    """A section heading of the form (h2; the page's h1 is COAL)."""
    return f'<h2 class="coal-sec{" coal-sec-first" if first else ""}">{text}</h2>'


def _code_choices():
    return [(f"{n}: {CODE_SHORT.get(n, name)}", n) for n, name in core.GENETIC_CODES.items()]


def _format_choices(features):
    return [(FORMAT_SHORT.get(f, core.OUTPUT_FORMATS.get(f, f)), f) for f in features.output_formats]


def _format_info(features):
    """No help line: the option names say it ("PAML (for codeml)"). Only when the engine writes FASTA alone does
    a line say why the box has one choice (Microsoft: explain an unavailable setting)."""
    if len(features.output_formats) <= 1:
        return "Only FASTA with this COAL version."
    return None


def _hint_parts(report):
    hints = report.get("hints") or []
    if not hints:
        return gr.update(value="", visible=False), gr.update(visible=False, variant="secondary"), None
    fix = core.code_fix(report)
    if fix:  # the banner says it all; the re-run is its primary action, under it in the same box
        return (gr.update(value="", visible=False),
                gr.update(value=f"Re-run with genetic code {fix['code']}", visible=True, variant="primary"),
                fix["code"])
    code = next((h.get("suggested_codontable") for h in hints if h.get("suggested_codontable")), None)
    # the tip's closing question ("Re-run with code 2?") is the button right under it
    text = "\n\n".join("**Tip:** " + _strip_rerun(h.get("text") or "", h.get("suggested_codontable"))
                         for h in hints)
    btn = (gr.update(value=f"Re-run with genetic code {code}", visible=True, variant="secondary") if code
           else gr.update(visible=False, variant="secondary"))
    return gr.update(value=text, visible=True), btn, code


def _strip_rerun(text: str, code) -> str:
    q = f" Re-run with code {code}?"
    return text[:-len(q)] if code and text.endswith(q) else text


def _results(report, output_path, report_path, coal_cmd, pal2nal_cmd, stderr, banner,
             view_html="", view_path=None, view_ok=False, preview=None, ctx=None):
    """banner: HTML from gui_html.
    preview (what the alignment preview drew, viewer.view's info) and ctx (what is needed to
    draw it again centred elsewhere); two more outputs: the run context and the "Centre on column" row."""
    hint_md, hint_btn, code = _hint_parts(report)
    # "Closest PAL2NAL v14 command" with its caveats (PAL2NAL cannot always do the same run)
    cmds = core.commands_text(coal_cmd, pal2nal_cmd, (report or {}).get("options") or {})
    return (
        gr.update(value=banner),
        hint_md,
        hint_btn,
        code,
        gr.update(value=output_path, visible=bool(output_path)),
        gr.update(value=report_path, visible=bool(report_path)),
        # links only if the view rendered, and only to cells in the preview
        gr.update(value=core.health_table_html(report, links=view_ok, preview=preview)),
        gr.update(value=cmds),
        gr.update(value=stderr or "(no messages)"),
        gr.update(value=view_html or VIEW_EMPTY),
        gr.update(value=view_path, visible=bool(view_path)),
        gr.update(value=core.methods_text(report) if report else "",
                  visible=bool(report and core.methods_text(report))),
        ctx,
        gr.update(visible=bool(preview and preview.get("preview") and preview.get("last", 0) - preview.get("first", 1)
                               + 1 < preview.get("ncols", 0))),
        gr.update(visible=True),  # the results area, hidden until the first result
    )


def _view_tail(preview) -> str:
    """A line under a preview, above its download button."""
    if not (preview and preview.get("preview")):
        return ""
    # one sentence: the preview line above the grid already says how to see everything
    return ('<p class="coal-small" id="coal-view-dl-note">The downloadable view is this preview, not the whole '
            'alignment.</p>')


VIEW_EMPTY = '<p class="coal-small">Drawn after a run.</p>'


def run_ui(protein, dna, fmt, code, nogap, nomismatch, protein_path="", dna_path="", strict=False, blockonly=False,
           preview_size="standard"):
    # preview_size: the Settings tab's "Alignment preview" (settings.PREVIEW_SIZES)
    # a typed/pasted path (local only) wins over an upload and is read in place, never copied
    try:
        protein, p_desc = core.resolve_input(protein, protein_path)
        dna, d_desc = core.resolve_input(dna, dna_path)
    except ValueError as e:
        return _results({}, None, None, "", "", str(e), gui_html.message_banner("failed", "Could not use that file",
                                                                                str(e)))
    if not protein or not dna:
        # says which file is missing
        missing = " and ".join(n for n, v in (("a protein alignment", protein), ("a coding DNA file", dna)) if not v)
        return _results({}, None, None, "", "", "", gui_html.message_banner(
            "warning", "Please add both files", text=f"Add {missing}, or try an example."))
    out_dir = core.new_run_dir()  # runs older than 24 h are pruned
    try:
        r = core.run_coal(protein, dna, fmt, int(code), bool(nogap), bool(nomismatch), out_dir=out_dir,
                          strict=bool(strict), blockonly=bool(blockonly))
    except ValueError as e:
        return _results({}, None, None, "", "", str(e), gui_html.message_banner("failed", "Could not run", str(e)))
    head = core.banner_head(r.ok, r.report.get("summary", {}), r.hints)  # see core.banner_head
    # the HTML summary banner (headline with counts, one tile per gene, the sentence and which files were used:
    # path box or upload)
    banner = gui_html.run_banner(head, r.report, core.summary_sentence(r.report),
                                 f"Files used: protein {p_desc}; DNA {d_desc}.")
    fix = core.code_fix(r.report)
    if fix:  # the genetic-code hint explains the run: say so first, with the fix
        banner = gui_html.code_fix_banner(fix, r.report, CODE_SHORT, f"Files used: protein {p_desc}; DNA {d_desc}.")
    view_html, view_path, view_ok = "", None, False
    info: dict = {}
    ctx = {"protein": protein, "output": r.output_path, "report": r.report_path, "dna": dna, "code": int(code),
           "out_dir": out_dir, "preview": preview_size}
    try:
        # fragment and page rendered once
        view_html, view_path = viewer.render_run(protein, r.output_path, r.report, dna, int(code), out_dir,
                                                 limits=settings.preview_limits(preview_size), info=info)
        view_html += _view_tail(info)
        view_ok = True
    except (OSError, ValueError) as e:
        # file names only, no server paths, in the message
        import html as _html
        view_html = ("<p>Could not draw the alignment: "
                     f"{_html.escape(core.public_error(e, [protein, dna, r.output_path, out_dir]))}</p>")
    return _results(r.report, r.output_path, r.report_path, r.coal_cmd, r.pal2nal_cmd, r.stderr, banner,
                    view_html, view_path, view_ok, preview=info if view_ok else None, ctx=ctx if view_ok else None)


def centre_ui(ctx, column):
    """Draw the preview again, centred on `column` (the "Centre the preview on column" box).
    Outputs: alignment view, its download button, health table."""
    import html as _html
    if not ctx:
        return gr.update(), gr.update(), gr.update()
    try:
        col = int(column)
    except (TypeError, ValueError):
        return (gr.update(value=f'<p class="coal-small">Type a column number first.</p>{VIEW_EMPTY}'), gr.update(),
                gr.update())
    try:
        import json
        with open(ctx["report"], encoding="utf-8") as fh:
            report = json.load(fh)
        info: dict = {}
        frag, path = viewer.render_run(ctx["protein"], ctx["output"], report, ctx["dna"], ctx["code"],
                                       ctx["out_dir"], limits=settings.preview_limits(ctx.get("preview")),
                                       centre=max(1, col), info=info)
    except (OSError, ValueError, KeyError, TypeError) as e:
        msg = core.public_error(e, [ctx.get("protein"), ctx.get("dna"), ctx.get("out_dir")])
        return (gr.update(value=f"<p>Could not draw the alignment: {_html.escape(msg)}</p>"), gr.update(),
                gr.update())
    return (gr.update(value=frag + _view_tail(info)), gr.update(value=path, visible=True),
            gr.update(value=core.health_table_html(report, links=True, preview=info)))


def example_form(name: str) -> tuple:
    """The form values an example runs with: genetic code, remove gap columns, remove mismatched
    columns, refuse probably-wrong genes, keep only marked columns (core.EXAMPLE_SETTINGS)."""
    e = core.EXAMPLE_SETTINGS[name]
    return e["codontable"], e["nogap"], e["nomismatch"], e["strict"], e["blockonly"]


def _clear_path(_upload):
    return ""


def show_sample():
    rep = core.load_sample_report()
    banner = gui_html.sample_banner(rep, core.summary_sentence(rep))
    return _results(rep, None, core.SAMPLE_REPORT, rep.get("command", ""), rep.get("pal2nal_equivalent", ""),
                    "", banner)


# ---------------------------------------------------------------- settings

SETTING_ORDER = ["codontable", "output_format", "nogap", "nomismatch", "strict", "preview"]


def _settings_where() -> str:
    if settings.stored_in_browser():
        return "Settings are kept in this web browser only."
    return f"Settings file: `{settings.settings_path()}`"


def _apply(data: dict, status: str):
    """Outputs for the settings controls, the main form, the status line and the browser store."""
    vals = [data[k] for k in SETTING_ORDER]
    form = [data["codontable"], data["output_format"], data["nogap"], data["nomismatch"], data["strict"]]
    return (*vals, *form, status, data)


def _now() -> str:
    import time
    return time.strftime("%H:%M")


def load_settings_ui(stored=None, formats=None):
    """demo.load: the saved settings (a file locally, the browser's store on the hosted Space) into the
    Settings tab and the main form."""
    if settings.stored_in_browser():
        data, note = settings.clean(stored, formats), ""
    else:
        data, note = settings.load(formats)
    return _apply(data, note)


def save_settings_ui(code, fmt, nogap, nomismatch, strict, preview, formats=None):
    """A settings control changed: save at once (no Save button) and apply to the main form."""
    data, err = settings.save({"codontable": code, "output_format": fmt, "nogap": nogap, "nomismatch": nomismatch,
                               "strict": strict, "preview": preview}, formats)
    return _apply(data, err or f"✔ Saved at {_now()}")


def reset_settings_ui(formats=None):
    data, err = settings.reset(formats)
    return _apply(data, err or f"✔ Back to the defaults ({_now()})")


def build(demo_mode: bool = False) -> gr.Blocks:
    """Layout, top to bottom: short header; the two files side by side (path boxes under a
    "Use a file path instead" disclosure; examples as one line under them); genetic code and output format, the
    four flags under "More options"; one primary button; results: summary banner (with the genetic-code tip and
    its re-run button), downloads, health table (OK genes in one collapsed row), alignment view, and one
    "Details" area for the methods paragraph, commands, raw messages and the glossary; "Centre the preview on
    column" under the view. All of that is the first of two tabs; the second, "Settings", holds the defaults a
    user sets once (saved at once, applied to the form)."""
    features = core.detect_features()
    formats = list(features.output_formats)
    saved, _note = settings.load(formats)  # the form starts with the saved defaults (locally)
    fmt0 = saved["output_format"] if saved["output_format"] in formats else formats[0]
    with gr.Blocks(title="COAL codon alignment", theme=_theme(), css=CSS, analytics_enabled=False) as demo:
        # the header is the mark and "COAL"; the where-the-data-goes note is the page's footer
        gr.HTML(gui_html.header_html(), elem_id="coal-header", padding=False)
        local = core.path_inputs_enabled()  # path boxes only when running locally
        # two tabs: the work, then its Settings (last and on the right, as Windows apps do); "Align" is a verb,
        # like the button (the header already says COAL)
        with gr.Tabs(elem_id="coal-tabs"):
            with gr.Tab("Align", id="align", elem_id="coal-align-tab"):
                # the form reads as sections, each under one heading: Files, Genetic code, Output,
                # then the button; Results below (docs/ui-principles.md)
                gr.HTML(_sec("Files", first=True), padding=False)
                with gr.Row(equal_height=True):
                    # no fixed height (a fixed one clips the drop zone); the formats are a quiet hint inside each
                    # drop zone (CSS); the labels are not numbered (the two files are a pair, not steps)
                    protein = gr.File(label="Protein alignment", type="filepath",
                                      file_types=[".fasta", ".fa", ".faa", ".aln", ".fas", ".clustal", ".txt"],
                                      elem_id="coal-protein")
                    dna = gr.File(label="Coding DNA", type="filepath",
                                  file_types=[".fasta", ".fa", ".fna", ".ffn", ".fas", ".txt"], elem_id="coal-dna")
                # one quiet line of link-style buttons directly under the drop zones (Squoosh: "Or try one of these:")
                with gr.Row(elem_id="coal-examples"):
                    gr.Markdown("Try an example:", elem_classes=["coal-ex-label"])
                    ex_buttons = {n: gr.Button(core.EXAMPLES[n][0], size="sm", variant="secondary")
                                  for n in core.available_examples()}
                    sample_btn = gr.Button("Sample report", size="sm", variant="secondary")
                # the path boxes, in the Files section under a closed disclosure (opened by the desktop wrapper's
                # ?protein=&dna= link)
                with gr.Accordion("Use file paths (large files)", open=False, visible=local,
                                  elem_id="coal-paths") as paths_box:
                    with gr.Row():
                        with gr.Column():
                            protein_path = gr.Textbox(label="Protein alignment path",
                                                      placeholder="/path/to/protein.aln.fasta", visible=local,
                                                      max_lines=1, elem_id="coal-protein-path")
                            protein_info = gr.Markdown(visible=local)
                        with gr.Column():
                            dna_path = gr.Textbox(label="Coding DNA path", placeholder="/path/to/dna.fasta",
                                                  visible=local, max_lines=1, elem_id="coal-dna-path")
                            dna_info = gr.Markdown(visible=local)

                # short option names, the wider box for the genetic code. The genetic code is the choice that
                # most often makes results wrong: it is a section of its own, first, its label drawn as the heading
                # (CSS), a larger value and one hint line. Output format is secondary: a narrower box under
                # "Output" (its label is the heading; the box keeps it as its accessible name)
                code = gr.Dropdown(_code_choices(), value=saved["codontable"], label="Genetic code", info=CODE_HINT,
                                   min_width=260, elem_id="coal-code")
                gr.HTML(_sec("Output"), padding=False)
                fmt = gr.Dropdown(_format_choices(features), value=fmt0, label="Output format", show_label=False,
                                  info=_format_info(features), min_width=220, elem_id="coal-format")
                # the four flags under "More options", stacked, plain label first, the flag at the end of the
                # one-line help text
                with gr.Accordion("More options", open=False, elem_id="coal-more"):
                    nogap = gr.Checkbox(value=saved["nogap"], label="Remove gap columns",
                                        info="Columns with a gap or stop codon. (-nogap)")
                    nomismatch = gr.Checkbox(value=saved["nomismatch"], label="Remove mismatched columns",
                                             info="Columns where a codon doesn't match its amino acid. (-nomismatch)")
                    strict = gr.Checkbox(value=saved["strict"], label="Refuse probably-wrong genes",
                                         info=STRICT_INFO)
                    blockonly = gr.Checkbox(label="Keep only '#'-marked columns",
                                            info="CLUSTAL input marked by Gblocks or similar. (-blockonly)")
                run_btn = gr.Button("Start alignment", variant="primary", size="lg", elem_id="coal-run")

                # hidden until there is a result; every result (run, example, sample report, missing-file
                # message) shows it (_results)
                with gr.Column(elem_id="coal-results", visible=bool(demo_mode)) as results_box:
                    gr.HTML(_sec("Results"), padding=False)
                    # banner, tip and re-run button in one column; when the genetic-code hint
                    # explains the run, CSS draws them as one box with the button as its primary action
                    with gr.Column(elem_id="coal-fixbox"):
                        banner = gr.HTML(gui_html.empty_banner(), elem_id="coal-banner", padding=False)
                        hint_md = gr.Markdown(visible=False, elem_id="coal-hint")
                        hint_btn = gr.Button("Re-run", visible=False, variant="secondary", size="md",
                                             elem_id="coal-hint-rerun")
                    hint_code = gr.State(None)
                    with gr.Row(elem_id="coal-downloads"):
                        # size md: the default lg is as big as the primary button
                        dl_out = gr.DownloadButton("Download codon alignment", visible=False, variant="secondary",
                                                   size="md")
                        dl_rep = gr.DownloadButton("Download report (JSON)", visible=False, variant="secondary",
                                                   size="md")
                    table = gr.HTML(elem_id="coal-health-box", padding=False)
                    with gr.Accordion("View alignment", open=True, elem_id="coal-view"):
                        view = gr.HTML(VIEW_EMPTY, padding=False)
                        # shown when the preview leaves columns out
                        with gr.Row(visible=False, elem_id="coal-centre") as centre_row:
                            centre_col = gr.Number(label="Centre the preview on column", precision=0, minimum=1,
                                                   elem_id="coal-centre-col")
                            centre_btn = gr.Button("Show these columns", variant="secondary",
                                                   elem_id="coal-centre-go")
                        dl_view = gr.DownloadButton("Download alignment view (HTML)", visible=False,
                                                    variant="secondary", size="md", elem_id="coal-view-dl")
                    # methods, commands, raw messages and the glossary in one area
                    with gr.Accordion("Methods, commands and glossary", open=False, elem_id="coal-details"):
                        methods = gr.Textbox(label="Methods paragraph", lines=4,
                                             max_lines=8, interactive=False, show_copy_button=True, visible=False)
                        cmds = gr.Code(label="Equivalent commands", language="shell", interactive=False)
                        raw = gr.Textbox(label="Raw COAL messages", lines=6, max_lines=20, interactive=False)
                        gr.Markdown("**Glossary**\n" + GLOSSARY)

            # Settings. Few settings (only what people set once), plain words, saved as soon as
            # they change (no Save button), "Reset to defaults", and where they are kept. See docs/ui-principles.md
            # No intro paragraph, no help lines; one-word group headers (Microsoft: "Accounts", not "Account
            # settings"); one "About" line
            with gr.Tab("Settings", id="settings", elem_id="coal-settings-tab"):
                with gr.Column(elem_id="coal-settings"):
                    gr.HTML(_sec("Defaults", first=True), padding=False)
                    with gr.Row():
                        s_code = gr.Dropdown(_code_choices(), value=saved["codontable"], label="Genetic code",
                                             scale=3, min_width=260, elem_id="coal-set-code")
                        s_fmt = gr.Dropdown(_format_choices(features), value=fmt0, label="Output format",
                                            info=_format_info(features), scale=2, min_width=220,
                                            elem_id="coal-set-format")
                    s_nogap = gr.Checkbox(value=saved["nogap"], label="Remove gap columns")
                    s_nomismatch = gr.Checkbox(value=saved["nomismatch"], label="Remove mismatched columns")
                    s_strict = gr.Checkbox(value=saved["strict"], label="Refuse probably-wrong genes")
                    gr.HTML(_sec("Preview"), padding=False)
                    s_preview = gr.Radio([(settings.PREVIEW_LABELS[k], k) for k in settings.PREVIEW_SIZES],
                                         value=saved["preview"], label="Largest preview (sequences × columns)",
                                         elem_id="coal-set-preview")
                    gr.HTML(_sec("About"), padding=False)
                    gr.Markdown(
                        _settings_where() + f"\n\nResults are deleted after {core.RUNS_MAX_AGE_HOURS} hours"
                        + (f" (`{core.RUNS_DIR}`)." if local else "."), elem_classes=["coal-small"])
                    s_status = gr.Markdown(_note, elem_id="coal-settings-status")
                    s_reset = gr.Button("Reset to defaults", variant="secondary", size="sm", elem_id="coal-set-reset")
        # the where-the-data-goes note, small, at the foot of the page, only when there is one (the hosted Space)
        if core.intro_text():
            gr.HTML(gui_html.footer_html(core.intro_text()), elem_id="coal-footer", padding=False)
        # per-browser store: used on the hosted Space only (locally the settings file is the store, because the
        # desktop app's web view keeps no browser storage between starts). The storage key and the "secret" are
        # fixed so the value survives a restart of the Space (Gradio's guide: a random secret is made per launch
        # otherwise); the data is not sensitive, the secret only obfuscates it.
        if hasattr(gr, "BrowserState"):
            s_store = gr.BrowserState(dict(settings.DEFAULTS), storage_key="coal-settings-v1",
                                      secret="coal-settings-v1-not-a-secret")
        else:
            s_store = gr.State(dict(settings.DEFAULTS))

        outputs = [banner, hint_md, hint_btn, hint_code, dl_out, dl_rep, table, cmds, raw, view, dl_view, methods]
        # the run context (for "Centre the preview on column") and that control's row
        run_ctx = gr.State(None)
        outputs += [run_ctx, centre_row]
        outputs += [results_box]
        inputs = [protein, dna, fmt, code, nogap, nomismatch, protein_path, dna_path, strict, blockonly, s_preview]
        run_btn.click(run_ui, inputs, outputs, show_api=False)
        protein_path.change(core.path_status, protein_path, protein_info, show_api=False)
        dna_path.change(core.path_status, dna_path, dna_info, show_api=False)
        # uploading a file clears its path box: a path wins over an upload, and "Open with" fills the boxes
        protein.upload(_clear_path, protein, protein_path, show_api=False)
        dna.upload(_clear_path, dna, dna_path, show_api=False)
        centre_btn.click(centre_ui, [run_ctx, centre_col], [view, dl_view, table], show_api=False)
        centre_col.submit(centre_ui, [run_ctx, centre_col], [view, dl_view, table], show_api=False)

        # Settings: "input" events only (the user's own changes), so loading or resetting them
        # does not save again. Each change is saved and applied to the main form at once.
        s_controls = [s_code, s_fmt, s_nogap, s_nomismatch, s_strict, s_preview]
        s_outputs = s_controls + [code, fmt, nogap, nomismatch, strict, s_status, s_store]

        def _save(c, f, g, m, st, pv):
            return save_settings_ui(c, f, g, m, st, pv, formats)

        for comp in s_controls:
            comp.input(_save, s_controls, s_outputs, show_api=False)
        s_reset.click(lambda: reset_settings_ui(formats), None, s_outputs, show_api=False)
        demo.load(lambda stored: load_settings_ui(stored, formats), s_store, s_outputs, show_api=False)

        if local:
            # Desktop wrappers: open http://127.0.0.1:7860/?protein=<path>&dna=<path> (URL-encoded) to fill the
            # path boxes.
            def fill_from_url(request: gr.Request):
                params = dict(request.query_params) if request is not None else {}
                p, d = core.paths_from_query(params)
                # also opens the path disclosure when the link filled a box
                return p, d, gr.update(open=bool(p or d))

            demo.load(fill_from_url, None, [protein_path, dna_path, paths_box], show_api=False)

        # an example also sets the genetic code and the four options to its own settings (core.EXAMPLE_SETTINGS),
        # shown in the form, before it runs, so it never inherits a code from an earlier run. A user's own files
        # keep the form as is.
        for name, btn in ex_buttons.items():
            p, d = core.example_paths(name)

            def load_example(p=p, d=d, name=name):
                return (p, d, "", "", *example_form(name))

            btn.click(load_example, None, [protein, dna, protein_path, dna_path, code, nogap, nomismatch, strict,
                                           blockonly], show_api=False).then(run_ui, inputs, outputs, show_api=False)

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
    core.prune_runs()  # runs older than 24 h are removed at startup
    demo = build(demo_mode)
    demo.launch(server_name=host, server_port=port, mcp_server=mcp, inbrowser=False,
                allowed_paths=[core.EXAMPLES_DIR, core.RUNS_DIR, os.path.dirname(core.SAMPLE_REPORT)],
                prevent_thread_lock=prevent_thread_lock, show_error=True)
    return demo
