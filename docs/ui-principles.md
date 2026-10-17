# COAL app: usability principles and what they mean here

Written for GUI round 2 (overnight/integration/CHANGES.md, rows GUI-8 to GUI-12). The owner's rule comes first: if
the experience isn't clean, then for a user the tool does not work. The recommendations below come from published
guidance and are listed in priority order. "Done" means it is in the app now.

## Recommendations

| # | Recommendation | Source | Status |
|---|---|---|---|
| 1 | At scale, say plainly that the grid is a preview: how many sequences and which columns are drawn, and that the downloaded files are complete. | NN/g heuristic 1, visibility of system status | Done (GUI-10): "Preview: 100 of 20,000 sequences, problems first (…); columns 1–1,000 of 12,000." |
| 2 | Show an overview of the whole alignment first, then the detail. Mark the part on screen with a box, as Jalview's overview window does. | Shneiderman 1996; Jalview overview window | Done (GUI-10): a sequence strip and a column strip (one small SVG each, from totals), with a box over the drawn columns. |
| 3 | Pick the rows to draw by what the user is there for: every problem first, then OK sequences from across the file for comparison. | Ten simple rules: be ready for data growth (rule 3); NN/g heuristic 8 | Done (GUI-10) |
| 4 | Send people to dedicated tools for the full alignment instead of drawing 100,000 rows in a web page. Use standard formats. | Ten simple rules, rule 4 (standard formats); AliView and Jalview are built for large alignments | Done (GUI-10): "open the codon alignment in Jalview or AliView" in the preview note and under the view. |
| 5 | Keep settings few. Use them only for things set once (defaults), and keep everyday choices on the main form. | Microsoft Windows app settings guidelines; Apple HIG (Settings) | Done (GUI-9): five defaults and a preview size. Genetic code and format stay on the form. |
| 6 | Apply a setting as soon as it changes, with no Save button. Confirm that it was saved, offer "Reset to defaults", and say where it is stored. | Microsoft settings guidelines; NN/g heuristic 1 and 3 (user control: a way back) | Done (GUI-9): "✔ Saved at 14:02 …", Reset to defaults, the path of settings.json (or "this browser" on the hosted copy). |
| 7 | Put Settings last, apart from the work, and label its groups with short headers. Don't nest disclosures more than two levels deep. | Microsoft settings guidelines; NN/g progressive disclosure | Done (GUI-9): a "Settings" tab on the right, with three short groups. |
| 8 | Name options in the user's words, short enough to read whole, with the explanation on the help line. | NN/g heuristics 2 and 6 (recognition rather than recall) | Done (GUI-8): "1: Standard (nuclear genes)". GUI-13: the option names carry the meaning, so the help lines are gone. |
| 9 | Never link somewhere that doesn't exist. A jump to a column outside the preview says so and tells you how to get there. | NN/g heuristics 1 and 9 | Done (GUI-10): no dead links; "Centre the preview on column" redraws around any column. |
| 10 | Keep long lists from pushing the main result off screen. Show the first few and put the rest behind a labelled disclosure. | NN/g progressive disclosure; heuristic 8 | Done (GUI-10): 8 problems are open above the grid, the rest under "N more …". |
| 11 | Expose only the parameters that matter, with conservative defaults. | Ten simple rules, rule 5 | Kept: advanced flags stay under "More options". Keep-going is always on and is described, not offered as a setting. |
| 12 | Errors and summaries: one summary at the top that links to each problem, and each message placed next to the thing it is about. | GOV.UK Design System error summary; NN/g error-message guidelines | Partly done (banner, then the problems-first health table with jump links). Not done: move keyboard focus to the banner after a run, and show "add this file" next to the empty drop zone. Gradio 5 has no supported focus or per-field error API, so both need custom JavaScript. |
| 13 | For very large views, draw only the rows on screen (virtual scrolling) and let a click on the overview move the view. | Jalview overview (click to centre); AliView (large alignments) | Not done: gr.HTML does not run scripts, so this needs a Gradio custom component. The column box plus "Centre the preview on column" covers the main need. |
| 14 | Get users started quickly with demo data. | Ten simple rules, rule 8 | Already there: the examples row and the demo sample. |
| 15 | Test with real users early, and watch where they hesitate. | Ten simple rules, rule 2; NN/g | Open: the owner's Windows run is the first round. A short scripted session with two biologists is the next step. |

## Sources

- Nielsen Norman Group, "10 Usability Heuristics for User Interface Design": https://www.nngroup.com/articles/ten-usability-heuristics/
- Nielsen Norman Group, "Progressive Disclosure": https://www.nngroup.com/articles/progressive-disclosure/
- Nielsen Norman Group, "Error-Message Guidelines": https://www.nngroup.com/articles/error-message-guidelines/
- GOV.UK Design System, "Error summary": https://design-system.service.gov.uk/components/error-summary/
- Microsoft, "Guidelines for app settings" (Windows apps): https://learn.microsoft.com/en-us/windows/apps/design/app-settings/guidelines-for-app-settings
- Apple Human Interface Guidelines, "Settings": https://developer.apple.com/design/human-interface-guidelines/settings (the page needs JavaScript, so only its summary could be read here)
- List M, Ebert P, Albrecht F (2017), "Ten simple rules for developing usable software in computational biology", PLOS Comput Biol 13(1): e1005265, https://doi.org/10.1371/journal.pcbi.1005265 (read via https://limswiki.org/index.php/Journal:Ten_simple_rules_for_developing_usable_software_in_computational_biology)
- Shneiderman B (1996), "The Eyes Have It: A Task by Data Type Taxonomy for Information Visualizations": https://drum.lib.umd.edu/handle/1903/5784
- Jalview help, "Overview window": https://www.jalview.org/help/html/features/overview.html
- Larsson A (2014), "AliView: a fast and lightweight alignment viewer and editor for large datasets", Bioinformatics, doi:10.1093/bioinformatics/btu531
- NCBI Insights, "Visualize and interpret alignment data with the Multiple Sequence Alignment Viewer" (2017): https://ncbiinsights.ncbi.nlm.nih.gov/2017/01/25/visualize-and-interpret-alignment-data-with-the-multiple-sequence-alignment-viewer/
- Gradio 5.49.1 source, `gradio/components/browser_state.py`, and the guide "State in Blocks" (BrowserState: the value lives in the browser's localStorage, encrypted with `secret`; a fixed `storage_key` and `secret` are needed for it to survive a restart): https://www.gradio.app/guides/state-in-blocks

## Round 3 baseline

Written for GUI round 3 (rows GUI-13 and GUI-14 in overnight/integration/CHANGES.md). The owner's feedback: "visual
clutter … the descriptions are too long — shorten or get rid of them, use existing software to make the decisions."
Six well-regarded tools are the baseline. For each: what COAL takes from it, and what to avoid.

| Tool | Convention taken | Applied in COAL |
|---|---|---|
| [Squoosh](https://squoosh.app/) (Google's image compressor: drop a file, pick options, get the result) | The start screen is the name and one drop target ("Drop OR Paste"), then "Or try one of these:" with sample files. The privacy promise ("Images never leave your device") sits below the work, not in the header. | Header is the mark and "COAL". "Try an example:" is one line of links under the drop zones. "Nothing leaves your computer" is a small footer line. |
| Microsoft Windows: [writing style](https://learn.microsoft.com/windows/uwp/design/style/writing-style) and [app settings guidelines](https://learn.microsoft.com/en-us/windows/apps/design/app-settings/guidelines-for-app-settings) | "Use no more words than you need"; labels without a closing period; a setting's description is optional, "add it only when a label alone isn't clear", and is used to explain why a setting is unavailable; one-word section headers ("Accounts", not "Account settings"); app information goes in an About section at the end of Settings. | No help lines under the dropdowns. The output-format box has a line only when the engine writes FASTA alone. Settings groups are "Defaults", "Preview" and "About"; the settings-file path and the 24-hour clean-up are the About line. |
| Apple HIG, [Writing](https://developer.apple.com/design/human-interface-guidelines/writing) | "Check each word to be sure it needs to be there. If you can use fewer words, do so." Describe what a setting does when it is on, and let people infer the opposite. Empty states: one clear next step. Buttons: a verb. | Checkbox help says what happens when ticked, in one line ("Columns with a gap or stop codon. (-nogap)"). The empty result area says "Results appear here."; the button above it already names the action. Tab "Align" (a verb), not "Codon alignment". |
| GitHub [Primer forms](https://primer.style/product/ui-patterns/forms/) (and GitHub Desktop, which uses it) | Labels of at most three words, in sentence case; a caption only when the label is not enough, "If the caption feels redundant, try removing it"; never a placeholder instead of a label; stack fields vertically by default. | "Protein alignment", "Coding DNA", "Genetic code", "Output format". The accepted formats are a quiet hint inside each drop zone, not part of the label. The four option checkboxes are stacked. |
| [Linear's UI redesign](https://linear.app/now/how-we-redesigned-the-linear-ui) | "Reduce visual noise, maintain visual alignment": align labels, icons and buttons on one edge; a neutral palette with the accent used sparingly. | No panel or border around each field (boxes within boxes): fields sit on the page and every label starts on the same left edge; disclosures are separated by one hairline. Ember stays on the one primary button and the selected tab. |
| [HandBrake](https://handbrake.fr/docs/en/1.7.0/workflow/open-video-source.html) (drop-in desktop tool) | Open a source by dropping a file on the window or with one button; after that, the window shows the file's name, nothing more. | Once a file is in, the drop zone shows its name (Gradio's file row); the run banner names both files once ("Files used: …"). |

Bioinformatics tools, for what to copy and what to avoid:

- [AliView](https://pmc.ncbi.nlm.nih.gov/articles/PMC4221126) (Larsson 2014) was built "with extensive focus on usability and intuitive handling" by following the conventions of software people already know (text editors, browsers). COAL keeps standard controls (drop zone, dropdowns, one primary button) rather than inventing its own.
- [Jalview's overview window](https://www.jalview.org/help/html/features/overview.html): the whole alignment small, with a box over what is on screen. Kept from round 2 (GUI-10).
- Galaxy tool forms, to avoid: users in a course were "overwhelmed by the choices/parameters in tools" ([Galaxy help forum, "Complexity of tools and user feedback"](https://help.galaxyproject.org/t/complexity-of-tools-and-user-feedback/13855)), and long parameter lists with the help text far below the form were raised on the developers' list ([galaxy-dev, "TopHat and other tools with too many options"](https://lists.galaxyproject.org/archives/list/galaxy-dev@lists.galaxyproject.org/thread/NBH7MA46FNPDS3EMKUIJVWHCJIPG2KGJ/)). COAL shows two choices; the four flags stay closed under "More options", each with one line next to it.

Rules that follow, for any later change: a new line of text on the main form needs a reason a label cannot carry;
numbers ("1.", "2.") only for real sequences of steps; one sentence per fact on screen (the banner, the table and
the alignment view must not say the same thing twice); explanations belong in the glossary, the docs or the About
line, not under a field.

## Round 4 polish

Written for v0.4.1 (rows POL-1 to POL-6 in overnight/integration/CHANGES.md). The owner, after running the desktop
app on Windows: cut the reassurance lines, give each section a clear hierarchy as on the website, emphasise the
genetic code, and replace "Results appear here." Three reference points, and what was taken from each:

| Reference | Convention | Applied in COAL |
|---|---|---|
| [Squoosh](https://squoosh.app/) | The options panel is grouped under short section headings ("Edit", "Compress"), with lighter option labels under them; there is no result panel until an image is open. | Headings "Files", "Genetic code", "Output", "Results" (h2, the site's h3 size), field labels one step lighter. The results area (heading, banner, view, methods) is hidden until there is a result. |
| GitHub [Primer forms](https://primer.style/product/ui-patterns/forms/) | Related fields grouped under one heading; field width fits the expected value; no caption unless the label is not enough. | One heading per group; the genetic-code box is at most 640 px and the output-format box 320 px wide; one hint line on the genetic code only. |
| MEGA (Kumar et al.) and other codon-aware tools | The genetic code is a first-class choice for coding data (MEGA asks for it when coding sequences are opened), not a detail next to the file format. | The genetic code is its own section, first after the files, its label drawn as the heading, a larger value, and one line: "Mitochondrial, bacterial and plastid genes need their own code." Output format is secondary, under "Output". |

Also: NN/g, [empty states in complex applications](https://www.nngroup.com/articles/empty-state-interface-design/):
an empty state should help people act. Here the form above already names the next step (the button), so a
placeholder under it added nothing; the area is simply absent until it has content.

Removed: the splash line "This computer only. No internet needed." and the footer "Nothing leaves your computer:
COAL runs here, on this machine." (reassurance, not something to act on; the facts stay in README.md,
desktop/README.md and on the site; the hosted Space keeps its note, because there files do go to a server). Gradio's
own footer link "Settings" (display language and theme for Gradio's labels only) is hidden, so COAL has one
settings place, its Settings tab; "Use via API or MCP" stays. The round 3 rows above that mention the footer note
and "Results appear here." describe that round.
