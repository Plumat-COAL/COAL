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
| 8 | Name options in the user's words, short enough to read whole, with the explanation on the help line. | NN/g heuristics 2 and 6 (recognition rather than recall) | Done (GUI-8): "1: Standard (nuclear genes)"; the help line lists the codes people need most (1, 2, 11). |
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
